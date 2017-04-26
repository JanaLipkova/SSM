/*
 *  DelayTauLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "DelayTauLeapingPhilippe.h"

DelayTauLeapingPhilippe::DelayTauLeapingPhilippe(Simulation * simulation):
LeapMethod(simulation)
{
}

DelayTauLeapingPhilippe::~DelayTauLeapingPhilippe()
{
}

void DelayTauLeapingPhilippe::solve()
{
	cout << "Delay-TauLeaping-Philippe..." << endl;

	int numDelayed = numberOfDelayedReactions();
	if (numDelayed == 0)
	{
		cout << "You are using the Delay-TauLeaping method but have not specified a reaction that has a delay." << endl;
		cout << "I shall now quit and give you some time to think about this." << endl; 
		return;
	}

	delayedReactionsIndices.resize( numberOfDelayedReactions() );
	
	setDelayedReactionsIndices();

	double a0			= 0.0;
	int delayIndex		= 0;
	
	Array<int, 1> numberToFireForEveryChannel(sbmlModel->getNumReactions());
	numberToFireForEveryChannel = 0;
	
	double tauAverage = 0.0;
	
	vector <double> queueTau;
	
	openAuxiliaryStream( (simulation->ModelName) + "-histogram-DTauLeaping.txt");
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations	= 0;
		timePoint			= 0;
		simulation->loadInitialConditions();

		delayedReactionsTimePoints.clear();
		delayedReactionsTimeIndices.clear();
		delayedReactionsTimeLeapLength.clear();
		
		double minValueDTDQ = 1.0;
		vector <double> meanValueDTDQ;
		double maxValueDTDQ = 0.0;
		queueTau.clear();
		
		vector <double> statistics;
		
		vector < double > timeSteps;
		
		while (t < tEnd)
		{
			saveData();
			
			// compute the propensities
			computePropensities(propensitiesVector);
	//exit(0);		
			a0 = blitz::sum(propensitiesVector);
			dt = computeTimeStep();
			
			if (dt <= (2.0/a0))
			{
				cout << "dt is too small: " << endl;
				cout << "dt: " << dt << endl;
				cout << 2.0/a0 << endl;
				//return;
			}
									
			numberToFireForEveryChannel = 0;
			delayIndex = 0;
			
			double stat = 0.0;
			while ( isDelayedReactionScheduled(t, dt, delayIndex) == true ) // fire delay that is scheduled within this timestep
			{
				double qt			= delayedReactionsTimePoints[delayIndex];
				double span			= queueTau[delayIndex];
				long long int kd	= delayedReactionsTimeLeapLength[delayIndex];
				
				double ratio = (t + dt - qt) / span;
				meanValueDTDQ.push_back( ratio );
				
				
				minValueDTDQ = min( minValueDTDQ, ratio );
				maxValueDTDQ = max( maxValueDTDQ, ratio );
				
				long int partialAmountFired = ignbin( kd,  min( ratio, 1.0 ) );
				
				stat += ((double)partialAmountFired);
				//stat += ( (((double)partialAmountFired) / ((double)kd) ) );
				
				fireReactionProposed( delayedReactionsTimeIndices[delayIndex] , partialAmountFired );
				acceptNewSpeciesValues();
				
				// do not remove the delayed reaction, but do give it a new queued dt and K
				queueTau[delayIndex]					   -= (t + dt - qt);     	// - update span
				delayedReactionsTimeLeapLength[delayIndex] -= partialAmountFired;  // - update kd
				delayedReactionsTimePoints[delayIndex]		= t + dt; // it's still up - update qt
				
				if (delayedReactionsTimeLeapLength[delayIndex] == 0)
				{
					removeDelayedReactionsTime(delayIndex);
					queueTau.erase( queueTau.begin() + delayIndex, queueTau.begin()+delayIndex+1);
				}
				
			}
			statistics.push_back( stat );
			
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
			{				
				numberToFireForEveryChannel( j )  = ignpoi( propensitiesVector(j)*dt );
			}
			
			for (int reactionIndex = 0; reactionIndex < numberToFireForEveryChannel.extent(firstDim); ++reactionIndex)
			{
				if ( isReactionDelayed(reactionIndex) == false ) // not delayed
				{
					fireReactionProposed(reactionIndex, numberToFireForEveryChannel(reactionIndex));
				}
				else if ( numberToFireForEveryChannel(reactionIndex) > 0 ) // queue the delay
				{
					setDelayedReactionsTime(reactionIndex, t, 0.0*dt,  numberToFireForEveryChannel(reactionIndex) );
					queueTau.push_back(dt);
				}
				else {}
			}
			
			if (isProposedNegative() == false)
			{
				acceptNewSpeciesValues();
				++numberOfIterations;
				
				t += dt;
				tauAverage += dt;
			}
			else
			{
				cout << "Negative species at time: " << t << endl;
				reloadProposedSpeciesValues();
			}
			timeSteps.push_back(dt);
		}

		saveData();
		cout << "Sample: " << samples << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
		
		/*
//		saveData();
		// statistics
//		writeData(outputFileName, samples);
//		zeroData();
//		cout << "		sample: " << samples << endl;	
		
//		writeSTLVector( statistics, sbmlModel->getName() + "time-steps-stat.txt");
//		writeSTLVector( timeSteps, sbmlModel->getName() + "time-steps-tau.txt");
		*/
		
/*		cout << "minValueDTDQ: " << minValueDTDQ << endl << endl;
		cout << "maxValueDTDQ: " << maxValueDTDQ << endl << endl;
		double mv = 0.0;
		for (int i = 0; i < meanValueDTDQ.size(); ++i)
		{
			mv += meanValueDTDQ[i];
		}
		mv /= meanValueDTDQ.size();
		cout << "meanValueDTDQ: " << mv << endl << endl;
		
		cout << "Sample: " << samples << endl;
		cout << "Average tau: " << tauAverage / ((double)numberOfIterations) << endl;*/
	}
	writeData(outputFileName);
	closeAuxiliaryStream();
}
