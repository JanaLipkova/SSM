/*
 *  DelayTauLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "DelayTauLeaping.h"

DelayTauLeaping::DelayTauLeaping(Simulation * simulation):
LeapMethod(simulation)
{
}

DelayTauLeaping::~DelayTauLeaping()
{
}

void DelayTauLeaping::solve()
{
	cout << "Delay-TauLeaping..." << endl;

/*	for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
	{
		SSMReaction * reaction = simulation->ssmReactionList[i];;
		reaction->toString();
	}
exit(0);*/

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
	
	vector < double > timeSteps;
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations	= 0;
		timePoint			= 0;
		simulation->loadInitialConditions();

		delayedReactionsTimePoints.clear();
		delayedReactionsTimeIndices.clear();
		delayedReactionsTimeLeapLength.clear();
		
		while (t < tEnd)
		{
			saveData();
			
			// compute the propensities
			computePropensities(propensitiesVector);
			//cout << "r4: " << propensitiesVector(3) << endl; exit(0);
			a0 = blitz::sum(propensitiesVector);
			//cout << "	a0: " << a0 << endl;
			dt = computeTimeStep();
			
			/*cout << "tau: " << dt << endl;
			cout << "a0: " << a0 << endl;
			cout << "r: " << 10.0 / a0 << endl;
			cout << "bool: " << (dt <= (10.0 / a0)) << endl << endl; */
			
			if (dt <= (5.0/a0))
			{
				cout << "dt is too small: " << endl;
				cout << "dt: " << dt << endl;
				cout << 5.0/a0 << endl;
				//return;
			}
									
			numberToFireForEveryChannel = 0;
			delayIndex = 0;
			
			while ( isDelayedReactionScheduled(t, dt, delayIndex) == true ) // fire delay that is scheduled within this timestep
			{
				/*
				cout << "delayIndex: " << delayIndex << endl; 
				cout << "delayedReactionsTimeIndices[delayIndex]: " <<  delayedReactionsTimeIndices[delayIndex] << endl; 
				cout << "delayedReactionsIndices( delayedReactionsTimeIndices[delayIndex] ): " << delayedReactionsIndices( delayedReactionsTimeIndices[delayIndex] ); exit(0);
				*/
				
				fireReactionProposed( delayedReactionsTimeIndices[delayIndex] , delayedReactionsTimeLeapLength[delayIndex] );	
				//fireReactionProposed( delayedReactionsIndices( delayedReactionsTimeIndices[delayIndex] ), delayedReactionsTimeLeapLength[delayIndex] );				
				// always fire the delayed reactions (nonconsuming reactions)
				acceptNewSpeciesValues();
				removeDelayedReactionsTime(delayIndex);
			}
			
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
					setDelayedReactionsTime(reactionIndex, t, 0.5*dt,  numberToFireForEveryChannel(reactionIndex) );
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
		cout << "Average tau: " << tauAverage / ((double)numberOfIterations) << endl;
		
		writeSTLVector( timeSteps, sbmlModel->getName() + "time-steps.txt");
		
		
	}
	writeData(outputFileName);
}


