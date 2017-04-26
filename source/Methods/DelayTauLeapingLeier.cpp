/*
 *  DelayTauLeapingLeier.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "DelayTauLeapingLeier.h"

DelayTauLeapingLeier::DelayTauLeapingLeier(Simulation * simulation):
LeapMethod(simulation)
{
}

DelayTauLeapingLeier::~DelayTauLeapingLeier()
{
}

void DelayTauLeapingLeier::solve()
{

	cout << "Delay-TauLeaping-LEIER..." << endl;

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
			a0 = blitz::sum(propensitiesVector);
			dt = computeTimeStep();
						
			numberToFireForEveryChannel = 0;			
			delayIndex = 0;
			
			while ( isDelayedReactionScheduled(t, dt, delayIndex) == true ) // fire delay that is scheduled within this timestep
			{
				fireReactionProposed( delayedReactionsTimeIndices[delayIndex], delayedReactionsTimeLeapLength[delayIndex] );
				//fireReactionProposed( delayedReactionsIndices( delayedReactionsTimeIndices[delayIndex] ), delayedReactionsTimeLeapLength[delayIndex] );				
				// always fire the delayed reactions (nonconsuming reactions)
				acceptNewSpeciesValues();
				removeDelayedReactionsTime(delayIndex);
			}
			
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
			{				
				numberToFireForEveryChannel( j ) = ignpoi( propensitiesVector(j)*dt );				
			}
			
			for (int reactionIndex = 0; reactionIndex < numberToFireForEveryChannel.extent(firstDim); ++reactionIndex)
			{
				if ( isReactionDelayed(reactionIndex) == false ) // not delayed
				{
					fireReactionProposed(reactionIndex, numberToFireForEveryChannel(reactionIndex));
				}
				else if ( numberToFireForEveryChannel(reactionIndex) > 0 ) // queue the delay
				{
					for (int i = 0; i < numberToFireForEveryChannel(reactionIndex); ++i)
					{
						setDelayedReactionsTime(reactionIndex, t, ranf()*dt, 1 );
					}
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
		}
		saveData();
		// statistics
//		writeData(outputFileName, samples);
//		zeroData();
		cout << "		sample: " << samples << endl;		
//		cout << "Sample: " << samples << endl;
//		cout << "Average tau: " << tauAverage / ((double)numberOfIterations) << endl;
	}
	writeData(outputFileName);
}
