/*
 *  DelaySSA.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/6/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "DelaySSA.h"

DelaySSA::DelaySSA(Simulation * simulation):
Method(simulation)
{
}

DelaySSA::~DelaySSA()
{
}

double DelaySSA::randomVariate(double a0)
{
	//return (genexp( (1.0 / a0) ));
	
	// select Exponential probability distribution
    boost::exponential_distribution<double> exp_dist(a0);
 
    // bind random number generator to distribution, forming a function
    boost::variate_generator<boost::mt19937&, boost::exponential_distribution<double> >  exponential_sampler(rng, exp_dist);
	return exponential_sampler();

}

void DelaySSA::solve()
{
	cout << "Delay-SSA-Cai-..." << endl;
	
	int numDelayed = numberOfDelayedReactions();
	if (numDelayed == 0)
	{
		cout << "You are using the Delay-SSA method but have not specified a reaction that has a delay." << endl;
		cout << "I shall now quit and give you some time to think about this." << endl; 
		return;
	}

	delayedReactionsIndices.resize( numberOfDelayedReactions() );
	
	setDelayedReactionsIndices();

	double a0 = 0.0;
	double r1;
	int reactionIndex = 0;
	double cummulative = 0.0;
	int delayIndex = 0;
	double tauAverage = 0.0;
	
	openAuxiliaryStream( (simulation->ModelName) + "-histogram-dssa.txt");
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		timePoint = 0;
		simulation->loadInitialConditions();

		delayedReactionsTimePoints.clear();
		delayedReactionsTimeIndices.clear();
		delayedReactionsTimeLeapLength.clear();
		
		while (t < tEnd)
		{ 
			saveData();
			
			delayIndex = 0;
			computePropensities(propensitiesVector);
			a0 = blitz::sum(propensitiesVector);
			dt = randomVariate(a0);
				
			while ( isDelayedReactionScheduled(t, dt, delayIndex) == true ) // fire delay that is scheduled within this timestep
			{
				//cout << "drti: " << delayedReactionsTimeIndices[delayIndex] << endl;
				//fireReaction( delayedReactionsIndices( delayedReactionsTimeIndices[delayIndex] ), 1 );
				fireReaction( delayedReactionsTimeIndices[delayIndex] , delayedReactionsTimeLeapLength[delayIndex] );
				t = delayedReactionsTimePoints[delayIndex];	
				removeDelayedReactionsTime(delayIndex);
				
				// rejection method
				computePropensities(propensitiesVector);
				a0 = blitz::sum(propensitiesVector);
				dt = randomVariate(a0);
			}
			
			r1 = ranf();
			reactionIndex = 0;
			cummulative = 0.0;
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
			{
				cummulative += propensitiesVector(j);
				if ( cummulative > a0*r1 )
				{
					reactionIndex = j;
					break;
				}
			}
			
			if ( isReactionDelayed(reactionIndex) == false ) // not delayed
			{
				fireReaction(reactionIndex, 1);
			}
			else // queue the delay
			{
				setDelayedReactionsTime(reactionIndex, t, dt, 1);
			}
			
			++numberOfIterations;
			t += dt;
			tauAverage += dt;
		}
		
//		saveData();
		// statistics
//		writeData(outputFileName, samples);
//		zeroData();
//		cout << "		s: " << samples << endl;
				
//		cout << "Sample: " << samples << endl;
		saveData();
		cout << "Sample: " << samples << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
		cout << "Average tau: " << tauAverage / ((double)numberOfIterations) << endl;
	}
	writeData(outputFileName);
	closeAuxiliaryStream();
}

