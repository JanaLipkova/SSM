/*
 *  TauLeapingCTS.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Roger Rossé on 02/12/12.
 *  Copyright 2012 Roger Rossé. All rights reserved.
 *
 */

#include "TauLeapingCT.h"

TauLeapingCT::TauLeapingCT(Simulation * simulation):
LeapMethod(simulation)
{
}

TauLeapingCT::~TauLeapingCT()
{
}

double TauLeapingCT::computeTimeStep()
{
	double tauPrime = simulation->Tau;
	return tauPrime;
}

void TauLeapingCT::solve()
{
	cout << "TauLeapingCT..." << endl;
	
	double aj;
	long int kj;
	double a0 = 0.0;
	double dt;
	bool isNegative = false;
	
	openAuxiliaryStream( (simulation->ModelName) + "-histogram-tauLeaping.txt");
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		timePoint = 0;
		simulation->loadInitialConditions();
		
		while (t < tEnd)
		{
			saveData();
			computePropensities(propensitiesVector, 0);
			a0 = blitz::sum(propensitiesVector);
	
			if (isNegative == false)
			{
				dt = computeTimeStep();
			}
			
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
			{				
				aj				= propensitiesVector(j);
				kj				= ignpoi( aj*dt );
				
				fireReactionProposed( j , kj );
			}

			if (isProposedNegative() == false)
			{
				acceptNewSpeciesValues();
				++numberOfIterations;
				t += dt;
				isNegative = false;
			}
			else
			{
				cout << "Negative species at time: " << t << endl;
				dt = dt * 0.5;
				reloadProposedSpeciesValues();
				isNegative = true;
			}			
		}
		
		saveData();
		cout << "Sample: " << samples << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
	}
	writeData(outputFileName);
	closeAuxiliaryStream();
}

