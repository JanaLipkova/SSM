/*
 *  SSALDM.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "SSALDM.h"

SSALDM::SSALDM(Simulation * simulation):
Method(simulation)
{
}

SSALDM::~SSALDM()
{
}

int SSALDM::center( int left, int right )
{
	double v = 0.5*((double)right + (double)left);
	return ((int)(v));
}

int SSALDM::binaryTreeSearch( Array< double , 1 > & cummulativeSum, double key, int left, int right, int & counter )
{	
	counter++;
			
	// check for a solution
	if ( key < cummulativeSum(right) && key >= cummulativeSum(right-1) )
	{ return right;}
	
	if ( left > 0)
	{
		if ( key < cummulativeSum(left) && key >= cummulativeSum(left-1) )
		{ return left;}
	}
	else
	{
		if ( key < cummulativeSum(left) )
		{ return left;}
	}

	int arrayCenter = center(left, right);
	if ( key > cummulativeSum(arrayCenter) )
	{ return binaryTreeSearch( cummulativeSum, key, arrayCenter, right, counter ); }
	else
	{ return binaryTreeSearch( cummulativeSum, key, left, arrayCenter, counter ); }
}

void SSALDM::solve()
{
	cout << "SSA-LDM..." << endl;
		
	Array< double , 1 > cummulativeSum;
	cummulativeSum.resize( propensitiesVector.extent(firstDim) );
	cummulativeSum = 0.0;	
	
	int		reactionIndex = 0;	
	int		counter       = -1;
	double	meanCounter   = 0.0;
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		timePoint = 0;
		simulation->loadInitialConditions();
		
		computePropensities		(propensitiesVector, 0); // O (M)
		
		while (t < tEnd)
		{
			saveData();
			
			//computePropensities		(propensitiesVector,	reactionIndex, reactionIndex+2);	// O (1)
			//computeCummulativeSum	(cummulativeSum,		reactionIndex);							// O (1/2 M)
			
			int minReactionIndex;
			computePropensitiesQuickly	( propensitiesVector, reactionIndex, minReactionIndex	);  // O (M ?)
			computeCummulativeSum		(cummulativeSum,      minReactionIndex					);	// O (1/2 M)
			
			double a0		= cummulativeSum( cummulativeSum.extent(firstDim)-1 );					// O (1)
			
			// refrain from dividing by zero
			if (a0 == 0.0)
			{ 
				t += dt;
				continue; 
			}
			
			dt				= (1.0/a0) * sgamma( (double)1.0 );
			double	r1		= ranf();
	
			// logarithmic direct method
			counter			= -1;
			double key		= a0*r1;
			
			reactionIndex	= binaryTreeSearch( cummulativeSum, key, 0, cummulativeSum.extent(firstDim)-1, counter ); // O (log M)
			meanCounter		+= ((double)counter);
						
			fireReaction(reactionIndex, 1);
			++numberOfIterations;
			t += dt;
		}
		cout << "Sample: " << samples << endl;
		cout << "Mean binary tree counter: " << meanCounter/((double)numberOfIterations) << endl;
		saveData();
	}
	writeData(outputFileName);
}

