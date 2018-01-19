/*
 *  SSA.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 *
 */

#include "SSA.h"

SSA::SSA(Simulation * simulation):
Method(simulation)
{ }

SSA::~SSA()
{ }

void SSA::solve()
{
	cout << "SSA..." << endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");

	double a0 = 0.0;
	double r1;
	int reactionIndex = 0;
	double cummulative = 0.0;
	double averNumberOfRealizations = 0.0;
	double genTime = 2100;

	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		timePoint = 0;
		//XXX
		whenToSave = t;
		zeroData();
		simulation->loadInitialConditions();



		#ifdef DEBUG_PRINT
			tempArray.resize(sbmlModel->getNumSpecies());
			myfile.open ("all-times.txt");

			myfile << t << "\t";
			tempArray = simulation->speciesValues(Range::all());
			for (int i = 0; i < tempArray.extent(firstDim); ++i)
			{
				myfile << tempArray(i) << "\t";
			}
			myfile << endl;
		#endif

		saveData();

		while (t < tEnd)
		{

#ifdef LacZLacY
    		// RNAP     = S(1) ~ N(35),3.5^2)
    		// Ribosome = S(9) ~ N(350,35^2)
    		simulation->speciesValues(1)  = 35;// gennor(35   * (1 + t/genTime), 3.5);
		simulation->speciesValues(9)  = 350;//gennor(350  * (1 + t/genTime),  35);
	     	//computePropensitiesGrowingVolume(propensitiesVector,t,genTime);
		computePropensities(propensitiesVector, 0);
#else
			computePropensities(propensitiesVector, 0);
#endif


			a0 = blitz::sum(propensitiesVector);
			dt = (1.0/a0) * sgamma( (double)1.0 );

			r1 = ranf();
			reactionIndex = 0;
			cummulative = 0.0;
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j){
				cummulative += propensitiesVector(j);
				if ( cummulative > a0*r1 ){
					reactionIndex = j;
					break;
				}
			}

			fireReaction(reactionIndex, 1);
			++numberOfIterations;

			//XXX
			t_old = t;
			t += dt;

			saveData();

			#ifdef DEBUG_PRINT
				myfile << min(t,tEnd) << "\t";
				if(t<tEnd)
					tempArray =  simulation->speciesValues(Range::all());
				else
					tempArray =  simulation->old_speciesValues(Range::all());

				for (int i = 0; i < tempArray.extent(firstDim); ++i){
					myfile << tempArray(i) << "\t";
				}
				myfile << endl;
		 	#endif

		}
        
        cout << "Sample: " << samples << endl;
        writeToAuxiliaryStream( simulation->speciesValues );
        averNumberOfRealizations += numberOfIterations;

		#ifdef DEBUG_PRINT
				myfile.close();
		#endif

	}

	writeData(outputFileName);
	closeAuxiliaryStream();

	cout << " Average number of Realizations in Gillespie SSA:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;

}































/*
{
	cout << "SSA..." << endl;

	double a0				= 0.0;
	double r1;
	int reactionIndex		= 0;
	double cummulative		= 0.0;

	Array< double , 1 > cummulativeSum;
	cummulativeSum.resize( propensitiesVector.extent(firstDim) );
	cummulativeSum			= 0.0;
	double	meanCounter		= 0.0;

	openAuxiliaryStream( (simulation->ModelName) + "-histogram-ssa.txt");
	vector<double> histogram;

	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations	= 0;
		timePoint			= 0;
		simulation->loadInitialConditions();

		computePropensities		(propensitiesVector, 0); // O (M)
		while (t < tEnd)
		{
			saveData();

			//computePropensities		(propensitiesVector,	reactionIndex);
			//computeCummulativeSum	(cummulativeSum,		reactionIndex);

			//computePropensities		(propensitiesVector,	reactionIndex, reactionIndex+2); // O (1)
			//computeCummulativeSum	(cummulativeSum,		reactionIndex);					 // O (1/2 M)

			int minReactionIndex;
			computePropensitiesQuickly	( propensitiesVector, reactionIndex, minReactionIndex	);  // O (M ?)
			computeCummulativeSum		(cummulativeSum,      minReactionIndex					);	// O (1/2 M)

			a0				= cummulativeSum( cummulativeSum.extent(firstDim)-1 );

			// refrain from dividing by zero
			if (a0 == 0.0)
			{
				t += dt;
				continue;
			}
			dt				= (1.0/a0) * sgamma( (double)1.0 );

			r1				= ranf();
			reactionIndex	= 0;
			cummulative		= 0.0;
			int counter		= -1;
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
			{
				++counter;
				cummulative += propensitiesVector(j);
				if ( cummulative > a0*r1 )
				{
					reactionIndex = j;
					break;
				}
			}

			meanCounter		+= ((double)counter);
			fireReaction(reactionIndex, 1);
			++numberOfIterations;
			t += dt;

		}
		saveData();
		cout << "Sample: " << samples << endl;
		//cout << "Mean SSA counter: " << meanCounter/((double)numberOfIterations) << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
	}
	writeData(outputFileName);
	closeAuxiliaryStream();
}
*/
