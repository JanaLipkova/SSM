/*
 *  SSA.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "SSA.h"

SSA::SSA(Simulation * simulation):
Method(simulation)
{
}

SSA::~SSA()
{
}



void SSA::_writeDiagnostic(FILE* myfile, int steps, double dt_sum)
{
	double aver_dt = dt_sum/ (double) steps;

	// time, dt, aver_dt, L,averL, #iterations
	if (myfile!=NULL)
	{
		fprintf(myfile, "%f  %e  %f \n", t, dt, steps  );
	}
}

void SSA::solve()
{
	cout << "SSA..." << endl;
#ifndef NDEBUG
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
#endif

	double a0 = 0.0;
	double r1;
	int reactionIndex = 0;
	double cummulative = 0.0;
	double averNumberOfRealizations = 0.0;



	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		timePoint = 0;
		zeroData();
		simulation->loadInitialConditions();


#ifndef NDIAGNOSTIC  //diagnostic
		FILE* myfile = fopen("SSA_Diag.txt", "w");
		double whenToWriteOffset = tEnd / numberOfFrames;
		double whenToWrite = whenToWriteOffset;

		int steps = 0;
		double dt_sum = 0.;
#endif


		while (t < tEnd)
		{
#ifndef NDEBUG
			saveData();
#endif
			computePropensities(propensitiesVector, 0);
			a0 = blitz::sum(propensitiesVector);

			dt = (1.0/a0) * sgamma( (double)1.0 );

			if (t+dt >= tEnd)
				break;

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

			fireReaction(reactionIndex, 1);
			++numberOfIterations;
			t += dt;

#ifndef NDIAGNOSTIC  // diagnostic

			steps++;
			dt_sum +=dt;

			if ( t >= ((double) (whenToWrite)))
			{
				_writeDiagnostic(myfile, steps,dt_sum);
				whenToWrite = whenToWrite + whenToWriteOffset;

				// reset counters
				steps = 0;
				dt_sum = 0.;
			}
#endif
		}

#ifndef NDEBUG
		saveData();

		cout << "Sample: " << samples << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
		//writeData(localOutputFileName,samples);
		averNumberOfRealizations += numberOfIterations;

		#ifndef NDIAGNOSTIC
		fclose(myfile);
		#endif
#endif

	}

#ifndef NDEBUG
	writeData(outputFileName);
	closeAuxiliaryStream();
	cout << " Average number of Realizations in Gillespie SSA:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
#endif
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
