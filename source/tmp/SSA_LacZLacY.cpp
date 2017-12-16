/*
 *  SSA_LacZLacY.cpp
 *  SSM_Xcode
 *
 *  Created by Lipkova on 9/3/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "SSA_LacZLacY.h"



SSA_LacZLacY::SSA_LacZLacY(Simulation * simulation):
Method(simulation)
{ }

SSA_LacZLacY::~SSA_LacZLacY()
{ }


void SSA_LacZLacY::_writeDiagnostic(FILE* myfile, int steps, double dt_sum)
{
	double aver_dt = dt_sum/ (double) steps;
	
	// time, dt, aver_dt, L,averL, #iterations
	if (myfile!=NULL)
		fprintf(myfile, "%f  %e  %f \n", t, dt, steps  );
}

void SSA_LacZLacY:: _writeTrajectories(FILE* myfile, vector<int> data)
{
	if (myfile!=NULL) 
	{
		for (int i=0; i < data.size() ; ++i) 
			fprintf(myfile, "%i ", data[i], " " );
		
		fprintf(myfile, "\n");
	}
}

void SSA_LacZLacY::solve()
{
	cout << "SSA LacZLacY..." << endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	
#ifdef TRAJECTORIES
	FILE* file_RbsLacY   = fopen("SSA_RbsLacY.txt"  , "w");
	FILE* file_TrLacZ2   = fopen("SSA_TrLacZ2.txt"  , "w");
	FILE* file_TrRbsLacY = fopen("SSA_TrRbsLacY.txt", "w");
	FILE* file_TrRbsLacZ = fopen("SSA_TrRbsLacZ.txt", "w");
#endif
	
	double a0 = 0.0;
	double r1;
	int reactionIndex = 0;
	double cummulative = 0.0;
	double averNumberOfRealizations = 0.0;
	
	double genTime = 2100.;   // generation time
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t					= simulation->StartTime;
		whenToSave			= simulation->StartTime;
		numberOfIterations  = 0;
		timePoint			= 0;
		zeroData();
		simulation->loadInitialConditions();
		
#ifdef TRAJECTORIES
		vector<int> RbsLacY;
		vector<int> TrLacZ2;
		vector<int> TrRbsLacY;
		vector<int> TrRbsLacZ;
		
		double whenToWriteTrajectories		 = simulation->StartTime;
		double whenToWriteTrajectoriesOffset = 1.;
#endif
		
		
		while (t < tEnd)
		{
			saveData();
			
			// set random variables values
			// RNAP     = S(1) ~ N(35),3.5^2)
			// Ribosome = S(9) ~ N(350,35^2)
			// + mean values of these pools growth with cell volume so that the concentrations of these molecules remain constant
			simulation->speciesValues(1)  = gennor(35  * (1 + t/genTime),  3.5);
			simulation->speciesValues(9)  = gennor(350 * (1 + t/genTime),   35);

			computePropensitiesGrowingVolume(propensitiesVector, t, genTime);
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
			
			
#ifdef TRAJECTORIES
			// write trajectories of RbsLacY=S(7), TrLacZ2=S(5),TrRbsLacY=S(13),TrRbsLacZ=S(12)
			if (t >= whenToWriteTrajectories )
			{
				RbsLacY.push_back(simulation->speciesValues(7));
				TrLacZ2.push_back(simulation->speciesValues(5));
				TrRbsLacY.push_back(simulation->speciesValues(13));
				TrRbsLacZ.push_back(simulation->speciesValues(12));
				
				whenToWriteTrajectories += whenToWriteTrajectoriesOffset;
			}
#endif
		}  // end of while loop
		
#ifdef TRAJECTORIES
		_writeTrajectories(file_RbsLacY,RbsLacY);
		_writeTrajectories(file_TrLacZ2,TrLacZ2);
		_writeTrajectories(file_TrRbsLacY,TrRbsLacY);
		_writeTrajectories(file_TrRbsLacZ,TrRbsLacZ);	
		cout << "Sample: " << samples << endl;
#endif
		
		saveData();
		cout << "Sample: " << samples << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
		//writeData(localOutputFileName,samples);
		averNumberOfRealizations += numberOfIterations;
	}
	
	writeData(outputFileName);
	closeAuxiliaryStream();
	cout << " Average number of Realizations in Gillespie SSA:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
	
#ifdef TRAJECTORIES
	fclose(file_RbsLacY);
	fclose(file_TrLacZ2);
	fclose(file_TrRbsLacY);
	fclose(file_TrRbsLacZ);
#endif
}
