/*
 *  TauLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "TauLeaping.h"

TauLeaping::TauLeaping(Simulation * simulation):
LeapMethod(simulation)
{
}

TauLeaping::~TauLeaping()
{
}

double TauLeaping::computeTimeStep()
{
	double epsilon	= simulation->Epsilon;
	
	int numberOfSpecies		= sbmlModel->getNumSpecies();
	Array<int, 1> hor			(numberOfSpecies);
	Array<int, 1> nuHor			(numberOfSpecies);
	Array<double, 1> muHat		(numberOfSpecies);
	Array<double, 1> sigmaHat2	(numberOfSpecies);
	Array<double, 1> varHat		(numberOfSpecies);
	hor = 0; nuHor = 0; muHat = 0.0; sigmaHat2 = 0.0;
	
	computeHor(hor, nuHor);
	computeMuHatSigmaHat2(muHat, sigmaHat2);
	
	double tau, taup,  epsi, epsixi, epsixisq;
	double xi;
	
	tau = HUGE_VAL;
	
	double a0 = (double)blitz::sum(propensitiesVector);
	for (int is = 0; is < numberOfSpecies; is++)
	{
		varHat(is) = sigmaHat2(is) - (1.0/a0) * muHat(is) * muHat(is);
	}
	
	for (int is = 0; is < numberOfSpecies; ++is)
	{
		taup = (HUGE_VALF*0.5);
		xi = (double)simulation->speciesValues(is);
		switch (hor(is)) {
			case 0:
				break;
			case 1:
				epsi = epsilon;
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabsf(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				break;
			case 2:
				if (nuHor(is) == 1)
					epsi = 0.5*epsilon;
				else
					epsi = epsilon*(xi-1.0)/(2.0*(xi-1.0)+1.0);
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabs(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				break;
			case 3:
				if (nuHor(is)==1)
					epsi = 0.3333333333*epsilon;
				else if (nuHor(is) == 2)
					epsi = epsilon*(xi-1)/(3.0*(xi-1)+1.5);
				else
					epsi = epsilon*(xi-1)*(xi-2)/(3.0*(xi-1)*(xi-2)+(xi-2)+2.0*(xi-1));
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabsf(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				break;
			default:
				break;
		}
	}
	
	return tau;	
}

// if proposed time step is too small,execute numberOfIterations steps of SSA   
void TauLeaping::_executeSSA(double& t, int SSAsteps)
{
	int count = 0.;
	double a0 = 0.;
	double tau;
	double r1;
	int reactionIndex = 0;
	double cummulative = 0.0;
	
	while (count < SSAsteps)
	{
		count++;
		computePropensities(propensitiesVector, 0);
		a0 = blitz::sum(propensitiesVector);
		tau = (1.0/a0) * sgamma( (double)1.0 );
		
		r1 = ranf();
		reactionIndex = -1;
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
		
		if (reactionIndex != -1)
		{
			fireReaction(reactionIndex, 1);
			t += tau;
		}
		else 
		{ 
			t = HUGE_VAL;
			break;
		}
	}
}

void TauLeaping::_writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum)
{
	double aver_dt = dt_sum/ (double) steps;
	double aver_L = L_sum/ (double) steps;
	
	// time, dt, aver_dt, L,averL, #iterations
	if (myfile!=NULL)
	{
		fprintf(myfile, "%f  %e  %e  %i  %f  %i \n", t, dt, aver_dt, L, aver_L, steps  );
	}
}

void TauLeaping::solve()
{
#ifndef NDEBUG 
	cout << "TauLeaping..." << endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	FILE* rejectionsfile = fopen("Rejections_T.txt", "w");
#endif
	
	double aj;
	long int kj;
	double a0						= 0.0;
	bool isNegative					= false;
	double averNumberOfRealizations = 0.0;
	double numberOfRejections		= 0.0;
	
	#ifndef NSSASTEP
	const int SSAfactor = 10;
	const int SSAsteps = 1;
	#endif

	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		numberOfRejections = 0.;
		timePoint = 0;
		zeroData();
		simulation->loadInitialConditions();
		isNegative = false;

		
#ifndef NDIAGNOSTIC  //diagnostic
		FILE* myfile = fopen("Diag_T.txt", "w");
		double whenToWriteOffset = tEnd / numberOfFrames;
		double whenToWrite = whenToWriteOffset;
		
		int iterSteps = 0;
		double dt_sum = 0.;
		long int L_sum = 0;
		long int L = 0;
		long int SSA_L = 0;
		double LAverage = 0.0;
#endif
		
		while (t < tEnd)
		{
			
#ifndef NDEBUG
			saveData();
#endif
			
			computePropensities(propensitiesVector, 0);
			a0 = blitz::sum(propensitiesVector);
			
			if (isNegative == false)
			{
				dt = computeTimeStep();
				//if (dt >= HUGE_VAL) {dt= tEnd*10; cout<<"stop"<<endl;	break;}
			}
			
#ifndef NSSASTEP
			if (dt <= SSAfactor * (1.0/a0) * sgamma( (double)1.0 ) ) 
			{
				_executeSSA(t, SSAsteps);
				#ifndef NDEBUG
				SSA_L += SSAsteps;
				#endif
			}
			else
#endif
			{				
				// sampling
				for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
				{				
					aj				= propensitiesVector(j);
					
					if (aj != 0) 
					{
						kj = ignpoi( aj*dt );
										
						fireReactionProposed( j , kj );
						#ifndef NDIAGNOSTIC
						L += kj;
						#endif
					}
				}
				//cout << endl;

				
				if (isProposedNegative() == false)
				{
					acceptNewSpeciesValues();
					++numberOfIterations;
					t += dt;
					isNegative = false;

#ifndef NDIAGNOSTIC  // diagnostic

					LAverage += (double) L;
					L=0;
					
					iterSteps = iterSteps + 1;
					L_sum  += L;
					dt_sum +=dt;
					
					if ( t >= ((double) (whenToWrite)))
					{
						_writeDiagnostic(myfile, L, iterSteps, L_sum, dt_sum);
						whenToWrite = whenToWrite + whenToWriteOffset;
						
						// reset counters
						iterSteps = 0;
						L_sum = 0.;
						dt_sum = 0.;
					}
#endif
				}
				else
				{
#ifndef NDEBUG			
					cout << "Negative species at time: " << t << endl;
					++numberOfRejections;
#endif
					dt = dt * 0.5;
					reloadProposedSpeciesValues();
					isNegative = true;
				}	
				
			} // end of if (tau <tauSSA )
		} // end of while (t< tEnd)
				
#ifndef NDEBUG

		saveData();
		
		cout << "Sample: " << samples << endl;
		#ifndef NSSASTEP
		cout<< "Number of SSA steps: "<< SSA_L / SSAsteps << endl;
		#endif
		
		#ifndef NDIAGNOSTIC
		cout << "Average L: " << LAverage/((double)numberOfIterations) << endl;
		fclose(myfile);
		#endif
		
		writeToAuxiliaryStream( simulation->speciesValues );
		averNumberOfRealizations += numberOfIterations;
		
		// report on negative population
		if (rejectionsfile!=NULL)
			fprintf(rejectionsfile, "%i %f %i \n", samples, numberOfRejections ,numberOfIterations);
#endif

	}
	
#ifndef NDEBUG
	writeData(outputFileName);
	closeAuxiliaryStream();
	cout << " Average number of Realizations in Tau-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
	
	fclose(rejectionsfile);

#endif
	
}

