/*
 *  TauLeaping_LacZLacY.cpp
 *  SSM_Xcode
 *
 *  Created by Lipkova on 9/4/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TauLeaping_LacZLacY.h"

// Constructor & destructor
TauLeaping_LacZLacY::TauLeaping_LacZLacY(Simulation* simulation):
LeapMethod(simulation)
{
}

TauLeaping_LacZLacY::~TauLeaping_LacZLacY()
{
}


void TauLeaping_LacZLacY::identifyReactions(vector<int>& criticalReactions, vector<int>& nonCriticalReactions)
{
	int Ncrit = 10;
	
	for (int ir = 0; ir < sbmlModel->getNumReactions(); ++ir)
	{
		long int lj = 2147483647; // MAXIMUM INTEGER
		SSMReaction * ssmReaction = simulation->ssmReactionList[ir];
		if (propensitiesVector(ir) > 0.0)
		{
			const vector<int> & changes		= ssmReaction->getChanges();
			const vector<int> & nuChanges	= ssmReaction->getNuChanges();
			
			for (int is = 0; is < changes.size() ; ++is)
			{
				if (nuChanges[is] > 0) break;
				lj = min(lj, -simulation->speciesValues(changes[is])/ nuChanges[is] );
			}
		}
		
		if (lj < Ncrit) 
			criticalReactions.push_back(ir);
		else
			nonCriticalReactions.push_back(ir);
	}
}

//****************************
//  compute time step
//****************************
double TauLeaping_LacZLacY::computeTimeStep(vector<int> criticalReactions, vector<int> nonCriticalReactions)
{
	
	double epsilon	        = simulation->Epsilon;	
	int numberOfSpecies		= sbmlModel->getNumSpecies();
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	Array<int, 1> hor			(numberOfSpecies);
	Array<int, 1> nuHor			(numberOfSpecies);
	Array<double, 1> muHat		(numberOfSpecies);
	Array<double, 1> sigmaHat2	(numberOfSpecies);
	Array<double, 1> varHat		(numberOfSpecies);
	
	hor = 0; nuHor = 0; muHat = 0.0; sigmaHat2 = 0.0;
	computeHor(hor, nuHor);
	
	//computeMuHatSigmaHat2
	if (criticalReactions.empty())
		LeapMethod::computeMuHatSigmaHat2(muHat, sigmaHat2);
	else  // exclude critical reactions 
		computeMuHatSigmaHat2(muHat, sigmaHat2,nonCriticalReactions);	
	
	
	double tau, taup,  epsi, epsixi, epsixisq;
	double xi;
	
    tau = HUGE_VAL;
	
	double a0 = (double)blitz::sum(propensitiesVector);
	for (int is = 0; is < numberOfSpecies; is++)
		varHat(is) = sigmaHat2(is) - (1.0/a0) * muHat(is) * muHat(is);
	
	
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
	
};
//****************************

//****************************
// time of critical reaction
//****************************
double TauLeaping_LacZLacY::computeTimeOfCritical(vector<int> criticalReactions, double& ac0) 
{
	double tau2;
	ac0 = 0.0;
//	
//	for (vector<int>::iterator it = criticalReactions.begin(); it!=criticalReactions.end() ; ++it)
//		ac0 += propensitiesVector(*it);
	
	if (criticalReactions.empty()) 
		tau2 = HUGE_VAL;
	else
	{
		tau2 = (1.0/ac0) * sgamma( (double)1.0 );
		
		for (vector<int>::iterator it = criticalReactions.begin(); it!=criticalReactions.end() ; ++it)
			ac0 += propensitiesVector(*it);
	}
	
	return tau2;
};
//****************************


//****************************
//  muHat, sigmaHat2 computation:
//****************************
// this method is overloaded from the LeapMethod class, since in adaptive tau leaping we exclude
// critical reactions from computation of explicit tau and additionaly we exclude reactions in PEC
// in computation of implicit tau
void TauLeaping_LacZLacY::computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2, vector<int> non_critical)
{
	int is, ir, ns, indx, nr;
	double tmpfloat;
	nr = sbmlModel->getNumReactions();
	
	for (int numbS = 0; numbS < sbmlModel->getNumSpecies(); ++numbS)
	{
		muHat(numbS) = 0.0;
		sigmaHat2(numbS) = 0.0;
	}
	
	// consider only non-critical reactions
	for (vector<int>::iterator ir = non_critical.begin(); ir!=non_critical.end(); ++ir)	
	{
		SSMReaction* ri = simulation->ssmReactionList[*ir];
		double  riPropensity = propensitiesVector(*ir);
		
		const vector<int> & changes = ri->getChanges();
		const vector<int> & nuChanges = ri->getNuChanges();
		
		ns = changes.size();
		for (is = 0; is < ns; is++ )
		{
			indx = changes[is];
			tmpfloat = nuChanges[is] * riPropensity;
			muHat(indx) += tmpfloat;
			sigmaHat2(indx) += nuChanges[is] * tmpfloat;
		}
	}
}
//****************************


//****************************
//     Execute SSA
//****************************
// if proposed time step is too small,execute numberOfIterations steps of SSA   
void TauLeaping_LacZLacY::_executeSSA(double& t, int SSAsteps, double genTime)
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
		computePropensitiesGrowingVolume(propensitiesVector, t, genTime);
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
			if (t > tEnd) 
				break;
		}
		else 
		{ 
			t = HUGE_VAL;
			break;
		}
	}
	
}
//****************************

//****************************
//   For given tau sample reactions
//****************************

//
void TauLeaping_LacZLacY::sampling(short int crit, vector<int> criticalReactions, vector<int> nonCriticalReactions, long int& L, double ac0)
{
	
#ifndef NDEBUG
	L = 0;
#endif
	
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	double aj;
	long int kj;
	
	if (criticalReactions.empty()) // no critical reactions, sample from all chanels
	{
		for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
		{				
			aj = propensitiesVector(j);
			
			if (aj != 0) 
			{
				kj = ignpoi( aj*dt );
				fireReactionProposed( j , kj );
				
#ifndef NDEBUG
				L += kj;
#endif
			}
		}
	}
	else  
	{
		// sample from non-critical
		int k = 0;
		for (vector<int>::iterator it=nonCriticalReactions.begin(); it!=nonCriticalReactions.end(); ++it) {
			aj = propensitiesVector(*it);
			if (aj != 0) {
				kj =  ignpoi(aj*dt);
				fireReactionProposed( *it , kj );
				#ifndef NDEBUG
				L += kj;
				#endif
			}
		}
	}
	
	// from all critical reactions sample one
	if (crit == 1)
	{
		double r1 = ranf();
		int reactionIndex = 0;
		vector<int>::iterator it = criticalReactions.begin();
		double cummulative = propensitiesVector(*it);
		
		while (cummulative < r1 * ac0)
		{
			++it;
			cummulative += propensitiesVector(*it);
		}
		
		fireReactionProposed( *it , 1 );
		
#ifndef NDEBUG		
		L++;
#endif
	}
}


//****************************
//   Diagnostics and Trajecotries 
//****************************
void TauLeaping_LacZLacY::_writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum)
{
	double aver_dt = dt_sum/ (double) steps;
	double aver_L = L_sum/ (double) steps;
	
	// time, dt, aver_dt, L,averL, #iterations
	if (myfile!=NULL)
	{
		fprintf(myfile, "%f  %e  %e  %i  %f  %f \n", t, dt, aver_dt, L, aver_L, steps  );
	}
}


void TauLeaping_LacZLacY:: _writeTrajectories(FILE* myfile, vector<int> data)
{
	if (myfile!=NULL) 
	{
		for (int i=0; i < data.size() ; ++i) 
			fprintf(myfile, "%i ", data[i], " " );
		
		fprintf(myfile, "\n");
	}
}
//****************************

//****************************
//  Solve it :-P
//****************************
void TauLeaping_LacZLacY::solve()
{
	
#ifndef NDEBUG 
	cout << "Tau Leaping with List of Critical reactions..." <<endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	FILE* rejectionsfile = fopen("Rejections_TNN.txt", "w");
#endif
	
#ifdef TRAJECTORIES
	FILE* file_RbsLacY   = fopen("T_RbsLacY.txt"  , "w");
	FILE* file_TrLacZ2   = fopen("T_TrLacZ2.txt"  , "w");
	FILE* file_TrRbsLacY = fopen("T_TrRbsLacY.txt", "w");
	FILE* file_TrRbsLacZ = fopen("T_TrRbsLacZ.txt", "w");
#endif
	
	/* initialize helping variables */
	double a0						= 0.0;
	double ac0						= 0.0;
	double dt2						= 0.; 
	long int L						= 0;	
	short int crit					= 0;
	bool isNegative					= false;
	double averNumberOfRealizations = 0.0;
	double numberOfRejections;
	
	double genTime = 2100.;   // generation time	
	
#ifndef NSSASTEP
	const int SSAfactor = 10;
	const int SSAsteps = 100;
#endif
	
	/* execute code for each realization */
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		whenToSave			= simulation->StartTime;
		numberOfIterations	= 0;
		numberOfRejections	= 0.;
		timePoint			= 0;
		zeroData();
		simulation->loadInitialConditions();
		isNegative			= false;
		
#ifdef TRAJECTORIES
		vector<int> RbsLacY;
		vector<int> TrLacZ2;
		vector<int> TrRbsLacY;
		vector<int> TrRbsLacZ;
		
		double whenToWriteTrajectories		 = simulation->StartTime;
		double whenToWriteTrajectoriesOffset = 1.;
#endif
		
#ifndef NDIAGNOSTIC  //diagnostic
		FILE* myfile = fopen("Diag_TNN.txt", "w");
		double whenToWriteOffset = tEnd / numberOfFrames;
		double whenToWrite = whenToWriteOffset;
		
		int steps = 0;
		double dt_sum = 0.;
		long int L_sum = 0;
		double LAverage = 0.0;
#endif
		
#ifndef NDEBUG
		long int SSA_L = 0;
#endif

		
		while (t < tEnd)
		{
			#ifndef NDEBUG
			saveData();
			#endif
			
			// set random variables values
			// RNAP     = S(1) ~ N(35),3.5^2)
			// Ribosome = S(9) ~ N(350,35^2)
			// + mean values of these pools growth with cell volume so that the concentrations of these molecules remain constant
			simulation->speciesValues(1)  = gennor(35  * (1 + t/genTime),  3.5);
			simulation->speciesValues(9)  = gennor(350 * (1 + t/genTime),   35);
			
			computePropensitiesGrowingVolume(propensitiesVector, t, genTime);
			a0 = blitz::sum(propensitiesVector);
			
			/* 1) get critical reactions */
			vector<int> criticalReactions;
			vector<int> nonCriticalReactions;
			identifyReactions(criticalReactions,nonCriticalReactions);
			crit = 0;						// number of critical reactions to be executed;
			
			/* 2) compute time step */
			if (isNegative == false)
			{
				dt = computeTimeStep(criticalReactions, nonCriticalReactions);
				if ( (t+dt) > tEnd) 
					dt = tEnd - t;
			}
			
#ifndef NSSASTEP
			if (dt <= SSAfactor * (1.0/a0) ) // * sgamma( (double)1.0 ) ) 
			{
				_executeSSA(t, SSAsteps, genTime);
				#ifndef NDEBUG 
				SSA_L += SSAsteps;
				#endif
				
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
			}
			else
#endif  // end of SSASTEP
			{
				dt2 = computeTimeOfCritical(criticalReactions, ac0);				
				if (dt2 < dt) {
					dt	 =	dt2;
					crit = 1;
				}
				
				sampling(crit, criticalReactions,nonCriticalReactions, L,ac0);
				
				/* 5) check proposed update and update the system or repeat simulaiton */
				if (isProposedNegative() == false)
				{
					acceptNewSpeciesValues();
					++numberOfIterations;
					t += dt;
					isNegative = false;
					
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
					
#ifndef NDIAGNOSTIC  // diagnostic
					LAverage += (double) L;
					steps++;
					L_sum  += L;
					dt_sum +=dt;
					
					//if ( t >= ((double) (whenToWrite)))
					{
						_writeDiagnostic(myfile, L, steps, L_sum, dt_sum);
						whenToWrite = whenToWrite + whenToWriteOffset;
						
						// reset counters
						steps = 0;
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
				
			}   // end of if (SSA step) 
		}  // end of while( t < tEnd)
		
#ifdef TRAJECTORIES
		_writeTrajectories(file_RbsLacY,RbsLacY);
		_writeTrajectories(file_TrLacZ2,TrLacZ2);
		_writeTrajectories(file_TrRbsLacY,TrRbsLacY);
		_writeTrajectories(file_TrRbsLacZ,TrRbsLacZ);
		
		cout << "Sample: " << samples << endl;
#endif
		
		
#ifndef NDEBUG
		saveData();		
		cout << "Sample: " << samples << endl;
	#ifndef NDIAGNOSTIC
		cout << "Average L: " << LAverage/((double)numberOfIterations) << endl;
		fclose(myfile);
		#ifndef NSSASTEP
		cout<< "Number of SSA steps: "<< SSA_L / SSAsteps << endl;
		#endif
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
	cout << " Average number of Realizations in Non Negative Tau-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
	
	fclose(rejectionsfile);
#endif
	
#ifdef TRAJECTORIES
	fclose(file_RbsLacY);
	fclose(file_TrLacZ2);
	fclose(file_TrRbsLacY);
	fclose(file_TrRbsLacZ);
#endif	
}