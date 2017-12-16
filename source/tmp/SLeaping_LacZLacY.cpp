/*
 *  SLeaping_LacZLacY.cpp
 *  SSM_Xcode
 *
 *  Created by Lipkova on 9/8/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "SLeaping_LacZLacY.h"


SLeaping_LacZLacY::SLeaping_LacZLacY(Simulation * simulation):
LeapMethod(simulation)
{
}

SLeaping_LacZLacY::~SLeaping_LacZLacY()
{
}


void SLeaping_LacZLacY::computePropensitiesGrowingVolume(double genTime)
{
	
	double volume	= 1. + t/genTime;
	double ivolume	= 1./volume;
	
	int nu;
	ParticleType x;
	ParticleType num, denom;
	
	int ir;		// the reaction index
	
	propensitiesVector = 0.0;
	
	vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;
	
	Reaction * sbmlreaction; // added maagm
	
	for (int ev = 0; ev < eventVector.size(); ++ev) // event index
	{
		eventVector[ev]->propensity = 0.0;
		ir							= eventVector[ev]->index;
		
		//added maagm
		SSMReaction* reaction		= ssmReactionList[ir];
		
		sbmlreaction = sbmlModel->getReaction(ir);
		KineticLaw * kineticLaw = sbmlreaction->getKineticLaw();
		Parameter * parameter = kineticLaw->getParameter(0);
		double rate = parameter->getValue();
		
		if (kineticLaw->getNumParameters() == 5)
		{
			int dependentSpecies = getDependentSpecies(ir);
			double dependentValue = (double)simulation->speciesValues(dependentSpecies);
			double h =					kineticLaw->getParameter(2)->getValue();
			double defaultProduction =	kineticLaw->getParameter(3)->getValue();
			double cHill =				kineticLaw->getParameter(4)->getValue();
			//propensitiesVector(ir) = defaultProduction + rate*hillFunction(cHill, dependentValue, h);
			reaction->setPropensity(defaultProduction + rate*hillFunction(cHill, dependentValue, h));
		}
		else
		{
			//original bbayati
			SSMReaction* reaction		= ssmReactionList[ir];
			vector <int>  reactants		= reaction->getReactants();
			vector <int>  nu_reactants	= reaction->getNuReactants();
			int order					= reaction->getOrder();
			
			reaction->setPropensity(reaction->getRate());
			
			for (int s = 0; s < reactants.size(); ++s)
			{
				nu		= nu_reactants[s];
				x		= simulation->speciesValues( reactants[s] );
				num		= x;
				denom	= nu;
				while ((--nu)>0)
				{
					denom	*= nu;
					num		*= (x - nu);
				}
				reaction->setPropensity( reaction->getPropensity()*((double)num/(double)denom) );
			}
			
			if (order == 2)
				reaction->setPropensity( reaction->getPropensity() * ivolume );
			
			if (order == 3) 
				reaction->setPropensity( reaction->getPropensity() * ivolume *ivolume);
			
			if (order > 3) 
			{
				std::cout<<"Aborting: Growing volume of reaction enviroment do not support reaction of order higher than 3, if you want it implement it"<<std::endl;
				std::abort();
			}	
		}
		
		propensitiesVector(ir)		= reaction->getPropensity();
		eventVector[ev]->propensity	= reaction->getPropensity();
	}
}



void SLeaping_LacZLacY::identifyReactions(vector<int>& criticalReactions, vector<int>& nonCriticalReactions)
{
	int Ncrit = 10;
	
	for (int ev=0; ev < eventVector.size(); ++ev)
	{
		long int lj = 2147483647; // MAXIMUM INTEGER
		int      ir = eventVector[ev]->index;
		SSMReaction * ssmReaction = simulation->ssmReactionList[ir];
		if (eventVector[ev]->propensity > 0.0)
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
			criticalReactions.push_back(ev);
		else
			nonCriticalReactions.push_back(ev);
	}
}

double SLeaping_LacZLacY::computeTimeStep(vector<int> criticalReactions, vector<int> nonCriticalReactions)
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

void SLeaping_LacZLacY::computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2, vector<int> nonCriticalReactions)
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
	for (vector<int>:: iterator ev=nonCriticalReactions.begin(); ev!=nonCriticalReactions.end(); ++ev)
	{
		SSMReaction* ri = simulation->ssmReactionList[eventVector[*ev]->index];
		double  riPropensity = eventVector[*ev]->propensity;
	
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

double SLeaping_LacZLacY::computeTimeOfCritical(vector<int> criticalReactions, double& ac0) 
{
	double tau2;
	ac0 = 0.0;
	
	for (vector<int>::iterator ev=criticalReactions.begin(); ev!=criticalReactions.end(); ++ev) 
		ac0 += eventVector[*ev]->propensity;
	
	if (criticalReactions.empty()) 
		tau2 = HUGE_VAL;
	else
		tau2 = (1.0/ac0) * sgamma( (double)1.0 );		
	
	return tau2;
};


void SLeaping_LacZLacY::sampling(long int Llocal, double a0, double ac0, vector<int> criticalReactions, vector<int>nonCriticalReactions, short int crit)
{
	long int numberOfReactions = sbmlModel->getNumReactions(); 
	
	if (criticalReactions.empty())
	{
		if ( Llocal > numberOfReactions )
		{
			double p = 0.0;
			double cummulative	= a0;
			long int k			= 0;
			
			for (int ev=0; ev < eventVector.size(); ++ev)
			{				
				if( (ev == numberOfReactions - 1 ) && (Llocal != 0) ){
					fireReactionProposed( eventVector[ev]->index, Llocal);	
					break;
				}
				
				cummulative		-= p;
				p				 = eventVector[ev]->propensity;
				if (p!=0)
				{
					k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
					Llocal			-= k;
					
					fireReactionProposed( eventVector[ev]->index , k);	
					if (Llocal == 0){ break; }
				}
			}
		}
//		else 
//		{
//			for (int s = 0; s<Llocal; s++){
//				double r1 = ranf();
//				int ev = 0;
//				double suma = eventVector[ev]->propensity;
//				
//				while ( suma < r1 * a0){
//					ev++;
//					suma += eventVector[ev]->propensity;
//				}
//				
//				fireReactionProposed(eventVector[ev]->index,1);
//			}
//		}
	}
	else  // sample just from non critical
	{
		if (Llocal > numberOfReactions)
		{
			double p			= 0.0;
			double cummulative	= a0-ac0;
			long int k			= 0;
			
			for (vector<int>:: iterator ev=nonCriticalReactions.begin(); ev!=nonCriticalReactions.end(); ++ev)
			{				
				if ( (*ev == nonCriticalReactions.size()-1) && (Llocal!=0) ) {
					fireReactionProposed( eventVector[*ev]->index, Llocal);	
					break;
				}
				
				cummulative		-= p;
				p				 = eventVector[*ev]->propensity;
				if (p!=0)
				{
					k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
					Llocal			-= k;
					
					fireReactionProposed( eventVector[*ev]->index , k);	
					if (Llocal == 0){ break; }
				}
			}
		}
//		else
//		{
//			double anc0 = a0 - ac0;
//			
//			for (int s=0; s< Llocal; ++s)
//			{
//				double r1 = ranf();
//				vector<int>::iterator ev = nonCriticalReactions.begin();
//				double suma = eventVector[*ev]->propensity;
//								
//				while (suma < r1 * anc0) {
//					++ev;
//					suma += eventVector[*ev]->propensity;
//				}
//				
//				fireReactionProposed(eventVector[*ev]->index,1);
//			}
//		}
	}

	
	if (crit == 1)
	{
		double r1 = ranf();
		int reactionIndex = 0;
		vector<int>::iterator it = criticalReactions.begin();
		double cummulative = eventVector[*it]->propensity;
		
		while (cummulative < r1 * ac0){
			++it;
			cummulative += eventVector[*it]->propensity;
		}
		
		fireReactionProposed( eventVector[*it]->index , 1 );
	}
}

void SLeaping_LacZLacY::_executeSSA(double& t, int SSAsteps, double genTime)
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
		computePropensitiesGrowingVolume(genTime);
		a0 = blitz::sum(propensitiesVector);
		tau = (1.0/a0) * sgamma( (double)1.0 );
		
		r1 = ranf();
		reactionIndex = -1;
		cummulative = 0.0;
		
		for (int ev=0; ev < eventVector.size(); ++ev)
		{
			cummulative += eventVector[ev]->propensity;
			if ( cummulative > a0*r1 )
			{
				reactionIndex = eventVector[ev]->index;
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

long int SLeaping_LacZLacY::computeLeapLength(double dt, double a0, double ac0, vector<int> criticalReactions)
{	
	if ( criticalReactions.empty() ) 
		return (long int)max( (long int)ignpoi(a0*dt), (long int)1) ; 
	else
		return (long int)max( (long int)ignpoi( (a0 - ac0) * dt), (long int)1) ; 
}


void SLeaping_LacZLacY:: _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum)
{
	double aver_dt = dt_sum/ (double) steps;
	double aver_L = L_sum/ (double) steps;
	
	// time, dt, aver_dt, L,averL, #iterations
	if (myfile!=NULL)
	{
		fprintf(myfile, "%f  %e  %e  %i  %f  %i \n", t, dt, aver_dt, L, aver_L, steps  );
	}
}

void SLeaping_LacZLacY:: _writeTrajectories(FILE* myfile, vector<int> data)
{
	if (myfile!=NULL) 
	{
		for (int i=0; i < data.size() ; ++i) 
			fprintf(myfile, "%i ", data[i], " " );
		
		fprintf(myfile, "\n");
	}
}

void SLeaping_LacZLacY::solve()
{
#ifndef NDEBUG	
	cout << "SLeaping LacZ/LacY..." << endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	FILE* rejectionsfile = fopen("Rejections_S.txt", "w");
#endif
	
#ifdef TRAJECTORIES
	FILE* file_RbsLacY   = fopen("S_RbsLacY.txt"  , "w");
	FILE* file_TrLacZ2   = fopen("S_TrLacZ2.txt"  , "w");
	FILE* file_TrRbsLacY = fopen("S_TrRbsLacY.txt", "w");
	FILE* file_TrRbsLacZ = fopen("S_TrRbsLacZ.txt", "w");
#endif
	
	double a0						= 0.0;
	double ac0						= 0.0;
	double dt2						= 0.; 
	long int L						= 1;
	short int crit					= 0;
	bool isNegative					= false;
	double averNumberOfRealizations = 0.0;
	double numberOfRejections = 0.0;
	
	double genTime = 2100.;   // generation time
	
#ifndef NSSASTEP
	const int SSAfactor = 10;
	const int SSAsteps = 100;
#endif
	
	for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
	{
		Event * e = new Event();
		e->index		= i;
		e->propensity	= 0.0;
		eventVector.push_back(e);
	}
	
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
		FILE* myfile = fopen("Diag_S.txt", "w");
		double whenToWriteOffset = tEnd / numberOfFrames;
		double whenToWrite = whenToWriteOffset;
		int iterSteps = 0;
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
			
			computePropensitiesGrowingVolume(genTime);
			a0 = blitz::sum(propensitiesVector);
			
			if (numberOfIterations % simulation->SortInterval == 0)
				sort(eventVector.begin(), eventVector.end(), EventSort());
			
			vector<int> criticalReactions;
			vector<int> nonCriticalReactions;
			identifyReactions(criticalReactions,nonCriticalReactions);
			crit = 0;
			
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
#endif
			{
				dt2 = computeTimeOfCritical(criticalReactions, ac0);
				if (dt2 < dt) {
					dt	 =	dt2;
					crit = 1;
				}
				
				L =  computeLeapLength(dt,a0,  ac0, criticalReactions);	
				sampling(L, a0, ac0, criticalReactions, nonCriticalReactions, crit);
				
				if (isProposedNegative() == false)
				{
					acceptNewSpeciesValues();
					++numberOfIterations;
					t += dt;
					isNegative = false;
				
					
			#ifdef TRAJECTORIES
					// write trajectories of RbsLacY=S(7), TrLacZ2=S(5),TrRbsLacY=S(13),TrRbsLacZ=S(12)
					if (t >= whenToWriteTrajectories ){
						RbsLacY.push_back(simulation->speciesValues(7));
						TrLacZ2.push_back(simulation->speciesValues(5));
						TrRbsLacY.push_back(simulation->speciesValues(13));
						TrRbsLacZ.push_back(simulation->speciesValues(12));
						
						whenToWriteTrajectories += whenToWriteTrajectoriesOffset;
					}
			#endif

			#ifndef NDIAGNOSTIC
					LAverage += (double)L;
					iterSteps = iterSteps + 1;
					L_sum  += L;
					dt_sum +=dt;
					
					if ( t >= ((double) (whenToWrite))){
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
					dt = dt *0.5;
					reloadProposedSpeciesValues();
					isNegative = true;
				}
			}  // end of SSA
		}     // end of while

		
#ifdef TRAJECTORIES
		_writeTrajectories(file_RbsLacY,RbsLacY);
		_writeTrajectories(file_TrLacZ2,TrLacZ2);
		_writeTrajectories(file_TrRbsLacY,TrRbsLacY);
		_writeTrajectories(file_TrRbsLacZ,TrRbsLacZ);	
		
		cout << "Sample: " << samples << endl;
#endif
		
#ifndef NDEBUG 		
		saveData();
		
#ifndef NDIAGNOSTIC
		cout << "Average L: " << LAverage/((double)numberOfIterations) << endl;
		fclose(myfile);
#endif
		writeToAuxiliaryStream( simulation->speciesValues );
		//writeData(localOutputFileName,samples);
		averNumberOfRealizations += numberOfIterations;
		
		// report on negative population
		if (rejectionsfile!=NULL)
			fprintf(rejectionsfile, "%i %f %i \n", samples, numberOfRejections ,numberOfIterations);
#endif
		
	}  // end of : for(samples, ...) 
	
#ifndef NDEBUG
	writeData(outputFileName);
	closeAuxiliaryStream();
	cout << " Average number of Realizations in R-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
	
	fclose(rejectionsfile);
#endif
	
#ifdef TRAJECTORIES
	fclose(file_RbsLacY);
	fclose(file_TrLacZ2);
	fclose(file_TrRbsLacY);
	fclose(file_TrRbsLacZ);
#endif
	
	for (int i = 0; i < eventVector.size(); ++i) { delete eventVector[i]; }
}






