/*
 *  SLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 4/25/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#include "SLeaping.h"

SLeaping::SLeaping(Simulation * simulation):
LeapMethod(simulation)
{ }

SLeaping::~SLeaping()
{ }

double SLeaping::computeTimeStep()
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
}

long int SLeaping::computeLeapLength(double dt, double a0)
{
	long int	L = (long int)max( (long int)ignpoi(a0*dt), (long int)1) ; 
	assert(L > 0);
	return L;
}

// this method is overloaded from the Methods class since R-Leaping needs to store both indices and 
// propensities (not just propensities).  These are located in the anonymous inner class called Event
void SLeaping::computePropensities()
{
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
			SSMReaction* reaction		= ssmReactionList[ir];
			vector <int>  reactants		= reaction->getReactants();
			vector <int>  nu_reactants	= reaction->getNuReactants();
			
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
		}
		propensitiesVector(ir)		= reaction->getPropensity();
		eventVector[ev]->propensity	= reaction->getPropensity();
	}
}




void SLeaping::sampling(long int Llocal, double a0)
{
    double p = 0.0;
    double cummulative	= a0;
    long int k			= 0;
    
    for (int j = 0; j < eventVector.size(); ++j){
        if( (j == eventVector.size() - 1 ) && (Llocal != 0) ){
            fireReactionProposed( eventVector[j]->index , Llocal);
            break;
        }
        
        cummulative		-= p;
        p				 = eventVector[j]->propensity;
        
        if (p!=0){
            k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
            Llocal				-= k;
            
            fireReactionProposed( eventVector[j]->index , k);
            if (Llocal == 0){ break; }
        }
    }
    
}

void SLeaping::_executeSSA(double& t, int SSAsteps)
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
		computePropensities();
		a0 = blitz::sum(propensitiesVector);
		tau = (1.0/a0) * sgamma( (double)1.0 );
		
		r1 = ranf();
		reactionIndex = -1;
		cummulative = 0.0;
		for (int j = 0; j < eventVector.size(); ++j)
		{
			cummulative += eventVector[j]->propensity;
			if ( cummulative > a0*r1 )
			{
				reactionIndex = eventVector[j]->index;
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

void SLeaping::	_writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum)
{
	double aver_dt = dt_sum/ (double) steps;
	double aver_L = L_sum/ (double) steps;
	
	if (myfile!=NULL)
		fprintf(myfile, "%f  %e  %e  %i  %f  %i \n", t, dt, aver_dt, L, aver_L, steps  );
}


void SLeaping::solve()
{
	cout << "SLeaping..." << endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	FILE* rejectionsfile = fopen("Rejections_S.txt", "w");
	
	double a0						= 0.0;
	long int L						= 1;
	bool isNegative					= false;
	double averNumberOfRealizations = 0.0;
	double numberOfRejections;
	
	#ifdef SSASTEP
	const int SSAfactor = 100;
	const int SSAsteps = 1;
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
		numberOfIterations		= 0;
		numberOfRejections		= 0.0;
		timePoint				= 0;
		zeroData();
		simulation->loadInitialConditions();
		L						=	1;
		isNegative				= false;
		
		while (t < tEnd) 
		{
        saveData();
			computePropensities();
			a0 = blitz::sum(propensitiesVector);
			
			if (numberOfIterations % simulation->SortInterval == 0)
				sort(eventVector.begin(), eventVector.end(), EventSort());
			
			if (isNegative == false){
				dt = computeTimeStep();
				if (dt >= HUGE_VAL) {dt= tEnd*10; cout<<"stop"<<endl;	break;}
			}
			
#ifdef SSASTEP
			if (dt <= SSAfactor *  (1.0/a0)  ) 
			{
				_executeSSA(t, SSAsteps);
				SSA_L += (double) SSAsteps;

			}
			else
#endif
			{
				L =  computeLeapLength(dt,a0);	
				sampling(L, a0);
				
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
					++numberOfRejections;
					dt = dt *0.5;
					reloadProposedSpeciesValues();
					isNegative = true;
				}
			}
		}
		
				saveData();

		cout << "Sample: " << samples << endl;
		#ifdef SSASTEP
		cout<< "Number of SSA steps: "<< SSA_L / SSAsteps << endl;
		#endif
		
		writeToAuxiliaryStream( simulation->speciesValues );
		averNumberOfRealizations += numberOfIterations;
		
		// report on negative population
		if (rejectionsfile!=NULL)
			fprintf(rejectionsfile, "%i %f %i \n", samples, numberOfRejections ,numberOfIterations);
		
	}
	

	writeData(outputFileName);
	closeAuxiliaryStream();
    fclose(rejectionsfile);

	cout << " Average number of Realizations in S-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
    
	for (int i = 0; i < eventVector.size(); ++i) { delete eventVector[i]; }
}


