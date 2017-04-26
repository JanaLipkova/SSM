/*
 *  RLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/27/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "RLeaping.h"

RLeaping::RLeaping(Simulation * simulation):
LeapMethod(simulation)
{
}

RLeaping::~RLeaping()
{
}


/*
 * returnval: L
 * constraint: L>0
 */
long int RLeaping::computeLeapLength()
{
	long int Lprime			 = 2147483647;//2^{31}-1
	long int L				 = 2147483647;
	
	double theta	= simulation->Theta;
	double epsilon	= simulation->Epsilon;
	
	// L prime first
	
	int numberOfSpecies		= sbmlModel->getNumSpecies();
	int numberOfReactions	= sbmlModel->getNumReactions(); 
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
	
	//cout << "a0 in computeLenghL: " << a0 << endl;
	
	for (int is = 0; is < numberOfSpecies; ++is)
	{
		taup = (HUGE_VALF*0.5);
		xi = (double)simulation->speciesValues(is);
		
		switch (hor(is)) {
			case 0:
				//cout << "case0: " << endl;
				break;
			case 1:
				//		cout << "case1: " << endl;
				
				epsi = epsilon;
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabsf(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				break;
			case 2:
				//		cout << "case2: " << endl;
				
				if (nuHor(is) == 1)
					epsi = 0.5*epsilon;
				else
					epsi = epsilon*(xi-1.0)/(2.0*(xi-1.0)+1.0);
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabs(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				
				//cout << "muHat(is): " << muHat(is) << endl;
				//cout << "varHat(is): " << varHat(is) << endl;
				
				break;
			case 3:
				//	cout << "case3: " << endl;
				
				if (nuHor(is)==1)
					epsi = 0.3333333333*epsilon;
				else if (nuHor(is) == 2)
					epsi = epsilon*(xi-1)/(3.0*(xi-1)+1.5);
				else
					epsi = epsilon*(xi-1)*(xi-2)/(3.0*(xi-1)*(xi-2)+(xi-2)+2.0*(xi-1));
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabs(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				break;
			default:
				break;
				
		}
		
		
		//cout << "tau * a0: " << (tau*a0) << endl;
		
	}
	Lprime = (long int)max((long int)(tau*a0), (long int)1);
	
	L = Lprime;
	//cout << "L before theta; " << L << endl;
	
	for (int ir = 0; ir < numberOfReactions; ++ir)
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
				//cout << "lj: " << lj << endl;
				//cout << "other: " << -simulation->speciesValues(changes[is])/ nuChanges[is] << endl;
				
			}
			long int propsedL = (long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj);
			if (propsedL < L && propsedL > 0)
			{
				L = propsedL;
			}
			//L = min( (long int)L, (long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj) );
			//cout << "	proposed	L: " << propsedL << endl;
			//cout << "	compared to L: " << L << endl;//(long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj) << endl;
		}
	}
	//cout << "L final: " << L << endl;
	//abort();
	
	assert(L > 0);
	
	return L;
}


// this method is overloaded from the Methods class since R-Leaping needs to store both indices and 
// propensities (not just propensities).  These are located in the anonymous inner class called Event
void RLeaping::computePropensities()
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
			//original bbayati
			
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
		//cout << "a" << ir << endl;;
		propensitiesVector(ir)		= reaction->getPropensity();
		eventVector[ev]->propensity	= reaction->getPropensity();
	}
}

void RLeaping::solve()
{
	cout << "RLeaping..." << endl;
	
	double a0			= 0.0;
	double p			= 0.0;
	double cummulative	= a0;
	long int k			= 0;
	long int Llocal		= 1;
	long int Lcurrent	= 1;
	bool isNegative		= false;
	
	for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
	{
		Event * e = new Event();
		e->index		= i;
		e->propensity	= 0.0;
		eventVector.push_back(e);
	}
	
	
	for (int noise = 0; noise < simulation->NumberOfNoiseLevels; noise++)
	{
		double noiseLevel = simulation->InitialNoise+noise*simulation->NoiseIncrement;
		int stringLength = outputFileName.length();
		string localOutputFileName = outputFileName.substr(0, stringLength-4);
		localOutputFileName = localOutputFileName + "-N"+ boost::lexical_cast<std::string>(noise) + "-S.txt";
		
		openAuxiliaryStream( (simulation->ModelName) + boost::lexical_cast<std::string>(noise) + "-histogram-rLeaping.txt");
		
		
		
		for (int samples = 0; samples < numberOfSamples; ++samples)
		{
			t = simulation->StartTime;
			numberOfIterations = 0;
			timePoint = 0;
			zeroData();
			simulation->loadInitialConditions(noiseLevel);
			Llocal = 1;
			Lcurrent = 1;
			isNegative = false;
			double LAverage = 0.0;
			
			double meanJ = 0.0;
			
			while (t < tEnd)
			{
				computePropensities();
				a0 = blitz::sum(propensitiesVector);
				
				// sort the list
				if (numberOfIterations % simulation->SortInterval == 0)
				{
					sort(eventVector.begin(), eventVector.end(), EventSort());
					cout << "			Sorting propensities." << endl;
				}
				
				if (isNegative == false)
				{ Lcurrent =  computeLeapLength(); }
				Llocal = Lcurrent;
				p = 0.0;
				cummulative	= a0;
				int j = 0;
				for (j = 0; j < eventVector.size(); ++j)
				{				
					cummulative		-= p;
					p				 = eventVector[j]->propensity;
					k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
					Llocal			-= k;
					
					fireReactionProposed( eventVector[j]->index , k);
					
					if (Llocal == 0){ break; }
				}
				meanJ += ((double)j);
				
				if (isProposedNegative() == false)
				{
					acceptNewSpeciesValues();
					++numberOfIterations;
					dt = (1.0/a0) * sgamma( (double)Lcurrent ); // Gamma ( L, 1.0 / a0 )
					t += dt;
					isNegative = false;
					LAverage += (double)Lcurrent;
				}
				else
				{
					cout << "Negative species at time: " << t << endl;
					Lcurrent = max( (int)(Lcurrent / 2.0), 1);
					cout << "Lcurrent = " << Lcurrent << endl;
					reloadProposedSpeciesValues();
					isNegative = true;
				}
				
			}
			saveData();
			cout << "Sample: " << samples << endl;
			cout << "Average L: " << LAverage/((double)numberOfIterations) << endl;
			cout << "Mean J: " << meanJ/((double)numberOfIterations) << endl;
			writeToAuxiliaryStream( simulation->speciesValues );
			writeData(localOutputFileName,samples);
		}
		//	writeData(outputFileName);
		closeAuxiliaryStream();
		
	}
	
	for (int i = 0; i < eventVector.size(); ++i) { delete eventVector[i]; }
}

