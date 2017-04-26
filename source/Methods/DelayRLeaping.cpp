/*
 *  DelayRLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/3/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "DelayRLeaping.h"

DelayRLeaping::DelayRLeaping(Simulation * simulation):
LeapMethod(simulation)
{
}

DelayRLeaping::~DelayRLeaping()
{
}
  #pragma mark  - RLeaping Methods - 
  
long int DelayRLeaping::computeLeapLength()
{
	long int Lprime			 = 2147483647;
	long int L				 = 2147483647;
	
	double theta	= simulation->Theta;
	double epsilon	= simulation->Epsilon;
		
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
				
		Lprime = (long int)max((long int)(tau*a0), (long int)1);
	}

	L = Lprime;
	
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
		}
	}
	assert(L > 0);
	return L;
}

  #pragma mark  - Delay Methods - 

void DelayRLeaping::computePropensities()
{
	propensitiesVector	= 0.0;
	//int lastIndex		= -1;
	Reaction * reaction;
	int reactionIndex;
	
	for (int ev = 0; ev < eventVector.size(); ++ev)
	{
		eventVector[ev]->propensity = 0.0;
		reactionIndex = eventVector[ev]->index;
		
		propensitiesVector(reactionIndex) = 0.0;
		
		reaction = sbmlModel->getReaction(reactionIndex);
		
		KineticLaw * kineticLaw = reaction->getKineticLaw();
		Parameter * parameter = kineticLaw->getParameter(0);
		double rate = parameter->getValue();
		
		if ( isReactionDelayed(reactionIndex) == true )
		{
			int dependentSpecies = getDependentSpecies(reactionIndex);
			double dependentValue = (double)simulation->speciesValues(dependentSpecies);
			double h =					kineticLaw->getParameter(2)->getValue();
			double defaultProduction =	kineticLaw->getParameter(3)->getValue();
			double cHill =				kineticLaw->getParameter(4)->getValue();
			propensitiesVector(reactionIndex) = defaultProduction + rate*hillFunction(cHill, dependentValue, h);
			//
//			int dependentSpecies = getDependentSpecies(reactionIndex);
//			double h =   kineticLaw->getParameter(2)->getValue();
//			double P0 =  kineticLaw->getParameter(3)->getValue();
//			propensitiesVector(reactionIndex)		= rate*hillFunction((double)simulation->speciesValues(dependentSpecies), P0, h);
			eventVector[ev]->propensity = defaultProduction + rate*hillFunction(cHill, dependentValue, h);//rate*hillFunction((double)simulation->speciesValues(dependentSpecies), P0, h);
		}
		else
		{
			propensitiesVector(reactionIndex)		= rate;
			eventVector[ev]->propensity = rate;
			if (kineticLaw->getFormula().substr(0, 8) == "Species:")
			{
				int dependentSpecies = getDependentSpecies(reactionIndex);
				propensitiesVector(reactionIndex)		*= (double)simulation->speciesValues(dependentSpecies);
				eventVector[ev]->propensity *= (double)simulation->speciesValues(dependentSpecies);
			}
			
			if ( reaction->getNumReactants() == 0 )
			{
			}
			else if ( reaction->getNumReactants() == 1 )
			{
				SpeciesReference * speciesReference = reaction->getReactant(0);
				string speciesName = speciesReference->getSpecies();
				int speciesIndex = simulation->getSpeciesIndex(speciesName);
				double value = (double)simulation->speciesValues(speciesIndex);
				
				propensitiesVector(reactionIndex)		*= value;
				eventVector[ev]->propensity *= value;
			}
			else
			{
				int nu;
				ParticleType x;
				ParticleType num, denom;
				
				int ir = reactionIndex;
				vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;
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
				propensitiesVector(reactionIndex) = reaction->getPropensity();
				eventVector[ev]->propensity       = reaction->getPropensity();
			}
		}
	} 
}

/*
{
	propensitiesVector	= 0.0;
	int lastIndex		= -1;
	Reaction * reaction;
	int reactionIndex;
	
	for (int ev = 0; ev < eventVector.size(); ++ev)
	{
		eventVector[ev]->propensity = 0.0;
		reactionIndex = eventVector[ev]->index;
		
		propensitiesVector(reactionIndex) = 0.0;

		reaction = sbmlModel->getReaction(reactionIndex);
		
		KineticLaw * kineticLaw = reaction->getKineticLaw();
		Parameter * parameter = kineticLaw->getParameter(0);
		double rate = parameter->getValue();

		if ( isReactionDelayed(reactionIndex) == true )
		{
			int dependentSpecies = getDependentSpecies(reactionIndex);
			double h =   kineticLaw->getParameter(2)->getValue();
			double P0 =  kineticLaw->getParameter(3)->getValue();
			propensitiesVector(reactionIndex)		= rate*hillFunction((double)simulation->speciesValues(dependentSpecies), P0, h);
			eventVector[ev]->propensity = rate*hillFunction((double)simulation->speciesValues(dependentSpecies), P0, h);
		}
		else
		{
			propensitiesVector(reactionIndex)		= rate;
			eventVector[ev]->propensity = rate;
			if (kineticLaw->getFormula().substr(0, 8) == "Species:")
			{
				int dependentSpecies = getDependentSpecies(reactionIndex);
				propensitiesVector(reactionIndex)		*= (double)simulation->speciesValues(dependentSpecies);
				eventVector[ev]->propensity *= (double)simulation->speciesValues(dependentSpecies);
			}
			
			if ( reaction->getNumReactants() == 0 )
			{
			}
			else if ( reaction->getNumReactants() == 1 )
			{
				SpeciesReference * speciesReference = reaction->getReactant(0);
				string speciesName = speciesReference->getSpecies();
				int speciesIndex = simulation->getSpeciesIndex(speciesName);
				double value = (double)simulation->speciesValues(speciesIndex);
			
				propensitiesVector(reactionIndex)		*= value;
				eventVector[ev]->propensity *= value;
			}
			else
			{
				SpeciesReference * speciesReference = reaction->getReactant(0);
				string speciesName = speciesReference->getSpecies();
				int speciesIndex = simulation->getSpeciesIndex(speciesName);
				double value = (double)simulation->speciesValues(speciesIndex);
				propensitiesVector(reactionIndex)		*= value;
				eventVector[ev]->propensity *= value;
				lastIndex = speciesIndex;
				
				for (int j = 1; j < reaction->getNumReactants(); ++j)
				{
					SpeciesReference * speciesReference = reaction->getReactant(j);
					string speciesName = speciesReference->getSpecies();
					int speciesIndex = simulation->getSpeciesIndex(speciesName);
					double value = (double)simulation->speciesValues(speciesIndex);
					
					if (speciesIndex == lastIndex)
					{
						propensitiesVector(reactionIndex)		*= (value - 1.0);
						eventVector[ev]->propensity *= (value - 1.0);
					}
					else
					{
						eventVector[ev]->propensity *= value;
						propensitiesVector(reactionIndex)		*= value;
					}
					lastIndex = speciesIndex;
				}
			}
		}
	} 
}
 */

double DelayRLeaping::randomVariate(double a0, long int L)
{
	return ((1.0/a0) * sgamma( (double)L )); // Gamma ( L, 1.0 / a0 )
}

void DelayRLeaping::solve()
{
	cout << "Delay-RLeaping..." << endl;
	
	int numDelayed = numberOfDelayedReactions();
	if (numDelayed == 0)
	{
		cout << "You are using the Delay-RLeaping method but have not specified a reaction that has a delay." << endl;
		cout << "I shall now quit and give you some time to think about this." << endl; 
		return;
	}

	delayedReactionsIndices.resize( numberOfDelayedReactions() );
	
	setDelayedReactionsIndices();

	double a0			= 0.0;
	int delayIndex		= 0;
	double p			= 0.0;
	double cummulative	= a0;
	long int k			= 0;
	long int Llocal		= 1;
	long int Lcurrent	= 1;
	
	Array<int, 1> numberToFireForEveryChannel(sbmlModel->getNumReactions());
	numberToFireForEveryChannel = 0;
	
	for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
	{
		Event * e = new Event();
		e->index		= i;
		e->propensity	= 0.0;
		eventVector.push_back(e);
	}
	
	double tauAverage = 0.0;
	
	vector <double> queueTau;
	
	openAuxiliaryStream( (simulation->ModelName) + "-histogram-DRLeaping.txt");
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations	= 0;
		timePoint			= 0;
		zeroData();

		//simulation->numberOfSamples();

		delayedReactionsTimePoints.clear();
		delayedReactionsTimeIndices.clear();
		delayedReactionsTimeLeapLength.clear();
		
		double minValueDTDQ = 1.0;
		vector <double> meanValueDTDQ;
		double maxValueDTDQ = 0.0;
		queueTau.clear();
		
		vector <double> timeSteps;
		
		
		while (t < tEnd)
		{
			saveData();
			
			// compute the propensities
			
			computePropensities();
			a0 = blitz::sum(propensitiesVector);
			Lcurrent =  computeLeapLength(); 
			dt = randomVariate(a0, Lcurrent);
			
			if (dt <= (5.0/a0))
			{
				cout << "dt is too small: " << endl;
				cout << "dt: " << dt << endl;
				cout << 5.0/a0 << endl;
			}
			
			numberToFireForEveryChannel = 0;
			
			while ( isDelayedReactionScheduled(t, dt, delayIndex) == true ) // fire delay that is scheduled within this timestep
			{
				double qt			= delayedReactionsTimePoints[delayIndex];
				double span			= queueTau[delayIndex];
				long long int kd	= delayedReactionsTimeLeapLength[delayIndex];
				
				double ratio = (t + dt - qt) / span;
				meanValueDTDQ.push_back( ratio );
				
				
				minValueDTDQ = min( minValueDTDQ, ratio );
				maxValueDTDQ = max( maxValueDTDQ, ratio );
				
				long int partialAmountFired = ignbin( kd,  min( ratio, 1.0 ) );
				
				fireReactionProposed( delayedReactionsTimeIndices[delayIndex] , partialAmountFired );
				acceptNewSpeciesValues();
				
				// do not remove the delayed reaction, but do give it a new queued dt and K
				queueTau[delayIndex]					   -= (t + dt - qt);     	// - update span
				delayedReactionsTimeLeapLength[delayIndex] -= partialAmountFired;  // - update kd
				delayedReactionsTimePoints[delayIndex]		= t + dt; // it's still up - update qt
				
				if (delayedReactionsTimeLeapLength[delayIndex] == 0)
				{
					removeDelayedReactionsTime(delayIndex);
					queueTau.erase( queueTau.begin() + delayIndex, queueTau.begin()+delayIndex+1);
				}
				
			}
			
			// sort the list
			
			if (numberOfIterations % simulation->SortInterval == 0)
			{
				sort(eventVector.begin(), eventVector.end(), EventSort());
				cout << "			Sorting propensities." << endl;
			}
			
			
			Llocal = Lcurrent;
			p = 0.0;
			cummulative	= a0;
			
			for (int j = 0; j < eventVector.size(); ++j)
			{				
				cummulative		-= p;
				p				 = eventVector[j]->propensity;
				k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
				Llocal			-= k;
				
				numberToFireForEveryChannel( eventVector[j]->index ) =  k;
				
				if (Llocal == 0){ break; }
			}
			
			for (int reactionIndex = 0; reactionIndex < numberToFireForEveryChannel.extent(firstDim); ++reactionIndex)
			{
				if ( isReactionDelayed(reactionIndex) == false ) // not delayed
				{
					fireReactionProposed(reactionIndex, numberToFireForEveryChannel(reactionIndex));
				}
				else if ( numberToFireForEveryChannel(reactionIndex) > 0 ) // queue the delay
				{
					setDelayedReactionsTime(reactionIndex, t, 0.0*dt, numberToFireForEveryChannel(reactionIndex) );
					queueTau.push_back(dt);
				}
				else {}
			}
			
			acceptNewSpeciesValues();
			++numberOfIterations;
				
			t += dt;
			timeSteps.push_back(dt);
		}
		
		saveData();
		writeToAuxiliaryStream( simulation->speciesValues );
		writeSTLVector( timeSteps, sbmlModel->getName() + "time-steps-r.txt");
		
		cout << "minValueDTDQ: " << minValueDTDQ << endl << endl;
		cout << "maxValueDTDQ: " << maxValueDTDQ << endl << endl;
		double mv = 0.0;
		for (int i = 0; i < meanValueDTDQ.size(); ++i)
		{
			mv += meanValueDTDQ[i];
		}
		mv /= meanValueDTDQ.size();
		cout << "meanValueDTDQ: " << mv << endl << endl;
		
		cout << "Sample: " << samples << endl;
		cout << "Average tau: " << tauAverage / ((double)numberOfIterations) << endl;
	}
	writeData(outputFileName);
	closeAuxiliaryStream();
	for (int i = 0; i < eventVector.size(); ++i)
	{
		delete eventVector[i];
	}
}
