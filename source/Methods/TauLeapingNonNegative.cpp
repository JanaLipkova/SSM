/*
 *  TauLeapingNonNegative.cpp
 *  
 *
 *  Created by Lipkova on 3/5/14.
 *  Copyright 2014 Jana Lipkova All rights reserved.
 *
 */


//****************************
// IDEA OF Non Negative TAU ALGORITHM:
//****************************

// 1. compute propensities
// 2. build list of critical reactions
// 3. compute tau1 excluding critical reactions
// 4. tau2 = time of next critical reactions
// 6. if (tau1 < tau2) => tau = tau1, and no crtical reaction will happen
//    else tau = tau2, => evolve non critical reactions normally, and 
//							the one critical reaction is simulated with SSA
//****************************

#include "TauLeapingNonNegative.h"

TauLeapingNonNegative::TauLeapingNonNegative(Simulation* simulation):
LeapMethod(simulation)
{ }

TauLeapingNonNegative::~TauLeapingNonNegative()
{ }


//****************************
//  build list of critical reactions
//****************************
vector<int> TauLeapingNonNegative::listOfCriticalReactions()
{	
	vector<int> critical_reactions;
	int Ncrit = 10;  // Cao et. al suggest Ncritical 10
	
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
			critical_reactions.push_back(ir);
	}

	return critical_reactions;
}
//****************************

//****************************
//  compute time step
//****************************
double TauLeapingNonNegative::computeTimeStep(vector<int> criticalReactions)
{
	vector<int>    non_critical;    // list of non critical reactions used to compute muHat, signaHat2 
	
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
	{
	    int j = 0;
		for (int ir = 0; ir < numberOfReactions; ++ir)
		{
			if (ir != criticalReactions[j] )   // non critical reactions
				non_critical.push_back(ir);
			else 
				j++; 
		}
		computeMuHatSigmaHat2(muHat, sigmaHat2,non_critical);	
	}
	
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
double TauLeapingNonNegative::computeTimeOfCritical(vector<int> criticalReactions, double& ac0) 
{
	double tau2;
	
	if (criticalReactions.empty()) 
	{
		tau2 = HUGE_VAL;
		ac0 = 0.0;
	}
	else
	{
		ac0 = 0.0;
		
		for (vector<int>::iterator it = criticalReactions.begin(); it!=criticalReactions.end() ; ++it)
			ac0 += propensitiesVector(*it);
		
		tau2 = (1.0/ac0) * sgamma( (double)1.0 );		
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
void TauLeapingNonNegative::computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2, vector<int> non_critical)
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
void TauLeapingNonNegative::executeSSA(double& t, int SSAsteps)
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
void TauLeapingNonNegative::sampling(short int crit, vector<int> criticalReactions, double ac0)
{
	
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	double aj;
	long int kj;
	
	if (criticalReactions.empty()) // no critical reactions, sample from all chanels
	{
		for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
		{				
			aj  = propensitiesVector(j);
            kj  = (aj == 0) ? 0. : ignpoi( aj*dt );
			fireReactionProposed( j , kj );
		}
	}
	else  
	{
		// sample from non-critical
		int k = 0;
		for (int j = 0; j < numberOfReactions; ++j)
		{
			if ( ( j!=criticalReactions[k] )||( k==criticalReactions.size() ) ) // sample non critical reactions
			{
				aj	= propensitiesVector(j);
                kj  = (aj == 0) ? 0. : ignpoi( aj*dt );
                fireReactionProposed( j , kj );
			}
			else  // do not sample critical reactions
				k++;
			
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
	}
}



//****************************

//****************************
//  Solve it :-P
//****************************
void TauLeapingNonNegative::solve()
{
    cout << "Tau Leaping with List of Critical reactions..." <<endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	
	/* initialize helping variables */
	double a0						= 0.0;
	double ac0						= 0.0;
	double dt2						= 0.; 
	short int crit					= 0;
	bool isNegative					= false;
	double averNumberOfRealizations = 0.0;
    vector<int> rejectionsVector(numberOfSamples);
	int numberOfRejections;
	const int SSAfactor = 10;
	const int SSAsteps = 100;
	
	/* execute code for each realization */
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations =  0;
		numberOfRejections	= 0.;
		timePoint			= 0;
		zeroData();
		simulation->loadInitialConditions();
		isNegative			= false;
		
		while (t < tEnd)
		{
			//saveData();
			computePropensities(propensitiesVector, 0);
			a0 = blitz::sum(propensitiesVector);
			
			/* 1) get critical reactions */
			vector<int> criticalReactions = listOfCriticalReactions();
			crit = 0;						// number of critical reactions to be executed;
			
			/* 2) compute time step */
			if (isNegative == false){
		          dt = computeTimeStep(criticalReactions);
			  if (dt >= HUGE_VAL) {dt= tEnd*10; cout<<"stop"<<endl;   break;}
			}
			
			if (dt <= SSAfactor * (1.0/a0) )
				executeSSA(t, SSAsteps);
			else
			{
				/* 3) compute time of critical reaction*/
				dt2 = computeTimeOfCritical(criticalReactions, ac0);
				
				if (dt2 < dt) {
					dt	 =	dt2;
					crit = 1;
				}
				
				/* 4) sampling */
				sampling(crit, criticalReactions,ac0);
								
				/* 5) check proposed update and update the system or repeat simulaiton */
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
					dt = dt * 0.5;
					reloadProposedSpeciesValues();
					isNegative = true;
				}
					
			}  
		}
		
        cout << "Sample: " << samples << endl;

		//saveData();
        rejectionsVector[samples] = numberOfRejections;
		writeToAuxiliaryStream( simulation->speciesValues );
		averNumberOfRealizations += numberOfIterations;
				
	}
	
	//writeData(outputFileName);
	closeAuxiliaryStream();

	cout << " Average number of Realizations in Non Negative Tau-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
    
    int rejectionSum = std::accumulate(rejectionsVector.begin(), rejectionsVector.end(), 0);
    std::cout<<"Negative species appeared in total:" << rejectionSum << " times" << std::endl;
}
