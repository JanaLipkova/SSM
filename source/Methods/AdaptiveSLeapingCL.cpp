/*
 *  AdaptiveSLeapingCL.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 7/3/13.
 *  Copyright 2013 CSE Lab. All rights reserved.
 *
 */



//****************************
// IDEA OF ADAPTIVE S ALGORITHM:
//****************************

// 1. compute propensities
// 2. build list of critical reactions
// 3. compute tau1 with PEC and type of stifness
// 4. compare with SSA and solve SSA
// 5. tau2 = time of next critical reactions
// 6. choose tau from tau1 and tau2, based on type of stiffness do sampling k_j
// 7. update state and time, and check for negative population

// Note: This S-leaping method use the control of negative species 
// introduced for Adaptive tau 
//****************************




#include "AdaptiveSLeapingCL.h"
//#define PEC

// Class Constructor
AdaptiveSLeapingCL::AdaptiveSLeapingCL(Simulation* simulation):
LeapMethod(simulation)
{
}
// Destuctor
AdaptiveSLeapingCL::~AdaptiveSLeapingCL()
{
}

//****************************
//  build list of critical reactions
//****************************
vector<int> AdaptiveSLeapingCL::listOfCriticalReactions()
{
	int numberOfSpecies		= sbmlModel->getNumSpecies();
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	
	vector<int> critical_reactions;
	int Ncrit = 10;  
	
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
				
			}
		}
		
		if (lj < Ncrit)
			critical_reactions.push_back(ir);
	}
	
	return critical_reactions;
}

//****************************
//  compute tau1 with PEC and type of stifness
//****************************
double AdaptiveSLeapingCL::computeAdaptiveTimeStep(vector<int> critical, int& type, double& tau_expl)
{
	double tauPrimeExp;			// explicit tau
	double tauPrimeImp;         // implicit tau
	double epsilon	= simulation->Epsilon;
	double delta = simulation->Delta;  
	
	//vector<int> reversible = simulation -> Reversible;  // indecies of reverisble reactions
	vector<int> rev;	
	rev.push_back(1);
	rev.push_back(2);
	
	list<int>    non_critical;    // list of reactions used to compute muHat, signaHat2 
	vector<int>  in_pec;          // list of reactions in PEC
	int Nstiff = 100;       // recomended in Cao at all, "The Adaptive Explicit-Implict Tau-Leaping Method with Automatic tau selection"   
	
	int numberOfSpecies		= sbmlModel->getNumSpecies();
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	
	
	Array<int, 1> hor			(numberOfSpecies);
	Array<int, 1> nuHor			(numberOfSpecies);
	Array<double, 1> muHat		(numberOfSpecies);
	Array<double, 1> sigmaHat2	(numberOfSpecies);
	Array<double, 1> varHat		(numberOfSpecies);
	
	hor = 0; nuHor = 0;
	computeHor(hor, nuHor);
	
	/* 1. STEP: EXPLICIT TAU */
	//-----------------------
	// 1. remove critical reactions
	// 2. compute muHat, sigmaHat2 for corresponding reactions
	// 3. compute explicit tau
	
    muHat = 0.0; sigmaHat2 = 0.0;
	
	if (critical.empty())
	{LeapMethod::computeMuHatSigmaHat2(muHat, sigmaHat2);
	}
	else 
	{
	    int j = 0;
		for (int ir = 0; ir < numberOfReactions; ++ir)
		{
			if (ir != critical[j] )   // non critical reactions
			{non_critical.push_back(ir);}
			else { j++; }
		}
		
		computeMuHatSigmaHat2(muHat, sigmaHat2,non_critical);	
	}
	
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
	tauPrimeExp = tau;
	tau_expl = tau;
	
#ifdef PEC
	/* 2. STEP IMPLICIT TAU */
	//------------------------
	// 1. find reactions in PEC
	// 2. build list of reactions considered in muHat,sigmaHat2 computation
	// 3. compute implicit tau
	
	double temp1;  // to store propenisty of forward reactions
	double temp2;  // to store propensity of backward reactions
	
	if (rev.empty()) 
	{
		tauPrimeImp = tauPrimeExp;
	}
	else 
	{	
		for (vector<int>::iterator it = rev.begin(); it !=rev.end(); it=it+2)
		{
			temp1 = propensitiesVector(*it);    // propensity of forward reaction
			temp2 = propensitiesVector(*it+1);  // propensity of backward reacrion
			
			if ( abs( temp1 - temp2 ) < delta* abs( temp1 + temp2 ) )
 			{
				in_pec.push_back(*it);
				in_pec.push_back(*it+1);
			}
		}
		
		if (in_pec.empty())
		{
			tauPrimeImp = tauPrimeExp;
		}
		else
		{
			if (critical.empty())
			{
				int k =0;
				for (int ir =0; ir < numberOfReactions ; ++ir)
				{
					if (ir != in_pec[k]) 
					{
						non_critical.push_back(ir);
					}
					else { k++; }
					
				}
			}
			else 
			{
				for (vector<int>::iterator it = in_pec.begin(); it !=in_pec.end(); ++it)
				{
					non_critical.remove(*it);
				}
			}
			
			muHat = 0.0; sigmaHat2 = 0.0;
			computeMuHatSigmaHat2(muHat, sigmaHat2,non_critical);
			
			double tau, taup,  epsi, epsixi, epsixisq;
			double xi;
			
			tau = HUGE_VAL;
			
			double a0 = (double)blitz::sum(propensitiesVector);
			for (int is = 0; is < numberOfSpecies; is++)
			{
				varHat(is) = sigmaHat2(is);// - (1.0/a0) * muHat(is) * muHat(is);
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
			
			tauPrimeImp = tau;
			
		}
	}
	
	
	/* 3. TAU AND STIFFNES */
	if (tauPrimeImp > Nstiff * tauPrimeExp)
	{
		type = 1;       // Stiff system take implicit tau
		tau = tauPrimeImp; 
	}
	else
	{
		type = 0;      // non-stiff system and explicit tau
		tau = tauPrimeExp;
	}
#else
	type = 0;
	tau = tauPrimeExp;
#endif
	
	return tau;
	
}
//****************************


//****************************
//  muHat, sigmaHat2 computation:
//****************************
// this method is overloaded from the LeapMethod class, since in adaptive tau leaping we exclude
// critical reactions from computation of explicit tau and additionaly we exclude reactions in PEC
// in computation of implicit tau
void AdaptiveSLeapingCL::computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2, list<int> non_critical)
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
	for (list<int>::iterator ir = non_critical.begin(); ir!=non_critical.end(); ++ir)	
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
// if previous step was SSA or explicit execute 100 steps of SSA
// if previous step was implicit execute 10 steps of SSA 
// number of steps recomended in Cao at all, "The Adaptive Explicit-Implict Tau-Leaping Method with Automatic tau selection"   
void AdaptiveSLeapingCL::execute_SSA(int& type, double& t, int& numberOfIterations)
{
	double a0 = 0.0;
	double dt;
	double r1;
	int reactionIndex = 0;
	double cummulative = 0.0;
	int steps;
	int count = 0;
	double time = 0.0;
	
	vector<int> k(sbmlModel->getNumReactions(),0);
	
	
	//check type of previus time step, type = 0 is for explicit or SSA, type = 1 is for implicit
	if (type == 0)
	{
		steps = 100;
		type = 0;      
	}
	else 
	{
		steps = 10;
		type = 0;     
	}
	
	while (count < steps)
	{
		count++;
		computePropensities();
		a0 = blitz::sum(propensitiesVector);
		
		dt = (1.0/a0) * sgamma( (double)1.0 );
		
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
			++numberOfIterations;
			time += dt;
		}
		else { 
			count = steps;
			t = HUGE_VAL;
			//cout << " end of simulation in SSA"<<endl;
		}
	}
	
	t += time;
}
//****************************

//****************************
//     time of next critical reactions
//****************************
double AdaptiveSLeapingCL::time_of_next_critical_reaction(vector<int> critical)
{
	double tau2;
	
	if (critical.empty()) 
	{ tau2 = HUGE_VAL;}
	else
	{
		double ac0 = 0.0;
		
		for (vector<int>::iterator it = critical.begin(); it!=critical.end() ; ++it)
		{
			ac0 += propensitiesVector(*it);
		}
		
		tau2 = (1.0/ac0) * sgamma( (double)1.0 );		
	}
	
	return tau2;
}
//****************************


//****************************************
// Final Tau, L and corresponding sampling
//****************************************
double AdaptiveSLeapingCL::tauL_sampling(double tau1, double tau2, double tau_expl, int& type,vector<int> critical)
{
	double tau;
	long int L;
	int crit_count = 0;  // count number of critical reactions that needs to be fired
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	
	
	// choose final tau and type of method
	if (tau1 < tau2) 
	{
		tau = tau1;
		crit_count = 0;
	}
	else
	{
		tau = tau2;
		crit_count = 1;   // one crit. reaction must be executer
		
		if ((type ==0)||(type==1 && tau2 < tau_expl) )
		{ type = 0; }
		else{ type = 1; }
	}
		
	// compute L w.r.t to type of method, and sample :)
	if (type == 0)
	{
		double a0 = (double)blitz::sum(propensitiesVector);
		
		if (critical.empty())
		{
			L = (long int)max( (long int)ignpoi(a0*tau), (long int)1);
			//long int L = (long int)max( (long int)(a0*tau), (long int)1);
			sampling(L,a0);
		}
		else 
		{
			double a0nc = 0.0; // sum of propensitirs of non-critical reactiosn
			double ac0 = 0.0;  // sum of propensites of critical reactions
			
			for (vector<int>::iterator it = critical.begin(); it != critical.end(); ++it){
				ac0 += propensitiesVector(*it);
			}
			
			a0nc = a0 - ac0;
			
			L = (long int)max( (long int)ignpoi(a0nc*tau), (long int)1);
			//long int L = (long int)max( (long int)(a0nc*tau), (long int)1);
			sampling_critical(L, critical, a0nc,crit_count,ac0);
		}
		
	}
	else
	{
		implicit_sampling(tau, critical); 		
	}
	
	return tau;
	
}
//****************************


/* "experience in this lab is when chaos is created every once in a while 
    and you can handle it" */




//***********************************
//       Sampling 
// Binomial sampling from reordered propensities
// no critical reactions considered
//***********************************
void AdaptiveSLeapingCL::sampling(long int L, double a0)
{
	
	double p = 0.0;
	double cummulative	= a0;
	long int k			= 0;
	
//	int j = 0;
//	for (j = 0; j < eventVector.size(); ++j)
//	{				
//		cummulative		-= p;
//		p				 = eventVector[j]->propensity;
//		k				 = ignbin(L, min(p/cummulative, 1.0) );
//		L				-= k;
//		
//		fireReactionProposed( eventVector[j]->index , k);			
//		if (L == 0){ break; }
//	}

	if ( L > sbmlModel->getNumReactions() )
	{
		double p = 0.0;
		double cummulative	= a0;
		long int k			= 0;
		
		int j = 0;
		for (j = 0; j < eventVector.size(); ++j)
		{				
			cummulative		-= p;
			p				 = eventVector[j]->propensity;
			k				 = ignbin(L, min(p/cummulative, 1.0) );
			L				-= k;
			
			fireReactionProposed( eventVector[j]->index , k);			
			if (L == 0){ break; }
		}
	}
	else {
		for (int s = 0; s<L; s++)
		{
			double r1 = ranf();
			int j = 0;
			double suma = propensitiesVector(j)/a0;
			
			while (r1 > suma)
			{
				j++;
				suma += propensitiesVector(j)/a0;
			}
			
			fireReactionProposed(j,1);
		}
	}
	
}
//***********************************


//****************************
//         Sampling
//  Binomial sampling from reordered propensities
//  excluding critical reactions
// a0nc = sum of propensities of non critical reactions
// ac0 = sum of propensities of critical reactions
//****************************
void AdaptiveSLeapingCL::sampling_critical(long int L,vector<int> critical, double a0nc, int crit_count, double ac0)
{	
	if (L > sbmlModel->getNumReactions() ) 
	{
		double p			= 0.0;
		double cummulative  = a0nc;
		long int k		    = 0;
		
		for (int j = 0; j < eventVector.size(); ++j)
		{
			if (!iscritical(j,critical) ) 
			{
				cummulative		-= p;
				p				 = eventVector[j]->propensity;
				k				 = ignbin(L, min(p/cummulative, 1.0) );
				L				-= k;
				fireReactionProposed( eventVector[j]->index , k);
				if (L == 0){ break; }
			}
		}
		
	}	
	else
	{
		for (int s=0; s<L; ++s)
		{
			double r1 = ranf();
			double suma = 0.0;
			int j = -1;
						
			while (r1 > suma)
			{
				j++;
				if (!iscritical(j,critical) ){
					suma += propensitiesVector(j)/a0nc;
				}
			}
			
			fireReactionProposed(j,1);
						
		}
	}
	
	
	// if tau = tau2, find index of critical reaction to be fired= reactionIndex	
	if (crit_count==1) 
	{
		int reactionIndex = -1;
		double r1 = ranf();
		double cummulative = 0.0;
		
		for (vector<int>::iterator it = critical.begin(); it != critical.end(); ++it) 
		{
			cummulative += propensitiesVector(*it);
			if ( cummulative > ac0*r1 )
			{
				reactionIndex = *it;
				//cout << "reactionIndex of cricial one " << reactionIndex << endl;
				break;
			}	
		}
		
		if (reactionIndex==-1) {
			cout << "index of critical reaction -1";
			abort();
		}
		
		fireReactionProposed(reactionIndex,1);
	}

}
//***********************************




//**************************************************
// Teste weather index j is in the list if criticla reactions
//**************************************************
bool AdaptiveSLeapingCL::iscritical(int j, vector<int> critical)
{
	bool crit = false;
	for (vector<int>::iterator it = critical.begin(); it!=critical.end(); ++it)
	{
		if (*it == j) {
			crit = true;
			break;
		}
	}
	
	return crit;
}
//**************************************************



void AdaptiveSLeapingCL::implicit_sampling(double tau, vector<int> critical)
{
	double aj;
	long int kj;
	int MaxNumberOfIterations = 100;
	
	int numberOfSpecies		= sbmlModel->getNumSpecies();
	int numberOfReactions	= sbmlModel->getNumReactions(); 
	
	vector<long int> k(numberOfReactions,0);  // sampled reactions,on j-th position is how many times j-th reaction should be fired
	vector<double> aL(numberOfReactions);
	vector<double> B(numberOfSpecies,0);
	vector<double> implicitPropenisty(numberOfReactions,0);   // propenisties in the roots of implicit system
	vector<double> roots(numberOfSpecies,0);                  // to store roots of implict system, denote X^ in literature
	
	
	// STEP 1: precompute k, aj/ao*L, L
	double a0  = (double)blitz::sum(propensitiesVector);
//	long int L = (long int)max( (long int)ignpoi(a0*tau), (long int)1);
	long int L = (long int)max( (long int)(a0*tau), (long int)1);

	
	// sample kj at X(t)
	double p = 0.0;
	double cummulative	= a0;
	long int kk			= 0;
	
	int j = 0;
	for (j = 0; j < eventVector.size(); ++j)
	{				
		cummulative		-= p;
		p				 = eventVector[j]->propensity;
		kk				 = ignbin(L, min(p/cummulative, 1.0) );
		L				-= kk;
		
		k[eventVector[j]->index] = kk;
		
		if (L == 0){ break; }
	}
	
	for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
	{
		aj = propensitiesVector(j);
//		aL[j] = (aj*L)/a0;
		aL[j] = aj*tau;

	}
	
	// build B
	for (int i = 0; i < numberOfSpecies; i++)
		B[i] = simulation->speciesValues(i);
	
	for (int j = 0; j < numberOfReactions; j++)
	{
		SSMReaction * ri = simulation->ssmReactionList[j];		
		const vector<int> & changes = ri->getChanges();
		const vector<int> & nuChanges = ri->getNuChanges();
		
		for (int s = 0; s < changes.size(); s++)
			B[changes[s]] +=nuChanges[s] * (k[j] - aL[j]);
	}
	
	//STEP 2: solve implicit system
	RootFinderJacobian::RootFinderSetUp(simulation, tau, MaxNumberOfIterations, B );
	RootFinderJacobian::find_roots(roots, implicitPropenisty, simulation->speciesValues);
	
	if (critical.empty())
	{
		for (int j=0; j < numberOfReactions; j++) 
		{
			kj = round( implicitPropenisty[j]*tau + k[j] - aL[j] );
			fireReactionProposed( j , kj );
		}
	}
	else
	{
		int m = 0;
		for (int j=0; j<numberOfReactions; j++)
		{
			if(j!= critical[m])
			{
				kj = round( implicitPropenisty[j]*tau + k[j] - aL[j] );
				fireReactionProposed( j , kj );
			}
		}
	}
}


//****************************
//   COMPUTE PROPENSITIES
//****************************
// this method is overloaded from the Methods class since R-Leaping needs to store both indices and 
// propensities (not just propensities).  These are located in the anonymous inner class called Event
void AdaptiveSLeapingCL::computePropensities()
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
//****************************






//****************************
//      SOLVE
//****************************
void AdaptiveSLeapingCL::solve()
{
	cout << "Adaptive S-Leaping with lict of critical reactions..." << endl;
	
	double a0 = 0.0;
	double tau;        // time step
	double tau1;       // 1st tau candidate (from non-critical reactions)
	double tau2;	   // 2nd tau candidate (from critical reactions)	
	double tau_expl;   // explicit tau	int type = 0;      // type = 0 for stiff system, type = 1 for non-stiff one
	int type = 0;      // type = 0 for stiff system, type = 1 for non-stiff one
	int typePrevios   = 0;
	bool isNegative = false;
	double averNumberOfRealizations = 0.0;
	double averNumberOfNegative = 0.0;
	vector<int> critical;           // list of critical reactions
	
	for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
	{
		Event * e = new Event();
		e->index		= i;
		e->propensity	= 0.0;
		eventVector.push_back(e);
	}
	
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");
	
	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations		= 0;
		timePoint				= 0;
		zeroData();
		int numberOfNegative	= 0;
		isNegative = false;

		simulation->loadInitialConditions();
		
		while (t < tEnd)
		{
			saveData();
			
			computePropensities();
			a0 = blitz::sum(propensitiesVector);
			typePrevios = type;
			
			// build list of critical reactions
			critical = listOfCriticalReactions();

			// sort the list
			if (numberOfIterations % simulation->SortInterval == 0)
				sort(eventVector.begin(), eventVector.end(), EventSort());
			
			
			if (isNegative == false)
			{
				tau1 = computeAdaptiveTimeStep(critical,type, tau_expl);
				if (tau1 > 2147483647)
				{
					t = tEnd;  // stoping criteria
				}
			}
			
			// if tau is less than SSA execute few SSA steps
			if (tau1  < (1.0/a0) * sgamma( (double)1.0 ) )
			{
				execute_SSA(type,t,numberOfIterations);
			}
			else 
			{
				tau2 = time_of_next_critical_reaction(critical);
				tau = tauL_sampling(tau1,tau2,tau_expl,type,critical);
				
				if ( isProposedNegative() == false)
				{
					acceptNewSpeciesValues();
					++numberOfIterations;
					t += tau;
					isNegative = false;
				}
				else 
				{
					////cout << " Negative species at time: "<< t << endl;
					tau1 = tau1 * 0.5;
					reloadProposedSpeciesValues();
					isNegative = true;
					numberOfNegative ++;
				}
			}

		}
		
		saveData();
		cout << "Sample: " << samples << endl;
		cout << " Number of Negative species " << numberOfNegative << endl;
		cout << "Number of realizations: " << numberOfIterations << endl;
		writeToAuxiliaryStream( simulation->speciesValues );
		averNumberOfRealizations += numberOfIterations;
		averNumberOfNegative += numberOfNegative;
	}
	
	writeData(outputFileName); 
	closeAuxiliaryStream();
	
	for (int i = 0; i < eventVector.size(); ++i) { delete eventVector[i]; }
	
	cout << " Average number of Realizations in Adaptive S-leaping CL:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;
	
	cout << " Average number of Negative Species in Adaptive S-leaping CL:" << endl;
	cout << averNumberOfNegative/numberOfSamples << endl;
}



















