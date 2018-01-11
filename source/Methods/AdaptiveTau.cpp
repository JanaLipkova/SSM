/*
 *  AdaptiveTau.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 5/2/13.
 *  Copyright 2013 CSElab . All rights reserved.
 *
 */


//****************************
// IDEA OF ADAPTIVE TAU ALGORITHM:
//****************************

// 1. compute propensities
// 2. build list of critical reactions
// 3. compute tau1 with PEC and type of stifness
// 4. compare with SSA and solve SSA
// 5. tau2 = time of next critical reactions
// 6. choose tau from tau1 and tau2, based on type of stiffness do sampling k_j
// 7. update state and time, and check for negative population
//****************************


#include "AdaptiveTau.h"
#include "RootFinderJacobian.h"




AdaptiveTau::AdaptiveTau(Simulation * simulation):
LeapMethod(simulation)
{ }

AdaptiveTau::~AdaptiveTau()
{ }

//****************************
//  build list of critical reactions
//****************************
vector<int> AdaptiveTau::listOfCriticalReactions()
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
double AdaptiveTau::computeTimeStep(vector<int> criticalReactions, int& type, int& crit)
{
	double tauPrimeExp;			// explicit tau
	double tauPrimeImp;         // implicit tau
	double tau1;				// tau for non-critical reactions
	double tau2;				// time of 1st critical reaction
	double epsilon	= simulation->Epsilon;
        double delta = 0.05; //simulation->Delta;

	vector<int> rev;  // vector of reversible reactions

#ifdef Dimerization
	rev.push_back(1);
	rev.push_back(2);
#elif defined(LacZLacY)
    rev.push_back(0 );
    rev.push_back(1 );
    rev.push_back(7 );
    rev.push_back(8 );
    rev.push_back(9 );
    rev.push_back(10);
#endif

	list<int>    non_critical;    // list of reactions used to compute muHat, signaHat2
	vector<int>  in_pec;          // list of reactions in PEC
	int Nstiff = 100;			  // recomended in Cao at all, "The Adaptive Explicit-Implict Tau-Leaping Method with Automatic tau selection"

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

	if (criticalReactions.empty())
		LeapMethod::computeMuHatSigmaHat2(muHat, sigmaHat2);
	else
	{
	    int j = 0;
		for (int ir = 0; ir < numberOfReactions; ++ir)
		{
			if (ir != criticalReactions[j] )   // non critical reactions
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

	/* 2. STEP IMPLICIT TAU */
	//------------------------
	// 1. find reactions in PEC
	// 2. build list of reactions considered in muHat,sigmaHat2 computation
	// 3. compute implicit tau

	double temp1;  // to store propenisty of forward reactions
	double temp2;  // to store propensity of backward reactions

	if (rev.empty())
		tauPrimeImp = tauPrimeExp;
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
			tauPrimeImp = tauPrimeExp;
		else
		{
			if (criticalReactions.empty())
			{
				int k =0;
				for (int ir =0; ir < numberOfReactions ; ++ir)
				{
					if (ir != in_pec[k])
						non_critical.push_back(ir);
					else { k++; }
				}
			}
			else
			{
				for (vector<int>::iterator it = in_pec.begin(); it !=in_pec.end(); ++it)
					non_critical.remove(*it);
			}

			muHat = 0.0; sigmaHat2 = 0.0;
			computeMuHatSigmaHat2(muHat, sigmaHat2,non_critical);

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

			tauPrimeImp = tau;
		}
	}


	/* 3. TAU AND STIFFNES */
	if (tauPrimeImp > Nstiff * tauPrimeExp)
	{
		type = 1;       // Stiff system take implicit tau
		tau1 = tauPrimeImp;
	}
	else
	{
		type = 0;      // non-stiff system and explicit tau
		tau1 = tauPrimeExp;
	}


	/* 3. TIME OF CRITICAL REACTION */
	//-------------------------------
	if (criticalReactions.empty())
	{ tau2 = HUGE_VAL;}
	else
	{
		double ac0 = 0.0;

		for (vector<int>::iterator it = criticalReactions.begin(); it!=criticalReactions.end() ; ++it)
			ac0 += propensitiesVector(*it);

		tau2 = (1.0/ac0) * sgamma( (double)1.0 );
	}


	/* 4. Choose final tau and number of critical reactions to be fired */
	//----------------------------------------------------------------------
	if (tau1 < tau2)
	{
		tau = tau1;
		crit = 0;
	}
	else
	{
		tau = tau2;
		crit = 1;

		if ((type ==0)||(type==1 && tau2 < tauPrimeExp) )
		{ type = 0; }
		else{ type = 1; }
	}

	return tau;

}
//****************************











//****************************
//  muHat, sigmaHat2 computation:
//****************************
// this method is overloaded from the LeapMethod class, since in adaptive tau leaping we exclude
// critical reactions from computation of explicit tau and additionaly we exclude reactions in PEC
// in computation of implicit tau
void AdaptiveTau::computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2, list<int> non_critical)
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
void AdaptiveTau::executeSSA(int& type, double& t, int& numberOfIterations)
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
	{using namespace RootFinderJacobian;
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
		computePropensities(propensitiesVector, 0);
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
			t = HUGE_VAL;;
		}
	}

	t += time;
}
//****************************













//****************************
//    Sampling
//****************************
void AdaptiveTau::sampling(double tau, int type,vector<int> criticalReactions, int crit)
{

	int numberOfReactions	= sbmlModel->getNumReactions();
	vector<long int> fire(numberOfReactions,0);   // to store temporal state of fireReactionProposed
	double aj;

	if (type == 0)  // explicit sampling from non-critical reactions
	{
		if ( criticalReactions.empty() ) {
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j){
				aj = propensitiesVector(j);
                fire[j] = (aj == 0) ? 0. : ignpoi( aj*tau );
			}
		}
		else
		{
			int k = 0;
			for (int j = 0; j < propensitiesVector.extent(firstDim); ++j){
				if ( j!=criticalReactions[k] ){
					aj	= propensitiesVector(j);
                    fire[j] = (aj == 0) ? 0. : ignpoi( aj*tau );
				}
				else{
					k++;
					fire[j]=0;
				}
			}
		}
	}
	else
		implicit_sampling(tau,criticalReactions,fire);


	// add one critical reaction
	if (crit == 1)
	{
		double ac0 = 0.0;   // sum of propensities of critical reactions

		for (vector<int>::iterator it = criticalReactions.begin(); it!=criticalReactions.end() ; ++it)
		 ac0 += propensitiesVector(*it);

		double r1 = ranf();
		int reactionIndex = 0;
		double cummulative = 0.0;

		for (vector<int>::iterator it = criticalReactions.begin(); it != criticalReactions.end(); ++it)
		{
			cummulative += propensitiesVector(*it);
			if ( cummulative > ac0*r1 )
			{
				reactionIndex = *it;
				break;
			}
		}

		fire[reactionIndex]=1;
	}

	// propose reaction to be fired
	for (int j=0; j<numberOfReactions; j++)
	{
		if (fire[j] != 0)
			fireReactionProposed( j , fire[j] );

		if (fire[j]<0)
		{
			cout << "Abort! Proposed k[j] is negative "<<endl;
			abort();
		}
	}
}
//****************************







//****************************
//      IMPLICI SAMPLING
//****************************
//  1. solve:
//      X(t+tau) = X(t) + sum_j vj*aj(X(t+tau))*tau + sum_j vj*[Poiss(aj(X(t)*tau) - aj(X(t))*tau]
//  2. denote solution of above system Xhat and find implicit sampling as:
//     kj = round[aj(Xhat) + Poiss(aj(X(t))*tau) - aj(X(t))*tau]

// notation used for simplification:
//          vector of pairs, k = where kj = Poiss(aj(X(t))*tau);  // length = numberOfReactions
//			vector aTau, aTau[j] = aj(X(t))*tau        // length = numberOfReactions
//          vector  B, B = X(t) + sum_j  vj*[Poiss(aj(X(t)*tau) - aj(X(t))*tau]
//          ... where rhs stands for all term that are fixed number in the rootfinder procedure, and thus can pre-computed

void AdaptiveTau::implicit_sampling( double tau, vector<int> critical, vector<long int>& fire)
{

	using namespace RootFinderJacobian;

    double aj;
	long int kj;
	int MaxNumberOfIterations = 100;

	int numberOfSpecies		= sbmlModel->getNumSpecies();
	int numberOfReactions	= sbmlModel->getNumReactions();

	vector<long int> k(numberOfReactions);  // sampled reactions,on j-th position is how many times j-th reaction should be fired
	vector<double> aTau(numberOfReactions);
	vector<double> B(numberOfSpecies,0);
	vector<double> implicitPropenisty(numberOfReactions,0);   // propenisties in the roots of implicit system
	vector<double> roots(numberOfSpecies,0);                  // to store roots of implict system, denote X^ in literature


	//STEP 1: precompute k, aTau, B
	for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
    {
		aj = propensitiesVector(j);
        k[j]  = (aj == 0) ? 0. : ignpoi( aj*tau );
		aTau[j] = aj*tau;
	}


	for (int i = 0; i < numberOfSpecies; i++)
		B[i] = simulation->speciesValues(i);


	for (int j = 0; j < numberOfReactions; j++)
	{
		SSMReaction * ri = simulation->ssmReactionList[j];
		const vector<int> & changes = ri->getChanges();
		const vector<int> & nuChanges = ri->getNuChanges();

		for (int s = 0; s < changes.size(); s++)
			B[changes[s]] +=nuChanges[s] * (k[j] - aTau[j]);
	}

	//STEP 2: solve implicit system
	using namespace RootFinderJacobian;
	RootFinderSetUp(simulation, tau, MaxNumberOfIterations, B );
	find_roots(roots, implicitPropenisty, simulation->speciesValues);


	// STEP 3: sample reaction channels from implicit stat: k_j(X(t+tau))
	if (critical.empty())
	{
		for (int j=0; j < numberOfReactions; j++)
		{
			kj = round(implicitPropenisty[j]*tau + k[j] - aTau[j]);
			fire[j]=kj;

			if (kj<0)
			{
				cout << "Abort in AdaptiveTau::implicit_sampling, in non critical section. Negative kj!!!"<<endl;
				//abort();
				fire[j] = 0;
			}
		}
	}
	else
	{
		int m = 0;
		for (int j=0; j<numberOfReactions; j++)
		{
			if(j!= critical[m])
			{
				kj = round(implicitPropenisty[j]*tau + k[j] - aTau[j]);
				fire[j]=kj;


				if (kj<0)
				{
					cout << "Abort in AdaptiveTau::implicit_sampling in critical section. Negative kj!!!"<<endl;
					abort();
				}
			}
			else
			{
				m++;
				fire[j] = 0;
			}
		}
	}

}
//****************************


void AdaptiveTau::executeSSA_lacZlacY(double& t, int SSAsteps, double genTime)
{
    int count = 0.;
    double a0 = 0.;
    double tau;
    double r1;
    int reactionIndex = 0;
    double cummulative = 0.0;

    while (count < SSAsteps)
    {
        computePropensitiesGrowingVolume(propensitiesVector,t,genTime);
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

        count++;

        // RNAP     = S(1) ~ N(35),3.5^2)
        // Ribosome = S(9) ~ N(350,35^2)
           simulation->speciesValues(1)  = gennor(35   * (1 + t/genTime), 3.5);
           simulation->speciesValues(9)  = gennor(350  * (1 + t/genTime),  35);
    }

}



//****************************
//      SOLVE
//****************************
void AdaptiveTau::solve()
{
	cout << "Adaptive Tau Leaping..." << endl;
    openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");

	double a0 = 0.0;
	double tau;
	int type = 0;      // type = 0 for stiff system, type = 1 for non-stiff one
	int typePrevios   = 0;
	bool isNegative = false;
	double averNumberOfRealizations = 0.0;
    	vector<int> rejectionsVector(numberOfSamples);
   	int numberOfRejections;
    	const int SSAfactor = 10;
    	const int SSAsteps = 100;
	double genTime = 2100;

	for (int samples = 0; samples < numberOfSamples; ++samples)
	{
		t = simulation->StartTime;
		numberOfIterations = 0;
		numberOfRejections = 0;
		timePoint = 0;
		simulation->loadInitialConditions();
        
        while (t < tEnd)
        {
            saveData();

#ifdef LacZLacY
            // RNAP     = S(1) ~ N(35),3.5^2)
            // Ribosome = S(9) ~ N(350,35^2)
            simulation->speciesValues(1)  = gennor(35   * (1 + t/genTime), 3.5);
            simulation->speciesValues(9)  = gennor(350  * (1 + t/genTime),  35);
            computePropensitiesGrowingVolume(propensitiesVector,t,genTime);
#else
            computePropensities(propensitiesVector, 0);
#endif
            a0 = blitz::sum(propensitiesVector);
            typePrevios = type;
            int crit = 0;  // number of critical reactions to be fired
            
            vector<int> criticalReactions = listOfCriticalReactions();
            
            if (isNegative == false){
                tau = computeTimeStep(criticalReactions,type,crit);
                if (tau > HUGE_VAL){t = tEnd; break;}  // stoping criteria
            }
            
            //if( dt <= SSAfactor * (1.0/a0) * sgamma( (double)1.0 ) )
	    //{
            //                 #ifdef LacZLacY
            //      executeSSA_lacZlacY(t, SSAsteps, genTime);
            //      #else
            // 		  executeSSA(type, t, numberOfIterations);
            //      #endif
	   // }
	   // else
           // {
                sampling(tau,type,criticalReactions,crit);
                
                if (isProposedNegative() == false){
                    acceptNewSpeciesValues();
                    ++numberOfIterations;
                    t += tau;
                    isNegative = false;
                }
                else{
                    tau = tau * 0.5;
                    reloadProposedSpeciesValues();
                    isNegative = true;
                    ++numberOfRejections;
                }
            //}
	}

        cout << "Sample: " << samples << endl;
        saveData();
        rejectionsVector[samples] = numberOfRejections;
	writeToAuxiliaryStream( simulation->speciesValues );
	averNumberOfRealizations += numberOfIterations;
       }

	writeData(outputFileName);
	closeAuxiliaryStream();

	cout << " Average number of Realizations in Adaptive Tau-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;

    int rejectionSum = std::accumulate(rejectionsVector.begin(), rejectionsVector.end(), 0);
    std::cout<<"Negative species appeared in total:" << rejectionSum << " times" << std::endl;
}
