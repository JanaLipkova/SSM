/*
 *  AdaptiveSLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 5/17/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#include "AdaptiveSLeaping.h"
#include "RootFinderJacobian.h"

// class constructor
AdaptiveSLeaping::AdaptiveSLeaping(Simulation* simulation):
LeapMethod(simulation)
{ }

// destructor
AdaptiveSLeaping::~AdaptiveSLeaping()
{ }


//****************************
//  compute tau1 with PEC and type of stifness
//****************************
double AdaptiveSLeaping::computeAdaptiveTimeStep(int& type)
{
    double tauPrimeExp;			// explicit tau
    double tauPrimeImp;         // implicit tau
    double epsilon	= simulation->Epsilon;
    double delta = 0.05; //simulation->Delta;
    
    //vector<int> reversible = simulation -> Reversible;  // indecies of reverisble reactions
    vector<int> rev;
    
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
    muHat = 0.0; sigmaHat2 = 0.0;
    LeapMethod::computeMuHatSigmaHat2(muHat, sigmaHat2);
    
    double tau, taup,  epsi, epsixi, epsixisq;
    double xi;
    tau = HUGE_VAL;
    
    double a0 = (double)blitz::sum(propensitiesVector);
    for (int is = 0; is < numberOfSpecies; is++)
    {varHat(is) = sigmaHat2(is) - (1.0/a0) * muHat(is) * muHat(is);}
    
    
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
    
    tauPrimeExp = tau;    // explicit tau
    
    
    
    /* 2. STEP IMPLICIT TAU */
    // 1. find reactions in PEC
    // 2. build list of reactions considered in muHat,sigmaHat2 computation
    // 3. compute implicit tau
    
    double temp1;  // to store propenisty of forward reactions
    double temp2;  // to store propensity of backward reactions
    
    if (rev.empty())
    { tauPrimeImp = tauPrimeExp;}
    else
    {
        for (vector<int>::iterator it = rev.begin(); it !=rev.end(); it=it+2)
        {
            temp1 = propensitiesVector(*it);    // propensity of forward reaction
            temp2 = propensitiesVector(*it+1);  // propensity of backward reacrion
            //cout <<"temp1="<<temp1<<" temp2 = "<<temp2<<endl;
            if ( abs( temp1 - temp2 ) < delta* abs( temp1 + temp2 ) )
            {
                in_pec.push_back(*it);
                in_pec.push_back(*it+1);
            }
        }
        
        if (in_pec.empty())
        { tauPrimeImp = tauPrimeExp; }
        else
        {
            int k =0;
            for (int ir =0; ir < numberOfReactions ; ++ir)
            {
                if (ir != in_pec[k])
                { non_critical.push_back(ir); }
                else { k++; }
            }
            
            muHat = 0.0; sigmaHat2 = 0.0;
            computeMuHatSigmaHat2(muHat, sigmaHat2,non_critical);
            
            double tau, taup,  epsi, epsixi, epsixisq;
            double xi;
            
            tau = HUGE_VAL;
            
            //			a0 = 0.0;
            //			for (list<int>::iterator it = non_critical.begin(); it != non_critical.end(); ++it)
            //				a0 += propensitiesVector(*it);
            //
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
    if (tauPrimeImp > Nstiff * tauPrimeExp){
        type = 1;       // Stiff system take implicit tau
        tau = tauPrimeImp;
    }
    else{
        type = 0;      // non-stiff system and explicit tau
        tau = tauPrimeExp;
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
void AdaptiveSLeaping::computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2, list<int> non_critical)
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
void AdaptiveSLeaping:: execute_SSA(int& type, double& t, int& numberOfIterations)
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
        //computePropensities(propensitiesVector,0);
        computePropensities();
        a0 = blitz::sum(propensitiesVector);
        
        dt = (1.0/a0) * sgamma( (double)1.0 );
        
        r1 = ranf();
        reactionIndex = -1;
        cummulative = 0.0;
        
        
        for(int ev = 0; ev < eventVector.size(); ++ev)
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
            ++numberOfIterations;
            time += dt;
        }
        else { count = steps; }
    }
    
    t += time;
}
//****************************

void AdaptiveSLeaping::sampling( double& tau, int type, double a0, vector<AdaptiveSLeaping::Event *>& eventVector)
{
    if (type == 0)
        explicit_sampling(tau, a0, eventVector);
    else
        implicit_sampling(tau, a0, eventVector);
    
}


void AdaptiveSLeaping::explicit_sampling(double& tau, double a0, vector<AdaptiveSLeaping::Event *>& eventVector)
{
  
    long int L =  (long int)max( (long int)ignpoi(a0*tau), (long int)1);
  
    long int Llocal = L;
    double r1;
    double suma = 0.0;
    
    double p = 0.0;
    double cummulative	= a0;
    long int k			= 0;
    
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
    
   tau =  (1.0/a0) * sgamma( (double)L );
}


void AdaptiveSLeaping::implicit_sampling(double& tau, double a0, vector<AdaptiveSLeaping::Event *>& eventVector)
{
    long int L;
    int ir;
    double aj;
    long int kj;
    int MaxNumberOfIterations = 100;
    
    int numberOfSpecies		= sbmlModel->getNumSpecies();
    int numberOfReactions	= sbmlModel->getNumReactions();
    
    vector<int>    k(numberOfReactions,0);        // sampled reactions,on j-th position is how many times j-th reaction should be fired
    vector<double> B(numberOfSpecies,0);         // rhs of the implicit system
    vector<double> tmp(numberOfReactions,0);      // aj * L / a0 term
    vector<double> implicitPropenisty(numberOfReactions,0);   // propenisties in the roots of implicit system
    vector<double> roots(numberOfSpecies,0);                  // to store roots of implict system, denote X^ in literature
    
    //	// STEP 1: precompute k, B
    long int Lexp = (long int)max( (long int)ignpoi(a0*tau), (long int)1);
    long int Llocal = Lexp;
    double r1;
    double suma = 0.0;
    
    double p = 0.0;
    double cummulative	= a0;
    
    for (int ev = 0; ev < eventVector.size(); ++ev)
    {
        int ir = eventVector[ev]->index;
        cummulative		-= p;
        p				 = eventVector[ev]->propensity;
        k[ir]			 = ignbin(Llocal, min(p/cummulative, 1.0) );
        Llocal			-= k[ir];
        
        if (Llocal == 0){ break; }
    }
    
    
    for (int i = 0; i < numberOfSpecies; i++)
        B[i] = simulation->speciesValues(i);
    
    for (int j = 0; j < numberOfReactions; j++)
    {
        SSMReaction * ri = simulation->ssmReactionList[j];
        const vector<int> & changes = ri->getChanges();
        const vector<int> & nuChanges = ri->getNuChanges();
        
        tmp[j] = propensitiesVector(j) * Lexp / a0;
        
        for (int s = 0; s < changes.size(); s++)
            B[changes[s]] +=nuChanges[s] * (k[j] - tmp[j]);
    }
    
    //STEP 2: solve implicit system
    RootFinderJacobian::RootFinderSetUp(simulation, tau, MaxNumberOfIterations, B );
    RootFinderJacobian::find_roots(roots, implicitPropenisty, simulation->speciesValues);
    
    // STEP 3: use implicitPropensities to propose reactions to be fired
    int L_imp = 0;
    double a0_imp = 0.;
    for (int j=0; j < numberOfReactions; j++)
    {
        int kj = round(implicitPropenisty[j]*tau + k[j] - tmp[j]);
        fireReactionProposed( j , kj);
        L_imp += kj;
        a0_imp += implicitPropenisty[j];
    }

   //cout<< "old L =" << Lexp << "newL="<<L_imp<<endl;
   double tau_old = tau;
   tau = (1.0/a0_imp) * sgamma( (double)L_imp ); 
   //cout<< "old tau =" << tau_old << "new tau="<<tau<<endl;
}

////********************************************
////   Compute L = number of all reactions to be fired
////********************************************
//long int AdaptiveSLeaping::computeLeapLength(double& tau,int& type, vector<AdaptiveSLeaping::Event *>& eventVector)
//{
//
//	long int L;   // final L
//	double theta	= simulation->Theta;
//
//	if (type==0)
//	{
//		double a0 = (double)blitz::sum(propensitiesVector);
//		L =  (long int)max( (long int)ignpoi(a0*tau), (long int)1);
//
//        // L with theta (control of negative species)
//        if (theta > 0.)
//        {
//            for (int ir = 0; ir < sbmlModel->getNumReactions(); ++ir)
//            {
//                long int lj = 2147483647; // MAXIMUM INTEGER
//                SSMReaction * ssmReaction = simulation->ssmReactionList[ir];
//                if (propensitiesVector(ir) > 0.0)
//                {
//                    const vector<int> & changes		= ssmReaction->getChanges();
//                    const vector<int> & nuChanges	= ssmReaction->getNuChanges();
//
//                    for (int is = 0; is < changes.size() ; ++is)
//                    {
//                        if (nuChanges[is] > 0) break;
//                        lj = min(lj, -simulation->speciesValues(changes[is])/ nuChanges[is] );
//                    }
//                    //long int propsedL = (long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj);
//                    if (lj < L && lj > 0)
//                    { L = lj; }
//                }
//            }
//        }
//	}
//	else  /* implicit L*/
//	{
//		L = compute_implicit_L(tau, eventVector, theta) ;
//		//double a0 = (double)blitz::sum(propensitiesVector);
//		//L =  (long int)max( (long int)(a0*tau), (long int)1);
//
//
//		vector<int> nopec(2,0);
//		nopec[0] = 0;
//		nopec[1] = 3;
//		vector<int> pec(2,0);
//		pec[0] = 1;
//		pec[1] = 2;
//
//		for (vector<int>::iterator ir=nopec.begin(); ir!=nopec.end(); ++ir)
//		{
//			long int lj = 2147483647.; // MAXIMUM INTEGER
//			SSMReaction * ssmReaction = simulation->ssmReactionList[*ir];
//			if (propensitiesVector(*ir) > 0.0)
//			{
//				const vector<int> & changes		= ssmReaction->getChanges();
//				const vector<int> & nuChanges	= ssmReaction->getNuChanges();
//
//				for (int is = 0; is < changes.size() ; ++is)
//				{
//					if (nuChanges[is] > 0) break;
//					lj = min(lj, (long int) (-simulation->speciesValues(changes[is])/ nuChanges[is]) );
//					//long int temp = -roots[changes[is]]/ nuChanges[is];
//					//cout << "temp = " << temp << endl;
//					//lj = min(lj, temp );
//
//				}
//				if (lj < L && lj > 0)
//				{ L = lj; }
//			}
//		}
//	}
//
//
//	/*---*/
////	double a0 = (double)blitz::sum(propensitiesVector);
////	L =  (long int)max( (long int)ignpoi(a0*tau), (long int)1);
////
////	vector<int> nopec(2,0);
////	nopec[0] = 0;
////	nopec[1] = 3;
////	vector<int> pec(2,0);
////	pec[0] = 1;
////	pec[1] = 2;
////
////	for (vector<int>::iterator ir=nopec.begin(); ir!=nopec.end(); ++ir)
////	{
////		long int lj = 2147483647.; // MAXIMUM INTEGER
////		SSMReaction * ssmReaction = simulation->ssmReactionList[*ir];
////		if (propensitiesVector(*ir) > 0.0)
////		{
////			const vector<int> & changes		= ssmReaction->getChanges();
////			const vector<int> & nuChanges	= ssmReaction->getNuChanges();
////
////			for (int is = 0; is < changes.size() ; ++is)
////			{
////				if (nuChanges[is] > 0) break;
////				lj = min(lj, (long int) (-simulation->speciesValues(changes[is])/ nuChanges[is]) );
////				//long int temp = -roots[changes[is]]/ nuChanges[is];
////				//cout << "temp = " << temp << endl;
////				//lj = min(lj, temp );
////
////			}
////			if (lj < L && lj > 0)
////			{ L = lj; }
////		}
////	}
//
////	long int lj = min( (long int) ((simulation->speciesValues(0)*0.5 + simulation->speciesValues(1)*1 )), L);
////	if (lj < L && lj >0) {
////		L=lj;}
//
//	double a0 = (double)blitz::sum(propensitiesVector);
//	tau = (1.0/a0) * sgamma( (double)L );
//
//	assert(L > 0);
//	return L;
//}
////****************************
//
//
////****************************
////      IMPLICI L
////****************************
//long int AdaptiveSLeaping::compute_implicit_L(double tau, vector<AdaptiveSLeaping::Event *>& eventVector, double theta)
//{
//	long int L;
//	int ir;
//	double aj;
//	long int kj;
//	int MaxNumberOfIterations = 100;
//
//	int numberOfSpecies		= sbmlModel->getNumSpecies();
//	int numberOfReactions	= sbmlModel->getNumReactions();
//
//	vector<int> k(numberOfReactions);  // sampled reactions,on j-th position is how many times j-th reaction should be fired
//	vector<double> aTau(numberOfReactions);
//	vector<double> B(numberOfSpecies,0);
//	vector<double> implicitPropenisty(numberOfReactions,0);   // propenisties in the roots of implicit system
//	vector<double> roots(numberOfSpecies,0);                  // to store roots of implict system, denote X^ in literature
//
//
//    //	// STEP 1: precompute k, aTau, B
//    for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
//    {
//        aj				= propensitiesVector(j);
//        if (aj == 0){
//            k[j] = 0;
//        } else {
//            k[j] = ignpoi( aj*tau );
//        }
//        aTau[j] = aj*tau;
//    }
//
//    for (int i = 0; i < numberOfSpecies; i++)
//        B[i] = simulation->speciesValues(i);
//
//    for (int j = 0; j < numberOfReactions; j++)
//    {
//        SSMReaction * ri = simulation->ssmReactionList[j];
//        const vector<int> & changes = ri->getChanges();
//        const vector<int> & nuChanges = ri->getNuChanges();
//
//        for (int s = 0; s < changes.size(); s++)
//            B[changes[s]] +=nuChanges[s] * (k[j] - aTau[j]);
//    }
//
//    //STEP 2: solve implicit system
//    RootFinderJacobian::RootFinderSetUp(simulation, tau, MaxNumberOfIterations, B );
//    RootFinderJacobian::find_roots(roots, implicitPropenisty, simulation->speciesValues);
//
//
//
//	// STEP 3: set propensities vector and eventVector to new implicit propensities need for impicit sampling
//	for (int ev = 0; ev < eventVector.size(); ++ev)
//	{
//		//eventVector[ev]->propensity = 0.0;
//		ir							= eventVector[ev]->index;
//		SSMReaction* reaction		= simulation->ssmReactionList[ir];
//
//		reaction->setPropensity(implicitPropenisty[ir]);
//		propensitiesVector(ir)		= implicitPropenisty[ir];//reaction->getPropensity();
//		eventVector[ev]->propensity = implicitPropenisty[ir];//reaction->getPropensity();
//	}
//
//	// STEP 4: compute implicit L
//	double a0 = (double)blitz::sum(propensitiesVector);
////	L =  (long int)max( (long int)ignpoi(a0*tau), (long int)1);
//	L =  (long int)max( (long int)(a0*tau), (long int)1);
//
//
//	// L with theta (control of negative species) with current state
//	vector<int> nopec(2,0);
//	nopec[0] = 0;
//	nopec[1] = 3;
//	vector<int> pec(2,0);
//	pec[0] = 1;
//	pec[1] = 2;
//
//	for (vector<int>::iterator ir=nopec.begin(); ir!=nopec.end(); ++ir)
//	{
//		long int lj = 2147483647.; // MAXIMUM INTEGER
//		SSMReaction * ssmReaction = simulation->ssmReactionList[*ir];
//		if (implicitPropenisty[*ir] > 0.0)
//		{
//			const vector<int> & changes		= ssmReaction->getChanges();
//			const vector<int> & nuChanges	= ssmReaction->getNuChanges();
//
//			for (int is = 0; is < changes.size() ; ++is)
//			{
//				if (nuChanges[is] > 0) break;
//				//lj = min(lj, (long int) (-simulation->speciesValues(changes[is])/ nuChanges[is]) );
//				long int temp = -roots[changes[is]]/ nuChanges[is];
//				//cout << "temp = " << temp << endl;
//				lj = min(lj, temp );
//
//			}
//			if (lj < L && lj > 0)
//			{ L = lj; }
//		}
//	}
//
//	// reactions in PEC
////	long int lj = min( (long int) ((roots[0]*0.5 + roots[1]*1 )), L);
////	if (lj < L && lj >0) {
////		L=lj;}
//
//
//
//
//
////	for (int ir = 0; ir < numberOfReactions; ++ir)
////	{
////		long int lj = 2147483647; // MAXIMUM INTEGER
////		SSMReaction * ssmReaction = simulation->ssmReactionList[ir];
////		if (propensitiesVector(ir) > 0.0)
////		{
////			const vector<int> & changes		= ssmReaction->getChanges();
////			const vector<int> & nuChanges	= ssmReaction->getNuChanges();
////
////			for (int is = 0; is < changes.size() ; ++is)
////			{
////				if (nuChanges[is] > 0) break;
////				lj = min(lj, -simulation->speciesValues(changes[is])/ nuChanges[is] );
////				//cout << "lj: " << lj << endl;
////				//cout << "other: " << -simulation->speciesValues(changes[is])/ nuChanges[is] << endl;
////
////			}
////			long int propsedL = lj;
////			if (propsedL < L && propsedL > 0)
////			{
////				L = propsedL;
////			}
////			//L = min( (long int)L, (long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj) );
////			//cout << "	proposed	L: " << propsedL << endl;
////			//cout << "	compared to L: " << L << endl;//(long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj) << endl;
////		}
////	}
//
////	for (int i = 0; i <3; i++) {
////		simulation->speciesValues(i) = roots[i];
////		simulation->proposedSpeciesValues(i) = roots[i];
////	}
//
////
////	cout << " L after control= " << L << endl;
////	cout << endl;
//
//	assert(L > 0);
//	return L;
//}
////****************************
//
//
//
////****************************
//// Binomial sampling from reordered propensities
////****************************
//void AdaptiveSLeaping::sampling(long int L, double a0)
//{
//    double Llocal = L;
//    double r1;
//    double suma = 0.0;
//
//    double p = 0.0;
//    double cummulative	= a0;
//    long int k			= 0;
//
//    int j = 0;
//    for (j = 0; j < eventVector.size(); ++j)
//    {
//        cummulative		-= p;
//        p				 = eventVector[j]->propensity;
//        k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
//        Llocal			-= k;
//
//        fireReactionProposed( eventVector[j]->index , k);
//
//        if (Llocal == 0){ break; }
//    }
//
//#ifdef BLABLA
//	double Llocal = Lcurrent;
//	double r1;
//	double suma = 0.0;
//	double a0 = (double)blitz::sum(propensitiesVector);
//
//
//
//	if ( Llocal > sbmlModel->getNumReactions() )
//	{
//		double p = 0.0;
//		double cummulative	= a0;
//		long int k			= 0;
//
//		int j = 0;
//		for (j = 0; j < eventVector.size(); ++j)
//		{
//			cummulative		-= p;
//			p				 = eventVector[j]->propensity;
//			k				 = ignbin(Llocal, min(p/cummulative, 1.0) );
//			Llocal			-= k;
//
//			fireReactionProposed( eventVector[j]->index , k);
//
//			if (Llocal == 0){ break; }
//		}
//	}
//	else {
//		for (int s = 0; s<Llocal; s++)
//		{
//			r1 = ranf();
//			int j = 0;
//			suma = propensitiesVector(j)/a0;
//
//			while (r1 > suma)
//			{
//				j++;
//				suma += propensitiesVector(j)/a0;
//			}
//
//			fireReactionProposed(j,1);
//		}
//	}
//
//#endif
//}
////****************************
//


//****************************
// this method is overloaded from the Methods class since R-Leaping needs to store both indices and
// propensities (not just propensities).  These are located in the anonymous inner class called Event
//****************************
void AdaptiveSLeaping::computePropensities()
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
        propensitiesVector(ir)		= reaction->getPropensity();
        eventVector[ev]->propensity	= reaction->getPropensity();
    }
}


void AdaptiveSLeaping::computePropensitiesGrowingVolume(Array< double , 1 > & propensitiesVector, double time, double genTime)
{

        double volume   = 1. + time/genTime;
        double ivolume  = 1./volume;

        int nu;
        ParticleType x;
        ParticleType num, denom;

        int ir;         // the reaction index

        propensitiesVector = 0.0;
        vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;
        Reaction * sbmlreaction; // added maagm

        for (int ev = 0; ev < eventVector.size(); ++ev) // event index
        {
                eventVector[ev]->propensity = 0.0;
                ir                           = eventVector[ev]->index;

                SSMReaction* reaction           = ssmReactionList[ir];
                vector <int>  reactants         = reaction->getReactants();
                vector <int>  nu_reactants      = reaction->getNuReactants();
                int order                       = reaction->getOrder();

                 reaction->setPropensity(reaction->getRate());

                 for (int s = 0; s < reactants.size(); ++s)
                 {
                       nu              = nu_reactants[s];
                       x               = simulation->speciesValues( reactants[s] );
                       num             = x;
                       denom   = nu;
                       while ((--nu)>0)
                       {
                            denom   *= nu;
                            num     *= (x - nu);
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


                propensitiesVector(ir)          = reaction->getPropensity();
                eventVector[ev]->propensity     = reaction->getPropensity();
        }
}
//****************************
//      SOLVE
//****************************
void AdaptiveSLeaping::solve()
{
    cout << "Adaptive S-Leaping..." << endl;
    openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");

    double a0 = 0.0;
    double tau;
    double tau_exp = 0.0;  // place holder for explicit time step, needed in the implicit solver
    int type = 0;      // type = 0 for stiff system, type = 1 for non-stiff one
    long int Llocal			= 1;
    long int Lcurrent		= 1;
    bool isNegative = false;
    double averNumberOfRealizations = 0.0;
    vector<int> rejectionsVector(numberOfSamples);
    int numberOfRejections;
    const int SSAfactor = 10;
    const int SSAsteps = 100;
    double genTime = 2100;
    
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
        numberOfRejections      = 0;
        timePoint				= 0;
        zeroData();
        isNegative = false;
        Llocal = 1;
        Lcurrent = 1;
        simulation->loadInitialConditions();
        
        while (t < tEnd)
        {
            saveData();
            computePropensities();
            a0 = blitz::sum(propensitiesVector);
            
            if (numberOfIterations % simulation->SortInterval == 0)
                sort(eventVector.begin(), eventVector.end(), EventSort());
            
            if (isNegative == false)
            {
                tau = computeAdaptiveTimeStep(type);
                if (tau > HUGE_VAL) {t = tEnd; break;}  // stoping criteria
            }
          
                sampling(tau, type,  a0, eventVector);

                if ( isProposedNegative() == false)
                {
                    acceptNewSpeciesValues();
                    ++numberOfIterations;
                    t += tau;
                    isNegative = false;
                }
                else
                {
                    tau = tau * 0.5;
                    reloadProposedSpeciesValues();
                    isNegative = true;
                    ++numberOfRejections;
                }
        }
        
        cout << "Sample: " << samples << endl;
        saveData();
        rejectionsVector[samples] = numberOfRejections;
        writeToAuxiliaryStream( simulation->speciesValues );
        averNumberOfRealizations += numberOfIterations;
    }
    
    writeData(outputFileName);
    closeAuxiliaryStream();
    
    cout << " Average number of Realizations in Adaptive S-leaping:" << endl;
    cout << averNumberOfRealizations/numberOfSamples << endl;
    
    int rejectionSum = std::accumulate(rejectionsVector.begin(), rejectionsVector.end(), 0);
    std::cout<<"Negative species appeared in total:" << rejectionSum << " times" << std::endl;
    
   for (int i = 0; i < eventVector.size(); ++i) { delete eventVector[i]; }
}
