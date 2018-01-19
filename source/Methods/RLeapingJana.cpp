/*
 *  RLeapingJana.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Lipkova on 4/24/13.
 *  Copyright 2013 CSElab. All rights reserved.
 *
 */


/* Improved RLeaping.cpp:
   - max #sampling in original RLeaping is M, not M-1 !!!
   - this is fixed in this version
 */

#include "RLeapingJana.h"

RLeapingJana::RLeapingJana(Simulation * simulation):
LeapMethod(simulation)
{ }

RLeapingJana::~RLeapingJana()
{ }

long int RLeapingJana::computeLeapLength()
{
	long int L				 = 2147483647;//2^{31}-1

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
				tau = min(tau,epsixi/fabs(muHat(is)));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/varHat(is));
				break;
			default:
				break;

		}
	}

	L = (long int)max((long int)(tau*a0), (long int)1);

	//   used to speed up simulation for theta != 0
/*    if(theta > 0 ){
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
                long int propsedL = (long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj);
                if (propsedL < L && propsedL > 0)
                    L = propsedL;

                L = min( (long int)L, (long int)((1.0-theta*(1.0-a0/propensitiesVector(ir)))*lj) );
            }
        }
    }
*/
	return L;
}

// this method is overloaded from the Methods class since R-Leaping needs to store both indices and
// propensities (not just propensities).  These are located in the anonymous inner class called Event
void RLeapingJana::computePropensities()
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


void RLeapingJana::computePropensitiesGrowingVolume(Array< double , 1 > & propensitiesVector, double time, double genTime)
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




//***********************************
//       Sampling
// Binomial sampling from reordered propensities
// no critical reactions considered
//***********************************
void RLeapingJana::sampling(long int L, double a0)
{
  	double p = 0.0;
    	double cummulative	= a0;
    	long int k			= 0;

    	for (int j = 0; j < eventVector.size(); ++j){
        	if( (j == eventVector.size() - 1 ) && (L != 0) ) // last reaction to be fires
        	{
            		fireReactionProposed( eventVector[j]->index , L);
            		break;
        	}	

        	cummulative -= p;
        	p = eventVector[j]->propensity;

        	if(p!=0)
        	{
            		k = ignbin(L, min(p/cummulative, 1.0) );
            		L -= k;

            		fireReactionProposed( eventVector[j]->index , k);
            		if (L == 0){ break; }
        	}
	}

}


void RLeapingJana::executeSSA(double& t, int SSAsteps)
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

void RLeapingJana::executeSSA_lacZlacY(double& t, int SSAsteps, double genTime)
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




//***********************************

void RLeapingJana::solve()
{
	cout << "RLeaping..." << endl;
	openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");

	double a0						= 0.0;
	long int Lcurrent				= 1;
	bool isNegative					= false;
	double averNumberOfRealizations = 0.0;
    vector<int> rejectionsVector(numberOfSamples);
	int numberOfRejections          = 0;
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
		numberOfIterations	= 0;
		numberOfRejections	= 0.;
		timePoint			= 0;
        whenToSave          = t;
		zeroData();
		simulation->loadInitialConditions();
		Lcurrent			= 1;
		isNegative			= false;
        
        #ifdef DEBUG_PRINT
            tempArray.resize(sbmlModel->getNumSpecies());
            myfile.open ("all-times.txt");
        
            myfile << t << "\t";
            tempArray = simulation->speciesValues(Range::all());
            for (int i = 0; i < tempArray.extent(firstDim); ++i){
                myfile << tempArray(i) << "\t";
            }
            myfile << endl;
        #endif
        
        saveData();

        while (t < tEnd)
        {
            #ifdef LacZLacY
                // RNAP     = S(1) ~ N(35),3.5^2)
                // Ribosome = S(9) ~ N(350,35^2)
                simulation->speciesValues(1)  = gennor(35   * (1 + t/genTime),  3.5);
                simulation->speciesValues(9)  = gennor(350  * (1 + t/genTime),   35);
                computePropensitiesGrowingVolume(propensitiesVector,t,genTime);
                #else
                computePropensities();
            #endif
            
            a0 = blitz::sum(propensitiesVector);
            
            if (numberOfIterations % simulation->SortInterval == 0)
                sort(eventVector.begin(), eventVector.end(), EventSort());
            
            if (isNegative == false)
                Lcurrent = computeLeapLength();

            sampling(Lcurrent, a0);
            
            if (isProposedNegative() == false)
            {
                acceptNewSpeciesValues();
                ++numberOfIterations;
                dt = (1.0/a0) * sgamma( (double)Lcurrent ); // Gamma ( L, 1.0 / a0 )
                t_old = t;
                t += dt;
                isNegative = false;
                saveData();
                
                #ifdef DEBUG_PRINT
                    myfile << min(t,tEnd) << "\t";
                    if(t<tEnd)
                        tempArray =  simulation->speciesValues(Range::all());
                    else
                        tempArray =  simulation->old_speciesValues(Range::all());
                
                    for (int i = 0; i < tempArray.extent(firstDim); ++i){
                        myfile << tempArray(i) << "\t";
                    }
                    myfile << endl;
                #endif
            }
            else
            {
                ++numberOfRejections;
                Lcurrent = max( (int)(Lcurrent*0.5), 1);
                reloadProposedSpeciesValues();
                isNegative = true;
            }
        }
        
        cout << "Sample: " << samples << endl;
        rejectionsVector[samples] = numberOfRejections;
        writeToAuxiliaryStream( simulation->speciesValues );
        averNumberOfRealizations += numberOfIterations;
        
        #ifdef DEBUG_PRINT
            myfile.close();
        #endif
    }

	writeData(outputFileName);
	closeAuxiliaryStream();
    
	cout << " Average number of Realizations in R-leaping:" << endl;
	cout << averNumberOfRealizations/numberOfSamples << endl;

    int rejectionSum = std::accumulate(rejectionsVector.begin(), rejectionsVector.end(), 0);
    std::cout<<"Negative species appeared in total:" << rejectionSum << " times" << std::endl;

	for (int i = 0; i < eventVector.size(); ++i) { delete eventVector[i]; }

}
