/*
 *  TauLeaping.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "TauLeaping.h"

TauLeaping::TauLeaping(Simulation * simulation):
LeapMethod(simulation)
{ }

TauLeaping::~TauLeaping()
{ }

double TauLeaping::computeTimeStep()
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
    
    return tau;
}



//****************************
//     Execute SSA
//****************************
// if proposed time step is too small,execute numberOfIterations steps of SSA
void TauLeaping::executeSSA(double& t, int SSAsteps)
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
            
            t_old = t;
            t += tau;
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

void TauLeaping::executeSSA_lacZlacY(double& t, int SSAsteps, double genTime)
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
            
            t_old = t;
            t += tau;
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




void TauLeaping::solve()
{
    cout << "TauLeaping..." << endl;
    openAuxiliaryStream( (simulation->ModelName) + "_histogram.txt");

    double aj;
    long int kj;
    double a0						= 0.0;
    bool isNegative					= false;
    double averNumberOfRealizations = 0.0;
    vector<int> rejectionsVector(numberOfSamples);
    int numberOfRejections;
    const int SSAfactor             = 10;
    const int SSAsteps              = 100;
    double genTime                  = 2100;   // generation time

    // create C++11 rng
    std::default_random_engine engine;
    std::poisson_distribution<int> pois_dist(4.1);
 
    for (int samples = 0; samples < numberOfSamples; ++samples)
    {
        t                   = simulation->StartTime;
        whenToSave          = t;
        numberOfIterations = 0;
        numberOfRejections = 0.;
        timePoint = 0;
        zeroData();
        simulation->loadInitialConditions();
        isNegative = false;
        
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
                simulation->speciesValues(1)  = 35;//gennor(35   * (1 + t/genTime), 3.5);
                simulation->speciesValues(9)  = 350;//gennor(350  * (1 + t/genTime),  35);
                //computePropensitiesGrowingVolume(propensitiesVector,t,genTime);
            	computePropensities(propensitiesVector, 0);
	    #else
                computePropensities(propensitiesVector, 0);
            #endif
            
            a0 = blitz::sum(propensitiesVector);
            
            if (isNegative == false)
                dt = computeTimeStep();

		 
                for (int j = 0; j < propensitiesVector.extent(firstDim); ++j)
                {
                    aj = propensitiesVector(j);
                    pois_dist = std::poisson_distribution<int>(aj*dt);
                    kj =  pois_dist(engine);
		    //kj = (aj == 0) ? 0. : ignpoi( aj*dt );
                    fireReactionProposed( j , kj );
                }

                if (isProposedNegative() == false)
                {
                    acceptNewSpeciesValues();
                    ++numberOfIterations;
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
                  //  cout << "Negative species at time: " << t << endl;
                    ++numberOfRejections;
                    dt = dt * 0.5;
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
    
    cout << " Average number of Realizations in Tau-leaping:" << endl;
    cout << averNumberOfRealizations/numberOfSamples << endl;
    
    int rejectionSum = std::accumulate(rejectionsVector.begin(), rejectionsVector.end(), 0);
    std::cout<<"Negative species appeared in total:" << rejectionSum << " times" << std::endl;
}

