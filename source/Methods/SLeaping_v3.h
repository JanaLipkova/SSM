/*
 *  SLeaping_v3.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 4/25/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#pragma once
#include "LeapMethod.h"

class SLeaping_v3: public LeapMethod
{
public:
    SLeaping_v3(Simulation * simulation);
    ~SLeaping_v3();

    // override the virtual method
    void solve();

private:
    double    computeTimeStep();
    void      computePropensities();   	// override the standard calculation of propensities
    void      computePropensitiesGrowingVolume(Array< double , 1 > & propensitiesVector, double time, double genTime);
    void      sampling(double& dt, double a0, long int L);
    //void      sampling(double& dt, double a0);
    void      executeSSA(double& t, int SSAsteps);
    void      executeSSA_lacZlacY(double& t, int SSAsteps, double genTime);

    // anonymous inner class, S-Leaping( same as R) needs to store the indices and propensities of reactions
    class Event
    {
    public:
        Event() {}
        ~Event() {}

        double	propensity;
        int		index;
    };

    // sorts the global reactions with respect to the largest propensities (i.e. descending order
    // of propensities - S-Leaping requirement
    class EventSort
    {
    public:
        bool operator() ( Event * a, Event * b)
        {
            return (a->propensity > b->propensity);
        }
    };

    vector<Event*> eventVector;

};
