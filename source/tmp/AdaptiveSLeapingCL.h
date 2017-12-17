/*
 *  AdaptiveSLeapingCL.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 7/3/13.
 *  Copyright 2013 CSE Lab. All rights reserved.
 *
 */


#pragma once

#include "LeapMethod.h"
#include "RootFinderJacobian.h"

class AdaptiveSLeapingCL : public LeapMethod 
{
public:
	AdaptiveSLeapingCL(Simulation * simulation);
	~AdaptiveSLeapingCL();
	
	// override the virtual method
	void solve();
	
	// anonymous inner class, R-Leaping needs to store the indices and propensities of reactions
	class Event
	{
	public:
		Event() {}
		~Event() {}
		
		double	propensity;
		int		index;
	};
	
	
	class EventSort 
	{
	public:
		bool operator() ( Event * a, Event * b) 
		{
			return (a->propensity > b->propensity);
		}
	};
	
	vector<Event*> eventVector;
	
private:
	
	vector<int> listOfCriticalReactions();	
	double computeAdaptiveTimeStep(vector<int> critical, int& type, double& tau_expl);
	
	// overwrite standard method for MuHat and SigmaHat2 computation,
	// now both parameters are computed only wrt given list of non-critical reactions
	// this is needed for implicit tau step computation
	void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::list<int> non_critical);
	
    void execute_SSA(int& type, double& t, int& numberOfIterations);
	double time_of_next_critical_reaction(vector<int> critical);
	
	double tauL_sampling(double tau1, double tau2, double tau_expl, int& type,vector<int> critical);
	void sampling(long int L,double a0);
	void sampling_critical(long int L,vector<int> critical,double a0nc, int crit_count,double ac0);	
	void implicit_sampling(double tau, vector<int> critical);

	bool iscritical(int j, vector<int> critical);

	// override the standard calculation of propensities
	void computePropensities();
	
	
};






