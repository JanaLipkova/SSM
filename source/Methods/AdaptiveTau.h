/*
 *  AdaptiveTau.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 5/2/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"
#include "RootFinderJacobian.h"

class AdaptiveTau : public LeapMethod 
{
public:
	AdaptiveTau(Simulation * simulation);
	~AdaptiveTau();
	
	// override the virtual method
	void solve();
	
private:
	
	vector<int> listOfCriticalReactions();

	double computeTimeStep(vector<int> criticalReactions, int& type, int& crit);

	// overwrite standard method for MuHat and SigmaHat2 computation,
	// now both parameters are computed only wrt given list of non-critical reactions
	// this is needed because of list of critical reactions 
	// and for the PEC in implicit tau step computation
	void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::list<int> non_critical);

	void execute_SSA(int& type, double& t, int& numberOfIterations);
	void sampling(double tau, int type, vector<int> criticalReactions, int crit);
	void implicit_sampling( double tau,vector<int> critical, vector<long int>& fire);
	
};

