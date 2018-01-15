/*
 *  AdaptiveSLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 5/17/13.
 *  Copyright 2013 CSE Lab. All rights reserved.
 *
 */

#pragma once
#include "LeapMethod.h"


class AdaptiveSLeaping : public LeapMethod
{
public:
	AdaptiveSLeaping(Simulation* simulation);
	~AdaptiveSLeaping();

	//overrude the virtual method
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

	double computeAdaptiveTimeStep(int& type);

	// overwrite standard method for MuHat and SigmaHat2 computation,
	// now both parameters are computed only wrt to reactions out of PEC
	void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::list<int> non_critical);
	void execute_SSA(int& type, double& t, int& numberOfIterations);

    void sampling( double& tau, int type, double a0, vector<AdaptiveSLeaping::Event *>& eventVector);
    void explicit_sampling(double& tau, double a0, vector<AdaptiveSLeaping::Event *>& eventVector);
    void implicit_sampling(double& tau, double a0, vector<AdaptiveSLeaping::Event *>& eventVector);

//	long int computeLeapLength(double& tau,int& type, vector<AdaptiveSLeaping::Event *>& eventVector);
//	long int compute_implicit_L(double tau, vector<AdaptiveSLeaping::Event *>& eventVector, double theta);
//	void sampling(long int Lcurrent);

	// override the standard calculation of propensities
	void computePropensities();
        void      computePropensitiesGrowingVolume(Array< double , 1 > & propensitiesVector, double time, double genTime);
};
