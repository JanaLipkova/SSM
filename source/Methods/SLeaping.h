/*
 *  SLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 4/25/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class SLeaping: public LeapMethod
{
public:
	SLeaping(Simulation * simulation);
	~SLeaping();
	
	// override the virtual method
	void solve();
	
private:
	double computeTimeStep();
	long int computeLeapLength(double dt, double a0);
	void sampling(long int L, double a0);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);
	void _executeSSA(double& t,  int SSAsteps);

	// override the standard calculation of propensities
	void computePropensities();
	
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