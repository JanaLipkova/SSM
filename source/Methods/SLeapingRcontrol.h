/*
 *  SLeapingRcontrol.h
 *  SSM_Xcode
 *
 *  Created by Lipkova on 8/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */



#pragma once

#include "LeapMethod.h"

class SLeapingRcontrol: public LeapMethod
{
public:
	SLeapingRcontrol(Simulation * simulation);
	~SLeapingRcontrol();
	
	// override the virtual method
	void solve();
	
private:
	double computeTimeStep();
	long int computeLeapLength(double& dt, double a0);
	void sampling(long int L, double a0);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);
	
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