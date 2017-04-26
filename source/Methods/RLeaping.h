/*
 *  RLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/27/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class RLeaping : public LeapMethod
{
public:
	RLeaping(Simulation * simulation);
	~RLeaping();
	
	// override the virtual method
	void solve();
private:
	long int computeLeapLength();
	
	// override the standard calculation of propensities
	void computePropensities();
	
	// anonymous inner class, R-Leaping needs to store the indices and propensities of reactions
	class Event
	{
		public:
			Event() {}
			~Event() {}
			
			double	propensity;
			int		index;
	};

	// sorts the global reactions with respect to the largest propensities (i.e. descending order 
	// of propensities - R-Leaping requirement
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

