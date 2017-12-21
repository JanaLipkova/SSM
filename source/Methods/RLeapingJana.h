/*
 *  RLeapingJana.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 4/24/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#pragma once
#include "LeapMethod.h"

class RLeapingJana : public LeapMethod
{
public:
	RLeapingJana(Simulation * simulation);
	~RLeapingJana();

	// override the virtual method
	void solve();
private:
	long int computeLeapLength();

	// override the standard calculation of propensities
	void computePropensities();
	void sampling(long int L, double a0);

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
