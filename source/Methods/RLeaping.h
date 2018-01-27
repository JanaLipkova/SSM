/*
 *  RLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by  Lipkova on 4/24/13.
 *  Copyright 2013  Lipkova. All rights reserved.
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
	void computePropensitiesGrowingVolume(Array< double , 1 > & propensitiesVector, double time, double genTime);
	void sampling(long int L, double a0);
	

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
