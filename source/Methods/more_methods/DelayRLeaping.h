/*
 *  DelayRLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/3/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class DelayRLeaping : public LeapMethod
{
public:
	DelayRLeaping(Simulation * simulation);
	~DelayRLeaping();
		
	// override the virtual method
	void solve();
private:

  #pragma mark  - RLeaping Methods - 
	long int computeLeapLength();
	
	//void computePropensities();

  #pragma mark  - RLeaping Inner Classes - 
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

  #pragma mark  - RLeaping Variables - 
	vector<Event*> eventVector;

  #pragma mark  - Delay Methods - 
  
	// override the standard calculation of propensities
	void computePropensities();


	double randomVariate(double a0, long int L);



};


