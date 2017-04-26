/*
 *  RLeaping_LacZLacY.h
 *  SSM_Xcode
 *
 *  Created by Lipkova on 9/4/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

/* LacZ/LacY system:
 = system with 22 reactions, 19 species
 = with liearly growing volume of reaction enviroment 
 => propensities of higner order reactions must be rescaled accordingly each step
 = amount of some reactants is drawn from random pools
 */


#pragma once

#include "LeapMethod.h"

class RLeaping_LacZLacY : public LeapMethod
{
public:
	RLeaping_LacZLacY(Simulation * simulation);
	~RLeaping_LacZLacY();
	
	// override the virtual method
	void solve();
private:
	long int computeLeapLength();
	
	// override the standard calculation of propensities
	void computePropensitiesGrowingVolume(double genTime);
	void sampling(long int L, double a0);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);
	void _writeTrajectories(FILE* myfile, vector<int> data);
	
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


