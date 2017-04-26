/*
 *  SLeaping_LacZLacY.h
 *  SSM_Xcode
 *
 *  Created by Lipkova on 9/8/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

/* LacZ/LacY system:
 = system with 22 reactions, 19 species
 = with linearly growing volume of reaction enviroment 
 => propensities of higner order reactions must be rescaled accordingly each step
 = amount of some reactants is drawn from random pools
 */


#pragma once
#include "LeapMethod.h"

class SLeaping_LacZLacY : public LeapMethod
{
public:
	SLeaping_LacZLacY(Simulation * simulation);
	~SLeaping_LacZLacY();
	
	// override the virtual method
	void solve();
private:
	
	void     identifyReactions(vector<int>& criticalReactions, vector<int>& nonCriticalReactions);
	double   computeTimeStep(vector<int> criticalReactions, vector<int> nonCriticalReactions);
	void	 computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::vector<int> nonCriticalReactions);	
	double   computeTimeOfCritical(vector<int> criticalReactions, double& ac0); 
	
	void	 _executeSSA(double& t, int SSAsteps, double genTime);
	long int computeLeapLength(double dt, double a0, double ac0, vector<int> criticalReactions);
	void     sampling(long int L, double a0, double ac0, vector<int> criticalReactions, vector<int> nonCriticalReactions, short int crit);

	// override the standard calculation of propensities
	void computePropensitiesGrowingVolume(double genTime);

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


