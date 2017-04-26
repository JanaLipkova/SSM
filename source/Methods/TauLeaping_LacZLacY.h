/*
 *  TauLeaping_LacZLacY.h
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

/* Non negative Tau Leaping method
 * = is an extension of the Tau Leaping by adding a list of critical reactions
 * = critical reactions are solver by SSA
 * = prevence appearance of negative population
 * = for more info see: Cao at al: "Avoiding negative populations in explicit Poisson tau-leaping"
 */

# pragma once
#include "LeapMethod.h"

class TauLeaping_LacZLacY : public LeapMethod
{
public:
	TauLeaping_LacZLacY(Simulation * simulation);
	~TauLeaping_LacZLacY();
	
	// override the virtual method
	void solve();
private:
	double computeTimeStep(vector<int> criticalReactions, vector<int> nonCriticalReactions);
	double computeTimeOfCritical(vector<int> criticalReactions, double& ac0); 
	void identifyReactions(vector<int>& criticalReactions, vector<int>& nonCriticalReactions);
	void _executeSSA(double& t, int SSAsteps, double genTime);
	void sampling(short int crit, vector<int> criticalReactions, vector<int> nonCriticalReactions, long int & L, double ac0);
	void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::vector<int> non_critical);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);
	void _writeTrajectories(FILE* myfile, vector<int> data);

};


