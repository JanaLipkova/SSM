/*
 *  TauLeapingNonNegative.h
 *  
 *
 *  Created by Lipkova on 3/5/14.
 *  Copyright 2014 Jana Lipkova. All rights reserved.
 *
 */

/* Non negative Tau Leaping method
 * = is an extension of the Tau Leaping by adding a list of critical reactions
 * = critical reactions are solver by SSA
 * = prevence appearance of negative population
 * = for more info see: Cao at al: "Avoiding negative populations in explicit Poisson tau-leaping"
 */

# pragma once
#include "LeapMethod.h"

class TauLeapingNonNegative : public LeapMethod
{
public:
	TauLeapingNonNegative(Simulation * simulation);
	~TauLeapingNonNegative();
	
	// override the virtual method
	void solve();
private:
	double computeTimeStep(vector<int> criticalReactions);
	double computeTimeOfCritical(vector<int> criticalReactions, double& ac0); 
	vector<int> listOfCriticalReactions();
	void _executeSSA(double& t, int SSAsteps);
	void sampling(short int crit, vector<int> criticalReactions, long int & L, double ac0);
	void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::vector<int> non_critical);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);

};


