/*
 *  SLeapingNonNegative.h
 *  SSM_Xcode
 *
 *  Created by Lipkova on 8/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "LeapMethod.h"

class SLeapingNonNegative: public LeapMethod
{
public:
	SLeapingNonNegative(Simulation * simulation);
	~SLeapingNonNegative();
	
	// override the virtual method
	void solve();
	
private:
	vector<int> listOfCriticalReactions();
	double computeTimeStep(vector<int> criticalReactions);
	double computeTimeOfCritical(vector<int> criticalReactions, double& ac0) ;
	
	long int computeLeapLength(double dt, double a0, double ac0, vector<int> criticalReactions);

	void sampling(long int L, double a0, double ac0, vector<int> criticalReactions, short int crit);
	void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2,std::vector<int> non_critical);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);
	void _executeSSA(double& t, int SSAsteps);
	
};