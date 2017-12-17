/*
 *  SSA_LacZLacY.h
 *  SSM_Xcode
 *
 *  Created by Lipkova on 9/3/14.
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
//#include "HeaderFiles.h"
#include "Method.h"

class SSA_LacZLacY : public Method
{
public:
	SSA_LacZLacY(Simulation * simulation);
	~SSA_LacZLacY();
	
	void _writeDiagnostic(FILE* myfile, int steps, double dt_sum);
	void _writeTrajectories(FILE* myfile, vector<int> data);

	// override the virtual method
	void solve();
private:
};

