/*
 *  TauLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class TauLeaping : public LeapMethod
{
public:
	TauLeaping(Simulation * simulation);
	~TauLeaping();

	void _executeSSA(double& t, int SSAsteps);
	void _writeDiagnostic(FILE* myfile, long int L, int steps, long int L_sum, double dt_sum);
	// override the virtual method
	void solve();
private:
	double computeTimeStep();
};

