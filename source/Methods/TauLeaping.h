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

	void solve();
private:
	double computeTimeStep();
    void   executeSSA(double& t, int SSAsteps);
    void   executeSSA_lacZlacY(double& t, int SSAsteps, double genTime);
};

