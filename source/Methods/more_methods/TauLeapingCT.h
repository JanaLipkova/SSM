/*
 *  TauLeapingCTS.h
 *  StochasticSimulationMethods
 *
 *  Created by Roger Rossé on 02/12/12.
 *  Copyright 2012 Roger Rossé. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class TauLeapingCT : public LeapMethod
{
public:
	TauLeapingCT(Simulation * simulation);
	~TauLeapingCT();
	
	// override the virtual method
	void solve();
private:
	double computeTimeStep();
};

