/*
 *  DelayTauLeapingLeier.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class DelayTauLeapingLeier : public LeapMethod
{
public:
	DelayTauLeapingLeier(Simulation * simulation);
	~DelayTauLeapingLeier();
		
	// override the virtual method
	void solve();
private:

};

