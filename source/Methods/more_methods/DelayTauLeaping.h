/*
 *  DelayTauLeaping.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class DelayTauLeaping : public LeapMethod
{
public:
	DelayTauLeaping(Simulation * simulation);
	~DelayTauLeaping();
		
	// override the virtual method
	void solve();
private:
	

};

