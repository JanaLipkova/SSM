/*
 *  DelayTauLeapingPhilippe.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 10/7/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

#include "LeapMethod.h"

class DelayTauLeapingPhilippe : public LeapMethod
{
public:
	DelayTauLeapingPhilippe(Simulation * simulation);
	~DelayTauLeapingPhilippe();
		
	// override the virtual method
	void solve();
private:
	

};

