/*
 *  SSA.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once
//#include "HeaderFiles.h"
#include "Method.h"

class SSA : public Method
{
public:
	SSA(Simulation * simulation);
	~SSA();
    
	void solve();
private:
};
