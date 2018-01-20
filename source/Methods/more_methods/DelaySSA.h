/*
 *  DelaySSA.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/6/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

//#include "HeaderFiles.h"
#include "Method.h"

class DelaySSA : public Method
{
public:
	DelaySSA(Simulation * simulation);
	~DelaySSA();
		
	// override the virtual method
	void solve();
private:
	
	
	double randomVariate(double a0);
	
	
	// Create a Mersenne twister random number generator
	boost::mt19937 rng;

};






