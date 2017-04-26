/*
 *  SSALDM.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once

//#include "HeaderFiles.h"
#include "Method.h"

class SSALDM : public Method
{
public:
	SSALDM(Simulation * simulation);
	~SSALDM();
	
	// override the virtual method
	void solve();
private:

	int center( int left, int right );
	int binaryTreeSearch( Array< double , 1 > & cummulativeSum, double key, int left, int right, int & counter );
};
