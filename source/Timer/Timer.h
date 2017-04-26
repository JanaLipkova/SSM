/*
 *  Timer.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */
 
#pragma once
#include <sys/time.h>
#include <stdlib.h>
//#include "HeaderFiles.h"

class Timer
{
	timeval tv_start, tv_stop;
public:
	void StartSW();
	void StopSW();
	double ReadSW();
};

