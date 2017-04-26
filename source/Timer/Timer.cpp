/*
 *  Timer.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */
 
#include "Timer.h"

void Timer::StartSW()
{
	gettimeofday(&tv_start, NULL);
}

void Timer::StopSW()
{
	gettimeofday(&tv_stop, NULL); 
	
}

double Timer::ReadSW()
{
	long isec, iusec;
	double time;
	isec = tv_stop.tv_sec - tv_start.tv_sec;
	iusec = tv_stop.tv_usec - tv_start.tv_usec;
	time = (double) isec + ( (double) iusec ) / 1000000.;
	return((double)time);
}



