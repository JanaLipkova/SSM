#ifndef MY_RAND_H
#define MY_RAND_H

#include "my_rand.h"
#include <stdlib.h>
#include <random>


namespace myrand{
	std::default_random_engine 			engine;
	std::poisson_distribution<int>  	pois_dist;
	std::gamma_distribution<double> 	gam_dist;
	std::binomial_distribution<int> 	bino_dist;
	std::uniform_real_distribution<double> 	unif_dist(0.0,1.0);
}

#endif
