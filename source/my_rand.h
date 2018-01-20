#ifndef MY_RAND_H
#define MY_RAND_H

#include <stdlib.h>

namespace myrand{

	extern std::default_random_engine 				engine;

	extern std::poisson_distribution<int>  			pois_dist;
	extern std::gamma_distribution<double> 			gam_dist;
	extern std::binomial_distribution<int> 			bino_dist;
	extern std::uniform_real_distribution<double> 	unif_dist;

}

#endif
