/*
 *  RootFinder.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 5/9/13.
 *  Copyright 2013 CSE lab. All rights reserved.
 *
 */

// AIM: Multidimensional root finidng using GNU library
// for more info see: http://www.gnu.org/software/gsl/manual/html_node/Overview-of-Multidimensional-Root-Finding.html

#pragma once

#include "../HeaderFiles.h"
#include "../Simulation.h"

namespace RootFinder
{
	void RootFinderSetUp(Simulation * simulation, double tau, int MaxNumberOfIterations, vector<double> B);
		
	void find_roots( vector<double>& output, vector<double>& impicitPropensities, Array<ParticleType,1> init_guess);
	void print_state (size_t iter, gsl_multiroot_fsolver * s, int N);
	
	double get_propensity(int j, vector< double> X);
	int implicit_system(const gsl_vector* x, void* params, gsl_vector* f );
};

