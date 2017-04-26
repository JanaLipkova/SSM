/*
 *  RootFinderJacobian.h
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 6/3/13.
 *  Copyright 2013 CSE lab. All rights reserved.
 *
 */

// AIM: Multidimensional root finidng using GNU library with Jacobian as input
// for more info see: http://www.gnu.org/software/gsl/manual/html_node/Overview-of-Multidimensional-Root-Finding.html

#pragma once

#include "../HeaderFiles.h"
#include "../Simulation.h"

namespace RootFinderJacobian
{
	void RootFinderSetUp(Simulation * simulation, double tau, int MaxNumberOfIterations, vector<double> B);
	
	void find_roots( vector<double>& output, vector<double>& impicitPropensities, Array<ParticleType,1> init_guess);
	void print_state (size_t iter, gsl_multiroot_fdfsolver * s, int N);
	
	double get_propensity(int j, vector< double> X);
	int implicit_system_f(const gsl_vector* x, void* params, gsl_vector* f );
	int implicit_system_df (const gsl_vector * x, void *params, gsl_matrix * J);
	int implicit_system_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);


};

