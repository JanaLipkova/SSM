 /*
 *  RootFinderJacobian.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 6/3/13.
 *  Copyright 2013 CSE Lab. All rights reserved.
 *
 */

#include "RootFinderJacobian.h"

namespace RootFinderJacobian
{
	Simulation * simulation;
	double tau;
	int numberOfSpecies;    
	int numberOfReactions;  
	int maxIter;    // max number of iteration for the solver
	double scale;
	bool isFirst = true; 

	struct parameters 
	{
		vector<double> par;
	};
	
	parameters B;
}


//************************
//  Set up parameters
//************************
void RootFinderJacobian::RootFinderSetUp(Simulation * simulation_, double tau_, int MaxNumberOfIterations, vector<double> B_ )
{
	RootFinderJacobian::simulation = simulation_;
	RootFinderJacobian::tau = tau_;
	RootFinderJacobian::maxIter = MaxNumberOfIterations;
	Model * sbmlModel = simulation->getSBMLModel();
	RootFinderJacobian::numberOfSpecies	    = sbmlModel->getNumSpecies();
	RootFinderJacobian::numberOfReactions	= sbmlModel->getNumReactions(); 
	RootFinderJacobian::scale = 1;
	
	B.par = B_;
	
//	for (int i=0; i<RootFinderJacobian::numberOfSpecies; i++)
//	{
//		B.par.push_back(B_[i]);
//	}
}

//************************
//   Propensity
// return propensity a_j for given state X
//************************
double RootFinderJacobian::get_propensity(int j, vector<double> X)
{
	int nu;
	//ParticleType xs;
	//ParticleType num, denom;
	double xs;
	double num, denom;
	
	SSMReaction * r = RootFinderJacobian::simulation->ssmReactionList[j];
	vector <int>  reactants		= r->getReactants();
	vector <int>  nu_reactants	= r->getNuReactants();
	
	double propens = r->getRate();
		
	for (int s = 0; s < reactants.size(); ++s)
	{
		nu = nu_reactants[s];
		xs = X[reactants[s]];
		num = xs;
		denom = nu;
		while ((--nu)>0)
		{
			denom *= nu;
			num *= (xs - nu);
		}
		
		propens = propens * ( (double)num / (double)denom ) ;
	}
	
	return propens;
	
}
//************************


//****************************************
// SET UP SYSTEM OF EQUATIONS - input for root finter solver
//****************************************
//  build system of equations whoose roots we are looking for
//   gsl_vector* x => init_gues, or x_old
//   void* params  => pointer to functions parametes, in this case NULL
//   f = f(x,params)  => functions whose roots we look for
//************************
int RootFinderJacobian::implicit_system_f(const gsl_vector* x, void* params, gsl_vector* f )
{
		
		/* 1. init_gues */
		vector<double> x_old(RootFinderJacobian::numberOfSpecies);
	
		for (int i=0; i<RootFinderJacobian::numberOfSpecies; i++)
			x_old[i] = (double)gsl_vector_get(x,i);
	
	
		/* 2. define system of equations 
		// Z(:) = X(t+tau) - sum_j(vj*aj(t=tau)*tau) - B(X(t)) = 0 */
		vector<double> z(RootFinderJacobian::numberOfSpecies,0);
		
		for (int i=0 ; i<RootFinderJacobian::numberOfSpecies; i++)
			z[i] = x_old[i] - RootFinderJacobian::B.par[i];

		
		for (int j = 0; j<RootFinderJacobian::numberOfReactions; j++)
		{
			SSMReaction * ri = simulation->ssmReactionList[j];		
			const vector<int> & changes = ri->getChanges();
			const vector<int> & nuChanges = ri->getNuChanges();
			
			for (int ns = 0; ns < changes.size(); ns++)
				z[changes[ns]] -= nuChanges[ns]* get_propensity(j,x_old) * (RootFinderJacobian::tau);

		}
	
	
	// security check
	for (int i = 0; i < RootFinderJacobian::numberOfSpecies; i++)
		if (isnan(z[i]) )
		{
			cout << "Root finder can't stand numbers bigger than e+10" << endl;
			abort();
		}
	
	// set system of equations for the solver
	for (int i = 0; i<RootFinderJacobian::numberOfSpecies; i++)
		gsl_vector_set(f,i,z[i]);
	
	return GSL_SUCCESS;	
}
//********************************


// !!!!!! Jacobian now hardcoded for the Decaying dimerization case !!!!!!!

//********************************
// Jacobian of the system
//********************************
int RootFinderJacobian::implicit_system_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
    
#ifdef BLABLA
	const double x0 = (double)gsl_vector_get(x,0);  // get species S1
	
	// get reaction rates
	vector<double> c(4,0);
	for (int j=0; j< RootFinderJacobian::numberOfReactions; ++j)
	{
		SSMReaction * r = RootFinderJacobian::simulation->ssmReactionList[j];
		c[j] = r->getRate();
	}
	c[1] = c[1]*0.5; // since it is second order reaction !!! NOT SURE ABOUT THIS

	// hardcode stochiometric coefficients for dec.dimerization
	vector<int> v0(3,0);	
	vector<int> v1(3,0);
	vector<int> v2(3,0);
	vector<int> v3(3,0);
	
	v0[0] = -1;
	
	v1[0] = -2;
	v1[1] =  1;
	
	v2[0] =  2;
	v2[1] = -1;
	
	v3[1] = -1;
	v3[2] =  1;
	
	
	 double t = RootFinderJacobian::tau;
	
	 // jacobian inputs, where dfij = @f_{i}\@x_{j}
	 double df00 = 1. - v0[0] * c[0] * t - 2. * x0 * v1[0] * c[1] * t + v1[0] * c[1] * t;
	 double df01 = - v2[0] * c[2] * t - v3[0] * c[3] * t;
	 double df02 = 0.;
	
	 double df10 = - v0[1] * c[0] * t - 2. * x0 * v1[1] * c[1] * t + v1[1] * c[1] * t;
	 double df11 = 1. - v2[1] * c[2] * t - v3[1] * c[3] * t;
	 double df12 = 0.;
	
	 double df20 = -v0[2] * c[0] * t - 2. * x0 * v1[2] * c[1] * t + v1[2] * c[1] * t;
	 double df21 = -v2[2] * c[2] * t - v3[2] * c[3] * t;
	 double df22 = 1.;
	
	vector<double> df(9,0);
	
	df[0] = df00;
	df[1] = df01;
	df[2] = df02;
	
	df[3] = df10;
	df[4] = df11;
	df[5] = df12;
	
	df[6] = df20;
	df[7] = df21;
	df[8] = df22;
	
		
	gsl_matrix_set (J, 0, 0, df[0]);
	gsl_matrix_set (J, 0, 1, df[1]);
	gsl_matrix_set (J, 0, 2, df[2]);

	gsl_matrix_set (J, 1, 0, df[3]);
	gsl_matrix_set (J, 1, 1, df[4]);
	gsl_matrix_set (J, 1, 2, df[5]);

	gsl_matrix_set (J, 2, 0, df[6]);
	gsl_matrix_set (J, 2, 1, df[7]);
	gsl_matrix_set (J, 2, 2, df[8]);
    
#endif
    
    
    // Jacobian from the file
    int M = RootFinderJacobian::numberOfReactions;
    int N = RootFinderJacobian::numberOfSpecies;

    vector<double> Jacobian(M*M, 0.);   // allocate jacobian into vector
    
    // get reaction rates
    vector<double> rates(M,0.);
    for (int j=0; j< RootFinderJacobian::numberOfReactions; ++j)
    {
        SSMReaction * r = RootFinderJacobian::simulation->ssmReactionList[j];
        rates[j] = r->getRate();
    }
    
    // get vector X with current state of the system
    vector<double> X(N,0.);   // doubles since Newton-Rhapson roots might not be integer
    for(int i=1; i<N; i++)
        X[i] = (double)gsl_vector_get(x,i);
    
    double tau = RootFinderJacobian::tau;
    computeJacobian(Jacobian,X,rates,tau);
    
    // fill in gsl_matrix with jacobian
    for(int ix=0; ix < N; ++ix)
        for(int iy=0; iy < N; ++iy)
            gsl_matrix_set (J, ix, iy, Jacobian[ix + N*iy]);

        
	return GSL_SUCCESS;
}
//********************************


//********************************
// launch both f and jacobian
//********************************
int RootFinderJacobian::implicit_system_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
	implicit_system_f (x, params, f);
	implicit_system_df (x, params, J);
	
	return GSL_SUCCESS;
}
//********************************



//********************************
//  main part, set up system 
//  and iterate to find roots
//********************************
void RootFinderJacobian::find_roots( vector<double>& output, vector<double>& implicitPropensities, Array<ParticleType,1> init_guess)
{

	int status;
	size_t iter = 0;
	
	//1. allocate pointers for solver:	
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;
	
	// 2. define system of equations 
	gsl_multiroot_function_fdf f = {&implicit_system_f,
		&implicit_system_df,
		&implicit_system_fdf,
		 RootFinderJacobian::numberOfSpecies, NULL};
	
	//3. set initial guess
		gsl_vector* x = gsl_vector_alloc(RootFinderJacobian::numberOfSpecies);
		for (int i = 0; i<RootFinderJacobian::numberOfSpecies; i++)
			gsl_vector_set(x,i,(double)init_guess(i) );

	
	//4. choose and allocate solver  newton
	T = gsl_multiroot_fdfsolver_hybridj;  // options: "newton", "hybridj", "hybridsj", "gnewton"
	s = gsl_multiroot_fdfsolver_alloc (T, RootFinderJacobian::numberOfSpecies);
	gsl_multiroot_fdfsolver_set (s, &f, x);
	//	print_state (iter, s, RootFinderJacobian::numberOfSpecies);
	
	
	/* 5. iterate to find roots */
	do
	{
		iter++;
		
		status = gsl_multiroot_fdfsolver_iterate (s);
		//print_state (iter, s, RootFinderJacobian::numberOfSpecies);
		
		if (status)
			break;
		
		status = gsl_multiroot_test_residual (s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < RootFinderJacobian::maxIter);
	
	//6. fill output vector
	for (int i=0; i<RootFinderJacobian::numberOfSpecies; i++)
		output[i] = (gsl_vector_get(s->x,i));
	
	
	for (int j=0; j< RootFinderJacobian::numberOfReactions; j++)
		implicitPropensities[j] = get_propensity(j,output);

	
	// 7. clean up
	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (x);
}
//********************************


//****************
// PRINT ROOTS AND FUNC.VALUE in founded roots
//****************
void RootFinderJacobian::print_state (size_t iter, gsl_multiroot_fdfsolver * s, int N)
{
	cout <<endl;
	cout << "iteration:" << iter << endl;
	cout << "roots " << " " << " f(root) " << endl;
	for (int i=0; i<N; i++)
	{
		cout<<setprecision(9) <<"x=" << gsl_vector_get(s->x, i) << " f(x)=" << gsl_vector_get(s->f,i) << endl;
	}
	
	cout << endl;
}
//********************************
