/*
 *  RootFinder.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Jana Lipkova on 5/9/13.
 *  Copyright 2013 Jana Lipkova. All rights reserved.
 *
 */

#include "RootFinder.h"

namespace RootFinder
{
	Simulation * simulation;
	double tau;
	int numberOfSpecies;    
	int numberOfReactions;  
	int maxIter;    // max number of iteration for the solver
	
	struct parameters 
	{
		vector<double> par;
	};
	
	parameters B;
}

void RootFinder::RootFinderSetUp(Simulation * simulation_, double tau_, int MaxNumberOfIterations, vector<double> B_ )
{
	RootFinder::simulation = simulation_;
	RootFinder::tau = tau_;
	RootFinder::maxIter = MaxNumberOfIterations;
	Model * sbmlModel = simulation->getSBMLModel();
	RootFinder::numberOfSpecies	    = sbmlModel->getNumSpecies();
	RootFinder::numberOfReactions	= sbmlModel->getNumReactions(); 
		
	for (int i=0; i<RootFinder::numberOfSpecies; i++)
	{
		B.par.push_back(B_[i]);
	}
}


//************************
//   Propensity
// return propensity a_j for given state X
//************************
double RootFinder::get_propensity(int j, vector<double> X)
{
	int nu;
	ParticleType xs;
	ParticleType num, denom;
	
	SSMReaction * r = RootFinder::simulation->ssmReactionList[j];
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
int RootFinder::implicit_system(const gsl_vector* x, void* params, gsl_vector* f )
{
	
	
	/* 1. init_gues */
	vector<double> x_old(RootFinder::numberOfSpecies);
	double scale = 1;
	vector<double> z(RootFinder::numberOfSpecies,0);

	
	for (int i=0; i<RootFinder::numberOfSpecies; i++)
	{
		x_old[i] = gsl_vector_get(x,i);
	}
	
	
	/* 2. define system of equations 
	 // Z(:) = X(t+tau) - sum_j(vj*aj(t=tau)*tau) - B(X(t)) = 0 */
	for (int i=0 ; i<RootFinder::numberOfSpecies; i++)
	{ 
		z[i] += x_old[i] - RootFinder::B.par[i];
	}
		
	for (int j = 0; j<RootFinder::numberOfReactions; j++)
	{
		SSMReaction * ri = simulation->ssmReactionList[j];		
		const vector<int> & changes = ri->getChanges();
		const vector<int> & nuChanges = ri->getNuChanges();
		
		for (int ns = 0; ns < changes.size(); ns++)
			z[changes[ns]] -= nuChanges[ns]* get_propensity(j,x_old) * (RootFinder::tau);
	}
//	cout << endl;
//	cout << " call one" << endl;
//	cout << endl;
	
// find the order of each z componenet and rescale but the same factor
	for (int i = 0; i < 3; i++)
	{
		//cout << "z["<<i<<"]="<<  z[i] << endl;
		int num = (int) abs(z[i]);
		double mag = 0;

		while (num > 0.1)
		{
			mag++;
			num = num*0.1;;
		}
		
		scale = max(scale, mag);
		cout << "scale="<<scale<<endl;
	}
	
	double temp = pow(10, scale);
	
	for (int i=0; i < 3; i++)
	{
		z[i] = z[i]/temp;
	}
	

	
	
	for (int i = 0; i < RootFinder::numberOfSpecies; i++)
		if (isnan(z[i]) )
		{
			cout << "Root finder can't stand numbers bigger than e+10" << endl;
			abort();
		}
	
	// set system of equations for the solver
	for (int i = 0; i<RootFinder::numberOfSpecies; i++)
		gsl_vector_set(f,i,z[i]);
	
	return GSL_SUCCESS;	
	
	
	
//	
//	/* 1. init_gues */
//	vector<double> x_old(RootFinder::numberOfSpecies);
//	
//	for (int i=0; i<RootFinder::numberOfSpecies; i++)
//	{
//		x_old[i] = gsl_vector_get(x,i);
////		cout << "x_old["<<i<<"]="<<x_old[i]<<endl;
//	}
//
//	
//	/* 2. define system of equations 
//	// Z(:) = X(t+tau) - sum_j(vj*aj(t=tau)*tau) - B(X(t)) = 0 */
//	vector<double> z(RootFinder::numberOfSpecies,0);
//	
//	for (int i=0 ; i<RootFinder::numberOfSpecies; i++)
//		z[i] = x_old[i] - RootFinder::B.par[i];
//
//	
//	for (int j = 0; j<RootFinder::numberOfReactions; j++)
//	{
//		SSMReaction * ri = simulation->ssmReactionList[j];		
//		const vector<int> & changes = ri->getChanges();
//		const vector<int> & nuChanges = ri->getNuChanges();
//		
//		for (int ns = 0; ns < changes.size(); ns++)
//			z[changes[ns]] -= nuChanges[ns]* get_propensity(j,x_old) * (RootFinder::tau);
//		
////		cout <<"propensity["<<j<<"]="<<get_propensity(j,x_old)<<endl;
////		cout << "aj*tau"<< get_propensity(j,x_old) * (RootFinder::tau)<<endl;
//	}
//		
//	
//	//----
//	
////	for (int i = 0; i < RootFinder::numberOfSpecies; i++)
////	{	cout << "z["<<i<<"]="<<z[i]<<endl;
////			z[i] = z[i]/(abs(round(RootFinder::B.par[i])));
////	}
//	//----
//	
//	for (int i = 0; i < RootFinder::numberOfSpecies; i++)
//         if (isnan(z[i]) )
//		 {
//			 cout << "Root finder can't stand numbers bigger than e+10" << endl;
//			 abort();
//		 }
//	
//	// set system of equations for the solver
//	for (int i = 0; i<RootFinder::numberOfSpecies; i++)
//		gsl_vector_set(f,i,z[i]);
//
//	return GSL_SUCCESS;	
}
//********************************


//********************************
//   ROOT FINDER METHOD
//********************************
void RootFinder::find_roots( vector<double>& output, vector<double>& implicitPropensities, Array<ParticleType,1> init_guess)
{
	int status;
	int iter = 0;
	
	//1. allocate pointers for solver:	
	const gsl_multiroot_fsolver_type* T ;    //returns a pointer to a newly allocated instance of a solver of type T for a system of n equations															
	gsl_multiroot_fsolver* s;                //allocat hybrid algorithm          
	
	// 2. define system of equations 
	gsl_multiroot_function f = {&implicit_system, RootFinder::numberOfSpecies, NULL};            // input in form { system of equations, dimenstion, pointer to parameters}
	
	vector<double> IC(3);
	
//	IC[0] = 2970.; 
//	IC[1] = 40155.;
//	IC[2] = 3445.;
//	
//	
//	/* small tau*/
	IC[0] = 4150.; 
	IC[1] = 39565;
	IC[2] = 3445.;
//	
	gsl_vector* x = gsl_vector_alloc(RootFinder::numberOfSpecies);
	for (int i = 0; i<RootFinder::numberOfSpecies; i++)
		gsl_vector_set(x,i,(double)IC[i] );
	
	//3. set initial guess
//	gsl_vector* x = gsl_vector_alloc(RootFinder::numberOfSpecies);
//	for (int i = 0; i<RootFinder::numberOfSpecies; i++)
//		gsl_vector_set(x,i,(double)init_guess(i) );

	//3. choose and allocate solver
	T = gsl_multiroot_fsolver_hybrids;           
	s = gsl_multiroot_fsolver_alloc (T, RootFinder::numberOfSpecies);   // allocat hybrid algorithm for 2 dim problem
	
	// 4. link solver with your system
	gsl_multiroot_fsolver_set (s, &f, x);
	
	print_state (iter, s, RootFinder::numberOfSpecies);
	
	/* 5. iterate to find roots */
	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);
		
		print_state (iter, s, RootFinder::numberOfSpecies);
		
		if (status)   /* check if solver is stuck */
			break;
		
		status =
		gsl_multiroot_test_residual (s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < RootFinder::maxIter);
	
	// fill the output vector
	for (int i=0; i<RootFinder::numberOfSpecies; i++)
		{
			output[i] = (gsl_vector_get(s->x,i));
			//cout<< "Output[i]=" << output[i] << endl;
		}
	
	for (int j=0; j< RootFinder::numberOfReactions; j++)
	{
		implicitPropensities[j] = get_propensity(j,output);
		//cout<<"propensity[j]"<<implicitPropensities[j]<< endl;
	}
	
	//clean up
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
	
}
//****************

//****************
// PRINT ROOTS AND FUNC.VALUE in founded roots
//****************
void RootFinder::print_state (size_t iter, gsl_multiroot_fsolver * s, int N)
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


