/*
 *  SSMReaction.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/27/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */
  
 /*  
	Example of the usage of this class
 
 Gray-Scott reactions:
 
 Reaction * r0 = new Reaction("U + 2V -> 3V", 2*(dc*dc)*amplification );
	Reaction * r1 = new Reaction("V -> 0", (F+kappa )*amplification);
	Reaction * r2 = new Reaction("U -> 0", (F)*amplification);
	Reaction * r3 = new Reaction("0 -> U", F/dc*amplification);
	
// U + 2V -> 3V
	r0->addReactant(0);
	r0->addReactant(1);
	r0->setLastReactantNu(2);
	r0->addProduct(1);
	r0->setLastProductNu(3);
	r0->finalizeReaction();

// V -> 0
	r1->addReactant(1);
	r1->finalizeReaction();

// U -> 0
	r2->addReactant(0);
	r2->finalizeReaction();

// 0 -> U
	r3->addProduct(0);
	r3->finalizeReaction();
 
 */
 
 
#pragma once
#include "HeaderFiles.h"

class SSMReaction
{
public:
	// Empty Constructor
	SSMReaction();
	
	// Overloaded constructors
	SSMReaction(string name);
	SSMReaction(string name, double rate);
	SSMReaction(string name, double rate, unsigned int nReactants,  unsigned int nProducts, ... );

	~SSMReaction();

	// SSMReaction member methods
	void		addReactant(int index);		// the index of the Species vector [0...N-1]
	void		setLastReactantNu(int nu);
	
	void		addProduct(int index);		// the index of the Species vector [0...N-1]
	void		setLastProductNu(int nu);
	
	void		finalizeReaction();	// this method calls 3 private methods:
									// sortSpecies(), computeChanges(), and computeOrder()

	// SSMReaction set methods
	void setPropensity(double _propensity){ propensity = _propensity; }
	void setRate(double rate);

	// SSMReaction get methods
	string		getName();
	int			getOrder();
	
	//access to info, as fast as possible
	int getChangesSize() const { return changes.size();}
	
	const vector<int> & getReactants() const { return reactants;}
	int getReactants(int index) const { return reactants[index]; }
	
	int getReactantsSize() const	{ return reactants.size();}
	const vector<int> & getNuReactants() const { return nuReactants;}
	
	int getNuReactants(int index) const { return nuReactants[index];}
	const vector<int> & getProducts() const { return products; }
	
	int getProducts(int index) const { return products[index];}
	int	getProductsSize() const { return products.size();}
	
	const vector<int> & getNuProducts() const { return nuProducts;}
	int	getNuProducts(int index)	const { return nuProducts[index];}
	
	const vector<int> & getChanges() const {return changes;}
	int getChanges(int index) const { return changes[index];}
	
	const vector<int> & getNuChanges() const{ return nuChanges;}
	int getNuChanges(int index) const { return nuChanges[index];}
	
	double getRate() const{ return rate;}
	double getPropensity() const {return propensity;}

	void		toString(); // for testing

private:
	// private methods
	void		sortSpecies();
	void		computeChanges();
	void		computeOrder();
	void		Init();

	// private variables
	string		name;			// the name of the SSMReaction
	
	vector<int> reactants;		// an STL vector of reactants
	vector<int> nuReactants;	// an STL vector of stoichiometric coefficients for the reactants
	vector<int> products;		// an STL vector of products
	vector<int> nuProducts;		// an STL vector of stoichiometric coefficients for the products
	vector<double> gradPropensity;	// an STL vector for the gradient of the propensity which is 
								// used by some of the leaping methods
	vector<int> changes;		// an STL vector of changes in the SSMReaction (usually denoted as 
								// v in literature)
	vector<int> nuChanges;		// an STL vector of stoichiometric coefficients for the changes
	
	double		rate;			// the rate of the SSMReaction (usually denoted by c in literature)
	double		propensity;		// the propensity of the SSMReaction (usually denoted a in literature)
	int			order;			// the order of the SSMReaction
	
	
};

