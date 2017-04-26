/*
 *  SSMReaction.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/27/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "SSMReaction.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// overloaded constructors
SSMReaction::SSMReaction()
{
	name = "";
	rate = 0.0;
	propensity = 0.0;
	order = 0;

}

SSMReaction::SSMReaction(string name)
{
	Init();
	this->name = name;
}

SSMReaction::SSMReaction(string name, double rate)
{
	Init();
	this->name = name;
	this->rate = rate;
}

SSMReaction::SSMReaction(string name, double rate, unsigned int nArgsReactant,  unsigned int nArgsProduct, ... )
{
	Init();
	this->name = name;
	this->rate = rate;
	
	va_list arglist;
	va_start(arglist,nArgsProduct);
	for(unsigned int i=0; i<nArgsReactant; i++) 
	{
		addReactant(va_arg(arglist, unsigned int));
		setLastReactantNu(1); //later we should change it
	}
	
	for(unsigned int i=0; i<nArgsProduct; i++) 
	{
		addProduct(va_arg(arglist, unsigned int));
		setLastProductNu(1);//later we should change it
	}
	
	va_end(arglist);
}

SSMReaction::~SSMReaction(){}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SSMReaction::addReactant(int index)
{
	reactants.push_back(index);
	nuReactants.push_back(1);
}

void SSMReaction::setLastReactantNu(int nu)
{
	int reactantVectorSize = reactants.size();
	nuReactants[ (reactantVectorSize - 1) ] = nu;
}

void SSMReaction::addProduct(int index)
{
	products.push_back(index);
	nuProducts.push_back(1);
}	

void SSMReaction::setLastProductNu(int nu)
{
	int productVectorSize = products.size();
	nuProducts[ (productVectorSize - 1) ] = nu;	
}

void SSMReaction::finalizeReaction()
{
	sortSpecies();
	computeChanges();
	computeOrder();
	
	// make room for the gradPropensity vector
	for (int i = 0; i < (int)reactants.size(); ++i)
	{
		gradPropensity.push_back(0);
	}
}	

void SSMReaction::setRate(double rate)
{
	this->rate = rate;
}

string SSMReaction::getName()
{
	return name;
}

int	SSMReaction::getOrder()
{
	return order;
}

void SSMReaction::toString()
{
	cout << "name: " << name << " propensity: " << propensity << endl;
	
	cout << "changes info: " << endl;
	cout << "changes: " << endl;
	for (int i = 0; i < changes.size(); ++i) {cout << changes[i] << " " ;} cout << endl;
	
	cout << "changes nu: " << endl;
	for (int i = 0; i < nuChanges.size(); ++i) {cout << nuChanges[i] << " " ;} cout << endl;
	
	cout << endl << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct idnu_s
{
	int id;
	int nu;
} idnu_t;	

// private methods
int	sortComparison(const void * a, const void * b)
{
	return(((idnu_t *)a)->id - ((idnu_t *)b)->id);
}

void SSMReaction::sortSpecies()
{
	/*
	 This function is called when a SSMReaction is finalized (end=SSMReaction).
	 It sorts the reactants and products by their index.
	 This will make a fast sparse matrix multiplication possible 
	 as needed by approximate leaping schemes
	 */
	unsigned int is;
	
	idnu_t * array = (idnu_t *)calloc(reactants.size(),sizeof(idnu_t));
	for (is=0; is < reactants.size(); ++is)
	{
		array[is].id = reactants[is];
		array[is].nu = nuReactants[is];
	}
	qsort(array,reactants.size(),sizeof(idnu_t),
		  sortComparison);
	for (is=0; is < reactants.size(); ++is)
	{
		reactants[is] = array[is].id; 
		nuReactants[is] = array[is].nu;
	}
	
	array = (idnu_t *)realloc(array, products.size() * sizeof(idnu_t));
	for (is=0;is<products.size();is++)
	{
		array[is].id = products[is];
		array[is].nu = nuProducts[is];
	}
	qsort(array,products.size(),sizeof(idnu_t),
		  sortComparison);
	for (is = 0; is < products.size(); ++is)
	{
		products[is] = array[is].id; 
		nuProducts[is] = array[is].nu;
	}
	
	free ( array );
	//for (is = 0; is < reactants.size(); ++is)
	//	cout << " " << reactants[is] << " " << nuReactants[is] << endl;
}

void SSMReaction::computeChanges()
{

	int is, js, ic;
	int changesNumber;
	
	// allocate space for the changes vectors
	for (unsigned int i = 0; i < (reactants.size() + products.size()); ++i)
	{
		changes.push_back(0);
		nuChanges.push_back(0);
	}
	
	
	
	// Merging the nu reactants and nu products into a 
	// Changes array	
	for (is = 0; is < reactants.size(); ++is)
	{
		changes[is] = reactants[is];
		nuChanges[is] = -nuReactants[is];
		//changes.push_back(reactants[is]);
		//nuChanges.push_back( -nuReactants[is] );
 	}
	
	
	//ic = reactants.size()+1;
	ic = reactants.size();
	is = 0; 
	js = 0;
	while ( (is < reactants.size()) && (js < products.size()) )
	{
		if (changes[is] > products[js])
		{
			changes[ic] = products[js];
			nuChanges[ic] = nuProducts[js];
			++ic;
			++js;
		}
		else if (changes[is] < products[js])
		{
			++is;
		}
		else
		{
			nuChanges[is] += nuProducts[js];
			++is;
			++js;
		}
	}
	for (; js < products.size(); ++js)
	{
		changes[ic] = products[js];
		nuChanges[ic] = nuProducts[js];
		++ic;
	}
	changesNumber = ic;
	//cout << "ic: " << ic << endl;
	// remove the zero elements of the changes and nuChanges vectors
	is = 0;
	//while (is < changesNumber)
	while (is < changes.size())
	{
		if (nuChanges[is] == 0)
		{
			nuChanges.erase( (nuChanges.begin()+is), (nuChanges.begin()+is+1) );
			changes.erase( (changes.begin()+is), (changes.begin()+is+1) );
			changesNumber -= 1;
			//cout << "is: " << is << endl;
		}
		else
		{
			++is;
		}
	}
}

void SSMReaction::computeOrder()
{
	order = 0;
	for (int i = 0; i < (int)reactants.size(); ++i)
	{
		order += nuReactants[i];
	}
}

void SSMReaction::Init()
{
	name = "";
	rate = 0.0;
	propensity = 0.0;
	order = 0;
}




