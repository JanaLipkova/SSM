/*
 *  StochasticSimulationMethodsGenerateXML.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

// C/C++, Blitz++, RanLib, SBML, and Boost
#include "HeaderFiles.h"


// "main" to generate the species an reactions for a SBML xml file
// 1-D Diffusion

long int	omega( double x, double nu, double t, double q )
{
	return (long int)((( ( x * exp ( - (x * x) / ( 1.0 + 4.0*nu*t ) )  )  / pow( (1.0 + 4.0*nu*t), (3.0/2.0) ) )*q) + 0.5);
}

void		printReaction( int Rn, int i, int j, double rate )
{
	cout << "<reaction name=\"R" << Rn << "\">" << endl;
	cout << "<listOfReactants>" << endl;
	cout << " <speciesReference species=\"S" << i << "\" />" << endl;
	cout << "</listOfReactants>" << endl;
	cout << "<listOfProducts>" << endl;
	cout << " <speciesReference species=\"S" << j << "\" />" << endl;
	cout << "</listOfProducts>" << endl;
	cout << "<kineticLaw formula=\"\">" << endl;
	cout << "<listOfParameters>" << endl;
	cout << "<parameter value=\"" << rate << "\" name=\"c" << i << "\" />" << endl;
	cout << "</listOfParameters>" << endl;
	cout << "</kineticLaw>" << endl;
	cout << "</reaction>" << endl;
	cout << endl;
}

void		printReaction( int Rn, int i, int j, double rate, string speciesName )
{
	cout << "<reaction name=\"R" << Rn << "\">" << endl;
	cout << "<listOfReactants>" << endl;
	cout << " <speciesReference species=\"" << speciesName << i << "\" />" << endl;
	cout << "</listOfReactants>" << endl;
	cout << "<listOfProducts>" << endl;
	cout << " <speciesReference species=\"" << speciesName << j << "\" />" << endl;
	cout << "</listOfProducts>" << endl;
	cout << "<kineticLaw formula=\"\">" << endl;
	cout << "<listOfParameters>" << endl;
	cout << "<parameter value=\"" << rate << "\" name=\"c" << i << "\" />" << endl;
	cout << "</listOfParameters>" << endl;
	cout << "</kineticLaw>" << endl;
	cout << "</reaction>" << endl;
	cout << endl;
}

void		printReactionFisher( int Rn, int i, int j, double rate, string speciesNameA, string speciesNameB )
{
	cout << "<reaction name=\"R" << Rn << "\">" << endl;
	cout << "<listOfReactants>" << endl;
	cout << " <speciesReference species=\"" << speciesNameA << i << "\" />" << endl;
	cout << " <speciesReference species=\"" << speciesNameB << i << "\" />" << endl;
	cout << "</listOfReactants>" << endl;
	cout << "<listOfProducts>" << endl;
	cout << " <speciesReference species=\"" << speciesNameA << j << "\"" << " stoichiometry=\"2\"" << " />" << endl;
	cout << "</listOfProducts>" << endl;
	cout << "<kineticLaw formula=\"\">" << endl;
	cout << "<listOfParameters>" << endl;
	cout << "<parameter value=\"" << rate << "\" name=\"c" << i << "\" />" << endl;
	cout << "</listOfParameters>" << endl;
	cout << "</kineticLaw>" << endl;
	cout << "</reaction>" << endl;
	cout << endl;
}

void		printReaction( int Rn, int i, double rate )
{
	cout << "<reaction name=\"R" << Rn << "\">" << endl;
	cout << "<listOfReactants>" << endl;
	cout << " <speciesReference species=\"S" << i << "\" />" << endl;
	cout << "</listOfReactants>" << endl;
	cout << "<listOfProducts>" << endl;
	cout << "</listOfProducts>" << endl;
	cout << "<kineticLaw formula=\"\">" << endl;
	cout << "<listOfParameters>" << endl;
	cout << "<parameter value=\"" << rate << "\" name=\"c" << i << "\" />" << endl;
	cout << "</listOfParameters>" << endl;
	cout << "</kineticLaw>" << endl;
	cout << "</reaction>" << endl;
	cout << endl;
}

int main (int argc, char * const argv[]) 
{
	// Fisher Model
	
	int numberOfCells	= 100;
	double q			= 100.0;
	int numberOfSpecies = numberOfCells;
	double h			= 1.0;
	double nu			= 1.0;
	double k			= 1.0/q;
	
	// species
	cout << "<listOfSpecies>" << endl;
	// species 0...n-1
	for (int i = 0; i < numberOfSpecies; ++i)
	{
		double x = 0.0;
		if (i <  ( numberOfCells/2 )  ) {x = 1.0;}
		cout << "<species compartment=\"comp1\" initialAmount=\"" << x*q << "\" name=\"A" << i << "\"  units=\"mole\">" << endl;
		cout << "</species>" << endl << endl;
	}
	for (int i = 0; i < numberOfSpecies; ++i)
	{
		double x = 1.0;
		if (i < ( numberOfCells/2 ) ) {x = 0.0;}
		cout << "<species compartment=\"comp1\" initialAmount=\"" << x*q << "\" name=\"B" << i << "\"  units=\"mole\">" << endl;
		cout << "</species>" << endl << endl;
	}
	cout << "</listOfSpecies>" << endl << endl;
	
	// diffusion reactions
	string species = "A";
	int counter = 0;
	{
		cout << "<listOfReactions>" << endl;
		
		double rate			= nu / (h*h);
		
		for (int i = 1; i < numberOfSpecies-1; ++i)
		{
			// left
			printReaction( counter++, i, i+1, rate, species );
			// right
			printReaction( counter++, i, i-1, rate, species );
		}
		
		// boundaries
		// neumann on the left
		printReaction( counter++, 0, 1, rate, species );
		printReaction( counter++, 0, 0, rate, species );
		
		// neumann on the right
		printReaction( counter++, numberOfSpecies-1, numberOfSpecies-2, rate, species );
		printReaction( counter++, numberOfSpecies-1, numberOfSpecies-1, rate, species );
	}
	species = "B";
	{		
		double rate			= nu / (h*h);
		
		for (int i = 1; i < numberOfSpecies-1; ++i)
		{
			// left
			printReaction( counter++, i, i+1, rate, species );
			// right
			printReaction( counter++, i, i-1, rate, species );
		}
		
		// boundaries
		// neumann on the left
		printReaction( counter++, 0, 1, rate, species );
		printReaction( counter++, 0, 0, rate, species );
		
		// neumann on the right
		printReaction( counter++, numberOfSpecies-1, numberOfSpecies-2, rate, species );
		printReaction( counter++, numberOfSpecies-1, numberOfSpecies-1, rate, species );
	}
	
	
	for (int i = 1; i < numberOfSpecies; ++i)
	{
		printReactionFisher( counter++, i, i, k, "A", "B" );
	}
	
	
	cout << "</listOfReactions>" << endl;

}

/*
int main (int argc, char * const argv[]) 
{
	// Diffusion Model

	int numberOfCells	= 50;
	double q			= 1000000.0;
	int numberOfSpecies = numberOfCells;
	double h			= 5.0 / ( (double)numberOfCells );
	double nu			= 1e-2;
	double t			= 0.0;
	
	// species
	cout << "<listOfSpecies>" << endl;
	// species 0...n-1
	for (int i = 0; i < numberOfSpecies; ++i)
	{
		double x = 0.5*h + ( ((double)i)*h );
		cout << "<species compartment=\"comp1\" initialAmount=\"" << omega( x, nu, t, q ) << "\" name=\"S" << i << "\"  units=\"mole\">" << endl;
		cout << "</species>" << endl << endl;
	}
	cout << "</listOfSpecies>" << endl << endl;
	
	// reactions
	cout << "<listOfReactions>" << endl;
	
	double rate			= nu / (h*h);
	
	// bulk reactions, S_{1}...S_{n-2}
	int counter = 0;
	for (int i = 1; i < numberOfSpecies-1; ++i)
	{
		// left
		printReaction( counter++, i, i+1, rate );
		// right
		printReaction( counter++, i, i-1, rate );
	}
	
	// boundaries
	// dirichlet on the left
	printReaction( counter++, 0, 1, rate );
	printReaction( counter++, 0, 2.0*rate );
	
	// neumann on the right
	printReaction( counter++, numberOfSpecies-1, numberOfSpecies-2, rate );
	printReaction( counter++, numberOfSpecies-1, numberOfSpecies-1, rate );
	
	cout << "</listOfReactions>" << endl;
}
*/

/*
// "main" to generate the species an reactions for a SBML xml file	
int main (int argc, char * const argv[]) 
{
	// Lienar Chain Model

	int numberToGenerate	= 1000;
	bool species			= false;
	// species
	
	if (species == true)
	{		
		cout << "<species compartment=\"comp1\" initialAmount=\"10000.0\" name=\"S" << 1 << "\"  units=\"mole\">" << endl;
		cout << "</species>" << endl << endl;
			
		for (int i = 2; i <= numberToGenerate; ++i)
		{
			cout << "<species compartment=\"comp1\" initialAmount=\"0.0\" name=\"S" << i << "\"  units=\"mole\">" << endl;
			cout << "</species>" << endl << endl;
		}
	}
	else
	{
	// reactions
	
		for (int i = 1; i <= numberToGenerate-1; ++i)
		{
			cout << "<reaction name=\"R" << i << "\">" << endl;
			cout << "<listOfReactants>" << endl;
			cout << " <speciesReference species=\"S" << i << "\" />" << endl;
			cout << "</listOfReactants>" << endl;
			cout << "<listOfProducts>" << endl;
			cout << " <speciesReference species=\"S" << i+1 << "\" />" << endl;
			cout << "</listOfProducts>" << endl;
			cout << "<kineticLaw formula=\"\">" << endl;
			cout << "<listOfParameters>" << endl;
			cout << "<parameter value=\"" << 1 << ".0\" name=\"c" << i << "\" />" << endl;
			cout << "</listOfParameters>" << endl;
			cout << "</kineticLaw>" << endl;
			cout << "</reaction>" << endl;
			
			cout << endl;
		}
	}
}

*/

/*
{
	int numberToGenerate = 10;
	
	
	// species
	for (int i = 1; i <= numberToGenerate; ++i)
	{
		cout << "<species compartment=\"comp1\" initialAmount=\"0.0\" name=\"S" << i << "\"  units=\"mole\">" << endl;
		cout << "</species>" << endl << endl;
	}
	
	// reactions
	for (int i = 1; i <= numberToGenerate; ++i)
	{
		cout << "<reaction name=\"R" << i << "\">" << endl;
		cout << "<listOfReactants>" << endl;
		
		cout << "</listOfReactants>" << endl;
		cout << "<listOfProducts>" << endl;
		cout << " <speciesReference species=\"S" << i << "\" />" << endl;
		cout << "</listOfProducts>" << endl;
		cout << "<kineticLaw formula=\"\">" << endl;
		cout << "<listOfParameters>" << endl;
		
		if ( i != numberToGenerate)
			cout << "<parameter value=\"" << 1 << ".0\" name=\"c" << i << "\" />" << endl;
		else
			cout << "<parameter value=\"" << 100 << ".0\" name=\"c" << i << "\" />" << endl;
		cout << "</listOfParameters>" << endl;
		cout << "</kineticLaw>" << endl;
		cout << "</reaction>" << endl;
		
		cout << endl;
	}
}
*/
