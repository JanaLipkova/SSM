/*
 *  StochasticSimulationMethods.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

// C/C++, Blitz++, SBML, and Boost


#include "HeaderFiles.h"

// Program Files
#include "SBMLReaderAndParser.h"
#include "Simulation.h"

#define PARTICLEDEFINED
#include "Methods/Method.h"

#include "Timer/Timer.h"

// Simulation Methods

// SSA
#include "Methods/SSA.h"

// R-Leaping
#include "Methods/RLeaping.h"

// S-Leaping and co.
#include "Methods/SLeaping.h"
#include "Methods/AdaptiveSLeaping.h"

// Tau-Leaping and co.
#include "Methods/TauLeaping.h"
#include "Methods/AdaptiveTau.h"


#include "my_rand.h"


int main (int argc, char * const argv[])
{
	// Check args for input file name
	// You must supply a SBML file
	if (argc < 2)
	{
		cout << endl << "Oh blimey.  Evidently something has gone askew old chap." << endl;
		cout << "Do check that you have supplied an input file via the command line." << endl;
		return EXIT_FAILURE;
	}

	string filename = argv[1];

	// next argument are RNG seeds
	if (argc >= 3){
		cout << "Random Number Generator Seeds Set." << endl;
		cout << "Seed : " << (long int) atof(argv[2]) << endl;
		myrand::engine.seed(atof(argv[2]));
	}
	else{
		myrand::engine.seed(1234);
	}



	// Parse SBML file
	SBMLReaderAndParser * sBMLReaderAndParser	= new SBMLReaderAndParser(filename);
	sBMLReaderAndParser->readAndParse();
	SBMLDocument * sbmlDocument					= sBMLReaderAndParser->getSBMLDocument();
	delete sBMLReaderAndParser;



	Simulation	* simulation	= new Simulation(sbmlDocument);
	Method		* method		= NULL;

	string selectedMethod = boost::algorithm::trim_copy( simulation->StochasticSimulationMethod );

	if ( ( selectedMethod == "SSA")	 || ( selectedMethod == "ESSA")	|| ( selectedMethod == "StochasticSimulationAlgorithm") )
		method	= new SSA(simulation);

 	else if ( ( selectedMethod == "TauLeap") || ( selectedMethod == "TauLeaping") || ( selectedMethod == "Tau-Leap") || ( selectedMethod == "Tau-Leaping") )
		method	= new TauLeaping (simulation);

	else if ( ( selectedMethod == "AdaptiveTauLeap") || ( selectedMethod == "AdaptiveTauLeaping") || ( selectedMethod == "AdaptiveTau")|| ( selectedMethod == "AdatpiveTau-Leaping") )
		method	= new AdaptiveTau (simulation);

	else if ( ( selectedMethod == "RLeap") || ( selectedMethod == "RLeaping") || ( selectedMethod == "R-Leap") || ( selectedMethod == "R-Leaping") )
		method	= new RLeaping (simulation);

    else if ( ( selectedMethod == "SLeap") || ( selectedMethod == "SLeaping") || ( selectedMethod == "S-Leap") || ( selectedMethod == "S-Leaping") )
    	method	= new SLeaping (simulation);

	else if ( ( selectedMethod == "AdaptiveSLeap") || ( selectedMethod == "AdaptiveSLeaping") || ( selectedMethod == "AdaptiveS")|| ( selectedMethod == "AdatpiveS-Leaping") )
		method	= new AdaptiveSLeaping (simulation);

	else
	{
		cout << "Oh blimey.  The simulation method that you requested, \"" << selectedMethod << "\", is not yet supported." << endl;
	}

	if ( method != NULL) {
		Timer timer;
		timer.StartSW();
		method->solve();
		timer.StopSW();
		cout << "Running time: " << timer.ReadSW() << endl;
	}

	delete method;
	delete simulation;

	delete sbmlDocument;

	return EXIT_SUCCESS;
}
