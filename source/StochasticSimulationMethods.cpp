/*
 *  StochasticSimulationMethods.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

// C/C++, Blitz++, RanLib, SBML, and Boost

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
#include "Methods/SSALDM.h"

// R-Leaping
#include "Methods/RLeaping.h"
#include "Methods/RLeapingJana.h"

// S-Leaping and co.
#include "Methods/SLeaping.h"

#include "Methods/SLeaping_v3.h"
#include "Methods/SLeaping_v4.h"
#include "Methods/SLeaping_v5.h"


#include "Methods/AdaptiveSLeaping.h"
//#include "Methods/AdaptiveSLeapingCL.h"
//#include "Methods/SLeapingNonNegative.h"
//#include "Methods/SLeapingRcontrol.h"

// Tau-Leaping and co.
#include "Methods/TauLeaping.h"
#include "Methods/TauLeapingCT.h"
#include "Methods/AdaptiveTau.h"
#include "Methods/TauLeapingNonNegative.h"

// Dealy methods
#include "Methods/DelaySSA.h"
#include "Methods/DelayRLeaping.h"
#include "Methods/DelayTauLeaping.h"
#include "Methods/DelayTauLeapingLeier.h"
#include "Methods/DelayTauLeapingPhilippe.h"

//// LacZ / LacY system
#include "Methods/SSA_LacZLacY.h"
//#include "Methods/RLeaping_LacZLacY.h"
//#include "Methods/TauLeaping_LacZLacY.h"
//#include "Methods/SLeaping_LacZLacY.h"

//#include "Methods/RLeapingGPU.h"


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

	// next 2 arguments are RNG seeds
	if (argc >= 4)
	{
		setall( (long int) atof(argv[2]), (long int) atof(argv[3]) ); // RNG seeds
		cout << "Random Number Generator Seeds Set." << endl;
		cout << "Seed 1: " << (long int) atof(argv[2]) << endl;
		cout << "Seed 2: " << (long int) atof(argv[3]) << endl << endl;
	}
	else // take generic seeds
	{
		setall(1056577, 112025); // RNG seeds
	}



	// Parse SBML file
	SBMLReaderAndParser * sBMLReaderAndParser	= new SBMLReaderAndParser(filename);
	sBMLReaderAndParser->readAndParse();
	SBMLDocument * sbmlDocument					= sBMLReaderAndParser->getSBMLDocument();

	//cout << "Are you content with the parsing of the SBML file? (y/n)" << endl;
	//string inputString;
	//cin  >> inputString;
	//string toSimulateOrNotToSimulate = boost::algorithm::trim_copy( inputString );

//	if ( (toSimulateOrNotToSimulate != "y") && (toSimulateOrNotToSimulate != "Y") && (toSimulateOrNotToSimulate != "yes") && (toSimulateOrNotToSimulate != "YES") )
//	{
//		cout << "Since you are discontent, the process will terminate.  Do review your SBML input file." << endl;
//		delete sBMLReaderAndParser;
//		return EXIT_SUCCESS;
//	}
	delete sBMLReaderAndParser;



	Simulation	* simulation	= new Simulation(sbmlDocument);
	Method		* method		= NULL;

	string selectedMethod = boost::algorithm::trim_copy( simulation->StochasticSimulationMethod );

	if ( ( selectedMethod == "SSA")				|| ( selectedMethod == "ESSA")				|| ( selectedMethod == "StochasticSimulationAlgorithm") )
	{ method	= new SSA		(simulation); }
    //else if ( ( selectedMethod == "SSALDM")		|| ( selectedMethod == "SSA-LDM")			|| ( selectedMethod == "SSA-LogDM")		|| ( selectedMethod == "SSA-LogDirectMethod") )
    //{ method	= new SSALDM	(simulation); }
	//else if ( ( selectedMethod == "DSSA")		|| ( selectedMethod == "DelaySSA")			|| ( selectedMethod == "D-SSA")			|| ( selectedMethod == "Delay-SSA") )
	//{ method	= new DelaySSA	(simulation); }
	//else if ( ( selectedMethod == "RLeap")		|| ( selectedMethod == "RLeaping")			|| ( selectedMethod == "R-Leap")		|| ( selectedMethod == "R-Leaping") )
	//{ method	= new RLeaping	(simulation); }
	// else if ( ( selectedMethod == "TauLeap")	|| ( selectedMethod == "TauLeaping")		|| ( selectedMethod == "Tau-Leap")		|| ( selectedMethod == "Tau-Leaping") )
	// { method	= new TauLeaping (simulation); }
	//else if ( ( selectedMethod == "DRLeap")		|| ( selectedMethod == "DelayRLeaping")		|| ( selectedMethod == "DR-Leap")		|| ( selectedMethod == "DR-Leaping") )
	// { method	= new DelayRLeaping (simulation); }
	// else if ( ( selectedMethod == "DTauLeap")	|| ( selectedMethod == "DelayTauLeaping")	|| ( selectedMethod == "DTau-Leap")		|| ( selectedMethod == "DTau-Leaping") )
	// { method	= new DelayTauLeapingPhilippe (simulation); }
	// else if ( ( selectedMethod == "DTauLeapLeier")	|| ( selectedMethod == "DelayTauLeapingLeier")	|| ( selectedMethod == "DTau-LeapLeier")		|| ( selectedMethod == "DTau-LeapingLeier") )
	// { method	= new DelayTauLeapingLeier (simulation); }
	// else if ( ( selectedMethod == "DRLeap")		|| ( selectedMethod == "DelayRLeaping")		|| ( selectedMethod == "DR-Leap")		|| ( selectedMethod == "DR-Leaping") )
	// { method	= new DelayRLeaping (simulation); }
	// else if ( ( selectedMethod == "TauLeapCT")	|| ( selectedMethod == "TauLeapingCT")		|| ( selectedMethod == "Tau-LeapCT")		|| ( selectedMethod == "Tau-LeapingCT") )
	// { method	= new TauLeapingCT (simulation); }
	// else if ( ( selectedMethod == "RLeapJana")		|| ( selectedMethod == "RLeapingJana")			|| ( selectedMethod == "R-LeapJana")		|| ( selectedMethod == "R-LeapingJana") )
	// { method	= new RLeapingJana (simulation); }
	//  else if ( ( selectedMethod == "SLeap")		|| ( selectedMethod == "SLeaping")			|| ( selectedMethod == "S-Leap")		|| ( selectedMethod == "S-Leaping") )
	//  { method	= new SLeaping (simulation); }
     else if ( ( selectedMethod == "SLeap_v3")		|| ( selectedMethod == "SLeaping_v3")			|| ( selectedMethod == "S-Leap_v3")		|| ( selectedMethod == "S-Leaping_v3") )
     { method	= new SLeaping_v3 (simulation); }
     else if ( ( selectedMethod == "SLeap_v4")		|| ( selectedMethod == "SLeaping_v4")			|| ( selectedMethod == "S-Leap_v4")		|| ( selectedMethod == "S-Leaping_v4") )
     { method	= new SLeaping_v4 (simulation); }
     else if ( ( selectedMethod == "SLeap_v5")		|| ( selectedMethod == "SLeaping_v5")			|| ( selectedMethod == "S-Leap_v5")		|| ( selectedMethod == "S-Leaping_v5") )
     { method	= new SLeaping_v5 (simulation); }
	else if ( ( selectedMethod == "AdaptiveTauLeap")|| ( selectedMethod == "AdaptiveTauLeaping")|| ( selectedMethod == "AdaptiveTau")|| ( selectedMethod == "AdatpiveTau-Leaping") )
	{ method	= new AdaptiveTau (simulation); }
	else if ( ( selectedMethod == "AdaptiveSLeap")|| ( selectedMethod == "AdaptiveSLeaping")|| ( selectedMethod == "AdaptiveS")|| ( selectedMethod == "AdatpiveS-Leaping") )
	{ method	= new AdaptiveSLeaping (simulation); }
	// //else if ( ( selectedMethod == "AdaptiveSLeapCL")|| ( selectedMethod == "AdaptiveSLeapingCL")|| ( selectedMethod == "AdaptiveSCL")|| ( selectedMethod == "AdatpiveS-LeapingCL") )
	// //{ method	= new AdaptiveSLeapingCL (simulation); }
	// else if ( ( selectedMethod == "TauLeapingNonNegative")|| ( selectedMethod == "TauLeapNonNegative")|| ( selectedMethod == "TauLeapingNN")|| ( selectedMethod == "TauLeaping-NN") )
	// { method	= new TauLeapingNonNegative (simulation); }
	//else if ( ( selectedMethod == "SSALacZLacY")		|| ( selectedMethod == "SSA-Lac")			|| ( selectedMethod == "SSALAC")			|| ( selectedMethod == "SSALac") )
	// { method	= new SSA_LacZLacY	(simulation); }
	//else if ( ( selectedMethod == "RLeapingLacZLacY")	|| ( selectedMethod == "RLeaping-Lac")		|| ( selectedMethod == "RLeapingLAC")		|| ( selectedMethod == "RLeapingLac") )
	//{ method	= new RLeaping_LacZLacY	(simulation); }
	// else if ( ( selectedMethod == "TauLeapingLacZLacY")	|| ( selectedMethod == "TauLeaping-Lac")		|| ( selectedMethod == "TauLeapingLAC")		|| ( selectedMethod == "TauLeapingLac") )
	// { method	= new TauLeaping_LacZLacY	(simulation); }
	//else if ( ( selectedMethod == "SLeapingLacZLacY")	|| ( selectedMethod == "SLeaping-Lac")		|| ( selectedMethod == "SLeapingLAC")		|| ( selectedMethod == "SLeapingLac") )
	//{ method	= new SLeaping_LacZLacY	(simulation); }

//	else if ( ( selectedMethod == "RLeapGPU") || ( selectedMethod == "RLeapingGPU") )
//	{ method	= new RLeapingGPU(simulation); }
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
