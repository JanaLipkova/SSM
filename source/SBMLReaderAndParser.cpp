/*
 *  SBMLReaderAndParser.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/2/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 *
 */

#include "SBMLReaderAndParser.h"

SBMLReaderAndParser::SBMLReaderAndParser(string filename)
{
	this->filename = filename;

	timeStartL1[0] = "<stochSim:TimeStart>";
	timeStartL1[1] = "</stochSim:TimeStart>";

	timeEndL1  [0] = "<stochSim:TimeEnd>";
	timeEndL1  [1] = "</stochSim:TimeEnd>";

	storeIntervalL1[0] = "<stochSim:StoreInterval>";
	storeIntervalL1[1] = "</stochSim:StoreInterval>";

        epsilonL1[0]       = "<stochSim:Epsilon>";
	epsilonL1[1]       = "</stochSim:Epsilon>";

	deltaL1[0]         = "<stochSim:Delta>";
        deltaL1[1]         = "</stochSim:Delta>";

        numberOfSamplesL1[0] = "<stochSim:NumberOfSamples>";
	numberOfSamplesL1[1] = "</stochSim:NumberOfSamples>";

	methodL1[0] = "<stochSim:Method>";
	methodL1[1] = "</stochSim:Method>";

	thetaL1[0] = "<stochSim:Theta>";
	thetaL1[1] = "</stochSim:Theta>";

	sortIntervalL1[0] = "<stochSim:SortInterval>";
	sortIntervalL1[1] = "</stochSim:SortInterval>";

	initialNoiseL1[0]			= "<stochSim:InitialNoise>";
	initialNoiseL1[1]			= "</stochSim:InitialNoise>";

	noiseIncrementL1[0]		= "<stochSim:NoiseIncrement>";
	noiseIncrementL1[1]		= "</stochSim:NoiseIncrement>";

	numberOfNoiseLevelsL1[0]	= "<stochSim:NumberOfNoiseLevels>";
	numberOfNoiseLevelsL1[1]	= "</stochSim:NumberOfNoiseLevels>";


	timeStartL2 			= "stochSim:TimeStart=\"";
	timeEndL2             	= "stochSim:TimeEnd=\"";
	storeIntervalL2       	= "stochSim:StoreInterval=\"";
	epsilonL2             = "stochSim:Epsilon=\"";
	deltaL2               = "stochSim:Delta=\"";
	numberOfSamplesL2     = "stochSim:NumberOfSamples=\"";
	methodL2              = "stochSim:Method=\"";
	thetaL2               = "stochSim:Theta=\"";
	sortIntervalL2				= "stochSim:SortInterval=\"";
	initialNoiseL2			  = "stochSim:InitialNoise=\"";
	noiseIncrementL2			= "stochSim:NoiseIncrement=\"";
	numberOfNoiseLevelsL2 = "stochSim:NumberOfNoiseLevels=\"";

}




SBMLReaderAndParser::~SBMLReaderAndParser()
{
}




string SBMLReaderAndParser::findSubstringAnnotationL1(string annotation, string s[2])
{
	int f1, f2;
	f1 = annotation.find(s[0]) + s[0].length();
	f2 = annotation.find(s[1]);
	string ret = annotation.substr(f1, f2-f1);
	boost::algorithm::trim(ret);
	return ret;
}



string SBMLReaderAndParser::findSubstringAnnotationL2(string annotation, string s)
{
	int f1;
	f1 = annotation.find(s) + s.length();
	string t (annotation,f1,annotation.length()-f1);
	f1 = t.find("\"");
	string r = t.substr(0, f1);
	boost::algorithm::trim(r);

	return r;
}







/**
 * readAndParse
 * collects information in SBML
 * and outputs it to cout
 *
 * User can validate output to see whether the SBML file could be parsed correctly
 *
 * TODO no. of compartments currently supported = 1
 */
int SBMLReaderAndParser::readAndParse()
{
	this->document = readSBML(filename.c_str());

	cout << "-----------" << endl;
	cout << "SBML: Level:"   << document->getLevel();
	cout << "  Version:" << document->getVersion() << endl;


	unsigned int errors = document->getNumErrors();
	if (errors > 0)
	{
		document->printErrors(cerr);
		return errors;
	}
	cout << endl;
	cout << "Filename: " << filename << endl;
	cout << "Error(s): " << errors  << endl;
	cout << endl;

	Model * model = document->getModel();
	cout << "Model name: " << model->getName() << endl << endl;



	// Stochastic Simulation Parameters are in the "Annotation" Field
	// string annotation = model->getAnnotationString(); //libSBML-4.0.1
	string annotation = model->getAnnotationString();


	if( document->getLevel() == 1){
		cout << "T start:                 "	<<  findSubstringAnnotationL1(annotation, timeStartL1)		<< endl;
		cout << "T final:                 "	<<  findSubstringAnnotationL1(annotation, timeEndL1)		<< endl;
		cout << "Simulation method:       "	<<  findSubstringAnnotationL1(annotation, methodL1)				<< endl;
		cout << "Number of samples:       "	<<  findSubstringAnnotationL1(annotation, numberOfSamplesL1)	<< endl;
		cout << endl;
		cout << "Store interval:          "	<<  findSubstringAnnotationL1(annotation, storeIntervalL1)	<< endl;
		cout << "Epsilon:                 "	<<  findSubstringAnnotationL1(annotation, epsilonL1)		<< endl;
		cout << "Delta:                   "     <<  findSubstringAnnotationL1(annotation, deltaL1)            << endl;
		cout << "Theta:                   "	<<  findSubstringAnnotationL1(annotation, thetaL1)			<< endl;
		cout << "Sort interval:           "	<<  findSubstringAnnotationL1(annotation, sortIntervalL1)		<< endl;
		cout << "Initial noise:           "	<<  findSubstringAnnotationL1(annotation, initialNoiseL1)		<< endl;
		cout << "Noise increment:         "	<<  findSubstringAnnotationL1(annotation, noiseIncrementL1)		<< endl;
		cout << "Number of noise levels:  "	<<  findSubstringAnnotationL1(annotation, numberOfNoiseLevelsL1)		<< endl;
		cout << endl;
	}
	else if( document->getLevel() == 2){
			cout << "T start:                 "	<<  findSubstringAnnotationL2(annotation, timeStartL2)		<< endl;
			cout << "T final:                 "	<<  findSubstringAnnotationL2(annotation, timeEndL2)		<< endl;
			cout << "Simulation method:       "	<<  findSubstringAnnotationL2(annotation, methodL2)				<< endl;
			cout << "Number of samples:       "	<<  findSubstringAnnotationL2(annotation, numberOfSamplesL2)	<< endl;
			cout << endl;
			cout << "Store interval:          "	<<  findSubstringAnnotationL2(annotation, storeIntervalL2)	<< endl;
			cout << "Epsilon:                 "	<<  findSubstringAnnotationL2(annotation, epsilonL2)		<< endl;
			cout << "Delta:                   "     <<  findSubstringAnnotationL2(annotation, deltaL2)            << endl;
			cout << "Theta:                   "	<<  findSubstringAnnotationL2(annotation, thetaL2)			<< endl;
			cout << "Sort interval:           "	<<  findSubstringAnnotationL2(annotation, sortIntervalL2)		<< endl;
			cout << "Initial noise:           "	<<  findSubstringAnnotationL2(annotation, initialNoiseL2)		<< endl;
			cout << "Noise increment:         "	<<  findSubstringAnnotationL2(annotation, noiseIncrementL2)		<< endl;
			cout << "Number of noise levels:  "	<<  findSubstringAnnotationL2(annotation, numberOfNoiseLevelsL2)		<< endl;
			cout << endl;
	}
	else{
		cout << "We work only with Level 1 and Level 2 SBML files. Exiting..." << endl ;
		exit(EXIT_FAILURE);
	}


cout << "-----------" << endl;



	//cout << "Number of compartments: " << model->getNumCompartments() << endl;
	cout << "Number of species  : " << model->getNumSpecies() << endl;
	cout << "Number of reactions: " << model->getNumReactions() << endl;
	cout << endl << endl;

	Species * species;
	cout << "	Species i: initial amount" << endl;
	for (int i = 0; i < model->getNumSpecies(); ++i)
	{
		cout << "			";
		species = model->getSpecies(i);
		cout << species->getName() << ":	" << species->getInitialAmount() << "	" << species->getName() << endl;
	}
	cout << endl << endl;



	Reaction * reaction;
	for (int i = 0; i < model->getNumReactions(); ++i)
	{
		reaction = model->getReaction(i);
		KineticLaw * kineticLaw = reaction->getKineticLaw();
		Parameter * parameter = kineticLaw->getParameter(0);
		double rate = parameter->getValue();

		cout << "	Reaction i: " << i << "	with rate: " << rate << endl;

		if (kineticLaw->getNumParameters() > 1)
		{
			cout << "		Note: reaction has more than 1 parameter => delayed." << endl;
			cout << "		delay time:			" << kineticLaw->getParameter(1)->getValue() << endl;
			cout << "		Hill coefficient:		" << kineticLaw->getParameter(2)->getValue() << endl;
			cout << "		P_0:					" << kineticLaw->getParameter(3)->getValue() << endl;
		}

		if (kineticLaw->getFormula().substr(0, 8) == "Species:")
		{ cout << "		Note: rate is dependent on: " << kineticLaw->getFormula() << endl;  } //cout << kineticLaw->getFormula().substr(8, 1) << endl;

		cout << "		Number of reactants: " << reaction->getNumReactants() << endl;
		cout << "		Number of products : " << reaction->getNumProducts() << endl;

		cout << "			";

		for (int j = 0; j < reaction->getNumReactants(); ++j)
		{
			SpeciesReference * speciesReference = reaction->getReactant(j);
			cout << speciesReference->getStoichiometry() << "  ";
			cout << speciesReference->getSpecies();
			if (j != reaction->getNumReactants()-1)
				cout << "	+	";
			else
				cout << "		->		";
		}
		if (reaction->getNumReactants() == 0)
		{
			cout << "0	" << "	->	   ";
		}

		for (int j = 0; j < reaction->getNumProducts(); ++j)
		{
			SpeciesReference * speciesReference = reaction->getProduct(j);
			cout << speciesReference->getStoichiometry() << "  ";
			cout << speciesReference->getSpecies();
			if (j != reaction->getNumProducts()-1)
				cout << "	+	";
		}
		if (reaction->getNumProducts() == 0)
		{
			cout << "0	";
		}

		cout << endl;
		cout << endl;

	}


	return EXIT_SUCCESS;
}

SBMLDocument * SBMLReaderAndParser::getSBMLDocument()
{
	return this->document;
}
