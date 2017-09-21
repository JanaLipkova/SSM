/*
 *  SBMLReaderAndParser.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/2/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "SBMLReaderAndParser.h"

SBMLReaderAndParser::SBMLReaderAndParser(string filename)
{
	this->filename = filename;

	timeStart[0] = "<stochSim:TimeStart>";
	timeStart[1] = "</stochSim:TimeStart>";

	timeEnd  [0] = "<stochSim:TimeEnd>";
	timeEnd  [1] = "</stochSim:TimeEnd>";

	storeInterval[0] = "<stochSim:StoreInterval>";
	storeInterval[1] = "</stochSim:StoreInterval>";

  epsilon[0]       = "<stochSim:Epsilon>";
	epsilon[1]       = "</stochSim:Epsilon>";

  numberOfSamples[0] = "<stochSim:NumberOfSamples>";
	numberOfSamples[1] = "</stochSim:NumberOfSamples>";

	method[0] = "<stochSim:Method>";
	method[1] = "</stochSim:Method>";

	theta[0] = "<stochSim:Theta>";
	theta[1] = "</stochSim:Theta>";

	sortInterval[0] = "<stochSim:SortInterval>";
	sortInterval[1] = "</stochSim:SortInterval>";

	initialNoise[0]			= "<stochSim:InitialNoise>";
	initialNoise[1]			= "</stochSim:InitialNoise>";

	noiseIncrement[0]		= "<stochSim:NoiseIncrement>";
	noiseIncrement[1]		= "</stochSim:NoiseIncrement>";

	numberOfNoiseLevels[0]	= "<stochSim:NumberOfNoiseLevels>";
	numberOfNoiseLevels[1]	= "</stochSim:NumberOfNoiseLevels>";
}

SBMLReaderAndParser::~SBMLReaderAndParser()
{
}

string SBMLReaderAndParser::findSubstringAnnotation(string annotation, string s [2])
{
	int f1, f2;
	f1 = annotation.find(s[0]) + s[0].length();
	f2 = annotation.find(s[1]);
	string ret = annotation.substr(f1, f2-f1);
	boost::algorithm::trim(ret);
	return ret;
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
	cout << "T start:				 "	<<  findSubstringAnnotation(annotation, timeStart)		<< endl;
	cout << "T final:				 "	<<  findSubstringAnnotation(annotation, timeEnd)		<< endl;
	cout << "Simulation method:		 "	<<  findSubstringAnnotation(annotation, method)				<< endl;
	cout << "Number of samples:		 "	<<  findSubstringAnnotation(annotation, numberOfSamples)	<< endl;
	cout << endl;
	cout << "Store interval:		 "	<<  findSubstringAnnotation(annotation, storeInterval)	<< endl;
	cout << "Epsilon:				 "	<<  findSubstringAnnotation(annotation, epsilon)		<< endl;
	cout << "Theta:					 "	<<  findSubstringAnnotation(annotation, theta)			<< endl;
	cout << "Sort interval:			 "	<<  findSubstringAnnotation(annotation, sortInterval)		<< endl;
	cout << "Initial noise:          "	<<  findSubstringAnnotation(annotation, initialNoise)		<< endl;
	cout << "Noise increment:        "	<<  findSubstringAnnotation(annotation, noiseIncrement)		<< endl;
	cout << "Number of noise levels: "	<<  findSubstringAnnotation(annotation, numberOfNoiseLevels)		<< endl;

	cout << endl;

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
		cout << species->getName() << ":	" << species->getInitialAmount() << "	" << species->getId() << endl;
		//cout << species->getName() << ":	" << species->getInitialAmount() << "	" << i << endl;
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
