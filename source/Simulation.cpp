/*
 *  Simulation.cpp
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#include "Simulation.h"

Simulation::Simulation(SBMLDocument * sbmlDocument)
{
	this->sbmlDocument	= sbmlDocument;
	this->sbmlModel		= sbmlDocument->getModel();

	ModelName			= sbmlModel->getName();

	speciesValues.resize(sbmlModel->getNumSpecies());
	proposedSpeciesValues.resize(sbmlModel->getNumSpecies());
	proposedSpeciesValues = (ParticleType)0.0;
	loadInitialConditions();


	if( sbmlModel->getLevel()==1 ){
			string timeStartL1			[2];
			string timeEndL1			[2];
			string storeIntervalL1		[2];
			string epsilonL1			[2];
			string numberOfSamplesL1	[2];
			string methodL1				[2];
			string thetaL1				[2];
			string sortIntervalL1		[2];
			string initialNoiseL1		[2];
			string noiseIncrementL1		[2];
			string numberOfNoiseLevelsL1[2];
			string tauL1				[2];
			string deltaL1				[2];


			timeStartL1[0]			= "<stochSim:TimeStart>";
			timeStartL1[1]			= "</stochSim:TimeStart>";

			timeEndL1[0]			= "<stochSim:TimeEnd>";
			timeEndL1[1]			= "</stochSim:TimeEnd>";

		    storeIntervalL1[0]		= "<stochSim:StoreInterval>";
			storeIntervalL1[1]		= "</stochSim:StoreInterval>";

		    epsilonL1[0]				= "<stochSim:Epsilon>";
			epsilonL1[1]				= "</stochSim:Epsilon>";

		    numberOfSamplesL1[0]		= "<stochSim:NumberOfSamples>";
			numberOfSamplesL1[1]		= "</stochSim:NumberOfSamples>";

			methodL1[0]				= "<stochSim:Method>";
			methodL1[1]				= "</stochSim:Method>";

			thetaL1[0]				= "<stochSim:Theta>";
			thetaL1[1]				= "</stochSim:Theta>";

			sortIntervalL1[0]			= "<stochSim:SortInterval>";
			sortIntervalL1[1]			= "</stochSim:SortInterval>";

			initialNoiseL1[0]			= "<stochSim:InitialNoise>";
			initialNoiseL1[1]			= "</stochSim:InitialNoise>";

			noiseIncrementL1[0]		= "<stochSim:NoiseIncrement>";
			noiseIncrementL1[1]		= "</stochSim:NoiseIncrement>";

			numberOfNoiseLevelsL1[0]	= "<stochSim:NumberOfNoiseLevels>";
			numberOfNoiseLevelsL1[1]	= "</stochSim:NumberOfNoiseLevels>";

			tauL1[0]					= "<stochSim:Tau>";
			tauL1[1]					= "</stochSim:Tau>";

			deltaL1[0]				= "<stochSim:Delta>";
			deltaL1[1]				= "</stochSim:Delta>";

			//string annotation = sbmlModel->getAnnotationString(); //libSBML-4.0.1	StartTime		= dFindSubstringAnnotation(annotation,		timeStart);
			string annotation = sbmlModel->getAnnotationString();

			StartTime = dFindSubstringAnnotationL1(annotation,	timeStartL1);

			EndTime			= dFindSubstringAnnotationL1(annotation,		timeEndL1);
			StoreInterval   = (int)dFindSubstringAnnotationL1(annotation,	storeIntervalL1);
			Epsilon			= dFindSubstringAnnotationL1(annotation,		epsilonL1);
			NumberOfSamples = (int)dFindSubstringAnnotationL1(annotation,	numberOfSamplesL1);

			StochasticSimulationMethod = sFindSubstringAnnotationL1(annotation, methodL1);

			SortInterval	= (int)dFindSubstringAnnotationL1(annotation,	sortIntervalL1);
			Theta			= dFindSubstringAnnotationL1(annotation,	thetaL1);

			InitialNoise= dFindSubstringAnnotationL1(annotation,	initialNoiseL1);
			NoiseIncrement = dFindSubstringAnnotationL1(annotation,	noiseIncrementL1);
			NumberOfNoiseLevels= (int)dFindSubstringAnnotationL1(annotation,	numberOfNoiseLevelsL1);

			Tau		= dFindSubstringAnnotationL1(annotation, tauL1);
			Delta	= dFindSubstringAnnotationL1(annotation, deltaL1);
}
else if( sbmlModel->getLevel()==2 ){


		string timeStartL2				= "stochSim:TimeStart=\"";
		string timeEndL2             	= "stochSim:TimeEnd=\"";
		string storeIntervalL2       	= "stochSim:StoreInterval=\"";
		string epsilonL2             	= "stochSim:Epsilon=\"";
		string numberOfSamplesL2     	= "stochSim:NumberOfSamples=\"";
		string methodL2              	= "stochSim:Method=\"";
		string thetaL2               	= "stochSim:Theta=\"";
		string sortIntervalL2			= "stochSim:SortInterval=\"";
		string initialNoiseL2			= "stochSim:InitialNoise=\"";
		string noiseIncrementL2			= "stochSim:NoiseIncrement=\"";
		string numberOfNoiseLevelsL2	= "stochSim:NumberOfNoiseLevels=\"";
		string tauL2					= "stochSim:Tau";
		string deltaL2					= "stochSim:Delta";


		string annotation = sbmlModel->getAnnotationString();

		StartTime = dFindSubstringAnnotationL2(annotation,		timeStartL2);

		EndTime			= 		dFindSubstringAnnotationL2(annotation,	timeEndL2);
		StoreInterval   = (int)	dFindSubstringAnnotationL2(annotation,	storeIntervalL2);
		Epsilon			= 		dFindSubstringAnnotationL2(annotation,	epsilonL2);
		NumberOfSamples = (int)	dFindSubstringAnnotationL2(annotation,	numberOfSamplesL2);

		StochasticSimulationMethod = sFindSubstringAnnotationL2(annotation, methodL2);

		SortInterval	= (int)	dFindSubstringAnnotationL2(annotation,	sortIntervalL2);
		Theta			= 		dFindSubstringAnnotationL2(annotation,	thetaL2);

		InitialNoise		= 		dFindSubstringAnnotationL2(annotation,	initialNoiseL2);
		NoiseIncrement 		= 		dFindSubstringAnnotationL2(annotation,	noiseIncrementL2);
		NumberOfNoiseLevels	= (int)	dFindSubstringAnnotationL2(annotation,	numberOfNoiseLevelsL2);

		Tau		= dFindSubstringAnnotationL2(annotation, tauL2);
		Delta	= dFindSubstringAnnotationL2(annotation, deltaL2);

}

	speciesEnsemble.resize(sbmlModel->getNumSpecies(), StoreInterval+1);
	speciesEnsemble = 0.0;

	// create SSMReaction objects
	for (int ir = 0; ir < sbmlModel->getNumReactions(); ++ir)
	{
		SSMReaction * ssmr = new SSMReaction();


		Reaction * sbmlReaction = sbmlModel->getReaction(ir);

		// save the reactants
		for (int j = 0; j < sbmlReaction->getNumReactants(); ++j)
		{
			SpeciesReference * speciesReference = sbmlReaction->getReactant(j);
			string speciesName = speciesReference->getSpecies();
			int speciesIndex = getSpeciesIndex(speciesName);
			ssmr->addReactant(speciesIndex);

			double stochiometricCoefficient = speciesReference->getStoichiometry();
			ssmr->setLastReactantNu((int)stochiometricCoefficient);

			//cout << j << "--" << speciesName << "---" << speciesIndex << endl ;

		}

		// save the products
		for (int j = 0; j < sbmlReaction->getNumProducts(); ++j)
		{
			SpeciesReference * speciesReference = sbmlReaction->getProduct(j);
			string speciesName = speciesReference->getSpecies();
			int speciesIndex = getSpeciesIndex(speciesName);
			ssmr->addProduct(speciesIndex);

			double stochiometricCoefficient = speciesReference->getStoichiometry();
			ssmr->setLastProductNu((int)stochiometricCoefficient);
		}

		KineticLaw * kineticLaw = sbmlReaction->getKineticLaw();

		Parameter * parameter = kineticLaw->getParameter(0);

		// save the rate
		ssmr->setRate( parameter->getValue() );


		ssmr->finalizeReaction();

		// for testing
		//ssmr->toString();
		ssmReactionList.push_back(ssmr);
	}


	// s	species
	// ir	reaction
	// r	reactants
	for (unsigned int s = 0; s < sbmlModel->getNumSpecies(); ++s)
	{
		// allocate a pointer to hold the reactions
		backPointers.push_back( vector<int>() );

		// find the reactions in which a particular species is involved
		for (unsigned int ir = 0; ir < sbmlModel->getNumReactions(); ++ir)
		{
			bool isInvolved = false;

			SSMReaction * reaction		= ssmReactionList[ir];
			vector <int>  reactants		= reaction->getReactants();

			// loop through the reactants of reaction "ir"
			for (unsigned int r = 0; r < reactants.size(); ++r)
			{
				if ( s == reactants[r] )
				{
					isInvolved = true;
					break;
				}
			}

			if (isInvolved == true)
			{
				backPointers[s].push_back(ir);
			}
		}

	}



	// for testing...
	/*
	cout << "Species back-pointers to reactions" << endl;
	for (unsigned int s = 0; s < sbmlModel->getNumSpecies(); ++s)
	{
		vector <int> & temp = backPointers[s];
		for (unsigned int j = 0; j < temp.size(); ++j)
		{
			cout << temp[j] << " ";
		}
		cout << endl;
	}
	exit(0);
	*/
}


Simulation::Simulation()
{
}


Simulation::~Simulation()
{
	for (int ir = 0; ir < ssmReactionList.size(); ++ir) { delete ssmReactionList[ir]; }
}


Model * Simulation::getSBMLModel()
{
	return this->sbmlModel;
}


int Simulation::getSpeciesIndex(string speciesName)
{
	iter = speciesData.begin();
	iter = speciesData.find(speciesName);
	return iter->second;
}




void Simulation::loadInitialConditions()
{
	for (int i = 0; i < sbmlModel->getNumSpecies(); ++i)
	{
		Species* sbmlSpecies = sbmlModel->getSpecies(i);
		speciesData.insert(pair<string, int>(sbmlSpecies->getName(), i));
		speciesValues(i) = (ParticleType)sbmlSpecies->getInitialAmount();
		proposedSpeciesValues(i) = (ParticleType)sbmlSpecies->getInitialAmount();
	}
}



void Simulation::loadInitialConditions(float noiseLevel)
{
	for (int i = 0; i < sbmlModel->getNumSpecies(); ++i)
	{
		Species* sbmlSpecies = sbmlModel->getSpecies(i);
		speciesData.insert(pair<string, int>(sbmlSpecies->getName(), i));
		double perturbedVal = (1.0+noiseLevel*2*((double)ranf()-0.5))*sbmlSpecies->getInitialAmount();
		speciesValues(i) = (ParticleType)perturbedVal;
		proposedSpeciesValues(i) = (ParticleType)perturbedVal;
	}
}






double Simulation::dFindSubstringAnnotationL1(string annotation, string s [2])
{
	int f1, f2;
	f1 = annotation.find(s[0]) + s[0].length();
	f2 = annotation.find(s[1]);
	string ret = annotation.substr(f1, f2-f1);
	return std::atof(ret.c_str());
}

string Simulation::sFindSubstringAnnotationL1(string annotation, string s [2])
{
	int f1, f2;
	f1 = annotation.find(s[0]) + s[0].length();
	f2 = annotation.find(s[1]);
	string ret = annotation.substr(f1, f2-f1);
	return ret;
}





double Simulation::dFindSubstringAnnotationL2(string annotation, string s)
{
	int f1;
	f1 = annotation.find(s) + s.length();
	string t (annotation,f1,annotation.length()-f1);
	f1 = t.find("\"");
	string ret = t.substr(0, f1);
	return std::atof(ret.c_str());
}

string Simulation::sFindSubstringAnnotationL2(string annotation, string s)
{
	int f1;
	f1 = annotation.find(s) + s.length();
	string t (annotation,f1,annotation.length()-f1);
	f1 = t.find("\"");
	string ret = t.substr(0, f1);
	return ret;
}
