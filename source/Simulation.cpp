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
	
	string timeStart			[2];
	string timeEnd				[2];
	string storeInterval		[2];
	string epsilon				[2];
	string numberOfSamples		[2];
	string method				[2];
	string theta				[2];
	string sortInterval			[2];
	string initialNoise			[2];
	string noiseIncrement		[2];
	string numberOfNoiseLevels	[2];
	string tau					[2];
	string delta				[2];

	
	timeStart[0]			= "<stochSim:TimeStart>";
	timeStart[1]			= "</stochSim:TimeStart>";
	
	timeEnd  [0]			= "<stochSim:TimeEnd>";
	timeEnd  [1]			= "</stochSim:TimeEnd>";
	
    storeInterval[0]		= "<stochSim:StoreInterval>";
	storeInterval[1]		= "</stochSim:StoreInterval>";
	
    epsilon[0]				= "<stochSim:Epsilon>";
	epsilon[1]				= "</stochSim:Epsilon>";
	
    numberOfSamples[0]		= "<stochSim:NumberOfSamples>";
	numberOfSamples[1]		= "</stochSim:NumberOfSamples>";	
	
	method[0]				= "<stochSim:Method>";
	method[1]				= "</stochSim:Method>";
	
	theta[0]				= "<stochSim:Theta>";
	theta[1]				= "</stochSim:Theta>";
	
	sortInterval[0]			= "<stochSim:SortInterval>";
	sortInterval[1]			= "</stochSim:SortInterval>";

	initialNoise[0]			= "<stochSim:InitialNoise>";
	initialNoise[1]			= "</stochSim:InitialNoise>";

	noiseIncrement[0]		= "<stochSim:NoiseIncrement>";
	noiseIncrement[1]		= "</stochSim:NoiseIncrement>";
	
	numberOfNoiseLevels[0]	= "<stochSim:NumberOfNoiseLevels>";
	numberOfNoiseLevels[1]	= "</stochSim:NumberOfNoiseLevels>";

	tau[0]					= "<stochSim:Tau>";
	tau[1]					= "</stochSim:Tau>";
	
	delta[0]				= "<stochSim:Delta>";
	delta[1]				= "</stochSim:Delta>";
	
	//string annotation = sbmlModel->getAnnotationString(); //libSBML-4.0.1	StartTime		= dFindSubstringAnnotation(annotation,		timeStart);
	string annotation = sbmlModel->getAnnotationString();
	
	StartTime = dFindSubstringAnnotation(annotation,		timeStart);

	EndTime			= dFindSubstringAnnotation(annotation,		timeEnd);
	StoreInterval   = (int)dFindSubstringAnnotation(annotation,	storeInterval);
	Epsilon			= dFindSubstringAnnotation(annotation,		epsilon);
	NumberOfSamples = (int)dFindSubstringAnnotation(annotation,	numberOfSamples);
	
	StochasticSimulationMethod = sFindSubstringAnnotation(annotation, method); 
	
	SortInterval	= (int)dFindSubstringAnnotation(annotation,	sortInterval);
	Theta			= dFindSubstringAnnotation(annotation,	theta);
	
	InitialNoise= dFindSubstringAnnotation(annotation,	initialNoise);
	NoiseIncrement = dFindSubstringAnnotation(annotation,	noiseIncrement);
	NumberOfNoiseLevels= (int)dFindSubstringAnnotation(annotation,	numberOfNoiseLevels);

	Tau			= dFindSubstringAnnotation(annotation,		tau);
	Delta			= dFindSubstringAnnotation(annotation,		delta);

	
	
	
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

double Simulation::dFindSubstringAnnotation(string annotation, string s [2])
{
	int f1, f2;
	f1 = annotation.find(s[0]) + s[0].length();
	f2 = annotation.find(s[1]);
	string ret = annotation.substr(f1, f2-f1);
	return std::atof(ret.c_str());
}

string Simulation::sFindSubstringAnnotation(string annotation, string s [2])
{
	int f1, f2;
	f1 = annotation.find(s[0]) + s[0].length();
	f2 = annotation.find(s[1]);
	string ret = annotation.substr(f1, f2-f1);
	return ret;
}

//void Simulation::parseAnnotation()
//{
//
//}



























