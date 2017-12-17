/*
 *  Simulation.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once
#include "HeaderFiles.h"
#include "SSMReaction.h"


class Simulation
{
public:

	// Constructors
	Simulation ();
	Simulation (SBMLDocument* sbmlDocument);
	~Simulation();
	
	// Public Members
	Array< ParticleType, 1 > speciesValues;
	Array< ParticleType, 1 > proposedSpeciesValues;
	Array<   double    , 2 > speciesEnsemble;
	
	vector<SSMReaction* >	ssmReactionList;
	vector < vector <int> > backPointers;  // matrix, (species) x (reactions in which the speices is involved)
	vector<int> Reversible;   // indecies of pairs of reversible reactions, i and i+1 is seen as reversible pair
	
	
	double StartTime, EndTime, Epsilon, Theta, InitialNoise, NoiseIncrement, Tau, Delta; // Simulation Parameters
	int    NumberOfSamples,NumberOfNoiseLevels, StoreInterval, SortInterval;
	
	string StochasticSimulationMethod;
	string ModelName;
	
	
	// Public Functions
	Model*	getSBMLModel();
	int		getSpeciesIndex(string speciesName);
	
	void	loadInitialConditions();
	void	loadInitialConditions(float noiseLevel);
	
	
	
private:
	
	map< string, int >				speciesData;
	map< string, int >::iterator	iter;
	
	SBMLDocument	*	sbmlDocument;
	Model			*	sbmlModel;
//	Reaction		*	sbmlReaction;
//	Species			*	sbmlSpecies;
	
	//
//	void		parseAnnotation();
	double		dFindSubstringAnnotation(string annotation, string s [2]);
	string		sFindSubstringAnnotation(string annotation, string s [2]);


	
};
