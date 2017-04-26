/*
 *  SBMLReaderAndParser.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/2/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

#pragma once
#include "HeaderFiles.h"

class SBMLReaderAndParser
{
public:

	SBMLReaderAndParser(string filename);
	~SBMLReaderAndParser();

	int readAndParse();
	SBMLDocument * getSBMLDocument();
private:

	string findSubstringAnnotation(string annotation, string s [2]);

	string filename;
	SBMLDocument * document;
	
	string timeStart		[2];
	string timeEnd			[2];
	string storeInterval	[2];
	string epsilon			[2];
	string delta			[2];
	string numberOfSamples	[2];
	string method			[2];
	string theta			[2];
	string sortInterval		[2];
	
	
	string initialNoise[2];
	
	string noiseIncrement[2];
	
	string numberOfNoiseLevels[2];
};
