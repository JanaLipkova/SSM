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

	string findSubstringAnnotationL1(string annotation, string s [2]);
	string findSubstringAnnotationL2(string annotation, string s);

	string filename;
	SBMLDocument * document;

	string timeStartL1			[2];
	string timeEndL1			[2];
	string storeIntervalL1		[2];
	string epsilonL1			[2];
	string deltaL1				[2];
	string numberOfSamplesL1	[2];
	string methodL1				[2];
	string thetaL1				[2];
	string sortIntervalL1		[2];
	string initialNoiseL1[2];
	string noiseIncrementL1[2];
	string numberOfNoiseLevelsL1[2];


	string timeStartL2;
	string timeEndL2;
	string storeIntervalL2;
	string epsilonL2;
	string deltaL2;
	string numberOfSamplesL2;
	string methodL2;
	string thetaL2;
	string sortIntervalL2;
	string initialNoiseL2;
	string noiseIncrementL2;
	string numberOfNoiseLevelsL2;




};
