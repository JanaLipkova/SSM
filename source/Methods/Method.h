/*
 *  Method.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */


#pragma once
#include "../HeaderFiles.h"
#include "../Simulation.h"
#include "../SSMReaction.h"

class Method
	{

	public:
		Method(Simulation * simulation)
		{
			this->simulation = simulation;
			this->sbmlModel  = simulation->getSBMLModel();
			propensitiesVector.resize(sbmlModel->getNumReactions());
			propensitiesVector = 0.0;

			t					= simulation->StartTime;
			t_old				= simulation->StartTime;
			tEnd				= simulation->EndTime;
			numberOfFrames		= simulation->StoreInterval;
			tDiff				= (tEnd - simulation->StartTime) / ((double)(numberOfFrames));
			whenToSave			= simulation->StartTime;
			numberOfSamples		= simulation->NumberOfSamples;
			outputFileName		= (simulation->ModelName) + "_Output.txt";
		}

		virtual ~Method(){};

		virtual void   solve() = 0;

	protected:


		double t;
		double t_old;
		double dt;
		double tEnd;

		int numberOfIterations; // # times the core of the SSA is executed
		int numberOfSamples; // # times the complete simulation is run

		int timePoint;		// measurement timepoint for statistical analysis of the ensemble of simulations
		double whenToSave;  // measurement timepoint for statistical analysis of the ensemble of simulation
		int numberOfFrames; // # measurements taken during a single simulation
		double tDiff;       // time difference between measurements

		string outputFileName;

		Simulation * simulation;
		Model      * sbmlModel;

		Array< double , 1 > propensitiesVector;

		// Delay variables
		Array<    int      , 1 > delayedReactionsIndices;

		// the queue
		vector<double> delayedReactionsTimePoints ;
		vector<int   > delayedReactionsTimeIndices;

		// additional part to the queue for the leaping scheme
		vector<int   > delayedReactionsTimeLeapLength;

		// auxiliary stream
		ofstream auxiliaryStream; // used for histograms

#ifdef DEBUG_PRINT
		Array<double, 1> tempArray;
		ofstream myfile;
#endif





#pragma mark - Standard SSA Methods -

		virtual void computeCummulativeSum(Array< double , 1 > & cummulativeSum, int startingIndex)
		{
			cummulativeSum(0) = propensitiesVector(0);

			int start;
			if ( startingIndex <= 0)
			{
				start  = 1;
			}
			else
			{
				start = startingIndex;
			}

			for (int j = start; j < propensitiesVector.extent(firstDim); ++j)
			{
				cummulativeSum(j) = cummulativeSum(j-1) + propensitiesVector(j);
			}
		}

		// calculation of propensities using reaction-pointers
		virtual void computePropensitiesQuickly( Array< double , 1 > & propensitiesVector, int reactionFired, int & minReactionIndex )
		{
			vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;

			// add the reactants and products of the previously fired reaction to a set
			set < int >					reactantsAndProducts;
			SSMReaction * reaction		= ssmReactionList[reactionFired];
			vector <int>  reactants		= reaction->getReactants();
			for (unsigned int r = 0; r < reactants.size(); ++r)
			{
				reactantsAndProducts.insert( reactants[r] );
			}

			vector <int>  products		= reaction->getProducts();
			for (unsigned int p = 0; p < products.size(); ++p)
			{
				reactantsAndProducts.insert( products[p] );
			}

			// tag the reactions that need to be updated and save the minimum value
			minReactionIndex = ssmReactionList.size() - 1;

			set < int >					reactionsToUpdate;
			set < int >::const_iterator position;
			for (position = reactantsAndProducts.begin(); position != reactantsAndProducts.end(); ++position)
			{
				int value = *position;
				vector <int> & temp = simulation->backPointers[value];
				for (unsigned int j = 0; j < temp.size(); ++j)
				{
					int reactionIndex = temp[j];
					reactionsToUpdate.insert( reactionIndex );

					if (reactionIndex < minReactionIndex)
					{
						minReactionIndex = reactionIndex;
					}
				}
			}


			// compute the propensities for a subset of the reactions
			int nu;
			ParticleType x;
			ParticleType num, denom;
			for (position = reactionsToUpdate.begin(); position != reactionsToUpdate.end(); ++position)
			{
				int ir = *position;

				SSMReaction* reaction		= ssmReactionList[ir];
				vector <int>  reactants		= reaction->getReactants();
				vector <int>  nu_reactants	= reaction->getNuReactants();

				reaction->setPropensity(reaction->getRate());

				for (int s = 0; s < reactants.size(); ++s)
				{
					nu		= nu_reactants[s];
					x		= simulation->speciesValues( reactants[s] );
					num		= x;
					denom	= nu;
					while ((--nu)>0)
					{
						denom	*= nu;
						num		*= (x - nu);
					}
					reaction->setPropensity( reaction->getPropensity()*((double)num/(double)denom) );
				}
				propensitiesVector(ir) = reaction->getPropensity();
			}
		}


		// standard calculation of propensities
		virtual void computePropensities(Array< double , 1 > & propensitiesVector, int startingValue, int  endingValue)
		{
			int nu;
			ParticleType x;
			ParticleType num, denom;

			vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;

			int start = startingValue;
			if (startingValue < 0)
			{
				start = 0;
			}

			if (endingValue >= ssmReactionList.size())
			{
				endingValue = ssmReactionList.size();
			}

			for (int ir = start; ir < endingValue; ++ir)
			{

				SSMReaction* reaction		= ssmReactionList[ir];
				vector <int>  reactants		= reaction->getReactants();
				vector <int>  nu_reactants	= reaction->getNuReactants();

				reaction->setPropensity(reaction->getRate());

				for (int s = 0; s < reactants.size(); ++s)
				{
					nu		= nu_reactants[s];
					x		= simulation->speciesValues( reactants[s] );
					num		= x;
					denom	= nu;
					while ((--nu)>0)
					{
						denom	*= nu;
						num		*= (x - nu);
					}
					reaction->setPropensity( reaction->getPropensity()*((double)num/(double)denom) );
				}
				propensitiesVector(ir) = reaction->getPropensity();
			}
		}


		virtual void computePropensities(Array< double , 1 > & propensitiesVector, int startingValue)
		{
			int nu;
			ParticleType x;
			ParticleType num, denom;

			vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;

			int start = startingValue;
			if (startingValue < 0)
			{
				start = 0;
			}

			Reaction * reaction;
			for (int ir = 0; ir < sbmlModel->getNumReactions(); ++ir)
			{

				reaction = sbmlModel->getReaction(ir);
				KineticLaw * kineticLaw = reaction->getKineticLaw();
				Parameter  * parameter   = kineticLaw->getParameter(0);
				double rate = parameter->getValue();

				SSMReaction* reaction		= ssmReactionList[ir];
				vector <int>  reactants		= reaction->getReactants();
				vector <int>  nu_reactants	= reaction->getNuReactants();



				reaction->setPropensity(reaction->getRate());


				for (int s = 0; s < reactants.size(); ++s)
				{
					nu		= nu_reactants[s];
					x		= simulation->speciesValues( reactants[s] );
					num		= x;
					denom	= nu;

					while ((--nu)>0)
					{
						denom	*= nu;
						num		*= (x - nu);
					}

					reaction->setPropensity( reaction->getPropensity()*((double)num/(double)denom) );
				}


				propensitiesVector(ir) = reaction->getPropensity();



			}

		}

		// for systems with linearly growing volume, i.e.
		//		V = 1 + time/genTime  // current time, generation time
		// propensities of hiher order reactions must be rescaled by current volume
		// for more info see: "A.M. Kriezek: Stocks: STOCHastic Kinetics Simulations ... , (2001)"
		virtual void computePropensitiesGrowingVolume(Array< double , 1 > & propensitiesVector, double time, double genTime)
		{
			double volume	= 1. + time/genTime;
			double ivolume	= 1./volume;

			int nu;
			ParticleType x;
			ParticleType num, denom;

			vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;

			for (int ir = 0; ir < ssmReactionList.size(); ++ir)
			{
				SSMReaction* reaction		= ssmReactionList[ir];
				vector <int>  reactants		= reaction->getReactants();
				vector <int>  nu_reactants	= reaction->getNuReactants();
				int order					= reaction->getOrder();

				reaction->setPropensity(reaction->getRate());

				for (int s = 0; s < reactants.size(); ++s)
				{
					nu		= nu_reactants[s];
					x		= simulation->speciesValues( reactants[s] );
					num		= x;
					denom	= nu;
					while ((--nu)>0)
					{
						denom	*= nu;
						num		*= (x - nu);
					}
					reaction->setPropensity( reaction->getPropensity()*((double)num/(double)denom) );
				}

				if (order == 2)
					reaction->setPropensity( reaction->getPropensity() * ivolume );

				if (order == 3)
					reaction->setPropensity( reaction->getPropensity() * ivolume *ivolume);

				if (order > 3)
				{
					std::cout<<"Aborting: Growing volume of reaction enviroment do not support reaction of order higher than 3, if you want it implement it"<<std::endl;
					std::abort();
				}

				propensitiesVector(ir) = reaction->getPropensity();
			}
		}

		// return false if there is a negative species
		virtual void fireReaction(int reactionIndex, int numberOfTimes)
		{

			SSMReaction * r = simulation->ssmReactionList[reactionIndex];
			vector <int> changes = r->getChanges();
			vector <int> nuChanges = r->getNuChanges();

			// XXX
			simulation->old_speciesValues = simulation->speciesValues;



			for (int i = 0; i < changes.size(); ++i)
			{
				simulation->speciesValues(changes[i]) += (nuChanges[i]*((ParticleType)(numberOfTimes)));
			}

			/*
			 Reaction * reaction = sbmlModel->getReaction(reactionIndex);
			 SpeciesReference * speciesReference;
			 string speciesName;
			 int speciesIndex;

			 for (int j = 0; j < reaction->getNumReactants(); ++j)
			 {
			 speciesReference = reaction->getReactant(j);
			 speciesName = speciesReference->getSpecies();
			 speciesIndex = simulation->getSpeciesIndex(speciesName);
			 simulation->speciesValues(speciesIndex) -= (ParticleType)(numberOfTimes);
			 }
			 for (int j = 0; j < reaction->getNumProducts(); ++j)
			 {
			 speciesReference = reaction->getProduct(j);
			 speciesName = speciesReference->getSpecies();
			 speciesIndex = simulation->getSpeciesIndex(speciesName);
			 simulation->speciesValues(speciesIndex) += (ParticleType)(numberOfTimes);
			 }
			 */
		}

#pragma mark - Negative Species Methods -
		virtual void fireReactionProposed(int reactionIndex, int numberOfTimes)
		{
			SSMReaction * r = simulation->ssmReactionList[reactionIndex];
			vector <int> changes = r->getChanges();
			vector <int> nuChanges = r->getNuChanges();


			for (int i = 0; i < changes.size(); ++i)
			{
				simulation->proposedSpeciesValues(changes[i]) += (nuChanges[i]*((ParticleType)(numberOfTimes)));
			}

			/*
			 Reaction * reaction = sbmlModel->getReaction(reactionIndex);
			 SpeciesReference * speciesReference;
			 string speciesName;
			 int speciesIndex;

			 for (int j = 0; j < reaction->getNumReactants(); ++j)
			 {
			 speciesReference = reaction->getReactant(j);
			 speciesName = speciesReference->getSpecies();
			 speciesIndex = simulation->getSpeciesIndex(speciesName);
			 simulation->proposedSpeciesValues(speciesIndex) -= (ParticleType)(numberOfTimes);
			 }
			 for (int j = 0; j < reaction->getNumProducts(); ++j)
			 {
			 speciesReference = reaction->getProduct(j);
			 speciesName = speciesReference->getSpecies();
			 speciesIndex = simulation->getSpeciesIndex(speciesName);
			 simulation->proposedSpeciesValues(speciesIndex) += (ParticleType)(numberOfTimes);
			 }
			 */
		}

		virtual bool isProposedNegative()
		{
			ParticleType minValue = blitz::min(simulation->proposedSpeciesValues);
			if (minValue < (ParticleType)0)
				return true;
			else
				return false;
		}

		virtual void acceptNewSpeciesValues()
		{
			// XXX
			simulation->old_speciesValues = simulation->speciesValues;

			simulation->speciesValues = simulation->proposedSpeciesValues;
			simulation->proposedSpeciesValues = simulation->speciesValues;
		}

		virtual void reloadProposedSpeciesValues()
		{
			simulation->proposedSpeciesValues = simulation->speciesValues;
		}

#pragma mark - Processing Methods -

		virtual void openAuxiliaryStream(string filename)
		{
			auxiliaryStream.open (filename.c_str(), ios::out);
		}

		virtual void writeToAuxiliaryStream( Array< ParticleType, 1 > & speciesValues )
		{
			for (int j = 0; j < speciesValues.extent(firstDim); ++j)
			{
				auxiliaryStream << speciesValues(j) << "\t";
			}
			auxiliaryStream << endl;
		}

		virtual void closeAuxiliaryStream()
		{
			auxiliaryStream.close();
		}





		void saveData()
		{
			// cout << "---------->" << t << "  --   "<< whenToSave << "  --   "<< t_old << endl ;
			if( t_old <= whenToSave && t > whenToSave){
				while( timePoint<=numberOfFrames && whenToSave <= t ){

					cout << "	 Saving data at time t = " << whenToSave << endl;
					simulation->speciesEnsemble(Range::all(), timePoint) += simulation->old_speciesValues;
					++timePoint;
					whenToSave += tDiff;

				}
			}

			// if (t >= whenToSave)
			// {
            //
			// 	simulation->speciesEnsemble(Range::all(), timePoint) += simulation->speciesValues;
			// 	++timePoint;
			// 	whenToSave += tDiff;
            //
			// 	cout << "	 Saving data at time t = " << t << endl;
			// }

		}



		virtual void writeData(string filename)
		{
			ofstream myfile;
			myfile.open (filename.c_str(), ios::out);
			myfile << scientific;
			myfile << setprecision (9);
			Array<double, 1> tempArray(sbmlModel->getNumSpecies());
			for (int j = 0; j < simulation->speciesEnsemble.extent(secondDim); ++j)
			{
				// write the time
				myfile << ((double)j)*tDiff << "\t";

				// write the data
				tempArray = simulation->speciesEnsemble(Range::all(), j);
				for (int i = 0; i < tempArray.extent(firstDim); ++i)
                {
					myfile << (tempArray(i) / (double)numberOfSamples) << "\t";
				}
				myfile << endl;
			}
			myfile.close();
		}

		virtual void writeData(string filename, int filenumber)
		{
			string fileNumber = boost::lexical_cast<std::string>(filenumber);

			ofstream myfile;
			int stringLength = filename.length();
			string fileToOpen = filename.substr(0, stringLength-4);
			fileToOpen = fileToOpen + fileNumber + ".txt";
			myfile.open (fileToOpen.c_str(), ios::out);
			Array<double, 1> tempArray(sbmlModel->getNumSpecies());
			for (int j = 0; j < simulation->speciesEnsemble.extent(secondDim); ++j)
			{
				// write the time
				myfile << ((double)j)*tDiff << "\t";

				// write the data
				for (int i = 0; i < tempArray.extent(firstDim); ++i)
				{
					myfile << (tempArray(i) ) << "\t";
				}
				myfile << endl;
			}
			myfile.close();
		}


		virtual void zeroData() // careful...
		{
			simulation->speciesEnsemble = (ParticleType)0.0;
		}

		virtual void writeSTLVector( vector <double> & v, string fileName)
		{
			FILE* myfile;
			cout << "---------------writing file: " << fileName << " ---------------" << endl;

			myfile = fopen(fileName.c_str(), "w");

			// write the data
			for(int i = 0; i < v.size(); ++i)
			{
				fprintf(myfile, "%d\t%e\n",  i, v[i] );
			}

			fclose(myfile);
		}


#pragma mark - Delay SSA Methods -

		// TODO
		// correct getFormula().substr() to use Species Names instead of Numbers...
		virtual int getDependentSpecies(int reactionIndex)
		{
			Reaction * reaction = sbmlModel->getReaction(reactionIndex);
			KineticLaw * kineticLaw = reaction->getKineticLaw();
			string speciesS = kineticLaw->getFormula().substr(8);
			double speciesD = std::atof(speciesS.c_str());
			return ((int)speciesD);
		}

		virtual int getDependentSpecies(string f)
		{
			string speciesS = f.substr(8);
			double speciesD = std::atof(speciesS.c_str());
			return ((int)speciesD);
		}

		virtual void computePropensities(Array< double , 1 > & propensitiesVector)
		{
			propensitiesVector = 0.0;
			//int lastIndex = -1;
			Reaction * reaction;
			for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
			{
				reaction = sbmlModel->getReaction(i);
				KineticLaw * kineticLaw = reaction->getKineticLaw();
				Parameter * parameter = kineticLaw->getParameter(0);
				double rate = parameter->getValue();

				if ( isReactionDelayed(i) == true && ( kineticLaw->getNumParameters() == 5) )
				{
					int dependentSpecies = getDependentSpecies(i);
					double dependentValue = (double)simulation->speciesValues(dependentSpecies);
					double h =					kineticLaw->getParameter(2)->getValue();
					double defaultProduction =	kineticLaw->getParameter(3)->getValue();
					double cHill =				kineticLaw->getParameter(4)->getValue();
					propensitiesVector(i) = defaultProduction + rate*hillFunction(cHill, dependentValue, h);

				}
				else if ( isReactionDelayed(i) == true && ( kineticLaw->getNumParameters() > 4) )
				{
					//cout << "cooking"<< endl;

					int dependentSpecies[3] = {0, 0, 0};
					double criticalValues[3] = {0.0, 0.0, 0.0};


					dependentSpecies[0] = getDependentSpecies( kineticLaw->getParameter(4)->getName() );
					criticalValues[0] = kineticLaw->getParameter(5)->getValue();

					dependentSpecies[1] = getDependentSpecies( kineticLaw->getParameter(6)->getName() );
					criticalValues[1] = kineticLaw->getParameter(7)->getValue();

					dependentSpecies[2] = getDependentSpecies( kineticLaw->getParameter(8)->getName() );
					criticalValues[2] = kineticLaw->getParameter(9)->getValue();


					double gamma[4] = {0.0, 0.0, 0.0, 0.0};
					for (int g = 0; g < 4; ++g)
					{
						gamma[g] = kineticLaw->getParameter(10+g)->getValue();;
					}



					propensitiesVector(i) = rate*hillFunctionCooking(
																	 ((double)simulation->speciesValues(dependentSpecies[0])), criticalValues[0],
																	 ((double)simulation->speciesValues(dependentSpecies[1])), criticalValues[1],
																	 ((double)simulation->speciesValues(dependentSpecies[2])), criticalValues[2],
																	 gamma );

					//propensitiesVector(i) = rate*hillFunction(((double)simulation->speciesValues(dependentSpecies)), P0, h);
				}
				else if ( isReactionDelayed(i) == true && ( kineticLaw->getNumParameters() <= 4) )
				{
					//cout << "standard delays" << endl;
					int dependentSpecies = getDependentSpecies(i);
					double h =   kineticLaw->getParameter(2)->getValue();
					double P0 =  kineticLaw->getParameter(3)->getValue();
					propensitiesVector(i) = rate*hillFunction((double)simulation->speciesValues(dependentSpecies), P0, h);
				}
				else
				{
					propensitiesVector(i) = rate;
					if (kineticLaw->getFormula().substr(0, 8) == "Species:")
					{
						int dependentSpecies = getDependentSpecies(i);
						propensitiesVector(i) *= (double)simulation->speciesValues(dependentSpecies);
					}

					if ( reaction->getNumReactants() == 0 )
					{
					}
					else if ( reaction->getNumReactants() == 1 )
					{
						SpeciesReference * speciesReference = reaction->getReactant(0);
						string speciesName = speciesReference->getSpecies();
						int speciesIndex = simulation->getSpeciesIndex(speciesName);
						double value = (double)simulation->speciesValues(speciesIndex);

						propensitiesVector(i) *= value;
					}
					else
					{
						int nu;
						ParticleType x;
						ParticleType num, denom;

						int ir = i;
						vector<SSMReaction* > ssmReactionList = simulation->ssmReactionList;
						SSMReaction* reaction		= ssmReactionList[ir];
						vector <int>  reactants		= reaction->getReactants();
						vector <int>  nu_reactants	= reaction->getNuReactants();

						reaction->setPropensity(reaction->getRate());

						for (int s = 0; s < reactants.size(); ++s)
						{
							nu		= nu_reactants[s];
							x		= simulation->speciesValues( reactants[s] );
							num		= x;
							denom	= nu;
							while ((--nu)>0)
							{
								denom	*= nu;
								num		*= (x - nu);
							}
							reaction->setPropensity( reaction->getPropensity()*((double)num/(double)denom) );
						}
						propensitiesVector(ir) = reaction->getPropensity();
					}
				}
			}
		}

		virtual double hillFunctionCooking(double Pher1, double Pher1Critical,
										   double Pher7, double Pher7Critical,
										   double NeighborDelta, double NeighborDeltaCritical,
										   double gamma[] )
		{
			double phidhat = NeighborDelta / NeighborDeltaCritical;
			double phiher1 = Pher1 / Pher1Critical;
			double phiher7 = Pher7 / Pher7Critical;

			double cooking = (gamma[0]) +
							 (gamma[1]*phidhat / (1.0 + phidhat)) +
							 (gamma[2] / (1.0 + phiher1*phiher7)) +
							 ( gamma[3]*(phidhat/(1.0 + phidhat))*(1.0 / (1.0 + phiher1*phiher7) ) ) ;

			return cooking;
		}

		/*
		 hillFunction(Ka, L, n)
		 Calculates a generic HillFunction
		                     L^n
		    h(Ka,L,n) = ––––––––––––
		                 Ka^n + L^n
		 */
		virtual double hillFunction(double Ka, double L, double n)
		{
//			if (h < 1.0e-6)
//			{
//				return Pt;
//			}
//			else
//				return ( 1.0 / ( 1.0 + pow( (Pt / P0) , h) )  );
			if (n == 1)
			{
				return L / (Ka + L);
			}
			else {
				double Ln = pow(L, n);
				return (Ln / (pow(Ka, n) + Ln));
			}
		}

		virtual int numberOfDelayedReactions()
		{
			int num = 0;
			Reaction * reaction;
			for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
			{
				reaction = sbmlModel->getReaction(i);
				KineticLaw * kineticLaw = reaction->getKineticLaw();
				if (kineticLaw->getNumParameters() > 1)
				{ ++num; }
			}
			return num;
		}

		virtual void setDelayedReactionsTime(int reactionIndex, double t, double dt, int k)
		{
			Reaction * reaction = sbmlModel->getReaction( reactionIndex );
			KineticLaw * kineticLaw = reaction->getKineticLaw();

			delayedReactionsTimePoints.push_back( t + dt + kineticLaw->getParameter(1)->getValue() );
			delayedReactionsTimeIndices.push_back( reactionIndex );
			delayedReactionsTimeLeapLength.push_back( k );
		}

		virtual void removeDelayedReactionsTime(int delayIndex)
		{
			delayedReactionsTimeIndices.erase( delayedReactionsTimeIndices.begin() + delayIndex, delayedReactionsTimeIndices.begin()+delayIndex+1);
			delayedReactionsTimePoints.erase( delayedReactionsTimePoints.begin() + delayIndex, delayedReactionsTimePoints.begin()+delayIndex+1);
			delayedReactionsTimeLeapLength.erase( delayedReactionsTimeLeapLength.begin() + delayIndex, delayedReactionsTimeLeapLength.begin()+delayIndex+1);
		}

		virtual double getDelayTime(int reactionIndex)
		{
			Reaction * reaction = sbmlModel->getReaction( reactionIndex );
			KineticLaw * kineticLaw = reaction->getKineticLaw();
			return (kineticLaw->getParameter(1)->getValue());
		}

		virtual void setDelayedReactionsIndices()
		{
			int num = 0;
			Reaction * reaction;
			for (int i = 0; i < sbmlModel->getNumReactions(); ++i)
			{
				reaction = sbmlModel->getReaction(i);
				KineticLaw * kineticLaw = reaction->getKineticLaw();
				if (kineticLaw->getNumParameters() > 1)
				{
					delayedReactionsIndices(num) = i;
					++num;
				}
			}
		}

		virtual bool isDelayedReactionScheduled(double t, double dt, int & delayIndex)
		{
			for (int i = 0; i < delayedReactionsTimePoints.size(); ++i)
			{
				if ( (delayedReactionsTimePoints[i] >= t) && (delayedReactionsTimePoints[i] < (t + dt) ))
				{
					delayIndex = i;
					return true;
				}
			}
			return false;
		}

		virtual bool isReactionDelayed(int reactionIndex)
		{
			for (int i = 0; i < delayedReactionsIndices.extent(firstDim); ++i)
			{
				if ( delayedReactionsIndices(i) == reactionIndex )
				{
					return true;
				}
			}
			return false;
		}


	};
