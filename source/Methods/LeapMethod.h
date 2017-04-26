/*
 *  RLeapMethod.h
 *  SSM
 *
 *  Created by Martin Maag on 15.06.2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once

#include "../HeaderFiles.h"
#include "Method.h"

class LeapMethod : public Method
{
public:
	LeapMethod(Simulation * simulation): Method(simulation)
	{}


		protected:
	
	

#pragma mark - Leaping Error Control - 	  
	virtual void computeHor(Array<int, 1> & hor, Array<int, 1> & nuHor)
	{
		for (int numbS = 0; numbS < sbmlModel->getNumSpecies(); ++numbS)
		{
			hor(numbS) = 0;
			nuHor(numbS) = 0;
		}
		
		int numberOfReactions = sbmlModel->getNumReactions();
		for (int ir = 0; ir < numberOfReactions; ++ir)
		{
			SSMReaction* ssmReaction = simulation->ssmReactionList[ir];
			
			const vector<int> & reactantsVector  = ssmReaction->getReactants();
			const vector<int> & nuReactantsVector = ssmReaction->getNuReactants();
			int order = 0;
			for (int is = 0; is < reactantsVector.size(); ++is)
			{
				order += nuReactantsVector[is];
			}
			for (int is = 0; is < reactantsVector.size(); ++is)
			{
				if (order > hor(reactantsVector[is]))
				{
					hor(reactantsVector[is]) = order;
					nuHor(reactantsVector[is]) = nuReactantsVector[is];
				}
			}
			
		}
	}
	
	virtual void computeMuHatSigmaHat2(Array<double, 1> & muHat, Array<double, 1> & sigmaHat2)
	{
		int is, ir, ns, indx, nr;
		double tmpfloat;
		nr = sbmlModel->getNumReactions();
		
		for (int numbS = 0; numbS < sbmlModel->getNumSpecies(); ++numbS)
		{
			muHat(numbS) = 0.0;
			sigmaHat2(numbS) = 0.0;
		}
		
		for (ir = 0; ir < nr; ++ir)
		{
			SSMReaction* ri = simulation->ssmReactionList[ir];
			double  riPropensity = propensitiesVector(ir);
			
			const vector<int> & changes = ri->getChanges();
			const vector<int> & nuChanges = ri->getNuChanges();
			
			ns = changes.size();
			for (is = 0; is < ns; is++ )
			{
				indx = changes[is];
				tmpfloat = nuChanges[is] * riPropensity;
				muHat(indx) += tmpfloat;
				sigmaHat2(indx) += nuChanges[is] * tmpfloat;
			}
		}
		
		//cout << "mmmu:"<<muHat << " sifma1 :"<<sigmaHat2 << endl;
	}
	/*{
	 int is, ir, ns, indx, nr;
	 double tmpfloat;
	 nr = sbmlModel->getNumReactions();
	 for (ir = 0; ir < nr; ++ir)
	 {
	 SSMReaction* ri = simulation->ssmReactionList[ir];
	 double  riPropensity = propensitiesVector(ir);
	 
	 const vector<int> & reactants = ri->getReactants();
	 const vector<int> & nuReactants = ri->getNuReactants();
	 
	 ns = reactants.size();
	 for (is = 0; is < ns; is++ )
	 {
	 indx = reactants[is];
	 tmpfloat = -nuReactants[is] * riPropensity;
	 muHat(indx) += tmpfloat;
	 sigmaHat2(indx) += -nuReactants[is] * tmpfloat;
	 }
	 }
	 }*/
	virtual double computeTimeStep()
	{
		double tauPrime;
		
		double epsilon	= simulation->Epsilon;
		
		int numberOfSpecies		= sbmlModel->getNumSpecies();
		Array<int, 1> hor			(numberOfSpecies);
		Array<int, 1> nuHor			(numberOfSpecies);
		Array<double, 1> muHat		(numberOfSpecies);
		Array<double, 1> sigmaHat2	(numberOfSpecies);
		Array<double, 1> varHat		(numberOfSpecies);
		hor = 0; nuHor = 0; muHat = 0.0; sigmaHat2 = 0.0;
		
		computeHor(hor, nuHor);
		computeMuHatSigmaHat2(muHat, sigmaHat2);
		
		double tau, taup,  epsi, epsixi, epsixisq;
		double xi;
		
		tau = HUGE_VAL;
		
		double a0 = (double)blitz::sum(propensitiesVector);
		
		//cout << "a0: " << a0 << endl;
		
		for (int is = 0; is < numberOfSpecies; is++)
		{
			varHat(is) = sigmaHat2(is) - (1.0/a0) * muHat(is) * muHat(is);
		}
		
		for (int is = 0; is < numberOfSpecies; ++is)
		{
			taup = (HUGE_VALF*0.5);
			xi = (double)simulation->speciesValues(is);
			switch (hor(is)) {
				case 0:
					break;
				case 1:
					epsi = epsilon;
					epsixi = epsi * xi;
					epsixi = max(epsixi,1.0);
					tau = min(tau,epsixi/fabsf(muHat(is)));
					epsixisq = epsixi*epsixi;
					tau = min(tau,epsixisq/varHat(is));
					
					//cout << "		species: " << is << "	tau: " << tau << endl << endl;
					
					
					break;
				case 2:
					if (nuHor(is) == 1)
						epsi = 0.5*epsilon;
					else
						epsi = epsilon*(xi-1.0)/(2.0*(xi-1.0)+1.0);
					epsixi = epsi * xi;
					epsixi = max(epsixi,1.0);
					tau = min(tau,epsixi/fabs(muHat(is)));
					epsixisq = epsixi*epsixi;
					tau = min(tau,epsixisq/varHat(is));
					break;
				case 3:
					if (nuHor(is)==1)
						epsi = 0.3333333333*epsilon;
					else if (nuHor(is) == 2)
						epsi = epsilon*(xi-1)/(3.0*(xi-1)+1.5);
					else
						epsi = epsilon*(xi-1)*(xi-2)/(3.0*(xi-1)*(xi-2)+(xi-2)+2.0*(xi-1));
					epsixi = epsi * xi;
					epsixi = max(epsixi,1.0);
					tau = min(tau,epsixi/fabsf(muHat(is)));
					epsixisq = epsixi*epsixi;
					tau = min(tau,epsixisq/varHat(is));
					break;
				default:
					break;
			}
		}
		tauPrime = tau;
		return tauPrime;
	}
	
	
};

