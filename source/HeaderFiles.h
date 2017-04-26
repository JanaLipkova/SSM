/*
 *  HeaderFiles.h
 *  StochasticSimulationMethods
 *
 *  Created by Basil Bayati on 5/5/08.
 *  Copyright 2008 Basil Bayati. All rights reserved.
 *
 */

// C/C++ Headers
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <cstring>
#include <string>



// Blitz++ Arrays
#include <blitz/array.h>
#include <random/gamma.h>

// STL Library
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <string>

// Random Number Library
#include "RNGLib/ranlib.h"

// Systems Biology Markup Language
//#include "sbml/SBMLTypes.h"
#include <sbml/SBMLTypes.h>

// Boost
#include <boost/algorithm/string/trim.hpp>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>

// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;
using namespace blitz;
//using namespace boost;
//LIBSBML_CPP_NAMESPACE_USE //libSBML-4.0.1

#ifndef PARTICLEDEFINED
	#define ParticleType long int
#endif

