#! /bin/csh -f
#
#  
#  
#########################################################################################################

# for SBML
# Own libraries in ssm-libs
set SBML_INC_DIR=$HOME/ssm-libs/libSBML-5.8.0-Source/include
set SBML_LIB_DIR=$HOME/ssm-libs/libSBML-5.8.0-Source/lib

set BLITZ_INC_DIR=$HOME/ssm-libs/blitz-0.10
set BLITZ_LIB_DIR=$HOME/ssm-libs/blitz-0.10/lib

# LRZ libraries in /lrz/sys/libraries
set BOOST_INC_DIR=${BOOST_INCDIR}
set BOOST_LIB_DIR=${BOOST_LIBDIR}

set GSL_INC_DIR=${GSL_INCDIR}
set GSL_LIB_DIR=${GSL_LIBDIR}

#echo BOOST_INC=${BOOST_INC}
#echo GSL=${GSL_INCDIR}



set CPPFLAGS=($SBML_LIB_DIR/libsbml.a -L$BLITZ_LIB_DIR -lblitz -lexpat -lbz2 -lz -L$GSL_LIB_DIR -lgsl -lgslcblas -lm )

echo --------------------------------------------------------------------------------------------------------------
echo --------------------------------------LINKING--ERRORS-AND-WARNINGS--------------------------------------------
echo --------------------------------------------------------------------------------------------------------------

g++ *.o $CPPFLAGS -o StochasticSimulationMethods

#g++ *.o $EXTRA_LIB/libsbml.a -lblitz -lexpat -lbz2 -lz -lgsl -lgslcblas -lm -o StochasticSimulationMethods

#rm *.o

echo --------------------------------------------------------------------------------------------------------------
echo --------------------------------------------------------------------------------------------------------------

echo COMPILE_SCRIPT_FINISHED

echo --------------------------------------------------------------------------------------------------------------
echo --------------------------------------------------------------------------------------------------------------

