#! /bin/bash -f


# modules
source /etc/profile.d/modules.sh
module load boost java/1.6
module unload gcc
module load python
module load gsl gcc blast

#libraries
export LD_LIBRARY_PATH=$HOME/ssm-libs/libSBML-5.8.0-Source/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/ssm-libs/blitz-0.10/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$BOOST_LIBDIR:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GSL_LIBDIR:$LD_LIBRARY_PATH

./ssm lacz_lacy_TAU_0.01.xml $1 $2
rm ssm
