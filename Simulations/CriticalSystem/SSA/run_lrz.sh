#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J SSA
#SBATCH --clusters=mpp2
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=28
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=05:00:00

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

#Threads 
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
program=ssm
system=*.xml

echo "In the directory: $PWD"
echo "Running program on $SLURM_NODES nodes, with $SLURM_CPUS_ON_NODE cores on node, each with $SLURM_CPUS_PER_TASK cores."

./$program  $system
