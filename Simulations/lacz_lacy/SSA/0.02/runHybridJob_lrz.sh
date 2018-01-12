#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J ex1
#SBATCH --clusters=mpp2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=00:02:00

source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default intel/16.0 mkl/11.3
module load gcc/4.9  mpi.ompi/1.10/gcc
module list

export LANG=C
export LC_ALL=C

export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH
echo "we use this mpicc:"
which mpicc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#srun --cpus-per-task=2 -n 2 ./ex1
#mpiexec --perhost=2 -n 2 ./ex1
#mpiexec -n 2 --perhost=2 ./ex1
mpirun -np 2 simulate_all 10
