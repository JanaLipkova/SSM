#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J Bsub
#SBATCH --clusters=mpp2
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=00:25:00

# Load torc codules
source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default intel/16.0 mkl/11.3
module load python
module load boost java/1.8 gsl blast gcc
#module load boost java/1.8
#module load python
#module load gsl blast
#module load gcc/4.9  mpi.ompi/1.10/gcc

module list

# LRZ enviroment set up
export LANG=C
export LC_ALL=C

#Setup paths
export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH
echo "we use this mpicc:"
which mpicc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "In the directory: $PWD"
echo "Running program with $SLURM_NTASKS MPI tasks, each with $SLURM_CPUS_PER_TASK threads."

#srun --cpus-per-task=2 -n 2 ./ex1
#mpiexec --perhost=2 -n 2 ./ex1
#mpiexec -n 2 --perhost=2 ./ex1

# number of ssm samples, ppn = mpi processes per host, max 28 on mpp2
SSMsamples=10024
ppn=28
mpirun -np $SLURM_NTASKS -env TORC_WORKERS 1 ./simulate_all $SSMsamples
