#! /bin/bash -f

module load python
./ssm lacz_lacy_SSA_0.02.xml $1 $2
rm ssm
