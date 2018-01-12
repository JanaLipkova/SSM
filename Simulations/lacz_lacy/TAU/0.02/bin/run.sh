#! /bin/bash -f

module load python
./ssm lacz_lacy_TAU_0.02.xml $1 $2
rm ssm
