#!/bin/bash


./ssm lacy_lacz_SSA_0.05.xml
mv  lacy_lacz_Output.txt  out_ssa.txt


./ssm lacy_lacz_AdaptiveSLeaping_0.05.xml;
mv  lacy_lacz_Output.txt  out_sleap.txt


./ssm lacy_lacz_AdaptiveTau_0.05.xml;
mv  lacy_lacz_Output.txt  out_tau.txt


./ssm lacy_lacz_RLeapJana_0.05.xml;
mv  lacy_lacz_Output.txt  out_rleap.txt
