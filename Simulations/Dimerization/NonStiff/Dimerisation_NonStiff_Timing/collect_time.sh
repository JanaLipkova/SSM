#!/bin/sh -f

inputFile=$1

for f in `find ${inputFile}/ -name output`; do tail -n1 $f | awk -F ' ' '{print $3}'; done 


