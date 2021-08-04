#!/bin/bash

numProc=2
addOpt=
if [ $# -gt 0 ]
then
    numProc=$1
fi
if [ $# -gt 1 ]
then
    shift
fi
#echo npis${numProc}
#echo optis$@
mpirun -np ${numProc} --oversubscribe ./main.exe $@ -ksp_monitor
