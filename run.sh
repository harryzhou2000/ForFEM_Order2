#!/bin/bash

numProc=2
additional=
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

#additional="-start_in_debugger ${additional}"
#additional="-ksp_monitor ${additional}"

mpirun -np ${numProc} --oversubscribe ./main.exe $@  ${additional}
