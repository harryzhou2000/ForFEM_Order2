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

#config that works for cooler:
# additional="-ksp_monitor -pc_type asm -sub_pc_type icc\
#     -ksp_type fcg 	-ksp_gmres_restart 120 	-ksp_gmres_cgs_refinement_type\
#     -ksp_gmres_cgs_refinement_type refine_ifneeded \
#     -ksp_max_it 1000 -ksp_rtol 0.5e-9 -ksp_converged_reason\
#     ${additional}"

additional="-ksp_monitor -pc_type asm -sub_pc_type icc\
    -ksp_type fcg 	-ksp_gmres_restart 120 	-ksp_gmres_cgs_refinement_type\
    -ksp_gmres_cgs_refinement_type refine_ifneeded \
    -ksp_max_it 1000 -ksp_rtol 0.5e-9 -ksp_converged_reason\
    ${additional}"

mpirun -np ${numProc}  ./main_cooler.exe $@  ${additional}
