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

additional="-ksp_monitor -ksp_type gmres -pc_type asm -sub_ksp_type preonly -sub_pc_type cholesky\
    -sub_pc_factor_levels 0 -sub_pc_factor_mat_solver_type superlu_dist\
    -ksp_gmres_restart 120 -ksp_gmres_cgs_refinement_type refine_ifneeded \
    -ksp_max_it 10000 -ksp_rtol 0.5e-5 -ksp_converged_reason\
    ${additional}"

# additional="-ksp_monitor -ksp_type gmres -pc_type asm -sub_ksp_type preonly -sub_pc_type icc\
#     -sub_pc_factor_levels 0 -sub_pc_factor_mat_solver_type petsc\
#     -ksp_gmres_restart 120 -ksp_gmres_cgs_refinement_type refine_ifneeded \
#     -ksp_max_it 10000 -ksp_rtol 0.5e-5 -ksp_converged_reason\
#     ${additional}"

    # -sub_ksp_max_it 100 -sub_ksp_rtol 0.5e-4 -sub_ksp_converged_reason\
additional="-eps_monitor -eps_converged_reason -eps_nev 5 -eps_ncv 20\
    -eps_type primme -st_pc_type cholesky -st_sub_pc_type cholesky \
    -eps_tol 1e-5 -eps_max_it 20000 \
    ${additional}"

# rqcg lobpcg primme krylovschur arpack arnoldi
# -st_factor_mat_solver_type superlu_dist\

#-eps_type gd  -st_type precond  \

#-eps_type primme -st_pc_type cholesky

# additional="-ksp_monitor -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist\
#     -ksp_gmres_restart 120 	-ksp_gmres_cgs_refinement_type\
#     -ksp_gmres_cgs_refinement_type refine_ifneeded \
#     -ksp_max_it 10000 -ksp_rtol 0.5e-9 -ksp_converged_reason\
#     ${additional}"

# mpirun -np ${numProc}  ./main_grid.exe $@  ${additional}
# mpirun -np ${numProc}  ./main_beam.exe $@  ${additional}
mpirun -np ${numProc}  ./main_cooler.exe $@  ${additional}
