#include <petsc/finclude/petscsys.h>
#include <slepc/finclude/slepcsys.h>
program main_cooler

    use globals
    use fem_order2
    use petscsys
    use slepcsys
    use common_utils
    use elastic_constitution
    implicit none
    integer ierr, rank, i
    real(8) start, end, localstart

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    FILEINP = "./mesh/cooler.neu"
    ! FILEINP = "mesh/COOLER_3_SMALL.neu"

    call initializeStatus
    call initializeLib
    ! 10-0.028
    !    1000

    call set_Unit_LTF(2.8e-3_8,1.0_8,1.0_8)
    ! call set_const_thremal_constitution(400.0_8, 390e6_8 )
    ! call set_const_elastic_constitution(1.1e5_8,0.37_8,8900e-12_8)
    ! call set_expansion_properties(293.15_8, 2.4e-5_8)
    call set_const_thremal_constitution(400.0_8, 390.0_8 )
    call set_const_elastic_constitution(1.1e11_8,0.37_8,8900.0_8)
    call set_expansion_properties(293.15_8, 2.4e-5_8)

    if_dynamic_elas = .true.
    if_dynamic_ther = .true.

    if(rank == 0 )then
        call readgfile
        call ReducePoints
        call getVolumes
        call output_plt_mesh("./out/cooler_0_mesh.plt", "goodstart")
        call output_plt_scalar("./out/cooler_VOL.plt", "goodstart",cell_volumes,"cell_volume", .true.)
    endif

    localstart=MPI_Wtime()
    call SetUpPartition
    if(rank == 0)then
        print*,'Partition Time:',MPI_Wtime()-localstart
    endif

    start = MPI_Wtime()

    !!!!!THERMAL SOLVE
    if(rank == 0) then
        call init_thermalBC_BLOCKED
        call set_thermalBC_BLOCKED(1, 0, 373.15_8, 0.0_8)
        call set_thermalBC_BLOCKED(2, 1, -10e3_8, 0.0_8)
    endif
    localstart=MPI_Wtime()
    call SetUpThermalBC_BLOCKED
    call SetUpThermal_InitializeObjects
    call SetUpThermalPara
    if(rank == 0)then
        print*,'Thermal Setup Time:',MPI_Wtime()-localstart
    endif
    call SolveThermal_Initialize
    localstart=MPI_Wtime()
    call SolveThermal
    if(rank == 0)then
        print*,'Thermal Solve Time:',MPI_Wtime()-localstart
    endif
    !call output_plt_thermal("./out/cooler_ther.plt", "goodstart")
    ! call SolveThermalMode_Initialize
    ! call SolveThermalMode
    ! do i = 1,nsolvedEigenTher
    !     call output_plt_thermal_mode("./out/cooler_ther_mode", "goodstart", i-1)
    ! enddo

    !!!!!ELASTIC SOLVE
    if(rank == 0) then
        call init_elasticBC_BLOCKED
        set_elasticBC_BLOCKED(1, 0, (/0.0_8,0.0_8,0.0_8/), reshape(eye3x3(),(/9/)) * 1e2)
        set_elasticBC_BLOCKED(2, 1, (/0.0_8,0.0_8,0.0_8/), &
                              (/50.0000000000000,50.0000000000000,0.,50.0000000000000,50.0000000000000,0.,0.,0.,0./)* 0.0_8)
    endif
    localstart=MPI_Wtime()
    call SetUpElasticBC_BLOCKED
    call SetUpElasticity_InitializeObjects
    call SetUpElasticityPara
    if(rank == 0)then
        print*,'Elastic Setup Time:',MPI_Wtime()-localstart
    endif
    ! call dumpAelasMelas('./out/cooler_elas_origin')
    ! call SolveElasticMode_Initialize
    ! ! elastic mode
    ! call SolveElasticMode
    ! do i = 1,nsolvedEigenElas
    !     call output_plt_elasticity_mode("./out/cooler_elas_mode", "goodstart", i-1)
    ! enddo

    call SolveElasticity_Initialize
    localstart=MPI_Wtime()
    call SolveElasticity
    if(rank == 0)then
        print*,'Elastic Solve Time:',MPI_Wtime()-localstart
    endif
    !call output_plt_elasticity("./out2_elas.plt", "goodstart")

    call GetElasticityUGradient
    call GetStrainStress
    call output_plt_elasticity_all("./out/cooler_elas_all.plt", "goodstart")

    print*,rank,'done'
    if(rank == 0)then
        print*,'TimeElapsed:',MPI_Wtime()-start
    endif

    call SlepcFinalize(ierr)
    call PetscFinalize(ierr)
    CHKERRA(ierr)
end program
