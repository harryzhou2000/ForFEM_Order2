#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>
#include <slepc/finclude/slepcsys.h>
program main_cooler

    use globals
    use fem_order2
    use petscsys
    use petscksp
    use slepcsys
    use common_utils
    use elastic_constitution
    implicit none
    integer ierr, rank, i
    real(8) start, end, localstart

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    FILEINP = "mesh/NOZZ_2_L.neu"
    call initializeStatus
    call initializeLib
    call set_const_thremal_constitution(397.0_8,390e6_8)
    call set_const_elastic_constitution(1.1e5_8,0.37_8,8900e-12_8)
    call set_expansion_properties(293.15_8, 2.4e-5_8)
    ! call set_copper_elastic_constitution
    if_dynamic_elas = .true.
    if_dynamic_ther = .true.

    if(rank == 0 )then
        call readgfile
        call ReducePoints
        call getVolumes
        
        call output_plt_mesh("./out/cooler_0_mesh.plt", "goodstart")
        print*,'here1'
        call output_plt_scalar("./out/cooler_VOL.plt", "goodstart",cell_volumes,"cell_volume", .true.)
        print*,'here2'
    endif

    localstart=MPI_Wtime()
    
    call SetUpPartition
    
    if(rank == 0)then
        print*,'Partition Time:',MPI_Wtime()-localstart
    endif

    start = MPI_Wtime()

    !!!!!THERMAL SOLVE
    if(rank == 0) then
        if(NBSETS<2) then
            print*,'nbset too small, need 2 at least'
        endif
        allocate(bcValueTher(NBSETS))
        allocate(bcTypeTher(NBSETS))
        allocate(bcValueTher2(NBSETS))
        bcValueTher2 = 0.000_8 ! h
        bcValueTher = 1.0_8  ! h*phi_b
        bcTypeTher(1) = 0
        bcValueTher(1) = 1073.15_8
        ! bcTypeTher(1) = 1
        ! bcValueTher(1) = 1e6*373.15_8
        ! bcValueTher2(1) = 1e6
        bcTypeTher(2) = 1
        bcValueTher(2) = 1e2*1e-3*473.15_8
        bcValueTher2(2) = 1e2*1e-3
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
    ! thermal mode
    ! call SolveThermalMode_Initialize
    ! call SolveThermalMode
    ! do i = 1,nsolvedEigenTher
    !     call output_plt_thermal_mode("./out/cooler_ther_mode", "goodstart", i-1)
    ! enddo
    call KSPDestroy(KSPther,ierr)
    call MatDestroy(Ather,ierr)
    call MatDestroy(Mther,ierr)

    !!!!!ELASTIC SOLVE
    if(rank == 0) then
        if(NBSETS<2) then
            print*,'nbset to small, need 2 at least'
        endif
        allocate(bcValueElas(NBSETS*3))
        allocate(bcTypeElas(NBSETS))
        allocate(bcValueElas2(NBSETS*9))
        bcValueElas2 = 0.000_8 ! h
        bcValueElas = 1.0_8  ! h*phi_b
        bcTypeElas = 1
        ! 1
        bcTypeElas(1) = 0
        bcValueElas(1*3-2) = 0.0_8
        bcValueElas(1*3-1) = 0.0_8
        bcValueElas(1*3-0) = 0.0_8
        ! 2
        bcTypeElas(2) = 1 ! 1 to use the mat below
        bcValueElas(2*3-2:2*3) = 0.0_8
        bcValueElas(2*3-2) = 0e-3_8
        bcValueElas2(2*9-8:2*9) = (/50.0000000000000,50.0000000000000,0.,50.0000000000000,50.0000000000000,0.,0.,0.,0./)
        bcValueElas2(2*9-8:2*9) = bcValueElas2(2*9-8:2*9) * 0.0_8
    endif
    localstart=MPI_Wtime()
    call SetUpElasticBC_BLOCKED
    call SetUpElasticity_InitializeObjects
    call SetUpElasticityPara
    if(rank == 0)then
        print*,'Elastic Setup Time:',MPI_Wtime()-localstart
    endif
    call SolveElasticity_Initialize
    localstart=MPI_Wtime()
    call SolveElasticity
    if(rank == 0)then
        print*,'Elastic Solve Time:',MPI_Wtime()-localstart
    endif
    ! elastic mode
    call SolveElasticMode_Initialize
    call SolveElasticMode
    do i = 1,nsolvedEigenElas
        call output_plt_elasticity_mode("./out/cooler_elas_mode", "goodstart", i-1)
    enddo

    ! strain stress
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
