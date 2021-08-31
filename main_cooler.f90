#include <petsc/finclude/petscsys.h>
program main_cooler


    use globals
    use fem_order2
    use petscsys
    use common_utils
    implicit none
    integer ierr, rank
    real(8) start, end

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); 
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr); 
    FILEINP = "./mesh/cooler.neu"
    call initializeStatus
    call initializeLib
    call InitializeConstitution

    if(rank == 0 )then
        call readgfile
        call ReducePoints
        call getVolumes
        call output_plt_mesh("./out/cooler_0_mesh.plt", "goodstart")
        call output_plt_scalar("./out/cooler_VOL.plt", "goodstart",cell_volumes,"cell_volume", .true.)
    endif

    call SetUpPartition

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
        bcValueTher(1) = 100.0_8
        bcTypeTher(2) = 1
        bcValueTher(2) = -10.0_8
        ! bcValueTher(2) = -1.0e-1_8*15
        ! bcValueTher2(2) = 1.0e-1_8

    endif
    call SetUpThermalBC_BLOCKED
    call SetUpThermal_InitializeObjects
    call SetUpThermalPara
    call SolveThermal_Initialize
    call SolveThermal
    !call output_plt_thermal("./out/cooler_ther.plt", "goodstart")

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
    call SetUpElasticBC_BLOCKED
    call SetUpElasticity_InitializeObjects
    call SetUpElasticityPara
    call SolveElasticity_Initialize
    call SolveElasticity
    !call output_plt_elasticity("./out2_elas.plt", "goodstart")
    call GetElasticityUGradient
    call GetStrainStress
    call output_plt_elasticity_all("./out/cooler_elas_all.plt", "goodstart")

    print*,rank,'done'
    if(rank == 0)then
        print*,'TimeElapsed:',MPI_Wtime()-start
    endif
    call PetscFinalize(ierr)
    CHKERRA(ierr)
end program
