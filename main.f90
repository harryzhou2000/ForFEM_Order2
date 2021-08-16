program main
#include <petsc/finclude/petscsys.h>

    use globals
    use fem_order2
    use petscsys
    use common_utils
    implicit none
    integer ierr, rank
    real(8) start, end

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); 
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr); 
    FILEINP = "./mark2_external.neu"
    call initializeStatus
    call initializeLib
    call InitializeConstitution

    if(rank == 0 )then
        call readgfile
        call ReducePoints
        call getVolumes
        call output_plt_mesh("./out2.plt", "goodstart")
        call output_plt_scalar("./out2_data1.plt", "goodstart",cell_volumes,"cell_volume", .true.)
    endif

    call SetUpPartition

    start = MPI_Wtime()

    !!!!!THERMAL SOLVE
    if(rank == 0) then
        if(NBSETS<6) then
            print*,'nbset to small, need 6 at least'
        endif
        allocate(bcValueTher(NBSETS))
        allocate(bcTypeTher(NBSETS))
        allocate(bcValueTher2(NBSETS))
        bcValueTher2 = 0.000_8 ! h
        bcValueTher = 1.0_8  ! h*phi_b
        bcTypeTher = 1
        bcValueTher(1) = 1
        bcValueTher(2) = -1.0e-1_8*15
        bcValueTher2(2) = 1.0e-1_8

        bcTypeTher(3) = 1
        bcValueTher(3) = 0
        bcTypeTher(4) = 1
        bcValueTher(4) = 0
        bcValueTher(5) = 0
        bcTypeTher(6) = 0
        bcValueTher(6) = 1.0_8
    endif
    call SetUpThermalBC_BLOCKED
    call SetUpThermal_InitializeObjects
    call SetUpThermalPara
    call SolveThermal_Initialize
    call SolveThermal
    !call output_plt_thermal("./out2_ther.plt", "goodstart")

    !!!!!ELASTIC SOLVE
    if(rank == 0) then
        if(NBSETS<6) then
            print*,'nbset to small, need 6 at least'
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
        bcValueElas(1*3-1) = myMinusNan()
        bcValueElas(1*3-0) = myMinusNan()
        ! 2
        bcTypeElas(2) = 0 ! 1 to use the mat below
        bcValueElas(2*3-2:2*3) = 0.0_8
        bcValueElas(2*3-2) = 0e-3_8
        bcValueElas2(2*9-8:2*9) = (/50.0000000000000,50.0000000000000,0.,50.0000000000000,50.0000000000000,0.,0.,0.,0./)
        ! 3 4
        bcTypeElas(3) = 1
        bcValueElas(3*3-2:3*3) = 0
        bcTypeElas(4) = 1
        bcValueElas(4*3-2:4*3) = 0
        ! 5
        bcTypeElas(5) = 0
        bcValueElas(5*3) = 0
        bcValueElas(5*3-2:5*3-1) = myMinusNan()
        ! 6
        bcTypeElas(6) = 1
        bcValueElas(6*3-2:6*3) = 0.0_8
    endif
    call SetUpElasticBC_BLOCKED
    call SetUpElasticity_InitializeObjects
    call SetUpElasticityPara
    call SolveElasticity_Initialize
    call SolveElasticity
    !call output_plt_elasticity("./out2_elas.plt", "goodstart")
    call GetElasticityUGradient
    call GetStrainStress
    call output_plt_elasticity_all("./out2_elas_all.plt", "goodstart")

    print*,rank,'done'
    if(rank == 0)then
        print*,'TimeElapsed:',MPI_Wtime()-start
    endif
    call PetscFinalize(ierr); 
    CHKERRA(ierr)
end program
