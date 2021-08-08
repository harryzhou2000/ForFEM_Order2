program main
!#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscsys.h>

    use globals
    use fem_order2
    !use petscksp
    use petscsys
    use common_utils
    implicit none
    real :: a, b, result
    integer :: c; 
    integer, dimension(5) :: someInts
    character(80) :: outfile
    character(10) :: outtitle = "goodstart"
    real(8) :: D(3,3) = reshape((/1, 2, 3, 5 ,5, 6, 1, 8, 9/),(/3,3/))
    integer ierr, rank

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr); 
    CHKERRA(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr); 
    CHKERRA(ierr)
    if(rank == 0) then
        print*,"This is a PETSC_NULL_CHARACTER: ",PETSC_NULL_CHARACTER
    endif

    FILEINP = "./mark2_external.neu"
    outfile = "./out2.plt"
    call initializeStatus
    call initializeLib
    call InitializeConstitution

    if(rank == 0 )then
        call readgfile
        call ReducePoints
        call getVolumes
        call output_plt_mesh(outfile, outtitle)
        call output_plt_scalar("./out2_data1.plt", "goodstart",cell_volumes,"cell_volume", .true.)

        a = -12.5
        b = 15.0
        result = a + b
        c = a
        someInts(1) = 1; 
        someInts(2:4) = (/1,2,3/)

        ! print*,D
        ! print*, matmul(directInverse3x3(D),D)
    endif

    call SetUpPartition

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
        bcValueTher(2) = -14.0_8
        bcValueTher2(2) = 14.0_8

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
    call output_plt_thermal("./out2_ther.plt", "goodstart")


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
        bcValueElas(1*3-2:1*3) = 0.0_8
        ! 2
        bcTypeElas(2) = 1
        bcValueElas(2*3-2:2*3) = 0.0_8
        bcValueElas(2*3-2) = 1e-3_8
        bcValueElas2(2*9-8:2*9) = reshape(eye3x3() * 0.5,(/9/))
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
    call output_plt_elasticity("./out2_elas.plt", "goodstart")

    print*,rank,'done'
    call PetscFinalize(ierr); 
    CHKERRA(ierr)
end program
