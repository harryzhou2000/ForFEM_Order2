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
    if(rank == 0) then
        allocate(bcValueTher(NBSETS))
        allocate(bcTypeTher(NBSETS))
        bcValueTher = 0.0_8
        bcTypeTher = 1
        bcTypeTher(1) = 0
        bcTypeTher(3) = 0
        bcTypeTher(4) = 0
        bcTypeTher(6) = 0
        bcValueTher(6) = 1.0_8
    endif
    call SetUpThermalBC_BLOCKED
    call SetUpThermal_InitializeObjects
    call SetUpThermal
    call SolveThermal_Initialize
    call SolveThermal
    call output_plt_thermal("./out2_ther.plt", "goodstart")

    if(rank == 0) then
        print *, 'The total is ', result , ' cis ', c
    endif

    print*,rank,'done'
    call PetscFinalize(ierr); 
    CHKERRA(ierr)
end program