#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdmplex.h>
!#include <petsc/finclude/petscsection.h>

program test
    use petscsys
    use petscdmplex
    !use petscsection
    implicit none
    DM dm
    integer ierr
    PetscInt verts(3)
    integer mpirank, mpisize
    PetscInt p0, pe, f0, fe, v0, ve
    PetscSection s
    integer i

!!!!!!!!!!
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,ierr)
    if(mpisize < 1) then
        print *, 'not enough procs'
        stop
    endif
    CHKERRA(ierr)
    !call DMPlexCreate(MPI_COMM_WORLD, dm, ierr)
    call DMPlexCreateGmshFromFile(MPI_COMM_WORLD, 'Plane_FEM.msh', PETSC_FALSE, dm, ierr)
    ! No interpolate for gmsh so only gets cell->vertice
    call PetscSectionCreate(MPI_COMM_WORLD, s, ierr)
    CHKERRA(ierr)

    ! if(mpirank == 0) then
    !     call DMPlexSetChart(dm, 0, 8, ierr)
    !     call DMPlexSetConeSize(dm,5,3,ierr) !
    !     call DMPlexSetConeSize(dm,6,3,ierr) !
    !     call DMPlexSetConeSize(dm,7,3,ierr) !

    !     call DMSetUp(dm,ierr)
    !     call DMPlexSetCone(dm, 5, (/0,1,4/), ierr) !
    !     call DMPlexSetCone(dm, 6, (/1,3,4/), ierr) !
    !     call DMPlexSetCone(dm, 7, (/1,3,2/), ierr) !

    !     call DMPlexSymmetrize(dm, ierr)
    !     call DMPlexStratify(dm, ierr); 
    !     CHKERRA(ierr)
    ! else
    !     call DMPlexSetChart(dm, 8, 8, ierr)

    !     call DMSetUp(dm,ierr)

    !     call DMPlexSymmetrize(dm, ierr)
    !     call DMPlexStratify(dm, ierr); 
    !     CHKERRA(ierr)
    ! endif

    call DMPlexGetChart(dm, p0, pe, ierr)
    call PetscSectionSetChart(s, p0, pe, ierr)
    call DMPlexGetHeightStratum(dm, 0, f0, fe, ierr)
    call DMPlexGetHeightStratum(dm, 1, v0, ve, ierr)
    do i = f0, fe-1
        call PetscSectionSetDof(s, i, 0, ierr)
    end do
    do i = v0, ve-1
        call PetscSectionSetDof(s, i, 4, ierr)
    end do
    call PetscSectionSetUp(s, ierr)

    print*,'rank', mpirank,f0,fe
    CHKERRA(ierr)

    call PetscFinalize(ierr); 
    CHKERRA(ierr)
end program
