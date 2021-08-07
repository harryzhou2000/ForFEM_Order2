#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdmplex.h>
#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscsection.h>

program test
    use petscsys
    use petscdmplex
    use petscvec
    !use petscsection
    implicit none
    DM dm
    integer ierr
    PetscInt verts(3)
    integer mpirank, mpisize
    PetscInt p0, pe, f0, fe, v0, ve
    PetscSection s
    integer i
    
    !!ghost vec test
    Vec Gvec
    PetscInt lo,hi
    PetscScalar val
    PetscScalar, pointer :: parray(:)

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

    
    if(mpirank == 0) then
        call VecCreateGhost(MPI_COMM_WORLD,5,PETSC_DECIDE,2,(/4,5/),Gvec, ierr)
        
    else
        call VecCreateGhost(MPI_COMM_WORLD,5,PETSC_DECIDE,2,(/1,0/),Gvec, ierr)
    endif

    call VecGetOwnershipRange(Gvec,lo, hi, ierr)
    do i = lo,hi-1
        val = i
        call VecSetValue(Gvec,i,val,INSERT_VALUES, ierr)
    enddo
    call VecAssemblyBegin(Gvec,ierr)
    call VecAssemblyEnd  (Gvec,ierr)
    call VecView(Gvec,PETSC_VIEWER_STDOUT_WORLD,ierr)

    call VecGhostUpdateBegin(Gvec,INSERT_VALUES,SCATTER_FORWARD, ierr);
    call VecGhostUpdateEnd  (Gvec,INSERT_VALUES,SCATTER_FORWARD, ierr);
    call VecGetArrayReadF90    (Gvec, parray ,ierr)
    print*, 'rank=',mpirank
    print*, parray(1:8)
    call VecRestoreArrayReadF90(Gvec, parray, ierr)


    call PetscFinalize(ierr); 
    CHKERRA(ierr)
end program
