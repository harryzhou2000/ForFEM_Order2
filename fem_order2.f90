!defining global variables for
!2nd order classsic isoparametric FEM
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscksp.h"

module fem_order2
    !#include "petsc/mpiuni/mpif.h"
    !#include "mpif.h"
    use petscvec
    use petscmat
    use petscsys
    use petscis
    use petscisdef
    use petscksp
    use mat_csr

    use common_utils

    !use mpi

    implicit none

    type fem_element
        integer(kind=4)           :: num_node ! M
        integer(kind=4)           :: num_intpoint ! R
        real(kind=8), allocatable :: coord_intpoint(:, :) !coords of integration point, 3xR
        real(kind=8), allocatable :: coord_intweight(:)  !weights of integration point, R
        real(kind=8), allocatable :: NImr(:, :) !N at intpoints M*R
        real(kind=8), allocatable :: dNIdLimr(:, :, :) !dNdL at intpoints 3*M*R

        integer(kind=4)              :: num_face       ! F
        integer(kind=4), allocatable :: face_size(:)   ! FS
        integer(kind=4), allocatable :: face_node(:,:) ! F*max(FS)
    end type

    type petscInt_vardimWrapper
        PetscInt, allocatable :: N(:)
    endtype

    type petscScalar_vardimWrapper
        PetscScalar, allocatable :: N(:)
    endtype

    type output_plt_value
        character(512) :: varName, varName2, varName3
        PetscScalar, pointer :: vardata(:)
        logical ifCell
    endtype

    type petscScalar_pointerWrapper
        PetscScalar, pointer :: P(:)
    endtype

    !!!! INPUT NEEDED:
    ! proc0:
    !   module global: CELL COORD bcread
    !   module fem_order2: bcValueT, bcValueT2, bcTypeT (where T = elas and ther)
    !

    !!!! other inputs:
    ! subrotine...

    !!!!!!!!!!!!!!!!!!!!!!!!
    !                      !
    !!!! partition data !!!!
    !                      !
    !!!!!!!!!!!!!!!!!!!!!!!!
    MatPartitioning partition
    Mat adjmat
    type(tIS) partitionIndex, rowsTwoSidedIndex, &
        partitionedNumberingIndex, partitionedNumberingIndexLoc

    type(tIs) partitionedNumberingIndexInversed, partitionedNumberingIndex3x

    PetscInt :: adjMaxWidth

    PetscInt :: localDOFs ! corresponding to rowsTwoSidedIndex, one dim per point
    PetscInt :: globalDOFs
    PetscInt indexLo ! just partitioned index range for local vectors
    PetscInt indexHi

    !!!!!!!!!!!!!!!!!!!!!!!!
    !                      !
    !!!! partition data !!!!
    !                      !
    !!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                           !
    !!!! cell partition data !!!!
    !                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! p0
    PetscInt, allocatable :: cellPartition(:)      !ncell
    PetscInt, allocatable :: cellLocalNumbering(:) !ncell
    ! all
    PetscInt nghost
    PetscInt, allocatable :: ghostingGlobal(:)     ! ghosting is based on local cells
    PetscInt, allocatable :: ghostingGlobal3x(:)   ! ghosting is based on local cells

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                          !
    !!!!cell partition data !!!!
    !                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !                        !
    !!!! partitioned mesh !!!!
    !                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    Vec localCoords ! ghosting is based on local cells
    type(csr) localCells
    PetscInt  nLocalCells
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !                        !
    !!!! partitioned mesh !!!!
    !                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                !
    !!!! solver and solution data !!!!
    !                                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    KSP  :: KSPelas
    Mat :: Aelas !stiffness for elasticity
    Vec :: Pelas !load Vector for elasticity
    Vec :: Uelas !solution Vector for elasticity
    Vec :: UGradientElas(9)
    Vec :: StrainElas(6)
    Vec :: StressElas(6)
    Vec :: VonMisesElas
    logical if_VMElas_alive, if_StrainElas_alive, if_StressElas_alive, &
        if_UGElas_alive
    !real(8), allocatable :: gatheredUelas(:)
    Vec dofFixElas
    PetscScalar, allocatable :: bcValueElas(:)  ! nbc * 3, 0p input
    PetscScalar, allocatable :: bcValueElas2(:) ! nbc * 9, 0p input
    PetscInt, allocatable :: bcTypeElas(:)      ! nbc    , 0p input
    !PetscInt, allocatable :: bcSetSizeElas(:)   ! nbc, calculated
    type(petscInt_vardimWrapper), allocatable :: bcDOFsElas(:) ! size is bc_related dofs
    type(petscScalar_vardimWrapper), allocatable :: bcVALsElas(:), bcVALsElas2(:)
    Vec dofFixElasDist

    KSP  :: KSPther
    Mat :: Ather !stiffness for thermal
    Vec :: Pther !load Vector for thermal
    Vec :: Uther !solution Vector for elasticity
    !real(8), allocatable :: gatheredUther(:)
    Vec dofFixTher
    PetscScalar, allocatable :: bcValueTher(:)  ! nbc, 0p input
    PetscScalar, allocatable :: bcValueTher2(:) ! nbc, 0p input
    PetscInt, allocatable :: bcTypeTher(:)      ! nbc, 0p input
    !PetscInt, allocatable :: bcSetSizeTher(:)   ! nbc, calculated
    type(petscInt_vardimWrapper), allocatable :: bcDOFsTher(:)
    type(petscScalar_vardimWrapper), allocatable :: bcVALsTher(:), bcVALsTher2(:)
    Vec dofFixTherDist
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                !
    !!!! solver and solution data !!!!
    !                                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!
    !                  !
    !!!! scatterers !!!!
    !                  !
    !!!!!!!!!!!!!!!!!!!!
    VecScatter procToZeroScatter
    VecScatter procToZeroScatter3x
    logical if_procToZeroScatter_alive, if_procToZeroScatter3x_alive
    !!!!!!!!!!!!!!!!!!!!
    !                  !
    !!!! scatterers !!!!
    !                  !
    !!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!
    !                      !
    !!!! auxiliary data !!!!
    !                      !
    !!!!!!!!!!!!!!!!!!!!!!!!
    !topology of mesh nodes, should be parallized or upgraded to PETSC's dmplex
    !serial for proc0
    type(csr), target :: AdjacencyCounter
    logical if_AdjacencyCounter_alive
    !node->nodes counter uncompressed
    !serial for proc0
    integer, allocatable :: AdjacencyNum0(:)
    !serial for proc0
    PetscInt, allocatable :: localCellSizes(:), localGhostSizes(:) ! size of comm
    !parallel
    PetscInt, allocatable :: localAdjacencyNum(:)
    PetscInt, allocatable :: ghostAdjacencyNum(:) !this is the ghosting in matrices
    integer, allocatable :: Ather_r_pos
    !serial for proc0
    real(8), allocatable :: cell_volumes(:)
    !!!!!!!!!!!!!!!!!!!!!!!!
    !                      !
    !!!! auxiliary data !!!!
    !                      !
    !!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !                        !
    !!!! FEM fundamentals !!!!
    !                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! constitutional:
    real(8) k_ther(3,3)
    integer, parameter:: nlib = 3, nlib_face = 2
    type(fem_element) elem_lib(nlib)
    type(fem_element) elem_lib_face(nlib_face)
    integer max_numnode,max_numintpoint
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !                        !
    !!!! FEM fundamentals !!!!
    !                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!

contains
    ! clear procToZeroScatter*
    subroutine initializeStatus
        PetscInt ierr, idim
        call VecScatterDestroy(procToZeroScatter, ierr)
        call VecScatterDestroy(procToZeroScatter3x, ierr)
        if_procToZeroScatter3x_alive = .false.
        if_procToZeroScatter_alive = .false.
        call Csr_DeleteMat(AdjacencyCounter)
        if_AdjacencyCounter_alive = .false.
        if(allocated(cell_volumes)) then
            deallocate(cell_volumes)
        endif

        if_UGElas_alive = .false.
        do idim = 1, 9
            call VecDestroy(UGradientElas(idim),ierr)
        enddo
        if_StrainElas_alive = .false.
        if_StressElas_alive = .false.
        do idim = 1,6
            call VecDestroy(StrainElas(idim), ierr)
            call VecDestroy(StressElas(idim), ierr)
        enddo
        if_VMElas_alive = .false.
        call VecDestroy(VonMisesElas, ierr)
    end subroutine

    subroutine initializeLib
        real(8) :: my_coord(3)

        !hammer 3rd tetra int
        real(8) :: tet_as(3,5)
        real(8) :: tet_ws(5)

        !hammer 5th tirangle int
        real(8), parameter :: tri_a1 = 0.0597158717_8
        real(8), parameter :: tri_b1 = 0.4701420641_8
        real(8), parameter :: tri_a2 = 0.7974269853_8
        real(8), parameter :: tri_b2 = 0.1012865073_8
        real(8), parameter :: tri_w0 = 0.225_8        * 0.5_8
        real(8), parameter :: tri_w1 = 0.1323941527_8 * 0.5_8
        real(8), parameter :: tri_w2 = 0.1259391805_8 * 0.5_8
        real(8) :: tri_as(2, 7), tri_ws(7)

        !gauss 4point int points
        real(8), parameter :: lin_a1 = 0.861136311594053_8
        real(8), parameter :: lin_a2 = 0.339981043584856_8
        real(8), parameter :: lin_w1 = 0.347854845137454_8
        real(8), parameter :: lin_w2 = 0.652145154862546_8
        real(8) :: lin_as(4), lin_ws(4)

        real(8) :: x, y, z !local coordinates

        integer i, j, k, l, intpoint_id

        !hammer tet3 organization
        tet_as(1,1) = 1.0_8/4.0_8
        tet_as(2,1) = 1.0_8/4.0_8
        tet_as(3,1) = 1.0_8/4.0_8
        tet_as(1,2) = 1.0_8/6.0_8
        tet_as(2,2) = 1.0_8/6.0_8
        tet_as(3,2) = 1.0_8/6.0_8
        tet_as(1,3) = 1.0_8/2.0_8
        tet_as(2,3) = 1.0_8/6.0_8
        tet_as(3,3) = 1.0_8/6.0_8
        tet_as(1,4) = 1.0_8/6.0_8
        tet_as(2,4) = 1.0_8/2.0_8
        tet_as(3,4) = 1.0_8/6.0_8
        tet_as(1,5) = 1.0_8/6.0_8
        tet_as(2,5) = 1.0_8/6.0_8
        tet_as(3,5) = 1.0_8/2.0_8
        tet_ws(1) = -4.0_8/ 5.0_8  /6.0_8
        tet_ws(2) =  9.0_8/20.0_8  /6.0_8
        tet_ws(3) =  9.0_8/20.0_8  /6.0_8
        tet_ws(4) =  9.0_8/20.0_8  /6.0_8
        tet_ws(5) =  9.0_8/20.0_8  /6.0_8

        !hammer tri5 points organization
        tri_as(1, 1) = 1.0/3.0; tri_as(2, 1) = 1.0/3.0
        tri_as(1, 2) = tri_a1; tri_as(2, 2) = tri_b1
        tri_as(1, 3) = tri_b1; tri_as(2, 3) = tri_a1
        tri_as(1, 4) = tri_b1; tri_as(2, 4) = tri_b1
        tri_as(1, 5) = tri_a2; tri_as(2, 5) = tri_b2
        tri_as(1, 6) = tri_b2; tri_as(2, 6) = tri_a2
        tri_as(1, 7) = tri_b2; tri_as(2, 7) = tri_b2
        tri_ws(1) = tri_w0
        tri_ws(2) = tri_w1; tri_ws(3) = tri_w1; tri_ws(4) = tri_w1
        tri_ws(5) = tri_w2; tri_ws(6) = tri_w2; tri_ws(7) = tri_w2

        !gauss points organization
        lin_as(1) = -lin_a1
        lin_as(2) = -lin_a2
        lin_as(3) =  lin_a2
        lin_as(4) =  lin_a1
        lin_ws(1) =  lin_w1
        lin_ws(2) =  lin_w2
        lin_ws(3) =  lin_w2
        lin_ws(4) =  lin_w1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                              !
        !    3D Element Definition     !
        !                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i = 1
        ! 15 node wedge
        elem_lib(i)%num_node = 15
        elem_lib(i)%num_intpoint = 7*4
        allocate (elem_lib(i)%coord_intpoint(3, elem_lib(i)%num_intpoint))
        allocate (elem_lib(i)%coord_intweight(elem_lib(i)%num_intpoint))

        do j = 1, 4
            do k = 1, 7
                my_coord(1) = tri_as(1, k)
                my_coord(2) = tri_as(2, k)
                my_coord(3) = lin_as(j)
                intpoint_id = k + (j - 1)*7
                elem_lib(i)%coord_intpoint(:, intpoint_id) = my_coord
                elem_lib(i)%coord_intweight(intpoint_id) = tri_ws(k)*lin_ws(j)
            end do
        end do
        ! shape functions
        allocate (elem_lib(i)%NImr(elem_lib(i)%num_node, elem_lib(i)%num_intpoint))
        allocate (elem_lib(i)%dNIdLimr(3, elem_lib(i)%num_node, elem_lib(i)%num_intpoint))
        do j = 1, elem_lib(i)%num_intpoint
            x = elem_lib(i)%coord_intpoint(1, j)
            y = elem_lib(i)%coord_intpoint(2, j)
            z = elem_lib(i)%coord_intpoint(3, j)

            !Ns
            elem_lib(i)%NImr(1,j)=-((z - 1.0)*(x + y - 1.0)*(2.0*x + 2.0*y + z))/2.0
            elem_lib(i)%NImr(2,j)=(x*(z - 1.0)*(z - 2.0*x + 2.0))/2.0
            elem_lib(i)%NImr(3,j)=(y*(z - 1.0)*(z - 2.0*y + 2.0))/2.0
            elem_lib(i)%NImr(4,j)=((z + 1.0)*(x + y - 1.0)*(2.0*x + 2.0*y - z))/2.0
            elem_lib(i)%NImr(5,j)=(x*(z + 1.0)*(2.0*x + z - 2.0))/2.0
            elem_lib(i)%NImr(6,j)=(y*(z + 1.0)*(2.0*y + z - 2.0))/2.0
            elem_lib(i)%NImr(7,j)=2.0*x*(z - 1.0)*(x + y - 1.0)
            elem_lib(i)%NImr(8,j)=-2.0*x*y*(z - 1.0)
            elem_lib(i)%NImr(9,j)=2.0*y*(z - 1.0)*(x + y - 1.0)
            elem_lib(i)%NImr(10,j)=-2.0*x*(z + 1.0)*(x + y - 1.0)
            elem_lib(i)%NImr(11,j)=2.0*x*y*(z + 1.0)
            elem_lib(i)%NImr(12,j)=-2.0*y*(z + 1.0)*(x + y - 1.0)
            elem_lib(i)%NImr(13,j)=(z**2.0 - 1.0)*(x + y - 1.0)
            elem_lib(i)%NImr(14,j)=-x*(z**2.0 - 1.0)
            elem_lib(i)%NImr(15,j)=-y*(z**2.0 - 1.0)

            !dirvatives of Ns
            elem_lib(i)%dNIdLimr(1,1,j)=-((z - 1.0)*(4.0*x + 4.0*y + z - 2.0))/2.0
            elem_lib(i)%dNIdLimr(1,2,j)=((z - 1.0)*(z - 4.0*x + 2.0))/2.0
            elem_lib(i)%dNIdLimr(1,3,j)=0.0
            elem_lib(i)%dNIdLimr(1,4,j)=((z + 1.0)*(4.0*x + 4.0*y - z - 2.0))/2.0
            elem_lib(i)%dNIdLimr(1,5,j)=((z + 1.0)*(4.0*x + z - 2.0))/2.0
            elem_lib(i)%dNIdLimr(1,6,j)=0.0
            elem_lib(i)%dNIdLimr(1,7,j)=2.0*(z - 1.0)*(2.0*x + y - 1.0)
            elem_lib(i)%dNIdLimr(1,8,j)=-2.0*y*(z - 1.0)
            elem_lib(i)%dNIdLimr(1,9,j)=2.0*y*(z - 1.0)
            elem_lib(i)%dNIdLimr(1,10,j)=-2.0*(z + 1.0)*(2.0*x + y - 1.0)
            elem_lib(i)%dNIdLimr(1,11,j)=2.0*y*(z + 1.0)
            elem_lib(i)%dNIdLimr(1,12,j)=-2.0*y*(z + 1.0)
            elem_lib(i)%dNIdLimr(1,13,j)=z**2.0 - 1.0
            elem_lib(i)%dNIdLimr(1,14,j)=1.0 - z**2.0
            elem_lib(i)%dNIdLimr(1,15,j)=0.0
            elem_lib(i)%dNIdLimr(2,1,j)=-((z - 1.0)*(4.0*x + 4.0*y + z - 2.0))/2.0
            elem_lib(i)%dNIdLimr(2,2,j)=0.0
            elem_lib(i)%dNIdLimr(2,3,j)=((z - 1.0)*(z - 4.0*y + 2.0))/2.0
            elem_lib(i)%dNIdLimr(2,4,j)=((z + 1.0)*(4.0*x + 4.0*y - z - 2.0))/2.0
            elem_lib(i)%dNIdLimr(2,5,j)=0.0
            elem_lib(i)%dNIdLimr(2,6,j)=((z + 1.0)*(4.0*y + z - 2.0))/2.0
            elem_lib(i)%dNIdLimr(2,7,j)=2.0*x*(z - 1.0)
            elem_lib(i)%dNIdLimr(2,8,j)=-2.0*x*(z - 1.0)
            elem_lib(i)%dNIdLimr(2,9,j)=2.0*(z - 1.0)*(x + 2.0*y - 1.0)
            elem_lib(i)%dNIdLimr(2,10,j)=-2.0*x*(z + 1.0)
            elem_lib(i)%dNIdLimr(2,11,j)=2.0*x*(z + 1.0)
            elem_lib(i)%dNIdLimr(2,12,j)=-2.0*(z + 1.0)*(x + 2.0*y - 1.0)
            elem_lib(i)%dNIdLimr(2,13,j)=z**2.0 - 1.0
            elem_lib(i)%dNIdLimr(2,14,j)=0.0
            elem_lib(i)%dNIdLimr(2,15,j)=1.0 - z**2.0
            elem_lib(i)%dNIdLimr(3,1,j)=-((x + y - 1.0)*(2.0*x + 2.0*y + 2.0*z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,2,j)=(x*(2.0*z - 2.0*x + 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,3,j)=(y*(2.0*z - 2.0*y + 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,4,j)=((x + y - 1.0)*(2.0*x + 2.0*y - 2.0*z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,5,j)=(x*(2.0*x + 2.0*z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,6,j)=(y*(2.0*y + 2.0*z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,7,j)=2.0*x*(x + y - 1.0)
            elem_lib(i)%dNIdLimr(3,8,j)=-2.0*x*y
            elem_lib(i)%dNIdLimr(3,9,j)=2.0*y*(x + y - 1.0)
            elem_lib(i)%dNIdLimr(3,10,j)=-2.0*x*(x + y - 1.0)
            elem_lib(i)%dNIdLimr(3,11,j)=2.0*x*y
            elem_lib(i)%dNIdLimr(3,12,j)=-2.0*y*(x + y - 1.0)
            elem_lib(i)%dNIdLimr(3,13,j)=2.0*z*(x + y - 1.0)
            elem_lib(i)%dNIdLimr(3,14,j)=-2.0*x*z
            elem_lib(i)%dNIdLimr(3,15,j)=-2.0*y*z
        end do

        !! set faces
        elem_lib(i)%num_face = 5
        allocate (elem_lib(i)%face_size(elem_lib(i)%num_face))
        elem_lib(i)%face_size = (/8,8,8,6,6/)
        allocate (elem_lib(i)%face_node(elem_lib(i)%num_face, 8))
        elem_lib(i)%face_node(1, 1:elem_lib(i)%face_size(1)) = (/1,7,2,14,5,10,4,13/)
        elem_lib(i)%face_node(2, 1:elem_lib(i)%face_size(2)) = (/2,8,3,15,6,11,5,14/)
        elem_lib(i)%face_node(3, 1:elem_lib(i)%face_size(3)) = (/3,9,1,13,4,12,6,15/)
        elem_lib(i)%face_node(4, 1:elem_lib(i)%face_size(4)) = (/1,9,3,8,2,7/)
        elem_lib(i)%face_node(5, 1:elem_lib(i)%face_size(5)) = (/4,10,5,11,6,12/)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i = 2
        ! 20 node brick
        ! intpoints and wight
        elem_lib(i)%num_node = 20
        elem_lib(i)%num_intpoint = 4*4*4
        allocate (elem_lib(i)%coord_intpoint(3, elem_lib(i)%num_intpoint))
        allocate (elem_lib(i)%coord_intweight(elem_lib(i)%num_intpoint))
        do j = 1, 4
            do k = 1, 4
                do l = 1, 4
                    my_coord(1) = lin_as(j)
                    my_coord(2) = lin_as(k)
                    my_coord(3) = lin_as(l)
                    intpoint_id = l + ((k - 1) + (j - 1)*4)*4
                    elem_lib(i)%coord_intpoint(:, intpoint_id) = my_coord
                    elem_lib(i)%coord_intweight(intpoint_id) = lin_ws(j)*lin_ws(k)*lin_ws(l)
                end do
            end do
        end do
        ! shape functions
        allocate (elem_lib(i)%NImr(elem_lib(i)%num_node, elem_lib(i)%num_intpoint))
        allocate (elem_lib(i)%dNIdLimr(3, elem_lib(i)%num_node, elem_lib(i)%num_intpoint))
        do j = 1, elem_lib(i)%num_intpoint
            x = elem_lib(i)%coord_intpoint(1, j)
            y = elem_lib(i)%coord_intpoint(2, j)
            z = elem_lib(i)%coord_intpoint(3, j)

            !Ns
            elem_lib(i)%NImr(1,j)=((x - 1.0)*(y - 1.0)*(z - 1.0)*(x + y + z + 2.0))/8.0
            elem_lib(i)%NImr(2,j)=-((x + 1.0)*(y - 1.0)*(z - 1.0)*(y - x + z + 2.0))/8.0
            elem_lib(i)%NImr(3,j)=-((x + 1.0)*(y + 1.0)*(z - 1.0)*(x + y - z - 2.0))/8.0
            elem_lib(i)%NImr(4,j)=-((x - 1.0)*(y + 1.0)*(z - 1.0)*(x - y + z + 2.0))/8.0
            elem_lib(i)%NImr(5,j)=-((x - 1.0)*(y - 1.0)*(z + 1.0)*(x + y - z + 2.0))/8.0
            elem_lib(i)%NImr(6,j)=-((x + 1.0)*(y - 1.0)*(z + 1.0)*(x - y + z - 2.0))/8.0
            elem_lib(i)%NImr(7,j)=((x + 1.0)*(y + 1.0)*(z + 1.0)*(x + y + z - 2.0))/8.0
            elem_lib(i)%NImr(8,j)=((x - 1.0)*(y + 1.0)*(z + 1.0)*(x - y - z + 2.0))/8.0
            elem_lib(i)%NImr(9,j)=-((x**2.0 - 1.0)*(y - 1.0)*(z - 1.0))/4.0
            elem_lib(i)%NImr(10,j)=((y**2.0 - 1.0)*(x + 1.0)*(z - 1.0))/4.0
            elem_lib(i)%NImr(11,j)=((x**2.0 - 1.0)*(y + 1.0)*(z - 1.0))/4.0
            elem_lib(i)%NImr(12,j)=-((y**2.0 - 1.0)*(x - 1.0)*(z - 1.0))/4.0
            elem_lib(i)%NImr(13,j)=((x**2.0 - 1.0)*(y - 1.0)*(z + 1.0))/4.0
            elem_lib(i)%NImr(14,j)=-((y**2.0 - 1.0)*(x + 1.0)*(z + 1.0))/4.0
            elem_lib(i)%NImr(15,j)=-((x**2.0 - 1.0)*(y + 1.0)*(z + 1.0))/4.0
            elem_lib(i)%NImr(16,j)=((y**2.0 - 1.0)*(x - 1.0)*(z + 1.0))/4.0
            elem_lib(i)%NImr(17,j)=-((z**2.0 - 1.0)*(x - 1.0)*(y - 1.0))/4.0
            elem_lib(i)%NImr(18,j)=((z**2.0 - 1.0)*(x + 1.0)*(y - 1.0))/4.0
            elem_lib(i)%NImr(19,j)=-((z**2.0 - 1.0)*(x + 1.0)*(y + 1.0))/4.0
            elem_lib(i)%NImr(20,j)=((z**2.0 - 1.0)*(x - 1.0)*(y + 1.0))/4.0

            !dirvatives of Ns
            elem_lib(i)%dNIdLimr(1,1,j)=((y - 1.0)*(z - 1.0)*(2.0*x + y + z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,2,j)=-((y - 1.0)*(z - 1.0)*(y - 2.0*x + z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,3,j)=-((y + 1.0)*(z - 1.0)*(2.0*x + y - z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,4,j)=-((y + 1.0)*(z - 1.0)*(2.0*x - y + z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,5,j)=-((y - 1.0)*(z + 1.0)*(2.0*x + y - z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,6,j)=-((y - 1.0)*(z + 1.0)*(2.0*x - y + z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,7,j)=((y + 1.0)*(z + 1.0)*(2.0*x + y + z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,8,j)=((y + 1.0)*(z + 1.0)*(2.0*x - y - z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(1,9,j)=-(x*(y - 1.0)*(z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(1,10,j)=((y**2.0 - 1.0)*(z - 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,11,j)=(x*(y + 1.0)*(z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(1,12,j)=-((y**2.0 - 1.0)*(z - 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,13,j)=(x*(y - 1.0)*(z + 1.0))/2.0
            elem_lib(i)%dNIdLimr(1,14,j)=-((y**2.0 - 1.0)*(z + 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,15,j)=-(x*(y + 1.0)*(z + 1.0))/2.0
            elem_lib(i)%dNIdLimr(1,16,j)=((y**2.0 - 1.0)*(z + 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,17,j)=-((z**2.0 - 1.0)*(y - 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,18,j)=((z**2.0 - 1.0)*(y - 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,19,j)=-((z**2.0 - 1.0)*(y + 1.0))/4.0
            elem_lib(i)%dNIdLimr(1,20,j)=((z**2.0 - 1.0)*(y + 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,1,j)=((x - 1.0)*(z - 1.0)*(x + 2.0*y + z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,2,j)=-((x + 1.0)*(z - 1.0)*(2.0*y - x + z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,3,j)=-((x + 1.0)*(z - 1.0)*(x + 2.0*y - z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,4,j)=-((x - 1.0)*(z - 1.0)*(x - 2.0*y + z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,5,j)=-((x - 1.0)*(z + 1.0)*(x + 2.0*y - z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,6,j)=-((x + 1.0)*(z + 1.0)*(x - 2.0*y + z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,7,j)=((x + 1.0)*(z + 1.0)*(x + 2.0*y + z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,8,j)=((x - 1.0)*(z + 1.0)*(x - 2.0*y - z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(2,9,j)=-((x**2.0 - 1.0)*(z - 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,10,j)=(y*(x + 1.0)*(z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(2,11,j)=((x**2.0 - 1.0)*(z - 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,12,j)=-(y*(x - 1.0)*(z - 1.0))/2.0
            elem_lib(i)%dNIdLimr(2,13,j)=((x**2.0 - 1.0)*(z + 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,14,j)=-(y*(x + 1.0)*(z + 1.0))/2.0
            elem_lib(i)%dNIdLimr(2,15,j)=-((x**2.0 - 1.0)*(z + 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,16,j)=(y*(x - 1.0)*(z + 1.0))/2.0
            elem_lib(i)%dNIdLimr(2,17,j)=-((z**2.0 - 1.0)*(x - 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,18,j)=((z**2.0 - 1.0)*(x + 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,19,j)=-((z**2.0 - 1.0)*(x + 1.0))/4.0
            elem_lib(i)%dNIdLimr(2,20,j)=((z**2.0 - 1.0)*(x - 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,1,j)=((x - 1.0)*(y - 1.0)*(x + y + 2.0*z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,2,j)=-((x + 1.0)*(y - 1.0)*(y - x + 2.0*z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,3,j)=-((x + 1.0)*(y + 1.0)*(x + y - 2.0*z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,4,j)=-((x - 1.0)*(y + 1.0)*(x - y + 2.0*z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,5,j)=-((x - 1.0)*(y - 1.0)*(x + y - 2.0*z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,6,j)=-((x + 1.0)*(y - 1.0)*(x - y + 2.0*z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,7,j)=((x + 1.0)*(y + 1.0)*(x + y + 2.0*z - 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,8,j)=((x - 1.0)*(y + 1.0)*(x - y - 2.0*z + 1.0))/8.0
            elem_lib(i)%dNIdLimr(3,9,j)=-((x**2.0 - 1.0)*(y - 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,10,j)=((y**2.0 - 1.0)*(x + 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,11,j)=((x**2.0 - 1.0)*(y + 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,12,j)=-((y**2.0 - 1.0)*(x - 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,13,j)=((x**2.0 - 1.0)*(y - 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,14,j)=-((y**2.0 - 1.0)*(x + 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,15,j)=-((x**2.0 - 1.0)*(y + 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,16,j)=((y**2.0 - 1.0)*(x - 1.0))/4.0
            elem_lib(i)%dNIdLimr(3,17,j)=-(z*(x - 1.0)*(y - 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,18,j)=(z*(x + 1.0)*(y - 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,19,j)=-(z*(x + 1.0)*(y + 1.0))/2.0
            elem_lib(i)%dNIdLimr(3,20,j)=(z*(x - 1.0)*(y + 1.0))/2.0
        end do

        !! set faces
        elem_lib(i)%num_face = 6
        allocate (elem_lib(i)%face_size(elem_lib(i)%num_face))
        elem_lib(i)%face_size = (/8,8,8,8,8,8/)
        allocate (elem_lib(i)%face_node(elem_lib(i)%num_face, 8))
        elem_lib(i)%face_node(1, 1:elem_lib(i)%face_size(1)) = (/1,12,4,11,3,10,2,9/)
        elem_lib(i)%face_node(2, 1:elem_lib(i)%face_size(2)) = (/4,20,8,15,7,19,3,11/)
        elem_lib(i)%face_node(3, 1:elem_lib(i)%face_size(3)) = (/8,16,5,13,6,14,7,15/)
        elem_lib(i)%face_node(4, 1:elem_lib(i)%face_size(4)) = (/5,17,1,9,2,18,6,13/)
        elem_lib(i)%face_node(5, 1:elem_lib(i)%face_size(5)) = (/4,12,1,17,5,16,8,20/)
        elem_lib(i)%face_node(6, 1:elem_lib(i)%face_size(6)) = (/2,10,3,19,7,14,6,18/)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i = 3
        ! 10 node tetra
        ! intpoints and wight
        elem_lib(i)%num_node = 10
        elem_lib(i)%num_intpoint = 5
        allocate (elem_lib(i)%coord_intpoint(3, elem_lib(i)%num_intpoint))
        allocate (elem_lib(i)%coord_intweight(elem_lib(i)%num_intpoint))
        elem_lib(i)%coord_intpoint = tet_as; 
        elem_lib(i)%coord_intweight = tet_ws
        ! shape functions
        allocate (elem_lib(i)%NImr(elem_lib(i)%num_node, elem_lib(i)%num_intpoint))
        allocate (elem_lib(i)%dNIdLimr(3, elem_lib(i)%num_node, elem_lib(i)%num_intpoint))
        do j = 1, elem_lib(i)%num_intpoint
            x = elem_lib(i)%coord_intpoint(1, j)
            y = elem_lib(i)%coord_intpoint(2, j)
            z = elem_lib(i)%coord_intpoint(3, j)

            !Ns
            elem_lib(i)%NImr(1,j)=(x + y + z - 1.0)*(2.0*x + 2.0*y + 2.0*z - 1.0)
            elem_lib(i)%NImr(2,j)=z*(2.0*z - 1.0)
            elem_lib(i)%NImr(3,j)=y*(2.0*y - 1.0)
            elem_lib(i)%NImr(4,j)=x*(2.0*x - 1.0)
            elem_lib(i)%NImr(5,j)=-4.0*z*(x + y + z - 1.0)
            elem_lib(i)%NImr(6,j)=4.0*y*z
            elem_lib(i)%NImr(7,j)=-4.0*y*(x + y + z - 1.0)
            elem_lib(i)%NImr(8,j)=-4.0*x*(x + y + z - 1.0)
            elem_lib(i)%NImr(9,j)=4.0*x*z
            elem_lib(i)%NImr(10,j)=4.0*x*y

            !dirvatives of Ns
            elem_lib(i)%dNIdLimr(1,1,j)=4.0*x + 4.0*y + 4.0*z - 3.0
            elem_lib(i)%dNIdLimr(1,2,j)=0.0
            elem_lib(i)%dNIdLimr(1,3,j)=0.0
            elem_lib(i)%dNIdLimr(1,4,j)=4.0*x - 1.0
            elem_lib(i)%dNIdLimr(1,5,j)=-4.0*z
            elem_lib(i)%dNIdLimr(1,6,j)=0.0
            elem_lib(i)%dNIdLimr(1,7,j)=-4.0*y
            elem_lib(i)%dNIdLimr(1,8,j)=4.0 - 4.0*y - 4.0*z - 8.0*x
            elem_lib(i)%dNIdLimr(1,9,j)=4.0*z
            elem_lib(i)%dNIdLimr(1,10,j)=4.0*y
            elem_lib(i)%dNIdLimr(2,1,j)=4.0*x + 4.0*y + 4.0*z - 3.0
            elem_lib(i)%dNIdLimr(2,2,j)=0.0
            elem_lib(i)%dNIdLimr(2,3,j)=4.0*y - 1.0
            elem_lib(i)%dNIdLimr(2,4,j)=0.0
            elem_lib(i)%dNIdLimr(2,5,j)=-4.0*z
            elem_lib(i)%dNIdLimr(2,6,j)=4.0*z
            elem_lib(i)%dNIdLimr(2,7,j)=4.0 - 8.0*y - 4.0*z - 4.0*x
            elem_lib(i)%dNIdLimr(2,8,j)=-4.0*x
            elem_lib(i)%dNIdLimr(2,9,j)=0.0
            elem_lib(i)%dNIdLimr(2,10,j)=4.0*x
            elem_lib(i)%dNIdLimr(3,1,j)=4.0*x + 4.0*y + 4.0*z - 3.0
            elem_lib(i)%dNIdLimr(3,2,j)=4.0*z - 1.0
            elem_lib(i)%dNIdLimr(3,3,j)=0.0
            elem_lib(i)%dNIdLimr(3,4,j)=0.0
            elem_lib(i)%dNIdLimr(3,5,j)=4.0 - 4.0*y - 8.0*z - 4.0*x
            elem_lib(i)%dNIdLimr(3,6,j)=4.0*y
            elem_lib(i)%dNIdLimr(3,7,j)=-4.0*y
            elem_lib(i)%dNIdLimr(3,8,j)=-4.0*x
            elem_lib(i)%dNIdLimr(3,9,j)=4.0*x
            elem_lib(i)%dNIdLimr(3,10,j)=0.0
        end do

        !! set faces
        elem_lib(i)%num_face = 4
        allocate (elem_lib(i)%face_size(elem_lib(i)%num_face))
        elem_lib(i)%face_size = (/6,6,6,6/)
        allocate (elem_lib(i)%face_node(elem_lib(i)%num_face, 8))
        elem_lib(i)%face_node(1, 1:elem_lib(i)%face_size(1)) = (/4,8,1,7,3,10/)
        elem_lib(i)%face_node(2, 1:elem_lib(i)%face_size(2)) = (/1,8,4,9,2,5/)
        elem_lib(i)%face_node(3, 1:elem_lib(i)%face_size(3)) = (/4,10,3,6,2,9/)
        elem_lib(i)%face_node(4, 1:elem_lib(i)%face_size(4)) = (/3,7,1,5,2,6/)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                              !
        !    2D Element Definition     !
        !                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i = 1
        ! 8 node rect
        ! intpoints and wight
        elem_lib_face(i)%num_node = 8
        elem_lib_face(i)%num_intpoint = 4*4
        allocate (elem_lib_face(i)%coord_intpoint(3, elem_lib_face(i)%num_intpoint))
        allocate (elem_lib_face(i)%coord_intweight(elem_lib_face(i)%num_intpoint))
        do j = 1,4
            do k = 1,4
                my_coord(1) = lin_as(j)
                my_coord(2) = lin_as(k)
                my_coord(3) = 0.0_8
                intpoint_id = k + (j - 1)*4
                elem_lib_face(i)%coord_intpoint(:, intpoint_id) = my_coord
                elem_lib_face(i)%coord_intweight(intpoint_id) = lin_ws(j)*lin_ws(k)
            enddo
        enddo
        ! shape functions
        allocate (elem_lib_face(i)%NImr(elem_lib_face(i)%num_node, elem_lib_face(i)%num_intpoint))
        allocate (elem_lib_face(i)%dNIdLimr(3, elem_lib_face(i)%num_node, elem_lib_face(i)%num_intpoint))
        do j = 1, elem_lib_face(i)%num_intpoint
            x = elem_lib_face(i)%coord_intpoint(1, j)
            y = elem_lib_face(i)%coord_intpoint(2, j)
            z = elem_lib_face(i)%coord_intpoint(3, j)

            !Ns
            elem_lib_face(i)%NImr(1,j)=-((x - 1.0)*(y - 1.0)*(x + y + 1.0))/4.0
            elem_lib_face(i)%NImr(2,j)=((x**2.0 - 1.0)*(y - 1.0))/2.0
            elem_lib_face(i)%NImr(3,j)=((x + 1.0)*(y - 1.0)*(y - x + 1.0))/4.0
            elem_lib_face(i)%NImr(4,j)=-((y**2.0 - 1.0)*(x + 1.0))/2.0
            elem_lib_face(i)%NImr(5,j)=((x + 1.0)*(y + 1.0)*(x + y - 1.0))/4.0
            elem_lib_face(i)%NImr(6,j)=-((x**2.0 - 1.0)*(y + 1.0))/2.0
            elem_lib_face(i)%NImr(7,j)=((x - 1.0)*(y + 1.0)*(x - y + 1.0))/4.0
            elem_lib_face(i)%NImr(8,j)=((y**2.0 - 1.0)*(x - 1.0))/2.0

            !dirvatives of Ns
            elem_lib_face(i)%dNIdLimr(1,1,j)=-((2.0*x + y)*(y - 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(1,2,j)=x*(y - 1.0)
            elem_lib_face(i)%dNIdLimr(1,3,j)=-((2.0*x - y)*(y - 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(1,4,j)=1.0/2.0 - y**2.0/2.0
            elem_lib_face(i)%dNIdLimr(1,5,j)=((2.0*x + y)*(y + 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(1,6,j)=-x*(y + 1.0)
            elem_lib_face(i)%dNIdLimr(1,7,j)=((2.0*x - y)*(y + 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(1,8,j)=y**2.0/2.0 - 1.0/2.0
            elem_lib_face(i)%dNIdLimr(2,1,j)=-((x + 2.0*y)*(x - 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(2,2,j)=x**2.0/2.0 - 1.0/2.0
            elem_lib_face(i)%dNIdLimr(2,3,j)=-((x - 2.0*y)*(x + 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(2,4,j)=-y*(x + 1.0)
            elem_lib_face(i)%dNIdLimr(2,5,j)=((x + 2.0*y)*(x + 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(2,6,j)=1.0/2.0 - x**2.0/2.0
            elem_lib_face(i)%dNIdLimr(2,7,j)=((x - 2.0*y)*(x - 1.0))/4.0
            elem_lib_face(i)%dNIdLimr(2,8,j)=y*(x - 1.0)
            elem_lib_face(i)%dNIdLimr(3,1,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,2,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,3,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,4,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,5,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,6,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,7,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,8,j)=0.0
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i = 2
        ! 6 node rect
        ! intpoints and wight
        elem_lib_face(i)%num_node = 6
        elem_lib_face(i)%num_intpoint = 7
        allocate (elem_lib_face(i)%coord_intpoint(3, elem_lib_face(i)%num_intpoint))
        allocate (elem_lib_face(i)%coord_intweight(elem_lib_face(i)%num_intpoint))
        elem_lib_face(i)%coord_intpoint(1:2, :) = tri_as
        elem_lib_face(i)%coord_intpoint(3,   :) = 0.0_8
        elem_lib_face(i)%coord_intweight = tri_ws
        ! shape functions
        allocate (elem_lib_face(i)%NImr(elem_lib_face(i)%num_node, elem_lib_face(i)%num_intpoint))
        allocate (elem_lib_face(i)%dNIdLimr(3, elem_lib_face(i)%num_node, elem_lib_face(i)%num_intpoint))
        do j = 1, elem_lib_face(i)%num_intpoint
            x = elem_lib_face(i)%coord_intpoint(1, j)
            y = elem_lib_face(i)%coord_intpoint(2, j)
            z = elem_lib_face(i)%coord_intpoint(3, j)

            !Ns
            elem_lib_face(i)%NImr(1,j)=(2.0*x + 2.0*y - 1.0)*(x + y - 1.0)
            elem_lib_face(i)%NImr(2,j)=-4.0*x*(x + y - 1.0)
            elem_lib_face(i)%NImr(3,j)=x*(2.0*x - 1.0)
            elem_lib_face(i)%NImr(4,j)=4.0*x*y
            elem_lib_face(i)%NImr(5,j)=y*(2.0*y - 1.0)
            elem_lib_face(i)%NImr(6,j)=-4.0*y*(x + y - 1.0)

            !dirvatives of Ns
            elem_lib_face(i)%dNIdLimr(1,1,j)=4.0*x + 4.0*y - 3.0
            elem_lib_face(i)%dNIdLimr(1,2,j)=4.0 - 4.0*y - 8.0*x
            elem_lib_face(i)%dNIdLimr(1,3,j)=4.0*x - 1.0
            elem_lib_face(i)%dNIdLimr(1,4,j)=4.0*y
            elem_lib_face(i)%dNIdLimr(1,5,j)=0.0
            elem_lib_face(i)%dNIdLimr(1,6,j)=-4.0*y
            elem_lib_face(i)%dNIdLimr(2,1,j)=4.0*x + 4.0*y - 3.0
            elem_lib_face(i)%dNIdLimr(2,2,j)=-4.0*x
            elem_lib_face(i)%dNIdLimr(2,3,j)=0.0
            elem_lib_face(i)%dNIdLimr(2,4,j)=4.0*x
            elem_lib_face(i)%dNIdLimr(2,5,j)=4.0*y - 1.0
            elem_lib_face(i)%dNIdLimr(2,6,j)=4.0 - 8.0*y - 4.0*x
            elem_lib_face(i)%dNIdLimr(3,1,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,2,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,3,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,4,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,5,j)=0.0
            elem_lib_face(i)%dNIdLimr(3,6,j)=0.0
        end do

        max_numintpoint = 0
        max_numnode = 0
        do i = 1,nlib
            max_numintpoint = max(max_numintpoint, elem_lib(i)%num_intpoint)
            max_numnode = max(max_numnode, elem_lib(i)%num_node)
        enddo
        do i = 1,nlib_face
            max_numintpoint = max(max_numintpoint, elem_lib_face(i)%num_intpoint)
            max_numnode = max(max_numnode, elem_lib_face(i)%num_node)
        enddo
    end subroutine

    subroutine InitializeConstitution
        k_ther = 0.0_8
        k_ther(1,1) = 1.0_8
        k_ther(2,2) = 1.0_8
        k_ther(3,3) = 1.0_8
        k_ther = k_ther * 110.0_8
    end subroutine

    !this is a tototally globals.mod-procedure
    subroutine ReducePoints
        use globals
        integer, allocatable :: newPointInd(:)
        logical, allocatable :: hasRefrence(:)
        real(8), allocatable :: newCoords(:,:)
        integer :: ncells, nverts, elem_id, num_node, i, j, newInd
        integer rank, siz, ierr
        !!!
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
        ncells = size(CELL)
        nverts = size(COORD,dim=2)
        if(rank == 0 .or. .true. )then
            allocate(hasRefrence(nverts),source=.false.)
            allocate(newPointInd(nverts),source=-1)
            do i = 1,ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                do j = 1,num_node
                    hasRefrence(CELL(i)%N(j)) = .true.
                enddo
            enddo
            ! print*,hasRefrence
            newInd = 0
            do i = 1, nverts
                if(hasRefrence(i)) then
                    newInd = newInd + 1
                    newPointInd(i) = newInd
                endif
            enddo
            allocate(newCoords(3,newInd))
            !wash the node indices of cells
            do i = 1, ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                do j = 1,num_node
                    CELL(i)%N(j) = newPointInd(CELL(i)%N(j))
                enddo
            enddo

            do i = 1, nverts
                if(hasRefrence(i)) then
                    newCoords(:,newPointInd(i)) = COORD(:,i)
                endif
            enddo
            deallocate(COORD)
            allocate(COORD(3,newInd))
            COORD = newCoords

            deallocate(newCoords)
            deallocate(newPointInd)
            deallocate(hasRefrence)
        endif
    end subroutine

    ! A seq routine, only 1 process should go
    subroutine getVolumes
        use globals
        integer i, ip, in, elem_id, num_intpoint, num_node
        real(8) :: elem_coord(27,3) !max num node is 27 ! warning
        real(8) :: Jacobi(3,3)
        real(8) :: detJacobi

        allocate(cell_volumes(size(CELL)))
        do i = 1, size(CELL)
            elem_id = getElemID(size(CELL(i)%N))
            cell_volumes(i) = 0.0_8
            num_intpoint = elem_lib(elem_id)%num_intpoint
            num_node = elem_lib(elem_id)%num_node
            if(num_node .ne. size(CELL(i)%N)) then
                print*,"getVolumes::Lib not matching allocated CELL(i)%N, i=",i
                stop
            end if
            do in = 1, num_node
                elem_coord(in,:) = (COORD(:,CELL(I)%N(in)))
            end do
            do ip = 1, num_intpoint
                Jacobi = matmul(elem_lib(elem_id)%dNIdLimr(:,:,ip),elem_coord(1:num_node,:))
                detJacobi = directDet3x3(Jacobi)
                cell_volumes(i) = cell_volumes(i) + detJacobi * elem_lib(elem_id)%coord_intweight(ip)
            end do
        end do
    end subroutine

    ! gadget
    function getElemID(CELLI) result(ID)
        integer, intent(in) :: CELLI
        integer :: ID
        select case((CELLI))
        case(15)
            ID = 1
        case(20)
            ID = 2
        case(10)
            ID = 3
        case default
            print *, "Error::getElemID::Current Num Node ",(CELLI)
            print *, "Not supported"
            stop
        end select
    end function

    ! gadget
    function getElemIDFace(FACEI) result(ID)
        integer, intent(in) :: FACEI
        integer :: ID
        select case((FACEI))
        case(8)
            ID = 1
        case(6)
            ID = 2
        case default
            print *, "Error::getFACEID::Current Num Node ",(FACEI)
            print *, "Not supported"
            stop
        end select
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                            !
    !!!!!!!!!!!Tecplot output functions!!!!!!!!!!!
    !                                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! A seq routine, only 1 process should go
    !https://tecplot.azureedge.net/products/360/current/360_data_format_guide.pdf
    subroutine output_plt_mesh(path, title)
        use globals
        implicit none
        character(*), intent(in) :: path
        integer(kind = 4) :: buffer_INT32
        character(512) :: headhead
        character(*),intent(in) :: title
        integer :: num_node_in_elem, i
        integer(kind = 4) :: aux_node(8)
        real(kind = 8) :: minCoord, maxCoord
        ! integer rank, ierr
        ! call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
        ! if (rank .ne. 0) then
        !     return
        ! end if

        if (.not. allocated(COORD)) then !using coord allocated to
            print *,"error::output_plt_mesh::coord data not allocated when calling output mesh"
            stop
        end if
        print *,"===Writing Mesh, num point=== ", size(COORD,2), " num elem ", size(CELL)

        close(IOUT2)
        open(IOUT2, file = path , Form = 'Unformatted' , Access = "Stream", STATUS = "Replace",  Action = "Write")
        headhead = "#!TDV112"
        write(IOUT2) headhead(1:8), 1_4, 1_4 !_ _ gridtype
        call write_binary_plt_string(title, IOUT2)
        write(IOUT2) 3_4
        headhead = "X"
        call write_binary_plt_string(headhead, IOUT2)
        headhead = "Y"
        call write_binary_plt_string(headhead, IOUT2)
        headhead = "Z"
        call write_binary_plt_string(headhead, IOUT2)

        !one zone
        write(IOUT2) 299.0_4 !zone
        headhead = "main_zone"
        call write_binary_plt_string(headhead, IOUT2)
        write(IOUT2) -1_4, -1_4, 0.0_8, -1_4, 5_4 !ParentZone, StrandID, SolutionTime, DZC, zoneType=FEBRICK
        write(IOUT2) 0_4, 0_4, 0_4 !specifyVarLocation=false
        buffer_INT32 = size(COORD,2)
        write(IOUT2) buffer_INT32 !numpoints
        buffer_INT32 = size(CELL)
        write(IOUT2) buffer_INT32, 0_4, 0_4, 0_4, 0 !numelems, IJKCelldim, no auxiliary data
        write(IOUT2) 357.0_4 !end of header

        write(IOUT2) 299.0_4 !zone
        write(IOUT2) 2_4, 2_4, 2_4 !print xyz as double
        write(IOUT2) 0_4, 0_4, -1_4 !nopassive nosharing nosharing
        minCoord = COORD(1, minloc(COORD(1,:),1))
        maxCoord = COORD(1, maxloc(COORD(1,:),1))
        write(IOUT2) minCoord, maxCoord
        minCoord = COORD(2, minloc(COORD(2,:),1))
        maxCoord = COORD(2, maxloc(COORD(2,:),1))
        write(IOUT2) minCoord, maxCoord
        minCoord = COORD(3, minloc(COORD(3,:),1))
        maxCoord = COORD(3, maxloc(COORD(3,:),1))
        write(IOUT2) minCoord, maxCoord
        write(IOUT2) COORD(1,:) !x
        write(IOUT2) COORD(2,:) !y
        write(IOUT2) COORD(3,:) !z

        do i = 1,size(CELL)
            num_node_in_elem = size(CELL(i)%N)
            select case(num_node_in_elem)
            case(15) ! 15 node wedge
                ! aux_node(1) = CELL(i)%N( 0+1)
                ! aux_node(2) = CELL(i)%N( 5+1)
                ! aux_node(3) = CELL(i)%N(14+1)
                ! aux_node(4) = CELL(i)%N( 9+1)
                ! aux_node(5) = CELL(i)%N( 2+1)
                ! aux_node(6) = CELL(i)%N( 2+1)
                ! aux_node(7) = CELL(i)%N(11+1)
                ! aux_node(8) = CELL(i)%N(11+1)
                aux_node(1:3) = CELL(i)%N(1:3)
                aux_node(4) = CELL(i)%N(3)
                aux_node(5:7) = CELL(i)%N(4:6)
                aux_node(8) = CELL(i)%N(6)

            case(20)
                ! aux_node(1) = CELL(i)%N( 0+1)
                ! aux_node(2) = CELL(i)%N(12+1)
                ! aux_node(3) = CELL(i)%N(14+1)
                ! aux_node(4) = CELL(i)%N( 2+1)
                ! aux_node(5) = CELL(i)%N( 5+1)
                ! aux_node(6) = CELL(i)%N(17+1)
                ! aux_node(7) = CELL(i)%N(19+1)
                ! aux_node(8) = CELL(i)%N( 7+1)
                aux_node = CELL(i)%N(1:8)

            case default
                print *, "Error::output_plt_mesh::Current Num Node ",num_node_in_elem
                print *, "Not supported"
                stop
            end select
            write(IOUT2) aux_node-1
            !print *, aux_node
        end do
        close(IOUT2)
        print *,"===Writing Mesh Successful==="
    end subroutine

    ! A seq routine, only 1 process should go
    ! TODO: adapt to PETSC_SCALAR
    subroutine output_plt_scalar(path, title, DATAin, DATAname, ifCell)
        use globals
        implicit none
        character(*), intent(in) :: path
        integer(kind = 4) :: buffer_INT32
        character(len=512) :: headhead
        character(*),intent(in) :: title
        !integer :: num_node_in_elem, i
        !integer(kind = 4) :: aux_node(8)
        real(kind = 8) :: minCoord, maxCoord
        real(kind = 8), intent(in) :: DATAin(:)
        character(*), intent(in) :: DATAname
        integer nmaxloc(1), nminloc(1)
        logical ifCell
        integer(4) CellInd
        ! integer rank, ierr
        ! call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
        ! if (rank .ne. 0) then
        !     return
        ! end if
        CellInd = 0_4
        if(ifCell) then
            CellInd = 1_4
        endif

        if (.not. allocated(COORD)) then !using coord allocated to
            print *,"error::output_plt_mesh::coord data not allocated when calling output mesh"
            stop
        end if
        print *,"===Writing Data === ",DATAname,", num point ", size(COORD,2), " num elem ", size(CELL)

        close(IOUT2)
        open(IOUT2, file = path , Form = 'Unformatted' , Access = "Stream", STATUS = "Replace",  Action = "Write")
        headhead = "#!TDV112"
        write(IOUT2) headhead(1:8), 1_4, 2_4 !_ _ solutiontype
        call write_binary_plt_string(title, IOUT2)
        write(IOUT2) 1_4 !only one data
        call write_binary_plt_string(DATAname, IOUT2)

        !one zone
        write(IOUT2) 299.0_4 !zone
        headhead = "main_zone"
        call write_binary_plt_string(headhead, IOUT2)
        write(IOUT2) -1_4, -1_4, 0.0_8, -1_4, 5_4 !ParentZone, StrandID, SolutionTime, DZC, zoneType=FEBRICK
        write(IOUT2) 1_4, CellInd !specifyVarLocation
        write(IOUT2) 0_4, 0_4
        buffer_INT32 = size(COORD,2)
        write(IOUT2) buffer_INT32 !numpoints
        buffer_INT32 = size(CELL)
        write(IOUT2) buffer_INT32, 0_4, 0_4, 0_4, 0 !numelems, IJKCelldim, no auxiliary data

        write(IOUT2) 357.0_4 !end of header

        write(IOUT2) 299.0_4 !zone
        write(IOUT2) 2_4 !print xyz as double
        write(IOUT2) 0_4, 0_4, -1_4 !nopassive nosharing nosharing
        nmaxloc = maxloc(DATAin)
        nminloc = minloc(DATAin)
        minCoord = DATAin(nminloc(1))
        maxCoord = DATAin(nmaxloc(1))
        write(IOUT2) minCoord, maxCoord
        write(IOUT2) DATAin
        if(ifCell) then
            if(size(DATAin) .ne. size(CELL)) then
                print*,"Error::output_plt_scalar::DATA size ",size(DATAin)," not equal to CELL ",size(CELL)
                stop
            end if
        else
            if(size(DATAin) .ne. size(COORD,2)) then
                print*,"Error::output_plt_scalar::DATA size ",size(DATAin)," not equal to VERT",size(COORD,2)
                stop
            end if
        endif

        close(IOUT2)
        print*,'===Writing Data done ==='
    end subroutine

    ! Input DATAin is {d1[1] d2[1] d3[1] d1[2] d2[2] d3[2] ...}
    ! TODO: adapt to PETSC_SCALAR
    subroutine output_plt_scalar3x(path, title, DATAin, DATAname1, DATAname2, DATAname3, ifCell)
        use globals
        implicit none
        character(*), intent(in) :: path
        integer(kind = 4) :: buffer_INT32
        character(len=512) :: headhead
        character(*),intent(in) :: title
        !integer :: num_node_in_elem, i
        !integer(kind = 4) :: aux_node(8)
        !real(kind = 8) :: minCoord, maxCoord
        real(kind = 8), intent(in) :: DATAin(:)
        character(*), intent(in) :: DATAname1, DATAname2, DATAname3
        integer  i
        logical ifCell
        integer(4) CellInd
        real(8) minmaxs(6)
        ! integer rank, ierr
        ! call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
        ! if (rank .ne. 0) then
        !     return
        ! end if
        CellInd = 0_4
        if(ifCell) then
            CellInd = 1_4
        endif

        if (.not. allocated(COORD)) then !using coord allocated to
            print *,"error::output_plt_mesh::coord data not allocated when calling output mesh"
            stop
        end if
        print *,"===Writing Data=== ",DATAname1,',', DATAname2,',', DATAname3,','&
            ,", num point ", size(COORD,2), " num elem ", size(CELL)

        close(IOUT2)
        open(IOUT2, file = path , Form = 'Unformatted' , Access = "Stream", STATUS = "Replace",  Action = "Write")
        headhead = "#!TDV112"
        write(IOUT2) headhead(1:8), 1_4, 2_4 !_ _ solutiontype
        call write_binary_plt_string(title, IOUT2)
        write(IOUT2) 3_4 !only 3 data
        call write_binary_plt_string(DATAname1, IOUT2)
        call write_binary_plt_string(DATAname2, IOUT2)
        call write_binary_plt_string(DATAname3, IOUT2)

        !one zone
        write(IOUT2) 299.0_4 !zone
        headhead = "main_zone"
        call write_binary_plt_string(headhead, IOUT2)
        write(IOUT2) -1_4, -1_4, 0.0_8, -1_4, 5_4 !ParentZone, StrandID, SolutionTime, DZC, zoneType=FEBRICK
        write(IOUT2) 1_4, CellInd, CellInd, CellInd  !specifyVarLocation
        write(IOUT2) 0_4, 0_4
        buffer_INT32 = size(COORD,2)
        write(IOUT2) buffer_INT32 !numpoints
        buffer_INT32 = size(CELL)
        write(IOUT2) buffer_INT32, 0_4, 0_4, 0_4, 0_4 !numelems, IJKCelldim, no auxiliary data

        write(IOUT2) 357.0_4 !end of header

        write(IOUT2) 299.0_4 !zone
        write(IOUT2) 2_4, 2_4, 2_4 !print xyz as double
        write(IOUT2) 0_4, 0_4, -1_4 !nopassive nosharing nosharing
        minmaxs = getMinMax3x(DATAin)
        write(IOUT2) minmaxs
        !print*,minmaxs
        do i = 1,(size(DATAin)/3)
            write(IOUT2) DATAin((i-1)*3+1)
        enddo
        do i = 1,(size(DATAin)/3)
            write(IOUT2) DATAin((i-1)*3+2)
        enddo
        do i = 1,(size(DATAin)/3)
            write(IOUT2) DATAin((i-1)*3+3)
        enddo
        ! write(IOUT2) DATAin
        if(ifCell) then
            if(size(DATAin) .ne. size(CELL) * 3) then
                print*,"Error::output_plt_scalar::DATA size ",size(DATAin)," not equal to CELL ",size(CELL)
                stop
            end if
        else
            if(size(DATAin) .ne. size(COORD,2) * 3) then
                print*,"Error::output_plt_scalar::DATA size ",size(DATAin)," not equal to VERT",size(COORD,2)
                stop
            end if
        endif

        close(IOUT2)
        print*,'===Writing Data done ==='
    end subroutine

    subroutine output_plt_scalars(path, title, dataset)
        use globals
        implicit none
        character(*), intent(in) :: path
        integer(kind = 4) :: buffer_INT32
        character(len=512) :: headhead
        character(*),intent(in) :: title
        type(output_plt_value), intent(in) ::dataset(:)
        integer  i,j
        logical, allocatable :: if3d(:)
        integer(4) PetscScalarKind, totalnumdata

        if(size(dataset)<=0)then
            return
        endif
        if(size(dataset(1)%vardata)<=0)then
            print*,'error::output_plt_scalars::dataset(1)%vardata missing'
            stop
        endif
        PetscScalarKind = kind(dataset(1)%vardata(1))
        if(PetscScalarKind .ne. 8)then
            print*,'error::output_plt_scalars::do not support data rather than double'
            stop
        endif
        if (.not. allocated(COORD)) then !using coord allocated to
            print *,"error::output_plt_scalars::coord data not allocated when calling output mesh"
            stop
        end if
        allocate(if3d(size(dataset)))
        totalnumdata = 0
        do j = 1,size(dataset)
            if(dataset(j)%ifCell) then
                if(size(dataset(j)%vardata) == size(CELL))then
                    if3d(j) = .false.
                    totalnumdata = totalnumdata + 1
                elseif(size(dataset(j)%vardata) == size(CELL) * 3) then
                    if3d(j) = .true.
                    totalnumdata = totalnumdata + 3
                else
                    print*, "error::output_plt_scalars::non-conforming size of data", size(CELL),&
                        'datasize: ', size(dataset(j)%vardata)
                    stop
                endif
            else
                if(size(dataset(j)%vardata) == size(COORD,2))then
                    if3d(j) = .false.
                    totalnumdata = totalnumdata + 1
                elseif(size(dataset(j)%vardata) == size(COORD,2) * 3) then
                    if3d(j) = .true.
                    totalnumdata = totalnumdata + 3
                else
                    print*, "error::output_plt_scalars::non-conforming size of data", size(COORD,2),&
                        'datasize: ', size(dataset(j)%vardata)
                    stop
                endif
            endif
        enddo
        print *,"===Writing Data Set=== size ",size(dataset)
        do j = 1,size(dataset)
            print*,trim(dataset(j)%varName)
            if(if3d(j))then
                print*,trim(dataset(j)%varName2)
                print*,trim(dataset(j)%varName3)
            endif
        enddo
        print *,"total", totalnumdata, ", num point ", size(COORD,2), " num elem ", size(CELL)

        close(IOUT2)
        open(IOUT2, file = path , Form = 'Unformatted' , Access = "Stream", STATUS = "Replace",  Action = "Write")
        headhead = "#!TDV112"
        write(IOUT2) headhead(1:8), 1_4, 2_4 !_ _ solutiontype
        call write_binary_plt_string(title, IOUT2)
        write(IOUT2) totalnumdata
        do j = 1,size(dataset)
            call write_binary_plt_string(dataset(j)%varName, IOUT2)
            if(if3d(j))then
                call write_binary_plt_string(dataset(j)%varName2, IOUT2)
                call write_binary_plt_string(dataset(j)%varName3, IOUT2)
            endif
        enddo

        !one zone
        write(IOUT2) 299.0_4 !zone
        headhead = "main_zone"
        call write_binary_plt_string(headhead, IOUT2)
        write(IOUT2) -1_4, -1_4, 0.0_8, -1_4, 5_4 !ParentZone, StrandID, SolutionTime, DZC, zoneType=FEBRICK
        write(IOUT2) 1_4  !specifyVarLocation
        do j = 1, size(dataset)
            if(if3d(j)) then
                if(dataset(j)%ifCell)then
                    write(IOUT2)  1_4, 1_4, 1_4
                else
                    write(IOUT2)  0_4, 0_4, 0_4
                endif
            else
                if(dataset(j)%ifCell)then
                    write(IOUT2)  1_4
                else
                    write(IOUT2)  0_4
                endif
            endif
        enddo
        write(IOUT2) 0_4, 0_4
        buffer_INT32 = size(COORD,2)
        write(IOUT2) buffer_INT32 !numpoints
        buffer_INT32 = size(CELL)
        write(IOUT2) buffer_INT32, 0_4, 0_4, 0_4, 0_4 !numelems, IJKCelldim, no auxiliary data

        write(IOUT2) 357.0_4 !end of header

        write(IOUT2) 299.0_4 !zone
        do j = 1, totalnumdata
            write(IOUT2) 2_4
        enddo
        write(IOUT2) 0_4, 0_4, -1_4 !nopassive nosharing nosharing
        do j = 1, size(dataset)
            if(if3d(j))then
                write(IOUT2) getMinMax3x(dataset(j)%vardata)
            else
                write(IOUT2) getMinMax(dataset(j)%vardata)
            endif
        enddo

        do j = 1, size(dataset)
            if(if3d(j)) then
                do i = 1,(size(dataset(j)%vardata)/3)
                    write(IOUT2) dataset(j)%vardata((i-1)*3+1)
                enddo
                do i = 1,(size(dataset(j)%vardata)/3)
                    write(IOUT2) dataset(j)%vardata((i-1)*3+2)
                enddo
                do i = 1,(size(dataset(j)%vardata)/3)
                    write(IOUT2) dataset(j)%vardata((i-1)*3+3)
                enddo
            else
                write(IOUT2) dataset(j)%vardata
            endif
        enddo
        close(IOUT2)
        print *,"===Writing Data Successful==="
    end subroutine

    subroutine write_binary_plt_string(title,IOUT2)
        implicit none
        character(*),intent(in) :: title
        integer ::IOUT2 ,i
        integer(kind = 4) :: buffer_INT32
        do i = 1, LEN(title)
            buffer_INT32 = ICHAR(title(i:i))
            if((buffer_INT32 == 0) .or. (buffer_INT32 == 32)) then
                exit
            end if
            write(IOUT2) buffer_INT32
        end do
        write(IOUT2) 0_4
    end subroutine

    !!!!!!!!!!!!!!!!! subroutines invoking PETSC
    ! builds adjacencyCounter, adjacendyNum
    subroutine SetUpPartition
        use globals
        use linear_set

        integer ncells, nverts,  elem_id, num_intpoint, num_node, current_row_start, current_ivert,&
            next_row_start, new_row_end, current_row_start_old, next_row_start_old
        integer rank, siz, ierr
        integer i,j,ir
        integer fillBuffer(27)
        integer, pointer:: tobesorted(:)
        integer, allocatable:: oldCols(:), oldRows(:)
        PetscInt, pointer :: parray(:), parray2(:)
        Mat m2, m3
        Vec v1,v2
        PetscInt msizem, msizen
        integer loc(1)
        PetscInt maxWidth
        PetscScalar, allocatable::cellMat(:,:)
        PetscInt,allocatable:: cellDOFs(:)
        real(8), allocatable :: point_partition(:)

        PetscInt,allocatable::numberingExtendedVal(:)
        PetscScalar, allocatable::tempScalarizedWidth(:)
        PetscInt, allocatable::tempScalarizedInsertionIndex(:)
        Vec widthVecZero
        Vec widthVecDist
        PetscScalar, pointer :: pwidth(:)
        PetscInt currentProc, numOffProc

        !!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)

        if(rank==0) then
            ncells = size(CELL)
            nverts = size(COORD,dim=2)
            globalDOFs = nverts
        endif
        call MPI_Bcast(globalDOFs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        !!! Start Deriving AdjacencyCounter
        if(rank == 0)then
            if(allocated(AdjacencyNum0)) then
                deallocate(AdjacencyNum0)
            endif
            allocate(AdjacencyNum0(nverts),source = 0)
            do i = 1,ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                do j = 1,num_node
                    AdjacencyNum0(CELL(i)%N(j)) = AdjacencyNum0(CELL(i)%N(j)) + num_node-1
                enddo
            enddo
            call Csr_SetNROW(AdjacencyCounter, nverts)
            call Csr_AllocateRowStart(AdjacencyCounter)
            AdjacencyCounter%rowStart(0+1) = 0
            do i = 1,nverts
                AdjacencyCounter%rowStart(i+1) = AdjacencyCounter%rowStart(i) + AdjacencyNum0(i)
            enddo
            call Csr_SetNNZ(AdjacencyCounter, AdjacencyCounter%rowStart(nverts+1))
            call Csr_AllocateColumn(AdjacencyCounter)

            AdjacencyNum0 = 0
            do i = 1,ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                do j = 1,num_node
                    fillBuffer(1:j-1) = CELL(i)%N(1:j-1)
                    fillBuffer(j:num_node-1) = CELL(i)%N(j+1:num_node)
                    current_ivert = CELL(i)%N(j)
                    current_row_start = AdjacencyCounter%rowStart(current_ivert) + AdjacencyNum0(current_ivert)
                    AdjacencyNum0(current_ivert) = AdjacencyNum0(current_ivert) + num_node - 1
                    AdjacencyCounter%column(current_row_start+1:current_row_start+num_node-1) &
                        = fillBuffer(1:num_node-1) - 1
                enddo
            enddo
            ! print*,'prev\',AdjacencyNum0
            tobesorted => AdjacencyCounter%column
            do i = 1,nverts
                current_row_start = AdjacencyCounter%rowStart(i)
                next_row_start = AdjacencyCounter%rowStart(i+1)
                call int_mergeSort(tobesorted, current_row_start+1, next_row_start+1)
                call int_reduceSorted(tobesorted,current_row_start+1,  next_row_start+1, new_row_end)
                AdjacencyNum0(i) = new_row_end - current_row_start - 1
            enddo
            ! print*,'after\',AdjacencyNum0
            ! print*,'rstart\', AdjacencyCounter%rowstart
            ! print*,'data\' ,AdjacencyCounter%column
            allocate(oldCols(size(AdjacencyCounter%column)))
            oldCols = AdjacencyCounter%column
            allocate(oldRows(size(AdjacencyCounter%rowStart)))
            oldRows = AdjacencyCounter%rowStart

            do i = 1,nverts
                AdjacencyCounter%rowStart(i+1) = AdjacencyCounter%rowStart(i) + AdjacencyNum0(i)
            enddo
            call Csr_SetNNZ(AdjacencyCounter, AdjacencyCounter%rowStart(nverts+1))
            call Csr_AllocateColumn(AdjacencyCounter)
            do i = 1,nverts
                current_row_start = AdjacencyCounter%rowStart(i)
                next_row_start    = AdjacencyCounter%rowStart(i+1)
                current_row_start_old = oldRows(i)
                AdjacencyCounter%column(current_row_start+1: next_row_start+1) = &
                    oldCols(current_row_start_old+1: AdjacencyNum0(i) + current_row_start_old)
            enddo
            ! print*,'rstart\', AdjacencyCounter%rowstart
            ! print*,'data\' ,AdjacencyCounter%column

            deallocate(oldCols)
            deallocate(oldRows)
        else
            call Csr_SetNROW(AdjacencyCounter,0)
            call Csr_AllocateRowStart(AdjacencyCounter)
            AdjacencyCounter%RowStart(1) = 0
            call Csr_SetNNZ(AdjacencyCounter,0)
            call Csr_AllocateColumn(AdjacencyCounter)
        endif

        call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INDEX, ierr)

        !!! Set MatAdj
        if(rank==0) then
            call MatCreateMPIAdj(MPI_COMM_WORLD,globalDOFs,globalDOFs,AdjacencyCounter%rowstart,AdjacencyCounter%column,&
                                 PETSC_NULL_INTEGER,adjMat, ierr)
        else
            call MatCreateMPIAdj(MPI_COMM_WORLD,0,         globalDOFs,AdjacencyCounter%rowstart,AdjacencyCounter%column,&
                                 PETSC_NULL_INTEGER,adjMat, ierr)
        endif
        !call MatView(adjMat,PETSC_VIEWER_STDOUT_WORLD,ierr)

        !!! Do Partition
        call MatPartitioningCreate(MPI_COMM_WORLD, partition, ierr)
        call MatPartitioningSetAdjacency(partition , adjmat, ierr)
        call MatPartitioningSetFromOptions(partition, ierr)
        call ISCreate(MPI_COMM_WORLD,partitionIndex,ierr)
        call MatPartitioningApply(partition,partitionIndex, ierr)
        call MatPartitioningDestroy(partition, ierr) !!!

        !!! get Two Sided
        call ISCreate(MPI_COMM_WORLD,rowsTwoSidedIndex,ierr)
        call ISBuildTwoSided(partitionIndex, PETSC_NULL_IS, rowsTwoSidedIndex, ierr)
        call ISGetLocalSize(rowsTwoSidedIndex,localDOFs,ierr)
        print*,rank,"LOCALDOFS",localDOFs

        !!! get Numbering
        call ISCreate(MPI_COMM_WORLD,partitionedNumberingIndex,ierr)
        call ISPartitioningToNumbering(partitionIndex, partitionedNumberingIndex, ierr)

        !!! get Numbering Extended
        call IsGetIndicesF90(partitionedNumberingIndex, parray2, ierr)
        if(rank == 0)then
            allocate(numberingExtendedVal(globalDOFs * 3))
            do i = 1, globalDOFs
                numberingExtendedVal((i-1) * 3 + 1) = parray2(i) * 3 + 0
                numberingExtendedVal((i-1) * 3 + 2) = parray2(i) * 3 + 1
                numberingExtendedVal((i-1) * 3 + 3) = parray2(i) * 3 + 2
            enddo
            call ISCreateGeneral(MPI_COMM_WORLD,globalDOFs*3,numberingExtendedVal,&
                                 PETSC_COPY_VALUES,partitionedNumberingIndex3x, ierr)
        else
            allocate(numberingExtendedVal(1))
            call ISCreateGeneral(MPI_COMM_WORLD,0,           numberingExtendedVal,&
                                 PETSC_COPY_VALUES,partitionedNumberingIndex3x, ierr)
        endif
        deallocate(numberingExtendedVal)
        call IsRestoreIndicesF90(partitionedNumberingIndex, parray2, ierr)

        !!! get Inumbering
        call ISCreate(MPI_COMM_WORLD,partitionedNumberingIndexInversed, ierr)
        call ISInvertPermutation(partitionedNumberingIndex,PETSC_DECIDE, partitionedNumberingIndexInversed,ierr)
        !partitionedNumberingIndexInversed is dataly same as rowsTwoSidedIndex

        !!! call ISSetBlockSize(partitionedNumberingIndex,123, ierr)

        !!!!CHECK
        ! call IsView(partitionedNumberingIndex3x, PETSC_VIEWER_STDOUT_WORLD, ierr)
        ! call MatCreateSubMatrix(adjMat,rowsTwoSidedIndex,rowsTwoSidedIndex,MAT_INITIAL_MATRIX,m2,ierr)
        ! !call MatPermute(adjMat,rowsTwoSidedIndex,rowsTwoSidedIndex,m3, ierr)
        ! !call MatView(m2, PETSC_VIEWER_STDOUT_WORLD, ierr)
        ! call MatGetLocalSize(m2, msizem, msizen, ierr)
        ! print*, rank, "m=", msizem, "n=", msizen
        call IsGetIndicesF90(partitionIndex, parray, ierr)
        if(rank ==0 ) then
            allocate(point_partition(globalDOFs))
            do i = 1,globalDOFs
                point_partition(i) = parray(i)
            enddo
            call output_plt_scalar('./out2_data1_part.plt', 'goodstart', point_partition, 'partition', .false.)
            deallocate(point_partition)
        endif
        call IsRestoreIndicesF90(partitionIndex, parray, ierr)
        !!!!CHECK

        !!! max Row Size
        if (rank == 0) then
            loc = maxloc(AdjacencyNum0)
            print*,"Max", AdjacencyNum0(loc(1))
            maxWidth = AdjacencyNum0(loc(1)) + 1
            loc = minloc(AdjacencyNum0)
            print*,"Min", AdjacencyNum0(loc(1)), size(AdjacencyNum0), loc
        endif
        call MPI_Bcast(maxWidth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adjMaxWidth = maxWidth !!!!!!!!
        !print*,rank, "MW=",maxWidth

        !!! create a mat for testing
        ! call MatCreate(MPI_COMM_WORLD, m3, ierr)
        ! call MatSetType(m3,MATMPIAIJ,ierr)
        ! call MatSetSizes(m3, localDOFs*3, localDOFs*3, globalDOFs*3, globalDOFs*3, ierr)
        ! call MatMPIAIJSetPreallocation(m3,maxWidth*3,PETSC_NULL_INTEGER,maxWidth*3,PETSC_NULL_INTEGER, ierr)
        ! call IsGetIndicesF90(partitionedNumberingIndex, parray, ierr) !parray is now partitionedNumberingIndex
        ! if(rank == 0)then
        !     do i = 1, ncells
        !         elem_id = getElemID(size(CELL(i)%N))
        !         !num_intpoint = elem_lib(elem_id)%num_intpoint
        !         num_node = elem_lib(elem_id)%num_node
        !         allocate(cellDOFs(num_node*3))
        !         allocate(cellMat(num_node*3,num_node*3))
        !         do j = 1, num_node
        !             cellDOFs(3*(j - 1) + 1) = parray(CELL(i)%N(j)-1) * 3 + 0
        !             cellDOFs(3*(j - 1) + 2) = parray(CELL(i)%N(j)-1) * 3 + 1
        !             cellDOFs(3*(j - 1) + 3) = parray(CELL(i)%N(j)-1) * 3 + 2
        !         enddo
        !         cellMat = 1
        !         call MatSetValues(m3, num_node*3, cellDOFs, num_node*3, cellDOFs, cellMat, ADD_VALUES, ierr)
        !         deallocate(cellDOFs)
        !         deallocate(cellMat)
        !     enddo
        ! endif
        ! call MatAssemblyBegin(m3, MAT_FINAL_ASSEMBLY, ierr)
        ! call MatAssemblyEnd(m3, MAT_FINAL_ASSEMBLY, ierr)
        ! !call MatView(m3, PETSC_VIEWER_STDOUT_WORLD, ierr)
        ! call VecCreateMPI(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3, v1, ierr)
        ! call VecCreateMPI(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3, v2, ierr)
        ! call VecSet(v1,1.0_8,ierr)
        ! call MatMult(m3, v1, v2, ierr)
        ! call VecView(v2,PETSC_VIEWER_STDOUT_WORLD, ierr)
        !!
        !!!!!!! distribute adjacency

        if(rank == 0)then
            allocate(tempScalarizedWidth(globalDOFs))
            allocate(tempScalarizedInsertionIndex(globalDOFs))
            do i = 1, globalDOFs
                tempScalarizedInsertionIndex(i) = i-1
            enddo
            tempScalarizedWidth = AdjacencyNum0
            call VecCreateMPI(MPI_COMM_WORLD, globalDOFs, globalDOFs, widthVecZero, ierr)
            call VecSetValues(widthVecZero,globalDOFs,tempScalarizedInsertionIndex,tempScalarizedWidth,INSERT_VALUES,ierr)
        else
            call VecCreateMPI(MPI_COMM_WORLD, 0         , globalDOFs, widthVecZero, ierr)
        endif
        call VecAssemblyBegin(widthVecZero, ierr)
        call VecAssemblyEnd(widthVecZero,ierr)
        !call VecView(widthVecZero,PETSC_VIEWER_STDOUT_WORLD, ierr)

        call VecCreateMPI(MPI_COMM_WORLD, localDOFs, globalDOFs, widthVecDist, ierr)
        call VecGetOwnershipRange(widthVecDist , indexLo, indexHi, ierr)
        call GatherVec1(widthVecDist, widthVecZero, .true. )

        call VecGetArrayReadF90(widthVecDist, pwidth, ierr)
        if(allocated(localAdjacencyNum)) then
            deallocate(localAdjacencyNum)
        endif
        allocate(localAdjacencyNum(localDOFs))
        localAdjacencyNum = nint(pwidth) + 1 ! nearest int, totalnum needs += 1
        call VecRestoreArrayReadF90(widthVecDist, pwidth, ierr)
        call IsGetIndicesF90(partitionIndex, parray, ierr)
        if(rank == 0)then
            do i = 1, globalDOFs
                currentProc = parray(i)
                tempScalarizedWidth(i) = 0._8
                do j = (AdjacencyCounter%rowStart(i)+1), (AdjacencyCounter%rowStart(i+1)-1+1)
                    if(currentProc .ne. parray((AdjacencyCounter%column(j)+1)) )then !if that vert is not on currentProc
                        tempScalarizedWidth(i) = tempScalarizedWidth(i) + 1._8
                    endif
                enddo
            enddo
            call VecSetValues(widthVecZero,globalDOFs,tempScalarizedInsertionIndex,tempScalarizedWidth,INSERT_VALUES,ierr)
            deallocate(tempScalarizedInsertionIndex)
            deallocate(tempScalarizedWidth)
        endif
        call VecAssemblyBegin(widthVecZero, ierr)
        call VecAssemblyEnd(widthVecZero,ierr)
        !call VecView(widthVecZero,PETSC_VIEWER_STDOUT_WORLD, ierr)
        call GatherVec1(widthVecDist, widthVecZero, .true. )
        call VecGetArrayReadF90(widthVecDist, pwidth, ierr)
        if(allocated(ghostAdjacencyNum)) then
            deallocate(ghostAdjacencyNum)
        endif
        allocate(ghostAdjacencyNum(localDOFs))
        print*,size(pwidth), localDOFs
        ghostAdjacencyNum = nint(pwidth) ! nearest int
        localAdjacencyNum = localAdjacencyNum - ghostAdjacencyNum ! save the ghosters
        call VecRestoreArrayReadF90(widthVecDist, pwidth, ierr)
        call VecDestroy(widthVecZero, ierr)
        call VecDestroy(widthVecDist, ierr)
        call IsRestoreIndicesF90(partitionIndex, parray, ierr)
        !!!!!!! distribute adjacency
        call DoCellPartition
        call DoGeomPartition
    end subroutine

    ! do not call
    subroutine DoCellPartition !creates ghostingGlobals and lcoalcells
        use globals
        use linear_set
        !!!!!!! cell partition and cell ghost deriving
        PetscInt ierr
        PetscInt ncells, nverts, i, j, ir, elem_id, num_node, current_row_start
        integer rank, siz
        PetscInt, pointer :: parray(:), parray2(:)
        type(petscInt_vardimWrapper), allocatable :: SendGhosts(:)
        type(Csr), allocatable :: SendCells(:)
        PetscInt ghostLocFound

        integer, allocatable :: sendreqs(:)
        integer recvreq(1)
        !integer, allocatable :: sendstatus(:)
        integer recvstatus(MPI_STATUS_SIZE),sendwaitindex,sendstatus(MPI_STATUS_SIZE)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, siz, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        call IsGetIndicesF90(partitionIndex, parray, ierr)
        if(rank == 0) then
            ncells = size(CELL)
            nverts = size(COORD,dim=2)
            if(allocated(cellPartition)) then
                deallocate(cellPartition)
            endif
            allocate(cellPartition(ncells))
            allocate(localCellSizes(siz), source = 0)
            allocate(SendCells(siz))
            allocate(localGhostSizes(siz), source = 0)
            allocate(SendGhosts(siz))
            do i = 1, ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                cellPartition(i) = parray(CELL(i)%N(mod(i,num_node) + 1))
                localCellSizes(cellPartition(i) + 1) = localCellSizes(cellPartition(i) + 1) + 1
                do j = 1, num_node
                    if(parray(CELL(i)%N(j)) .ne. cellPartition(i))then
                        localGhostSizes(cellPartition(i) + 1) = localGhostSizes(cellPartition(i) + 1) + 1
                    endif
                enddo
            enddo
            do ir = 0 , (siz-1)
                call Csr_SetNROW(SendCells(ir+1), localCellSizes(ir+1))
                call Csr_AllocateRowStart(SendCells(ir+1))
                allocate(SendGhosts(ir+1)%N(localGhostSizes(ir+1)))
            enddo
            localGhostSizes = 0
            localCellSizes = 0
            do i = 1, ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                ir = cellPartition(i)
                localCellSizes(ir + 1) = localCellSizes(ir + 1) + 1
                SendCells(ir + 1)%rowStart(localCellSizes(ir + 1) + 1) = &
                    SendCells(ir + 1)%rowStart(localCellSizes(ir + 1)) + num_node
                do j = 1, num_node
                    if(parray(CELL(i)%N(j)) .ne. cellPartition(i))then
                        localGhostSizes(cellPartition(i) + 1) = localGhostSizes(cellPartition(i) + 1) + 1
                        SendGhosts(cellPartition(i) + 1)%N(localGhostSizes(cellPartition(i) + 1)) = &
                            CELL(i)%N(j)
                    endif
                enddo
            enddo
            do ir = 0 , (siz-1)
                call Csr_SetNNZ(SendCells(ir+1), SendCells(ir+1)%rowStart(SendCells(ir+1)%nrow + 1))
                call Csr_AllocateColumn(SendCells(ir+1))
                call int_mergeSort(SendGhosts(ir+1)%N, 1, localGhostSizes(ir+1) + 1)
                call int_reduceSorted(SendGhosts(ir+1)%N, 1, size(SendGhosts(ir+1)%N)+1,&
                                      localGhostSizes(ir+1))
            enddo
        endif
        call IsGetIndicesF90(partitionedNumberingIndex, parray2, ierr)
        if(rank == 0) then
            localCellSizes = 0
            do i = 1, ncells
                elem_id = getElemID(size(CELL(i)%N))
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                ir = cellPartition(i)
                localCellSizes(ir + 1) = localCellSizes(ir + 1) + 1
                current_row_start = SendCells(ir+1)%rowStart(localCellSizes(ir + 1))
                !next_row_start    = SendCells(ir+1)%rowStart(localCellSizes(ir + 1) + 1)
                do j = 1, num_node
                    if(parray(CELL(i)%N(j)) .ne. cellPartition(i))then
                        if(int_searchBinary(SendGhosts(ir+1)%N, 1, localGhostSizes(ir+1), &
                                            CELL(i)%N(j), ghostLocFound)) then
                            SendCells(ir+1)%column(current_row_start + j) = -ghostLocFound !!! use minus to distinguish from global sets
                        else
                            print*, CELL(i)%N(j),'ghost search failed',SendGhosts(ir+1)%N
                            stop
                        endif
                    else
                        SendCells(ir+1)%column(current_row_start + j) = &
                            parray2(CELL(i)%N(j))
                    endif
                enddo
            enddo
            do ir = 0 , (siz-1)
                !print*,ir,":",SendGhosts(ir+1)%N
                do i = 1, (localGhostSizes(ir+1) - 1) ! localGhostSizes is now 1 larger than top
                    SendGhosts(ir+1)%N(i) = parray2(SendGhosts(ir+1)%N(i))
                enddo
            enddo
            !print*,'rank', rank, 'HERE'
        endif
        ! send nrow
        if(rank == 0) then
            allocate(sendreqs(siz))
            do ir = 0, (siz - 1)
                call MPI_Isend(SendCells(ir+1)%nrow, 1, MPI_INTEGER, ir, 1231, &
                               MPI_COMM_WORLD, sendreqs(ir+1), ierr)
            enddo
        endif
        call MPI_Irecv(localCells%nrow, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1231, &
                       MPI_COMM_WORLD, recvreq(1), ierr)
        if(rank == 0) then
            call MPI_Waitall(siz, sendreqs, MPI_STATUSES_IGNORE, ierr)
        endif
        call MPI_Wait(recvreq(1), recvstatus, ierr)
        call Csr_AllocateRowStart(localCells)
        ! send rowStarts
        if(rank == 0) then
            do ir = 0, (siz - 1)
                call MPI_Isend(SendCells(ir+1)%rowStart, SendCells(ir+1)%nrow+1, MPI_INTEGER, ir, 1232, &
                               MPI_COMM_WORLD, sendreqs(ir+1), ierr)
            enddo
        endif
        call MPI_Irecv(localCells%rowStart, localCells%nrow+1, MPI_INTEGER, MPI_ANY_SOURCE, 1232, &
                       MPI_COMM_WORLD, recvreq(1), ierr)
        if(rank == 0) then
            call MPI_Waitall(siz, sendreqs, MPI_STATUSES_IGNORE, ierr)
        endif
        call MPI_Wait(recvreq(1), recvstatus, ierr)
        ! send nnz
        if(rank == 0) then
            do ir = 0, (siz - 1)
                call MPI_Isend(SendCells(ir+1)%nnz, 1, MPI_INTEGER, ir, 1233, &
                               MPI_COMM_WORLD, sendreqs(ir+1), ierr)
            enddo
        endif
        call MPI_Irecv(localCells%nnz, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1233, &
                       MPI_COMM_WORLD, recvreq(1), ierr)
        if(rank == 0) then
            call MPI_Waitall(siz, sendreqs, MPI_STATUSES_IGNORE, ierr)
        endif
        call MPI_Wait(recvreq(1), recvstatus, ierr)
        call Csr_AllocateColumn(localCells)
        ! send cols
        if(rank == 0) then
            do ir = 0, (siz - 1)
                call MPI_Isend(SendCells(ir+1)%column, SendCells(ir+1)%nnz, MPI_INTEGER, ir, 1234, &
                               MPI_COMM_WORLD, sendreqs(ir+1), ierr)
            enddo
        endif
        call MPI_Irecv(localCells%column, localCells%nnz, MPI_INTEGER, MPI_ANY_SOURCE, 1234, &
                       MPI_COMM_WORLD, recvreq(1), ierr)
        if(rank == 0) then
            call MPI_Waitall(siz, sendreqs, MPI_STATUSES_IGNORE, ierr)
        endif
        call MPI_Wait(recvreq(1), recvstatus, ierr)
        ! send nghosts
        if(rank == 0) then
            do ir = 0, (siz - 1)
                call MPI_Isend(localGhostSizes(ir+1), 1, MPI_INTEGER, ir, 1235, &
                               MPI_COMM_WORLD, sendreqs(ir+1), ierr)
            enddo
        endif
        call MPI_Irecv(nghost, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1235, &
                       MPI_COMM_WORLD, recvreq(1), ierr)
        if(rank == 0) then
            call MPI_Waitall(siz, sendreqs, MPI_STATUSES_IGNORE, ierr)
        endif
        call MPI_Wait(recvreq(1), recvstatus, ierr)
        nghost = nghost - 1 ! mark that the localGhostSizes has never been corrcected after reducing
        print*,'ng',nghost
        ! send ghosts
        allocate(ghostingGlobal(nghost))
        if(rank == 0) then
            do ir = 0, (siz - 1)
                call MPI_Isend(SendGhosts(ir+1)%N, localGhostSizes(ir+1) - 1, MPI_INTEGER, ir, 1236, &
                               MPI_COMM_WORLD, sendreqs(ir+1), ierr)
            enddo
        endif
        call MPI_Irecv(ghostingGlobal, nghost, MPI_INTEGER, MPI_ANY_SOURCE, 1236, &
                       MPI_COMM_WORLD, recvreq(1), ierr)
        if(rank == 0) then
            call MPI_Waitall(siz, sendreqs, MPI_STATUSES_IGNORE, ierr)
        endif
        call MPI_Wait(recvreq(1), recvstatus, ierr)

        allocate(ghostingGlobal3x(nghost*3))
        do i = 1,nghost
            ghostingGlobal3x((i-1)*3+1) = ghostingGlobal(i) * 3 + 0
            ghostingGlobal3x((i-1)*3+2) = ghostingGlobal(i) * 3 + 1
            ghostingGlobal3x((i-1)*3+3) = ghostingGlobal(i) * 3 + 2
        enddo

        !print*,'rank',indexLo,indexHi
        do i = 1, localCells%nnz
            if (localCells%column(i) >= 0 )then
                if(localCells%column(i)>= indexHi .or. localCells%column(i) < indexLo)then
                    print*, 'illegal localCells%column(i) ! '
                    stop
                endif
                localCells%column(i) = localCells%column(i) - indexLo
            else ! is ghost
                localCells%column(i) = -1 - localCells%column(i) + localDOFs
            endif
        enddo

        if(rank == 0) then
            !print*,localCells%column
        endif

        ! do i = 1, localCells%nnz
        !     if (localCells%column(i) < 0 .or. localCells%column(i) >= localDOFs)then
        !         print*,rank,'not local!!'
        !     endif
        !     if (localCells%column(i) < 0 .or. localCells%column(i) >= localDOFs + nghost)then
        !         print*,rank,'not legal!!'
        !     endif
        ! enddo
        nlocalCells = localCells%nrow

        if(rank == 0) then
            do i = 1, siz
                deallocate(SendGhosts(i)%N)
                call Csr_DeleteMat(SendCells(i))
            enddo
            deallocate(SendCells)
            deallocate(SendGhosts)
            deallocate(sendreqs)
        endif

        call IsRestoreIndicesF90(partitionIndex,            parray,  ierr)
        call IsRestoreIndicesF90(partitionedNumberingIndex, parray2, ierr)
    end subroutine

    ! do not call
    subroutine DoGeomPartition ! creates localcoords(needs DoCellPartition) and do first ghost updating
        use globals
        PetscInt ierr, i
        integer rank, siz
        Vec zerocoord

        call MPI_COMM_SIZE(MPI_COMM_WORLD, siz, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        if(rank == 0) then
            call VecCreateMPI(MPI_COMM_WORLD, globalDOFs*3, globalDOFs*3, zerocoord, ierr)
            do i = 1, globalDOFs
                call VecSetValues(zerocoord,3,(/(i-1)*3+0, (i-1)*3+1, (i-1)*3+2/), COORD(:,i),INSERT_VALUES,ierr)
            enddo
        else
            call VecCreateMPI(MPI_COMM_WORLD, 0         ,   globalDOFs*3, zerocoord, ierr)
        endif
        call VecAssemblyBegin(zerocoord,ierr)
        call VecAssemblyEnd  (zerocoord,ierr)
        call VecCreateGhost(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3,nghost*3,ghostingGlobal3x ,localCoords, ierr)
        call GatherVec3(localCoords, zerocoord, .true. )
        !call VecView(localCoords,PETSC_VIEWER_STDOUT_WORLD, ierr)
        call VecDestroy(zerocoord, ierr)
        call VecGhostUpdateBegin(localCoords,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGhostUpdateEnd  (localCoords,INSERT_VALUES,SCATTER_FORWARD, ierr); 
    end subroutine

    ! do not call
    subroutine DoBCPartition ! create local bcvals

    end subroutine

    !after readgrid initializeLib
    subroutine SetUpThermalBC_BLOCKED
        use globals
        use linear_set
        integer :: rank, siz
        real(8) minus1
        PetscInt i, j, ncells, nverts, elem_id, num_node, num_intpoint, ip, in, in2, &
            ie, ibc, ifc, facesize, ierr, maxWidth
        PetscInt faceNodes(8)
        PetscInt bcDofsCount

        !!! startup
        maxWidth = adjMaxWidth
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
        call VecDestroy(dofFixTher, ierr)
        if(allocated(bcDOFsTher))then
            call deallocate_petscInt_vardimWrapper(bcDOFsTher)
            deallocate(bcDOFsTher)
        endif
        if(allocated(bcVALsTher))then
            call deallocate_petscScalar_vardimWrapper(bcVALsTher)
            deallocate(bcVALsTher)
        endif
        if(allocated(bcVALsTher2))then
            call deallocate_petscScalar_vardimWrapper(bcVALsTher2)
            deallocate(bcVALsTher2)
        endif
        !print*,rank,'here'
        if(rank==0) then
            allocate(bcDOFsTher(NBSETS))
            allocate(bcVALsTher(NBSETS))
            allocate(bcVALsTher2(NBSETS))
            ncells = size(CELL)
            nverts = size(COORD,dim=2)
            print*,nverts,'bcnverts'
            !!! mark set dofs
            call VecCreateMPI(MPI_COMM_WORLD, globalDOFs, globalDOFs, dofFixTher, ierr)
            minus1 = -1
            call VecSet(dofFixTher ,sqrt(minus1), ierr)
            do ibc = 1, NBSETS
                if(bcTypeTher(ibc) == 0) then
                    !print*,'fixedbc', ibc, bname(ibc),bcValueTher(ibc)
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        !print*,ibc,bcread(ibc)%ELEM_ID(i),&
                        !    bcread(ibc)%ELEM_FACE(i)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        faceNodes(1:facesize) = elem_lib(elem_id)%face_node(ifc, 1:facesize)
                        do j = 1,facesize
                            faceNodes(j) = CELL(ie)%N(faceNodes(j))
                            call VecSetValue(dofFixTher, faceNodes(j)-1, bcValueTher(ibc),INSERT_VALUES, ierr)
                        enddo
                    enddo
                else if(bcTypeTher(ibc) == 1) then ! linear flow condition
                    bcDofsCount = 0
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        bcDofsCount = bcDofsCount + facesize
                    enddo
                    allocate(bcDOFsTher(ibc)%N(1*bcDofsCount))
                    bcDofsCount = 0
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        bcDofsCount = bcDofsCount + facesize
                        faceNodes(1:facesize) = elem_lib(elem_id)%face_node(ifc, 1:facesize)
                        do j = 1,facesize
                            faceNodes(j) = CELL(ie)%N(faceNodes(j))
                        enddo
                        bcDOFsTher(ibc)%N(bcDofsCount-facesize+1:bcDofsCount) = faceNodes(1:facesize)
                    enddo
                    call int_mergeSort(bcDOFsTher(ibc)%N,1,bcDofsCount+1)
                    call int_reduceSorted(bcDOFsTher(ibc)%N,1,size(bcDOFsTher(ibc)%N),bcDofsCount)
                    allocate(bcVALsTher(ibc) %N(bcDofsCount-1))
                    allocate(bcVALsTher2(ibc)%N(bcDofsCount-1))
                    bcVALsTher(ibc)%N = bcValueTher(ibc)
                    bcVALsTher2(ibc)%N = bcValueTher2(ibc)
                endif
            enddo
        else
            call VecCreateMPI(MPI_COMM_WORLD,      0, globalDOFs, dofFixTher, ierr)
            call VecSet(dofFixTher ,sqrt(minus1), ierr) ! it is collective! must also call
        endif

        call VecAssemblyBegin(dofFixTher,ierr)
        call VecAssemblyEnd  (dofFixTher,ierr)
        call VecCreateGhost(MPI_COMM_WORLD, localDOFs, globalDOFs,nghost,ghostingGlobal ,dofFixTherDist, ierr)
        call GatherVec1(dofFixTherDist, dofFixTher, .true. )
        call VecGhostUpdateBegin(dofFixTherDist,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGhostUpdateEnd  (dofFixTherDist,INSERT_VALUES,SCATTER_FORWARD, ierr); 
    end subroutine

    subroutine SetUpThermal_InitializeObjects ! To create mat and vecs
        PetscInt ierr
        PetscInt maxWidth
        maxWidth = adjMaxWidth
        call MatDestroy(Ather,ierr)
        call VecDestroy(Uther,ierr)
        call VecDestroy(Pther,ierr)
        call MatCreate(MPI_COMM_WORLD, Ather, ierr)
        call MatSetType(Ather,MATMPIAIJ,ierr)
        call MatSetSizes(Ather, localDOFs, localDOFs, globalDOFs, globalDOFs, ierr)
        call MatMPIAIJSetPreallocation(Ather,maxWidth*1,localAdjacencyNum,maxWidth*1,ghostAdjacencyNum, ierr)
        call MatSetOption(Ather,MAT_ROW_ORIENTED,PETSC_FALSE,ierr) !!! Fortran is column major
        call VecCreateMPI(MPI_COMM_WORLD, localDOFs*1, globalDOFs*1, Pther, ierr)
        call VecCreateGhost(MPI_COMM_WORLD, localDOFs*1, globalDOFs*1,nghost, ghostingGlobal, Uther, ierr)
        CHKERRA(ierr)
    end subroutine

    subroutine SetUpThermalPara
        use globals
        use linear_set
        PetscInt ierr
        PetscInt maxWidth
        PetscInt, pointer :: parray(:)
        PetscScalar, pointer :: pfix(:), pfix0(:)
        PetscScalar, pointer :: pcoord(:)
        PetscInt celllo,cellhi
        integer :: rank, siz
        PetscInt i, j, ncells, nverts, elem_id,  face_id, num_node, num_intpoint, ip, in, in2, &
            ie, ibc, ifc, facesize
        PetscInt faceNodes(max_numnode)
        PetscInt :: cellDOFs(max_numnode)
        PetscScalar :: cellVec(max_numnode)
        PetscScalar, allocatable::cellMat(:,:)
        PetscScalar   cellData(max_numnode), cellData2(max_numnode), &
            pointData(max_numintpoint), pointData2(max_numintpoint)
        PetscScalar :: elem_coord(max_numnode,3), Jacobi(3,3), invJacobi(3,3)  !max num node is 27 ! warning
        PetscScalar :: dNjdxi_ij(3,max_numnode), k__mul__dNjdxi_ij(3,max_numnode)
        PetscScalar  nodeFix, nodeDiag, detJacobi, dsscale
        PetscInt boundPos
        logical touchedMat

        !!! startup
        maxWidth = adjMaxWidth
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
        ncells = size(CELL)
        nverts = size(COORD,dim=2)
        call IsGetIndicesF90(partitionedNumberingIndex, parray, ierr) !parray is now partitionedNumberingIndex
        call VecGetArrayReadF90(dofFixTherDist, pfix, ierr)
        call VecGetArrayReadF90(dofFixTher    , pfix0,ierr)
        call VecGetArrayReadF90(localCoords,    pcoord, ierr)

        call VecSet(Pther,0.0_8,ierr)
        !print*,rank,nlocalCells

        ! parallel A assemble
        do i = 1, nLocalCells
            celllo = localCells%rowStart(i)
            cellhi = localCells%rowStart(i+1)
            elem_id = getElemID(cellhi - celllo)
            num_intpoint = elem_lib(elem_id)%num_intpoint
            num_node = elem_lib(elem_id)%num_node
            allocate(cellMat(num_node*1,num_node*1))
            do j = 1, num_node
                if(localCells%column(celllo+j) < localDOFs) then
                    cellDOFs(1*(j - 1) + 1) = (localCells%column(celllo + j) + indexLo) * 1 + 0
                else
                    cellDOFs(1*(j - 1) + 1) = ghostingGlobal(localCells%column(celllo + j) - localDOFs + 1) * 1 + 0
                endif
            enddo
            !!! Integrate the cellMat
            ! get the coords
            do in = 1, num_node
                elem_coord(in,:) = pcoord((localCells%column(celllo+in))*3+1:(localCells%column(celllo+in))*3+3)
            end do
            cellMat = 0.0_8
            do ip = 1, num_intpoint
                Jacobi = matmul(elem_lib(elem_id)%dNIdLimr(:,:,ip),elem_coord(1:num_node,:)) !
                detJacobi = directDet3x3(Jacobi)
                if(detJacobi <= 0.0_8) then
                    print*,'Jacobi minus, mesh distortion'
                    stop
                endif
                invJacobi = directInverse3x3(Jacobi)
                dNjdxi_ij(:,1:num_node) = matmul(invJacobi,elem_lib(elem_id)%dNIdLimr(:,:,ip))
                k__mul__dNjdxi_ij(:,1:num_node) = matmul(k_ther, dNjdxi_ij)
                cellMat = cellMat + matmul(transpose(dNjdxi_ij(:,1:num_node)), k__mul__dNjdxi_ij(:,1:num_node))&
                          * elem_lib(elem_id)%coord_intweight(ip) * detJacobi
            end do

            touchedMat = .false.
            do in = 1, num_node
                if(.not. isnan(pfix(localCells%column(celllo+in)+1))) then
                    nodeFix = pfix(localCells%column(celllo+in)+1)
                    nodeDiag = cellMat(in,in)
                    cellVec(1:num_node) = -cellMat(:,in)
                    do in2 = 1, num_node
                        if(.not. isnan(pfix(localCells%column(celllo+in2)+1))) then
                            cellVec(in2) = 0.0_8 ! do not add rhs for other fixed dofs
                        endif
                    enddo
                    cellVec(in) = nodeDiag
                    cellVec = cellVec * nodeFix
                    call VecSetValues(Pther,num_node,cellDOFs,cellVec,ADD_VALUES,ierr)
                    cellMat(in,:) = 0.0_8
                    cellMat(:,in) = 0.0_8
                    cellMat(in,in) = nodeDiag
                    touchedMat = .true.
                endif
            enddo
            ! if(touchedMat) then
            !     print*, "CELL",i
            !     do in2 = 1, num_node
            !         print*, cellMat(in2,:)
            !     enddo
            ! endif
            ! if(rank == 0) then
            !     print*,'CELLDOFS',cellDOFs
            ! endif
            call MatSetValues(Ather, num_node*1, cellDOFs, num_node*1, cellDOFs, cellMat, ADD_VALUES, ierr)
            !print*,i
            !print*,cellMat
            deallocate(cellMat)
        enddo

        !!!! this bc set is now only serial
        if(rank == 0)then
            do ibc = 1, NBSETS
                if(bcTypeTher(ibc) == 1) then! linear flow condition
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        faceNodes(1:facesize) = elem_lib(elem_id)%face_node(ifc, 1:facesize)
                        face_id = getElemIDFace(facesize)
                        num_intpoint = elem_lib_face(face_id)%num_intpoint
                        allocate(cellMat(facesize*1,facesize*1))
                        do j = 1,facesize
                            faceNodes(j) = CELL(ie)%N(faceNodes(j))
                            elem_coord(j,:) = COORD(:,faceNodes(j))
                            cellDOFs(j) = parray(faceNodes(j))
                            !print*,'DEBUG',bcDOFsTher(ibc)%N
                            if( int_searchBinary(bcDOFsTher(ibc)%N, 1, size(bcVALsTher(ibc)%N)+1, &
                                                 faceNodes(j), boundPos))then
                                cellData(j) = bcVALsTher(ibc)%N(boundPos)
                                cellData2(j) = bcVALsTher2(ibc)%N(boundPos)
                            else
                                print*,'bc search failed ', ibc, i, j
                                stop
                            endif
                        enddo
                        pointData(1:num_intpoint) = matmul(cellData(1:facesize), elem_lib_face(face_id)%NImr)
                        pointData2(1:num_intpoint) = matmul(cellData2(1:facesize), elem_lib_face(face_id)%NImr)
                        cellMat = 0.0_8
                        cellVec = 0.0_8
                        !print*,cellData
                        do ip = 1, num_intpoint
                            Jacobi = matmul(elem_lib_face(face_id)%dNIdLimr(:,:,ip),elem_coord(1:facesize,:))
                            dsscale = norm2(crossProduct(Jacobi(1,:), Jacobi(2,:)))
                            if(dsscale <= 0.0_8) then
                                print*,'surface Jacobi minus, mesh distortion'
                                stop
                            endif
                            cellMat = cellMat + pointData2(ip) * dsscale * elem_lib_face(face_id)%coord_intweight(ip) * &
                                      matmul(&
                                      elem_lib_face(face_id)%NImr(:,ip:ip),&
                                      transpose(elem_lib_face(face_id)%NImr(:,ip:ip)))
                            cellVec(1:facesize) = cellVec(1:facesize) + &
                                                  pointData(ip)  * dsscale * elem_lib_face(face_id)%coord_intweight(ip) * &
                                                  elem_lib_face(face_id)%NImr(:,ip)
                        enddo

                        do j = 1, facesize
                            if(.not. isnan(pfix0(faceNodes(j)))) then
                                cellVec(j) = 0.0_8 ! do not add rhs for fixed dofs
                                cellMat(j,:) = 0.0_8 ! do not add mat for fixed dofs
                                cellMat(:,j) = 0.0_8
                            endif

                        enddo
                        call VecSetValues(Pther, facesize*1, cellDOFs, cellVec,ADD_VALUES,ierr)
                        call MatSetValues(Ather, facesize*1, cellDOFs, facesize*1, cellDOFs, cellMat, ADD_VALUES, ierr)

                        deallocate(cellMat)
                    enddo
                endif
            enddo
        endif

        call MatAssemblyBegin(Ather, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(Ather, MAT_FINAL_ASSEMBLY, ierr)
        call VecAssemblyBegin(Pther, ierr)
        call VecAssemblyEnd(Pther,ierr)
        call VecRestoreArrayReadF90(dofFixTherDist, pfix, ierr)
        call VecRestoreArrayReadF90(localCoords,    pcoord, ierr)
        call VecRestoreArrayReadF90(dofFixTher    , pfix0,ierr)
        call IsRestoreIndicesF90(partitionedNumberingIndex, parray, ierr)
        CHKERRA(ierr)
    end subroutine

    !after readgrid initializeLib
    subroutine SetUpElasticBC_BLOCKED
        use globals
        use linear_set
        integer :: rank, siz
        real(8) minus1
        PetscInt i, j, ncells, nverts, elem_id, num_node, num_intpoint, ip, in, in2, &
            ie, ibc, ifc, facesize, ierr, maxWidth
        PetscInt faceNodes(8)
        PetscInt bcDofsCount

        !!! startup
        maxWidth = adjMaxWidth
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
        call VecDestroy(dofFixElas, ierr)
        if(allocated(bcDOFsElas))then
            call deallocate_petscInt_vardimWrapper(bcDOFsElas)
            deallocate(bcDOFsElas)
        endif
        if(allocated(bcVALsElas))then
            call deallocate_petscScalar_vardimWrapper(bcVALsElas)
            deallocate(bcVALsElas)
        endif
        if(allocated(bcVALsElas2))then
            call deallocate_petscScalar_vardimWrapper(bcVALsElas2)
            deallocate(bcVALsElas2)
        endif
        !print*,rank,'here'
        if(rank==0) then
            allocate(bcDOFsElas(NBSETS))
            allocate(bcVALsElas(NBSETS))
            allocate(bcVALsElas2(NBSETS))
            ncells = size(CELL)
            nverts = size(COORD,dim=2)
            print*,nverts,'bcnverts'
            !!! mark set dofs
            call VecCreateMPI(MPI_COMM_WORLD, globalDOFs*3, globalDOFs*3, dofFixElas, ierr)
            minus1 = -1
            call VecSet(dofFixElas ,sqrt(minus1), ierr)
            do ibc = 1, NBSETS
                if(bcTypeElas(ibc) == 0) then
                    !print*,'fixedbc', ibc, bname(ibc),bcValueElas(ibc)
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        !print*,ibc,bcread(ibc)%ELEM_ID(i),&
                        !    bcread(ibc)%ELEM_FACE(i)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        faceNodes(1:facesize) = elem_lib(elem_id)%face_node(ifc, 1:facesize)
                        do j = 1,facesize
                            faceNodes(j) = CELL(ie)%N(faceNodes(j))
                            call VecSetValues(dofFixElas,3, (/faceNodes(j)*3-3,faceNodes(j)*3-2,faceNodes(j)*3-1/),&
                                              bcValueElas(ibc*3-2:ibc*3),INSERT_VALUES, ierr)
                        enddo
                    enddo
                else if(bcTypeElas(ibc) == 1) then ! linear flow condition
                    bcDofsCount = 0
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        bcDofsCount = bcDofsCount + facesize
                    enddo
                    allocate(bcDOFsElas(ibc)%N(bcDofsCount))
                    bcDofsCount = 0
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        bcDofsCount = bcDofsCount + facesize
                        faceNodes(1:facesize) = elem_lib(elem_id)%face_node(ifc, 1:facesize)
                        do j = 1,facesize
                            faceNodes(j) = CELL(ie)%N(faceNodes(j))
                            bcDOFsElas(ibc)%N(bcDofsCount-facesize+j) = faceNodes(j)
                        enddo
                    enddo
                    call int_mergeSort(bcDOFsElas(ibc)%N,1,bcDofsCount+1)
                    call int_reduceSorted(bcDOFsElas(ibc)%N,1,size(bcDOFsElas(ibc)%N),bcDofsCount)
                    allocate(bcVALsElas(ibc) %N(3*bcDofsCount-1))
                    allocate(bcVALsElas2(ibc)%N(9*(bcDofsCount-1))) !note that bcVALsElas2 is 9-sized blocking
                    do i = 1, (bcDofsCount-1)
                        bcVALsElas(ibc)%N(i*3-2:i*3) = bcValueElas(ibc*3-2:ibc*3)
                        bcVALsElas2(ibc)%N(i*9-8:i*9) = bcValueElas2(ibc*9-8:ibc*9)
                    enddo
                endif
            enddo
        else
            call VecCreateMPI(MPI_COMM_WORLD,      0, globalDOFs*3, dofFixElas, ierr)
            call VecSet(dofFixElas ,sqrt(minus1), ierr) ! it is collective! must also call
        endif
        call VecAssemblyBegin(dofFixElas,ierr)
        call VecAssemblyEnd  (dofFixElas,ierr)
        call VecCreateGhost(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3,nghost*3,ghostingGlobal3x ,dofFixElasDist, ierr)
        call GatherVec3(dofFixElasDist, dofFixElas, .true. )
        call VecGhostUpdateBegin(dofFixElasDist,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGhostUpdateEnd  (dofFixElasDist,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        !call VecView(dofFixElasDist,PETSC_VIEWER_STDOUT_WORLD, ierr)
    end subroutine

    subroutine SetUpElasticity_InitializeObjects ! To create mat and vecs
        PetscInt ierr
        PetscInt maxWidth,i
        PetscInt,allocatable :: localAdjacencyNum3x(:),ghostAdjacencyNum3x(:)
        maxWidth = adjMaxWidth
        allocate(localAdjacencyNum3x(3*localDOFs))
        allocate(ghostAdjacencyNum3x(3*localDOFs))
        do i = 1, localDOFs
            localAdjacencyNum3x(i*3-2:i*3) = localAdjacencyNum(i) * 3
            ghostAdjacencyNum3x(i*3-2:i*3) = ghostAdjacencyNum(i) * 3
        enddo
        call MatDestroy(Aelas,ierr)
        call VecDestroy(Uelas,ierr)
        call VecDestroy(Pelas,ierr)
        call MatCreate(MPI_COMM_WORLD, Aelas, ierr)
        call MatSetType(Aelas,MATMPIAIJ,ierr)
        call MatSetSizes(Aelas, localDOFs*3, localDOFs*3, globalDOFs*3, globalDOFs*3, ierr)
        call MatMPIAIJSetPreallocation(Aelas,maxWidth*3,localAdjacencyNum3x,maxWidth*3,ghostAdjacencyNum3x, ierr)
        call MatSetOption(Aelas,MAT_ROW_ORIENTED,PETSC_FALSE,ierr) !!! Fortran is column major
        call VecCreateMPI(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3, Pelas, ierr)
        call VecCreateGhost(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3,nghost*3,ghostingGlobal3x ,Uelas, ierr)
        deallocate(localAdjacencyNum3x)
        deallocate(ghostAdjacencyNum3x)
        CHKERRA(ierr)
    end subroutine

    subroutine SetUpElasticityPara
        use globals
        use linear_set
        use elastic_constitution
        PetscInt ierr
        PetscInt maxWidth
        PetscInt, pointer :: parray(:)
        PetscScalar, pointer :: pfix(:), pfix0(:)
        PetscScalar, pointer :: pcoord(:), pphi(:)
        PetscInt celllo,cellhi
        integer :: rank, siz
        PetscInt i, j, ncells, nverts, elem_id,  face_id, num_node, num_intpoint, ip, in, in2, &
            ie, ibc, ifc, facesize, idim
        PetscInt faceNodes(8)
        PetscInt cellDOFs(max_numnode*3)
        PetscScalar  cellVec(max_numnode*3),localFix(max_numnode*3)
        PetscScalar, allocatable::cellMat(:,:)
        PetscScalar cellData(max_numnode*3), cellData2(max_numnode*9), &
            pointData(max_numintpoint*3), pointData2(max_numintpoint*9), &
            cellDataTemp(max_numnode), pointDataTemp(max_numintpoint)
        PetscScalar :: elem_coord(max_numnode,3), Jacobi(3,3), invJacobi(3,3)  !max num node is 27 ! warning
        PetscScalar :: dNjdxi_ij(3,max_numnode), b__mul__dNjdxi_ij(6,max_numnode*3), &
            dNjdxi_ij_expand(9,max_numnode*3), NImr_expand(3,max_numnode*3), &
            b_dir2strain(6,9), unitbulkstrain(6), btd(max_numnode*3,6)
        ! the dNjdxi_ij_expand is [mat],   [mat]*[u1v1w1u2v2w2...]' = [dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz]'
        ! the b is geometry relation , b*dir = strain, dir = [dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz]'
        !                              strain = [sxx syy szz syz szx sxy]'
        PetscScalar  nodeFix, nodeDiag, detJacobi, dsscale, E, nu
        PetscInt boundPos
        logical touchedMat
        real(8) evreal(max_numnode*3), evimag(max_numnode*3), work(12*max_numnode)

        call VecSet(Pelas,0.0_8,ierr)

        !!! startup
        maxWidth = adjMaxWidth
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
        ncells = size(CELL)
        nverts = size(COORD,dim=2)
        call IsGetIndicesF90(partitionedNumberingIndex, parray, ierr) !parray is now partitionedNumberingIndex
        call VecGetArrayReadF90(dofFixElasDist, pfix, ierr)
        call VecGetArrayReadF90(dofFixElas    , pfix0,ierr)
        call VecGetArrayReadF90(localCoords,    pcoord, ierr)
        call VecGhostUpdateBegin(Uther,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGhostUpdateEnd  (Uther,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGetArrayReadF90(Uther,   pphi, ierr)
        !   print*, pphi

        dNjdxi_ij_expand = 0.0_8
        b_dir2strain = geometryRelation()
        ! set geometry relation, should dump to constitution
        unitbulkstrain = getUnitBulkStrain()

        !print*,rank,nlocalCells

        ! parallel A assemble
        do i = 1, nLocalCells
            celllo = localCells%rowStart(i)
            cellhi = localCells%rowStart(i+1)
            elem_id = getElemID(cellhi - celllo)
            num_intpoint = elem_lib(elem_id)%num_intpoint
            num_node = elem_lib(elem_id)%num_node
            allocate(cellMat(num_node*3,num_node*3))
            do j = 1, num_node
                if(localCells%column(celllo+j) < localDOFs) then
                    cellDOFs(3*(j-1)+1) = (localCells%column(celllo + j) + indexLo) * 3 + 0
                    cellDOFs(3*(j-1)+2) = (localCells%column(celllo + j) + indexLo) * 3 + 1
                    cellDOFs(3*(j-1)+3) = (localCells%column(celllo + j) + indexLo) * 3 + 2
                else
                    cellDOFs(3*(j-1)+1) = ghostingGlobal(localCells%column(celllo + j)-localDOFs+1) * 3 + 0
                    cellDOFs(3*(j-1)+2) = ghostingGlobal(localCells%column(celllo + j)-localDOFs+1) * 3 + 1
                    cellDOFs(3*(j-1)+3) = ghostingGlobal(localCells%column(celllo + j)-localDOFs+1) * 3 + 2
                endif
                !local dofs are x1 y1 z1 x2 y2 z2 ...
                cellData(j) = pphi(localCells%column(celllo+j) + 1) ! this is temprature phi
                localfix(j*3-2:j*3) = pfix(localCells%column(celllo+j)*3+1 : localCells%column(celllo+j)*3+3)
            enddo
            !print*,localFix
            pointData(1:num_intpoint) = matmul(cellData(1:num_node), elem_lib(elem_id)%NImr)
            !!! Integrate the cellMat
            ! get the coords
            do in = 1, num_node
                elem_coord(in,:) = pcoord((localCells%column(celllo+in))*3+1:(localCells%column(celllo+in))*3+3)
            end do
            cellMat = 0.0_8
            cellVec = 0.0_8
            do ip = 1, num_intpoint
                Jacobi = matmul(elem_lib(elem_id)%dNIdLimr(:,:,ip),elem_coord(1:num_node,:)) !
                detJacobi = directDet3x3(Jacobi)
                if(detJacobi <= 0.0_8) then
                    print*,'Jacobi minus, mesh distortion'
                    stop
                endif
                invJacobi = directInverse3x3(Jacobi)
                dNjdxi_ij(:,1:num_node) = matmul(invJacobi,elem_lib(elem_id)%dNIdLimr(:,:,ip))
                do in = 1, num_node
                    dNjdxi_ij_expand(1:3,in*3-2) = dNjdxi_ij(:,in)
                    dNjdxi_ij_expand(4:6,in*3-1) = dNjdxi_ij(:,in)
                    dNjdxi_ij_expand(7:9,in*3-0) = dNjdxi_ij(:,in)
                enddo
                b__mul__dNjdxi_ij(:,1:num_node*3) = matmul(b_dir2strain, dNjdxi_ij_expand(:,1:num_node*3))
                call tempToENU_0(pointData(ip), E, nu) ! use temp-realted E and nu
                btd(1:num_node*3,:) = matmul(transpose(b__mul__dNjdxi_ij(:,1:num_node*3)), constitutionMat(E,nu))

                cellMat = cellMat + matmul(btd (1:num_node*3,:) , b__mul__dNjdxi_ij(:,1:num_node*3)) &
                          * elem_lib(elem_id)%coord_intweight(ip) * detJacobi
                cellVec(1:num_node*3) = cellVec(1:num_node*3) + &
                                        matmul(btd(1:num_node*3,:), expandRate(pointData(ip))*unitbulkstrain)&
                                        * elem_lib(elem_id)%coord_intweight(ip) * detJacobi
            end do
            !print*,pointData

            do in2 = 1, num_node*3
                if(.not. isnan(localfix(in2))) then
                    cellVec(in2) = 0.0_8 ! do not add rhs for fixed dofs
                endif
            enddo

            call VecSetValues(PElas,num_node*3,cellDOFs,cellVec,ADD_VALUES,ierr)
            touchedMat = .false.
            do in = 1, num_node*3
                if(.not. isnan(localfix(in))) then
                    nodeFix = localfix(in)
                    !print*,nodeFix,mod(in,3)
                    nodeDiag = cellMat(in,in)
                    cellVec(1:num_node*3) = -cellMat(:,in)
                    do in2 = 1, num_node*3
                        if(.not. isnan(localfix(in2))) then
                            cellVec(in2) = 0.0_8 ! do not add rhs for fixed dofs
                        endif
                    enddo
                    cellVec(in) = nodeDiag
                    cellVec = cellVec * nodeFix
                    call VecSetValues(PElas,num_node*3,cellDOFs,cellVec,ADD_VALUES,ierr)
                    cellMat(in,:) = 0.0_8
                    cellMat(:,in) = 0.0_8
                    cellMat(in,in) = nodeDiag
                    touchedMat = .true.
                endif
            enddo
            ! if(touchedMat) then
            !     print*, "CELL",i
            !     do in2 = 1, num_node
            !         print*, cellMat(in2,:)
            !     enddo
            ! endif
            ! if(rank == 0) then
            !     print*,'CELLDOFS',cellDOFs
            ! endif
            ! print*,'Sym check',norm2(cellMat-transpose(cellMat))
            ! call dgeev        ('N',&
            !                    'N',&
            !                    num_node*3,&
            !                    cellMat,&
            !                    num_node*3,&
            !                    evreal,&
            !                    evimag,&
            !                    evreal,&
            !                    num_node*3,&
            !                    evreal,&
            !                    num_node*3,&
            !                    work,&
            !                    12*num_node,&
            !                    ierr &
            !                    )
            ! print*,'Eig check',evimag(1:num_node*3)
            call MatSetValues(AElas, num_node*3, cellDOFs, num_node*3, cellDOFs, cellMat, ADD_VALUES, ierr)
            !print*,i
            !print*,cellMat
            deallocate(cellMat)
        enddo

        !!!! this bc set is now only serial
        NImr_expand = 0.0_8
        if(rank == 0)then
            do ibc = 1, NBSETS
                if(bcTypeElas(ibc) == 1) then! linear flow condition
                    do i = 1, size(bcread(ibc)%ELEM_ID)
                        ie = bcread(ibc)%ELEM_ID(i)
                        ifc = bcread(ibc)%ELEM_FACE(i)
                        elem_id = getElemID(size(CELL(ie)%N))
                        facesize = elem_lib(elem_id)%face_size(ifc)
                        faceNodes(1:facesize) = elem_lib(elem_id)%face_node(ifc, 1:facesize)
                        face_id = getElemIDFace(facesize)
                        num_intpoint = elem_lib_face(face_id)%num_intpoint
                        allocate(cellMat(facesize*3,facesize*3))
                        do j = 1,facesize
                            faceNodes(j) = CELL(ie)%N(faceNodes(j))
                            elem_coord(j,:) = COORD(:,faceNodes(j))
                            cellDOFs(j*3-2) = parray(faceNodes(j))*3+0
                            cellDOFs(j*3-1) = parray(faceNodes(j))*3+1
                            cellDOFs(j*3-0) = parray(faceNodes(j))*3+2

                            !print*,'DEBUG',bcDOFsElas(ibc)%N
                            if( int_searchBinary(bcDOFsElas(ibc)%N, 1, size(bcVALsElas(ibc)%N)/3+1, &
                                                 faceNodes(j), boundPos))then
                                cellData(j*3-2:j*3) = bcVALsElas(ibc)%N(boundPos*3-2:boundPos*3)
                                cellData2(j*9-8:j*9) = bcVALsElas2(ibc)%N(boundPos*9-8:boundPos*9)
                            else
                                print*,'bc search failed ', ibc, i, j
                                stop
                            endif
                            localFix(j*3-2:j*3) = pfix(faceNodes(j)*3-2:faceNodes(j)*3)
                        enddo
                        do idim = 1,3
                            do j = 1,facesize
                                cellDataTemp(j) = cellData(j*3-3+idim)
                            enddo
                            pointDataTemp(1:num_intpoint) = matmul(cellDataTemp(1:facesize), elem_lib_face(face_id)%NImr)
                            do j = 1,num_intpoint
                                pointData(j*3-3+idim) = pointDataTemp(j)
                            enddo
                        enddo
                        do idim = 1,9
                            do j = 1,facesize
                                cellDataTemp(j) = cellData2(j*9-9+idim)
                            enddo
                            pointDataTemp(1:num_intpoint) = matmul(cellDataTemp(1:facesize), elem_lib_face(face_id)%NImr)
                            do j = 1,num_intpoint
                                pointData2(j*9-9+idim) = pointDataTemp(j)
                            enddo
                        enddo
                        cellMat = 0.0_8
                        cellVec = 0.0_8
                        !print*,cellData
                        do ip = 1, num_intpoint
                            Jacobi = matmul(elem_lib_face(face_id)%dNIdLimr(:,:,ip),elem_coord(1:facesize,:))
                            dsscale = norm2(crossProduct(Jacobi(1,:), Jacobi(2,:)))
                            if(dsscale <= 0.0_8) then
                                print*,'surface Jacobi minus, mesh distortion'
                                stop
                            endif
                            do j = 1,facesize
                                NImr_expand(1, j*3-2) = elem_lib_face(face_id)%NImr(j,ip)
                                NImr_expand(2, j*3-1) = elem_lib_face(face_id)%NImr(j,ip)
                                NImr_expand(3, j*3-0) = elem_lib_face(face_id)%NImr(j,ip)
                            enddo
                            cellMat = cellMat + dsscale * elem_lib_face(face_id)%coord_intweight(ip) * &
                                      matmul(&
                                      matmul(transpose(NImr_expand(:,1:3*facesize)),reshape(pointData2(ip*9-8:ip*9),(/3,3/))),&
                                      NImr_expand(:,1:3*facesize))
                            cellVec(1:facesize*3) = cellVec(1:facesize*3) + dsscale * &
                                                    elem_lib_face(face_id)%coord_intweight(ip)&
                                                    * &
                                                    matmul(pointData(ip*3-2:ip*3), NImr_expand(:,1:3*facesize))
                            !print*,pointData(ip*3-2:ip*3)
                        enddo

                        do j = 1, facesize*3
                            if(.not. isnan(localFix(j))) then
                                cellVec(j) = 0.0_8 ! do not add rhs for fixed dofs
                                cellMat(j,:) = 0.0_8 ! do not add mat for fixed dofs
                                cellMat(:,j) = 0.0_8
                            endif
                        enddo
                        call VecSetValues(Pelas, facesize*3, cellDOFs, cellVec,ADD_VALUES,ierr)
                        call MatSetValues(Aelas, facesize*3, cellDOFs, facesize*3, cellDOFs, cellMat, ADD_VALUES, ierr)
                        deallocate(cellMat)
                    enddo
                endif
            enddo
        endif

        call MatAssemblyBegin(AElas, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(AElas, MAT_FINAL_ASSEMBLY, ierr)
        call VecAssemblyBegin(PElas, ierr)
        call VecAssemblyEnd(PElas,ierr)
        call VecRestoreArrayReadF90(dofFixElasDist, pfix, ierr)
        call VecRestoreArrayReadF90(localCoords,    pcoord, ierr)
        call VecRestoreArrayReadF90(dofFixElas    , pfix0,ierr)
        call VecRestoreArrayReadF90(Uther,   pphi, ierr)
        call IsRestoreIndicesF90(partitionedNumberingIndex, parray, ierr)
        CHKERRA(ierr)
    end subroutine

    subroutine ClearThermal

    end subroutine

    subroutine ClearElasticity

    end subroutine

    subroutine SolveThermal_Initialize !after setting up
        PetscInt ierr
        PetscBool flag
        !call MatView(Ather,PETSC_VIEWER_STDOUT_WORLD, ierr)
        call MatIsSymmetric(Ather,1e-10_8,flag, ierr)
        print*,"Thermal Symmetric Check: ",flag
        call KSPCreate(MPI_COMM_WORLD, KSPther, ierr)
        call KSPSetOperators(KSPther, Ather, Ather, ierr)
        call KSPSetFromOptions(KSPther, ierr)
    end subroutine

    subroutine SolveThermal
        PetscInt ierr
        call KSPSolve(KSPther, Pther, Uther, ierr)
        !call VecView(Pther, PETSC_VIEWER_STDOUT_WORLD, ierr)
    end subroutine

    subroutine SolveElasticity_Initialize !after setting up
        PetscInt ierr
        PetscBool flag
        call MatIsSymmetric(Aelas,1e-10_8,flag, ierr)
        print*,"Elastic Symmetric Check: ",flag
        call KSPCreate(MPI_COMM_WORLD, KSPelas, ierr)
        call KSPSetOperators(KSPelas, Aelas, Aelas, ierr)
        call KSPSetFromOptions(KSPelas, ierr)
    end subroutine

    subroutine SolveElasticity
        PetscInt ierr
        call KSPSolve(KSPelas, Pelas, Uelas, ierr)
        !call VecView(Uelas, PETSC_VIEWER_STDOUT_WORLD, ierr)
    end subroutine

    subroutine GetElasticityUGradient ! uses lapack
        PetscInt ierr,i,j,idim,in,ip,celllo, cellhi, elem_id, num_intpoint, num_node
        PetscInt cellDOFs(max_numnode)
        PetscScalar cellData(max_numnode,3), elem_coord(max_numnode,3), &
            pointGradient(3,3), Jacobi(3,3), invJacobi(3,3),&
            pointGradientVecs(9,max(max_numintpoint, max_numnode)), NIMRT(max_numintpoint,max_numnode),&
            pointGradientSols(max(max_numintpoint, max_numnode),9),&
            nodeOnes(max_numnode),&
            nodeSolvedSingularValues(max(max_numintpoint, max_numnode))
        !PetscScalar GradientAverage(9)
        PetscScalar,allocatable :: workarray(:)
        PetscScalar,pointer :: pU(:), pCoord(:)
        integer rankout, lapackinfo, workarraysize

        Vec v2cadjcount

        !don't forget the ghost queries don't come along automatically
        call VecGhostUpdateBegin(Uelas,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGhostUpdateEnd  (Uelas,INSERT_VALUES,SCATTER_FORWARD, ierr); 
        call VecGetArrayReadF90(Uelas,pU,ierr)
        call VecGetArrayReadF90(localCoords,pCoord,ierr)
        call VecCreateMPI(MPI_COMM_WORLD, localDOFs, globalDOFs,v2cadjcount, ierr)
        call VecSet(v2cadjcount, 0.0_8, ierr)

        do idim = 1, 9
            call VecCreateMPI(MPI_COMM_WORLD, localDOFs, globalDOFs ,UGradientElas(idim), ierr)
            call VecSet(v2cadjcount, 0.0_8, ierr)
        enddo

        allocate(workarray(4096))
        nodeOnes = 1
        do i = 1, nLocalCells
            celllo = localCells%rowStart(i)
            cellhi = localCells%rowStart(i+1)
            elem_id = getElemID(cellhi - celllo)
            num_intpoint = elem_lib(elem_id)%num_intpoint
            num_node = elem_lib(elem_id)%num_node
            do j = 1, num_node
                if(localCells%column(celllo+j) < localDOFs) then
                    cellDOFs(1*(j-1)+1) = (localCells%column(celllo + j) + indexLo) * 1 + 0
                else
                    cellDOFs(1*(j-1)+1) = ghostingGlobal(localCells%column(celllo + j)-localDOFs+1) * 1 + 0
                endif
                !local dofs are p1 p2 p3 ...
                cellData(j,:) = pU(localCells%column(celllo+j) * 3+1:localCells%column(celllo+j) * 3+3)
                elem_coord(j,:) = pcoord((localCells%column(celllo+j))*3+1:(localCells%column(celllo+j))*3+3)
                !cell data is now disps
            enddo
            !print*,cellData

            do ip = 1, num_intpoint
                Jacobi = matmul(elem_lib(elem_id)%dNIdLimr(:,:,ip),elem_coord(1:num_node,:))
                invJacobi = directInverse3x3(Jacobi)
                pointGradient = matmul(invJacobi,matmul(elem_lib(elem_id)%dNIdLimr(:,:,ip),cellData(1:num_node,:)))
                !pointGradient is duj/dxi
                pointGradientVecs(:,ip) = reshape((pointGradient),(/9/))

            enddo
            pointGradientSols = transpose(pointGradientVecs)
            !!!!!!!!!!!!lapack
            call dgelss(num_intpoint, num_node, 9, NIMRT, max_numintpoint,&
                        pointGradientSols, size(pointGradientSols,1), nodeSolvedSingularValues,&
                        -1.0_8, rankout, workarray,-1, lapackinfo)
            if(lapackinfo .ne. 0)then
                print*, 'GetElasticityUGradient::lapack dgelss QueryFailed'
                stop
            endif
            workarraysize = nint(workarray(1))
            if(workarraysize>size(workarray)) then
                deallocate(workarray)
                allocate(workarray(workarraysize))
            endif
            NIMRT(1:num_intpoint, 1:num_node) = transpose(elem_lib(elem_id)%NImr)
            call dgelss(num_intpoint, num_node, 9, NIMRT, max_numintpoint,&
                        pointGradientSols, size(pointGradientSols,1), nodeSolvedSingularValues,&
                        -1.0_8, rankout, workarray,workarraysize, lapackinfo)
            if(lapackinfo .ne. 0)then
                print*, 'GetElasticityUGradient::lapack dgelss Failed'
                stop
            endif

            call VecSetValues(v2cadjcount, num_node, cellDOFs, nodeOnes, ADD_VALUES,ierr)
            do idim = 1,9
                !GradientAverage(idim) = sum(pointGradientVecs(idim,:))/num_intpoint
                call VecSetValues(UGradientElas(idim), num_node, cellDOFs, pointGradientSols(:,idim), ADD_VALUES,ierr)
            enddo
            !print*,GradientAverage

        enddo
        call VecAssemblyBegin(v2cadjcount, ierr)
        do idim = 1,9
            call VecAssemblyBegin(UGradientElas(idim),ierr)
        enddo
        call VecAssemblyEnd(v2cadjcount, ierr)
        do idim = 1,9
            call VecAssemblyEnd(UGradientElas(idim),ierr)
        enddo
        do idim = 1,9
            call VecPointWiseDivide(UGradientElas(idim),UGradientElas(idim),v2cadjcount,ierr)
        enddo

        deallocate(workarray)
        call VecRestoreArrayReadF90(Uelas,pU,ierr)
        call VecRestoreArrayReadF90(localCoords,pCoord,ierr)
        !call VecView(v2cadjcount,PETSC_VIEWER_STDOUT_WORLD, ierr)
        call VecDestroy(v2cadjcount, ierr)

    end subroutine

    subroutine GetStrainStress
        use elastic_constitution
        PetscInt ierr, idim, i
        PetscScalar b_dir2strain(6,9), strain(6), stress(6), E, nu, grad(9), ubulkStrain(6)
        type(petscScalar_pointerWrapper) :: pVM, pGrad(9), pStrain(6), pStress(6), pphi

        b_dir2strain = geometryRelation()
        ubulkStrain = getUnitBulkStrain()
        if(.not. if_VMElas_alive) then
            call VecCreateMPI(MPI_COMM_WORLD, localDOFs, globalDOFs ,VonMisesElas, ierr)
        endif
        if_VMElas_alive = .true.
        if(.not. if_StrainElas_alive)then
            do idim = 1,6
                call VecCreateMPI(MPI_COMM_WORLD, localDOFs, globalDOFs ,StrainElas(idim), ierr)
            enddo
        endif
        if_StrainElas_alive = .true.
        if(.not. if_StressElas_alive) then
            do idim = 1,6
                call VecCreateMPI(MPI_COMM_WORLD, localDOFs, globalDOFs ,StressElas(idim), ierr)
            enddo
        endif
        if_StressElas_alive = .true.

        call VecGetArrayF90(VonMisesElas, pVM%P, ierr)
        call VecGetArrayReadF90(Pther, pphi%P, ierr)
        do idim = 1,9
            call VecGetArrayReadF90(UGradientElas(idim), pGrad(idim)%P, ierr)
        enddo
        do idim = 1,6
            call VecGetArrayF90(StrainElas(idim), pStrain(idim)%P, ierr)
            call VecGetArrayF90(StressElas(idim), pStress(idim)%P, ierr)
        enddo
        do i = 1, localDOFs
            do idim = 1,9
                grad(idim) = pGrad(idim)%P(i)
            enddo
            strain = matmul(b_dir2strain, grad)
            do idim = 1,6
                pStrain(idim)%P(i) = strain(idim)
            enddo
            call tempToENU_0(pphi%P(i), E, nu)
            stress = matmul(constitutionMat(E,nu),strain-ubulkStrain*expandRate(pphi%P(i)))
            do idim = 1,6
                pStress(idim)%P(i) = stress(idim)
            enddo
            pVM%P(i) = getVonMises(stress)
        enddo
        call VecRestoreArrayF90(VonMisesElas, pVM%P, ierr)
        call VecRestoreArrayReadF90(Pther, pphi%P, ierr)
        do idim = 1,9
            call VecRestoreArrayReadF90(UGradientElas(idim), pGrad(idim)%P, ierr)
        enddo
        do idim = 1,6
            call VecRestoreArrayF90(StrainElas(idim), pStrain(idim)%P, ierr)
            call VecRestoreArrayF90(StressElas(idim), pStress(idim)%P, ierr)
        enddo
        call output_plt_vert_vec1('./out2_debug.plt','goodstartWhat', VonMisesElas, 'VMdebug')
    end subroutine

    ! when reverse, zero->dist, else dist->zero
    subroutine GatherVec1(vecdist, veczero, ifreverse)
        Vec vecdist
        Vec veczero
        ! Vec vecdist2
        PetscInt ierr
        logical ifreverse
        if(.not. (if_procToZeroScatter_alive)) then
            call VecScatterCreate(vecdist, partitionedNumberingIndex, veczero,PETSC_NULL_IS,&
                                  procToZeroScatter, ierr)
            if_procToZeroScatter_alive = .true.
        endif
        if(ifreverse) then
            call VecScatterBegin(procToZeroScatter,veczero,vecdist,INSERT_VALUES,SCATTER_REVERSE,ierr)
            call VecScatterEnd(procToZeroScatter,veczero,vecdist,INSERT_VALUES,SCATTER_REVERSE,ierr)
        else
            call VecScatterBegin(procToZeroScatter,vecdist,veczero,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(procToZeroScatter,vecdist,veczero,INSERT_VALUES,SCATTER_FORWARD,ierr)
        endif
        ! ! cheking if reverse is good
        ! call VecDuplicate(vecdist, vecdist2, ierr)
        ! call VecScatterBegin(procToZeroScatter,veczero,vecdist2,INSERT_VALUES,SCATTER_REVERSE,ierr)
        ! call VecScatterEnd(procToZeroScatter,veczero,vecdist2,INSERT_VALUES,SCATTER_REVERSE,ierr)
        ! call VecAXPY(vecdist2,-1.0_8,vecdist,ierr)
        ! call VecView(vecdist2, PETSC_VIEWER_STDOUT_WORLD, ierr)
    end subroutine

    ! when reverse, zero->dist, else dist->zero
    subroutine GatherVec3(vecdist, veczero, ifreverse)
        Vec vecdist
        Vec veczero
        PetscInt ierr
        logical ifreverse
        if(.not. (if_procToZeroScatter3x_alive)) then
            call VecScatterCreate(vecdist, partitionedNumberingIndex3x, veczero,PETSC_NULL_IS,&
                                  procToZeroScatter3x, ierr)
            if_procToZeroScatter3x_alive = .true.
        endif
        if(ifreverse) then
            call VecScatterBegin(procToZeroScatter3x,veczero,vecdist,INSERT_VALUES,SCATTER_REVERSE,ierr)
            call VecScatterEnd(procToZeroScatter3x,veczero,vecdist,INSERT_VALUES,SCATTER_REVERSE,ierr)
        else
            call VecScatterBegin(procToZeroScatter3x,vecdist,veczero,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(procToZeroScatter3x,vecdist,veczero,INSERT_VALUES,SCATTER_FORWARD,ierr)
        endif
    end subroutine

    subroutine output_plt_thermal(path, title)
        character(*), intent(in) :: path, title
        call output_plt_vert_vec1(path,title,Uther,'phi')
    end subroutine

    subroutine output_plt_elasticity(path, title)
        character(*), intent(in) :: path, title
        call output_plt_vert_vec3(path,title,Uelas,'u','v','w')
    end subroutine

    subroutine output_plt_elasticity_all(path, title)
        PetscInt ierr, rank, idim
        character(*), intent(in) :: path, title
        type(output_plt_value) outArray(1+6+6+2)
        Vec gathered(1+6+6+2)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        if(rank == 0) then
            call VecCreateMPI(MPI_COMM_WORLD,globalDOFs*3,PETSC_DETERMINE, gathered(1), ierr)
            do idim = 2,(1+6+6+2)
                call VecCreateMPI(MPI_COMM_WORLD,globalDOFs,PETSC_DETERMINE, gathered(idim), ierr)
            enddo
        else
            call VecCreateMPI(MPI_COMM_WORLD,           0,PETSC_DETERMINE, gathered(1), ierr)
            do idim = 2,(1+6+6+2)
                call VecCreateMPI(MPI_COMM_WORLD,         0,PETSC_DETERMINE, gathered(idim), ierr)
            enddo
        endif
        call GatherVec3(Uelas, gathered(1), .false.)
        do idim = 1, 6
            call GatherVec1(StrainElas(idim), gathered(1+idim), .false.)
        enddo
        do idim = 1, 6
            call GatherVec1(StressElas(idim), gathered(1+6+idim), .false.)
        enddo
        call GatherVec1(VonMisesElas, gathered(1+6+6+1), .false.)
        call GatherVec1(Uther, gathered(1+6+6+2), .false.)
        outArray(1)%varName="u"
        outArray(1)%varName2="v"
        outArray(1)%varName3="w"
        outArray(1+1)%varName="e11"
        outArray(1+2)%varName="e22"
        outArray(1+3)%varName="e33"
        outArray(1+4)%varName="e23"
        outArray(1+5)%varName="e31"
        outArray(1+6)%varName="e12"
        outArray(1+6+1)%varName="s11"
        outArray(1+6+2)%varName="s22"
        outArray(1+6+3)%varName="s33"
        outArray(1+6+4)%varName="s23"
        outArray(1+6+5)%varName="s31"
        outArray(1+6+6)%varName="s12"
        outArray(1+6+6+1)%varName="VM"
        outArray(1+6+6+2)%varName="phi"
        do idim = 1,(1+6+6+2)
            call VecGetArrayReadF90(gathered(idim), outArray(idim)%vardata, ierr)
            outArray(idim)%ifCell = .false.
        enddo
        if(rank == 0)then
            call output_plt_scalars(path, title, outArray)
        endif
        do idim = 1,(1+6+6+2)
            call VecRestoreArrayReadF90(gathered(idim), outArray(idim)%vardata, ierr)
            call VecDestroy(gathered(idim), ierr)
        enddo
    end subroutine

    subroutine output_plt_vert_vec3(path, title, outvec,name1,name2,name3)
        Vec veczero, outvec
        PetscInt ierr, rank
        PetscScalar, pointer :: gathered(:)
        character(*), intent(in) :: path, title, name1,name2,name3

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        if(rank == 0) then
            call VecCreateMPI(MPI_COMM_WORLD,globalDOFs*3,PETSC_DETERMINE, veczero, ierr)
        else
            call VecCreateMPI(MPI_COMM_WORLD,           0,PETSC_DETERMINE, veczero, ierr)
        endif

        call GatherVec3(outvec, veczero, .false.)
        call VecGetArrayReadF90(veczero, gathered, ierr)
        if(rank == 0) then
            call output_plt_scalar3x(path, title, gathered, name1, name2, name3, .false.)
        endif
        call VecRestoreArrayReadF90(veczero, gathered, ierr)
        call VecDestroy(veczero, ierr)
    end subroutine

    subroutine output_plt_vert_vec1(path, title, outvec,name)
        Vec veczero, outvec
        PetscInt ierr, rank
        PetscScalar, pointer :: gathered(:)
        character(*), intent(in) :: path, title,name

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        if(rank == 0) then
            call VecCreateMPI(MPI_COMM_WORLD,globalDOFs,PETSC_DETERMINE, veczero, ierr)
        else
            call VecCreateMPI(MPI_COMM_WORLD,         0,PETSC_DETERMINE, veczero, ierr)
        endif

        call GatherVec1(outvec, veczero, .false.)
        call VecGetArrayReadF90(veczero, gathered, ierr)
        ! gatheredUther = gathered
        if(rank == 0) then
            call output_plt_scalar(path, title, gathered, name, .false.)
        endif
        call VecRestoreArrayReadF90(veczero, gathered, ierr)
        call VecDestroy(veczero, ierr)
        ! if(rank == 0) then
        !     print*,'UtherGather:: ', allocated(gatheredUther)
        ! endif
    end subroutine

    subroutine deallocate_petscInt_vardimWrapper(obj)
        type(petscInt_vardimWrapper) :: obj(:)
        PetscInt i
        do i = 1, size(obj)
            if(allocated(obj(i)%N)) then
                deallocate(obj(i)%N)
            endif
        enddo
    end subroutine

    subroutine deallocate_petscScalar_vardimWrapper(obj)
        type(petscScalar_vardimWrapper) :: obj(:)
        PetscInt i
        do i = 1, size(obj)
            if(allocated(obj(i)%N)) then
                deallocate(obj(i)%N)
            endif
        enddo
    end subroutine

    !!!!!!!!!!deprecated!!!!!!!!!!
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !!!!!!!!!!deprecated!!!!!!!!!!
    subroutine SetUpThermal ! deprecated
        use globals
        PetscInt ierr
        PetscInt maxWidth
        PetscInt, pointer :: parray(:)
        PetscScalar, pointer :: pfix(:)
        integer :: rank, siz
        PetscInt i, j, ncells, nverts, elem_id, num_node, num_intpoint, ip, in, in2, &
            ie, ibc, ifc, facesize
        PetscInt faceNodes(8)
        PetscInt, allocatable::cellDOFs(:)
        PetscScalar, allocatable::cellMat(:,:), cellVec(:)
        real(8) :: elem_coord(27,3), Jacobi(3,3), invJacobi(3,3)  !max num node is 27 ! warning
        real(8) :: dNjdxi_ij(3,27), k__mul__dNjdxi_ij(3,27)
        real(8) minus1, nodeFix, nodeDiag
        logical touchedMat

        !
        ! Vec testvec, testvec2
        ! VecScatter procToZeroScatter
        ! PetscInt Gsize, csize
        ! PetscScalar cval

        !!! startup
        maxWidth = adjMaxWidth
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
        ncells = size(CELL)
        nverts = size(COORD,dim=2)
        call IsGetIndicesF90(partitionedNumberingIndex, parray, ierr) !parray is now partitionedNumberingIndex
        call VecGetArrayReadF90(dofFixTher, pfix, ierr)

        call VecSet(Pther,0.0_8,ierr)

        if(rank == 0)then
            do i = 1, ncells
                elem_id = getElemID(size(CELL(i)%N))
                num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                allocate(cellDOFs(num_node*1))
                allocate(cellMat(num_node*1,num_node*1))
                allocate(cellVec(num_node*1))
                do j = 1, num_node
                    cellDOFs(1*(j - 1) + 1) = parray(CELL(i)%N(j)) * 1 + 0
                enddo
                !!! Integrate the cellMat
                ! get the coords
                do in = 1, num_node
                    elem_coord(in,:) = COORD(:,CELL(I)%N(in))
                end do
                cellMat = 0.0_8
                do ip = 1, num_intpoint
                    Jacobi = matmul(elem_lib(elem_id)%dNIdLimr(:,:,ip),elem_coord(1:num_node,:)) !
                    invJacobi = directInverse3x3(Jacobi)
                    dNjdxi_ij(:,1:num_node) = matmul(invJacobi,elem_lib(elem_id)%dNIdLimr(:,:,ip))
                    k__mul__dNjdxi_ij(:,1:num_node) = matmul(k_ther, dNjdxi_ij)
                    cellMat = cellMat + matmul(transpose(dNjdxi_ij(:,1:num_node)), k__mul__dNjdxi_ij(:,1:num_node))
                end do

                touchedMat = .false.
                do in = 1, num_node
                    if(.not. isnan(pfix(CELL(I)%N(in)))) then
                        nodeFix = pfix(CELL(I)%N(in))
                        nodeDiag = cellMat(in,in)
                        cellVec = -cellMat(:,in)
                        do in2 = 1, num_node
                            if(.not. isnan(pfix(CELL(I)%N(in2)))) then
                                cellVec(in2) = 0.0_8 ! do not add rhs for other fixed dofs
                            endif
                        enddo
                        cellVec(in) = nodeDiag
                        cellVec = cellVec * nodeFix
                        call VecSetValues(Pther,num_node,cellDOFs,cellVec,ADD_VALUES,ierr)
                        cellMat(in,:) = 0.0_8
                        cellMat(:,in) = 0.0_8
                        cellMat(in,in) = nodeDiag
                        ! print*,"Fixing",CELL(I)%N(in), nodeFix
                        touchedMat = .true.

                    endif
                enddo

                ! if(touchedMat) then
                !     print*, "CELL",i
                !     do in2 = 1, num_node
                !         print*, cellMat(in2,:)
                !     enddo
                ! endif

                !!!
                call MatSetValues(Ather, num_node*1, cellDOFs, num_node*1, cellDOFs, cellMat, ADD_VALUES, ierr)
                !print*,i
                !print*,cellMat
                !print*,'CELLDOFS',cellDOFs
                deallocate(cellDOFs)
                deallocate(cellMat)
                deallocate(cellVec)
            enddo
        endif
        call VecRestoreArrayReadF90(dofFixTher, pfix, ierr)

        if(rank == 0)then
            do ibc = 1, NBSETS
                do i = 1, size(bcread(ibc)%ELEM_ID)

                enddo
            enddo
        endif

        call MatAssemblyBegin(Ather, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(Ather, MAT_FINAL_ASSEMBLY, ierr)
        call VecAssemblyBegin(Pther, ierr)
        call VecAssemblyEnd(Pther,ierr)
        !call MatView(Ather, PETSC_VIEWER_STDOUT_WORLD, ierr)
        !call VecView(Pther, PETSC_VIEWER_STDOUT_WORLD, ierr)

        ! call VecDuplicate(Pther, testvec, ierr)
        ! call VecGetSize(testvec, Gsize, ierr)
        ! if(rank == 0)then
        !     do i = 0,(Gsize-1)
        !         csize = i
        !         cval = i
        !         call VecSetValue(testvec, csize, cval, INSERT_VALUES, ierr)
        !     enddo
        ! endif
        ! call VecAssemblyBegin(testvec, ierr)
        ! call VecAssemblyEnd(testvec,ierr)

        ! if(rank == 0) then
        !     call VecCreateMPI(MPI_COMM_WORLD,Gsize,PETSC_DETERMINE, testvec2, ierr)
        ! else
        !     call VecCreateMPI(MPI_COMM_WORLD,0,PETSC_DETERMINE, testvec2, ierr)
        ! endif
        ! call VecScatterCreate(testvec, partitionedNumberingIndex, testvec2,PETSC_NULL_IS,&
        !                       procToZeroScatter, ierr)
        ! call VecScatterBegin(procToZeroScatter,testvec,testvec2,INSERT_VALUES,SCATTER_FORWARD,ierr)
        ! call VecScatterEnd(procToZeroScatter,testvec,testvec2,INSERT_VALUES,SCATTER_FORWARD,ierr)
        ! call VecScatterDestroy(procToZeroScatter, ierr)
        ! call VecView(testvec2, PETSC_VIEWER_STDOUT_WORLD, ierr)

        CHKERRA(ierr)
    end subroutine

    subroutine SetUpElasticity
    end subroutine

end module
