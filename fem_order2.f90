!defining global variables for
!2nd order classsic isoparametric FEM

module fem_order2

#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
!#include "petsc/mpiuni/mpif.h"
!#include "mpif.h"
    use petscvec
    use petscmat
    use petscsys
    use petscis
    use petscisdef
    use mat_csr

    !use mpi

    !use mpi
    implicit none
    MatPartitioning partition
    Mat adjmat
    type(tIS) partitionIndex, rowsTwoSidedIndex, &
        partitionedNumberingIndex, partitionedNumberingIndexLoc

    type(tIs) partitionedNumberingIndexInversed

    PetscInt :: localDOFs ! corresponding to rowsTwoSidedIndex, one dim per point
    PetscInt :: globalDOFs

    type(tMat) :: Aelas !stiffness for elasticity
    type(tVec) :: Pelas !load Vector for elasticity

    type(tMat) :: Ather !stiffness for thermal
    type(tVec) :: Pther !load Vector for thermal

    !topology of mesh nodes, should be parallized or upgraded to PETSC's dmplex
    !serial for proc0
    type(csr), target :: AdjacencyCounter
    !node->nodes counter uncompressed
    !serial for proc0
    integer, allocatable :: AdjacencyNum0(:)

    !using proc0, we assemble A and P

    integer, allocatable :: Ather_r_pos

    type fem_element
        integer(kind=4)           :: num_node ! M
        integer(kind=4)           :: num_intpoint ! R
        real(kind=8), allocatable :: coord_intpoint(:, :) !coords of integration point, 3xR
        real(kind=8), allocatable :: coord_intweight(:)  !weights of integration point, R
        real(kind=8), allocatable :: NImr(:, :) !N at intpoints M*R
        real(kind=8), allocatable :: dNIdLimr(:, :, :) !dNdL at intpoints 3*M*R
    end type

    integer, parameter:: nlib = 3
    type(fem_element) elem_lib(nlib)

    real(8), allocatable :: cell_volumes(:)
    real(8), allocatable :: point_partition(:)

contains
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
        tet_ws(4) =  9.0_8/20.0_8  /6.0_4
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
        lin_as(3) = lin_a2
        lin_as(4) = lin_a1
        lin_ws(1) = lin_w1
        lin_ws(2) = lin_w2
        lin_ws(3) = lin_w2
        lin_ws(4) = lin_w1

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

    end subroutine

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
                elem_id = getElemID(CELL(i)%N)
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
            do i = 1, ncells
                elem_id = getElemID(CELL(i)%N)
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
    subroutine getVolumes()
        use globals
        integer i, ip, in, elem_id, num_intpoint, num_node
        real(8) :: elem_coord(27,3) !max num node is 27 ! warning
        real(8) :: Jacobi(3,3)
        real(8) :: detJacobi

        allocate(cell_volumes(size(CELL)))
        do i = 1, size(CELL)
            elem_id = getElemID(CELL(i)%N)
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

    function directInverse3x3(A) result(AI)
        real(kind=8), intent(in)  :: A(3, 3)
        real(kind=8) :: AI(3, 3)
        real(kind=8) :: detA
        detA = A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) - A(1, 2)*A(2, 1)*A(3, 3) &
               + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1)
        if (abs(detA) < tiny(detA)*16) then
            print *, "Error===directInverse3x3===singular matrix, det = ", detA
            stop
        end if
        AI(1, 1) = (A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))/detA
        AI(1, 2) = -(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))/detA
        AI(1, 3) = (A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))/detA
        AI(2, 1) = -(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))/detA
        AI(2, 2) = (A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))/detA
        AI(2, 3) = -(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))/detA
        AI(3, 1) = (A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))/detA
        AI(3, 2) = -(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))/detA
        AI(3, 3) = (A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))/detA
    end function

    function directDet3x3(A) result(detA)
        real(kind=8), intent(in)  :: A(3, 3)
        real(kind=8) :: detA
        detA = A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) - A(1, 2)*A(2, 1)*A(3, 3) &
               + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1)
    end function

    function getElemID(CELLI) result(ID)
        integer, allocatable, intent(in) :: CELLI(:)
        integer :: ID
        select case(size(CELLI))
        case(15)
            ID = 1
        case(20)
            ID = 2
        case(10)
            ID = 3
        case default
            print *, "Error::getElemID::Current Num Node ",size(CELLI)
            print *, "Not supported"
            stop
        end select

    end function

    !!!!Tecplot output

    ! A seq routine, only 1 process should go
    !https://tecplot.azureedge.net/products/360/current/360_data_format_guide.pdf
    subroutine output_plt_mesh(path, title)
        use globals

        implicit none
        character(80), intent(in) :: path
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
        print *,"===Writing Mesh, num point ", size(COORD,2), " num elem ", size(CELL)

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
    end subroutine

    ! A seq routine, only 1 process should go
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
        print *,"===Writing Data ",DATAname,", num point ", size(COORD,2), " num elem ", size(CELL)

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

! interface
!     subroutine output_plt_mesh(path, COORD)
!     use iso_c_binding
!     real(kind=c_float), allocatable :: coord(:,:)

!     end subroutine output_plt_mesh
! end interface

!!!!!!!!!!!!!!!!! subroutines invoking PETSC

    subroutine SetUpPartition
        use globals
        use linear_set

        integer ncells, nverts,  elem_id, num_intpoint, num_node, current_row_start, current_ivert,&
            next_row_start, new_row_end, current_row_start_old, next_row_start_old
        integer rank, siz, ierr
        integer i,j
        integer fillBuffer(27)
        integer, pointer:: tobesorted(:)
        integer, allocatable:: oldCols(:), oldRows(:)
        PetscInt, pointer :: parray(:)
        Mat m2, m3
        Vec v1,v2
        PetscInt msizem, msizen
        integer loc(1)
        PetscInt maxWidth
        PetscScalar, allocatable::cellMat(:,:)
        PetscInt,allocatable:: cellDOFs(:)

        !!!!!!!!!!!!!!!!!!!!!!!!!
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)

        ncells = size(CELL)
        nverts = size(COORD,dim=2)
        if(rank==0) then
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
                elem_id = getElemID(CELL(i)%N)
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
                elem_id = getElemID(CELL(i)%N)
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

        !!! get Inumbering
        call ISCreate(MPI_COMM_WORLD,partitionedNumberingIndexInversed, ierr)
        call ISInvertPermutation(partitionedNumberingIndex,PETSC_DECIDE, partitionedNumberingIndexInversed,ierr)

        !!! call ISSetBlockSize(partitionedNumberingIndex,123, ierr)

        !!!!CHECK
        !call IsView(partitionedNumberingIndex, PETSC_VIEWER_STDOUT_WORLD, ierr)
        ! call MatCreateSubMatrix(adjMat,rowsTwoSidedIndex,rowsTwoSidedIndex,MAT_INITIAL_MATRIX,m2,ierr)
        ! !call MatPermute(adjMat,rowsTwoSidedIndex,rowsTwoSidedIndex,m3, ierr)
        ! !call MatView(m2, PETSC_VIEWER_STDOUT_WORLD, ierr)
        ! call MatGetLocalSize(m2, msizem, msizen, ierr)
        ! print*, rank, "m=", msizem, "n=", msizen
        call IsGetIndicesF90(partitionIndex, parray, ierr) ! Fortran90 call does not need restoring
        point_partition = parray
        if(rank ==0 ) then
            call output_plt_scalar('./out2_data1_part.plt', 'goodstart', point_partition, 'partition', .false.)
        endif
        !!!!CHECK

        !!! max Row Size
        if (rank == 0) then
            loc = maxloc(AdjacencyNum0)
            print*,"Max", AdjacencyNum0(loc(1))
            maxWidth = AdjacencyNum0(loc(1))
            loc = minloc(AdjacencyNum0)
            print*,"Min", AdjacencyNum0(loc(1)), size(AdjacencyNum0), loc
        endif
        call MPI_Bcast(maxWidth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        !print*,rank, "MW=",maxWidth

        !!! create a mat for testeing
        call MatCreate(MPI_COMM_WORLD, m3, ierr)
        call MatSetType(m3,MATMPIAIJ,ierr)
        call MatSetSizes(m3, localDOFs*3, localDOFs*3, globalDOFs*3, globalDOFs*3, ierr)
        call MatMPIAIJSetPreallocation(m3,maxWidth*3,PETSC_NULL_INTEGER,maxWidth*3,PETSC_NULL_INTEGER, ierr)
        call IsGetIndicesF90(partitionedNumberingIndex, parray, ierr) !parray is now partitionedNumberingIndex
        if(rank == 0)then
            do i = 1, ncells
                elem_id = getElemID(CELL(i)%N)
                !num_intpoint = elem_lib(elem_id)%num_intpoint
                num_node = elem_lib(elem_id)%num_node
                allocate(cellDOFs(num_node*3))
                allocate(cellMat(num_node*3,num_node*3))
                do j = 1, num_node
                    cellDOFs(3*(j - 1) + 1) = parray(CELL(i)%N(j)-1) * 3 + 0
                    cellDOFs(3*(j - 1) + 2) = parray(CELL(i)%N(j)-1) * 3 + 1
                    cellDOFs(3*(j - 1) + 3) = parray(CELL(i)%N(j)-1) * 3 + 2
                enddo
                cellMat = 1
                call MatSetValues(m3, num_node*3, cellDOFs, num_node*3, cellDOFs, cellMat, ADD_VALUES, ierr)
                deallocate(cellDOFs)
                deallocate(cellMat)
            enddo
        endif
        call MatAssemblyBegin(m3, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(m3, MAT_FINAL_ASSEMBLY, ierr)
        !call MatView(m3, PETSC_VIEWER_STDOUT_WORLD, ierr)
        call VecCreateMPI(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3, v1, ierr)
        call VecCreateMPI(MPI_COMM_WORLD, localDOFs*3, globalDOFs*3, v2, ierr)
        call VecSet(v1,1.0_8,ierr)
        call MatMult(m3, v1, v2, ierr)
        call VecView(v2,PETSC_VIEWER_STDOUT_WORLD, ierr)
        !!!

    end subroutine

    !after readgrid initializeLib
    subroutine SetUpThermal
        PetscInt ierr

        call MatCreate(PETSC_COMM_WORLD, ATher ,ierr)
        CHKERRA(ierr)
    end subroutine

    !after readgrid initializeLib
    subroutine SetUpElasticity

    end subroutine

    subroutine ClearThermal

    end subroutine

    subroutine ClearElasticity

    end subroutine

    subroutine ConstructThermalSeq

    end subroutine

    subroutine ConstructElasticitySeq

    end subroutine

    subroutine SolveThermal

    end subroutine

    subroutine SolveElasticity

    end subroutine

end module
