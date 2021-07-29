! TO NOT use para_csr in fem
! 0. get a nodal partition and local->global indice
! 1. distribute local cell->node info (none local node is -n as local index)
! 2. inquire nodal values: position bc source
! 2. determine DIAGONAL part and none-DIAGONAL part size=> petsc solve
! 3. solve the gradient: use dgels or dgelss solve least-square nodal gradient
!    3.1 average nodal values (push other nodes)




 

module mat_csr
    implicit none

    type Csr
        integer nrow, nnz !number of local rows and non-zeros, set with subroutines
        integer globalRowStart !where row 0 maps to the global row
        integer, allocatable :: rowStart(:)  ! set with subroutines
        integer, allocatable :: column(:)    ! set with subroutines
        real(8), allocatable :: value(:)     ! set with subroutines
        ! integer, allocatable :: OrowStart(:) ! set with subroutines
        ! integer, allocatable :: Ocolumn(:)   ! set with subroutines
    end type

contains
    ! resets the matrix to nnrow rows
    subroutine Csr_SetNROW(A,nnrow)
        type(Csr),intent(inout) :: A
        integer, intent(in) :: nnrow
        A%nrow = nnrow
    end subroutine

    subroutine Csr_SetNNZ(A,nnnz)
        type(Csr),intent(inout) :: A
        integer, intent(in) :: nnnz
        A%nnz = nnnz
    end subroutine

    subroutine Csr_AllocateRowStart(A)
        type(Csr),intent(inout) :: A
        if (allocated(A%rowStart)) then
            deallocate(A%rowStart)
        endif
        allocate(A%rowStart(A%nrow+1))
    end subroutine

    subroutine Csr_AllocateColumn(A)
        type(Csr),intent(inout) :: A
        if (allocated(A%column)) then
            deallocate(A%column)
        endif
        allocate(A%column(A%nnz))
    end subroutine

    subroutine Csr_AllocateValue(A)
        type(Csr),intent(inout) :: A
        if (allocated(A%value)) then
            deallocate(A%value)
        endif
        allocate(A%value(A%nnz))
    end subroutine

    subroutine Csr_DeleteMat(A)
        type(Csr),intent(inout) :: A
        if (allocated(A%value)) then
            deallocate(A%value)
        endif
        if (allocated(A%column)) then
            deallocate(A%column)
        endif
        if (allocated(A%rowStart)) then
            deallocate(A%rowStart)
        endif
    end subroutine

end module

