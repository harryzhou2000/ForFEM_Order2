module globals
    implicit none
!
    CHARACTER*80 FILEINP1,FILEINP2,FILEINP,FILEOUT
!   units definition
    integer :: IIN = 11     ! unit for input file
    integer :: IOUT = 12    ! unit for output file
    integer :: IOUT2 = 13   ! unit for output file
    integer :: TEMC = 14    ! unit for tem file
    integer :: FAT = 15     ! unit for hierachy
    INTEGER :: SPLIT_ACTION
!   variables definition
!   GAMBIT
    CHARACTER(79) :: GFILE ! I/O
    INTEGER :: NUMNP,NELEM,NINNER,NBOUND,NNELEM,NNBOUND
    INTEGER :: NGRPS,NBSETS,NDFCD,NDFVL
!   Mesh Parameters
    INTEGER :: NNODE,NEDGE,NFACE,NCELL,NFACE3I,NFACE4I,NFACE3B,NFACE4B
    INTEGER :: NTRI,NQUAD,NTET,NHEX,NPRIS,NPRY,NFN,NVN,NVNC,NNP
    integer,allocatable :: cell_type(:),cell_nnode(:),ITYPE(:),NENTRY(:),NVALUES(:)
    INTEGER,ALLOCATABLE :: IBCODE(:),NNENTRY(:)
    CHARACTER(len=5)  :: MESHTYPE
    CHARACTER(32) :: bname(10)

!   Define new type
    TYPE VAR_DIM
        INTEGER,ALLOCATABLE :: N(:)
    END TYPE VAR_DIM

    TYPE FACE_S
        INTEGER :: CE1,CE2
        INTEGER,ALLOCATABLE :: N(:)
    END TYPE FACE_S

    TYPE BC_READ
        integer,allocatable :: ELEM_ID(:)
        integer,allocatable :: ELEM_TYPE(:)
        integer,allocatable :: ELEM_FACE(:)
        integer :: BC
    end type

    type EDGE3Point
        integer :: L1(2)
        integer :: L2(2)
    end type

    type EDGE3F
        integer :: F(2)
    end type

    type FACE3F
        integer :: P(4)
    end type

!   EDGE DEFINATION
    INTEGER,PARAMETER :: EDGE_TET(12) = &
                         reshape((/1,2,2,3,3,1,1,4,2,4,3,4/),(/12/))
    INTEGER,PARAMETER :: EDGE_TETC(6) = (/5,6,7,8,9,10/)
    INTEGER,PARAMETER :: EDGE_HEX(24) = &
                         reshape((/1,2,2,3,3,4,4,1,5,1,6,2,7,3,8,4,&
                                   5,6,6,7,7,8,8,5/),(/24/))
    INTEGER,PARAMETER :: EDGE_HEXC(12) = &
                         (/9,10,11,12,13,14,15,16,17,18,19,20/)
    INTEGER,PARAMETER :: EDGE_WED(18) = &
                         reshape((/1,2,2,3,3,1,1,4,2,5,3,6,4,5,5,6,6,4/),(/18/))
    INTEGER,PARAMETER :: EDGE_WEDC(9) = &
                         (/7,8,9,10,11,12,13,14,15/)
    INTEGER,PARAMETER :: EDGE_PRY(16) = &
                         reshape((/1,2,2,3,3,4,4,1,1,5,2,5,3,5,4,5/),(/16/))
    INTEGER,PARAMETER :: EDGE_PRYC(8) = &
                         (/6,7,8,9,10,11,12,13/)

    real(8),allocatable :: COORD(:,:)
    real(8),allocatable :: NEW_COORD(:,:)
    type(VAR_DIM),allocatable,target :: CELL(:)
    type(VAR_DIM),allocatable,target :: NEW_CELL(:)
    type(VAR_DIM),allocatable,target :: FINCELL(:)
    type(VAR_DIM),allocatable,target :: CELLSNODE(:)
    type(VAR_DIM),allocatable,target :: EDGE(:)
    type(VAR_DIM),allocatable,target :: CELLEDGE(:)
    type(VAR_DIM),allocatable,target :: EDGENODE(:)
    type(VAR_DIM),allocatable,target :: CELLFACE(:)
    type(VAR_DIM),allocatable,target :: CELLTEMFACE(:)
    type(VAR_DIM),allocatable,target :: SUBELEM(:)
    TYPE(VAR_DIM),ALLOCATABLE,target :: ZONE(:)
    TYPE(VAR_DIM),ALLOCATABLE,target :: NZONE(:)
    type(VAR_DIM),ALLOCATABLE,target :: GRPIT(:)
    type(VAR_DIM),ALLOCATABLE,target :: NGRPIT(:)
    type(VAR_DIM),ALLOCATABLE,TARGET :: CELLTYPE(:)
    type(BC_READ) :: bcread(10),new_bcread(10)

    type(FACE_S),allocatable :: FACE(:)
    integer,allocatable :: EDGEF(:)
    integer,allocatable :: FACEF(:)
    type(EDGE3Point),allocatable :: EDGE3(:)
    type(EDGE3F),allocatable :: EDGEF3(:)
    type(FACE3F),allocatable :: FACEF3(:)
    character(32),allocatable:: zone_label(:)
    integer :: kindn(8)
!   time record
    real,allocatable :: TIM(:)
    integer :: NITER,level

end module
