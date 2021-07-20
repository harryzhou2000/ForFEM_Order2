subroutine readgfile
    !subroutine input params
    !fileinp, string, as import file name
    use globals
    implicit none
    integer :: cety,cenode,I,m,J

    OPEN(IIN   , FILE = trim(fileinp),  STATUS = "OLD")
    read(IIN,'(6/6(1X,I9)//)') NUMNP,NELEM,NGRPS,NBSETS,NDFCD,NDFVL

    allocate(COORD(3,NUMNP))
    allocate(CELL(NELEM))
    allocate(CELLTYPE(NELEM))

    do I = 1,NUMNP
        read(IIN,'(10X,3E20.11)') COORD(1:3,I)
    enddo

    read(IIN,*)
    read(IIN,*)

    do I = 1,NELEM
        read(IIN,fmt=300,advance='no') cety,cenode
        call noderead(cety,cenode,i)
    enddo
300 format(9X,I2,1X,I2,1X)

    !   read zones
    allocate(ZONE(NGRPS)) !zone elem set
    allocate(zone_label(NGRPS)) !zone label strings
    allocate(GRPIT(NGRPS)) !the basic info of zones

    do I = 1,NGRPS
        allocate(GRPIT(I)%N(4))
        read(IIN,400) GRPIT(I)%N(1:4)

        read(IIN,'(A32)') zone_label(I)
        read(IIN,*)
        allocate(ZONE(I)%N(GRPIT(I)%N(2)))
        read(IIN,'(10I8)') (ZONE(I)%N(J),J = 1,GRPIT(I)%N(2))
    enddo
    !read bsets
    ALLOCATE(ITYPE(NBSETS))
    ALLOCATE(NENTRY(NBSETS))
    ALLOCATE(NVALUES(NBSETS))
    ALLOCATE(IBCODE(NBSETS))

    do m = 1,NBSETS
        read(IIN,500) bname(m),itype(m),nentry(m),nvalues(m),ibcode(m)

        allocate(bcread(m)%elem_id(nentry(m)))
        allocate(bcread(m)%elem_type(nentry(m)))
        allocate(bcread(m)%elem_face(nentry(m)))

        do i = 1,nentry(m)
            read(IIN,600) bcread(m)%elem_id(i),bcread(m)%elem_type(i),bcread(m)%elem_face(i)
        enddo
    enddo
    !
    ! ... output the bad cell
    !

    close(IIN)
400 format(2/7X,I10,11X,I10,11X,I10,9X,I10)
600 format(I10,I5,I5)
500 format(2/,a32,4i8)
end subroutine

subroutine noderead(cety,cenode,i)
    use globals
    implicit none
    integer :: cety,cenode,i,tem(27),rtem(27), ctype, cnodes, solver_Nnode

    ctype = cety
    cnodes = cenode
    read(IIN,fmt=200) tem(1:cenode)

200 format(7I8:/(15X,7I8:))

    ! ... HEX
    IF( CTYPE == 4 ) THEN ! HEX
        IF(CNODES==8) THEN
            RTEM(1) = TEM(1)
            RTEM(2) = TEM(5)
            RTEM(3) = TEM(6)
            RTEM(4) = TEM(2)
            RTEM(5) = TEM(3)
            RTEM(6) = TEM(7)
            RTEM(7) = TEM(8)
            RTEM(8) = TEM(4)

            SOLVER_NNODE = 8

        ELSEIF(CNODES==20) THEN

            ! CORNER POINT
            RTEM(1) = TEM( 1)
            RTEM(2) = TEM(13)
            RTEM(3) = TEM(15)
            RTEM(4) = TEM( 3)
            RTEM(5) = TEM( 6)
            RTEM(6) = TEM(18)
            RTEM(7) = TEM(20)
            RTEM(8) = TEM( 8)

            !EDGE CENTER POINT
            RTEM(10) = TEM(14)
            RTEM(12) = TEM( 2)
            RTEM(16) = TEM( 7)
            RTEM(14) = TEM(19)
            RTEM(18) = TEM(16)
            RTEM(19) = TEM(17)
            RTEM(20) = TEM( 5)
            RTEM(17) = TEM( 4)
            RTEM( 9) = TEM( 9)
            RTEM(11) = TEM(10)
            RTEM(15) = TEM(12)
            RTEM(13) = TEM(11)

            SOLVER_NNODE = 20

        ELSEIF(CNODES==27) THEN
            ! CORNER POINT
            RTEM(1) = TEM( 1)
            RTEM(2) = TEM(19)
            RTEM(3) = TEM(21)
            RTEM(4) = TEM( 3)
            RTEM(5) = TEM( 7)
            RTEM(6) = TEM(25)
            RTEM(7) = TEM(27)
            RTEM(8) = TEM( 9)

            ! EDGE CENTER POINT
            RTEM(10) = TEM(20)
            RTEM(12) = TEM( 2)
            RTEM(16) = TEM( 8)
            RTEM(14) = TEM(26)
            RTEM(18) = TEM(22)
            RTEM(19) = TEM(24)
            RTEM(20) = TEM( 6)
            RTEM(17) = TEM( 4)
            RTEM( 9) = TEM(10)
            RTEM(11) = TEM(12)
            RTEM(15) = TEM(18)
            RTEM(13) = TEM(16)

            SOLVER_NNODE = 20

            ! STORE SURFACE NODES FOR ADAPTIVITY
            RTEM(21) = TEM(14)
            RTEM(22) = TEM( 5)
            RTEM(23) = TEM(23)
            RTEM(24) = TEM(11)
            RTEM(25) = TEM(17)
            RTEM(26) = TEM(13)
            RTEM(27) = TEM(15)

        ENDIF

        ! ... TET
    ELSEIF( CTYPE == 6 ) THEN ! TET
        IF(CNODES==4) THEN
            RTEM(1) = TEM(1)
            RTEM(2) = TEM(4)
            RTEM(3) = TEM(3)
            RTEM(4) = TEM(2)
            SOLVER_NNODE = 4

        ELSEIF(CNODES==10) THEN
            RTEM(1) = TEM( 1)
            RTEM(2) = TEM(10)
            RTEM(3) = TEM( 6)
            RTEM(4) = TEM( 3)

            RTEM(5) = TEM( 7)
            RTEM(6) = TEM( 9)
            RTEM(7) = TEM( 4)
            RTEM(8) = TEM( 2)
            RTEM(9) = TEM( 8)
            RTEM(10) = TEM(5)
            SOLVER_NNODE = 10
        ENDIF

! ... WEDGE
    ELSEIF( CTYPE == 5 ) THEN ! WEDGE
        IF(CNODES==6) THEN
            RTEM( 1) = TEM( 1)
            RTEM( 2) = TEM( 2)
            RTEM( 3) = TEM( 3)
            RTEM( 4) = TEM( 4)
            RTEM( 5) = TEM( 5)
            RTEM( 6) = TEM( 6)

            SOLVER_NNODE = 6

        ELSEIF(CNODES==15) THEN
            RTEM( 1) = TEM( 1)
            RTEM( 2) = TEM( 3)
            RTEM( 3) = TEM( 6)
            RTEM( 4) = TEM(10)
            RTEM( 5) = TEM(12)
            RTEM( 6) = TEM(15)

            RTEM( 7) = TEM( 2)
            RTEM( 8) = TEM( 5)
            RTEM( 9) = TEM( 4)
            RTEM(10) = TEM(11)
            RTEM(11) = TEM(14)
            RTEM(12) = TEM(13)
            RTEM(13) = TEM( 7)
            RTEM(14) = TEM( 8)
            RTEM(15) = TEM( 9)

            SOLVER_NNODE = 15

        ELSEIF(CNODES==18) THEN
            RTEM( 1) = TEM( 1)
            RTEM( 2) = TEM( 3)
            RTEM( 3) = TEM( 6)
            RTEM( 4) = TEM(13)
            RTEM( 5) = TEM(15)
            RTEM( 6) = TEM(18)

            RTEM( 7) = TEM( 2)
            RTEM( 8) = TEM( 5)
            RTEM( 9) = TEM( 4)
            RTEM(10) = TEM(14)
            RTEM(11) = TEM(17)
            RTEM(12) = TEM(16)
            RTEM(13) = TEM( 7)
            RTEM(14) = TEM( 9)
            RTEM(15) = TEM(12)

            SOLVER_NNODE = 15

            ! STORE SURFACE NODES FOR ADAPTIVITY
            RTEM(16) = TEM( 8)
            RTEM(17) = TEM(10)
            RTEM(18) = TEM(11)

        ENDIF

! ... PYRAMAID
    ELSEIF( CTYPE == 7 ) THEN ! PRYMAID

        IF(CNODES==5) THEN
            RTEM( 1) = TEM( 1)
            RTEM( 2) = TEM( 2)
            RTEM( 3) = TEM( 4)
            RTEM( 4) = TEM( 3)
            RTEM( 5) = TEM( 5)

            SOLVER_NNODE = 5

        ELSEIF(CNODES==13) THEN
            RTEM( 1) = TEM( 1)
            RTEM( 2) = TEM( 3)
            RTEM( 3) = TEM( 8)
            RTEM( 4) = TEM( 6)
            RTEM( 5) = TEM(13)

            RTEM( 6) = TEM( 2)
            RTEM( 7) = TEM( 5)
            RTEM( 8) = TEM( 7)
            RTEM( 9) = TEM( 4)
            RTEM(10) = TEM( 9)
            RTEM(11) = TEM(10)
            RTEM(12) = TEM(12)
            RTEM(13) = TEM(11)

            SOLVER_NNODE = 13

        ELSEIF(CNODES==14) THEN
            RTEM( 1) = TEM( 1)
            RTEM( 2) = TEM( 3)
            RTEM( 3) = TEM( 9)
            RTEM( 4) = TEM( 7)
            RTEM( 5) = TEM(14)

            RTEM( 6) = TEM( 2)
            RTEM( 7) = TEM( 6)
            RTEM( 8) = TEM( 8)
            RTEM( 9) = TEM( 4)
            RTEM(10) = TEM(10)
            RTEM(11) = TEM(11)
            RTEM(12) = TEM(13)
            RTEM(13) = TEM(12)

            SOLVER_NNODE = 13
            ! STORE SURFACE NODES FOR ADAPTIVITY
            RTEM(14) = TEM(5)
        ENDIF
    ENDIF
    allocate(CELL(i)%N(SOLVER_NNODE))
    CELL(i)%N = RTEM(1:SOLVER_NNODE)
! ...
    return

end subroutine

subroutine writeelemg(I)
    use globals
    implicit none
    integer :: cetype,cenode,NV,I

    NV = size(CELL(I)%N)
    if (NV .eq. 4) then
        cetype = 6
        cenode = 4
        write(IOUT) cetype,cenode
300     format(8X,2I3,1X)
        write(IOUT) CELL(i)%N(1:4)
200     format(7I10:/(15X,7I10:))
    elseif (NV .eq. 10) then
        cetype = 6
        cenode = 10
        write(IOUT) cetype,cenode
        write(IOUT) CELL(i)%N(1:10)
    elseif (NV .eq. 5) then
        cetype = 7
        cenode = 5
        write(IOUT) cetype,cenode
        write(IOUT) CELL(I)%N(1:5)
    elseif (NV .eq. 13) then
        cetype = 7
        cenode = 13
        write(IOUT) cetype,cenode
        write(IOUT) CELL(I)%N(1:13)
    elseif(NV .eq. 6) then
        cetype = 5
        cenode = 6
        write(IOUT) cetype,cenode
        write(IOUT) CELL(I)%N(1:6)
    elseif (NV .eq. 15) then
        cetype = 5
        cenode = 15
        write(IOUT) cetype,cenode
        write(IOUT) CELL(I)%N(1:15)
    elseif (NV .eq. 8) then
        cetype = 4
        cenode = 8
        write(IOUT) cetype,cenode
        write(IOUT) CELL(I)%N(1:8)
    elseif (NV .eq. 20) then
        cetype = 4
        cenode = 20
        write(IOUT) cetype,cenode
        write(IOUT) CELL(I)%N(1:20)
    endif

    return

end subroutine

subroutine writebdface(m,I)
    use globals
    implicit none
    integer :: m,I,TY,FA

    TY = bcread(m)%elem_type(I)
    FA = bcread(m)%elem_face(I)
    if (TY .eq. 4) then
        call facetranhex(m,I,FA)
    elseif (TY .eq. 5) then
        call facetranwed(m,I,FA)
    elseif (TY .eq. 6) then
        call facetrantet(m,I,FA)
    elseif (TY .eq. 7) then
        call facetranpyr(m,I,FA)
    endif

    return

end subroutine

subroutine facetranhex(m,I,FA)
    use globals
    implicit none
    integer :: FA,m,I

    select case (FA)
    case (1)
        write(IOUT) bcread(m)%elem_id(I),4,1
    case (2)
        write(IOUT) bcread(m)%elem_id(I),4,3
    case(3)
        write(IOUT) bcread(m)%elem_id(I),4,6
    case(4)
        write(IOUT) bcread(m)%elem_id(I),4,5
    case(5)
        write(IOUT) bcread(m)%elem_id(I),4,4
    case(6)
        write(IOUT) bcread(m)%elem_id(I),4,2
    end select

    return

end subroutine

subroutine facetranwed(m,I,FA)
    use globals
    implicit none
    integer :: m,I,FA

    select case(FA)
    case(1)
        write(IOUT) bcread(m)%elem_id(I),5,4
    case(2)
        write(IOUT) bcread(m)%elem_id(I),5,3
    case(3)
        write(IOUT) bcread(m)%elem_id(I),5,2
    case(4)
        write(IOUT) bcread(m)%elem_id(I),5,1
    case(5)
        write(IOUT) bcread(m)%elem_id(I),5,5
    end select

    return

end subroutine

subroutine facetrantet(m,I,FA)
    use globals
    implicit none
    integer :: m,I,FA

    select case(FA)
    case(1)
        write(IOUT) bcread(m)%elem_id(I),6,2
    case(2)
        write(IOUT) bcread(m)%elem_id(I),6,4
    case(3)
        write(IOUT) bcread(m)%elem_id(I),6,3
    case(4)
        write(IOUT) bcread(m)%elem_id(I),6,1
    end select

    return

end subroutine

subroutine facetranpyr(m,I,FA)
    use globals
    implicit none
    integer :: m,I,FA

    select case(FA)
    case(1)
        write(IOUT) bcread(m)%elem_id(I),7,1
    case(2)
        write(IOUT) bcread(m)%elem_id(I),7,2
    case(3)
        write(IOUT) bcread(m)%elem_id(I),7,3
    case(4)
        write(IOUT) bcread(m)%elem_id(I),7,4
    case(5)
        write(IOUT) bcread(m)%elem_id(I),7,5
    end select

    return

end subroutine

subroutine initialize
    use globals
    implicit none
    integer :: m

    deallocate(COORD)
    deallocate(CELL)
    deallocate(ZONE)
    deallocate(zone_label)
    deallocate(GRPIT)
    deALLOCATE(ITYPE)
    deALLOCATE(NENTRY)
    deALLOCATE(NVALUES)
    deALLOCATE(IBCODE)
    do m = 1,NBSETS
        deallocate(bcread(m)%elem_id  )
        deallocate(bcread(m)%elem_type)
        deallocate(bcread(m)%elem_face)
    enddo

end subroutine


