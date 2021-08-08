program test_set
    use linear_set
    implicit none
    type(lset) myset
    integer,allocatable,target ::seq(:), seq3(:)
    integer rseedsize, i,j,nhi
    integer,allocatable::rseed(:)
    integer,pointer::seq2(:)
    real rout
    integer,allocatable :: array(:)
    integer,allocatable :: anotherarray(:)
    integer place
    logical found
    real vec1(4), vec2(3)

    allocate(seq(3000))
    seq2 => seq
    call random_seed(size = rseedsize)
    allocate(rseed(rseedsize))
    rseed = 124
    call random_seed(put = rseed)

    do j = 1,100
        do i=1, 3000
            call random_number(rout)
            seq(i) = floor(rout*55)
        enddo
        call int_mergeSort(seq,1,3001)
        if(.not. int_checkSorted(seq2,1,3001)) then
            print*,"BAD ", j
        endif
    end do

    call int_reduceSorted(seq2,1,3001,nhi)
    ! print*,nhi
    ! print*,seq2(1:nhi-1)

    allocate(seq3(10))
    seq3 = (/1,3,4,1,2,4,5,8,11,9/)
    call int_mergeSort(seq3,1, 11)
    call int_reduceSorted(seq3,1,11,nhi)
    print*,'seq3',seq3(1:nhi-1)
    found = int_searchBinary(seq3,1,nhi,9,place)
    print*,found,place
    


    call lsetCreate(myset,10)
    do j = 1,1000
        call random_number(rout)
        call lsetPush(myset,floor(rout*1e4))

    end do
    !print*,myset%d(1:myset%siz)

    allocate(       array(1000000))
    allocate(anotherarray(1000000))
    do i = 1,10
        anotherarray = array
    end do

    vec1 = (/1,2,3,4/)
    vec2 = (/5,6,7/)
    

end program test_set


