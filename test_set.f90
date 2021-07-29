program test_set
    use linear_set
    implicit none
    type(lset) myset
    integer,allocatable,target ::seq(:)
    integer rseedsize, i,j,nhi
    integer,allocatable::rseed(:)
    integer,pointer::seq2(:)
    real rout
    integer,allocatable :: array(:)
    integer,allocatable :: anotherarray(:)

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
        call int_mergeSort(seq2,1,3001)
        if(.not. int_checkSorted(seq2,1,3001)) then
            print*,"BAD ", j
        endif
    end do

    call int_reduceSorted(seq2,1,3001,nhi)
    print*,nhi
    print*,seq2(1:nhi-1)


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

end program test_set


