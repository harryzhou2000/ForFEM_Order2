program test_set
    use linear_set
    implicit none
    type(lset) myset
    integer,allocatable ::seq(:)

    call lsetCreate(myset, 5)
    call lsetPush(myset, 1)
    call lsetPush(myset, 2)
    call lsetPush(myset, 3)
    call lsetPush(myset, 4)
    call lsetPush(myset, 5)
    call lsetPush(myset, 6)
    call lsetPush(myset, 6)
    call lsetPush(myset, 8)
    call lsetPush(myset, 6)
    call lsetPush(myset, 8)
    call lsetPush(myset, 1)
    call lsetPush(myset, 7)
    !print*,myset%d

    allocate(seq(17),source=[&
             14,5,1,3,&
             141,15,121,13,&
             1311,13,1,4,&
             155,14,5,1414,&
             77&
             ])
    call int_mergeSort(seq,1,18)
    print*,seq

end program test_set
