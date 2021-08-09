program test_set
    use linear_set
    use common_utils
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
    real(8) Mat(4,6),Mat0(4,6)
    real(8) x(6), b(6),s(4), work(4096)
    integer rankout,infoout,iworksize, iwork(44)

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

    Mat = reshape((/0.8147,0.9058,0.127,0.9134,0.6324,0.09754,0.2785,0.5469,0.9575,0.9649,&
    0.1576,0.9706,0.9572,0.4854,0.8003,0.1419,0.4218,0.9157,0.7922,0.9595,0.6557,0.03571,0.8491,0.934/),(/4,6/))
    Mat0 = Mat
    b = 1
    call printMat(Mat)
    print*,b
    call dgelss(4, 6, 1, Mat, 4, b, 6, S, -1.0_8,rankout, work,-1, infoout)
    print*,'query'
    print*,work(1)
    print*,'done'
    call dgelss(4, 6, 1, Mat, 4, b, 6, S, -1.0_8,rankout, work,4096, infoout)
    print*,b
    !print*,matmul(reshape(transpose(Mat),(/4,6/)),b)
    print*,matmul(Mat0,b)
    print*,S

end program test_set


