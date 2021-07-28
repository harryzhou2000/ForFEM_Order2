!linear set

module linear_set
    implicit none

    type lset
        integer, pointer ::d(:)
        integer(4) :: siz
    end type

contains

    subroutine lsetCreate (S, cap0)
        implicit none
        type(lset) S
        integer(4) cap0
        allocate(S%d(cap0))
        S%siz=0
    end subroutine

    subroutine lsetPush (S,newElem)
        implicit none
        type(lset) S
        integer newElem, i
        integer, pointer :: newd(:)
        integer(4) oldCap

        do i = 1, S%siz
            if(S%d(i) == newElem) then
                return
            endif
        enddo
        oldCap = size(S%d,dim=1)
        S%siz = S%siz + 1
        if( S%siz > oldCap ) then
            allocate(newd(oldCap*2))
            newd(1:oldCap) = S%d
            newd(S%siz) = newElem
            deallocate(S%d)
            S%d => newd
        else
            S%d(S%siz) = newElem
        endif

    end subroutine

    subroutine int_mergeSort(Seq,lo,hi)
        integer,allocatable :: Seq(:)
        integer lo, hi
        integer,allocatable :: Useq(:)
        if(.not. lo < hi) then
            return
        endif
        allocate(Useq((lo + hi)/2-lo))
        call int_mergeSort_Rec(Seq, lo, hi, Useq)
        deallocate(Useq)
    end subroutine

    recursive subroutine int_mergeSort_Rec(Seq, lo, hi, Useq) ![lo,hi)
        integer,allocatable :: Seq(:)
        integer lo, hi, mid
        integer,allocatable :: Useq(:)

        if(.not. lo < hi - 1) then
            return
        endif
        mid = (lo + hi)/2
        call int_mergeSort_Rec(Seq, lo, mid, Useq)
        call int_mergeSort_Rec(Seq, mid, hi, Useq)
        call int_merge_inplace(Seq, lo,mid,hi,Useq)
    end subroutine

    subroutine int_merge_inplace(Seq,lo,mid,hi,Useq)
        integer,allocatable :: Seq(:)
        integer,allocatable :: Useq(:)
        integer lo, mid, hi, lop, midp, newp
        Useq(lo:mid-1) = Seq(lo:mid-1)
        lop = lo
        midp = mid
        newp = lo
        do while(lop < mid .and. midp < hi)
            if(Seq(midp) < Useq(lop))then
                Seq(newp) = Seq(midp)
                newp = newp + 1
                midp = midp + 1
            else
                Seq(newp) = Useq(lop)
                newp = newp + 1
                lop = lop + 1
            endif
        enddo
        do while(lop < mid)
            Seq(newp) = Useq(lop)
            newp = newp + 1
            lop = lop + 1
        enddo
        do while(midp < hi)
            Seq(newp) = Seq(midp)
            newp = newp + 1
            midp = midp + 1
        enddo
    end subroutine

end module
