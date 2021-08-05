!linear set
!int is to be decided
#define int integer

module linear_set
    implicit none

    type lset
        int, pointer ::d(:)
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
        int newElem, i
        int, pointer :: newd(:)
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
        int:: Seq(:)
        int lo, hi
        int,pointer :: Useq(:)
        if(.not. lo < hi) then
            return
        endif
        allocate(Useq((lo + hi)/2-lo+1))
        call int_mergeSort_Rec(Seq, lo, hi, Useq)
        deallocate(Useq)
    end subroutine

    recursive subroutine int_mergeSort_Rec(Seq, lo, hi, Useq) ![lo,hi)
        int :: Seq(:)
        int lo, hi, mid
        int :: Useq(:)

        if(.not. lo < hi - 1) then
            return
        endif
        mid = (lo + hi)/2
        call int_mergeSort_Rec(Seq, lo, mid, Useq)
        call int_mergeSort_Rec(Seq, mid, hi, Useq)
        call int_merge_inplace(Seq, lo,mid,hi,Useq)
    end subroutine

    subroutine int_merge_inplace(Seq,lo,mid,hi,Useq)
        int :: Seq(:)
        int :: Useq(:)
        int lo, mid, hi, lop, midp, newp
        Useq(1:mid-lo) = Seq(lo:mid-1)
        lop = lo
        midp = mid
        newp = lo
        do while(lop < mid .and. midp < hi)
            if(Seq(midp) < Useq(lop-lo+1))then
                Seq(newp) = Seq(midp)
                newp = newp + 1
                midp = midp + 1
            else
                Seq(newp) = Useq(lop-lo+1)
                newp = newp + 1
                lop = lop + 1
            endif
        enddo
        do while(lop < mid)
            Seq(newp) = Useq(lop-lo+1)
            newp = newp + 1
            lop = lop + 1
        enddo
        do while(midp < hi)
            Seq(newp) = Seq(midp)
            newp = newp + 1
            midp = midp + 1
        enddo
    end subroutine

    function int_checkSorted(Seq,lo,hi) result(res) ! true if sorted
        int Seq(:)
        int lo,hi,i
        logical res
        res = .true.
        do i = lo+1,hi-2
            if (Seq(lo) > Seq(lo + 1)) then
                res = .false.
            end if
        end do
    end function

    subroutine int_reduceSorted(Seq,lo,hi,nhi)
        int :: Seq(:)
        int lo, hi
        int pwrite, psee
        int nhi
        pwrite = lo + 1
        psee = lo + 1

        do while(psee < hi)
            if(Seq(psee-1)==Seq(psee)) then
                psee = psee + 1
                cycle
            else
                Seq(pwrite)=Seq(psee)
                psee = psee + 1
                pwrite = pwrite + 1
                cycle
            endif
        enddo
        nhi = pwrite
        !fill null
        do while(pwrite < hi)
            Seq(pwrite) = huge(1_4)
            pwrite = pwrite + 1
        enddo
    end subroutine


    

end module

#undef int