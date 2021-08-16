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

    !unused
    subroutine lsetCreate (S, cap0)
        implicit none
        type(lset) S
        integer(4) cap0
        allocate(S%d(cap0))
        S%siz=0
    end subroutine

    !unused
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

    ! main entrance for int mergesort
    ! sort Seq's interval [lo,hi)
    subroutine int_mergeSort(Seq,lo,hi)
        int:: Seq(:)
        int,intent(in) ::  lo, hi
        int,pointer :: Useq(:)
        if(.not. lo < hi) then
            return
        endif
        allocate(Useq((lo + hi)/2-lo+1))
        call int_mergeSort_Rec(Seq, lo, hi, Useq)
        deallocate(Useq)
    end subroutine

    ! recursed function called by int_mergeSort
    recursive subroutine int_mergeSort_Rec(Seq, lo, hi, Useq) ![lo,hi)
        int :: Seq(:)
        int,intent(in) ::  lo, hi
        int mid
        int :: Useq(:)

        if(.not. lo < hi - 1) then
            return
        endif
        mid = lo + (hi - lo) / 2
        call int_mergeSort_Rec(Seq, lo, mid, Useq)
        call int_mergeSort_Rec(Seq, mid, hi, Useq)
        call int_merge_inplace(Seq, lo,mid,hi,Useq)
    end subroutine

    ! merge [lo,mid) with [mid,hi) with assistance of Useq
    subroutine int_merge_inplace(Seq,lo,mid,hi,Useq)
        int :: Seq(:)
        int :: Useq(:)
        int, intent(in) :: lo, mid, hi
        int lop, midp, newp
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

    ! check if [lo,hi) is sorted
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

    ! reduce dulplicates int a sorted [lo hi)
    ! return with the [lo,nhi) valid interval
    subroutine int_reduceSorted(Seq,lo,hi,nhi)
        int :: Seq(:)
        int,intent(in) :: lo, hi
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

    ! search tar within [loin,hiin), mid is the final place
    ! success is true if the final place equals tar
    function int_searchBinary(Seq,loin,hiin,tar,mid) result(success)
        int :: Seq(:)
        int, intent(in) :: loin, hiin ,tar
        int lo,hi,mid
        logical success
        !!
        lo = loin
        hi = hiin
        do while(lo < hi - 1)
            mid = lo + (hi - lo) / 2
            if(Seq(mid) > tar)then
                hi = mid
            else
                lo = mid
            endif
        enddo
        mid = lo
        success = .false.
        if(Seq(mid) == tar) then
            success = .true.
        endif
    end function
end module

#undef int