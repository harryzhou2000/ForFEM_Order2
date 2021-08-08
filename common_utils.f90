
#define scalar real(kind=8)
#define index integer(kind=4)

module common_utils
implicit none

contains
    function directInverse3x3(A) result(AI)
        scalar, intent(in)  :: A(3, 3)
        scalar :: AI(3, 3)
        scalar :: detA
        detA = A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) - A(1, 2)*A(2, 1)*A(3, 3) &
               + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1)
        if (abs(detA) < tiny(detA)*16) then
            print *, "Error::directInverse3x3::singular matrix, det = ", detA
            stop
        end if
        AI(1, 1) = (A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))/detA
        AI(1, 2) = -(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))/detA
        AI(1, 3) = (A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))/detA
        AI(2, 1) = -(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))/detA
        AI(2, 2) = (A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))/detA
        AI(2, 3) = -(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))/detA
        AI(3, 1) = (A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))/detA
        AI(3, 2) = -(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))/detA
        AI(3, 3) = (A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))/detA
    end function

    function directDet3x3(A) result(detA)
        scalar, intent(in)  :: A(3, 3)
        scalar :: detA
        detA = A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) - A(1, 2)*A(2, 1)*A(3, 3) &
               + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1)
    end function

    function crossProduct(a,b) result(c)
        scalar a(3),b(3),c(3)
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function

    function getMinMax3x(DATAin) result(minmax)
        index:: siz !size is 1/3 size of DATAin 
        scalar, intent(in)  ::  DATAin(:)
        scalar ::  minmax(6)
        index i
        minmax(1) =  huge(minmax(1))
        minmax(2) = -huge(minmax(1))
        minmax(3) =  huge(minmax(1))
        minmax(4) = -huge(minmax(1))
        minmax(5) =  huge(minmax(1))
        minmax(6) = -huge(minmax(1))
        siz = size(DATAin)
        if(mod(siz,3) .ne. 0)then
            print*,"ERROR::common_utils::getMinMax3x, Datasize not 3x but", siz
            stop
        endif
        siz = siz/3

        do i = 1,siz
            if(DATAin((i-1) * 3 + 1) < minmax(1)) then
                minmax(1) = DATAin((i-1) * 3 + 1)
            endif
            if(DATAin((i-1) * 3 + 1) > minmax(2)) then
                minmax(2) = DATAin((i-1) * 3 + 1)
            endif
            if(DATAin((i-1) * 3 + 2) < minmax(3)) then
                minmax(3) = DATAin((i-1) * 3 + 2)
            endif
            if(DATAin((i-1) * 3 + 2) > minmax(4)) then
                minmax(4) = DATAin((i-1) * 3 + 2)
            endif
            if(DATAin((i-1) * 3 + 3) < minmax(5)) then
                minmax(5) = DATAin((i-1) * 3 + 3)
            endif
            if(DATAin((i-1) * 3 + 3) > minmax(6)) then
                minmax(6) = DATAin((i-1) * 3 + 3)
            endif
        enddo

    end function

    function myMinusNan() result(res)
        scalar res, minus1
        minus1 = -1
        res = sqrt(minus1)
    end function
    


end module common_utils

#undef scalar
#undef index