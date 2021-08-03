

module common_utils
implicit none

contains
    function directInverse3x3(A) result(AI)
        real(kind=8), intent(in)  :: A(3, 3)
        real(kind=8) :: AI(3, 3)
        real(kind=8) :: detA
        detA = A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) - A(1, 2)*A(2, 1)*A(3, 3) &
               + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1)
        if (abs(detA) < tiny(detA)*16) then
            print *, "Error===directInverse3x3===singular matrix, det = ", detA
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
        real(kind=8), intent(in)  :: A(3, 3)
        real(kind=8) :: detA
        detA = A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) - A(1, 2)*A(2, 1)*A(3, 3) &
               + A(1, 2)*A(2, 3)*A(3, 1) + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1)
    end function

    function crossProduct(a,b) result(c)
        real(8) a(3),b(3),c(3)
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function
    


end module common_utils