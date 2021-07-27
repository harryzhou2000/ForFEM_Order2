module common_utils
implicit none

contains
    subroutine print_real(a)
        real, intent(in) :: a
        write(*,)
500 format(E20)
    end subroutine print_real


end module common_utils