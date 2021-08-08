#define Scalar real(8)

module elastic_constitution
    implicit none

contains
    function constitutionMat(E,nu) result(ret)
        Scalar E, nu
        Scalar ret(0:5,0:5)
        ret = 0.0_8
        ret(0, 0) = 1.0
        ret(1, 1) = 1.0
        ret(2, 2) = 1.0
        ret(5, 5) = (1.0 - 2.0 * nu) / ((1.0 - nu))! be ware of gamma12 and e12
        ret(4, 4) = ret(5, 5)
        ret(3, 3) = ret(4, 4)
        ret(2, 1) = nu / (1.0 - nu)
        ret(2, 0) = ret(2, 1)
        ret(1, 0) = ret(2, 0)
        ret(1, 2) = ret(1, 0)
        ret(0, 2) = ret(1, 2)
        ret(0, 1) = ret(0, 2)
        ret = ret * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E
    end function

    subroutine tempToENU_0(phi,E,nu)
        Scalar phi
        Scalar E,nu
        E = 1.0_8
        nu = 0.3_8
    end subroutine

    function expandRate(phi) result(eb)
        Scalar eb, phi
        eb = (phi - 0.0) * 1e-3
    end function

end module

#undef Scalar
