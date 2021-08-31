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
        eb = (phi - 0.0) * 1.3e-5
    end function

    function geometryRelation() result(b_dir2strain)
        Scalar b_dir2strain(6,9)
        b_dir2strain = 0.0_8
        b_dir2strain(1,1) = 1.0_8
        b_dir2strain(2,5) = 1.0_8
        b_dir2strain(3,9) = 1.0_8
        b_dir2strain(4,6) = 0.5_8
        b_dir2strain(4,8) = 0.5_8
        b_dir2strain(5,7) = 0.5_8
        b_dir2strain(5,3) = 0.5_8
        b_dir2strain(6,2) = 0.5_8
        b_dir2strain(6,4) = 0.5_8
    end function

    function getUnitBulkStrain() result(unitbulkstrain)
        Scalar unitbulkstrain(6)
        unitbulkstrain = 0.0_8
        unitbulkstrain(1) = 1.0_8
        unitbulkstrain(2) = 1.0_8
        unitbulkstrain(3) = 1.0_8
    end function

    function getVonMises(stress) result(VM)
        Scalar VM, stress(6)
        VM = sqrt((&
        (stress(1)-stress(2))**2+(stress(2)-stress(3))**2+(stress(3)-stress(1))**2+&
        6*(stress(4)**2+stress(5)**2+stress(6)**2)&
        )/2.0_8)

    end function

end module

#undef Scalar
