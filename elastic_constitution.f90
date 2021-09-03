#define Scalar real(8)
#define Index integer(8)

module elastic_constitution
    use common_utils
    implicit none
    Scalar E_poly(20), nu_poly(20)
    Index E_poly_siz,nu_poly_siz! order of polynomial is poly_siz-1
    Scalar rho_uni

    Scalar baseTemp, expansionRate 
    Scalar k_ther(3,3), heat_cap_uni

contains

    subroutine set_const_elastic_constitution(E,nu,rho)
        Scalar,intent(in) :: E,nu,rho
        E_poly_siz=1
        E_poly(1:1) = E
        nu_poly_siz=1
        nu_poly(1:1) = nu
        rho_uni = rho
    end subroutine

    subroutine set_const_thremal_constitution(k, c)
        Scalar,intent(in) :: k,c
        k_ther = 0.0_8
        k_ther(1,1) = k
        k_ther(2,2) = k
        k_ther(3,3) = k
        heat_cap_uni = c
    end subroutine

    subroutine set_copper_elastic_constitution
        E_poly_siz=5
        E_poly(1:5) = (/139.059524545458_8, -0.0134613936001288_8, &
                        -0.000120250376148280_8, 1.91953636792926e-07_8 ,-1.20519845938051e-10_8/)
        E_poly(1:5) = E_poly(1:5) * 1e3_8
        nu_poly_siz=1
        nu_poly(1:1) = 0.37_8
        rho_uni = 8900e-12_8
    end subroutine

    subroutine set_expansion_properties(baseTemprature, materialExpasionRate)
        Scalar,intent(in) :: baseTemprature,materialExpasionRate
        baseTemp = baseTemprature
        expansionRate = materialExpasionRate
    end subroutine

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
        E = singlePoly(E_poly(1:E_poly_siz),phi)
        nu = singlePoly(nu_poly(1:nu_poly_siz),phi)
    end subroutine

    function expandRate(phi) result(eb)
        Scalar eb, phi
        eb = (phi - baseTemp) * expansionRate
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
