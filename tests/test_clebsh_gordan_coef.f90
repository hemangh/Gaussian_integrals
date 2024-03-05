!***********************************************************************************************************************************
! clebsch_gordan_test_program
!
! Unit test program for the clebsch_gordan_mod module
!
!***********************************************************************************************************************************

program clebsch_gordan_test_program
    use clebsch_gordan_mod

    implicit none
    integer      :: j1, j2, j3, m1, m2
    real(kind=8) ::  cg_coefficient

    ! Set test case 1
    j1 = 2
    j2 = 1
    j3 = 3
    m1 = 1
    m2 = -1

    ! Calculate Clebsch-Gordan coefficient
    cg_coefficient = calculate_cg_coefficients(j1, j2, j3, m1, m2)

    ! Display result for test case 1
    print *, "Test Case 1:"
    print *, "j1 =", j1, "   j2 =", j2, "   j3 =", j3, "   m1 =", m1, "   m2 =", m2
    print *, "Clebsch-Gordan Coefficient (C) =", cg_coefficient
    print *, "-------------------------"
    if (cg_coefficient /= 0.44721359549995787d0)then
        print*,"case 1 failed"
        print*, cg_coefficient , ' v.s. ', 0.44721359549995787d0 
        stop 1
    end if

    ! Set test case 2
    j1 = 5
    j2 = 4
    j3 = 1
    m1 = 0
    m2 = 0

    ! Calculate Clebsch-Gordan coefficient
    cg_coefficient = calculate_cg_coefficients(j1, j2, j3, m1, m2)

    ! Display result for test case 2
    print *, "Test Case 2:"
    print *, "j1 =", j1, "   j2 =", j2, "   j3 =", j3, "   m1 =", m1, "   m2 =", m2
    print *, "Clebsch-Gordan Coefficient (C) =", cg_coefficient
    print *, "-------------------------"
    if (cg_coefficient /= 0.38924947208076155d0)then
        print*,"case 1 failed"
        print*, cg_coefficient , ' v.s. ', 0.38924947208076155d0 
        stop 1
    end if

    ! Test Case 3
    j1 = 7
    j2 = 3
    j3 = 9
    m1 = 2
    m2 = -1
    cg_coefficient = calculate_cg_coefficients(j1, j2, j3, m1, m2)
    ! Display result for test case 3
    print *, "Test Case 3:"
    print *, "j1 =", j1, "   j2 =", j2, "   j3 =", j3, "   m1 =", m1, "   m2 =", m2
    print *, "Clebsch-Gordan Coefficient (C) =", cg_coefficient
    print *, "-------------------------"
    if (cg_coefficient /= 0.52549385424538797d0)then
        print*,"case 1 failed"
        print*, cg_coefficient , ' v.s. ', 0.52549385424538797d0 
        stop 1
    end if
end program clebsch_gordan_test_program
