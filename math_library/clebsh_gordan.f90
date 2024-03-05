!***********************************************************************************************************************************
! clebsch_gordan_module
!
! Module providing Clebsch-Gordan coefficients
!
! Original Author: David G. Simpson, NASA Goddard Space Flight Center, Greenbelt, Maryland, 20771
! Date: April 20, 2005
!
!***********************************************************************************************************************************

module clebsch_gordan_mod
    use factorials_mod, only: factorial, compute_factorials

    implicit none
  
    private
  
  
    public :: calculate_cg_coefficients
  
    contains
  

!***********************************************************************************************************************************
  ! calculate_cg_coefficients
  !
  ! Subroutine to calculate Clebsch-Gordan coefficients.
  !
  ! Assumptions:
  !   - Assumes factorial is a precalculated array of factorials up to 100
  !
  !***********************************************************************************************************************************

    function calculate_cg_coefficients( j1,  j2,  j3, m1, m2 ) result(cg)
        integer, intent(in) ::  j1,  j2,  j3, m1, m2
        real(kind=8) :: cg
        integer :: k, m3
        real(kind=8) :: sum_k, term

        m3 =   m1 +   m2

    !
    !     Check for conditions that give cg = 0.
    !

    if ( ( j3 .lt. abs( j1- j2)) .or.  &
        ( j3 .gt. ( j1+ j2))    .or.  &
        (ABS(M1) .gt.  j1)    .or.  &
        (ABS(M2) .gt.  j2)    .or.  &
        (ABS(M3) .gt.  j3)) then
            cg = 0.0d0
        else
        ! (Factorial function calls need to be replaced with the actual values)
        call compute_factorials

        cg = sqrt((j3+j3 + 1) / factorial( j1 + j2 + j3 + 1))
        cg =   cg * sqrt(factorial(j1 + j2 -j3) * factorial(j2 + j3 - j1) &
                          * factorial( j3 + j1 -j2))
        cg =   cg * sqrt(factorial(j1 + m1) * factorial(j1 - m1) * factorial(j2 + m2) &
                        * factorial(j2 -m2) * factorial(j3 + m3) * factorial(j3 - m3))
      sum_k = 0.0d0
      do k = 0, 99
        if (   j1 +    j2 -    j3 - k .lt. 0.0d0) cycle
        if (   j3 -    j1 -   m2 + k .lt. 0.0d0) cycle
        if (   j3 -    j2 +   m1 + k .lt. 0.0d0) cycle
        if (   j1 -   m1 - k .lt. 0.0d0) cycle
        if (   j2 +   m2 - k .lt. 0.0d0) cycle
        term = factorial( (   j1 +    j2 -    j3 - k)) * factorial( (   j3 -    j1 -   m2 + k)) &
             * factorial( (   j3 -    j2 +   m1 + k)) * factorial( (   j1 -   m1 - k)) &
             * factorial( (   j2 +   m2 - k)) * factorial(k)
        if (mod(k, 2) .eq. 1) term = -term
        sum_k = sum_k + 1.0d0 / term
      end do
        cg =   cg * sum_k
    end if
    end function calculate_cg_coefficients
  
  end module clebsch_gordan_mod
  