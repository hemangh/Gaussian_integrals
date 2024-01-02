module Gaussian_expansion_coefficient_Ylm
use precision
use factorials_mod, only: binomial_coeff
use conversions, only: cart2sph
    implicit none

    contains
    !
    ! (x-Rx)^l * (y-Ry)^m * (z-Rz)^n = \sum_{ijk} A^{lmn}_{ijk} r^{i+j+k}
    !
    function cartesian_to_spherical_prefactor_coeff (x, Rx, l, y,Ry,m, z, Rz, n) result(sum_ijk)
        real(idp):: x, y, z
        real(idp):: Rx,Ry, Rz
        integer  :: l, m, n
        real(idp) :: sum_ijk
        

        integer :: i, j, k, r_pow
        real(idp) :: a, b, c
        real(idp) :: term
        real(idp), allocatable :: coeff(:,:)
        real(idp), allocatable :: coeff_vector(:)

        real(idp):: rthetaphi(3)
        real(idp) :: r, theta, phi

        allocate(coeff(0:max(l,m,n),3))
        coeff = 0._idp

        allocate(coeff_vector(0:sum([l,m,n])))
        coeff = 0._idp
        ! convert x, y, z to r, thea, phi (e.g., spherical coordinate)
        rthetaphi = cart2sph([x,y,z])
        r = rthetaphi(1)
        theta = rthetaphi(2)
        phi = rthetaphi(3)
        a = sin(theta) * cos(phi)
        b = sin(theta) * sin(phi)
        c = cos(theta)

        !(x-Rx)^l => \sum_{k=0}^l binomial_coeff(l,k) * (a*r)^{l-k} * Rx^k 
        do k = 0, l
            coeff(k,1) = binomial_coeff(l,k) * Rx**k * a**(l-k)
        end do

        !(y-Ry)^m => \sum_{k=0}^m binomial_coeff(m,k) * (b*r)^{m-k} * Ry^k 
        do k = 0, m
            coeff(k,2) = binomial_coeff(m,k) * Ry**k * b**(m-k)
        end do

        !(z-Rz)^n => \sum_{k=0}^n binomial_coeff(n,k) * (c*r)^{n-k} * Rz^k 
        do k = 0, n
            coeff(k,3) = binomial_coeff(n,k) * Rz**k * c**(n-k)
        end do

        ! now multiply each term from x, y and z into every other term:
        do i = 0, l
            do j = 0, m
                do k = 0, n
                    r_pow = i+j+k
                    term = coeff(i,1) * coeff(j,2) * coeff(k, 3) * r**r_pow
                    coeff_vector(r_pow) = coeff_vector(r_pow) + term
                end do
            end do
        end do

        sum_ijk = sum(coeff_vector)

    end function


end module