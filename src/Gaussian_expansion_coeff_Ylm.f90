module Gaussian_expansion_coefficient_Ylm
use precision, only: idp
use factorials_mod, only: binomial_coeff
use conversions, only: cart2sph
use special_functions, only: ikna
use sphericalHarmonics, only: real_spherical_harmonics
    implicit none
    integer, parameter :: max_l = 4
    ! Coefficient of x^i y^j z^k = B^{LM}_{ijk} r^{i+j+k} Y_{LM}   
    ! (L, M, i, j, k)                                    
    real(idp) :: xyz_YLM_coefficient(0:max_l, -max_l:max_l, 0:max_l, 0:max_l, 0:max_l) 

    ! xyz_YLM_coefficient(0, 0, 0, 0, 0) = 0.28209479177387814d0 ! Y_00

    ! xyz_YLM_coefficient(1, -1, 0, 1, 0) = 0.4886025119029199d0  ! Y_1m1
    ! xyz_YLM_coefficient(1, 0, 0, 0, 0)  = 0.4886025119029199d0  ! Y_10
    ! xyz_YLM_coefficient(1, 1, 1, 0, 0)  = 0.4886025119029199d0  ! Y_11

    ! xyz_YLM_coefficient(2, -2, 1, 1, 0) = 0.31539156525252005d0 ! Y_2m2
    ! xyz_YLM_coefficient(2, -1, 0, 1, 0) = 0.5462742152960396d0  ! Y_2m1
    ! xyz_YLM_coefficient(2, 0, 1, 1, 0)  = 0.5462742152960396d0   ! Y_20
    ! xyz_YLM_coefficient(2, 1, 1, 1, 0)  = 0.5462742152960396d0   ! Y_21
    ! xyz_YLM_coefficient(2, 2, 1, 1, 0)  = 0.5462742152960396d0   ! Y_22
  
    ! xyz_YLM_coefficient(3, -3, 0, 3, 0) = 0.13895681407382172d0 ! Y_3m3
    ! xyz_YLM_coefficient(3, -2, 1, 1, 0) = 0.5900435899266439d0  ! Y_3m2
    ! xyz_YLM_coefficient(3, -1, 0, 1, 0) = 0.28449475869314084d0 ! Y_3m1
    ! xyz_YLM_coefficient(3, 0, 1, 1, 0) = 0.5900435899266439d0   ! Y_30
    ! xyz_YLM_coefficient(3, 1, 1, 1, 0) = 0.28449475869314084d0  ! Y_31
    ! xyz_YLM_coefficient(3, 2, 1, 1, 0) = 0.5900435899266439d0   ! Y_32
    ! xyz_YLM_coefficient(3, 3, 0, 3, 0) = 0.13895681407382172d0 ! Y_33
  
    ! xyz_YLM_coefficient(4, -4, 3, 1, 0) = 0.4886025119029199d0  ! Y_4m4
    ! xyz_YLM_coefficient(4, -3, 0, 3, 0) = 0.46557123656573955d0 ! Y_4m3
    ! xyz_YLM_coefficient(4, -2, 1, 1, 0) = 0.12593510574971699d0 ! Y_4m2
    ! xyz_YLM_coefficient(4, -1, 0, 1, 0) = 0.12593510574971699d0 ! Y_4m1
    ! xyz_YLM_coefficient(4, 0, 4, 0, 0) = 0.119514472455632d0   ! Y_40
    ! xyz_YLM_coefficient(4, 1, 0, 3, 0) = 0.12593510574971699d0 ! Y_41
    ! xyz_YLM_coefficient(4, 2, 2, 2, 0) = 0.059857729514531774d0 ! Y_42
    ! xyz_YLM_coefficient(4, 3, 0, 3, 0) = 0.12593510574971699d0 ! Y_43
    ! xyz_YLM_coefficient(4, 4, 4, 0, 0) = 0.02370600766752924d0 ! Y_44

    contains
    !***********************************************************************************************
    ! (x-Rx)^l * (y-Ry)^m * (z-Rz)^n = \sum_{ijk} A^{lmn}_{ijk} x^{i}*y^{j}*z^{k}
    !
    ! Input: 
    !        Cartesian coordinate of primitive Gaussian's center: Rxyz_coord => [Rx, Ry, Rz] (double)
    !        power of x,y,z terms of the Gaussian primitive, in the non-local (off-center) coordinate
    !                               :  pwr(x,y,z) => [l, m, n] (integer)
    ! Output: matrix coefficients (A^{lmn}_{ijk}) where i, j, k are powers of x, y, z in the local coordinate
    !                               : coeff(0:max(l,m,n),3)  (double)
    !************************************************************************************************
    subroutine localcartesian_coeff ( Rxyz_coord, pwr, coeff)

       
        real(idp), intent(in) :: Rxyz_coord(3) 
        real(idp), intent(in) :: pwr(3)

        real(idp), allocatable, intent(out) :: coeff(:,:) !result

        real(idp):: Rx,Ry, Rz
        integer  :: l, m, n
        
        integer :: k

        l = pwr(1)
        m = pwr(2)
        n = pwr(3)
        allocate(coeff(0:max(l,m,n),3))
        coeff = 0._idp

        !(x-Rx)^l => \sum_{k=0}^l binomial_coeff(l,k) * (x)^{l-k} * Rx^k 
        do k = 0, l
            coeff(l-k,1) = binomial_coeff(l,k) * Rxyz_coord(1)**k !* x**(l-k)
        end do

        !(y-Ry)^m => \sum_{k=0}^m binomial_coeff(m,k) * (y)^{m-k} * Ry^k 
        do k = 0, m
            coeff(m-k,2) = binomial_coeff(m,k) * Rxyz_coord(2)**k !* y**(m-k)
        end do

        !(z-Rz)^n => \sum_{k=0}^n binomial_coeff(n,k) * (z)^{n-k} * Rz^k 
        do k = 0, n
            coeff(n-k,3) = binomial_coeff(n,k) * Rxyz_coord(3)**k !* z**(n-k)
        end do

    end subroutine localcartesian_coeff
    ! 
    ! \sum_{l,m} exp(-2*a*RA*r)* i_l(2*a*RA*r) * Y_{lm}(\hat{RA})
    ! 
    function sum_over_Bessel_i_Ylm(a, RA, ThetaA,PhiA, r, lmax) result(sum2)
        real(idp) :: a, RA, ThetaA,PhiA, r
        integer   :: lmax
        real(idp) :: sum2

        integer :: nm
        real(idp) , dimension(0:lmax) :: bi, di, bk, dk
        real(idp), allocatable :: Ylm(:,:,:), theta(:), phi(:)
        real(idp) :: x
        integer :: l, m

        sum2 = 0._idp
        x = 2._idp * a * RA * r
        allocate(theta(1), phi(1))
        theta(1) = ThetaA
        phi(1)   = PhiA
        call ikna(lmax, x, nm, bi, di, bk, dk)
        call real_spherical_harmonics(Ylm,theta,phi,1,lmax) !allocated according to phi/theta dim & lmax in the routine
        do l = 0, lmax
            do m = -l, l
                sum2 = sum2 + exp(-x) * bi(l) * Ylm(1,l,m) ! * B^{lm}_{ijk} * Ylm(r) * 
            end do
        end do

        end function sum_over_Bessel_i_Ylm

        !
        ! \sum_{l3, m3} \int Y_{l1,m1} Y_{l2, m2} Y_{l3, m3}
        !
        ! function sum_over_Clebsh_Gordon_Constants(l1, l2, l3, m1, m2, m3) result(sum3)
            
            ! integer, intent(in) :: l, m1, l2, m2, l3, m3
            ! real(idp) :: sum3
            
            ! end function sum_over_Clebsh_Gordon_Constants

end module