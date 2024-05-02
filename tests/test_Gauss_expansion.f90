program test_Gauss_expansion
    use precision, only: idp
    use constants, only : pi
    use special_functions, only: ikna
    use sphericalHarmonics, only: real_spherical_harmonics
    use Gaussian_expansion_coefficient_Ylm, only: expansion_primitiveG_Ylm
    implicit none

    integer :: l
    integer :: m
    integer :: nm
    real(idp) :: a ,r, ra, x
    real(idp) :: mod_bessel
    real(idp) , dimension(0:1) :: bi, di, bk, dk
    real(idp), allocatable :: Ylm(:,:,:)
    real(idp) :: C_10
    real(idp), allocatable :: coeff(:,:) 
    real(idp)  :: grid_point(1,3)

    ! Test Case: 1
    ! Test against the known expression for s function: G_{000} (r)
    ! C^{lm}_{000a} = 4 pi (2a/pi) ^ {3/4} exp(-a(r-R)^2) exp(-2arR) * i_l(2arR) * Y_lm(R)

    ! take a fix r, a, R, lm:
    l  = 1
    m  = 0
    ra = 5.d0
    r  = 1.d0
    a  = 1.d0

    x = 2.d0 * r * ra * a

    

    call ikna(l, x, nm, bi, di, bk, dk)
    mod_bessel = bi(l) * exp(-x)
    
    print*, "modified bessel:", bi(l), mod_bessel

    call real_spherical_harmonics(Ylm,[0.d0],[0.d0],1,1)
 
    
    print*, "real sphY_l=1,m=0:", Ylm(1,1,0)

    C_10 = 4.d0 * pi * (2.d0 * a/pi) ** (3/4) * exp(-a*(r-R)**2) * mod_bessel * Ylm(1,1,0)

    print*, "analytical C_10(r=1.0):", C_10 

    grid_point = reshape([r,0.d0,0.d0],(/1,3/))
    print*, grid_point(1,1), grid_point(1,2), grid_point(1,3)
    call expansion_primitiveG_Ylm(l, grid_point,[0,0,0],[0.d0,0.d0,ra],a,coeff)
    print*, "my C_10(r=1.0):", coeff(1,0) , coeff(1,1)
    
    
end program test_Gauss_expansion