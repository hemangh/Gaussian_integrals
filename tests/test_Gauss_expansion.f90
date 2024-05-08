program test_Gauss_expansion
    use precision, only: idp
    use constants, only : pi
    use special_functions, only: ikna
    use sphericalHarmonics, only: real_spherical_harmonics
    use Gaussian_expansion_coefficient_Ylm
    implicit none

    integer :: l
    integer :: m
    integer :: nm
    real(idp) :: a ,r, ra, x
    real(idp) :: mod_bessel
    real(idp) , dimension(0:1) :: bi, di, bk, dk
    real(idp), allocatable :: Ylm(:,:,:)
    real(idp) :: C_lm
    real(idp), allocatable :: coeff(:,:,:) 
    real(idp)  :: grid_point(1,3), RAsph_coord(3), Rxyz_coord(3)

    ! Test Case: 1
    ! Test against the known expression for s function: G_{000} (r)
    ! C^{lm}_{000a} = 4 pi (2a/pi) ^ {3/4} exp(-a(r-R)^2) exp(-2arR) * i_l(2arR) * Y_lm(R)

    ! take a fix r, a, R, lm:
    l  = 1
    !m  = 0
    
    a  = 1.d0

    grid_point = reshape([0.d0,0.d0,r],(/1,3/))
    print*, grid_point(1,1), grid_point(1,2), grid_point(1,3)
    r = sqrt(dot_product(grid_point(1,:),grid_point(1,:)))

    Rxyz_coord = [0.d0,sqrt(.5d0),sqrt(.5d0)]

    RAsph_coord = cart2sph(Rxyz_coord)
    ra = RAsph_coord(1)

    x = 2.d0 * r * ra * a

    call ikna(l, x, nm, bi, di, bk, dk)
    mod_bessel = bi(l) * exp(-x)
    
    print*, "modified bessel:", bi(l), mod_bessel

    call real_spherical_harmonics(Ylm,RAsph_coord(2),RAsph_coord(3),1,1)
    

    call expansion_primitiveG_Ylm(l, grid_point,[0,0,0],Rxyz_coord,a,coeff)

    do m = -l , l

        print*, "m:", m
    
        print'("real sphY_l=",I0,",m=",I0,":",E15.8,", B_",I0,":",E15.8)', l,m,Ylm(1,l,m), l , mod_bessel

        C_lm = 4.d0 * pi * (2.d0 * a/pi) ** (3/4) * exp(-a*(r-ra)**2) * mod_bessel * Ylm(1,l,m)

        print'( "analytical C_",I0,",",I0,"(r=1.0):",E15.8)',l,m, C_lm 

        
        print*, "norm:", normalization(a,[0,0,0])
        print*, "sum:", mod_bessel * Ylm(1,l,m), "exp:", exp(-a*(r-ra)**2) 

        print'("my C_",I0,",",I0,"(r=1.0):", E15.8)',l,m, coeff(1,l,m) 

        if(coeff(1,l,m) /= C_lm) print*, "Error:", coeff(1,l,m), "v.s" , C_lm, 2.d0*coeff(1,l,m)
        

    end do
    
    
    
    
    
end program test_Gauss_expansion