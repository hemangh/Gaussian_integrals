program test_Gauss_expansion
    use precision, only: idp
    use constants, only : pi
    use special_functions, only: ikna, sphi
    use sphericalHarmonics, only: real_spherical_harmonics
    use lebedev_quadrature, only: Lebedev_type
    use Gaussian_expansion_coefficient_Ylm
    implicit none

    integer :: l, il, i
    integer :: m
    integer :: nm
    integer :: rule
    integer , parameter :: num_pt = 10
    
    real(idp) :: a ,r, ra, x
    real(idp) :: mod_bessel
    real(idp) , allocatable, dimension(:) :: bi, di, bk, dk
    real(idp), allocatable :: Ylm(:,:,:)
    real(idp) :: C_lm
    real(idp) :: rand_fun
    real(idp), allocatable :: coeff(:,:,:) , rand_grid_pt(:,:), r2(:)
    real(idp), allocatable :: pxgrid(:), sgrid(:)
    real(idp)  :: grid_point(num_pt,3), RAsph_coord(3), Rxyz_coord(3), rand_sphgridpt(3)
    real(idp)  :: sum_g, scale
    real(idp)  :: Y20, Y11, Y21, Y10, Y22
    
    type(Lebedev_type) :: leb_quad

    integer, parameter :: lmax = 4

    REAL(kind = idp ), DIMENSION(0:2*lmax) :: si
    REAL(kind = idp ), DIMENSION(0:2*lmax) :: dsi

    ! Test Case: 1
    ! Test against the known expression for s function: G_{000} (r)
    ! C^{lm}_{000a} = 4 pi (2a/pi) ^ {3/4} exp(-a(r-R)^2) exp(-2arR) * i_l(2arR) * Y_lm(R)

    ! take a fix r, a, R, lm:
    
    a  = 0.5d0
 ! scale the grid points:
    scale = 1.0_idp

    allocate(bi(0:2*lmax), di(0:2*lmax),bk(0:2*lmax),dk(0:2*lmax))

    call generate_random_points(num_pt, grid_point)

    !print*, grid_point(1,1), grid_point(1,2), grid_point(1,3)

    

    ! Rxyz_coord = [0.d0,sqrt(.5d0),sqrt(.5d0)]
    Rxyz_coord = [0.25d0,sqrt(.5d0),0.d0]

    RAsph_coord = cart2sph(Rxyz_coord)
    print'("Spherical Coordinates of R_A: r_A=",F5.1, ", theta_A=", F8.3, " ,Phi_A=", F8.3 )', RAsph_coord(:)
    print*,RAsph_coord(1)
    ra = RAsph_coord(1)

    call real_spherical_harmonics(Ylm,[RAsph_coord(2)],[RAsph_coord(3)],1,lmax)
    do il = 0, lmax
        do m=-il,il
    print*, ylm(:,il,m), il, m
        end do
    end do
    Y11 = Ylm(1,1,1)
    Y10 = Ylm(1,1,1)
    Y20 = Ylm(1,2,0)
    Y21 = Ylm(1,2,1)
    Y22 = Ylm(1,2,2)
    
!     ! S-func expansion:
!     call expansion_primitiveG_Ylm(lmax, grid_point,[0,0,0],Rxyz_coord,a,coeff)

!     do i = 1, num_pt

!         r = sqrt(dot_product(grid_point(i,:),grid_point(i,:)))

!         x = 2.d0 * r * ra * a

!         call ikna(4, x, nm, bi, di, bk, dk)

!         call sphi(4,x,nm,si,dsi)

!         ! print*, bi(0), si(0), sinh(x)/x 
!         ! print*, bi(1), si(1), (x*cosh(x)-sinh(x))/x**2 
!         ! print*, bi(2), si(2), ((x**2+3._idp)*sinh(x)-3._idp*x*sinh(x))/x**3
!         ! print*, bi(3), si(3), ((x**3+15._idp*x)*cosh(x)-(6._idp*x**2+15_idp)*sinh(x))/x**4
!         ! print*, bi(4), si(4), ((x**4+45._idp*x**2+105._idp)*sinh(x)-(10._idp*x**3+105_idp*x)*cosh(x))/x**5
        

!         do l = 0,lmax

!         mod_bessel = si(l) * exp(-x)
    
!         print*, "modified bessel:", si(l), mod_bessel

        

!         do m = -l , l
    

!         print*, "m:", m
    
!        print'("real sphY_l=",I0,",m=",I0,":",E15.8,", B_",I0,":",E15.8)', l,m,Ylm(1,l,m), l , mod_bessel

!         C_lm = 4.d0 * pi * (2.d0 * a/pi) ** (3/4) * exp(-a*(r-ra)**2) * mod_bessel * Ylm(1,l,m)

!        print'( "analytical C_",I0,",",I0,"(r=",F3.1,"):", E15.8)',l,m,r, C_lm 

        
!         print*, "norm:", 4.d0 * pi*normalization(a,[0,0,0]),4.d0 * pi * (2.d0 * a/pi) ** (3/4)
!         !print*, "sum:", mod_bessel * Ylm(1,l,m), "exp:", exp(-a*(r-ra)**2) 

       

!         print'("my C_",I0,",",I0,"(r=",F3.1,"):", E15.8)',l,m,r, coeff(i,l,m) 

!          if(coeff(i,l,m) /= C_lm) print*, "Error:", coeff(i,l,m), "v.s" , C_lm ,'diff:',abs(C_lm-coeff(i,l,m)) ! 2.d0*coeff(1,l,m)
        
        
!     end do !m

! end do !l

! end do !point



!     ! ! make random grid points:
!     ! allocate(rand_grid_pt(num_pt, 3))

!     ! call generate_random_points(num_pt, rand_grid_pt)

    ! make a lebedev grid of just order 5:
    rule = 65
    call leb_quad%get_rule_points(rule)
    print*, 'precision of quadrature:',leb_quad%precision, "number of grid points:", leb_quad%num_points
    if (allocated(rand_grid_pt)) deallocate(rand_grid_pt)
    if (allocated(sgrid)) deallocate(sgrid)
    if (allocated(pxgrid)) deallocate(pxgrid)

    allocate(rand_grid_pt, source=leb_quad%cart_coord)
    allocate(sgrid(leb_quad%num_points))
    allocate(pxgrid(leb_quad%num_points))

   
    do i = 1, leb_quad%num_points
    
     rand_grid_pt(i,:) = scale * leb_quad%cart_coord(i,:)

    end do

    call real_spherical_harmonics(Ylm,leb_quad%theta,leb_quad%phi,leb_quad%num_points,lmax)

!     ! Check s function:
!     call expansion_primitiveG_Ylm(lmax =lmax, grid_points =rand_grid_pt,powers =[0,0,0],&
!                                      Rxyz_coord = Rxyz_coord, alpha = a,coeff = coeff)

!     sgrid = s(leb_quad%num_points,rand_grid_pt, Rxyz_coord,a)
    
!     print*, "s:"
!     do il = 0, lmax
!         do m = -il, il

!             sum_g = sum(4._idp*pi*sgrid(:) * Ylm(:,il,m)*leb_quad%w(:))

!             ! check: are all clm(r) the same?
!             print'("C_",I0,",",I0,"(r=",F3.2,"):", 3E15.8)',il,m,scale, coeff(1,il,m), &
!             coeff(leb_quad%num_points/2,il,m), coeff(leb_quad%num_points,il,m)

!             print'("Numeric: C_",I0,",",I0,"(r=",F3.2,"):", 3E15.8)',il,m,scale,  sum_g
!             if(abs(sum_g)>=1.e-12_idp)print*,"ratio:", coeff(1,il,m)/sum_g

! !         sum_g = sum_g + coeff(i,il,m) * Ylm(1,il,m) 

!         end do !l
!     end do !m
    
    
    
    ! Choose px function to expand:

    ! call expansion_primitiveG_Ylm(lmax =lmax, grid_points =rand_grid_pt,powers =[1,0,0],&
    !                                 Rxyz_coord = Rxyz_coord, alpha = a,coeff = coeff)
    
    ! pxgrid = px(leb_quad%num_points,rand_grid_pt, Rxyz_coord,a)
!     allocate(r2(leb_quad%num_points))
!     r2(:) = &
!    (leb_quad%cart_coord(:,1) - Rxyz_coord(1)) **2 + (leb_quad%cart_coord(:,2) - Rxyz_coord(2)) **2 + &
!    (leb_quad%cart_coord(:,3) - Rxyz_coord(3)) **2
    
!    print*,"Explicit sum on grid: (C_1,1:)", &
!    4._idp*pi*normalization(a,[1,0,0])*&
!    sum((leb_quad%cart_coord(:,1) - Rxyz_coord(1)) * exp(-a*r2(:)) * Ylm(:,1,1) * leb_quad%w(:))
    

    ! ! do i = 1, num_pt
    ! i = 1                                 
    !     rand_fun = px(rand_grid_pt(i:i,:), Rxyz_coord,a)

    !     rand_sphgridpt = cart2sph(rand_grid_pt(i,3))

    !     call real_spherical_harmonics(Ylm,rand_sphgridpt(2),rand_sphgridpt(3),1,l)

    !     sum_g = 0._idp

    ! print*,"px:"

    !     do il = 0, lmax
    !         do m = -il, il

   
    !             ! print*, sum(4.d0*pi* Ylm(:,il,m) * Ylm(:,il,m)*leb_quad%w(:))

    !             sum_g = sum(4.d0*pi*pxgrid(:) * Ylm(:,il,m)*leb_quad%w(:))

    !             ! check: are all clm(r) the same?
    !             print'("C_",I0,",",I0,"(r=",F3.1,"):", 1E15.8)',il,m,scale, coeff(1,il,m)!, &
    !             !coeff(leb_quad%num_points/2,il,m), coeff(leb_quad%num_points,il,m)

    !             print'("Numeric: C_",I0,",",I0,"(r=",F3.1,"):", 3E15.8)',il,m,scale, sum_g
    !             if(abs(sum_g)>=1.e-12_idp)print*,"ratio:", coeff(1,il,m)/sum_g

    ! !         sum_g = sum_g + coeff(i,il,m) * Ylm(1,il,m) 

    !         end do !l
    !     end do !m


        print*,"dxy:"

        ! Choose dxy function to expand:

        call expansion_primitiveG_Ylm(lmax =lmax, grid_points =rand_grid_pt,powers =[1,1,0],&
        Rxyz_coord = Rxyz_coord, alpha = a,coeff = coeff)

        pxgrid = dxy(leb_quad%num_points,rand_grid_pt, Rxyz_coord,a)

        do il = 0, lmax
            do m = -il, il

                ! il = 0 ; m =0
                ! print*, sum(4.d0*pi* Ylm(:,il,m) * Ylm(:,il,m)*leb_quad%w(:))

                sum_g = sum(4.d0*pi*pxgrid(:) * Ylm(:,il,m)*leb_quad%w(:))

                ! check: are all clm(r) the same?
                print'("C_",I0,",",I0,"(r=",F3.1,"):", 1E15.8)',il,m,scale, coeff(1,il,m)!, &
                !coeff(leb_quad%num_points/2,il,m), coeff(leb_quad%num_points,il,m)

                print'("Numeric: C_",I0,",",I0,"(r=",F3.1,"):", 3E15.8)',il,m,scale, sum_g
                if(abs(sum_g)>=1.e-12_idp)print*,"ratio:", coeff(1,il,m)/sum_g

    !         sum_g = sum_g + coeff(i,il,m) * Ylm(1,il,m) 

            end do !l
        end do !m

        ! print*, "B2-2_110 (Analytic):", 2.0_idp * sqrt(pi/15.0_idp), "B2-2_110 (Calculated):",xyz_YLM_coefficient(2,-2,1,1,0)

    !     print*, sum_g, rand_fun

    ! end do !i

        ! print*, "B11_100 (Analytic):", sqrt(4._idp * pi/3._idp), "B11_100 (Calculated):",xyz_YLM_coefficient(1,1,1,0,0)

        ! x = 2.d0 * scale * ra * a
        ! call sphi(lmax,x,nm,si,dsi)

        ! ! print*, si(0), si(1), si(2)
        ! ! print*, Y11, Rxyz_coord(1); stop

        ! mod_bessel = si(0) *  exp(-x)

        ! C_lm = 4.d0 * pi * normalization(a,[1,0,0]) * exp(-a*(scale-ra)**2) * scale* &
        ! exp(-x) * (si(0)  * 0.25_idp/pi + &
        ! si(1) * Y10  * integrals_3rspH(1,1,1,1,0,1) + &
        ! si(2) * Y20  * integrals_3rspH(1,2,1,1,0,1) + &
        ! si(2) * Y22  * integrals_3rspH(1,2,1,1,2,1))&
        ! * xyz_YLM_coefficient(1,1,1,0,0) 
        
        
        ! print*,"explicitly calculated C_1,1 (using fortran functions):", C_lm
        ! print*,"explicitly calculated C_1,1 (using mathematica numbers):", &
        ! 4.d0 * pi * (2._idp*a/pi)**(3/4) * sqrt(4._idp*a) * exp(-a*(scale*scale+ra*ra)) *( scale* &
        ! (si(0) * 0.25_idp/pi   -si(2) * Y20  *0.12615662610100803_idp &
        ! +si(2)*0.21850968611841617* Y22) * sqrt(4._idp*pi/3._idp) &
        ! - Rxyz_coord(1)* si(1)*Y11*xyz_YLM_coefficient(0,0,0,0,0)*0.5_idp/sqrt(pi))

        ! ! print*, "point by point comparison"


        ! ! do i = 1, leb_quad%num_points

        ! !     print*,pxgrid(i)* Ylm(i,il,m), C_lm , pxgrid(i)/ C_lm 

        ! ! end do


        ! ! print*,sum_over_Clebsh_Gordon_Constants(1,0,1,1,0,1)
        ! ! print*, sum_over_Clebsh_Gordon_Constants(1,1,1,0,1,1)
        ! ! print*, sum_over_Clebsh_Gordon_Constants(2,1,1,0,1,1), sum_over_Clebsh_Gordon_Constants(1,2,1,1,0,1)
        ! ! print*, calculate_cg_coefficients(1,1,1,0,0), calculate_cg_coefficients(1,1,1,0,1)

        ! print*, "comparison of numerical integral of 3 Ylms"
        ! print*, "==========================================="
        ! print*, "      CG sums", "               ",   "Lebedev sums"
    
        ! print*, sum_over_Clebsh_Gordon_Constants(1,0,1,1,0,1),&
        ! sum(4.d0*pi* Ylm(:,0,0) * Ylm(:,1,1)*Ylm(:,1,1)*leb_quad%w(:)), "<l=1,m=1|l'=0,m'=0|L=1,M=1>"
        ! print*, sum_over_Clebsh_Gordon_Constants(1,1,1,0,1,1),&
        ! sum(4.d0*pi* Ylm(:,1,0) * Ylm(:,1,1)*Ylm(:,1,1)*leb_quad%w(:)), "<l=1,m=1|l'=1,m'=0|L=1,M=1>"
        ! print*, sum_over_Clebsh_Gordon_Constants(2,1,1,0,1,1), &
        ! sum(4.d0*pi* Ylm(:,2,0) * Ylm(:,1,1)*Ylm(:,1,1)*leb_quad%w(:)), "<l=1,m=1|l'=2,m'=0|L=1,M=1>"
        ! print*, sum_over_Clebsh_Gordon_Constants(2,1,1,2,-1,1), &
        ! sum(4.d0*pi* Ylm(:,2,2) * Ylm(:,1,1)*Ylm(:,1,-1)*leb_quad%w(:)), "<l=1,m=1|l'=2,m'=2|L=1,M=-1>"
        ! print*, sum_over_Clebsh_Gordon_Constants(2,1,1,1,0,1), &
        ! sum(4.d0*pi* Ylm(:,2,1) * Ylm(:,1,1)*Ylm(:,1,0)*leb_quad%w(:)), "<l=1,m=1|l'=2,m'=1|L=1,M=0>"
        ! print*, sum_over_Clebsh_Gordon_Constants(2,1,1,-2,1,1), &
        ! sum(4.d0*pi* Ylm(:,2,-2) * Ylm(:,1,1)*Ylm(:,1,1)*leb_quad%w(:)), "<l=1,m=1|l'=2,m'=-2|L=1,M=1>"
        ! print*, sum_over_Clebsh_Gordon_Constants(2,1,1,2,1,1), &
        ! sum(4.d0*pi* Ylm(:,2,2) * Ylm(:,1,1)*Ylm(:,1,1)*leb_quad%w(:)), "<l=1,m=1|l'=2,m'=2|L=1,M=1>"

        ! print*,calculate_cg_coefficients(1,1,2,0,0), calculate_cg_coefficients(1,1,2,1,1), &
        ! sqrt((2*1+1) * (2*1+1)/4._idp/ pi/ (2*2+1) ) * calculate_cg_coefficients(1,1,2,0,0)&
        ! * calculate_cg_coefficients(1,1,2,1,1)* pi /4.d0
                

    contains 

    subroutine generate_random_points(num_pt, rand_gridpt)
        integer, intent(in) :: num_pt
        real(idp), intent(out) :: rand_gridpt(num_pt, 3)
        integer :: i, seed_size
        real(idp) :: seed_value
        integer, allocatable :: seed(:)
        
        ! Set the seed for reproducibility
        seed_value = 12345.0_idp
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed = floor(seed_value * (/(i, i=1, seed_size)/))
    
        call random_seed(put=seed)
    
        ! Generate random points
        do i = 1, num_pt
          call random_number(rand_gridpt(i, 1))
          rand_gridpt(i, 1) = rand_gridpt(i, 1) * 2.0_idp
          call random_number(rand_gridpt(i, 2))
          rand_gridpt(i, 2) = rand_gridpt(i, 2) * 2.0_idp
          call random_number(rand_gridpt(i, 3))
          rand_gridpt(i, 3) = rand_gridpt(i, 3) * 2.0_idp
        end do
      end subroutine generate_random_points


    function s(num_pt,cart_coord, center_coord, alpha) result(sval)

        integer   :: num_pt
        real(idp) :: cart_coord(num_pt,3), center_coord(3)
        real(idp) :: alpha
        real(idp) :: sval(num_pt)

        real(idp) :: new_coord(num_pt,3)
        real(idp) :: r2(num_pt)
        
        do i = 1, num_pt
            new_coord(i,1:3) = cart_coord(i,1:3) - center_coord(1:3)
            r2(i) = dot_product(new_coord(i,:),new_coord(i,:))
        end do

        sval(:) = normalization(alpha,[0,0,0]) * exp(-alpha*r2(:))

    end function

    function px(num_pt,cart_coord, center_coord, alpha) result(pxval)

        integer   :: num_pt
        real(idp) :: cart_coord(num_pt,3), center_coord(3)
        real(idp) :: alpha
        real(idp) :: pxval(num_pt)

        real(idp) :: new_coord(num_pt,3)
        real(idp) :: r2(num_pt)
        real(idp) :: norm

        norm =  normalization(alpha,[1,0,0])
        
        do i = 1, num_pt
            new_coord(i,1) = cart_coord(i,1) - center_coord(1)
            new_coord(i,2) = cart_coord(i,2) - center_coord(2)
            new_coord(i,3) = cart_coord(i,3) - center_coord(3)
            ! r2(i) = dot_product(new_coord(i,:),new_coord(i,:))
            r2(i) = new_coord(i,1) **2 + new_coord(i,2) **2 + new_coord(i,3) **2
            pxval(i) = norm*(new_coord(i,1)) * exp(-alpha*r2(i))
        end do

        

        ! pxval(:) = normalization(alpha,[1,0,0])*(new_coord(:,1)) * exp(-alpha*r2(:))

    end function

    function py(num_pt,cart_coord, center_coord, alpha) result(pyval)

        integer   :: num_pt
        real(idp) :: cart_coord(num_pt,3), center_coord(3)
        real(idp) :: alpha
        real(idp) :: pyval(num_pt)

        real(idp) :: new_coord(num_pt,3)
        real(idp) :: r2(num_pt)
        
        do i = 1, num_pt
            new_coord(i,1:3) = cart_coord(i,1:3) - center_coord(1:3)
            r2(i) = dot_product(new_coord(i,:),new_coord(i,:))
        end do

        

        pyval(:) = normalization(alpha,[0,1,0])*(new_coord(:,2)) * exp(-alpha*r2(:))

    end function


    function pz(num_pt,cart_coord, center_coord, alpha) result(pzval)

        integer   :: num_pt
        real(idp) :: cart_coord(num_pt,3), center_coord(3)
        real(idp) :: alpha
        real(idp) :: pzval(num_pt)

        real(idp) :: new_coord(num_pt,3)
        real(idp) :: r2(num_pt)
        
        do i = 1, num_pt
            new_coord(i,1:3) = cart_coord(i,1:3) - center_coord(1:3)
            r2(i) = dot_product(new_coord(i,:),new_coord(i,:))
        end do

        

        pzval(:) = normalization(alpha,[0,0,1])*(new_coord(:,3)) * exp(-alpha*r2(:))

    end function

    function dxy(num_pt,cart_coord, center_coord, alpha) result(dval)
    
        integer   :: num_pt
        real(idp) :: cart_coord(num_pt,3), center_coord(3)
        real(idp) :: alpha
        real(idp) :: dval(num_pt)

        real(idp) :: new_coord(num_pt,3)
        real(idp) :: r2(num_pt)
        
        do i = 1, num_pt
            new_coord(i,1:3) = cart_coord(i,1:3) - center_coord(1:3)
            r2(i) = dot_product(new_coord(i,:),new_coord(i,:))
        end do

        dval(:) = normalization(alpha,[1,1,0])*(new_coord(:,1)) * (new_coord(:,2))* exp(-alpha*r2(:))

    end function
    
end program test_Gauss_expansion