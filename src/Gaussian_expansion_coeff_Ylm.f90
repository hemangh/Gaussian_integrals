module Gaussian_expansion_coefficient_Ylm
use precision, only: idp
use constants, only: pi
use factorials_mod, only: binomial_coeff
use conversions, only: cart2sph
use special_functions, only: sphi
use sphericalHarmonics, only: real_spherical_harmonics, YLM_2_XYZ
use clebsch_gordan_mod, only: calculate_cg_coefficients
use matrix_operations, only: invs
use factorials_mod, only: compute_factorials, factorial_2n_minus_1
use integrals_3rspH_mod, only: read_integrals_3rspH
    implicit none
    integer, parameter :: max_l =5 
    ! Coefficient of x^i y^j z^k = B^{LM}_{ijk} r^{i+j+k} Y_{LM}   
    ! (L, M, i, j, k)                                    
    real(idp) :: xyz_YLM_coefficient(0:max_l, -max_l:max_l, 0:max_l, 0:max_l, 0:max_l) 
    logical :: set_xyz_YLM = .false.
    ! integrals <Y_l1,m1| Y_l2,m2 | Y_l3,m3>, allocated as (l1,l2,l3,m1,m2,m3)
    real(idp), allocatable :: integrals_3rspH(:,:,:,:,:,:)

   

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
        integer, intent(in)   :: pwr(3)

        real(idp), allocatable, intent(out) :: coeff(:,:) !result

        real(idp):: Rx,Ry, Rz
        integer  :: l, m, n
        
        integer :: k

        l = pwr(1)
        m = pwr(2)
        n = pwr(3)
        k = sum(pwr)
        allocate(coeff(0:k,3))
        coeff = 0._idp

        !(x-Rx)^l => \sum_{k=0}^l binomial_coeff(l,k) * (x)^{l-k} * Rx^k 
        do k = 0, l
            coeff(l-k,1) = binomial_coeff(l,k) * (-Rxyz_coord(1))**k !* x**(l-k)
        end do

        !(y-Ry)^m => \sum_{k=0}^m binomial_coeff(m,k) * (y)^{m-k} * Ry^k 
        do k = 0, m
            coeff(m-k,2) = binomial_coeff(m,k) * (-Rxyz_coord(2))**k !* y**(m-k)
        end do

        !(z-Rz)^n => \sum_{k=0}^n binomial_coeff(n,k) * (z)^{n-k} * Rz^k 
        do k = 0, n
            coeff(n-k,3) = binomial_coeff(n,k) * (-Rxyz_coord(3))**k !* z**(n-k)
        end do

    end subroutine localcartesian_coeff

    function normalization(alpha, pwr) result(norm)

        real(idp) :: alpha
        integer :: pwr(3)

        real(idp) :: norm

        integer :: pwrsum

        pwrsum = sum (pwr)

        call compute_factorials

        norm = (2.d0 * alpha / pi) ** (3/4)
        norm = norm * &
        sqrt((4.d0*alpha)**pwrsum/factorial_2n_minus_1(pwr(1))/factorial_2n_minus_1(pwr(2))/factorial_2n_minus_1(pwr(3)))

    end function

    ! 
    ! \sum_{l,m} exp(-2*a*RA*r)* i_l(2*a*RA*r) * Y_{lm}(\hat{RA})
    ! 
    ! function expansion_exp_Bessel_i_Ylm(a, r, RAsph_coord, ylm, lmin, lmax) result(term)
    !     real(idp) :: a, r, RAsph_coord(3)
    !     real(idp) :: ylm
    !     integer   :: lmin, lmax

    !     real(idp) :: term
        
    !     integer :: nm
    !     real(idp) , dimension(0:lmax) :: bi, di, bk, dk
    !     real(idp) :: x
    !     integer :: l, m

    !     term = 0._idp
    !     x = 2._idp * a * r * RAsph_coord(1)
        
    !     call ikna(lmax, x, nm, bi, di, bk, dk)
    !     ! call real_spherical_harmonics(Ylm,RAsph_coord(2),RAsph_coord(3),1,lmax) !allocated according to phi/theta dim & lmax in the routine
         
    !     term = term + exp(-x) * bi(l) * Ylm ! * B^{lm}_{ijk} * Ylm(r) * 
  

    !     end function expansion_exp_Bessel_i_Ylm

        !
        ! \sum_{l3, m3} \int Y_{l1,m1} * Y_{l2, m2} * Y_{l3, m3}
        !
         function sum_over_Clebsh_Gordon_Constants(l1, l2, l3, m1, m2, m3) result(sum3)
            
            integer, intent(in) :: l1, m1, l2, m2, l3, m3
            real(idp) :: sum3

            sum3 = 0._idp


            if(m3 == m1+m2 .and. l3 <= abs(l1+l2) .and. l3 >= abs(l1-l2)) then
            
                sum3 = (2*l1+1) * (2*l2+1)/4._idp/ pi/ (2*l3+1)
                sum3 = sqrt(sum3) / 2._idp /pi  ! => this is real spherical harmonics extra factor
                sum3 = sum3 * calculate_cg_coefficients(l1,l2,l3,0,0) * calculate_cg_coefficients(l1,l2,l3,m1,m2) 
                if(m1>0 .and. m2>0) then
                    sum3 = sum3 * pi /2.d0
                else if(m1>0 .and. m2<0 .and. m3<0) then
                    sum3 = sum3 * pi/2.d0
                else if(m1<0 .and. m3>0 .and. m3<0) then 
                    sum3 = sum3 * pi/2.d0
                else if (m1==0 .and. m2==0) then
                    sum3 = sum3 * 2.d0 * pi
                elseif ( m1==m3 .or.  m2==m3) then
                    sum3 = sum3 * 2.d0* pi
                else
                    sum3 = 0._idp
                end if

            endif
            
            
         end function sum_over_Clebsh_Gordon_Constants

        function sum_3(i,j,k, l1,l2,m1,m2,lmax) result(sum3)
            integer :: i,j,k
            integer :: l1,l2,m1,m2,lmax
            real(idp) :: sum3

            integer :: l, m
            real(idp) , allocatable :: coeff_b(:,:,:,:)
            
            
            sum3 = 0._idp
            do l = 0, lmax
                do m = 1, 2*l+1
                    !calculate
                    call get_coeff_b(l, coeff_b)
                    sum3 = sum3 +  coeff_b(i,j,k, m) *sum_over_Clebsh_Gordon_Constants(l1,l2,l, m1, m2, m)
                end do
            end do

        end function sum_3

        subroutine get_coeff_b(l, coeff_BLM_ijk)

                integer, intent(in) :: l
                real(idp), intent(out), allocatable :: coeff_BLM_ijk (:,:,:,:)

                real(idp), allocatable :: coeff(:,:) , invCoef(:,:)
                integer, allocatable   :: pwr(:,:)
                integer :: i,j,k, m
                integer :: count, mcount, max_pwr
                
                
                call YLM_2_XYZ(l,pwr,coeff)
                max_pwr = (l+1) * (l+2) /2
                if (allocated(coeff_BLM_ijk)) deallocate(coeff_BLM_ijk)
                allocate(coeff_BLM_ijk(0:l, 0:l, 0:l, -l:l))
                invCoef = invs(coeff)
                do count = 1, max_pwr
                    i = pwr(1,count)
                    j = pwr(2,count)
                    k = pwr(3,count)
                    mcount = 0
                    do m = -l, l
                        mcount = mcount+1
                        coeff_BLM_ijk(i,j,k, m) = invCoef(count,mcount)
                        xyz_YLM_coefficient(l,m,i,j,k) = coeff_BLM_ijk(i,j,k, m)
                    end do
                end do

        end subroutine



    !from here the input is considered:
        ! G_{p,q,r,a} = N_{p,q,r,a} (x-Rx)^p (y-Ry)^q (z-Rz)^r exp[-a (r-R_A)^2]
        ! 1) (x-Rx)^p (y-Ry)^q (z-Rz)^r = SUM_{i,j,k} A^{p,q,r}_{i,j,k} x^i y^j z^i
        ! 2) x^i y^j z^i =  r^{gamma} SUM_{L,M} B^{LM}_{i,j,k} * Y_{lm}, where gamma = i + j+ k and L = {0,2,...,gamma}
        !    if gamma is even and L = {1,3,..,gamma} if gamma is odd 
        ! 3) for each l and L pair find l' so that |L-l'| <= l <= |L + l'|
        ! 4) M + m' = m 
        ! 5) sum over modified spherical Bessels with index l' and Y_{l',m'} (R_A) there is a factor of 2 if m != 0 
        ! 7) sum over L, M index of B^{LM}_{i,j,k} * Clebsh_Gordan(lm,LM|l'm') * Clebsh_Gordan(l0,L0|l'0) 
        ! 8) if m >=0 M >=0 and m'>=0 else if m<0 either M>=0 and m'<0 or M<0 and m'>=0 
        ! 

    subroutine expansion_primitiveG_Ylm(lmax, grid_points, powers, Rxyz_coord, alpha,coeff)
        integer, intent(in) :: lmax
        real(idp), intent(in) :: grid_points(:,:) !(num of points, 1:3 (x,y,z)) cartesian coord
        integer, intent(in) :: powers(3) ! (x,y,z)
        real(idp), intent(in) :: Rxyz_coord(3) !(Rx, Ry, Rz)
        real(idp), intent(in) :: alpha ! exponent of Gaussian
        real(idp), intent(out), allocatable :: coeff(:,:,:)  !(num_gridpts,0:lmax,-lmax:lmax)

        
        integer :: i, j, k, l, m
        integer :: gamma
        integer :: igrid, num_gridpts
        integer :: nm
        integer :: l1, l2_min, l2_max, l2, m1, m2
        real(idp), allocatable :: rthetaphi(:,:)  !(num of points, 1:3 (r,theta,phi)) spherical coord
        real(idp), allocatable :: A_pqr_ijk(:,:)  ! allocated as (0:max(powers),1:3)
        real(idp), allocatable :: coeff_BLM_ijk (:,:,:,:) !(i,j,k,m)
        real(idp), allocatable :: Ylm(:,:,:)
        real(idp), allocatable :: bessel_argument(:), mod_bessel(:,:)
        real(idp) , allocatable :: outer_sum_grid(:), inner_outer_sum_grid(:)
        real(idp) :: inner_sum, A_coeff

        real(idp) ::RAsph_coord(3)
        
        integer:: lmaximum ! max(lmax, max_l)
        REAL(kind = idp ), allocatable :: si (:)  ! modified spherical Bessel function i_l(x) (0:lmaximum)
        REAL(kind = idp ), allocatable :: dsi (:) ! first derivative of modified spherical Bessel function i_l(x) [not used]
        real(idp) :: sum
        
        if (.not. set_xyz_YLM) then
            xyz_YLM_coefficient = 0._idp
            do l = 0, lmax
                 call get_coeff_b(l, coeff_BLM_ijk) ! This has to change.
             end do
             set_xyz_YLM = .true.
        end if
        ! print*,xyz_YLM_coefficient(0,0,0,0,0)
        ! print*,xyz_YLM_coefficient(1,0,0,0,1)
        ! print*,xyz_YLM_coefficient(1,1,1,0,0)
        ! print*,xyz_YLM_coefficient(1,-1,0,1,0)

        ! print*,xyz_YLM_coefficient(2,0,0,0,2)
        ! print*,xyz_YLM_coefficient(2,1,1,0,1)
        ! print*,xyz_YLM_coefficient(2,-1,0,1,1)
        ! print*,xyz_YLM_coefficient(2,2,0,2,0), xyz_YLM_coefficient(2,2,2,0,0)
        ! print*,xyz_YLM_coefficient(2,-2,1,1,0)
        ! step 0: convert cartesian to spherical

        ! Read all the integrals of 3 real spherical harmonics
        ! Goes l1, l2, l3, m1, m2, m3
        call read_integrals_3rspH(integrals_3rspH, "../../input/integrals_3rspH.txt")

        num_gridpts = size(grid_points,1)
        allocate(rthetaphi(num_gridpts,3))
        do igrid = 1, num_gridpts
            rthetaphi (igrid,:)= cart2sph(grid_points(igrid,:))
        end do

        ! same conversion for the center of the Gaussian in the off-center grid
        RAsph_coord = cart2sph(Rxyz_coord)

        ! allocate the variable arrays
        lmaximum = max(max_l,lmax)
        allocate(coeff(num_gridpts,0:lmax,-lmax:lmax), outer_sum_grid(num_gridpts))
        allocate(inner_outer_sum_grid(num_gridpts))
        allocate(bessel_argument(num_gridpts))
        allocate(si(0:lmaximum), dsi(0:lmaximum)) 
        allocate(mod_bessel(num_gridpts,0:lmaximum))

        ! some preperations:
        ! find the modified Bessel functions over grid points:
        bessel_argument(:) = 2._idp * alpha * rthetaphi (:,1) * RAsph_coord(1)
            
        do igrid = 1, num_gridpts
            call sphi(lmaximum,bessel_argument(igrid),nm,si,dsi)
            do l = 0, lmaximum
                mod_bessel(igrid,l) = si(l) * exp(-bessel_argument(igrid))
            end do
        end do

        ! find Ylm of center point over ls
        call real_spherical_harmonics(Ylm,RAsph_coord(2),RAsph_coord(3),1,lmaximum)


        ! step 1: find A coeff
        call localcartesian_coeff ( Rxyz_coord, powers, A_pqr_ijk)

        gamma = sum(powers)

        coeff = 0._idp
        ! step 2: iterate over l => this is C^{lm}_{apqr}(r) = <G_{pqr} (a, x,y,z, R) | Y_{lm} (omega)> 
       
        do l = 0, lmax
        
           do m = -l, l
        ! l =0; m=0
            outer_sum_grid = 0._idp
            ! step 3: find i, j, k combinations 
            do i = 0,gamma 
                do j = 0,gamma-i
                   do  k = 0,gamma - i - j
                    
                    A_coeff = A_pqr_ijk(i,1) * A_pqr_ijk(j,2) * A_pqr_ijk(k,3)
                    if ( A_coeff /= 0._idp) then
                        ! find all l prime => l2
                        inner_outer_sum_grid = 0._idp
                        do l2 = 0, max_l
                            do m2 = -l2, l2
                                ! inner_sum = 0._idp
                            
                                ! inner sum over LM
                                inner_sum = sum_inner_most(l,l2,m,m2,i,j,k)
                            
                                inner_outer_sum_grid(:) =  inner_outer_sum_grid(:) &
                                + mod_bessel(:,l2) * Ylm(1,l2,m2) * inner_sum

                                ! print'("real sphY_l=",I0,",m=",I0,":",E15.8,", B_",I0,":",E15.8)', &
                                ! l2,m2,Ylm(1,l2,m2),l, mod_bessel(2,l2)
                            
                            end do !m2
                            
                        end do !l2

                        outer_sum_grid(:) = outer_sum_grid(:) + A_coeff&
                        * rthetaphi(:,1) ** (i+j+k) * inner_outer_sum_grid(:)
                    end if


                    end do !k
                end do !j
            end do !i
            
        
            !print*, "norm:", normalization(alpha,powers)
            
            !print*, "sum:", outer_sum_grid(:), "exp:",exp(-alpha*(rthetaphi(:,1)-RAsph_coord(1))**2)
            coeff(:,l,m) = 4.d0 * pi* normalization(alpha, powers) &
            * exp(-alpha*(rthetaphi(:,1)-RAsph_coord(1))**2) * outer_sum_grid(:)
            end do !m
        end do !l

    end subroutine



    function sum_inner_most(l,l1,m,m1, i,j,k) result(suml2)

        integer:: l, l1, m, m1, i,j,k

        real(idp):: suml2

        integer :: l2, m2, l2_min, gamma
        !real(idp), allocatable :: coeff_BLM_ijk(:,:,:,:)

        gamma = i + j + k

        

        if (gamma > max_l) stop "increse max_l parameter!"
        
        suml2 = 0.d0
        
        l2_min = mod(gamma,2)


        do l2 = l2_min, gamma, 2

            !call get_coeff_b(l2, coeff_BLM_ijk) ! This has to change.

            if (l2 == 0) then

                m2 = 0
                ! suml2 = suml2 + coeff_BLM_ijk(i,j,k,m2) * sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
                suml2 = suml2 + xyz_YLM_coefficient(l2,m2,i,j,k) * integrals_3rspH(l1, l2, l, m1, m2, m)
 

            else

                do m2 = -l2, l2

                    ! suml2 = suml2 + coeff_BLM_ijk(i,j,k,m2) * sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
                    suml2 = suml2 + xyz_YLM_coefficient(l2,m2,i,j,k) * integrals_3rspH(l1, l2, l, m1, m2, m)

                end do !m2

            end if

        end do !l2

    end function

    ! function sum_l1(l,m, sum2) result(suml1)

    !     integer:: l, m
    !     real(idp) :: sum2
    !     real(idp) :: suml1

    !     integer :: l1, l1_min, l1_max

    !     sum_l1 = 0.d0

    !     do l1 = 0, max_l

    !         if (l1 == 0) then

    !             m1 = 0
    !             suml1 = suml1 + sum_l2(l,m,l1,m1,coef_bl1, gamma)
 

    !         else

    !             do m2 = -l2, l2

    !                 suml2 = sum_l2 + coef_bl1 * sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)

    !             end do !m2

    !         end if

    !     end do !l2

    ! end function

end module