module Gaussian_expansion_coefficient_Ylm
use precision, only: idp
use constants, only: pi
use factorials_mod, only: binomial_coeff
use conversions, only: cart2sph
use special_functions, only: ikna
use sphericalHarmonics, only: real_spherical_harmonics, YLM_2_XYZ
use clebsch_gordan_mod, only: calculate_cg_coefficients
use matrix_operations, only: invs
use factorials_mod, only: compute_factorials, factorial_2n_minus_1
    implicit none
    integer, parameter :: max_l = 2
    ! Coefficient of x^i y^j z^k = B^{LM}_{ijk} r^{i+j+k} Y_{LM}   
    ! (L, M, i, j, k)                                    
    ! real(idp) :: xyz_YLM_coefficient(0:max_l, -max_l:max_l, 0:max_l, 0:max_l, 0:max_l) 

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
        integer, intent(in)   :: pwr(3)

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
    function expansion_exp_Bessel_i_Ylm(a, r, RAsph_coord, ylm, lmin, lmax) result(term)
        real(idp) :: a, r, RAsph_coord(3)
        real(idp) :: ylm
        integer   :: lmin, lmax

        real(idp) :: term
        
        integer :: nm
        real(idp) , dimension(0:lmax) :: bi, di, bk, dk
        real(idp) :: x
        integer :: l, m

        term = 0._idp
        x = 2._idp * a * r * RAsph_coord(1)
        
        call ikna(lmax, x, nm, bi, di, bk, dk)
        ! call real_spherical_harmonics(Ylm,RAsph_coord(2),RAsph_coord(3),1,lmax) !allocated according to phi/theta dim & lmax in the routine
         
        term = term + exp(-x) * bi(l) * Ylm ! * B^{lm}_{ijk} * Ylm(r) * 
  

        end function expansion_exp_Bessel_i_Ylm

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
                    sum3 = sum3 * pi
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
        ! 5) sum over modified spherical Bessels with index l' and Y_{l',m'} (R_A)
        ! 7) sum over L, M index of B^{LM}_{i,j,k} * Clebsh_Gordan(lm,LM|l'm') * Clebsh_Gordan(l0,L0|l'0) 
        ! 8) if m >=0 M >=0 and m'>=0 else if m<0 either M>=0 and m'<0 or M<0 and m'>=0 
        ! 

    subroutine expansion_primitiveG_Ylm(lmax, grid_points, powers, Rxyz_coord, alpha,coeff)
        integer, intent(in) :: lmax
        real(idp), intent(in) :: grid_points(:,:) !(num of points, 1:3 (x,y,z)) cartesian coord
        integer, intent(in) :: powers(3) ! (x,y,z)
        real(idp), intent(in) :: Rxyz_coord(3) !(Rx, Ry, Rz)
        real(idp), intent(in) :: alpha ! exponent of Gaussian
        real(idp), intent(out), allocatable :: coeff(:,:)  !(num_gridpts,0:lmax)

        
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
        real(idp) :: inner_sum

        real(idp) ::RAsph_coord(3)
        real(idp) , dimension(0:lmax) :: bi, di, bk, dk
        real(idp) :: sum
        
        ! step 0: convert cartesian to spherical
        num_gridpts = size(grid_points,1)
        allocate(rthetaphi(num_gridpts,3))
        do igrid = 1, num_gridpts
            rthetaphi (igrid,:)= cart2sph(grid_points(igrid,:))
        end do

        ! same conversion for the center of the Gaussian in the off-center grid
        RAsph_coord = cart2sph(Rxyz_coord)

        ! allocate the variable arrays
        allocate(coeff(num_gridpts,0:lmax), outer_sum_grid(num_gridpts))
        allocate(inner_outer_sum_grid(num_gridpts))
        allocate(bessel_argument(num_gridpts))
        allocate(mod_bessel(num_gridpts,0:max(max_l,lmax)))

        ! some preperations:
        ! find the modified Bessel functions over grid points:
        bessel_argument(:) = 2._idp * alpha * rthetaphi (:,1) * RAsph_coord(1)
            
        do igrid = 1, num_gridpts
            call ikna(max(max_l,lmax), bessel_argument(igrid), nm, bi, di, bk, dk)
            do l = 0, max(max_l,lmax)
                mod_bessel(igrid,l) = bi(l) * exp(-bessel_argument(igrid))
            end do
        end do

        ! find Ylm of center point over ls
        call real_spherical_harmonics(Ylm,RAsph_coord(2),RAsph_coord(3),1,max(max_l,lmax))


        ! step 1: find A coeff
        call localcartesian_coeff ( Rxyz_coord, powers, A_pqr_ijk)

        ! step 2: iterate over l => this is C^{lm}_{apqr}(r) = <G_{pqr} (a, x,y,z, R) | Y_{lm} (omega)> 
        ! do l = 0, lmax
        l = 1
        m = 1
            outer_sum_grid = 0._idp
            ! step 3: find i, j, k combinations 
            do i = 0,size(A_pqr_ijk(:,1))-1
                do j = 0,size(A_pqr_ijk(:,2))-1
                    do k = 0,size(A_pqr_ijk(:,3))-1
                        gamma = i + j + k

                        ! do m = -l, l
                            ! find all L 
                            inner_outer_sum_grid = 0._idp
                            do l2 = 0, max_l
                                do m2 = -l2, l2
                                    inner_sum = 0._idp
                                    ! do l1 = gamma, 0, -2
                                       
                                    !     if(l1<0) exit
                                        
                                    !     ! step 4: find B coeff
                                    !     call get_coeff_b(l1, coeff_BLM_ijk)

                                    !     !find all l' => Sum_{l2,m2} Exp(-2arR) i_l2(2arR) Y_{l2,m2} (\hat{R})
                                    !     l2_min = abs(l-l1)
                                    !     l2_max = l+l1

                                    !     if (l2 <= l2_max .and. l2 <= l2_min) then
                                    !         ! Adjust loop bounds for m2
                                    !         if (l2 == 0) then
                                    !             m2 = 0
                                    !             ! Adjust loop bounds for m1
                                    !             if (l1 == 0) then
                                    !                 m1 = 0
                                    !                 if (m1 + m2 == m) then
                                    !                     sum = sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
                                    !                     inner_sum_grid(:) = inner_sum_grid(:) + &
                                    !                      coeff_BLM_ijk(i, j, k, m) * sum
                                    !                 end if
                                    !             else
                                    !                 do m1 = -l1, l1
                                    !                     if (m1 + m2 == m) then
                                    !                         sum = sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
                                    !                         inner_sum_grid(:) = inner_sum_grid(:) + &
                                    !                          coeff_BLM_ijk(i, j, k, m) * sum
                                    !                     end if
                                    !                 end do
                                    !             end if
                                    !         else
                                    !             do m2 = -l2, l2
                                    !                 ! Adjust loop bounds for m1
                                    !                 if (l1 == 0) then
                                    !                     m1 = 0
                                    !                     if (m1 + m2 == m) then
                                    !                         sum = sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
                                    !                         inner_sum_grid(:) = inner_sum_grid(:) + &
                                    !                         coeff_BLM_ijk(i, j, k, m) * sum
                                    !                     end if
                                    !                 else
                                    !                     do m1 = -l1, l1
                                    !                         if (m1 + m2 == m) then
                                    !                             sum = sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
                                    !                             inner_sum_grid(:) = inner_sum_grid(:) + &
                                    !                              coeff_BLM_ijk(i, j, k, m) * sum
                                    !                         end if
                                    !                     end do
                                    !                 end if
                                    !             end do
                                    !         end if
                                    !     end if
                                        

                                        
                                    ! end do !l1
                                
                                    inner_sum = sum_inner_most(l,l2,m,m2,i,j,k)
                                
                                    inner_outer_sum_grid(:) =  inner_outer_sum_grid(:) &
                                    + mod_bessel(:,l2) * Ylm(1,l2,m2) * inner_sum
                                
                                end do !m2
                                
                            end do !l2

                            outer_sum_grid(:) = outer_sum_grid(:) + A_pqr_ijk(i,1) * A_pqr_ijk(j,2) * A_pqr_ijk(k,3) &
                            * rthetaphi(:,1) ** gamma * inner_outer_sum_grid(:)

                        ! end do !m

                    end do !k
                end do !j
            end do !i
            print*, "norm:", normalization(alpha,powers)
            print*, "sum:", outer_sum_grid(:), "exp:",exp(-alpha*(rthetaphi(:,1)-RAsph_coord(1))**2)
            coeff(:,l) = 4.d0 * pi* normalization(alpha, powers) &
            * exp(-alpha*(rthetaphi(:,1)-RAsph_coord(1))**2) * outer_sum_grid(:)

        ! end do !l

    end subroutine



    function sum_inner_most(l,l1,m,m1, i,j,k) result(suml2)

        integer:: l, l1, m, m1, i,j,k

        real(idp):: suml2

        integer :: l2, m2, l2_min, gamma
        real(idp), allocatable :: coeff_BLM_ijk(:,:,:,:)

        gamma = i + j + k

        

        if (gamma > max_l) stop "increse max_l parameter!"
        
        suml2 = 0.d0
        
        l2_min = mod(gamma,2)


        do l2 = l2_min, gamma, 2

            call get_coeff_b(l2, coeff_BLM_ijk)

            if (l2 == 0) then

                m2 = 0
                suml2 = suml2 + coeff_BLM_ijk(i,j,k,m2) * sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)
 

            else

                do m2 = -l2, l2

                    suml2 = suml2 + coeff_BLM_ijk(i,j,k,m2) * sum_over_Clebsh_Gordon_Constants(l1, l2, l, m1, m2, m)

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