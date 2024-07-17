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
    ! Orbitals up to L=4 (f) can be used. This value can be changed by use case
    integer, parameter :: max_orbital_l = 4
    ! Coefficient of x^i y^j z^k = B^{LM}_{ijk} r^{i+j+k} Y_{LM}   
    ! (L, M, i, j, k)                                    
    real(idp), allocatable  :: xyz_YLM_coefficient(:,:,:,:,:)!(0:max_l, -max_l:max_l, 0:max_l, 0:max_l, 0:max_l) 
    logical :: set_xyz_YLM = .false. !key to prevent multiple computation of the xyz_YLM_coefficient
    
    ! integrals <Y_l1,m1| Y_l2,m2 | Y_l3,m3>, allocated as (l1,l2,l3,m1,m2,m3)
    real(idp), allocatable :: integrals_3rspH(:,:,:,:,:,:)
    logical :: set_3Ylm_integrals = .false. !key to prevent multiple computation of integrals_3rspH

   

    contains
!***********************************************************************************************
! SUBROUTINE: localcartesian_coeff
! PURPOSE:
!   This subroutine calculates the coefficients for transforming a Cartesian Gaussian function 
!   from a non-local (off-center) coordinate system to a local (centered) coordinate system.
!   This is essential in quantum chemistry when working with molecular orbitals centered 
!   on different atoms.
!
! MATHEMATICAL BASIS:
!   The transformation is based on the binomial expansion:
!      (x - Rx)^l * (y - Ry)^m * (z - Rz)^n = ∑_{i,j,k} A^{lmn}_{ijk} x^i * y^j * z^k
!   where:
!     * Rx, Ry, Rz are the coordinates of the Gaussian's center in the non-local frame.
!     * l, m, n are the powers of x, y, z in the non-local Gaussian.
!     * A^{lmn}_{ijk} are the coefficients to be determined.
!
! INPUT:
!   Rxyz_coord (real(idp)): An array of length 3 containing the Cartesian coordinates 
!                          [Rx, Ry, Rz] of the Gaussian's center in the non-local frame.
!   pwr        (integer): An array of length 3 containing the powers [l, m, n] 
!                          of x, y, z in the non-local Gaussian.
!
! OUTPUT:
!   coeff      (real(idp)): A dynamically allocated 2D array containing the calculated 
!                          coefficients A^{lmn}_{ijk}. The first dimension ranges from 0 to 
!                          the sum of powers (l + m + n), and the second dimension has size 3.
!
! EXTERNAL DEPENDENCIES:
!   This subroutine assumes the existence of a function `binomial_coeff` to calculate binomial coefficients.
!***********************************************************************************************

    subroutine localcartesian_coeff(Rxyz_coord, pwr, coeff)

        ! Declare variables with explicit kinds for better numerical precision
        real(idp), intent(in) :: Rxyz_coord(3) 
        integer, intent(in) :: pwr(3)
        real(idp), allocatable, intent(out) :: coeff(:,:) 
    
        real(idp) :: Rx, Ry, Rz  ! Local variables for clarity
        integer :: l, m, n, k     ! Loop counters and power variables
    
        ! Extract powers from input array for readability
        l = pwr(1)
        m = pwr(2)
        n = pwr(3)
    
        ! Calculate the total sum of powers
        k = sum(pwr)
    
        ! Allocate the coefficient array
        allocate(coeff(0:k, 3))
    
        ! Initialize the coefficient array to zeros
        coeff = 0._idp 
    
        ! Calculate coefficients for the x-term expansion (x - Rx)^l
        do k = 0, l
            coeff(l-k, 1) = binomial_coeff(l, k) * (-Rxyz_coord(1))**k  ! * x**(l-k) is implicit
        end do
    
        ! Calculate coefficients for the y-term expansion (y - Ry)^m
        do k = 0, m
            coeff(m-k, 2) = binomial_coeff(m, k) * (-Rxyz_coord(2))**k  ! * y**(m-k) is implicit
        end do
    
        ! Calculate coefficients for the z-term expansion (z - Rz)^n
        do k = 0, n
            coeff(n-k, 3) = binomial_coeff(n, k) * (-Rxyz_coord(3))**k  ! * z**(n-k) is implicit
        end do
    
    end subroutine localcartesian_coeff
    
!***************************************************************************************!
! FUNCTION: normalization
! PURPOSE:
!   Calculates the normalization factor for a Cartesian Gaussian function. 
!
! INPUT:
!   alpha (real(idp)): The exponent of the Gaussian function (determines width).
!   pwr   (integer): An array containing the powers (nx, ny, nz) in the Gaussian.
!
! OUTPUT:
!   norm  (real(idp)): The calculated normalization factor.
!
! EXTERNAL DEPENDENCIES:
!   This function assumes the existence of an external subroutine `compute_factorials` 
!   and a function `factorial_2n_minus_1` which are used for factorial calculations.
!***************************************************************************************!
    function normalization(alpha, pwr) result(norm)

        ! Declare variables with explicit kinds for precision
        real(idp) :: alpha
        integer :: pwr(3)
        real(idp) :: norm
    
        integer :: pwrsum  ! Temporary variable to hold the sum of powers
    
        ! Calculate the sum of powers (nx + ny + nz)
        pwrsum = sum(pwr) 
    
        ! Compute the required factorials. We're assuming the 
        ! existence of this external function to avoid redundancy.
        call compute_factorials 
    
        ! Calculate the normalization factor based on the formula
        ! The formula is derived from the normalization condition of the Gaussian integral.
        norm = (2.d0 * alpha / pi) ** (3/4) 
        norm = norm * &
        sqrt((4.d0*alpha)**pwrsum / factorial_2n_minus_1(pwr(1)) / factorial_2n_minus_1(pwr(2)) / factorial_2n_minus_1(pwr(3)))
    
    end function



 !***********************************************************************************************
! SUBROUTINE: get_coeff_b
! PURPOSE:
!   Calculates the transformation coefficients (B coefficients) used for converting between 
!   real spherical harmonics (Ylm) and Cartesian basis functions (x^i * y^j * z^k). 
!
! METHODOLOGY:
!   1. Call `YLM_2_XYZ`: Obtains a matrix (`coeff`) relating spherical harmonics (Ylm) to
!      Cartesian terms (x^i * y^j * z^k). This A matrix is inverse of B coefficients.
!   2. Invert Coefficients: Compute the inverse of the `coeff` matrix to get the B coefficients
!      (`invCoef`). These are the transformation coefficients we seek.
!   3. Populate `coeff_BLM_ijk`: Store the B coefficients in the 4D array `coeff_BLM_ijk`, indexed
!      by the Cartesian powers i, j, k, and the magnetic quantum number m.
!   4. Update Global Storage: Store the B coefficients in the global array `xyz_YLM_coefficient` 
!      for potential use elsewhere in the program.
!
! INPUT:
!   l (integer): The angular momentum quantum number (determines the maximum values of i, j, k, and m).
!
! OUTPUT:
!   coeff_BLM_ijk (real(idp)): A 4D array containing the calculated B coefficients, where:
!      - First dimension: Power of x (i)
!      - Second dimension: Power of y (j)
!      - Third dimension: Power of z (k)
!      - Fourth dimension: Magnetic quantum number (m)
!
! EXTERNAL DEPENDENCIES:
!   - YLM_2_XYZ: A subroutine to calculate the A coefficients (spherical to Cartesian).
!   - invs: A function to invert a matrix.
!   - xyz_YLM_coefficient (module variable): A global array to store the B coefficients for wider use.
!
! NOTE: 
!   This code assumes the existence of a module variable `xyz_YLM_coefficient` for storing the 
!   B coefficients globally. If this is not the case, modifications are needed to handle the storage
!   of the calculated coefficients appropriately.
!***********************************************************************************************

subroutine get_coeff_b(l, coeff_BLM_ijk)

    integer, intent(in) :: l
    real(idp), intent(out), allocatable :: coeff_BLM_ijk(:,:,:,:)

    real(idp), allocatable :: coeff(:,:), invCoef(:,:)  ! A coefficients, B coefficients (inverse)
    integer, allocatable :: pwr(:,:)                    ! Powers of x, y, z for each Ylm
    integer :: i, j, k, m                               ! Loop indices and magnetic quantum number
    integer :: count, mcount, max_pwr                   ! Loop counters and maximum power size

    ! 1. Get A coefficients (spherical harmonic to Cartesian)
    call YLM_2_XYZ(l, pwr, coeff)  
    max_pwr = (l + 1) * (l + 2) / 2 

    ! 2. Invert to get B coefficients (Cartesian to spherical harmonic)
    if (allocated(coeff_BLM_ijk)) deallocate(coeff_BLM_ijk)
    allocate(coeff_BLM_ijk(0:l, 0:l, 0:l, -l:l))
    invCoef = invs(coeff) 

    ! 3. Store B coefficients in the output array and the global array
    do count = 1, max_pwr
        i = pwr(1, count)
        j = pwr(2, count)
        k = pwr(3, count)
        mcount = 0
        do m = -l, l
            mcount = mcount + 1
            coeff_BLM_ijk(i, j, k, m) = invCoef(count, mcount)
            xyz_YLM_coefficient(l, m, i, j, k) = coeff_BLM_ijk(i, j, k, m)
        end do
    end do

end subroutine get_coeff_b


!***********************************************************************************************
! Subroutine: expansion_primitiveG_Ylm
! Purpose:
!   This subroutine expands a primitive Cartesian Gaussian function centered at Rxyz_coord into a series
!   of real spherical harmonics (Ylm) evaluated at specified grid points.
!
! Mathematical Background: Reference: H. LE ROUZO, IJQC, Vol. 64, 647-653 (1997)
!   The input Gaussian is of the form: G_{p,q,r,a} = N_{p,q,r,a} (x-Rx)^p (y-Ry)^q (z-Rz)^r exp[-a (r-R_A)^2]
!   The expansion follows these steps:
!   1. Expand the Cartesian terms using binomial theorem: (x-Rx)^p (y-Ry)^q (z-Rz)^r = SUM_{i,j,k} A^{p,q,r}_{i,j,k} x^i y^j z^k
!   2. Express Cartesian terms in spherical coordinates: x^i y^j z^k = r^{gamma} SUM_{L,M} B^{LM}_{i,j,k} * Y_{lm}, 
!      where gamma = i + j + k
!   3. For each l and L, find l' such that |L-l'| <= l <= |L + l'|
!   4. Enforce M + m' = m
!   5. Sum over modified spherical Bessel functions with index l' and Y_{l',m'} (R_A)
!   6. Sum over L, M indices of B^{LM}_{i,j,k} * Clebsh_Gordan(lm,LM|l'm') * Clebsh_Gordan(l0,L0|l'0)
!   7. Apply conditions on m, M, and m' based on sign of m
!
! Input:
!   lmax         (integer) : Maximum spherical harmonic degree (l) to include in the expansion
!   grid_points  (real(idp)): Array of Cartesian coordinates (x,y,z) where to evaluate the expansion
!   powers       (integer) : Powers (p,q,r) of the Cartesian Gaussian
!   Rxyz_coord   (real(idp)): Cartesian coordinates (Rx, Ry, Rz) of the Gaussian center
!   alpha        (real(idp)): Exponent of the Gaussian function
!
! Output:
!   coeff        (real(idp)): Expansion coefficients for each grid point and each (l,m) pair
!***********************************************************************************************

    subroutine expansion_primitiveG_Ylm(lmax, grid_points, powers, Rxyz_coord, alpha,coeff)
        integer, intent(inout) :: lmax
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
        integer :: max_l, int_lmax
        integer, save :: upper_l
        


        real(idp), allocatable :: rthetaphi(:,:)  !(num of points, 1:3 (r,theta,phi)) spherical coord
        real(idp), allocatable :: A_pqr_ijk(:,:)  ! allocated as (0:max(powers),1:3)
        real(idp), allocatable :: coeff_BLM_ijk (:,:,:,:) !(i,j,k,m)
        real(idp), allocatable :: Ylm(:,:,:)
        real(idp), allocatable :: bessel_argument(:), mod_bessel(:,:)
        real(idp) , allocatable :: outer_sum_grid(:), inner_outer_sum_grid(:)
        real(idp) :: inner_sum, A_coeff

        real(idp) ::RAsph_coord(3)
        
        REAL(kind = idp ), allocatable :: si (:)  ! modified spherical Bessel function i_l(x) (0:max_l)
        REAL(kind = idp ), allocatable :: dsi (:) ! first derivative of modified spherical Bessel function i_l(x) [not used]
        real(idp) :: sum
        


                ! Read all the integrals of 3 real spherical harmonics
        ! Goes l1, l2, l3, m1, m2, m3
        if (.not. set_3Ylm_integrals) then

            call read_integrals_3rspH(integrals_3rspH, "../../input/integrals_3rspH.txt", int_lmax)

            set_3Ylm_integrals = .true.
            upper_l =  int_lmax - max_orbital_l
            print'(A,x,I0,x,A)', "Note: only lmax of up to l= ", upper_l, "can be used."

        end if

        if(lmax > upper_l) lmax = upper_l

        max_l = lmax + max_orbital_l 
        
        if (.not. set_xyz_YLM) then
            allocate(xyz_YLM_coefficient(0:max_l, -max_l:max_l, 0:max_l, 0:max_l, 0:max_l))
            xyz_YLM_coefficient = 0._idp
            do l = 0, lmax
                 call get_coeff_b(l, coeff_BLM_ijk) 
             end do
             set_xyz_YLM = .true.
        end if
 

        ! step 0: convert cartesian to spherical
        num_gridpts = size(grid_points,1)
        allocate(rthetaphi(num_gridpts,3))
        do igrid = 1, num_gridpts
            rthetaphi (igrid,:)= cart2sph(grid_points(igrid,:))
        end do

        ! same conversion for the center of the Gaussian in the off-center grid
        RAsph_coord = cart2sph(Rxyz_coord)

        ! allocate the variable arrays
        
        allocate(coeff(num_gridpts,0:lmax,-lmax:lmax), outer_sum_grid(num_gridpts))
        allocate(inner_outer_sum_grid(num_gridpts))
        allocate(bessel_argument(num_gridpts))
        allocate(si(0:max_l), dsi(0:max_l)) 
        allocate(mod_bessel(num_gridpts,0:max_l))
        gamma = sum(powers)


        ! some preperations:
        ! find the modified Bessel functions over grid points:
        bessel_argument(:) = 2._idp * alpha * rthetaphi (:,1) * RAsph_coord(1)
            
        do igrid = 1, num_gridpts
            call sphi(max_l,bessel_argument(igrid),nm,si,dsi)
            do l = 0, max_l
                mod_bessel(igrid,l) = si(l) * exp(-bessel_argument(igrid))
            end do
        end do

        ! find Ylm of center point over ls
        call real_spherical_harmonics(Ylm,RAsph_coord(2),RAsph_coord(3),1,max_l)


        ! step 1: find A coeff
        call localcartesian_coeff ( Rxyz_coord, powers, A_pqr_ijk)

        

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
                            
                                ! inner sum over LM
                                inner_sum = sum_inner_most(l,l2,m,m2,i,j,k)
                            
                                inner_outer_sum_grid(:) =  inner_outer_sum_grid(:) &
                                + mod_bessel(:,l2) * Ylm(1,l2,m2) * inner_sum

                            
                            end do !m2
                            
                        end do !l2

                        outer_sum_grid(:) = outer_sum_grid(:) + A_coeff&
                        * rthetaphi(:,1) ** (i+j+k) * inner_outer_sum_grid(:)
                    end if


                    end do !k
                end do !j
            end do !i
            
            coeff(:,l,m) = 4.d0 * pi* normalization(alpha, powers) &
            * exp(-alpha*(rthetaphi(:,1)-RAsph_coord(1))**2) * outer_sum_grid(:)
            end do !m
        end do !l

    end subroutine



!***********************************************************************************************
! FUNCTION: sum_inner_most
! PURPOSE:
!   Calculates an inner sum within the expansion of a Cartesian Gaussian into spherical harmonics. 
!   This sum involves iterating over intermediate angular momentum (l2) and magnetic quantum number 
!   (m2), combining the B coefficients (for Cartesian to spherical harmonic transformation), and 
!   the precomputed integrals of spherical harmonics.
!
! CONTEXT:
!   This function is called within the main loop of the `expansion_primitiveG_Ylm` subroutine. 
!   It contributes to the overall calculation of the expansion coefficients for each spherical 
!   harmonic (l, m) at each grid point.
!
! INPUT:
!   l     (integer): The target spherical harmonic degree.
!   l1    (integer): An intermediate angular momentum value (used in the Clebsch-Gordan coefficients).
!   m     (integer): The target spherical harmonic order.
!   m1    (integer): An intermediate magnetic quantum number (used in the Clebsch-Gordan coefficients).
!   i, j, k (integer): The powers of x, y, z in the Cartesian Gaussian term.
!
! OUTPUT:
!   suml2 (real(idp)): The calculated inner sum value.
!
! EXTERNAL DEPENDENCIES:
!   - xyz_YLM_coefficient (module variable): A 4D array storing precomputed B coefficients.
!   - integrals_3rspH (module variable): A 6D array storing precomputed integrals of spherical harmonics.
!   - max_l (module variable): Maximum angular momentum for precomputed B coefficients.
!***********************************************************************************************

    function sum_inner_most(l, l1, m, m1, i, j, k) result(suml2)

        integer :: l, l1, m, m1, i, j, k
        real(idp) :: suml2
    
        integer :: l2, m2, l2_min, gamma ! Loop variables and total angular momentum
    
        ! Calculate total angular momentum (gamma = i + j + k)
        gamma = i + j + k
        
        ! Initialize the sum
        suml2 = 0.d0
        
        ! Determine the minimum l2 based on parity of gamma (even or odd)
        l2_min = mod(gamma, 2)
    
        ! Iterate over possible l2 values
        do l2 = l2_min, gamma, 2   ! Step by 2 due to selection rules
                ! Loop over possible m2 values
                do m2 = -l2, l2
                    suml2 = suml2 + xyz_YLM_coefficient(l2, m2, i, j, k) * integrals_3rspH(l1, l2, l, m1, m2, m)
                end do 
        end do 
    
    end function sum_inner_most
    
!***********************************************************************************************
! Subroutine: expansion_primitiveS_Ylm
! Purpose:
!   This subroutine expands a primitive S-type (spherically symmetric) Gaussian function, centered
!   at Rxyz_coord, into a series of real spherical harmonics (Ylm). The expansion is evaluated at
!   specified grid points.
!
! Mathematical Background:
!   The input S-type Gaussian is of the form: 
!       G_{0,0,0,a} = N_{0,0,0,a} exp[-a (r - R)^2]
!   The expansion follows these steps:
!   1. Express the S-type Gaussian in spherical coordinates.
!   2. Compute modified spherical Bessel functions (si) for each grid point.
!   3. Calculate real spherical harmonics (Ylm) for the Gaussian center.
!   4. For each spherical harmonic (l, m):
!      - Combine the Gaussian, ms Bessel function, and spherical harmonic to get the expansion coefficient.
!
! Input:
!   lmax         (integer) : Maximum spherical harmonic degree (l) to include in the expansion
!   grid_points  (real(idp)): Array of Cartesian coordinates (x,y,z) where to evaluate the expansion
!   Rxyz_coord   (real(idp)): Cartesian coordinates (Rx, Ry, Rz) of the Gaussian's center
!   alpha        (real(idp)): Exponent of the Gaussian function
!
! Output:
!   coeff        (real(idp)): Expansion coefficients for each grid point and each (l,m) pair
!***********************************************************************************************

    subroutine expansion_primitiveS_Ylm(lmax, grid_points, Rxyz_coord, alpha,coeff)

        ! Input arguments with intent and kind specification
        integer, intent(in) :: lmax
        real(kind=idp), intent(in) :: grid_points(:,:) 
        real(kind=idp), intent(in) :: Rxyz_coord(3) 
        real(kind=idp), intent(in) :: alpha 
    
        ! Output argument (to be allocated inside)
        real(kind=idp), intent(out), allocatable :: coeff(:,:,:)
    
        ! Local variables with modern Fortran declarations
        integer :: num_gridpts, i, l, m, nm  
        real(kind=idp), dimension(3) :: RAsph_coord
        real(kind=idp) :: ra, r, x, mod_bessel
        real(kind=idp), dimension(:),     allocatable :: si, dsi
        real(kind=idp), dimension(:,:,:), allocatable :: Ylm
    
        ! Constants : pi is taken from the module constants
        ! real(kind=idp), parameter :: pi = 3.14159265358979323846_idp
    
        num_gridpts = size(grid_points, 1)
    
        ! Convert Cartesian coordinates of the Gaussian center to spherical coordinates
        RAsph_coord = cart2sph(Rxyz_coord)
        ra = RAsph_coord(1)  ! Extract the radial distance
    
        ! Calculate real spherical harmonics for the Gaussian center
        call real_spherical_harmonics(Ylm, [RAsph_coord(2)], [RAsph_coord(3)], 1, lmax)
    
        ! Allocate memory for the expansion coefficients
        allocate(coeff(num_gridpts, 0:lmax, -lmax:lmax))

        !Allocate memory for the modified bessel (si) and its derivative (dsi)
        allocate(si(0:lmax), dsi(0:lmax))
    
        ! Iterate over all grid points
        do i = 1, num_gridpts
            ! Calculate the radial distance of the grid point from the origin
            r = sqrt(dot_product(grid_points(i,:), grid_points(i,:)))
            ! Calculate the argument for the modified spherical Bessel function
            x = 2.0_idp * r * ra * alpha
    
            ! Calculate modified spherical Bessel functions up to order 4 (sufficient for S-type Gaussian)
            call sphi(lmax, x, nm, si, dsi)
    
            ! Iterate over spherical harmonic degrees (l)
            do l = 0, lmax
                ! Calculate the modified spherical Bessel function for the current l
                mod_bessel = si(l) * exp(-x)
    
                ! Iterate over spherical harmonic orders (m)
                do m = -l, l
                    ! Calculate the expansion coefficient for the current (l, m) pair
                    coeff(:, l, m) = 4.0_idp * pi * (2.0_idp * alpha / pi)**(3.0_idp / 4.0_idp) &
                        * exp(-alpha * (r - ra)**2) * mod_bessel * Ylm(1, l, m)
                end do
            end do
        end do
    
    end subroutine expansion_primitiveS_Ylm
    

        


end module


    ! !
    ! ! \sum_{l3, m3} \int Y_{l1,m1} * Y_{l2, m2} * Y_{l3, m3}
    ! !
    ! function sum_over_Clebsh_Gordon_Constants(l1, l2, l3, m1, m2, m3) result(sum3)
    
    ! integer, intent(in) :: l1, m1, l2, m2, l3, m3
    ! real(idp) :: sum3

    ! sum3 = 0._idp


    ! if(m3 == m1+m2 .and. l3 <= abs(l1+l2) .and. l3 >= abs(l1-l2)) then
    
    !     sum3 = (2*l1+1) * (2*l2+1)/4._idp/ pi/ (2*l3+1)
    !     sum3 = sqrt(sum3) / 2._idp /pi  ! => this is real spherical harmonics extra factor
    !     sum3 = sum3 * calculate_cg_coefficients(l1,l2,l3,0,0) * calculate_cg_coefficients(l1,l2,l3,m1,m2) 
    !     if(m1>0 .and. m2>0) then
    !         sum3 = sum3 * pi /2.d0
    !     else if(m1>0 .and. m2<0 .and. m3<0) then
    !         sum3 = sum3 * pi/2.d0
    !     else if(m1<0 .and. m3>0 .and. m3<0) then 
    !         sum3 = sum3 * pi/2.d0
    !     else if (m1==0 .and. m2==0) then
    !         sum3 = sum3 * 2.d0 * pi
    !     elseif ( m1==m3 .or.  m2==m3) then
    !         sum3 = sum3 * 2.d0* pi
    !     else
    !         sum3 = 0._idp
    !     end if

    ! endif
    
    
    ! end function sum_over_Clebsh_Gordon_Constants

    !     function sum_3(i,j,k, l1,l2,m1,m2,lmax) result(sum3)
    !         integer :: i,j,k
    !         integer :: l1,l2,m1,m2,lmax
    !         real(idp) :: sum3

    !         integer :: l, m
    !         real(idp) , allocatable :: coeff_b(:,:,:,:)
            
            
    !         sum3 = 0._idp
    !         do l = 0, lmax
    !             do m = 1, 2*l+1
    !                 !calculate
    !                 call get_coeff_b(l, coeff_b)
    !                 sum3 = sum3 +  coeff_b(i,j,k, m) *sum_over_Clebsh_Gordon_Constants(l1,l2,l, m1, m2, m)
    !             end do
    !         end do

    !     end function sum_3