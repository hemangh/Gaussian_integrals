program test_sphHY
    use sphericalHarmonics, only: YLM_2_XYZ, real_spherical_harmonics, associated_legendre
    use matrix_operations, only: invs
    use constants, only : pi
    use precision, only:idp
    implicit none
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pwr
    !! powers of the Cartesian Gaussian in the conventional
    !! ordering allocated with dimensions(1:3,1:(l+1)*(l+2)/2)
    REAL(idp), ALLOCATABLE, DIMENSION(:,:) :: coeff
    !! conversion coefficients between Cartesian and Spherical
    !! Gaussians, dimension(-l:l,1:(l+1)*(l+2)/2)
    REAL(idp), ALLOCATABLE, DIMENSION(:,:) :: invCoef(:,:)
    !! conversion coefficients between Spherical and Cartesian
    !! Gaussians, dimension(1:(l+1)*(l+2)/2, -l:l)
    REAL (idp), allocatable, Dimension(:,:,:) :: Ylm
    real(idp),dimension(:,:),allocatable           :: plm
    integer :: i
    !**********************************************************************
    !        pwr(1:3 -->[x,y,z],1:(l+1)*(l+2)/2) powers of the 
    !                   Cartesian Gaussian in the conventional 
    !                   ordering described in the icart function.
    !**********************************************************************
    !        doubles: coeff(1:2l+1 --> [-m,m],1:(l+1)*(l+2)/2 )
    !                   the coefficient of transformation of a Cartesian
    !                   to a Spherical Gaussian.
    !***********************************************************************
    ! l = 0
    call YLM_2_XYZ(0,pwr,coeff)
    print*, coeff
    !===================================================
    !l = 1
    call YLM_2_XYZ(1,pwr,coeff)
    print*, "Coeff matrix l=1:"
    call print_matrix(coeff)
    call print_matrix_int(pwr)

    !check the individual coefficients:
    if(coeff(3,1) /= sqrt(3._idp/4._idp/pi)) call exit(1)
    invCoef = invs(coeff)
    print*,"inverse matrix l=1"
    call print_matrix(invCoef)
    call print_matrix(matmul(coeff,invCoef))
    print*, sqrt(4._idp * pi/3._idp)
    if(coeff(1,3)/= sqrt(4._idp * pi/3._idp)) call exit(1)
    
    !=======================================================

    ! l = 2

    call YLM_2_XYZ(2,pwr, coeff)
    ! do i = 1, size(pwr,2)
    !     print*, pwr(:,i), coeff(:,i)
    ! end do

    
    print*, "Coeff matrix l=2:"
    call print_matrix(coeff)
    call print_matrix_int(pwr)

    ! Now invert the transformation matrix and print it 
    ! print*, size(coeff,1), size(coeff,2)
    
    invCoef = invs(coeff)
    print*, "(psudo))inverse:"
    call print_matrix(invCoef)
    print*, "check inverse:"
    call print_matrix(matmul(coeff,invCoef))

        
    call YLM_2_XYZ(3, pwr, coeff)

    call print_matrix_int(pwr)
    print*, "Coeff matrix l=3:"
    call print_matrix(coeff)

    invCoef = invs(coeff)
    print*, "(psudo))inverse:"
    call print_matrix(invCoef)
    print*, "check inverse:"
    call print_matrix(matmul(coeff,invCoef))

    
    ! call YLM_2_XYZ(4, pwr, coeff)


    ! call print_matrix_int(pwr)
    ! print*, "Coeff matrix l=4:"
    ! call print_matrix(coeff)

    ! Test Case 2:

    call real_spherical_harmonics(ylm,[(pi/4._idp)],[pi/4._idp] ,1,2)
    allocate(plm(1,0:2))
    call associated_legendre(plm,[cos(pi/4._idp)],1,2,0)
    print*, plm(1,0), plm(1,1)*0.5d0*sqrt((2.d0*real(1)+1.d0)/pi), plm(1,2)
    ! print*, ylm(1,0,0)

    print*, ylm(1,1,0), ylm(1,1,1), Ylm(1,1,-1)
    print*, YLM(1,1,2), YLM(1,1,-2)

    contains

    subroutine print_matrix(A)
        real(8), intent(in) :: A(:,:)
        integer :: m, n
        integer :: i, j
        m = size(A,1)
        n = size(A,2)
        do i = 1, m
            do j = 1, n
                ! Adjust format specifiers as needed for desired output
                write(*, '(F10.6, 1X)', advance='no') A(i, j)
            end do
            print *, "" ! Newline after each row
        end do
    end subroutine print_matrix

    subroutine print_matrix_int(A)
        integer, intent(in) :: A(:,:)
        integer :: m, n
        integer :: i, j
        m = size(A,1)
        n = size(A,2)
        do i = 1, m
            do j = 1, n
                ! Adjust format specifiers as needed for desired output
                write(*, '(I10, 1X)', advance='no') A(i, j)
            end do
            print *, "" ! Newline after each row
        end do
    end subroutine print_matrix_int

    
end program