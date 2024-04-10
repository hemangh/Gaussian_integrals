program test_sphHY
    use sphericalHarmonics, only: YLM_2_XYZ
    use matrix_operations, only: invs
    use precision, only:idp
    implicit none
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pwr
    !! powers of the Cartesian Gaussian in the conventional
    !! ordering allocated with dimensions(1:3,1:(l+1)*(l+2)/2)
    REAL(idp), ALLOCATABLE, DIMENSION(:,:) :: coeff, invCoef(:,:)
    !! conversion coefficients between Cartesian and Spherical
    !! Gaussians, dimension(-l:l,1:(l+1)*(l+2)/2)
    integer :: i

    call YLM_2_XYZ(0,pwr,coeff)
    print*, coeff

    call YLM_2_XYZ(1,pwr,coeff)
    print*, "Coeff matrix l=1:"
    call print_matrix(coeff)
    call print_matrix_int(pwr)

    call YLM_2_XYZ(2,pwr, coeff)
    ! do i = 1, size(pwr,2)
    !     print*, pwr(:,i), coeff(:,i)
    ! end do

    call print_matrix_int(pwr)
    print*, "Coeff matrix l=2:"
    call print_matrix(coeff)

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