program matrix_operation_tests
    use matrix_operations
    implicit none
    integer :: n  ! Size of square matrices
    real(8), allocatable :: A(:,:), invA(:,:)
    logical :: result
    integer :: i

    ! Test 1: Identity Matrix Test
    n = 3
    allocate(A(n,n), invA(n,n))
    A = 0.0
    do i = 1, n
        A(i,i) = 1.0
    end do
    invA = invs(A)
    result = matrix_equal(matmul(A, invA), identity(n))
    print *, "Identity Matrix Test Passed:", result

    ! Test 2: Random Square Matrix Test
    deallocate(invA, A)
    call random_matrix(A, 3)
    
    ! invA = invs(A)
    ! call print_matrix(A)
    invA = invs(A)
    ! print*,  matmul( invA, A)
    result = matrix_equal(matmul(A, invA), identity(3))
    print *, "Random Square Matrix Test Passed:", result

    ! Test 3: Singular Matrix Test
    A = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], [n, n]) ! Singular matrix
    invA = invs(A)
    print *, "Singular Matrix Test Passed:", result

    ! Test 4: Orthogonal Matrix Test
    deallocate(A, invA)
    allocate(A(2,2))
    A = reshape([0., 1., -1., 0.], [2, 2]) ! Orthogonal matrix
    invA = invs(A)
    result = matrix_equal(matmul(A, invA), identity(2))
    print *, "Orthogonal Matrix Test Passed:", result

    ! Test 5: Ill-conditioned Matrix Test
    A = reshape([1., 0.999, 0.999, 0.998], [2, 2]) ! Ill-conditioned matrix
    invA = invs(A)
    result = matrix_equal(matmul(A, invA), identity(2))
    print *, "Ill-conditioned Matrix Test Passed:", result

    ! Test 6: Large Matrix Test
    n = 500
    deallocate(A, invA)
    allocate(A(n,n), invA(n,n))
    call random_matrix(A, n) ! Large random matrix
    invA = invs(A)
    result = matrix_equal(matmul(A, invA), identity(n))
    print *, "Large Matrix Test Passed:", result

    ! Test 7: Rectangular Matrix Psudoinverse
    deallocate(A, invA)
    allocate(A(3,2), invA(2,3))
    A = reshape([1.d0, 0.d0, 1.d0, -1.d0, 1.d0, 0.d0],[3,2])
    
    print*, "A (before):"
    call print_matrix(A)
    invA = invs(A)
    print*, "A (after):"
    call print_matrix(A)
    print*, "invsA*3:"
    call print_matrix(invA*3.d0)
    print*, "Ainv*A:"
    call print_matrix(matmul(invA,A))
    result = matrix_equal(matmul( invA, A), identity(2))
    print *, "Rectangular Matrix Psudoinverse Passed:", result


    print*, "A*Ainv:"
    call print_matrix(matmul(A,invA))

    contains

    ! Function to check if matrices are equal within a tolerance
    logical function matrix_equal(A, B)
        real(8), intent(in) :: A(:,:), B(:,:)
        real(8), parameter :: tolerance = 1.0E-10
        integer :: i, j

        matrix_equal = .true.
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                if (abs(A(i,j) - B(i,j)) > tolerance) then
                    matrix_equal = .false.
                    return
                end if
            end do
        end do
    end function matrix_equal

! Function to generate a random matrix
    subroutine random_matrix(A,  n)
        real(8), allocatable, intent(out) :: A(:,:)
        integer, intent(in) ::  n
        integer :: i, j
        if (allocated(A)) deallocate(A)
        allocate(A(n, n))
    
        do i = 1, n
            do j = 1, n
                call random_number(A(i, j))
            end do
        end do
    end subroutine random_matrix

    ! Function to calculate the identity matrix
    function identity(m) result(identity_matrix)
        integer, intent(in) :: m
        
        integer :: i, j
        real(8), allocatable :: identity_matrix(:,:)
        
        allocate(identity_matrix(m,m))
        identity_matrix = 0.0
        do i = 1, m
            identity_matrix(i,i) = 1.0
        end do
 
    end function identity

    




end program