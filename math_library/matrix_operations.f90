module matrix_operations
    use precision, only: idp
    implicit none
    
contains

    
      ! -- Returns the inverse of a general square matrix A 
      function invs(A) result(Ainv)
        
        real(idp),intent(in) :: A(:,:)
        real(idp), allocatable :: Ainv(:,:)!, Acopy(:,:), C(:,:)
        real(idp), allocatable :: work(:) !size(A,2))            ! work array for LAPACK
        integer              :: m,n,info, lwork, i
        integer,allocatable  ::ipiv(:)     ! pivot indices
        real(idp), allocatable :: s(:), u(:, :), vt(:, :)
        real(idp), allocatable :: s_dagger(:), ut(:, :), v(:, :)
    
        
        m = size(A,1) ! rows
        n = size(A,2) ! columns
        print*, m,'x',n
        allocate(Ainv(m,n))
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A


        if (m /= n) then
          print*, "none square matrix does not have an inverse. Returning psudo-inverse if possible."
          !     Compute the singular values and left and right singular vectors
          !     of A = u * s * vt
          !lda = m
          !ldu = m
          !ldvt = n
          
          lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))
          
          allocate(work(lwork))
          allocate(s(min(m,n)))
          allocate(u(m,m))
          allocate(vt(n,n))
          Call dgesvd('A', 'S', m, n, Ainv, m, s, u, m, vt, n, work, lwork, info)

          if (info/=0) Then
            Write (*, '(1X, A, I4)') 'Failure in DGESVD. INFO =', info
            return
          end If
          deallocate(Ainv) ! it is overwritten by dgesvd
          allocate(Ainv(n,m))
          ! Need At = v * s_dagger * ut
          allocate(ut, source=u)
          allocate(v, source =vt)
          allocate(s_dagger, source=s)
          ut = transpose(u)
          v  = transpose(vt)
          s_dagger(:) = 1.d0 / s(:)
          print*, "U:"
          call print_matrix(u)
          print*, "S:"
          print*, s
          print*, "VT:"
          call print_matrix(vt)
          Ainv = 0.d0
          do i = 1, size(s_dagger)
            Ainv(i,i)= s_dagger(i)
          end do
          ! allocate(Acopy, source=A)
          ! Acopy = 0.d0
          ! do i = 1, size(s)
          !   Acopy(i,i) = s(i)
          ! end do

          ! print*, "S_dagger"
          ! call print_matrix(Ainv)

          ! print*, "S"
          ! call print_matrix(Acopy)
          ! print*, "UT"
          ! call print_matrix(ut)
          

          ! print*, "Test"
          ! call print_matrix(matmul(u,ut))
          ! call print_matrix(matmul(vt,v))
          ! call print_matrix(matmul(Acopy,Ainv))
          ! call print_matrix(matmul(Ainv,Acopy))


          Ainv = matmul(Ainv, ut)
          ! Acopy = matmul(u,Acopy)
          ! call print_matrix(matmul(Acopy,Ainv))
          ! allocate(C(size(Acopy,1),size(Acopy,1)))
          ! call dgemm('N','N', size(Acopy,1), size(Ainv,2), size(Acopy,2), &
          !  1.d0, Acopy,size(Acopy,1) ,Ainv, size(Acopy,2),0.d0, C, size(Acopy,1) )
          ! call print_matrix(C)
          ! call print_matrix(matmul(Ainv, Acopy))
          ! print*, "S_dagger * UT"
          ! call print_matrix(Ainv)
          Ainv = matmul(v, Ainv) 

          return
        end if

        allocate(ipiv(n))
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n,n,Ainv,n,ipiv,info)
        if (info.ne.0) then
          print*, 'Matrix is numerically singular!'
          return
        end if
        allocate(work(size(A,2)))
        ! SGETRI computes the inverse of a matrix using the LU factorization
        ! computed by SGETRF.
        call DGETRI(n,Ainv,n,ipiv,work,n,info)
        if (info.ne.0) then
          print*, 'Matrix inversion failed!'
          return
        end if
      end function invs

      function matrix_transpose(A) result(At)

        real(idp) :: A(:,:)
        real(idp), allocatable :: At(:,:)

        allocate(At(size(A,1), size(A,2)))

        

      end function
    
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

end module matrix_operations