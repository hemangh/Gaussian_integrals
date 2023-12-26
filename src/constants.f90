module constants
  use precision, only : idp
  implicit none

  integer, private, parameter  :: maxfac = 30 
  real(idp), parameter  :: pi =  atan(1.0_idp) * 4.0_idp 
  real(idp), dimension(0:maxfac) :: factorial, double_factorial
  real(idp), dimension(-1:maxfac) :: factorial_2n_minus_1

contains

  subroutine factorials
    integer:: i
    logical:: called=.false.
    !implicit integer (a-z)
    !real *8 factorial, factorial_2n_minus_1
    !dimension factorial(0:maxfac), factorial_2n_minus_1(0:maxfac)

    if(called)then
      !write(iout,*) "*****factorial is already called******"
      return
    else
      write(*,*) "first time factorials is called"
      !----------------------------------------------------------------------c
      !               calculate factorials & double factorials               c
      !----------------------------------------------------------------------c

      factorial(0)=1.d+00
      factorial(1)=1.d+00
      double_factorial(0) = 1.d+00
      double_factorial(1) = 1.d+00
      if (maxfac.gt.1) then
        do 10 i=2,maxfac
          factorial(i)=i*factorial(i-1)
          double_factorial(i) = i*double_factorial(i-2)
          10     continue
        endif
        !----------------------------------------------------------------------c
        !           calculate (2*m-1) double factorial                         c
        !----------------------------------------------------------------------c
        factorial_2n_minus_1(-1) = 1.d+00
        factorial_2n_minus_1(0)=1.d+00
        factorial_2n_minus_1(1)=1.d+00
        factorial_2n_minus_1(2)=3.d+00
        if (maxfac.gt.2) then
          do 20 i=3,maxfac
            factorial_2n_minus_1(i)=(i+i-1)*factorial_2n_minus_1(i-1)
            20    continue
          endif
          called = .true.
          return
        end if
      end subroutine factorials

      function binomial_coeff(n, k) result(bico)
        integer :: n, k
        real(idp) :: bico
        bico = exp(factorial(n)-factorial(k)-factorial(n-k))
      end function

    end module constants
