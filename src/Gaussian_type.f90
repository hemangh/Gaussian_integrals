module Gaussian_type
  use precision, only : idp
  use constants, only : pi, factorials, ddfct
  implicit none

  type contrct_Gaussian
    integer                  :: num
    real(idp) , dimension(3) :: origin
    integer   , dimension(3) :: power
    real(idp) , allocatable  :: exps(:)
    real(idp) , allocatable  :: coefs(:) 
    real(idp) , allocatable  :: norm(:)
  contains
    procedure :: initialize
    procedure :: normalize
  end type
contains
  subroutine initialize(self, origin, power, exps, coefs)
  class(contrct_Gaussian) :: self

    real(idp) , dimension(3), intent(in) :: origin
    integer   , dimension(3), intent(in) :: power
    real(idp) , intent(in)  :: exps(:)
    real(idp) , intent(in)  :: coefs(:) 

    integer     :: num

    num = size(exps)
    !check size of inputs:
    if(num < 1) stop 'size of exps is not correct'
    if(num /= size(coefs)) stop 'size of coefs is not compatible with exps'

    self%num = num
    self%origin = origin
    self%power  = power
    allocate(self%exps, source = exps)
    allocate(self%coefs, source = coefs)
    allocate(self%norm(num))
    call self%normalize()
  end subroutine

  subroutine  normalize(self)
  class(contrct_Gaussian) :: self
    real(idp) :: norm
    real(idp) :: prefactor
    integer   :: l(3)
    integer   :: lmax
    integer   :: pow
    integer   :: dbl_fact_prod
    integer   :: ia, ib
    integer   :: num_exp

    call factorials()

    l(:) = self%power(:)
    lmax = sum (l)
    pow  = 2*lmax + 1.5
    dbl_fact_prod = ddfct(2*l(1) - 1) * ddfct(2*l(2) - 1) * ddfct(2*l(3) - 1)

    self%norm(:) = 2.d0 ** pow * & 
      self%exps(:) ** pow / dbl_fact_prod
    self%norm = sqrt(self%norm)

    prefactor = pi ** 1.5 * dbl_fact_prod / 2.d0 ** lmax

    norm = 0.d0
    num_exp = self%num
    do ia = 1, num_exp
      do ib = 1, num_exp
        norm = norm + self%norm(ia) * self%norm(ib) * self%coefs(ia) * self%coefs(ib) &
          / (self%coefs(ia) + self%coefs(ib) ) ** (lmax + 1.5)        
      end do
    end do

    norm = norm * prefactor
    norm = norm ** (-0.5)

    do ia = 1, num_exp
      self%coefs(ia) = self%coefs(ia) * norm
    end do

  end subroutine

end module Gaussian_type
