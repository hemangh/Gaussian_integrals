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
    procedure :: normalize
  end type
contains
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
