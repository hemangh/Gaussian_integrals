module nuclear_attrraction_integral
    use precision, only : idp
    use constants, only : pi, factorial, double_factorial, factorials
    use contracted_Gaussian_type
    use Gaussian_overlap, only : Hermite_Gaussian_coefficients
    use boys_function, only: boys
    implicit none
    
contains
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
! Returns the Coulomb Hermite integrals (R^0_{tuv})
! Returns a real(kind=8) val.
! Arguments:
! t,u,v:
! order of Coulomb Hermite derivative in x,y,z
! (see defs in Helgaker and Taylor)
! p:
! sum of exponents of two Gaussians
! PCx,y,z: Cartesian vector distance between Gaussian
! composite center P and nuclear center C
! RPC:
! Distance between P and C
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
function Couloumb_Hermite_integral_t0 (t,u,v,p,PCx,PCy,PCz,RPC) result(val)
    integer, intent(in) ::  t, u, v
    real(idp), intent(in) :: p, PCx, PCy, PCz, RPC
    real(idp) :: val

    real(idp), external :: Boys_func
    
    real(idp) :: x, prefac, b
    integer   :: i, j, k, order

    call factorials
    x = p * RPC * RPC
    val = 0.d0
    do i = 0, t/2
        do j = 0, u/2
            do k = 0, v/2
                prefac = factorial(t) * factorial(u) * factorial(v) &
                /double_factorial(2*i)/double_factorial(2*j)/double_factorial(2*k)/ &
                factorial(t-2*i)/factorial(u-2*j)/factorial(v-2*k)
                order = t+u+v-i-j-k
                b = boys(order,x)
                val = val +  PCx ** (t-2*i)*PCy**(u-2*j)*PCz**(v-2*k) * &
                prefac * b *(-2.d0*p) ** order
            end do
        end do
    end do

end function

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
! Returns the Coulomb auxiliary Hermite integrals
! Returns a real(kind=8) val.
! Arguments:
! t,u,v:
! order of Coulomb Hermite derivative in x,y,z
! (see defs in Helgaker and Taylor)
! p:
! sum of exponents of two Gaussians
! n:
! order of Boys function
! PCx,y,z: Cartesian vector distance between Gaussian
! composite center P and nuclear center C
! RPC:
! Distance between P and C
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
recursive function Couloumb_Hermite_integral (t,u,v,n,p,PCx,PCy,PCz,RPC) result(val)
    integer, intent(in) :: t, u, v, n
    real(idp), intent(in) :: p, PCx, PCy, PCz, RPC
    real(idp) :: val

    real(idp), external :: Boys_func
    
    real(idp) :: x

    x = p * RPC * RPC
    val = 0.0
    if (t == 0 .and. u == 0 .and. v == 0) then
        val = (-2.0d0 * p) ** n * Boys_func(n, x)
    elseif (t == 0 .and. u == 0) then
        if (v > 1) then
        val = (v - 1) * Couloumb_Hermite_integral(t, u, v - 2, n + 1, p, PCx, PCy, PCz, RPC)
        end if
        val = val + PCz * Couloumb_Hermite_integral(t, u, v - 1, n + 1, p, PCx, PCy, PCz, RPC)
    elseif (t == 0) then
        if (u > 1) then
        val = (u - 1) * Couloumb_Hermite_integral(t, u - 2, v, n + 1, p, PCx, PCy, PCz, RPC)
        end if
        val = val + PCy * Couloumb_Hermite_integral(t, u - 1, v, n + 1, p, PCx, PCy, PCz, RPC)
    else
        if (t > 1) then
        val = (t - 1) * Couloumb_Hermite_integral(t - 2, u, v, n + 1, p, PCx, PCy, PCz, RPC)
        end if
        val = val + PCx * Couloumb_Hermite_integral(t - 1, u, v, n + 1, p, PCx, PCy, PCz, RPC)
    end if
end function
    
function gaussian_product_center(exp_a, origin_a, exp_b, origin_b) result(product_cen)
    real(idp) :: exp_a, origin_a(3), exp_b, origin_b(3), product_cen(3)
    product_cen = (exp_a*origin_a+exp_b*origin_b)/(exp_a+exp_b)
end function

!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
! Evaluates kinetic energy integral between two Gaussians
! Output: result, real(8) integral result.
! Input:
! exp_a:
! orbital exponent on Gaussian 'a' (e.g. alpha in the text)
! exp_b:
! orbital exponent on Gaussian 'b' (e.g. beta in the text)
! pow_a: int 3-d array containing orbital angular momentum (e.g. (1,0,0))
! for Gaussian 'a'
! pow_b: int 3-d array containing orbital angular momentum for Gaussian 'b'
! origin_a:
! 3-d array containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
! origin_b:
! 3-d array containing origin of Gaussian 'b'
! C:
! 3-d array containing origin of nuclear center 'C'
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
function primitive_Gaussian_nuclear_attraction(exp_a, pow_a, origin_a, exp_b, pow_b, origin_b, C) result(val)
    
    real(idp) :: exp_a, exp_b
    integer   :: pow_a(3), pow_b(3)
    real(idp) :: origin_a(3), origin_b(3), C(3)
    real(idp) :: val
  
    integer :: l1, m1, n1, l2, m2, n2, t, u, v
    real(idp) :: gamma, RPC, P(3)
    real(idp) :: E1, E2, E3, RVal
  
    l1 = pow_a(1)
    m1 = pow_a(2)
    n1 = pow_a(3)
    l2 = pow_b(1)
    m2 = pow_b(2)
    n2 = pow_b(3)
  
    gamma = exp_a + exp_b
    P = gaussian_product_center(exp_a, origin_a, exp_b, origin_b)
    RPC = sqrt((P(1)-C(1))**2 + (P(2)-C(2))**2 + (P(3)-C(3))**2)
    val = 0.0
  
    do t = 0, l1 + l2
      do u = 0, m1 + m2
        do v = 0, n1 + n2
          E1 = Hermite_Gaussian_coefficients(l1, l2, t, origin_a(1) - origin_b(1), exp_a, exp_b)
          E2 = Hermite_Gaussian_coefficients(m1, m2, u, origin_a(2) - origin_b(2), exp_a, exp_b)
          E3 = Hermite_Gaussian_coefficients(n1, n2, v, origin_a(3) - origin_b(3), exp_a, exp_b)
        !  RVal = Couloumb_Hermite_integral(t, u, v, 0, gamma, P(1) - C(1), P(2) - C(2), P(3) - C(3), RPC)
          RVal = Couloumb_Hermite_integral_t0( t, u, v, gamma, P(1) - C(1), P(2) - C(2), P(3) - C(3), RPC)
          val = val + E1 * E2 * E3 * RVal
        end do
      end do
    end do
  
    val = val * 2 * pi / gamma
    
end function primitive_Gaussian_nuclear_attraction
  
  !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
  !! Evaluates nuclear attraction between two contracted Gaussians
  !!    Output (real=8) :: integral
  !!    Input: a, b (type(contrct_Gaussian))
  !!    a: contracted Gaussian 'a', contrct_Gaussian object
  !!    b: contracted Gaussian 'b', contrct_Gaussian object
  !!    C: 3-d (real=8) array, center of nucleus
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
  function contracted_Gaussian_nuclear_attraction(a,b,C) result(integral)
    
    type(contrct_Gaussian) :: a, b
    real(idp) :: C(3)
    real(idp) :: integral

    integer :: ia, ib
    integral = 0.d0
    do ia = 1, a%num
        do ib = 1, b%num
            integral = integral + a%norm(ia) * b%norm(ib) * a%coefs(ia) * b%coefs(ib) &
                    * primitive_Gaussian_nuclear_attraction (a%exps(ia), a%power,a%origin, b%exps(ib), b%power,b%origin, C)
        end do
    end do

end function

end module nuclear_attrraction_integral
