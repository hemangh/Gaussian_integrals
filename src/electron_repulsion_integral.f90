module electron_repulsion_integral
    use precision, only : idp
    use constants, only : pi
    use Gaussian_overlap, only : Hermite_Gaussian_coefficients
    use nuclear_attrraction_integral, only : Couloumb_Hermite_integral
    use contracted_Gaussian_type, only : contrct_Gaussian
    implicit none
    
contains
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
! Evaluates electron repulsion integral between four primitive Gaussians
! Returns a real(kind=8) result.
! exp_a,exp_b,exp_c,exp_d:
! orbital exponent on Gaussian 'a','b','c','d'
! pow_a,pow_b
! pow_c,pow_d: int 3-D arrays containing orbital angular momentum
! for Gaussian 'a','b','c','d', respectively
! origin_a,origin_b,origin_c,origin_d:
! 3-D arrays containing origin of Gaussian 'a','b','c','d'
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
function primitive_Gaussian_electron_repulsion(exp_a, pow_a, origin_a, exp_b, pow_b, origin_b, &
                                   exp_c, pow_c, origin_c, exp_d, pow_d, origin_d) result(result)
    
    real(kind=8) :: exp_a, exp_b, exp_c, exp_d
    integer, dimension(3) :: pow_a, pow_b, pow_c, pow_d
    real(kind=8), dimension(3) :: origin_a(3), origin_b(3), origin_c(3), origin_d(3)
    real(kind=8) :: result
    
    real(kind=8) :: gamma, lambda, alpha, RPQ, val
    integer :: t, u, v, tau, nu, phi
    real(kind=8) :: HermiteCoeff1, HermiteCoeff2, HermiteCoeff3, HermiteCoeff4, HermiteCoeff5, HermiteCoeff6
    real(kind=8) :: P(3), Q(3)
    
! Extracting components from pow_a, pow_b, pow_c, and pow_d
    integer :: l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4
    l1 = pow_a(1)
    m1 = pow_a(2)
    n1 = pow_a(3)
    l2 = pow_b(1)
    m2 = pow_b(2)
    n2 = pow_b(3)
    l3 = pow_c(1)
    m3 = pow_c(2)
    n3 = pow_c(3)
    l4 = pow_d(1)
    m4 = pow_d(2)
    n4 = pow_d(3)
    
    gamma = exp_a + exp_b
    lambda = exp_c + exp_d
    alpha = (gamma * lambda) / (gamma + lambda)
    
    ! Calculate Gaussian product center
    P = 0.5d0 * ((exp_a * origin_a) + (exp_b * origin_b)) / (exp_a + exp_b)
    Q = 0.5d0 * ((exp_c * origin_c) + (exp_d * origin_d)) / (exp_c + exp_d)
    RPQ = sqrt((P(1) - Q(1))**2 + (P(2) - Q(2))**2 + (P(3) - Q(3))**2)
    
    val = 0.0
    
    do t = 0, l1 + l2
      do u = 0, m1 + m2
        do v = 0, n1 + n2
          do tau = 0, l3 + l4
            do nu = 0, m3 + m4
              do phi = 0, n3 + n4
                HermiteCoeff1 = Hermite_Gaussian_coefficients(l1, l2, t, origin_a(1) - origin_b(1), exp_a, exp_b)
                HermiteCoeff2 = Hermite_Gaussian_coefficients(m1, m2, u, origin_a(2) - origin_b(2), exp_a, exp_b)
                HermiteCoeff3 = Hermite_Gaussian_coefficients(n1, n2, v, origin_a(3) - origin_b(3), exp_a, exp_b)
                HermiteCoeff4 = Hermite_Gaussian_coefficients(l3, l4, tau, origin_c(1) - origin_d(1), exp_c, exp_d)
                HermiteCoeff5 = Hermite_Gaussian_coefficients(m3, m4, nu, origin_c(2) - origin_d(2), exp_c, exp_d)
                HermiteCoeff6 = Hermite_Gaussian_coefficients(n3, n4, phi, origin_c(3) - origin_d(3), exp_c, exp_d)
                
                val = val + HermiteCoeff1 * HermiteCoeff2 * HermiteCoeff3 * HermiteCoeff4 * HermiteCoeff5 * HermiteCoeff6 * &
                      (-1.0d0)**(tau + nu + phi) * &
                       Couloumb_Hermite_integral(t + tau, u + nu, v + phi, 0, alpha, P(1) - Q(1), P(2) - Q(2), P(3) - Q(3), RPQ)
              end do
            end do
          end do
        end do
      end do
    end do
    
    result = val * 2.0d0 * (pi ** 2.5d0) / (gamma * lambda * sqrt(gamma + lambda))
    
  end function primitive_Gaussian_electron_repulsion
  
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
  !! Evaluates electron repulsion between four contracted Gaussians
  !!    Output (real=8) :: integral
  !!    Input: a, b (type(contrct_Gaussian))
  !!    a, b, c, d: contracted Gaussian 'a', 'b', 'c', 'd'.
  !!    contrct_Gaussian object
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
  function contracted_Gaussian_electron_repulsion_integral (a,b,c,d) result(integral)

    type(contrct_Gaussian) :: a, b, c, d
    real(idp) :: integral

    integer :: ia, ib, ic, id

    integral = 0.d0

    do ia = 1, a%num
        do ib = 1, b%num
            do ic = 1, c%num
                do id = 1, d%num
            integral = integral + a%norm(ia)* b%norm(ib) * a%coefs(ia) * b%coefs(ib) * &
            primitive_Gaussian_electron_repulsion (a%exps(ia), a%power, a%origin, b%exps(ib), b%power, b%origin , &
                                                   c%exps(ic), c%power, c%origin, d%exps(id), d%power, d%origin)
                end do
            end do
        end do
    end do

end function contracted_Gaussian_electron_repulsion_integral

    
end module electron_repulsion_integral