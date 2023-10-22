module nuclear_attrraction_integral
    use precision, only : idp
    use constants, only : pi
    use Gaussian_type
    use Gaussian_overlap, only : Hermite_Gaussian_coefficients
    use boys_function, only: boys
    implicit none
    
contains

recursive function Couloumb_Hermite_integral (t,u,v,n,p,PCx,PCy,PCz,RPC) result(val)
    integer, intent(in) :: t, u, v, n
    real(idp), intent(in) :: p, PCx, PCy, PCz, RPC
    real(idp) :: val
    
    real(idp) :: x

    x = p * RPC * RPC
    val = 0.0
    if (t == 0 .and. u == 0 .and. v == 0) then
        val = (-2.0d0 * p) ** n * boys(n, x)
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
subroutine primitive_Gaussian_nuclear_attraction(exp_a, pow_a, origin_a, exp_b, pow_b, origin_b, C, result)
    
    real(idp), intent(in) :: exp_a, exp_b
    integer, intent(in) :: pow_a(3), pow_b(3)
    real(idp), intent(in) :: origin_a(3), origin_b(3), C(3)
    real(idp), intent(out) :: result
  
    integer :: l1, m1, n1, l2, m2, n2, t, u, v
    real(idp) :: gamma, RPC, val, P(3)
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
          RVal = Couloumb_Hermite_integral(t, u, v, 0, gamma, P(1) - C(1), P(2) - C(2), P(3) - C(3), RPC)
          val = val + E1 * E2 * E3 * RVal
        end do
      end do
    end do
  
    val = val * 2 * pi / gamma
    result = val
  end subroutine primitive_Gaussian_nuclear_attraction
  

end module nuclear_attrraction_integral
