module kinetic_integral
    use precision, only: idp
    use Gaussian_type
    use Gaussian_overlap
    implicit none
    
contains
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !! Evaluates kinetic energy integral between two Gaussians
! Output: a real=8.
! Inputs:
! exp_a:
! real(8) orbital exponent on Gaussian 'a' (e.g. alpha in the text)
! exp_b:
! real(8) orbital exponent on Gaussian 'b' (e.g. beta in the text)
! pow_a: int array(3) containing orbital angular momentum (e.g. (1,0,0))
! for Gaussian 'a'
! pow_b: int array(3) containing orbital angular momentum for Gaussian 'b'
! origin_a:
! real(8) array containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
! origin_b:
! real(8) array containing origin of Gaussian 'b'
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function primitive_Gaussians_kinetic_integral (exp_a, pow_a, origin_a, exp_b, pow_b, origin_b ) result(integral)
        
        real(idp) :: exp_a, exp_b
        real(idp), dimension(3) :: origin_a, origin_b
        integer  , dimension(3) :: pow_a, pow_b     
        real(idp) :: integral

        real(idp) :: overlap0
        real(idp) :: overlap1_p, overlap2_p, overlap3_p
        real(idp) :: overlap1_m, overlap2_m, overlap3_m
        real(idp) :: term0, term1, term2
        integer   :: lmax_a, lmax_b
        
        lmax_a = sum (pow_a)
        lmax_b = sum (pow_b)
        overlap0 = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,pow_b,origin_b)
        overlap1_p = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,[pow_b(1)+2,pow_b(2), pow_b(3) ],origin_b)
        overlap2_p = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,[pow_b(1),pow_b(2)+2, pow_b(3) ],origin_b)
        overlap3_p = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,[pow_b(1),pow_b(2), pow_b(3)+2 ],origin_b)

        overlap1_m = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,[pow_b(1)-2,pow_b(2), pow_b(3) ],origin_b)
        overlap2_m = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,[pow_b(1),pow_b(2)-2, pow_b(3) ],origin_b)
        overlap3_m = primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,[pow_b(1),pow_b(2), pow_b(3)-2 ],origin_b)

        term0 = exp_b*(2.d0*lmax_b+3)*primitive_Gaussian_Overlap_Integral(exp_a,pow_a,origin_a,exp_b,pow_b,origin_b)
        term1 = -2.d0*exp_b**2*(overlap1_p + overlap2_p + overlap3_p)
        term2 = -0.5d0*(pow_b(1)*(pow_b(1)-1)*overlap1_m + pow_b(2)*(pow_b(2)-1)*overlap2_m + pow_b(3)*(pow_b(3)-1)*overlap3_m)

        integral = term0 + term1 + term2

    end function
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    ! '''Evaluates kinetic energy between two contracted Gaussians
    ! Output: real(8) integral.
    ! Inputs:
    ! a: contracted Gaussian 'a', contrct_Gaussian type
    ! b: contracted Gaussian 'b', contrct_Gaussian type
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    function contracted_Gaussian_kinetic_integral (a,b) result(integral)

        type(contrct_Gaussian) :: a, b
        real(idp) :: integral

        integer :: ia, ib

        integral = 0.d0

        do ia = 1, a%num
            do ib = 1, b%num
                integral = integral + a%norm(ia)* b%norm(ib) * a%coefs(ia) * b%coefs(ib) * &
                primitive_Gaussians_kinetic_integral(a%exps(ia), a%power, a%origin, b%exps(ib), b%power, b%origin)
            end do
        end do

    end function

end module kinetic_integral
    


