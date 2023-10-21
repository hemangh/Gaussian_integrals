module Gaussian_overlap
    use precision, only : idp, i4
    use constants, only : pi
    use Gaussian_type
    implicit none
    
    

contains
    
    recursive function Hermite_Gaussian_coefficients(l1,l2,num_nodes,distance,exp1,exp2) result(E)
    ! Recursive definition of Hermite Gaussian expansion coefficients.
    ! Returns a real double.
    ! l1,l2: orbital angular momentum number on Gaussian '2' and '1'
    ! num_nodes: number nodes in Hermite (depends on type of integral,
    ! e.g. always zero for overlap integrals)
    ! distance: distance between origins of Gaussian '1' and '2'
    ! exp1: orbital exponent on Gaussian '1' (e.g. alpha in the text)
    ! exp2: orbital exponent on Gaussian '2' (e.g. beta in the text)
    
    integer(i4) :: l1,l2
    integer(i4) :: num_nodes
    real(idp)   :: distance
    real(idp)   :: exp1, exp2

    real(idp) :: E
    real(idp) :: p, q
    
    p = exp1 + exp2
    q = exp1*exp2/p
    if ((num_nodes < 0) .or. (num_nodes > (l1 + l2))) then
    ! out of bounds for num_nodes
        E = 0.d0
    elseif (l1 == l2 .and. l2 == num_nodes .and. num_nodes == 0) then
    ! base case
        E = exp(-q*distance*distance) ! K_AB
    elseif (l2 == 0) then
    ! ! decrement index l1
        E = &
            (1.d0/(2.d0*p))  * Hermite_Gaussian_coefficients(l1-1,l2,num_nodes-1,distance,exp1,exp2) - &
            (q*distance/exp1)* Hermite_Gaussian_coefficients(l1-1,l2,num_nodes,distance,exp1,exp2) + &
            (num_nodes+1)*Hermite_Gaussian_coefficients(l1-1,l2,num_nodes+1,distance,exp1,exp2)
    else
    ! decrement index l2
        E = &
             (1/(2*p))*Hermite_Gaussian_coefficients(l1,l2-1,num_nodes-1,distance,exp1,exp2)   + &
            (q*distance/exp2)*Hermite_Gaussian_coefficients(l1,l2-1,num_nodes,distance,exp1,exp2)+ &
            (num_nodes+1)*Hermite_Gaussian_coefficients(l1,l2-1,num_nodes+1,distance,exp1,exp2)
    endif

    end function Hermite_Gaussian_coefficients


    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    !! Evaluation of overlap of two primitive Gaussians on output
    !! Input : a, b (real=8) :: exponents of the two Gaussians
    !! lmn1,lmn2 (integer)   :: array dimensioned 3, powers or orbital angular momentum of
    !!                       :: the two Gaussians in x, y, z
    !! A, B (real=8)         :: array dimensioned 3, coordinates of the origins of the two 
    !!                          Gaussians
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!! 
    function Analytic_Gaussian_Overlap_Integral (a, lmn1, coord_A, b, lmn2, coord_B) result (overlap)

        real(idp) :: a, b, coord_A(3), coord_B(3)
        integer   :: lmn1(3), lmn2(3)
        real(idp) :: overlap
        
        real(idp) :: S(3)
        integer   :: i
        
        do i = 1, 3

            S(i) = Hermite_Gaussian_coefficients (lmn1(i), lmn2(i),0, abs(coord_A(i)-coord_B(i)), a, b)

        end do    
        
        overlap = product(S) * (sqrt (pi/(a + b)))** 3

        end function
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    !! Evaluates overlap between two contracted Gaussians
    !!    Output (real=8) :: overlap
    !!    Input: a, b (type(contrct_Gaussian))
    !!    a: contracted Gaussian 'a', contrct_Gaussian object
    !!    b: contracted Gaussian 'b', contrct_Gaussian object
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    function contracted_Gaussian_overlap(a,b) result(overlap)
    
        type(contrct_Gaussian) :: a, b
        real(idp) :: overlap

        integer :: ia, ib
        overlap = 0.d0
        do ia = 1, a%num
            do ib = 1, b%num
                overlap = overlap + a%norm(ia) * b%norm(ib) * a%coefs(ia) * b%coefs(ib) &
                                  * Analytic_Gaussian_Overlap_Integral(a%exps(ia), a%power,a%origin, b%exps(ib), b%power,b%origin)
            end do
        end do

    end function

end module Gaussian_overlap