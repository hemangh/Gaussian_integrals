module Gaussian_overlap
    use precision, only : idp, i4
    use constants, only : pi
    use contracted_Gaussian_type, only: contrct_Gaussian

    implicit none
    
    

contains
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    ! Recursive definition of Hermite Gaussian expansion coefficients.
    ! Returns a real double.
    ! l1,l2: orbital angular momentum number on Gaussian '2' and '1'
    ! num_nodes: number nodes in Hermite (depends on type of integral,
    ! e.g. always zero for overlap integrals)
    ! distance: distance between origins of Gaussian '1' and '2'
    ! exp1: orbital exponent on Gaussian '1' (e.g. alpha in the text)
    ! exp2: orbital exponent on Gaussian '2' (e.g. beta in the text)
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
    recursive function Hermite_Gaussian_coefficients(l1,l2,num_nodes,distance,exp1,exp2) result(E)
    
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
             (1.d0/(2.d0*p))*Hermite_Gaussian_coefficients(l1,l2-1,num_nodes-1,distance,exp1,exp2)   + &
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
    function primitive_Gaussian_Overlap_Integral (a, lmn1, coord_A, b, lmn2, coord_B) result (overlap)

        real(idp) :: a, b, coord_A(3), coord_B(3)
        integer   :: lmn1(3), lmn2(3)
        real(idp) :: overlap
        
        real(idp) :: S(3)
        integer   :: i
        
        do i = 1, 3

            S(i) = Hermite_Gaussian_coefficients (lmn1(i), lmn2(i),0, coord_A(i)-coord_B(i), a, b)

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
        real(idp) :: overlap, sab

        integer :: ia, ib
        overlap = 0.d0
        do ia = 1, a%num
            do ib = 1, b%num
                sab = primitive_Gaussian_Overlap_Integral(a%exps(ia), a%power,a%origin, b%exps(ib), b%power,b%origin)
                overlap = overlap + a%norm(ia) * b%norm(ib) * a%coefs(ia) * b%coefs(ib) &
                * sab
                !* primitive_Gaussian_Overlap_Integral(a%exps(ia), a%power,a%origin, b%exps(ib), b%power,b%origin)
                ! *OV(a%exps(ia), a%power,a%origin, b%exps(ib), b%power,b%origin)
            end do
        end do

    end function

    function OV(alpha, LA, RA, beta, LB, RB) result(Overlap)
        real(kind=8), intent(in) :: alpha, beta
        real(kind=8), dimension(3), intent(in) :: RA, RB
        integer, dimension(3), intent(in) :: LA, LB
        real(kind=8) :: Overlap, EAB
        real(kind=8), dimension(:,:,:) , allocatable :: s
        integer :: i, a, b, max_LA, max_LB
        
        ! Initialize s(i, a, 0)
        max_LA = max(LA(1),LA(2),LA(3))
        max_LB = max(LB(1),LB(2),LB(3))
        allocate(s(1:3, 0:max_LA+2, 0:max_LB+2))
        do i = 1, 3
          s(i, 0, 0) = 1.0d0
          s(i, 1, 0) = -(RA(i) - (alpha*RA(i) + beta*RB(i)) / (alpha + beta))
          do a = 2, LA(i)
            s(i, a, 0) = -(RA(i) - (alpha*RA(i) + beta*RB(i)) / (alpha + beta)) * s(i, a - 1, 0) &
                         + ((a - 1) / (2.0d0 * (alpha + beta))) * s(i, a - 2, 0)
          end do
        end do
        
        ! Transfer equation
        do i = 1, 3
          do b = 0, LB(i)
            do a = 0, LA(i)
              s(i, a, b) = s(i, a + 1, b - 1) + (RA(i) - RB(i)) * s(i, a, b - 1)
            end do
          end do
        end do
        
        ! Compute EAB and overlap
        EAB = exp(-alpha * beta / (alpha + beta) * sum((RA - RB) ** 2))
        Overlap = EAB * (pi / (alpha + beta)) ** 1.5d0 * s(1, LA(1), LB(1)) * s(2, LA(2), LB(2)) * s(3, LA(3), LB(3))
      end function OV

end module Gaussian_overlap