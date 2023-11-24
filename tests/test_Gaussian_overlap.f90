program test
    use Gaussian_overlap
    use contracted_Gaussian_type, only : contrct_Gaussian
    ! use basisfunctions, only: contrct_Gaussian
    implicit none
    
    real(8) :: alpha1, alpha2, s(5)
    integer, dimension(3) :: pow1, pow2
    real(8), dimension(3) :: coord_1, coord_2
    real(8), allocatable  :: exps(:), coeff(:), HGcoeff(:,:,:)
    type(contrct_Gaussian) :: a, b
    type(contrct_Gaussian) :: g(7)
    type(contrct_Gaussian) :: f1(37), f2(37)

    real(8) :: analytic_ans, norm_contr
    integer :: i, j
    
    open(unit=7, file="test_log.txt")

    alpha1 = 1.d0
    alpha2 = 1.d0

    pow1 = [0, 0 , 0]
    pow2 = [0, 0,  0]

    coord_1 = [0.d0, 0.d0, 0.5d0]
    coord_2 = [0.d0, 0.d0,-0.5d0]
    
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,*) "Analytical overlap integrals for various primitive Gaussian functions"
    write(7,'(A,Es22.15E2)' ) "s and s primitives,  1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 0 , 0]
    pow2 = [1, 0,  0]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "s and px primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 0 , 0]
    pow2 = [0, 1,  0]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "s and py primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 0 , 0]
    pow2 = [0, 0,  1]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "s and pz primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [1, 0 , 0]
    pow2 = [0, 0,  1]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "px and pz primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 1 , 0]
    pow2 = [0, 0,  1]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "py and pz primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 0 , 1]
    pow2 = [0, 0,  1]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "pz and pz primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 1 , 0]
    pow2 = [0, 1,  0]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "py and py primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [1, 0 , 0]
    pow2 = [1, 0,  0]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "px and px primitives, 1 bohr apart alpha = 1:", analytic_ans

    pow1 = [0, 0 , 1]
    pow2 = [0, 0,  1]
    coord_1 = [0.d0, 0.d0, 0.05d0]
    coord_2 = [0.d0, 0.d0,-0.05d0]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "pz and pz primitives, 0.1 bohr apart, alpha = 1:", analytic_ans

    pow1 = [0, 1 , 0]
    pow2 = [0, 1,  0]
    
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "py and py primitives, 0.1 bohr apart, alpha = 1:", analytic_ans

    pow1 = [1, 0 , 0]
    pow2 = [1, 0,  0]
    
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,'(A,Es22.15E2)' ) "px and px primitives, 0.1 bohr apart, alpha = 1:", analytic_ans

    !! Specific tests:
    !! s & s overlap, H2O tests:
    !! orb. 1
    print*, "orbital 1: O 1s"
    allocate(exps(3), coeff(3))
    exps = [130.709320, 23.808861, 6.44360830] 
    coeff= [0.15432897, 0.53532814, 0.44463454]
    pow1 = [0, 0, 0]
    coord_1 = [0.000000000000 , -0.143225816552 ,  0.000000000000]
    
    call g(1)%initialize(coord_1, pow1, exps, coeff) 
    !! check each term for normalization:
    do i = 1, size(exps)
        analytic_ans = primitive_Gaussian_Overlap_Integral (exps(i),pow1,coord_1, exps(i),pow1,coord_1)
        print*, "primitive ",i,":",analytic_ans * g(1)%norm(i) ** 2
    end do

    norm_contr = contracted_Gaussian_overlap (g(1),g(1))  
    print*, "contracted basis norm:", norm_contr 
    ! orb. 2
    print*, "orbital 2: O 2s"
    exps = [5.0331513_idp, 1.1695961_idp, 0.3803890_idp] 
    coeff= [-0.09996723_idp, 0.39951283_idp, 0.70011547_idp]
    pow1 = [0, 0, 0]
    coord_1 = [0.000000000000_idp , -0.143225816552_idp ,  0.000000000000_idp]
     
    call g(2)%initialize(coord_1, pow1, exps, coeff)  

    analytic_ans = contracted_Gaussian_overlap (g(1),g(2))   
    print*, "overlap 1&2:", analytic_ans
    print*, "diff:", 0.236703936510848d0 - analytic_ans
    !orb. 3
    print*, "orbital 3: O 2px" 
    pow2 = [1,0,0]
    
    call g(3)%initialize(coord_1, pow2, exps, coeff)
    analytic_ans = contracted_Gaussian_overlap (g(1),g(3))   
    print*, "overlap 1&3:", analytic_ans
    !orb. 4
    print*, "orbital 4: O 2py"
    pow2 = [0,1,0]

    call g(4)%initialize(coord_1, pow2, exps, coeff)
    analytic_ans = contracted_Gaussian_overlap (g(1),g(4))   
    print*, "overlap 1&4:", analytic_ans
    !orb. 5
    print*, "orbital 5: O 2pz"
    pow2 = [0,0,1]
    
    call g(5)%initialize(coord_1, pow2, exps, coeff)
    analytic_ans = contracted_Gaussian_overlap (g(1),g(5))   
    print*, "overlap 1&5:", analytic_ans
    
    !orb. 6
    print*, "orbital 6: H1 1s"
    exps = [3.42525091d0, 0.62391373d0, 0.16885540d0] 
    coeff= [0.15432897d0, 0.53532814d0, 0.44463454d0]
    pow1 = [0, 0, 0]
    coord_2 = [1.638036840407d0 ,  1.136548822547d0 , -0.000000000000d0]
    
    
    call g(6)%initialize(coord_2, pow1, exps, coeff)
    analytic_ans = contracted_Gaussian_overlap (g(1),g(6))
    print*, "overlap 1&6:", analytic_ans
    !orb. 7
    print*, "orbital 7: H2 1s"
    coord_2 = [-1.638036840407d0 ,  1.136548822547d0 , -0.000000000000d0]
    
    call g(7)%initialize(coord_2, pow1, exps, coeff)
    analytic_ans = contracted_Gaussian_overlap (g(1),g(7))
    print*, "overlap 1&7:", analytic_ans
    print*, "diff:", 0.038405599785757d0 - analytic_ans
 
    ! !orb. 2
    ! exps = [5.0331513d0, 1.1695961d0, 0.3803890d0] 
    ! coeff= [-0.09996723d0, 0.39951283d0, 0.70011547d0]
    ! pow1 = [0, 0, 0]
    ! coord_1 = [0.000000000000d0 , -0.143225816552d0 ,  0.000000000000d0]
    ! call a%destroy() 
    ! call a%initialize(coord_1, pow1, exps, coeff)  

    analytic_ans = contracted_Gaussian_overlap (g(2),g(7))   
    print*, "overlap 7&2:", analytic_ans
    print*, "diff:", 0.386138840478250d0 - analytic_ans

    ! !orb. 3
    ! pow2 = [1,0,0]
    ! call a%destroy()
    ! call a%initialize(coord_1, pow2, exps, coeff)
    norm_contr = contracted_Gaussian_overlap (g(3),g(3))  
    print*, "contracted basis norm:", norm_contr 
    analytic_ans = contracted_Gaussian_overlap (g(3),g(7))   
    print*, "overlap 7&3:", analytic_ans
    print*, "diff:", -0.268438243716457d0 - analytic_ans
    ! !orb. 4
    ! pow2 = [0,1,0]
    ! call a%destroy()
    ! call a%initialize(coord_1, pow2, exps, coeff)
    norm_contr = contracted_Gaussian_overlap (g(4),g(4))  
    print*, "contracted basis norm:", norm_contr 
    analytic_ans = contracted_Gaussian_overlap (g(4),g(7))   
    print*, "overlap 7&4:", analytic_ans
    analytic_ans = contracted_Gaussian_overlap (g(7),g(4))   
    print*, "overlap 4&7:", analytic_ans
    print*, "diff:", 0.209726941420375d0 - analytic_ans
    ! !orb. 5
    ! pow2 = [0,0,1]
    ! call a%destroy()
    ! call a%initialize(coord_1, pow2, exps, coeff)
    norm_contr = contracted_Gaussian_overlap (g(5),g(5))  
    print*, "contracted basis norm:", norm_contr 
    analytic_ans = contracted_Gaussian_overlap (g(5),g(7))   
    print*, "overlap 7&5:", analytic_ans
    print*, "diff:", -0.000000000000000d0 - analytic_ans

    do i = 1, 7
        do j = 1, i
            write(7,'(x,2(I0,2x), F17.15)') j, i, contracted_Gaussian_overlap(g(j),g(i))
        end do
    end do


    print*, Hermite_Gaussian_coefficients(1,0,0,1.d0,1.d0,1.d0) / exp(-0.5d0)
    print*, Hermite_Gaussian_coefficients(0,1,0,1.d0,1.d0,1.d0) / exp(-0.5d0)

    print*, Hermite_Gaussian_coefficients(1,0,0,1.d0,0.25d0,0.75d0) / exp(-0.1875d0)
    
    print*, primitive_Gaussian_Overlap_Integral(1.d0,[1,0,0],[1.d0,0.d0,0.d0],&
    1.d0,[0,0,0],[0.d0,0.d0,0.d0])
    print*, primitive_Gaussian_Overlap_Integral(1.d0,[0,0,0],[1.d0,0.d0,0.d0],&
    1.d0,[1,0,0],[0.d0,0.d0,0.d0])

    coord_1 = [0.d0 ,0.d0,  .7d0]
    coord_2 = [0.d0 ,0.d0, -.7d0]
    s(1) = Hermite_Gaussian_coefficients(1,0,0,1.4d0,1.d3,1.d1)
    s(2) = Hermite_Gaussian_coefficients(0,0,0,0.d0,1.d3,1.d1)
    s(3) = Hermite_Gaussian_coefficients(0,0,0,0.d0,1.d3,1.d1)
    print*, s(1)
    print*, s(2)
    print*, s(3)
    print*, (sqrt (pi/(1.d3 + 1.d1)))** 3 * s(1) * s(2) * s(3)
    print*
    
    ! type=s
    ! 1000.00    1.00
    call f1(1)%initialize(coord_1,[0,0,0], [1.d3],[1.d0])
    call f2(1)%initialize(coord_2,[0,0,0], [1.d3],[1.d0])
    ! print'(F17.15)', contracted_Gaussian_overlap (f(1),f(1)) 
    !  type=s
    ! 100.00     1.00
    call f1(2)%initialize(coord_1,[0,0,0], [1.d2],[1.d0])
    call f2(2)%initialize(coord_2,[0,0,0], [1.d2],[1.d0])
    !  type=s
    ! 10.00      1.00
    call f1(3)%initialize(coord_1,[0,0,0], [1.d1],[1.d0])
    call f2(3)%initialize(coord_2,[0,0,0], [1.d1],[1.d0])
    !  type=s
    !  1.00      1.00
    call f1(4)%initialize(coord_1,[0,0,0], [1.d0],[1.d0])
    call f2(4)%initialize(coord_2,[0,0,0], [1.d0],[1.d0])
    !  type=s
    !   .100     1.00
    call f1(5)%initialize(coord_1,[0,0,0], [1.d-1],[1.d0])
    call f2(5)%initialize(coord_2,[0,0,0], [1.d-1],[1.d0])
    !  type=s
    !   .01      1.00  
    call f1(6)%initialize(coord_1,[0,0,0], [1.d-2],[1.d0])
    call f2(6)%initialize(coord_2,[0,0,0], [1.d-2],[1.d0])
    !  type=s
    !   .001     1.00
    call f1(7)%initialize(coord_1,[0,0,0], [1.d-3],[1.d0])
    call f2(7)%initialize(coord_2,[0,0,0], [1.d-3],[1.d0])
    !  type=p
    !   10.00    1.00
    call f1(8)%initialize(coord_1,[1,0,0], [1.d1],[1.d0])
    call f2(8)%initialize(coord_2,[1,0,0], [1.d1],[1.d0])
    call f1(9)%initialize(coord_1,[0,1,0], [1.d1],[1.d0])
    call f2(9)%initialize(coord_2,[0,1,0], [1.d1],[1.d0])
    call f1(10)%initialize(coord_1,[0,0,1], [1.d1],[1.d0])
    call f2(10)%initialize(coord_2,[0,0,1], [1.d1],[1.d0])
    !  type=p
    !   1.00     1.00

    call f1(11)%initialize(coord_1,[1,0,0], [1.d0],[1.d0])
    call f2(11)%initialize(coord_2,[1,0,0], [1.d0],[1.d0])
    call f1(12)%initialize(coord_1,[0,1,0], [1.d0],[1.d0])
    call f2(12)%initialize(coord_2,[0,1,0], [1.d0],[1.d0])
    call f1(13)%initialize(coord_1,[0,0,1], [1.d0],[1.d0])
    call f2(13)%initialize(coord_2,[0,0,1], [1.d0],[1.d0])
    !  type=p
    !    .10     1.00

    call f1(14)%initialize(coord_1,[1,0,0], [1.d-1],[1.d0])
    call f2(14)%initialize(coord_2,[1,0,0], [1.d-1],[1.d0])
    call f1(15)%initialize(coord_1,[0,1,0], [1.d-1],[1.d0])
    call f2(15)%initialize(coord_2,[0,1,0], [1.d-1],[1.d0])
    call f1(16)%initialize(coord_1,[0,0,1], [1.d-1],[1.d0])
    call f2(16)%initialize(coord_2,[0,0,1], [1.d-1],[1.d0])
    !  type=p
    !    .01     1.00
    call f1(17)%initialize(coord_1,[1,0,0], [1.d-2],[1.d0])
    call f2(17)%initialize(coord_2,[1,0,0], [1.d-2],[1.d0])
    call f1(18)%initialize(coord_1,[0,1,0], [1.d-2],[1.d0])
    call f2(18)%initialize(coord_2,[0,1,0], [1.d-2],[1.d0])
    call f1(19)%initialize(coord_1,[0,0,1], [1.d-2],[1.d0])
    call f2(19)%initialize(coord_2,[0,0,1], [1.d-2],[1.d0])
    !  type=d
    !   5.00     1.00
    call f1(20)%initialize(coord_1,[2,0,0], [5.d0],[1.d0])
    call f2(20)%initialize(coord_2,[2,0,0], [5.d0],[1.d0])
    call f1(21)%initialize(coord_1,[0,2,0], [5.d0],[1.d0])
    call f2(21)%initialize(coord_2,[0,2,0], [5.d0],[1.d0])
    call f1(22)%initialize(coord_1,[0,0,2], [5.d0],[1.d0])
    call f2(22)%initialize(coord_2,[0,0,2], [5.d0],[1.d0])
    call f1(23)%initialize(coord_1,[1,1,0], [5.d0],[1.d0])
    call f2(23)%initialize(coord_2,[1,1,0], [5.d0],[1.d0])
    call f1(24)%initialize(coord_1,[1,0,1], [5.d0],[1.d0])
    call f2(24)%initialize(coord_2,[1,0,1], [5.d0],[1.d0])
    call f1(25)%initialize(coord_1,[0,1,1], [5.d0],[1.d0])
    call f2(25)%initialize(coord_2,[0,1,1], [5.d0],[1.d0])
    !  type=d
    !   1.00     1.00
    call f1(26)%initialize(coord_1,[2,0,0], [1.d0],[1.d0])
    call f2(26)%initialize(coord_2,[2,0,0], [1.d0],[1.d0])
    call f1(27)%initialize(coord_1,[0,2,0], [1.d0],[1.d0])
    call f2(27)%initialize(coord_2,[0,2,0], [1.d0],[1.d0])
    call f1(28)%initialize(coord_1,[0,0,2], [1.d0],[1.d0])
    call f2(28)%initialize(coord_2,[0,0,2], [1.d0],[1.d0])
    call f1(29)%initialize(coord_1,[1,1,0], [1.d0],[1.d0])
    call f2(29)%initialize(coord_2,[1,1,0], [1.d0],[1.d0])
    call f1(30)%initialize(coord_1,[1,0,1], [1.d0],[1.d0])
    call f2(30)%initialize(coord_2,[1,0,1], [1.d0],[1.d0])
    call f1(31)%initialize(coord_1,[0,1,1], [1.d0],[1.d0])
    call f2(31)%initialize(coord_2,[0,1,1], [1.d0],[1.d0])
    !  type=d
    !    .10     1.00
    call f1(32)%initialize(coord_1,[2,0,0], [1.d-1],[1.d0])
    call f2(32)%initialize(coord_2,[2,0,0], [1.d-1],[1.d0])
    call f1(33)%initialize(coord_1,[0,2,0], [1.d-1],[1.d0])
    call f2(33)%initialize(coord_2,[0,2,0], [1.d-1],[1.d0])
    call f1(34)%initialize(coord_1,[0,0,2], [1.d-1],[1.d0])
    call f2(34)%initialize(coord_2,[0,0,2], [1.d-1],[1.d0])
    call f1(35)%initialize(coord_1,[1,1,0], [1.d-1],[1.d0])
    call f2(35)%initialize(coord_2,[1,1,0], [1.d-1],[1.d0])
    call f1(36)%initialize(coord_1,[1,0,1], [1.d-1],[1.d0])
    call f2(36)%initialize(coord_2,[1,0,1], [1.d-1],[1.d0])
    call f1(37)%initialize(coord_1,[0,1,1], [1.d-1],[1.d0])
    call f2(37)%initialize(coord_2,[0,1,1], [1.d-1],[1.d0])
    print'(7(Es15.8))', contracted_Gaussian_overlap (f1(1),f1(1)) ,&
    contracted_Gaussian_overlap (f1(1),f1(2)) ,&
    contracted_Gaussian_overlap (f1(1),f1(3)) ,&
    contracted_Gaussian_overlap (f1(1),f1(4)) ,&
    contracted_Gaussian_overlap (f1(1),f1(5)) ,&
    contracted_Gaussian_overlap (f1(1),f1(6)) ,&
    contracted_Gaussian_overlap (f1(1),f1(7)) !,&

    print'(7(Es15.8))', contracted_Gaussian_overlap (f1(1),f2(1)) ,&
    contracted_Gaussian_overlap (f1(1),f2(2)) ,&
    contracted_Gaussian_overlap (f1(1),f2(3)) ,&
    contracted_Gaussian_overlap (f1(1),f2(4)) ,&
    contracted_Gaussian_overlap (f1(1),f2(5)) ,&
    contracted_Gaussian_overlap (f1(1),f2(6)) ,&
    contracted_Gaussian_overlap (f1(1),f2(7)) !,&

    print'(7(E15.8))',contracted_Gaussian_overlap (f1(1),f2(10)) , &
    contracted_Gaussian_overlap (f1(2),f2(10)) ,&
    contracted_Gaussian_overlap (f1(3),f2(10)) ,&
    contracted_Gaussian_overlap (f1(4),f2(10)), &
    contracted_Gaussian_overlap (f1(5),f2(10)), &
    contracted_Gaussian_overlap (f1(6),f2(10)), &
    contracted_Gaussian_overlap (f1(7),f2(10))
    
   print*, contracted_Gaussian_overlap (f1(19),f1(19))
   print*, contracted_Gaussian_overlap (f1(32),f1(32)) 

    s(:) = 0.d0
    do i = 1, 5
        do j = 1, i
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i)) / sqrt(0.2d0/3.d0)
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 21, 25
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 26, 30
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 31, 35
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))/ sqrt(0.2d0)
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 36, 37
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do

    do i = 1, 5
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i)) / sqrt(0.2d0)
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 21, 25
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 26, 30
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 31, 35
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    do i = 36, 37
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        end do
        write(7,'(I0,x,5(E15.8))') i,s(:j-1)
    end do
    close(7)
    allocate(HGcoeff(0:1, 0:1, 0:1))
    call Hermite_Gaussian_coefficients_subroutine(1,1.d0,1.d0, 1.d0, HGcoeff)
    print*, HGcoeff
end program test