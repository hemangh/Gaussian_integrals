program test
    use Gaussian_overlap
    use contracted_Gaussian_type, only : contrct_Gaussian
    ! use basisfunctions, only: contrct_Gaussian
    implicit none
    
    real(8) :: alpha1, alpha2
    integer, dimension(3) :: pow1, pow2
    real(8), dimension(3) :: coord_1, coord_2
    real(8), allocatable  :: exps(:), coeff(:)
    type(contrct_Gaussian) :: a, b

    real(8) :: analytic_ans, norm_contr
    integer :: i
    
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
    allocate(exps(3), coeff(3))
    exps = [130.709320, 23.808861, 6.44360830] 
    coeff= [0.15432897, 0.53532814, 0.44463454]
    pow1 = [0, 0, 0]
    coord_1 = [0.000000000000 , -0.143225816552 ,  0.000000000000]
    
    call a%initialize(coord_1, pow1, exps, coeff) 
    !! check each term for normalization:
    do i = 1, size(exps)
        analytic_ans = primitive_Gaussian_Overlap_Integral (exps(i),pow1,coord_1, exps(i),pow1,coord_1)
        print*, "primitive ",i,":",analytic_ans * a%norm(i) ** 2
    end do

    norm_contr = contracted_Gaussian_overlap (a,a)  
    print*, "contracted basis norm:", norm_contr 
    close(7)
    
end program test