program test_simpleH2
    use Gaussian_overlap, only: primitive_Gaussian_Overlap_Integral
    use kinetic_integral, only: primitive_Gaussians_kinetic_integral
    use nuclear_attrraction_integral, only: primitive_Gaussian_nuclear_attraction
    use electron_repulsion_integral, only: primitive_Gaussian_electron_repulsion
    implicit none
    real(8) :: alpha1, alpha2
    integer, dimension(3) :: pow1, pow2
    real(8), dimension(3) :: coord_1, coord_2
    real(8) :: analytic_ans
    
    open(unit=7, file="simpleH2test_log.txt")

    alpha1 = 5.d-1
    alpha2 = 5.d-1

    pow1 = [0, 0 , 0]
    pow2 = [0, 0,  0]

    coord_1 = [0.d0, 0.d0, 0.5d0]
    coord_2 = [0.d0, 0.d0,-0.5d0]
    analytic_ans = primitive_Gaussian_Overlap_Integral (alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,*) "Analytical overlap integrals for 2 primitive Gaussian functions"
    write(7,'(A,Es22.15E2)' ) "s and s primitives,  1 bohr apart alpha = 1:", analytic_ans
    
    analytic_ans = primitive_Gaussians_kinetic_integral(alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,*) "Analytical kinetic integrals for 2 primitive Gaussian functions"
    write(7,'(A,Es22.15E2)' ) "s and s primitives,  1 bohr apart alpha = 1:", analytic_ans

    analytic_ans = primitive_Gaussian_nuclear_attraction(alpha1, pow1, coord_1, alpha2, pow2, coord_2, coord_1)
    write(7,*) "Analytical nuclear attraction integrals for 2 primitive Gaussian functions"
    write(7,'(A,Es22.15E2)' ) "s and s primitives,  1 bohr apart alpha = 1, nuc=1:", analytic_ans

    analytic_ans = primitive_Gaussian_nuclear_attraction(alpha1, pow1, coord_1, alpha2, pow2, coord_2, coord_2)
    write(7,*) "Analytical nuclear attraction integrals for 2 primitive Gaussian functions"
    write(7,'(A,Es22.15E2)' ) "s and s primitives,  1 bohr apart alpha = 1, nuc = 2:", analytic_ans

    analytic_ans = primitive_Gaussian_electron_repulsion(alpha1, pow1, coord_1, alpha2, pow2, coord_2, &
                                                         alpha1, pow1, coord_1, alpha2, pow2, coord_2)
    write(7,*) "Analytical electron repulsion integrals for 2 primitive Gaussian functions"
    write(7,'(A,Es22.15E2)' ) "s and s primitives,  1 bohr apart alpha = 1:", analytic_ans

    close(7)

end program