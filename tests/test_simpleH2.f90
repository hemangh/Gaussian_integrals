program test_simpleH2
    use Gaussian_overlap!, only: primitive_Gaussian_Overlap_Integral
    use kinetic_integral, only: primitive_Gaussians_kinetic_integral
    use nuclear_attrraction_integral, only: primitive_Gaussian_nuclear_attraction
    use electron_repulsion_integral, only: primitive_Gaussian_electron_repulsion
    use reader_mod

    implicit none
    real(8) :: alpha1, alpha2
    real(8) :: expon(1), coef(1)
    integer, dimension(3) :: pow1, pow2
    real(8), dimension(3) :: coord_1, coord_2
    real(8) :: analytic_ans, s(5)
    type(contrct_Gaussian) , allocatable :: basis(:,:)
    integer :: i, j
    character(len=:), allocatable :: symb

    
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


    open(unit=7, file="compH2_log.txt")

        call read_basis_format1 ("../../input/inpsmall.txt", basis)

        do i = 9, 11
             s(1) = contracted_Gaussian_overlap (basis(i,1),basis(i,1))
             s(2) = contracted_Gaussian_overlap(basis(i,2),basis(i,2))
             basis(i,1)%norm(:) = 1.d0/sqrt(s(1) / 3.d0) * basis(i,1)%norm(:)
             basis(i,2)%norm(:) = 1.d0/sqrt(s(2) / 3.d0) * basis(i,2)%norm(:)
             print*, contracted_Gaussian_overlap(basis(i,2),basis(i,2))
        end do 

        do i = 15, 17
            s(1) = contracted_Gaussian_overlap (basis(i,1),basis(i,1))
            s(2) = contracted_Gaussian_overlap(basis(i,2),basis(i,2))
            basis(i,1)%norm(:) = 1.d0/sqrt(s(1) / 3.d0) * basis(i,1)%norm(:)
            basis(i,2)%norm(:) = 1.d0/sqrt(s(2) / 3.d0) * basis(i,2)%norm(:)
            print*, contracted_Gaussian_overlap(basis(i,2),basis(i,2))
       end do 
       write(7, *) "++++++++++++++++++++++++++++over atom 1+++++++++++++++++++++++++++++++++++"
        write(7,*) "        s1-h1           s2-h1          x1-h1        y1-h1         z1-h1"
        do i = 1, 5
            do j = 1, i
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
            end do
            symb = pow_to_symb(basis(i, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 6, 10
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
            end do
            symb = pow_to_symb(basis(i, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 11, 15
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
            end do
            symb = pow_to_symb(basis(i, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 16, 20
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
            end do
            symb = pow_to_symb(basis(i, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        write(7, *) "++++++++++++++++++++++++++++over atom 2+++++++++++++++++++++++++++++++++++"
        write(7,*) "        s1-h1           s2-h1          x1-h1        y1-h1         z1-h1"
        do i = 1, 5
            do j = 1, i
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
            end do
            symb = pow_to_symb(basis(i, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 6, 10
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
            end do
            symb = pow_to_symb(basis(i, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 11, 15
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
            end do
            symb = pow_to_symb(basis(i, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 16, 20
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
            end do
            symb = pow_to_symb(basis(i, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do

        write(7, *) "++++++++++++++++++++++++++++over atom 1+++++++++++++++++++++++++++++++++++"
        write(7,*) "        x2-h1           y2-h1          z2-h1        xx1-h1         yy1-h1"
        do i = 1, 5
            do j = 1, i
                s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,1))
            end do
            symb = pow_to_symb(basis(i+5, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 6, 10
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,1))
            end do
            symb = pow_to_symb(basis(i+5, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 11, 15
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,1))
            end do
            symb = pow_to_symb(basis(i+5, 1)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do

        write(7, *) "++++++++++++++++++++++++++++over atom 2+++++++++++++++++++++++++++++++++++"
        write(7,*) "        x2-h1           y2-h1          z2-h1        xx1-h1         yy1-h1"
        do i = 1, 5
            do j = 1, i
                s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,2))
            end do
            symb = pow_to_symb(basis(i+5, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 6, 10
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,2))
            end do
            symb = pow_to_symb(basis(i+5, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do
        do i = 11, 15
            do j = 1, 5
                s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,2))
            end do
            symb = pow_to_symb(basis(i+5, 2)%power)
            write(7,'(A2,x,5(E15.8))') symb,s(:j-1)
            deallocate(symb)
        end do

        ! do i = 21, 25
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 26, 30
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 31, 35
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f1(i))/ sqrt(0.2d0)
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 36, 37
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f1(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
    
        ! do i = 1, 5
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 6, 10
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 11, 15
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 16, 20
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i)) / sqrt(0.2d0)
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 21, 25
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 26, 30
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 31, 35
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do
        ! do i = 36, 37
        !     do j = 1, 5
        !         s(j) = contracted_Gaussian_overlap (f1(j),f2(i))
        !     end do
        !     write(7,'(I0,x,5(E15.8))') i,s(:j-1)
        ! end do




    close(7)

    contains

    function pow_to_symb(pow) result(symb)

        integer , dimension(3) :: pow
        character(len=:), allocatable :: symb
        
        integer :: lmax

        lmax = sum(pow)

        select case(lmax)
        case(0)  
            symb = 's'
        case(1)
            if(pow(1)== 1) symb = 'x'
            if(pow(2)== 1) symb = 'y'
            if(pow(3)== 1) symb = 'z'
        case(2)
            if(pow(1) == 2) symb = 'xx'
            if(pow(2) == 2) symb = 'yy'
            if(pow(3) == 2) symb = 'zz'
            if(pow(1) == 1 .and. pow(2) == 1) symb = 'xy'
            if(pow(1) == 1 .and. pow(3) == 1) symb = 'xz'
            if(pow(2) == 1 .and. pow(3) == 1) symb = 'yz'
        case default
            symb = '-'
        
        end select

    end function 


end program


