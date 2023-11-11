program test_primitivesN2
    use Gaussian_overlap
    use nuclear_attrraction_integral
    use precision, only: idp
    use reader_mod
    implicit none
    


    type(contrct_Gaussian) , allocatable :: basis(:,:)
    integer :: i, j, basis_size, num_atoms
    character(len=:), allocatable :: symb
    real(8) :: analytic_ans, s(5)
    real(8), allocatable ::  mat_pot(:,:)

    open(unit=7, file="compN2_log.txt")

    call read_basis_format2 ("../../input/inpsmallN2.txt", basis)

    do i = 5, 7
            s(1) = contracted_Gaussian_overlap (basis(i,1),basis(i,1))
            s(2) = contracted_Gaussian_overlap(basis(i,2),basis(i,2))
            basis(i,1)%norm(:) = 1.d0/sqrt(s(1) / 3.d0) * basis(i,1)%norm(:)
            basis(i,2)%norm(:) = 1.d0/sqrt(s(2) / 3.d0) * basis(i,2)%norm(:)
            !print*, contracted_Gaussian_overlap(basis(i,2),basis(i,2))
    end do 

    do i = 11, 13
        s(1) = contracted_Gaussian_overlap (basis(i,1),basis(i,1))
        s(2) = contracted_Gaussian_overlap(basis(i,2),basis(i,2))
        basis(i,1)%norm(:) = 1.d0/sqrt(s(1)/15.d0 ) * basis(i,1)%norm(:)
        basis(i,2)%norm(:) = 1.d0/sqrt(s(2)/15.d0 ) * basis(i,2)%norm(:)
        print*, contracted_Gaussian_overlap(basis(i,1),basis(i,1))
        print*, contracted_Gaussian_overlap(basis(i,2),basis(i,2))
    end do 

    do i = 14, 19
        s(1) = contracted_Gaussian_overlap (basis(i,1),basis(i,1))
        s(2) = contracted_Gaussian_overlap(basis(i,2),basis(i,2))
        basis(i,1)%norm(:) = 1.d0/sqrt(s(1)/3.d0 ) * basis(i,1)%norm(:)
        basis(i,2)%norm(:) = 1.d0/sqrt(s(2)/3.d0 ) * basis(i,2)%norm(:)
        print*, contracted_Gaussian_overlap(basis(i,1),basis(i,1))
        print*, contracted_Gaussian_overlap(basis(i,2),basis(i,2))
    end do 

    print*, contracted_Gaussian_overlap(basis(i,1),basis(i,1))
    print*, contracted_Gaussian_overlap(basis(i,2),basis(i,2))

    write(7,*)  "                           Overlap Integrals                               "
    write(7, *) "++++++++++++++++++++++++++++over atom 1+++++++++++++++++++++++++++++++++++"
    write(7,*) "        s1-h1           x1-h1          y1-h1        z1-h1         xx1-h1"
    do i = 1, 5
        do j = 1, i
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,1))
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    write(7, *) "++++++++++++++++++++++++++++over atom 2+++++++++++++++++++++++++++++++++++"
    write(7,*) "        s1-h1           x1-h1          y1-h1        z1-h1         xx1-h1"
    do i = 1, 5
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do

    write(7, *) "++++++++++++++++++++++++++++over atom 1+++++++++++++++++++++++++++++++++++"
    write(7,*) "        yy1-h1           zz1-h1          xy1-h1        xz1-h1         yz1-h1"
    do i = 1, 5
        do j = 1, i
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,1))
        end do
        symb = pow_to_symb(basis(i+5, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,1))
        end do
        symb = pow_to_symb(basis(i+5, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i+5,1))
        end do
        symb = pow_to_symb(basis(i+5, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    write(7,*)
    write(7, *) "++++++++++++++++++++++++++++over atom 2+++++++++++++++++++++++++++++++++++"
    write(7,*) "        yy1-h1           zz1-h1          xy1-h1        xz1-h1         yz1-h1"
    do i = 1, 5
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = contracted_Gaussian_overlap (basis(j+5,1),basis(i,2))
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do

    print*, -7.d0* (contracted_Gaussian_nuclear_attraction(basis(1,1), basis(1,1),basis(1,1)%origin) + &
    contracted_Gaussian_nuclear_attraction(basis(1,1), basis(1,1),basis(1,2)%origin))
    
    num_atoms = size(basis,2)
    basis_size = size(basis, 1)
    
    allocate(mat_pot(basis_size * num_atoms,basis_size * num_atoms))
    mat_pot = 0.d0
    
    do i = 1, basis_size
        do j = i, basis_size
            mat_pot(i,j) = -7.d0* (contracted_Gaussian_nuclear_attraction(basis(i,1), basis(j,1),basis(1,1)%origin) + &
            contracted_Gaussian_nuclear_attraction(basis(i,1), basis(j,1),basis(1,2)%origin))
            mat_pot(j,i) = mat_pot(i,j)
        end do
    end do
    
    do i =  1, basis_size
        do j = 1, basis_size
            mat_pot(i,basis_size + j) = -7.d0* &
            (contracted_Gaussian_nuclear_attraction(basis(i,1), basis(j,2),basis(1,1)%origin) + &
            contracted_Gaussian_nuclear_attraction(basis(i,1), basis(j,2),basis(1,2)%origin))
            mat_pot(basis_size +j,i) = - mat_pot(i,basis_size +j)
        end do
    end do
    print*, -7.d0* &
    (contracted_Gaussian_nuclear_attraction(basis(5,1), basis(6,2),basis(1,1)%origin) + &
    contracted_Gaussian_nuclear_attraction(basis(5,1), basis(6,2),basis(1,2)%origin))
    ! call matrix_print(mat_pot)

    write(7,*)  "                 Couloumb  Potential energy Integrals                               "
    write(7, *) "++++++++++++++++++++++++++++over atom 1+++++++++++++++++++++++++++++++++++"
    write(7,*) "        s1-h1           x1-h1          y1-h1        z1-h1         xx1-h1"
    do i = 1, 5
        do j = 1, i
            s(j) = mat_pot(i,j)
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = mat_pot(i,j)
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = mat_pot(i,j)
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = mat_pot(i,j)
        end do
        symb = pow_to_symb(basis(i, 1)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do

    write(7, *) "++++++++++++++++++++++++++++over atom 2+++++++++++++++++++++++++++++++++++"
    write(7,*) "        s1-h1           x1-h1          y1-h1        z1-h1         xx1-h1"
    do i = 1, 5
        do j = 1, 5
            s(j) = mat_pot(i,j+basis_size)
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 6, 10
        do j = 1, 5
            s(j) = mat_pot(i,j+basis_size)
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 11, 15
        do j = 1, 5
            s(j) = mat_pot(i,j+basis_size)
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
    do i = 16, 20
        do j = 1, 5
            s(j) = mat_pot(i,j+basis_size)
        end do
        symb = pow_to_symb(basis(i, 2)%power)
        write(7,'(A3,x,5(E15.8))') symb,s(:j-1)
        deallocate(symb)
    end do
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
        case(3)
            if(pow(1) == 3) symb = 'xxx'
            if(pow(2) == 3) symb = 'yyy'
            if(pow(3) == 3) symb = 'zzz'
            if(pow(1) == 2 .and. pow(2) == 1) symb = 'xxy'
            if(pow(1) == 2 .and. pow(3) == 1) symb = 'xxz'
            if(pow(2) == 2.and. pow(3) == 1) symb = 'yyz'
            if(pow(2) == 2.and. pow(1) == 1) symb = 'yyx'
            if(pow(2) == 1.and. pow(3) == 2) symb = 'yzz'
            if(pow(1) == 1.and. pow(3) == 2) symb = 'xzz'
            if(pow(1) == 1.and. pow(2) == 1 .and. pow(3) == 1) symb = 'xyz'
        case default
            symb = '-'
        
        end select

    end function 

    subroutine matrix_print(mat, width)
        real(8), intent(in) :: mat(:,:)
        integer, intent(in) :: width
        integer :: i, j
        integer :: mat_size
        

        mat_size = size(mat, 1)
        print'(A)', "  "
        do i = 1, mat_size/ 2
            do j = i, width
                !print'(A3, )',
            end do
        end do    

        end subroutine
end program test_primitivesN2
