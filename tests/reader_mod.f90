module reader_mod
    use precision, only : idp
    use contracted_Gaussian_type, only: contrct_Gaussian
    use atom_type, only: atom
    implicit none
    

contains
    subroutine read_basis_format1(input_file, basis)
        character(len=*) , intent(in) :: input_file
        type(contrct_Gaussian), intent(out) , allocatable :: basis(:,:)
        type(atom) , allocatable:: atoms(:)
        integer :: num_atoms , num_orb,line_num,atm_num, i, j, k
        integer , dimension(3):: pow1, pow2
        real(idp) , allocatable:: coord(:,:)
        character(len=1) :: orb_char
        real(idp) :: exponent, coeff
        
        
        logical :: inputfile_exists
        !Check the file exists:    
        inquire(file=trim(input_file),exist=inputfile_exists)
        if(.not. inputfile_exists) then
            print*, "input file" // input_file //" is missing!"
            stop
        else 
        open(unit=8,file=input_file)
        endif
        
        read(8,*) num_atoms
        allocate(coord(3,num_atoms))
        allocate(atoms(num_atoms))
        do i = 1, num_atoms
            read(8,*) coord(1,i), coord(2, i), coord(3,i)
            print*, coord(:,i)
        end do
        
        read(8,*) line_num, num_orb
        print* , line_num, num_orb
        allocate(basis(num_orb, num_atoms))
        j = 0
        do i = 1, line_num

            read(8, *) orb_char
            read(8, *) exponent, coeff
            print*, orb_char, exponent, coeff

            if(orb_char == 's')then
                pow1 = [0,0,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
            end if

            if(orb_char == 'p')then
                do k = 1, 3
                    pow1 = [0,0,0]
                    pow1(k) = 1
                    j = j + 1
                    print*, j
                    call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                    call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                end do
            end if

            if(orb_char == 'd')then
                !xx, yy, zz
                do k = 1, 3
                    pow1 = [0,0,0]
                    pow1(k) = 2
                    j = j + 1
                    print*, j
                    call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                    call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                end do
                !xy
                pow1 = [1,1,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !xz
                pow1 = [1,0,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !yz
                pow1 = [0,1,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])


            end if
        
            if(orb_char == 'f')then
                !xxx, yyy, zzz
                do k = 1, 3
                    pow1 = [0,0,0]
                    pow1(k) = 3
                    j = j + 1
                    print*, j
                    call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                    call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                end do
                !xyy
                pow1 = [1,2,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !xxy
                pow1 = [2,1,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !xxz
                pow1 = [2,0,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !xzz
                pow1 = [1,0,2]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !yzz
                pow1 = [0,1,2]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !yyz
                pow1 = [0,2,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
                !xyz
                pow1 = [1,1,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, [exponent], [coeff])
                call basis(j,2)%initialize(coord(:,2), pow1, [exponent], [coeff])
            end if

        end do

        close(8)


    end subroutine


    subroutine read_basis_format2(input_file, basis)
        character(len=*) , intent(in) :: input_file
        type(contrct_Gaussian), intent(out) , allocatable :: basis(:,:)
        type(atom) , allocatable:: atoms(:)
        integer :: num_atoms , num_orb,line_num,atm_num, i, j, k,l, num_exp
        integer , dimension(3):: pow1, pow2
        real(idp) , allocatable:: coord(:,:)
        character(len=1) :: orb_char
        real(idp), allocatable :: exponent(:), coeff(:)
        
        
        logical :: inputfile_exists
        !Check the file exists:    
        inquire(file=trim(input_file),exist=inputfile_exists)
        if(.not. inputfile_exists) then
            print*, "input file" // input_file //" is missing!"
            stop
        else 
        open(unit=8,file=input_file)
        endif
        ! read number of atoms
        read(8,*) num_atoms
        allocate(coord(3,num_atoms))
        allocate(atoms(num_atoms))
        !read atoms origin coordinates
        do i = 1, num_atoms
            read(8,*) coord(1,i), coord(2, i), coord(3,i)
            print*, coord(:,i)
        end do
        
        read(8,*) line_num, num_orb
        print* , line_num, num_orb
        allocate(basis(num_orb, num_atoms))
        j = 0
        do i = 1, line_num

            read(8, *) orb_char, num_exp
            allocate(exponent(num_exp), coeff(num_exp))
            do l = 1, num_exp
                read(8, *) exponent(l), coeff(l)
            end do
            print*, orb_char, num_exp, exponent, coeff

            if(orb_char == 's')then
                pow1 = [0,0,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
            end if

            if(orb_char == 'p')then
                do k = 1, 3
                    pow1 = [0,0,0]
                    pow1(k) = 1
                    j = j + 1
                    print*, j
                    call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                    call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                end do
            end if

            if(orb_char == 'd')then
                !xx, yy, zz
                do k = 1, 3
                    pow1 = [0,0,0]
                    pow1(k) = 2
                    j = j + 1
                    print*, j
                    call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                    call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                end do
                !xy
                pow1 = [1,1,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !xz
                pow1 = [1,0,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !yz
                pow1 = [0,1,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)


            end if
        
            if(orb_char == 'f')then
                !xxx, yyy, zzz
                do k = 1, 3
                    pow1 = [0,0,0]
                    pow1(k) = 3
                    j = j + 1
                    print*, j
                    call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                    call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                end do
                !xyy
                pow1 = [1,2,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !xxy
                pow1 = [2,1,0]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !xxz
                pow1 = [2,0,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !xzz
                pow1 = [1,0,2]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !yzz
                pow1 = [0,1,2]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !yyz
                pow1 = [0,2,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
                !xyz
                pow1 = [1,1,1]
                j = j + 1
                print*, j
                call basis(j,1)%initialize(coord(:,1), pow1, exponent, coeff)
                call basis(j,2)%initialize(coord(:,2), pow1, exponent, coeff)
            end if
            deallocate(exponent, coeff)
        end do

        close(8)


    end subroutine read_basis_format2

end module reader_mod