module molden_parser
    use precision, only: idp
    use contracted_Gaussian_type, only: contrct_Gaussian
    use string_array_mod
    implicit none
    
contains
    subroutine read_molden_file (moldeninp,atm_contrct_Gaussian)
        character(len=*)                   , intent(in)  :: moldeninp
        type(contrct_Gaussian), allocatable, intent(out) :: atm_contrct_Gaussian(:)
        
        integer, parameter   :: max_num = 20

        integer :: iostat, linesum, indx, count
        integer :: i,j,k, number_atoms, atm_number, num_primitives
        integer :: saved_num_prim, saved_type
        integer, allocatable :: atmnum(:), charge(:)
        

        real(idp),allocatable :: xyzcoord(:,:), expn(:), coeff(:)
        real(idp),allocatable :: expn_saved(:), saved_coeff(:,:)
        real(idp)             :: energy
        real(idp)             :: occupancy
        
        character(len=15)::aline
        character(len=7)::atoms_key
        character(len=5) :: gto_key
        character(len=4) :: mo_key
        character(len=1)::bracket, orb_symb
        character(len=6)::dummy
        character(len=4) :: symkey, enkey
        character(len=5) :: spinkey, spin
        character(len=6) :: occkey, symm
        character(len=1) :: saved_symb

        type(string_array_type),allocatable::symb(:)

        open(unit=18, file=moldeninp, status='old', iostat=iostat)
        if ( iostat /= 0 ) stop 'molden file is not opened'

        read(18,*)
        read(18,'(A7)') atoms_key
        if (atoms_key .ne. '[Atoms]') stop 'molden file [Atoms] section is missing'
        !determine how many lines to the next section
        number_atoms = 0
        do 
           read(18,'(A1)') bracket
           if ( bracket .eq. '[' ) exit
           number_atoms = number_atoms + 1
        end do

        do i = 1, number_atoms+1
            backspace (18)
        end do
        
        allocate(symb(number_atoms))
        allocate(atmnum(number_atoms),charge(number_atoms)) 
        allocate(xyzcoord(3,number_atoms))
        ! number_atoms here is really the number of atoms
        do  i = 1, number_atoms   
            read(18,*) dummy, atmnum(i), charge(i), xyzcoord(:,i)
            symb(i)%string = trim(dummy)
            print*, symb(i)%string, atmnum(i), charge(i), xyzcoord(:,i)
        end do
        
        !initialize number of atoms numer type with orbitals
        allocate(atm_contrct_Gaussian(number_atoms))

        read(18,*) gto_key
        if (gto_key .ne. '[GTO]') stop 'molden file [GTO] section is missing'
        ! determine how many lines are in GTO section
        linesum = 1
        do i = 1, number_atoms
            k = 0
            read(18,*) atm_number
            if (atm_number /= i) stop 'reading GTOS ran into a problem'

            saved_num_prim = -1
            saved_symb     = 'x'
            saved_type     = 0
            expn_saved     = [0.d0]            

            do 
                linesum = linesum + 1
                read(18,'(A)')aline
                if(index(aline,'               ')>0) then 
                    exit
                else
                   backspace(18) 
                   read(18,*) orb_symb, num_primitives
                   print*, orb_symb, num_primitives
                   k = k + 1

                   if (orb_symb /= saved_symb .or. num_primitives /= saved_num_prim) then
                        
                        if(allocated(expn)) deallocate(expn)
                        allocate(expn(num_primitives))
                        saved_num_prim = num_primitives
                        saved_symb  = orb_symb
                        

                   endif           

                   do j = 1, num_primitives
                      read(18,*) expn(j)
                   end do
                   if (size(expn_saved) == num_primitives)then

                    if(all(expn_saved .eq. expn) ) then
                       
                    else
                        
                        expn_saved = expn
                        saved_type  = saved_type + 1

                    endif

                   else  
                 
                    saved_type  = saved_type + 1
                    deallocate(expn_saved)
                    expn_saved = expn

                   endif      

                   linesum = linesum + num_primitives
                
                end if


            end do 
            
            print*, i , k, linesum
            ! Initialize atom i with saved_type number of grouped orbitals
            call atm_contrct_Gaussian(i)%initialize(saved_type,xyzcoord(:,i))

            !print*, atm(i)%center
        end do
        linesum = linesum + 1
        print*, linesum
        do i = 1, linesum
            backspace(18)
        end do
        ! read(18,'(A)')aline
        ! print*, aline    
        do i = 1, number_atoms
            k = 0
            count = 0
            read(18,*) atm_number
            if (atm_number /= i) stop 'reading GTOS ran into a problem'
            
            saved_num_prim = -1
            saved_symb     = 'x'
            saved_type     = 0
            expn_saved     = [0.d0]
            do 
                read(18,'(A)')aline
                if(index(aline,'               ')>0) then 
                    exit
                else
                   backspace(18) 
                   read(18,*) orb_symb, num_primitives
                   
                   if (orb_symb /= saved_symb .or. num_primitives /= saved_num_prim) then
                       if (count /= 0) then
                        call atm(i)%initialize_gto_data_per_atom(saved_type,saved_num_prim,trim(saved_symb), count)
                        call atm(i)%input_gauss_data_per_orbital_type(saved_type,saved_num_prim,count, expn_saved,&
                                                                      saved_coeff(1:saved_num_prim,1:count))
                       endif
                       if(allocated(expn)) deallocate(expn)
                       if(allocated(coeff)) deallocate(coeff)
                       allocate(expn(num_primitives),coeff(num_primitives))
                       if (allocated(saved_coeff))then
                          deallocate(saved_coeff)
                       endif   
                       allocate(saved_coeff(num_primitives, max_num))
                       saved_num_prim = num_primitives
                       saved_symb  = orb_symb

                       count = 0
                   endif  

                   k = k + 1
                   !call atm(i)%gto(k)%initialize(num_primitives,trim(orb_symb))
                   !call atm(i)%initialize_gto_data_per_atom(k,num_primitives, trim(orb_symb))
                   
                   do j = 1, num_primitives
                      read(18,*) expn(j),coeff(j)
                      !atm(i)%gto(k)%expon(j), atm(i)%gto(k)%cont_coeff(j)
                   end do
 
                   if (size(expn_saved) == num_primitives)then

                    if(all(expn_saved .ne. expn) ) then
                       saved_type  = saved_type + 1
                       deallocate(expn_saved)
                       expn_saved = expn

                    else
                        
                        count = count + 1
                        saved_coeff(:,count) = coeff(:)

                    endif

                   else

                        saved_type  = saved_type + 1
                        deallocate(expn_saved)
                        expn_saved = expn
                        count = count + 1
                        saved_coeff(:,count) = coeff(:)

                   endif
                   
                   
                   !call atm(i)%input_gauss_data_per_orbital(k,num_primitives,expn,coeff,saved_type)
                
                end if


            end do 
            
            call atm(i)%initialize_gto_data_per_atom(saved_type,saved_num_prim,trim(saved_symb), count)
            call atm(i)%input_gauss_data_per_orbital_type(saved_type,saved_num_prim,count, expn_saved, &
                                                          saved_coeff(1:saved_num_prim,1:count))
        end do 

        
    end subroutine read_molden_file

end module molden_parser