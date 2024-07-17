module integrals_3rspH_mod

    use, intrinsic :: iso_fortran_env, only: input_unit, error_unit 
    
    implicit none
    private

    public read_integrals_3rspH
    
    integer, parameter :: default_lmax = 20

    contains

    subroutine read_integrals_3rspH(integrals_3rspH, filename, lmax)
    
        ! Argument declarations with explicit intent and types
        real(kind=8), allocatable, intent(inout) :: integrals_3rspH(:,:,:,:,:,:)
        character(*), intent(in) :: filename  ! Adjust character length in calling program
        integer,      intent(out) :: lmax
    
        ! Internal variables
        integer :: num_lines = 0, line_indx, l1, l2, l3, m1, m2, m3
        integer :: ios
        integer, dimension(6) :: expected_shape
        logical :: exists
    
        ! File existence check
        inquire(file=filename, exist=exists)
        if (.not. exists) then
            write(error_unit,'(A)') "Error: File not found: " // trim(filename) 
            stop  ! Stop the program with an error message
        end if
    
        ! Open file in read-only mode
        open(unit=input_unit, file=filename, action="read", status="old", iostat=ios)
        if (ios /= 0) then  ! If ios is not zero, there was an error opening
            write(error_unit,'(A)') "Error opening file: " // trim(filename)
            stop  ! Stop the program with an error message
        end if
    
        ! Read the first line to get lmax
        read(input_unit,*, iostat=ios) lmax  ! Directly read lmax (no need for 'line' here)
        if (ios /= 0) then  ! If ios is not zero, there was an error reading
            write(error_unit,'(A)') "Error reading lmax from file: " // trim(filename)
            stop 
        end if
    
        ! Count the remaining lines (those containing integrals)
        num_lines = 0
        do
            read(input_unit,*, iostat=ios)  ! Ignore line content, just count
            if (ios /= 0) exit
            num_lines = num_lines + 1
        end do
        rewind(input_unit)  ! Go back to the beginning of the file
    
        ! Set expected shape (explicitly)
        expected_shape = [0, 0, 0, -lmax, -lmax, -lmax]  ! Lower bounds
        expected_shape(1:6) = expected_shape(1:6) + lmax + 1  ! Upper bounds
    
        ! Check existing allocation and reallocate if needed
        if (allocated(integrals_3rspH)) then
            if (any(shape(integrals_3rspH) /= expected_shape)) then
                deallocate(integrals_3rspH)  ! Deallocate if wrong shape
            end if
        end if
    
        ! Allocate if not allocated or deallocated above
        if (.not. allocated(integrals_3rspH)) then
            allocate(integrals_3rspH(0:lmax, 0:lmax, 0:lmax, -lmax:lmax, -lmax:lmax, -lmax:lmax))
        end if
        integrals_3rspH = 0.0d0  ! Initialize all elements to zero (double precision)
    
        ! Skip the first line with lmax
        read(input_unit,*)
        ! Read Integrals from File (Corrected)
        do line_indx = 1, num_lines
            ! Read and immediately assign values to variables
            read(input_unit,*, iostat=ios) l1, l2, l3, m1, m2, m3, integrals_3rspH(l1, l2, l3, m1, m2, m3)
        end do  ! End of line reading loop
    
        ! Close File 
        close(input_unit) 
    
    end subroutine read_integrals_3rspH

end module integrals_3rspH_mod
