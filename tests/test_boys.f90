program boys_test
    use boys_function, only: boys
    implicit none
    ! double precision :: Boys_func
    integer :: n, io_error, nv
    double precision :: t, res, tv, mv, deps
    
    write(*,*) "-------------------------------------------------"// &
    "--------------------------------------------------------------------"
    write(*,'(A5,5A22)') "N", "X", "calculated value", "reference value",  &
     "difference", "rel. difference"
    write(*,*) "-------------------------------------------"// &
    "--------------------------------------------------------------------------"
    
    open(unit=20,file='../../input/benchmark.values',status='old',action='read',iostat=io_error)
    if ( io_error .eq. 0 ) then
    
            do
                    read(20,*,iostat=io_error) nv, tv, mv
                    if ( io_error .ne. 0 ) then
                            exit
                    end if
                    
                    res = boys(nv, tv)
                    if ( mv .gt. 0.0 ) then
                            write(*,'(I5,5E22.14)') nv, tv, res, mv, res-mv, (res-mv)/mv
                    else
                            write(*,'(I5,4E22.14, A22)') nv, tv, res, mv, res-mv, "--"
                    end if
            end do
    else
            write(*,*) "Error! Can't open 'benchmark.values'!"
    end if
    
    close(unit=20)
    
    write(*,*) "-----------------------------------------------------------"// &
    "----------------------------------------------------------"
            
end program
