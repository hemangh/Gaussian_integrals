program test_integral_3rspH
    use precision, only: idp
    use integrals_3rspH_mod
    implicit none

    real(idp), allocatable :: integrals(:,:,:,:,:,:)
    character(len=100) :: filename = "../../input/integrals_3rspH.txt"  
    integer :: lmax

    call read_integrals_3rspH(integrals, filename, lmax)

    ! Now use the 'integrals' array for further calculations

end program test_integral_3rspH
