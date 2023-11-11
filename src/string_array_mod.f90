module string_array_mod
    implicit none

    type string_array_type 
      character (len=:), allocatable :: string
    end type 

end module    