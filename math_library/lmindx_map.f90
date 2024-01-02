module lmindx_map
  implicit none
   private
   public make_lm_indx, indx_to_lm, indx_range, lookup_l
   public cnvt_indx_to_lm
   public kappamax_2_lmax, lmax_2_kappamax

  type one_to_two_indx
    integer  :: l
    integer  :: m
  end type



  type(one_to_two_indx), allocatable :: indxlm_series(:)

contains

  !>given l and m the angular and magnetic quantum numbers the
  !*function gives back a unique indx or index number for the pair of lm.
  function make_lm_indx(l,m) result(indx)

    integer              :: l,m
    integer              :: indx

      indx = l **2 + 1 + l - m  

  end function make_lm_indx

  function indx_range(l) result (karray)
    
    integer              :: l
    integer              :: karray(2)

    karray(1) = l * l + 1 
    karray(2) = l * l + 1 + 2 * l
  
  end function  

  function kappamax_2_lmax(kappamax) result(lmax)
    
    integer :: kappamax
    integer :: lmax 

    lmax =  int(sqrt(real(kappamax))) - 1 

  end function

  function lmax_2_kappamax(lmax) result(kappamax)
    
    integer :: kappamax
    integer :: lmax 

    kappamax = (lmax + 1)*(lmax + 1) 

  end function

  !> given a maximum angular momentum number,
  !* a map of individual {l,m} series is 
  !* is constructed as type one_to_two_indx
  !* described above. This is a common parameter
  !* to any program that uses this module
  !* the purpose is then that later, given 
  !* any oe index number kappa, the (l,m) pair
  !* is extracted
  subroutine maplmindx(lmax)
    integer, intent(in) :: lmax
    integer             :: max_num
    integer             :: k, l, m
    if(allocated(indxlm_series))then
      deallocate(indxlm_series)
    end if
    max_num = (lmax + 1)**2
    allocate(indxlm_series(1:max_num))
    k = 1
    indxlm_series(k)%l=0
    indxlm_series(k)%m=0

    do l=1, lmax
      do m = -l, l
        k = k + 1
        indxlm_series(k)%l = l
        indxlm_series(k)%m = m
      end do
    end do

  end subroutine maplmindx


  !>given indx or index number for the pair of lm the
  !*function gives back a unique pair l and m, 
  !*the angular and magnetic quantum numbers,
  !* in a two dimensional integer array lm:
  !* lm(1) = l ; lm(2) = m
  function indx_to_lm(indx,lmax) result(lm)
    integer, dimension(2)::lm
    integer              :: indx
    integer, optional    :: lmax

    if(.not. allocated(indxlm_series)) then
      if(present(lmax))then
        call maplmindx(lmax)
      else
        !call lnkerr('lmax optional argument must be provided the first time to make a map in lmindx_map.f90'); stop
        stop 'lmax optional argument must be provided the first time to make a map in lmindx_map.f90'
      end if
    else
      if(present(lmax))then
        call maplmindx(lmax)  
      endif    
    end if

    lm(1)=indxlm_series(indx)%l
    lm(2)=indxlm_series(indx)%m
  end function indx_to_lm


  !>given indx or index number for the pair of lm the
  !*function gives back a unique pair l and m, 
  !*the angular and magnetic quantum numbers
  !*if called with lmax argument, it maps lmax to the indx
  !*in the first place. Therefore, the subroutine has to be called 
  !*with lmax input the first time or when lmax changes, but subsequently
  !*it is not needed. 
  subroutine cnvt_indx_to_lm(indx,l,m,lmax)
    
    integer, intent(in)  :: indx
    integer, intent(out) :: l,m
    integer, optional    :: lmax
    if(.not. allocated(indxlm_series)) then
      if(present(lmax))then
        call maplmindx(lmax)
      else
        !call lnkerr('lmax optional argument must be provided the first time to make a map in lmindx_map.f90'); stop
        stop 'lmax optional argument must be provided the first time to make a map in lmindx_map.f90'
      end if
    else
      if(present(lmax))then
        call maplmindx(lmax)  
      endif   
    end if
    l=indxlm_series(indx)%l
    m=indxlm_series(indx)%m
  end subroutine cnvt_indx_to_lm


  !>given lmax and lambda number, all different combinations of
  !* |l_1 - l_2|<= lambda <= l_1+l_2 are listed in a two dimensional array
  subroutine list_l1_l2_for_lambda(lambda,lmax,l12combo)

    integer, intent(in) :: lambda
    integer, intent(in) :: lmax
    integer             :: l1, l2
    integer             :: lcombomax
    integer, allocatable, intent(out):: l12combo(:,:)
      
      lcombomax = 0
      do l1=0,lmax
        do l2=0,lmax

          if(abs(l2-l1)<=lambda .and. l1+l2>=lambda)then
            lcombomax = lcombomax + 1
          end if
        
        end do
      end do   
      
      allocate(l12combo(2,lcombomax)) 

           
      lcombomax = 0
      do l1=0,lmax
        do l2=0,lmax

          if(abs(l2-l1)<=lambda .and. l1+l2>=lambda)then
            lcombomax = lcombomax + 1
            l12combo(1,lcombomax) = l1
            l12combo(2,lcombomax) = l2

          end if
        
        end do
      end do   

  end subroutine list_l1_l2_for_lambda

  function lookup_l(symb) result(l)
    character(len=1) :: symb
    integer :: l

    select case(symb)

      case('s') 

        l = 0
      
      case('p')

        l = 1

      case('d')
        
        l = 2

      case('f')
        
        l = 3
      
      case('g')

        l = 4

      case ('h')
        
        l = 5

      case default
      
         stop 'no such orbital symbol'

      end select   

    end function lookup_l

end module
