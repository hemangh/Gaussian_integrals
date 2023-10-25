module constants
  use precision, only : idp
  implicit none

  integer, private, parameter  :: maxfac = 30 
  real(idp), parameter  :: pi =  atan(1.0_idp) * 4.0_idp 
  real(idp), dimension(0:maxfac) :: dfct
  real(idp), dimension(-1:maxfac) :: ddfct

contains

  subroutine factorials
    integer:: i
    logical:: called=.false.
    !implicit integer (a-z)
    !real *8 dfct, ddfct
    !dimension dfct(0:maxfac), ddfct(0:maxfac)

    if(called)then
      !write(iout,*) "*****factorial is already called******"
      return
    else
      write(*,*) "first time factorials is called"
      !----------------------------------------------------------------------c
      !               calculate factorials                                   c
      !----------------------------------------------------------------------c

      dfct(0)=1.d+00
      dfct(1)=1.d+00
      if (maxfac.gt.1) then
        do 10 i=2,maxfac
          dfct(i)=i*dfct(i-1)
          10     continue
        endif
        !----------------------------------------------------------------------c
        !           calculate (2*m-1) double factorial                         c
        !----------------------------------------------------------------------c
        ddfct(-1) = 1.d+00
        ddfct(0)=1.d+00
        ddfct(1)=1.d+00
        ddfct(2)=3.d+00
        if (maxfac.gt.2) then
          do 20 i=3,maxfac
            ddfct(i)=(i+i-1)*ddfct(i-1)
            20    continue
          endif
          called = .true.
          return
        end if
      end subroutine factorials

    end module constants
