module boys_function

       implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2021 The Numericus Group, LLC.  All rights reserved.
! You may use, distribute, and modify this code ONLY under the terms
! of the TNG IP license.
!
! You should have received a copy of the TNG IP license with this
! file.  If not, please write to: The Numericus Group, LLC,
! 2525 Arapahoe Ave. E4-431, Boulder, CO 80302.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   NAME:dboysfun12
!
!   DESC: Computes values of the Boys function for n=0,1,...  12
!         for a real argument
!
!        Input: x  --- argument, real *8, x >= 0
!
!        Output: vals  --- values of the Boys function, n = 0,1,...  12
!
!
!   DATE: MAY 2021
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       subroutine dboysfun12(x,vals)
        implicit none
        character, parameter :: name*12="dboysfun12: "
 !
        real *8, parameter :: tol= 1.0d-03
        real *8, parameter :: sqrtpio2=0.886226925452758014d0
        real *8, parameter :: t(0:  11)=[ &
                         0.20000000000000000D+01,&
                         0.66666666666666663D+00,&
                         0.40000000000000002D+00,&
                         0.28571428571428570D+00,&
                         0.22222222222222221D+00,&
                         0.18181818181818182D+00,&
                         0.15384615384615385D+00,&
                         0.13333333333333333D+00,&
                         0.11764705882352941D+00,&
                         0.10526315789473684D+00,&
                         0.95238095238095233D-01,&
                         0.86956521739130432D-01]
        complex *16, parameter :: zz(1:  10)=[ &
                       (  0.64304020652330500D+01,  0.18243694739308491D+02),&
                       (  0.64304020652330500D+01, -0.18243694739308491D+02),&
                       ( -0.12572081889410178D+01,  0.14121366415342502D+02),&
                       ( -0.12572081889410178D+01, -0.14121366415342502D+02),&
                       ( -0.54103079551670268D+01,  0.10457909575828442D+02),&
                       ( -0.54103079551670268D+01, -0.10457909575828442D+02),&
                       ( -0.78720025594983341D+01,  0.69309284623985663D+01),&
                       ( -0.78720025594983341D+01, -0.69309284623985663D+01),&
                       ( -0.92069621609035313D+01,  0.34559308619699376D+01),&
                       ( -0.92069621609035313D+01, -0.34559308619699376D+01)]
        complex *16, parameter :: fact(1:  10)=[ &
                       (  0.13249210991966042D-02,  0.91787356295447745D-03),&
                       (  0.13249210991966042D-02, -0.91787356295447745D-03),&
                       (  0.55545905103006735D-01, -0.35151540664451613D+01),&
                       (  0.55545905103006735D-01,  0.35151540664451613D+01),&
                       ( -0.11456407675096416D+03,  0.19213789620924834D+03),&
                       ( -0.11456407675096416D+03, -0.19213789620924834D+03),&
                       (  0.20915556220686653D+04, -0.15825742912360638D+04),&
                       (  0.20915556220686653D+04,  0.15825742912360638D+04),&
                       ( -0.94779394228935325D+04,  0.30814443710192086D+04),&
                       ( -0.94779394228935325D+04, -0.30814443710192086D+04)]
        complex *16, parameter :: ww(1:  10)=[ &
                       ( -0.83418049867878959D-08, -0.70958810331788253D-08),&
                       ( -0.83418050437598581D-08,  0.70958810084577824D-08),&
                       (  0.82436739552884774D-07, -0.27704117936134414D-06),&
                       (  0.82436739547688584D-07,  0.27704117938414886D-06),&
                       (  0.19838416382728666D-05,  0.78321058613942770D-06),&
                       (  0.19838416382681279D-05, -0.78321058613180811D-06),&
                       ( -0.47372729839268780D-05,  0.58076919074212929D-05),&
                       ( -0.47372729839287016D-05, -0.58076919074154416D-05),&
                       ( -0.68186014282131608D-05, -0.13515261354290787D-04),&
                       ( -0.68186014282138385D-05,  0.13515261354295612D-04)]
        real *8, parameter :: rzz(1:   1)=[ &
                        -0.96321934290343840D+01]
        real *8, parameter :: rfact(1:   1)=[ &
                         0.15247844519077540D+05]
       real *8, parameter :: rww(1:   1)=[ &
                         0.18995875677635889D-04]
        real *8, intent(in)  :: x
        real *8, intent(out) :: vals(0:  12)
        real *8    :: y, yy,  rtmp
        real *8    ::   p, q, tmp
        integer *4 :: n, k
 !
           y = exp(-x)
 !
           if (abs(x).ge. 0.45425955121971775D+01) then
           yy      = sqrt(x)
           vals(0) = sqrtpio2*erf(yy)/yy
          yy = y/2.0d0
           do n = 1,  12
           vals(n) = ((n -0.5d0)*vals(n-1) - yy)/x 
           enddo
           return
           endif
 !
           rtmp = 0
           do k = 1,   10,2
           rtmp = rtmp + ww(k)*(1.0d0 - fact(k)*y)/(x + zz(k))
           enddo
 !
           tmp = 0
           do k = 1,    1
           if (abs(x + rzz(k)).ge.tol) then 
           tmp = tmp + rww(k)*(1.0d0 - rfact(k)*y)/(x + rzz(k))
           else
           q = x+rzz(k)
           p = 1.0d0 - q/2.0d0 + q**2/6.0d0 - q**3/24.0d0 + q**4/120.0d0
           tmp = tmp +rww(k)*p
           endif
           enddo
 !
           vals(  12) = 2*rtmp+tmp
           yy = y/2.0d0
           do n =   11,0,-1
           vals(n) = (x*vals(n+1)+yy)*t(n)
           enddo
 !
      return
      end subroutine dboysfun12
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! END: dboysfun12
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      function boys(n,x) result(val)
       integer :: n
       real(kind=8) :: x, val
       real(kind=8) :: vals(0:12)
       call dboysfun12(x, vals)
       if(n > 12 .or. n < 0 ) stop 'wrong order given to the boys function 0<=n<=12'
       val = vals(n)
      end function

end module