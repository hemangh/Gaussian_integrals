    program test_lebedev
        use lebedev_quadrature
        use sphericalHarmonics, only: real_spherical_harmonics
        implicit none
        integer, parameter :: idp = selected_real_kind(15, 307)  ! Double precision kind
        real(idp), parameter :: pi = atan(1.0_idp) * 4.0_idp  

        type(Lebedev_type)     :: leb_ang
        integer           :: rule, lmax, l, m
        real(idp),dimension(:,:,:), allocatable        :: ylm   
        real(idp) :: integral
        rule = 5
        call leb_ang%get_rule_points(rule)
      lmax = leb_ang%precision/2
      print*, ' lmax : ', lmax
      ! Mapping Ylms up to lmax on the leedev grid
      call real_spherical_harmonics(ylm,leb_ang%theta,leb_ang%phi,leb_ang%num_points,lmax)
      do l = 0, lmax
          do m = -l, l
          ! check to see if the integral of Ylm * Ylm is one
          integral = sum(4.d0*pi*leb_ang%w(:)* ylm(:,l,m) ** 2)
          if(abs(integral-1._idp)>1.d-12)then
            print'(2(A,I0),F15.8)','l=',l, 'm=',m, integral
            call exit(1)
          else
            print *, 'Test passed: l = ', l, ', m = ', m
          end if
        end do
      end do

    end program