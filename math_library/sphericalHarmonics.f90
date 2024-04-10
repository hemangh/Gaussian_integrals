module sphericalHarmonics
  use precision, only: idp
  use constants, only: pi
  use factorials_mod, only: compute_factorials, factorial, binomial_coeff
  use lmindx_map

  implicit none
  

contains



      !!************************************************************************************
      subroutine associated_legendre (plm,x,npt,lmax,m)
        !!for a fixed m finds plm for (l=m,lmax)
        integer                               :: i, l
        integer                               :: lmax, m, npt
        real(idp)                             :: x(:)
        real(idp),dimension(npt)              :: somx2
        real(idp)                             :: norm, fact
        real(idp),dimension(npt,m:lmax)       :: plm
        !----------------------------------------------------------------------
        !           initialize the arrays and do checks
        !----------------------------------------------------------------------

        plm = 0.d0

        if(m.lt.0.d0 .or. m .gt. lmax) then
          write(*,*) 'bad l and m arguments in plm'
          stop
        endif

        do i = 1, npt
          if(abs(x(i)).gt.1.d0)then
            write(*,*) 'bad grid (x) arguments in plm, x must be &
            &  cos_theta'
            stop
          end if
        end do

        !----------------------------------------------------------------------
        !          compute p^m_m = (-1)^m (1-x^2)^(m/2) (2m-1)!!
        !----------------------------------------------------------------------
        plm(:,m)=1.d0
        if (m > 0) then
          somx2(:) = sqrt((1.d0 - x(:)) * (1.d0 + x(:)))
          fact = 1.d0
          do i= 1, m
            plm(:,m) =  - plm(:,m) * fact * somx2(:)
            fact = fact + 2.d0
          end do
        end if
        !----------------------------------------------------------------------
        !          compute p^m_{m+1} = x (2 m + 1) p^m_m
        !----------------------------------------------------------------------
        if(m /= lmax)then
          plm(:,m + 1) = x(:) * (2 * m + 1) * plm(:,m)
          !----------------------------------------------------------------------
          !           start recursion with plm(m,m) and plm(m+1,m)
          !                      and recur upward
          !            p^m_l (l-m) = x (2l-1)p^m_{l-1} - (l+m-1) p^m_{l-2}
          !----------------------------------------------------------------------
          if(m+1 /= lmax)then
            do l = m + 2, lmax
              plm(:,l)=(x(:)*(2 * l - 1)*plm(:,l-1)-(l+m-1)*plm(:,l-2))/(l - m)
            end do
          endif
        endif
        !---------  normalized--------------------------------------------------

        !           do l = m, lmax
        !                norm= sqrt(0.5d0 * (2 *l + 1) * factorial(l - m )/factorial(l + m))
        !                plm(:,l) = norm * plm (:,l)
        !           end do

      end subroutine associated_legendre

      !***************************************************************************************
      ! Here negative is always associated with sin(m*phi) and positive with
      ! cos(m*phi).
      ! The spherical harmonics thus produced are known as real spherical harmonics
      ! or tesseral spherical harmonics.
      !*****************************************************************************************

      subroutine real_spherical_harmonics(ylm,theta,phi,npt,lmax)
        !!for a given lmax finds all real ylm (l=0,lmax) and (m=-lmax,lmax)
        integer                                        :: i
        integer                                        :: m,l
        integer                                        :: lmax, npt
        real(idp),dimension(npt)                       :: theta
        real(idp),dimension(npt)                       :: phi
        real(idp)                                      :: norm
        real(idp), allocatable                         :: x(:)
        real(idp)                                      :: fact
        real(idp)                                      :: condon_fac
        real(idp),dimension(:,:),allocatable           :: plm
        real(idp),dimension(:,:,:), allocatable        :: ylm
        if(allocated(ylm)) deallocate(ylm)
        allocate(ylm(npt,0:lmax,-lmax:lmax))
        !****it is only called once to make the factorials
        call compute_factorials
        allocate(x(npt))
        x(:) = cos(theta(:))
        !****loop over positive and zero m values
        do m=0,lmax
          !***find all non-normalized plm's for a fixed value of m
          allocate(plm(npt,m:lmax))
          call associated_legendre (plm,x,npt,lmax,m)

          !***loop over l values for a fixed m
          do l= m, lmax
            !**find condon-shortley phase factor of (-1)**m

            if(mod(m,2)==0)then
              condon_fac=1.d0
            else
              condon_fac=-1.d0
            endif
            !***normalization factor of ylm's
            if (m==0)then
              !fact=sqrt((2*l+1)*factorial(l)/(4.d0*pi)/factorial(l))
              fact=0.5d0*sqrt((2.d0*real(l)+1.d0)/pi)
              !write(iout,*)'fact=',fact, pi
              !write(iout,*)'l:',l
              !write(iout,*)'m:',m
              do i=1,npt
                ylm(i,l,0)=fact*plm(i,l)
                !write(iout,*)'ylm:',fact,plm(i,l),l,m,ylm(i,l,m)
              end do
            else
              !fact=condon_fac*sqrt((2*l+1)*factorial(l-m)/(2.d0*pi)/factorial(l+m))
              fact=condon_fac*sqrt((2.d0*l+1)*factorial(l-m)/(2.d0*pi)/factorial(l+m))
              ! write(iout,*)'fact=',fact
              ! write(iout,*)'l:',l
              ! write(iout,*)'m:',m
              ! write(iout,*)'(l-m)!:',factorial(l-m)!factorial(l-m)
              ! write(iout,*)'(l+m)!:',factorial(l+m)!factorial(l+m)
              ! write(iout,*)'fact:',fact
              !***evaluate ylm's on the quadrature points on the surface of a
              !sphere
              do i=1,npt
                !            write(iout,*)'plm:',plm(i,l)
                !            write(iout,*)'cos(phi):',cos(m*phi(i))
                !            write(iout,*)'sin(phi):',sin(m*phi(i))
                ylm(i,l, m) = fact * plm(i,l) * cos(m*phi(i)) !positive m part
                ylm(i,l,-m) = fact * plm(i,l) * sin(m*phi(i)) !negative m  part

                !write(iout,*)'ylm:',l,m,ylm(i,l,m)
                !write(iout,*)'ylm:',l,-m,ylm(i,l,-m)
              end do

            end if
          end do !l
          deallocate(plm)
        end do !m

      end subroutine real_spherical_harmonics


        !! Everything below probably needs to be moved to a Gaussian making module

    !=========== icart ===============================================================!
      ! cartesian descending ordering: idx = (L+1-i)*(j+k)/2 + k + 1
        function icart(i,j,k) result(indx)

          integer :: indx
          integer :: i,j,k
          integer :: L
      
          L = i + j + k
          indx = ((L+1-i)*(j+k)/2) + k + 1
      
          end function icart 
      
      
      !=========== Nlm ===============================================================!
          real(idp) function compute_Nlm(l,m)
            
          integer, intent(in) :: l, m
          integer :: m2
          real(idp) :: denom,fac
      
          m2 = abs(m);
      
          !compute N
          denom = 1.0;
          if (m == 0) then 
               denom = 2.0
          end if
          fac = 2.0 * factorial(l+m2) * factorial(l-m2) / denom
          compute_Nlm = sqrt(fac) / ( (2**m2)* factorial(l) )
      
          end function compute_Nlm 
      
      !=========== Ctuvlm ===============================================================!
          real(idp) function Ctuvlm(l,m,t,u,v,vm)
      
      
          integer, intent(in) :: l, m, t, u, v, vm
          integer :: m2, sgn
          real(idp) :: denom
      
          m2 = abs(m);
      
          ! compute C
          if (mod(t+v,2) == 0) then 
              sgn = 1;
          else 
              sgn = -1;
          end if
      
          denom = 4**t
          Ctuvlm = sgn * binomial_coeff(l,t) * binomial_coeff(l-t,m2+t) * binomial_coeff(t,u) * binomial_coeff(m2,vm) / denom
      
          end function Ctuvlm

    !=========== spherical trans for shell: ===============================================================!
     subroutine YLM_2_XYZ(l,  pwr, coeff)
      !
      ! YLM_2_XYZ is a set of codes to analytically compute the expansion
      ! of a set of real solid spherical harmonics in terms of cartesian
      ! gaussians.  These coefficients are used to transform cartesian integrals
      ! to spherical harmonic integrals.  The convention used for the order of the
      ! harmonics are descending in m.  More critically, the order of the primitives
      ! is needed for the transformation.  The powers are explicity computed as are the
      ! cartesian symbols.
      !
      !  Input: integer : l, the total power of the Cartesian Gaussian
      !                   as well as the degree of the Spherical Harmonics
      !                   part of the Spherical Gaussian
      !  Output: integer: pwr(1:3 -->[x,y,z],1:(l+1)*(l+2)/2) powers of the 
      !                   Cartesian Gaussian in the conventional 
      !                   ordering described in the icart function.
      !                   Notice that pwr is allocatable and is deallocated 
      !                   at every call if it is already allocated 
      !          doubles: coeff(1:2l+1 --> [-m,m],1:(l+1)*(l+2)/2 )
      !                   the coefficient of transformation of a Cartesian
      !                   to a Spherical Gaussian, Notice that coeff is
      !                   allocatable and is deallocated 
      !                   at every call if it is already allocated
        integer  , intent(in)                  :: l
        integer, allocatable, intent(out)      :: pwr(:,:)   ! on output it is allocated as (2l+1, (l+2)*(l+1)/2)
        real(idp), allocatable, intent(out)    :: coeff(:,:) ! on output it is allocated as (3 (x,y,z) , (l+2)*(l+1)/2)
        
        integer                                :: n_t
        integer                                :: n_v
        integer                                :: vv
        integer                                :: m
        integer                                :: mm
        integer                                :: num
        integer                                :: t
        integer                                :: u
        integer                                :: v
        integer                                :: i
        integer                                :: j
        integer                                :: k
        integer                                :: c1
        integer                                :: count
        real(idp)                              :: Nlm
        real(idp)                              :: Ctuvlm2
        real(idp)                              :: fac
                      
      !initializations--------------------------------------------------     
        call compute_factorials !****it is only called once to make the factorials
        if (allocated(coeff)) deallocate(coeff)
        if (allocated(pwr)) deallocate(pwr)
        fac = sqrt((2*l+1)/pi)*.5d0
        
        n_v = 1
        vv = 0
        num = (l+1) * (l+2 ) / 2
        allocate(coeff(2*l+1,num))
        allocate(pwr(3,num))
        Nlm = 1.0
        ! write(*,*) 'Entering YLM_2_XYZ for l = ',l
        count = 1
       DO m = -l, l
        ! DO loop = 1, 2*l+1
        !    m=sp%m_list(loop)
           mm = abs(m)
      !                      compute Nlm
           if (mm > 0) then 
               Nlm = compute_Nlm(l,m)
           else 
               Nlm = 1.0
           end if
      !     write(iout,*) l, Nlm
      !                      find n_t
           n_t = (l-mm)/2
      !                      find n_v
           if (m >= 0) then 
               n_v = mm/2 
           else 
               n_v = (mm-1)/2
           end if
      ! spherical idx: for descending order
      !sidx = ipure2(l2, m); 
      ! spherical idx: for Psi4 order
      ! loop over tuv
           DO t=0,n_t
              DO v=0,n_v
                 vv = 2*v
                 if (m < 0) then 
                     vv = vv + 1
                 end if
                 DO u=0,t
                    Ctuvlm2 = Ctuvlm(l, m, t, u, v, vv)
                    j = 2*(u+v)
                    if (m < 0) then
                        j = j + 1
                    end if
                    i = (2*t) + mm - j 
                    k = l - i - j
                    c1 = icart(i,j,k)
                    pwr(1,c1) = i
                    pwr(2,c1) = j
                    pwr(3,c1) = k
                    coeff(count,c1) = Ctuvlm2 * Nlm * fac
                    
                 END DO
              END DO
           END DO
           count = count + 1
      !
        END DO
        ! call print_matrix_int(pwr)
        ! call print_matrix(coeff)
        
     
        end subroutine YLM_2_XYZ

    end module sphericalHarmonics
