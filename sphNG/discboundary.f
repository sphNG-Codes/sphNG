      subroutine discboundary(ibound)

      INCLUDE 'idim'

      INCLUDE 'COMMONS/boundheight'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/cgas'

      REAL kappa
      REAL*4 rhocold
      REAL*8 inverf

c      OPEN (unit=36,file='zbound.txt',access='append')

c FM: produce array of the boundary height for radii from rmin to rcyl
      nsteps = 1000
      deltar = (rcyl-rmind)/(1.*nsteps)

      G = 1.0d0
      Mstar_init = 1.0d0
      R_o = 1.0d0
      R_out = 25.0d0

      radius = rmind
      do i=1, nsteps+1

         cv = 1.5*Rg/(gmw*uergg)         
         rhocold = 1.0E-15/udens  !convert from cgs to code units

         IF (ibound.EQ.101) boundtempl = (gmw/(gamma*Rg/uergg))
     &        *hoverr**2*G*Mstar_init/R_o*(R_o/radius)**0.5
         radius = 1.0

         IF (ibound.EQ.102) boundtempl = gmw*hoverr**2/((Rg/uergg)*
     &        gamma*radius)

         ucold = boundtempl*cv
         icounter = 0
         do
            cv = GETCV(rhocold,ucold)
            T_iterative = ucold/cv
            icounter = icounter+1

            IF (abs(boundtempl-T_iterative) .le. 1.0) THEN
               GOTO 102
            ELSE
               fraction = T_iterative/boundtempl
               ucold = ucold/fraction
            ENDIF

            IF (icounter .gt. 1000000) THEN
               GOTO 101
            ENDIF
         enddo

 101     write(*,*)'icount exceeded: ',boundtempl,T_iterative,icounter 
         STOP
 102     continue
         
         kappa = getkappa(ucold, cv, rhocold)

         IF (ibound.EQ.101) sigma_init = hoverr/(2*pi)*Mstar_init/
     &        (R_o*R_out**3)**0.25*(R_o/radius)
         IF (ibound.EQ.102) sigma_init = 75.0*udist**2/umass*
     &        signorm/radius**0.5

         x = 1.0-2.0/(kappa*sigma_init)
         zoverh(i) = sqrt(2.0)*abs(inverf(x))

c         write(36,12345) radius, zoverh(i), sigma_init,boundtempl,
c     &        x,inverf(x),kappa,zoverh(i)*0.05*radius
12345    FORMAT(8(2X,1PE12.5))

         radius = rmind + (i-1)*deltar

      END DO
c      close(36)

      RETURN
      END

c-------------------

      FUNCTION inverf(x)

c***************************************
c                                      c
c This function calculates the inverse c
c error function of a variable x       c 
c                                      c
c**************************************c


      IMPLICIT NONE
      real*8 inverf, inverfnew, x, pi, c_func, difference
      integer k, kmax
      parameter (kmax=1000000)

c DEFINE VARIABLES ETC HERE

      inverf = 0.0
      pi = acos(-1.0)

c      write(*,*)'x = ',x

c Sum the Maclaurins series for erf
      do k=0, kmax 

         inverfnew = inverf+c_func(k)/(2*k+1)*(sqrt(pi)/2*x)**(2*k+1)
c         write(*,*) 'inverf, inverfnew = ', inverf, inverfnew

         if (inverfnew-inverf .LT. 1e-7) then
            difference = inverfnew-inverf
            inverf = inverfnew
c            write(*,*) 'k, inverf = ', k, inverf
c            write(*,*) 'difference = ', difference
            go to 100
         else
            inverf = inverfnew
         end if

      end do

c      write(*,*)'kmax reached: increase kmax'
      stop

 100    continue


      RETURN
      END

c-------------------

      FUNCTION c_func(k)

      integer k, kmax, m
      parameter (kmax=1000000)
      real*8 c_func
      real, save :: cdummy(0:kmax)

      cdummy(k) = 0.0

      cdummy(0) = 1.0

      if (k .gt. 0) then
         do m=0, k-1
            cdummy(k)= cdummy(k) + cdummy(m)*cdummy(k-1-m)/(real(m+1))
     &           /(real(2*m+1))
         end do
      end if

      c_func = cdummy(k)

      RETURN
      END
