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
      R_o = 1.0d0
      R_out = 25.0d0
      sigma_init = 0.

      radius = rmind
      do i=1, nsteps+1

         cv = 1.5*Rg/(gmw*uergg)         
         rhocold = 1.0E-15*udist**3/umass  !convert cgs to code units

         IF (ibound.EQ.101) boundtempl = (gmw/(Rg/uergg))*
     &        centralmass*hoverr**2*G/R_o*(R_o/radius)**0.5

         IF ((ibound.EQ.102 .OR. ibound.EQ.103) .AND. use_tprof) THEN
            boundtempl = gmw*centralmass*hoverr**2*radius**
     &           (tprof+1)/((Rg/uergg)*radius)
         ELSEIF (ibound.EQ.102 .OR. ibound.EQ.103) THEN
            boundtempl = gmw*centralmass*hoverr**2/((Rg/uergg)*
     &           radius)
         ELSE
            print *, 'You should not be in discboundary.f'
            CALL quit(0) 
         ENDIF

         ucold = getu(rhocold, boundtempl)
         cv = getcv(rhocold,ucold)

         kappa = getkappa(ucold, cv, rhocold)

         IF (ibound.EQ.101) sigma_init = hoverr/(2*pi)*centralmass/
     &        (R_o*R_out**3)**0.25*(R_o/radius)
         IF (ibound.EQ.102.OR. ibound.EQ.103) 
     &        sigma_init = (75.0*udist**2/umass)*
     &        signorm/radius**0.5

         x = 1.0-2.0/(kappa*sigma_init)
         zoverh(i) = 0.0
         IF (x .gt. 0.0 .AND. x.lt. 1.0) THEN
            zoverh(i) = sqrt(2.0)*abs(inverf(x))
         ENDIF
         zoverh(i) = MAX(1.5, zoverh(i))
c--This should be in H as above, NOT the value of H as below.
c   zoverh(i) = MAX(hoverr*radius**(0.5*(tprof+3))*1.5, zoverh(i))

c      write(36,12345) radius, zoverh(i), sigma_init,boundtempl,
c     &  x,kappa,zoverh(i)*hoverr*radius**(0.5*(tprof+3)),ucold/cv
c12345    FORMAT(8(2X,1PE12.5))

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
      integer k, kmax, sign
      parameter (kmax=1000001)

      sign = 1
      IF (x.ge.1.0 .OR. x.le.-1.0) THEN
         print *, 'ERROR: Invalid x value for InvErf'
         CALL quit(0) 
      ENDIF
      IF (x .lt. 0.0) then
         sign = -1
      ENDIF

      inverf = 0.0
      pi = acos(-1.0)
      
c Sum the Maclaurins series for erf
      do k=0, kmax 
         inverfnew = inverf+c_func(k)/(2*k+1)*(abs(x)*sqrt(pi)/2.0)
     &        **(2*k+1)
c         write(*,*) 'inverf, inverfnew = ', inverf, inverfnew

         if (inverfnew-inverf .LT. 1e-7) then
            difference = inverfnew-inverf
            inverf = inverfnew
            go to 100
         else
            inverf = inverfnew
         end if

      end do

      stop

 100  continue

      inverf = inverf*real(sign)

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
