      SUBROUTINE bonnorebert(ximax)
c
c--This subrountine generates rho and M_enc as a function of radius for a
c     Bonnor-Ebert sphere.  Note:  A King model -> a Bonnor-Ebert sphere
c     in the limit that Phi(0)/sigma^2 -> Infinity (Binney and Tremaine, 
c     p. 233).
c
      IMPLICIT NONE

      INCLUDE 'idim'

      REAL ximax, bemasstotal, berad, bemaxrad, befac,
     &     xi, phi, func, dxi, dfunc, central_density,
     &     containedmass, dmass, conmassnext, dphi, rho
      INTEGER loop,i,j

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/bonnortbl'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/rbnd'

      REAL*4 berho(nbonnor)

 40   bemasstotal = 1.0
      berad = 1.0
      WRITE (*,*) 'Generating Bonnor-Ebert density profile...'
c
c-- Generate the Bonnor-Ebert profile up to the required ximax
c
      loop = 100000000
      xi = 0.0
      phi = 0.0
      func = 0.0
      containedmass = 0.0
      dmass = 1.0E-3
      conmassnext = dmass
      dxi = 1.0D+03/DBLE(loop)
      dfunc = (-EXP(phi))*dxi
      bonnor_radmass(1,1) = 0.
      bonnor_radmass(2,1) = 0.
      berho(1) = 1.
      j = 2
      DO i = 1, loop
         xi = i*dxi
         func = func + dfunc
         dphi =  func*dxi
         phi = phi + dphi
         dfunc = (-EXP(phi) - 2.0*func/xi)*dxi
         rho = EXP(phi)
         containedmass = containedmass + 4.0*PI*xi*xi*rho*dxi
         IF (containedmass.GE.conmassnext) THEN
            bonnor_radmass(2,j) = containedmass
            bonnor_radmass(1,j) = xi
            berho(j) = rho
            IF(xi.GT.ximax) GOTO 200
            conmassnext = conmassnext + dmass
            j = j + 1
            IF (j.GT.nbonnor) THEN
               WRITE (*,*) 'ERROR - Bonnor Ebert table too small'
               CALL quit(0)
            ENDIF
         ENDIF
      END DO
      WRITE(*,*) 'ERROR: xi is too large for loop ',xi,loop,ximax
      CALL quit(0)

 200  CONTINUE
      ibelast = j - 1
      WRITE(*,*) 'Ratio rho_central/rho_out (critical value 14.1)',
     &     1.0/rho
      WRITE(*,*) 'Sphere extends to xi = ', xi, ' table ',ibelast
      WRITE(*,*) 'Value of func at xi_max = ',func
      central_density = -ximax/(func)/(4*pi)
      WRITE(*,98001) central_density/rmax**3
98001 FORMAT(' Value of central density (code units) = ',
     &     1PE12.5,'*(mass in code units)')
      WRITE(*,98002) -1./(func*ximax)*gmw/Rg*(gg*umass/udist)/rmax
98002 FORMAT(' Equilibrium temperature = '1PE10.3,
     &     '*(mass in code units)')

      bemaxrad = bonnor_radmass(1,ibelast)
      befac = berad/bemaxrad 
c
c-- Scale the masses and radii appropriately.  Also write output file
c     that gives the density profile and cumulative mass profile.
c
      OPEN(13,FILE='BonnorEbert.txt')
      DO j = 1, ibelast
         bonnor_radmass(1,j) = bonnor_radmass(1,j) * befac
         bonnor_radmass(2,j) = bonnor_radmass(2,j) / 
     &        bonnor_radmass(2,ibelast)*bemasstotal

         WRITE (13,99001) bonnor_radmass(1,j),bonnor_radmass(2,j),
     &        berho(j)
99001    FORMAT(3(1PE12.5,1X))
      END DO
      CLOSE(13)

      RETURN

      END
