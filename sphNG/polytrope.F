      SUBROUTINE polytrope(poly_index,ximax)
c
c--This subrountine generates rho and M_enc as a function of radius for a
c     polytropic sphere.  It is very similar to bonnorebert.f which
c     does the same thing for a Bonnor-Ebert sphere.
c

#ifdef MPIALL
#define IMPLICITNONE
#include "mpi_sup.h"
#else
      IMPLICIT NONE
#endif
      INCLUDE 'idim'
#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
#endif

      REAL poly_index

      REAL ximax, bemasstotal, berad, bemaxrad, befac,
     &     xi, phi, func, dxi, dfunc, 
     &     containedmass, dmass, conmassnext, dphi, rho
      INTEGER loop,i,j

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/bonnortbl'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/logun'

      REAL*4 berho(nbonnor)
      REAL*4 rhoc(nbonnor),cs2(nbonnor)

      bemasstotal = 1.0
      berad = 1.0
      WRITE (*,*) 'Generating Polytrope density profile...'
c
c-- Generate the Polytrope profile up to the required ximax
c
      loop = 100000000
      xi = 0.0
      phi = 1.0
      func = 0.0
      containedmass = 0.0
      dmass = 1.0E-3 
      conmassnext = dmass
      dxi = 1.0D+03/DBLE(loop)
      dfunc = -(phi**poly_index)*dxi
      bonnor_radmass(1,1) = 0.
      bonnor_radmass(2,1) = 0.
      berho(1) = 1.
      j = 2
      DO i = 1, loop
         xi = i*dxi
         func = func + dfunc
         dphi =  func*dxi
         phi = phi + dphi
         dfunc = (-phi**poly_index - 2.0*func/xi)*dxi
         rho = phi**poly_index
         containedmass = containedmass + 4.0*PI*xi*xi*rho*dxi
cWRITE(*,*) 'func=',xi*xi*func,'xi=',xi
cWRITE(*,*) 'containedmass=',containedmass,'phi=',phi
         IF(xi.GT.ximax) GOTO 200
         IF (containedmass.GE.conmassnext) THEN
            bonnor_radmass(2,j) = containedmass
            bonnor_radmass(1,j) = xi
            berho(j) = rho
             rhoc(j) = ximax**3/(4.0*PI*xi*xi*func)*udens
              cs2(j) = ximax/(xi*xi*func)*(gg*umass/udist)
     &               * (gmw/Rg)
            conmassnext = conmassnext + dmass
            j = j + 1
            IF (j.GT.nbonnor) THEN
               WRITE (*,*) 'ERROR - polytrope table too small'
               CALL quit(0)
            ENDIF
         ENDIF
      END DO
      WRITE(*,*) 'ERROR: xi is too large for loop ',xi,loop,ximax,
     & "conmassnext=",conmassnext
      CALL quit(0)

 200  CONTINUE
      ibelast = j - 1
      WRITE(*,*) 'Ratio rho_central/rho_out (critical value 14.1)',
     &     1.0/rho
      WRITE(*,*) 'Sphere extends to xi = ', xi, ' table ',ibelast
      WRITE(*,*) 'Value of func at xi_max = ',func
      central_density = -ximax/(func)/(4*pi)
      WRITE(*,*) 'Value of central density (code units) = ',
     &     central_density
      poly_const =  4.0*pi/(poly_index+1)/ximax**2 *
     &     central_density**(1.0-1.0/poly_index)
      WRITE(*,*) 'Value of poly_const = ',poly_const

      bemaxrad = bonnor_radmass(1,ibelast)
      befac = berad/bemaxrad 
      WRITE(*,*) "rmax=",bemaxrad,"temp=",cs2(ibelast)
c
c--Scale the masses and radii appropriately.
c
      DO j = 1, ibelast
         bonnor_radmass(1,j) = bonnor_radmass(1,j) * befac
         bonnor_radmass(2,j) = bonnor_radmass(2,j) /
     &        bonnor_radmass(2,ibelast)*bemasstotal
      END DO
c
c--Write output file that gives the density profile and
c     cumulative mass profile.
c
#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
      OPEN(idisk3,FILE='Polytrope.txt',
     &   STATUS='replace',ACTION='write')
      DO j = 1, ibelast
         WRITE (idisk3,99001) bonnor_radmass(1,j),bonnor_radmass(2,j),
     &        berho(j),rhoc(j),cs2(j)
99001    FORMAT(5(1PE12.5,1X))
      END DO
      CLOSE(idisk3)
#ifdef MPI
      ENDIF
#endif

      RETURN
      END
