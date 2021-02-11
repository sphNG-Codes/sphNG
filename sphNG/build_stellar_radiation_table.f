      SUBROUTINE build_stellar_radiation_table
c
c--Calculate integral of Q_v * B_v (absorption times Planck function)
c     and the integral of B_v (i.e. without absorption), but both with
c     extinction for the (gray) stellar radiation field from sink
c     particles.  Builds a lookup table with extinction values
c     from A_V = 10^-3 to A_V = ~1200 in log-space, and for stellar
c     temperatures T=1000 to ~55,000 K.
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/stellarradiation'

      xSTARfac = EXP(xSTARlogfac)
      dlogv = LOG(1.1)

      DO k = 1, nstellar_T
         T = 1000.0 * xSTARfac**(k-1)

         A_V_init = 1.0E-3
         DO j = 1, nstellar_Av
            A_V = A_V_init*xSTARfac**(j-1)

            x_integral = 0.
            x_integral2 = 0.
            DO i = 1, 1000
               v = (1.1**i)
               v = v / 1.0E+10
               xJ = planck(v,T)
               IF (xJ.GT.1e-40 .AND. v.GT.0.1) THEN
                  Qtot = Qv(v,Qabs)
                  Qtot550 = Qv(c/0.0000550,Qabs550)
                  quantity = EXP(-A_V/1.086*Qtot/
     &                 Qtot550 ) * xJ*v*dlogv
                  x_integral = x_integral + quantity * Qabs
                  x_integral2 = x_integral2 + quantity
               ENDIF
            END DO
c
c--NOTE: The integral over the Planck function is sigma_B*T^4/pi, while
c     the star's luminosity = 4*pi*R^2 * sigma_B*T^4
c     This can checked by commenting out the absorption terms in the 
c     integral.
c
            stellar_AvT_table(1,j,k) = x_integral * pi
            stellar_AvT_table(2,j,k) = x_integral2 * pi

c            WRITE (40,*) T,A_V,stellar_AvT_table(1,j,k)
         END DO
      END DO

c      STOP

      RETURN
      END

c-----------------------------------------------------------

      SUBROUTINE build_photoionisation_heating_rate_table
c
c--Heating rate due to hydrogen ionisation by UV photons.
c     Heating rate equation from, e.g. Kannan et al. (2019), eqn. 57
c     Cross-section is from Verner et al. (1996)
c     Depends on stellar spectrum, so tabulate based on stellar
c     blackbody temperature, from T=1000 to ~55,000 K.
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/stellarradiation'

      xSTARfac = EXP(xSTARlogfac)
      dlogv = LOG(1.1)

      DO k = 1, nstellar_T
         T = 1000.0 * xSTARfac**(k-1)

         x_integral = 0.
         x_integral_norm = 0.
         x_integral_2 = 0.
         x_integral_ionisingflux = 0.
         x_integral_N_ionphotons = 0.
c
c--Start from 13.6 eV photons (v=13.6 eV / h)
c
         v_ion = 13.6 * eleccharge / planckconst
         v = v_ion
         DO i = 1, 400
            xJ = planck(v,T)
            IF (xJ.GT.1e-40 .AND. v.GT.0.1) THEN
               sigma = hydrogen_cross_section(v)

               energyfac = planckconst*(v - v_ion)
               x_integral = x_integral + xJ*sigma*energyfac*v*dlogv

               x_integral_2 = x_integral_2 + 
     &              xJ*sigma*v_ion*planckconst*v*dlogv
               x_integral_norm = x_integral_norm + xJ*sigma*v*dlogv

               x_integral_ionisingflux = x_integral_ionisingflux +
     &              xJ*v*dlogv
               x_integral_N_ionphotons = x_integral_N_ionphotons +
     &              xJ*v*dlogv/(planckconst*v)

c               IF (k.EQ.88) WRITE (41,33001) v,sigma,!energyfac,
c     &              xJ*sigma*energyfac*v*dlogv,
c     &              xJ*sigma*v_ion*planckconst*v*dlogv
c     &              xJ*v*dlogv,planckconst*v*xJ*v*dlogv

            ENDIF
            v = v * 1.1
         END DO

         photoionisation_heating_T_table(k) = x_integral * pi

         IF (x_integral_norm.EQ.0.) x_integral_norm = pi

c         WRITE (90,33001) T,photoionisation_heating_T_table(k),
c     &        x_integral_2/x_integral_norm, 
c     &        photoionisation_heating_T_table(k)/x_integral_norm,
c     &        x_integral_N_ionphotons *pi
c     &        * 4.0*pi*(7.0 * 7.0E+10)**2
33001    FORMAT(5(1PE12.5,1X))

      END DO

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION hydrogen_cross_section(v)
c
c--Hydrogen photoionisation cross-section is from Verner et al. (1996)
c
c     Units are cm^2 (i.e. cm^2/H)
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT NONE

      REAL*8, PARAMETER :: P = 2.963D0
      REAL*8, PARAMETER :: ya = 3.288D1

      REAL*8 hydrogen_cross_section,v,E,x,y

      INCLUDE 'COMMONS/physcon'

      E = planckconst * v / eleccharge

      x = E/0.4298
      y = x
      hydrogen_cross_section = 5.475E-14 * (x - 1.0)**2 * 
     &     y**(0.5*P-5.5) * (1.0 + SQRT(y/ya))**(-P)

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION photoionisation_heating(Tstar)
c
c--Returns the photoionisation heating rate for a star with surface
c     temperature Tstar.
c     The heating rate is in erg/s per H.  But need to multiply
c     by the ratio of (Rstar/distance)^2 to get
c     the local heating rate in units of erg/s per H.
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/stellarradiation'

      REAL*8 photoionisation_heating,Tstar,xtemperature
      INTEGER itemperature

      IF (Tstar.LE.3000.) THEN
         photoionisation_heating = 0.
      ELSE
         xtemperature = LOG(Tstar/1000.0)/xSTARlogfac + 1.5
         itemperature = MIN(nstellar_T-2, INT(xtemperature)) + 1
         IF (xtemperature.GT.nstellar_T) xtemperature = nstellar_T
         
         photoionisation_heating = 
     &        photoionisation_heating_T_table(itemperature) + 
     &        (photoionisation_heating_T_table(itemperature + 1) - 
     &        photoionisation_heating_T_table(itemperature)) *
     &        (xtemperature - itemperature + 1)
      ENDIF
      
      RETURN
      END
