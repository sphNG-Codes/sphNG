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
