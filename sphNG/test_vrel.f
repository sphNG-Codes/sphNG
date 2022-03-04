      SUBROUTINE test_vrel
c************************************************************
c                                                           *
c  This routine can be called from step_P.F to estimate the *
c     relative grain velocities for a stand disc model.     *
c     It is not used by the SPH code itself.                *
c     It would need to be added into the Makefile if used.  *
c                                                           *
c     Matthew R. Bate, 15 December 2021                     *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/HY09accel'

      REAL omega

      idiff = 8
      ipart = 1
      radius = 0.005
      DO k = 0, 10

c         DO i = idiff + 1,HY09_ndust_bins
c            j = i - idiff
         DO i = k + 1,HY09_ndust_bins
            j = i - k

            temperature = 10.0
            rhoi = 1.0E-14/udens
            vthermal = SQRT(8.0*boltzmannk*temperature/(pi*gmw*mH))
c
c--Pre-compute coefficient for vrel from pressure gradient
c
         vrel_press_coeff = HY09_dust_density/(rhoi*udens*vthermal)*
     &        SQRT(gas_accel(1,ipart)**2 + gas_accel(2,ipart)**2 +
     &        gas_accel(3,ipart)**2)*udist/(utime**2)
c
c--Pre-compute coefficient for radial drift
c
         rgrav_accel = SQRT(grav_accel(1,ipart)**2 +
     &           grav_accel(2,ipart)**2 +
     &           grav_accel(3,ipart)**2)*udist/(utime)**2
         IF (rgrav_accel.EQ.0.) THEN
            WRITE (*,*) 'ERROR - rgrav_accel.EQ.0.'
            CALL quit(0)
         ENDIF
         dlnPdlnr = -1.0
         eta = - 0.5*(0.05)**2*dlnPdlnr
         vrel_rad_drift_coeff = -2.0*HY09_dust_density/
     &        (rhoi*udens*vthermal)*eta*rgrav_accel
c
c--Pre-compute quantities for turbulence relative velocities
c
         omega = SQRT(rgrav_accel/(radius*udist+solarr))
         vg2 = 1.0E-03*pi/8.0*vthermal**2

            vrel = HY09_vreldust(ipart,i,j,temperature,radius,rhoi,
     &           vthermal,
     &        vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &        vrel_settling_coeff,
     &        Stokes_k,Stokes_j,size_k,.TRUE.,.TRUE.)
            WRITE (99,100) HY09binsizes(i),vrel,Stokes_k,Stokes_j,size_k
 100        FORMAT(5(1PE12.5,1X))
         END DO
         radius = radius*2.0
      END DO

      STOP

      RETURN
      END
