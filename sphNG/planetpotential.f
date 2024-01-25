      SUBROUTINE planetpotential (i, px, py, pz, imigrate, rorbitmax,
     &     pmrate, time)

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/xforce'

      INTEGER i, imigrate
      REAL px, py, pz
      REAL rorbitmax, pmrate
      REAL time
      REAL rorbit, ang, angend

c--Allow for migration of potential at prescribed rate.
      IF (imigrate.EQ.1 .AND. time.LT.rorbitmax) THEN
         rorbit = planetsemiaxis(i) + (pmrate*time)
         ang = (2.*planetsemiaxis(i)*SQRT(xmass/planetsemiaxis(i)**3) - 
     &        2.*(pmrate*time+planetsemiaxis(i))*SQRT(xmass/(pmrate*
     &        time + planetsemiaxis(i))**3))/pmrate
      ELSEIF (imigrate.EQ.1 .AND. time.GE.rorbitmax) THEN
         rorbit = planetsemiaxis(i) + (pmrate*rorbitmax)
c--This only needs to be calculated once really, consider moving.
         angend = (2.*planetsemiaxis(i)*SQRT(xmass/planetsemiaxis(i)**3) 
     &     - 2.*(pmrate*rorbitmax+planetsemiaxis(i))*SQRT(xmass/(pmrate*
     &        rorbitmax + planetsemiaxis(i))**3))/pmrate
         
         ang = angend + (time-rorbitmax)/sqrt(rorbit**3/xmass)
      ELSE
         rorbit = planetsemiaxis(i)
         ang = time*sqrt(xmass/rorbit**3)
      ENDIF

      px = rorbit*COS(ang)
      py = rorbit*SIN(ang)
      pz = 0.0

      RETURN

      END SUBROUTINE planetpotential
c
c--------------------------------------------------------------------
c
c--Initial radius expansion factor for planets with surfaces.
c
c     For planets modelled using potentials it makes the effective
c     planet radius large (0.01 in code units) initially, which
c     decays exponentially based on the orbital period to allow 
c     gas to settle around the planet potential.
c
c     For planets modelled using sink particles with surfaces, makes 
c     the planets 10x larger initially (factor set in COMMONS/xforce)
c     and the factor pradfac_sinks that is used after that reduces 
c     exponentially at a rate set in derivi.
c
      FUNCTION pradfac(i,time)

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/physcon'

      INTEGER i
      REAL time
      REAL pradfac
      REAL period_over2pi

      IF (i.EQ.0) THEN  ! used for sink particles with surfaces
         IF (time.EQ.0.) THEN
            pradfac = pradfac_sinks_init
         ELSE
            pradfac = pradfac_sinks ! value evolved in derivi
         ENDIF
      ELSEIF (i.GT.numplanet) THEN
         WRITE (*,*) 'ERROR - pradfac: i.GT.numplanet'
         CALL quit(0)
      ELSE ! used for planets modelled using potentials
         period_over2pi = SQRT(planetsemiaxis(i)**3/xmass)
         IF (time.GT.40.0*period_over2pi) THEN
            pradfac = 1.
         ELSE
            pradfac =(planetradius(i) + 0.01*EXP(-time/period_over2pi))/
     &           planetradius(i)
         ENDIF
      ENDIF

      RETURN

      END FUNCTION pradfac
