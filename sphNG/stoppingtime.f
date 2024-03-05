      FUNCTION stoppingtime(ipart,rhoi,dustvari)
c************************************************************
c                                                           *
c  One-fluid dust:  This function returns the stopping      *
c   time of particle ipart as given in Ballabio et al. 2018 *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'

      REAL stoppingtime 
      REAL graindens, s, f, rhoi
      REAL dustvari
      REAL s2, epsilon
      REAL thermal_velocity
      REAL vthermal

c      REAL, PARAMETER :: eps_thresh = 0.01005

      INTEGER ipart
      INTEGER K

c      LOGICAL switch
c
c - Set stopping time switch
c
c      switch = .FALSE.
c
c--Grain density of silicate
c
      graindens = 3./udens  
c
c--Size of particle
c - 0.1 micron
c      s = 1E-05/udist
c - 1 micron
c      s = 1E-04/udist
c - 3 micron
c      s = 3E-04/udist
c - 10 micron
c      s = 1E-03/udist
c - 30 micron
c      s = 3E-03/udist
c - 100 micron
c      s = 1E-02/udist
c - 300 micron
c      s = 3E-02/udist
c - 1 mm
c      s = 1E-01/udist
c - 1 cm
      s = 1.0/udist
c - thermal velocity 
c
      vthermal = thermal_velocity(ipart)
c--Supersonic correction factor not used, see PHANTOM 2018 paper
c      f = 1.
c
c - Switch for stopping time constraints
c
c      s2 = dustvari * dustvari
c      epsilon = s2 / (s2 + 1.)
c
c - Turning stopping time limiter on if dustfrac falls below a certain
c   threshold - eps_thresh
c
c      IF (epsilon.GT.eps_thresh) THEN
c
c - This is the Ballabio et al. 2018 formualtion with switch
c
c         stoppingtime = DMIN1(graindens*s*SQRT(pi*gamma/8.)/rhoi
c     &                     /vsound(ipart),1000*xyzmh(5,ipart)/
c     &                     vsound(ipart))
c
c         stoppingtime = DMIN1(graindens*s*SQRT(pi*gamma/8.)/rhoi
c     &                     /vsound(ipart),xyzmh(5,ipart)/vsound(ipart))
c         PRINT *,"Stopping time = ",stoppingtime
c
c - This is Ballabio et al. 2018 formulation without switch
c
c      ELSE
c         stoppingtime = graindens*s*SQRT(pi*gamma/8.)/rhoi/vsound(ipart)
c         print*, 'stoppingtime = ',stoppingtime
c
c - Stopping time limiter for Loren-Aguilar & Bate stopping time
c
c      stoppingtime = DMIN1(graindens*s/vthermal/rhoi,xyzmh(5,ipart)/
c     &                          vsound(ipart))
c      ENDIF
c
c - Stopping time Bate & Loren-Aguilar 2017
c
c
c
      stoppingtime = graindens * s / vthermal / rhoi
c
c - calculate dust fraction      
c      s2 = dustvari*dustvari
c      epsilon = s2 / (s2 + 1.)
c - set drag coefficient
c      K = 1000
c - calculate stopping time
c      stoppingtime = epsilon * (1.0 - epsilon) * rhoi / K

c      stoppingtime = 0.005

      RETURN
      END
