      SUBROUTINE wrinsph
c************************************************************
c                                                           *
c  This subroutine writes a new file inname from parameters *
c     taken from setpart                                    *
c                                                           *
c************************************************************
#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'idim'

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/crpart'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/varmhd'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine wrinsph')
c
c--Open input file
c
#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
      OPEN (iterm, FILE=inname)
c
c--Determine options for evolution run
c
c--Write name of run
c
      WRITE (iterm, 99002) namenextrun
99002 FORMAT (A20)
c
c--Write name of file containing physical input
c
      WRITE (iterm, 99003) file1
99003 FORMAT (A7)
      WRITE (iterm, 99003) varsta
      IF (imhd.EQ.idim) WRITE (iterm, 10001) varmhd
10001 FORMAT(A4)
c
c--Write options
c
      WRITE (iterm, 99006) encal
99006 FORMAT (A1)
      IF (encal.EQ.'r') THEN
         WRITE (iterm, 99008) tolerance, boundtemp
99008    FORMAT(2X,1F6.4,1X,1PE12.5,'  Implicit radiative transfer',
     &        ' tolerance boundtemp')
      ENDIF

      WRITE (iterm, 89000) initialptm
89000 FORMAT(2X,I2,'  Point Masses Initially')
      WRITE (iterm, 89005) iaccevol
89005 FORMAT(A1,'  Variable:Roche,Sep/Fixed Accretion Radii (v/s/f)')
      IF (iaccevol.EQ.'v' .OR. iaccevol.EQ.'s')
     &     WRITE (iterm, 89010) accfac
89010 FORMAT(2X,1PE12.5,'  Fraction of Roche lobe size')
      WRITE (iterm, 89015) iptmass
89015 FORMAT(2X,I2,'  Point Mass Creation')
      WRITE (iterm, 88000) igrp
88000 FORMAT(2X,I2,'  Pressure')
      WRITE (iterm, 88001) igphi
88001 FORMAT(2X,I2,'  Gravity')
      WRITE (iterm, 88002) ifsvi,alpha,beta
88002 FORMAT(2X,I2,1X,1F6.3,1X,1F6.3,'  Artificial Viscosity')
      WRITE (iterm, 88003) ifcor
88003 FORMAT(2X,I2,'  Coriolis Forces')
      WRITE (iterm, 88004) ichoc
88004 FORMAT(2X,I2,'  Heating from Shocks')
      WRITE (iterm, 88005) iener
88005 FORMAT(2X,I2,'  PdV')
      WRITE (iterm, 88006) damp
88006 FORMAT(2X,1PE12.5,'  Damping')
      WRITE (iterm, 88007) ibound
88007 FORMAT(1X,I3,'  Boundry')
      WRITE (iterm, 88008) iexf
88008 FORMAT(2X,I2,'  External Forces')
      WRITE (iterm, 88009) iexpan
88009 FORMAT(2X,I2,'  Expansion')
      WRITE (iterm, 88010) nstep
88010 FORMAT(2X,I2,'  Binary Dump Every N Max Time Steps')
      WRITE (iterm, 88017) nfullstep
88017 FORMAT(2X,I2,'  Full Binary Dump Every N Binary Dumps')
      WRITE (iterm, 89020) iptoutnum
89020 FORMAT(2X,I4,'  N Ptmass Dumps Per Max Timestep') 
      WRITE (iterm, 88011) tol, tolptm, tolh
88011  FORMAT(2X,3(1PE12.5,1X),'Tolerance (Gas, Ptmass, Smooth Len)')
      WRITE (iterm, 88012) ipos
88012  FORMAT(2X,I3,'  File Position')
      WRITE (iterm, 88013) tmax
88013  FORMAT(2X,1PE12.5,'  Max Time (min) (Status Written)')
      WRITE (iterm, 88014) tstop  
88014  FORMAT(2X,1PE12.5,'  Max Dynamic Time')
      WRITE (iterm, 88015) dtmax
88015  FORMAT(2X,1PE14.7,'  Max Timestep')
      WRITE (iterm, 88016) dtini
88016  FORMAT(2X,1PE12.5,'  Initial Timestep')

      IF (ifcor.NE.0) THEN
         IF (job(1:9).EQ.'evolution') THEN
            omeg0r = omeg0/utime
         ELSE
            omeg0r = omeg0
         ENDIF
         WRITE (iterm, 88020) omeg0r
88020    FORMAT(2X,1PE12.5,'  Omega')
      ENDIF

      IF (iexpan.NE.0) THEN
         WRITE (iterm, 88025) vexpan
88025    FORMAT(2X,1PE12.5,'  Expansion Velocity')
       ENDIF

      IF (ibound.EQ.7) THEN
         WRITE (iterm, 88035) hmaximum
88035    FORMAT(2X,1PE12.5,'  Maximum h (0 to disable)')
         WRITE (iterm, 88030) pext
88030    FORMAT(2X,1PE12.5,'  Const. External Pressure')
      ENDIF

      IF (ibound.EQ.8) THEN
         WRITE (iterm, 88036) deadbound
88036    FORMAT(2X,1PE12.5,'  Dead Particle Boundary')
         WRITE (iterm, 88037) fractan, fracradial, nstop, nfastd
88037    FORMAT(2X,2(1PE12.5,1X),I6,1X,I6,
     &        '  New Particle Velocities And Nstop')
      ENDIF
      IF (ibound/10.EQ.9) THEN
         WRITE (iterm, 88036) deadbound
         WRITE (iterm, 88038) fractan, fracradial, nshell, rshell
88038    FORMAT(2X,2(1PE12.5,1X),I6,1X,1PE12.5,
     &        '  New Part. Vels., Nshell, and Rshell')
      ENDIF
      IF (ibound.EQ.100) THEN
         WRITE (iterm, 88039) variation,phibound,flowrate,signorm,
     &        hoverr,partm
88039        FORMAT(2X,6(1PE12.5,1X),
     &        '  Rad variation, Phi, flowrate, sigma, HoverR, partmass')
      ENDIF

      IF (iexf.EQ.5 .OR. iexf.EQ.6 .OR. iexf.EQ.7) THEN
         WRITE (iterm, 88040) xmass
88040    FORMAT(2X,1PE12.5,'  External Forces Mass')
      ENDIF
      IF (iexf.EQ.7) THEN
         WRITE (iterm, 88043) planetmass
88043    FORMAT(2X,1PE12.5,'  Planet Mass')
      ENDIF

      IF (iptmass.NE.0 .OR. nptmass.NE.0) THEN
         WRITE (iterm, 88042) hacc
88042    FORMAT(2X,1PE12.5,'  Outer Accretion Radius')
         WRITE (iterm, 88044) haccall
88044    FORMAT(2X,1PE12.5,'  Inner Accretion Radius')
      ENDIF

      IF (iptmass.GE.1) THEN
         WRITE (iterm, 88046) radcrit
88046    FORMAT(2X,1PE12.5,'  Critical Radius')
         WRITE (iterm, 88048) ptmcrit
88048    FORMAT(2X,1PE12.5,'  Critical Density')

      ENDIF

      rzero = 0.0
      rminus1 = -1.0
      IF (ibound.EQ.1 .OR. ibound.EQ.3 .OR. ibound.EQ.8 .OR. 
     &                                           ibound/10.EQ.9) THEN 
         WRITE (iterm, 88050) rmind, rmax,
     &                        xmin, xmax, ymin, ymax, zmin, zmax
      ELSEIF (ibound.EQ.2) THEN
         WRITE (iterm, 88050) rmind, rcyl, 
     &                        xmin, xmax, ymin, ymax, zmin, zmax
      ELSE
         IF (rmax.EQ.rcyl) THEN
            WRITE (iterm, 88050) rminus1, rmax,
     &                        xmin, xmax, ymin, ymax, zmin, zmax
         ELSE
            WRITE (iterm, 88050) rmind, rcyl, 
     &                        xmin, xmax, ymin, ymax, zmin, zmax
         ENDIF
      ENDIF
88050  FORMAT(8(1X,1PE12.5))
c
c--Check for consistency
c
      CALL chekopt

      IF (idebug(1:7).EQ.'wrinsph') THEN
         WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc, iener,
     &                         ibound, damp, varsta
99004    FORMAT (1X, 7(I2,1X), 2(E12.5,1X), 1X, A7)
         WRITE (iprint, 99005) file1, ipos, nstep
99005    FORMAT (1X, A7, 1X, I4, 1X, I4)
      ENDIF

      CLOSE (iterm)

#ifdef MPI
      ENDIF
#endif

      RETURN
      END
