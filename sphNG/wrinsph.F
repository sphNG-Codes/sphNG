      SUBROUTINE wrinsph
c************************************************************
c                                                           *
c  This subroutine writes a new file inname from parameters *
c     taken from setpart                                    *
c                                                           *
c************************************************************
#ifdef MPIALL
#include "mpi_sup.h"
#endif
      INCLUDE 'idim'

#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

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
      INCLUDE 'COMMONS/prdrag'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/grain'
      INCLUDE 'COMMONS/tstop'
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
#ifdef MPICOPY
      IF (iproc.EQ.0) THEN
#endif
      OPEN (iterm, FILE=inname, STATUS='replace', ACTION='write')
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
      IF (imhd.EQ.idim) THEN
         WRITE (iterm, 10001) varmhd
      ELSE
         WRITE (iterm, 10001) 'Brho'
      ENDIF
10001 FORMAT(A4)
c
c--Write options
c
      WRITE (iterm, 99006) encal
99006 FORMAT (A1)
      IF (encal.EQ.'r') THEN
         WRITE (iterm, 99008) tolerance_rt, boundtemp, bounddens,
     &        opdenom, iopmodel, metallicity
99008    FORMAT(2X,4(1PE12.5, 1X),I2,1X,1PE12.5,
     &    ' Implicit radiative transfer',
     &    ' tolerance, boundtemp, bounddens, opacity divisor,'
     &    ' opacity model, metallicity')
      ENDIF
      IF (idustRT.GT.0) THEN
         IF (iatomic) THEN
            WRITE (iterm, 99006) 'a'
         ELSE
            WRITE (iterm, 99006) 'm'
         ENDIF
      ENDIF

#ifdef MPI
      WRITE (iterm, 99016) mpitype, xlocationmpi, ylocationmpi,
     &     zlocationmpi, densmpi
99016 FORMAT (A1,4(1X,1PE16.9),'  Type of MPI decomposition ',
     &     '(s)phere, (d)isc, (c)artesian with origin')
#else
      WRITE (iterm, 99017)
99017 FORMAT ('d  0.0 0.0 0.0 0.0  Type of MPI decomposition ',
     &     '(s)phere, (d)isc, (c)artesian')
#endif

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
88001 FORMAT(2X,I2,'  Self-gravity (1=on, 0=off, -1=no back-reaction ',
     &     'on sinks from gas,etc)')
      IF (imhd.EQ.idim) THEN
         WRITE (iterm, 87002) ifsvi,alpha,beta,alphamin(2),alphamax(2)
         WRITE (iterm, 87003) ndivBsubcycles, overclean_fac, iresist, 
     &        etamhd
      ELSE
         WRITE (iterm, 87002) ifsvi,alpha,beta, 0., 0.
         WRITE (iterm, 87003) 0, 1.0, 0, 0.
      ENDIF
87002 FORMAT(2X,I2,4(1X,1F6.3),' Artificial Viscosity & Resistivity')
87003 FORMAT(2X,I2,1X,1F6.3,1X,I2,1X,1F6.3,'  MHD: DivB-Subcycles,',
     &     ' DivB-Overclean, Ohmic Resistivity ',
     &     '(1=const,2,3=func,4=NICIL; Negative is implicit),eta')
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
      WRITE (iterm, 88012) irec 
88012  FORMAT(2X,I4,'  Binary Dump Counter')
      WRITE (iterm, 88013) tmax
88013  FORMAT(2X,1PE12.5,'  Max Time (min) (Status Written)')
      WRITE (iterm, 88014) tstop  
88014  FORMAT(2X,1PE12.5,'  Max Dynamic Time')
      WRITE (iterm, 88015) dtmax
88015  FORMAT(2X,1PE14.7,'  Max Timestep')
      WRITE (iterm, 88016) dtini
88016  FORMAT(2X,1PE12.5,'  Initial Timestep')

      IF (ifcor.NE.0) THEN
         omeg0r = omeg0/utime
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
     &        hoverr,partm,centralmass
88039    FORMAT(2X,7(1PE11.5,1X),
     &        '  Rad variation, Phi, flowrate, sigma, HoverR, partmass,'
     &        ,' centralmass (for u)')
      ENDIF

      IF (ibound.EQ.102 .OR. ibound.EQ.103 .OR. ibound.EQ.105) THEN
88071    FORMAT(3X, L1, ' t-profile switch')
         WRITE (iterm, 88071) use_tprof
         IF (use_tprof) WRITE (iterm, 88072) tprof
88072    FORMAT(2X,1PE12.5, '  Temperature profile')
      ENDIF

      IF (ibound.EQ.104 .OR. ibound.EQ.105) THEN
         WRITE (iterm, 88041) variation,phibound,signorm,
     &        hoverr,partm,centralmass
88041    FORMAT(2X,6(1PE12.5,1X),
     &        '  Rad variation, Phi, sigma, HoverR, partmass,',
     &        ' centralmass (for u)')
      ENDIF

      IF (ibound.EQ.102 .OR. ibound.EQ.103) THEN
         WRITE (iterm, 88073) signorm,hoverr,partm,centralmass
88073    FORMAT(2X,4(1PE12.5,1X),'  sigma, HoverR, partmass',
     &        'centralmass (for u)')
         WRITE (iterm, 11074) imigrate
11074    FORMAT(2X, I2, ' prescribed migration')
         IF (imigrate.EQ.1) THEN
            WRITE (iterm, 11075) pmrate/(1.49598E13/3.1536E13)/
     &           (utime/udist)
11075       FORMAT(2X, 1PE12.5, ' applied migration rate (AU/Myr)')
            WRITE (iterm, 11076) rorbitmax*pmrate + rorbit_orig
11076       FORMAT(2X, 1PE12.5, ' stopping radius (code units)')
         ENDIF
      ENDIF

      IF (iexf.EQ.5 .OR. iexf.EQ.6 .OR. iexf.EQ.7 .OR. iexf.EQ.11) THEN
         WRITE (iterm, 88040) xmass
88040    FORMAT(2X,1PE12.5,'  External Forces Mass')
      ENDIF
      IF (iexf.EQ.7 .OR. 
     &     (iexf.EQ.5 .AND. (ibound.EQ.102 .OR. ibound.EQ.103)) ) THEN
         WRITE (iterm, 88100) numplanet
88100    FORMAT(2X,I2,'  Number of planet potentials')
         DO i = 1, numplanet
            WRITE (iterm, 88105) planetmass(i), i
         END DO
88105    FORMAT(2X,1PE15.8,'  Planet Mass ',I2)
         DO i = 1, numplanet
            WRITE (iterm, 88110) planetsemiaxis(i), i
         END DO
88110    FORMAT(2X,1PE15.8,'  Planet Semi-Major Axis ',I2)
         DO i = 1, numplanet
            WRITE (iterm, 88115) planetradius(i), i
         END DO
88115    FORMAT(2X,1PE12.5,'  Planet Radius ',I2)
         IF (iexf.EQ.7) THEN
            WRITE (iterm, 88200) iandrea
88200       FORMAT(2X,I2,'  Include planetesimal heating')
            WRITE (iterm, 88205) eccentricity_p
88205       FORMAT(2X,1PE12.5,'  Planetesimal eccentricity max')
            WRITE (iterm, 88210) inclination_p
88210       FORMAT(2X,1PE12.5,'  Planetesimal inclination max')
         ENDIF
      ENDIF

      IF ((ibound.EQ.102 .OR. ibound.EQ.103) .AND. initialptm.EQ.5) THEN
88300    FORMAT(2X,1PE12.5,'  Uniform Timestep radius')
         WRITE (iterm, 88300) uniformtslim
      ENDIF

      IF (iptmass.NE.0 .OR. nptmasstot.NE.0) THEN
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

      IF (encal.EQ.'c') THEN
         WRITE (iterm, 88049) coolingrate
88049    FORMAT(2X,1PE12.5,'  Cooling Rate for encal.EQ.c')
      ENDIF

      rzero = 0.0
      rminus1 = -1.0
      IF (ibound.EQ.1 .OR. ibound.EQ.3 .OR. ibound.EQ.8 .OR. 
     &                                           ibound/10.EQ.9) THEN 
         WRITE (iterm, 88050) rmind, rmax,
     &                        xmin, xmax, ymin, ymax, zmin, zmax
      ELSEIF (ibound.EQ.2 .OR. ibound.EQ.104) THEN
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
88050 FORMAT(8(1X,1PE12.5))
c
c--Parameters for planetesimals and dust
c
      IF (iplanetesimal_info.EQ.'p' .OR. iplanetesimal_info.EQ.'g'
     &     .OR. idustFluid.NE.0) THEN
         IF (idustFluid.NE.0) iplanetesimal_info = 'd'
         WRITE (iterm, 99006) iplanetesimal_info
11000    FORMAT (A1,'     Planetesimal or dust grain info (p/g/d)')
         IF (iplanetesimal_info.EQ.'p' .OR. iplanetesimal_info.EQ.'g')
     &        THEN
            WRITE (iterm, 11100) gasdrag
11100       FORMAT(2X,L1, '   Explicit planetesimal-gas drag (T/F)')
         ENDIF
         IF (idustFluid.NE.0 .OR. gasdrag) THEN
            WRITE (iterm, 11110) idragscheme
11110       FORMAT(2X,I2, '  Drag scheme (0,1,2,3,4)')
            IF (idragscheme.EQ.4) THEN
               WRITE (iterm, 11115) K_code
11115          FORMAT(2X,1PE12.5,'   Fixed drag coefficient')
            ENDIF
         ENDIF
         IF (gasdrag) THEN
            WRITE (iterm, 11120) prswitch
11120       FORMAT(2X, L1, '   PR-drag & rad. pres. switch')
            IF (prswitch) WRITE (iterm, 11130) pbeta
11130       FORMAT(2X,1PE12.5, ' Grain beta (Burns et al. 1979)')
         ENDIF
      ENDIF

      IF (iplanetesimal_info.EQ.'p' .AND. idragscheme.NE.4) THEN
         WRITE (iterm, 11140) n_planetesimal_types
11140    FORMAT(2X, I2, '  Number of planetesimal types')
         DO i = 1, n_planetesimal_types
            r_planetesimal = r_planetesimals(i)*udist
            rho_planetesimal = rho_planetesimals(i)*udens
            WRITE (iterm, 11150) r_planetesimal,rho_planetesimal
11150       FORMAT (4X,2(1X,1PE12.5),
     &           ' planetesimal radius (cm), density (g/cm^3)')
         END DO
      ENDIF

      IF (iplanetesimal_info.EQ.'g') THEN
         WRITE (iterm, 11160) igrowthmodel
11160    FORMAT(2X, I1, '   Growth model ',
     &        '(1:monodisperse 2:polydisperse)')
         WRITE (iterm, 11170) vfrag*udist/utime
11170    FORMAT(2X, 1PE12.5, ' Dust fragmentation velocity (cm/s)',
     &        '(negative value to switch frag off)')
         WRITE (iterm, 11180) rgrain_min*udist
11180    FORMAT(2X, 1PE12.5, ' Minimum grain size (cm)',
     &        '(negative value to remove minimum size)')
      ENDIF
c
c--Check for consistency
c
      CALL chekopt

      IF (idebug(1:7).EQ.'wrinsph') THEN
         WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc, iener,
     &                         ibound, damp, varsta
99004    FORMAT (1X, 7(I2,1X), 2(E12.5,1X), 1X, A7)
         WRITE (iprint, 99005) file1, irec, nstep
99005    FORMAT (1X, A7, 1X, I4, 1X, I4)
      ENDIF

      CLOSE (iterm)

#ifdef MPI
      ENDIF
#endif
#ifdef MPICOPY
      ENDIF
#endif

      RETURN
      END
