      SUBROUTINE header(where)
c************************************************************
c                                                           *
c  This routine writes on first page of listing the value   *
c     of all variable defined at the start of the run.      *
c                                                           *
c************************************************************
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/vargam'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/physeos'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/h2'

      REAL*8 angmom

      CHARACTER*35 var
      CHARACTER*7 where, where2

      DATA where2/'header'/

      nparttot = npart
      nactivetot = nactive
      fmas1tot = fmas1
      fmas2tot = fmas2
      n1tot = n1
      n2tot = n2
      nptmasstot = nptmass
#ifdef MPI
      CALL MPI_REDUCE(npart,nparttot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(n1,n1tot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(n2,n2tot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(fmas1,fmas1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(fmas2,fmas2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)

      CALL MPI_ALLREDUCE(nptmass,nptmasstot,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
#endif
#ifdef MPIALL
      IF (iproc.EQ.0) THEN
#endif
c
c--Scaling factors
c
      CALL scaling(gt, rscale, drdt, dlnrdt)
c
c--Write units
c
      angmom = umass*dble(udist**2)/dble(utime)
      velo = udist/utime
      WRITE (iprint, 99001) umass, udist, udens, utime, velo, uergg,
     &                      angmom,umagfd
99001 FORMAT (//, ' The computations are done in the following units',
     &        /, ' units of :  mass       :', 1PD12.4,
     &        '   distance    :', 1PE12.4, /,
     &        '             density    :', 1PE12.4, '   time        :',
     &        1PE12.4, /, '             velocity   :', 1PE12.4,
     &        '   energy/mass :', 1PE12.4, /,
     &        '             ang. mom.  :', 1PE12.4,
     &        '   mag field   :', 1PE12.4, //)
c
c--Write options
c
      IF (varsta.EQ.'entropy') THEN
         var = 'specific entropy'
      ELSE
         var = 'specific internal energy'
      ENDIF
      WRITE (iprint, 99002) var
99002 FORMAT (' Variable of state used : density and ', A35, /)
c
c--Mass fractions
c
      WRITE (iprint, 99003) nparttot, nactivetot, (fmas1tot + fmas2tot),
     &        n1tot, fmas1tot, n2tot, fmas2tot
99003 FORMAT (/, ' Total number of particles used : ', I8, /,
     &        ' Number of active particles     : ', I8, /,
     &        ' Total mass                     :', 1PE12.4, /,
     &        ' Distribution : object 1 :  number of particles : ', I8,
     &        /, '                            mass of object      : ',
     &        1PE12.4, /,
     &        '                object 2 :  number of particles : ', I8,
     &        /, '                            mass of object      : ',
     &        1PE12.4, /)

c
c--Calculate tff
c
      tcomp = SQRT((3 * pi) / (32 * rhozero))
      tff = tcomp * utime
      WRITE (iprint, 99004) tff, tcomp
99004 FORMAT(/,' The free fall time is : ',1PE14.7,/,
     &            ' in computational units: ',1PE14.7)

      IF (where(1:6).NE.'newrun') THEN
c
c--Options of the code
c
         WRITE (iprint, 99005) encal, igrp, igphi, ifsvi, ifcor, 
     &           ibound, iexf, iener, ichoc, iexpan, damp
  
99005    FORMAT (/, ' The following options were set :', /,
     &           ' energy calculation     : ', A1, /,
     &           ' pressure gradients     : ', I2,
     &           '   self-gravity     : ', I2, /,
     &           ' artificial viscosity   : ', I2,
     &           '   coriolis         : ', I2, /,
     &           ' boundary type          : ', I3,
     &           '   external force   : ', I2, /,
     &           ' energy conservation    : ', I2,
     &           '   shock heating    : ', I2, /,
     &           ' homologous exp.        : ', I2,
     &           '   general damping  : ', F3.1, /)
c
c--Gravity calculations
c
         IF (igphi.EQ.1) THEN
            IF (igrape.EQ.1) THEN
               WRITE (iprint, 98001)
98001          FORMAT (' Gravity calculated using: GRAPE board')
            ELSEIF (igrape.EQ.0) THEN
               WRITE (iprint, 98002) 
98002          FORMAT (' Gravity calculated using: Binary Tree')
            ELSE
               CALL error(where2, 1)
            ENDIF
 
            IF (isoft.EQ.1) THEN
               WRITE (iprint, 98003) psoft
98003          FORMAT ('                           Softening = ',
     &              1PE12.5)
            ELSEIF (isoft.EQ.0) THEN
               WRITE (iprint, 98004)
98004          FORMAT ('                           Kernel Softening')
            ELSE
               CALL error(where2, 2)
            ENDIF

            IF (igrape.EQ.0) THEN
               WRITE (iprint, 98303) ipartialrevtree
98303          FORMAT (' Partial tree revision: ',I1,/)
            ENDIF
         ENDIF
c
c--Rotation of expansion velocity
c
         IF (iexpan.NE.0) THEN
            WRITE (iprint, 99006) drdt*udist/utime
99006       FORMAT (/,
     &           ' Calculations done in a frame homologously expanding'
     &           , ' at ', 1PE12.5, ' cm/s', /)
         ENDIF
         IF (ifcor.GT.0) THEN
            WRITE (iprint, 99007) omeg0/utime
99007       FORMAT (/,
     &           ' Calculations done in a frame in uniform rotation',
     &           ' at ', 1PE12.5, ' /s ', /)
         ENDIF
c
c--Boundaries
c
         IF (ibound.EQ.1) THEN
            WRITE (iprint, 99008) xmin, xmax, ymin, ymax, zmin, zmax
99008       FORMAT (/,' Boundary type : reflective ', /,
     &              ' position :  cartesian    :  xmin : ', F7.3,
     &              '  xmax :', F7.3, /,
     &              '                             ymin : ', F7.3,
     &              '  ymax :', F7.3, /,
     &              '                             zmin : ', F7.3,
     &              '  zmax :', F7.3, /)
         ELSEIF (ibound.EQ.2) THEN
            WRITE (iprint, 99009) rcyl, zmin, zmax
99009       FORMAT(/,' Boundary type : cylindrical reflective ',/,
     &                 '           rcyl : ', F7.3, /,
     &                 '           zmin : ', F7.3, /,
     &                 '           zmax : ', F7.3, /)
            IF (rmind.NE.0) WRITE (iprint, 99010) rmind
99010       FORMAT (/,' accretion disk inner boundary :',F7.3,/)
         ELSEIF (ibound.EQ.3) THEN 
            WRITE (iprint, 99011) rmax
99011       FORMAT (/,' Boundary type : spherical reflective ',/,
     &                  '          rmax : ', F7.3, /)  
         ELSEIF (ibound.EQ.7) THEN 
            WRITE (iprint,99501) pext
99501       FORMAT (/,' Boundary type: constant pressure',/,
     &           '          pressure: ',1PE12.4,/)
         ELSEIF (ibound.EQ.8 .OR. ibound/10.EQ.9) THEN 
            WRITE (iprint,99502) deadbound, specang, fractan, fracradial
99502       FORMAT (/,' Boundary type: dead particle',/,
     &           '            radius: ',1PE12.4,/,
     &           '   spec. ang. mom.: ',1PE12.4,/,
     &           '  frac. tangential: ',1PE12.4,/,
     &           '  frac. radial    : ',1PE12.4,/)
         ELSEIF (ibound.EQ.100) THEN 
            WRITE (iprint,99533) variation,phibound,flowrate,
     &           signorm,hoverr,planetmass,rplanet
99533       FORMAT (/,' Boundary type: planet/disc  ',/,
     &           '  radius variation: ',1PE12.4,/,
     &           '  phi boundary    : ',1PE12.4,/,
     &           '  flowrate        : ',1PE12.4,/,
     &           '  sigma normal.   : ',1PE12.4,/,
     &           '  disc H/R        : ',1PE12.4,/,
     &           '  planet mass     : ',1PE12.4,/,
     &           '  planet radius   : ',1PE12.4,/)
         ELSEIF (ibound.EQ.102) THEN
            IF (initialptm.EQ.5) THEN
               WRITE (iprint,99535) rcyl,rmind,
     &              signorm,hoverr,xyzmh(4,listpm(1)),rplanet,
     &              uniformtslim
99535          FORMAT (/,' Boundary type: whole disc  ',/,
     &              '  outer radius    : ',1PE12.4,/,
     &              '  inner radius    : ',1PE12.4,/,
     &              '  sigma normal.   : ',1PE12.4,/,
     &              '  disc H/R        : ',1PE12.4,/,
     &              '  planet mass     : ',1PE12.4,/,
     &              '  planet radius   : ',1PE12.4,/,
     &              '  uniform ts rad  : ',1PE12.4,/)
            ELSEIF (initialptm.EQ.0) THEN
               WRITE (iprint,99534) rcyl,rmind,
     &              signorm,hoverr,planetmass,rplanet
            ELSE
               WRITE (iprint,99534) rcyl,rmind,
     &              signorm,hoverr,xyzmh(4,listpm(1)),rplanet
99534          FORMAT (/,' Boundary type: whole disc  ',/,
     &              '  outer radius    : ',1PE12.4,/,
     &              '  inner radius    : ',1PE12.4,/,
     &              '  sigma normal.   : ',1PE12.4,/,
     &              '  disc H/R        : ',1PE12.4,/,
     &              '  planet mass     : ',1PE12.4,/,
     &              '  planet radius   : ',1PE12.4,/)
            ENDIF
         ENDIF
c
c--Critical densities for variable gamma
c
         WRITE (iprint, 99504) gamma
99504    FORMAT (' gamma                : ', 1PE12.3)
         WRITE (iprint, 99503) gmw
99503    FORMAT (/,' Mean molecular weight: ', 1PE12.3, /)

         IF (encal.EQ.'v') THEN
            rhocrt = rhocrit * udens
            rhocrt2 = rhocrit2 * udens
            rhocrt3 = rhocrit3 * udens
            WRITE (iprint, 99012) gam, rhocrt, gam, gamdh, rhocrt2, 
     &             gamdh, gamah, rhocrt3
99012       FORMAT (/,' Critical densities for changing gamma ', /,
     &      '           from 1.00    to ',G8.3,': ', 1PE12.3, /,
     &      '           from ',G8.3,'to ',G8.3,': ', 1PE12.3, /,
     &      '           from ',G8.3,'to ',G8.3,': ', 1PE12.3, //)
         ELSEIF (encal.EQ.'x') THEN
            rp1 = rhophys1 * udens
            rc1 = rhochange1 * udens
            rc2 = rhochange2 * udens
            WRITE (iprint, 99040) gamphys1, rp1, gamphys1, 
     &             gamphys2, rc1, gamphys2, gamphys3, rc2
99040       FORMAT (/,' Physical equation of state ', /,
     &      '           from 1.00    to ',G8.3,': ', 1PE12.3, /,
     &      '           from ',G8.3,'to ',G8.3,': ', 1PE12.3, /,
     &      '           from ',G8.3,'to ',G8.3,': ', 1PE12.3, //)
         ENDIF
c
c--Print out constants used for integration
c
         WRITE (iprint, 99013) alpha, beta, acc, tol, tolptm, tolh
99013    FORMAT (/,' Numerical constants used in this run :', /,
     &           ' artificial viscosity  alpha     : ', 1PE12.3, /,
     &           '                       beta      : ', 1PE12.3, /,
     &           ' binary tree accuracy param.     : ', 1PE12.3, /,
     &           ' RK2 tolerance - gas             : ', 1PE12.3, /,
     &           '               - point masses    : ', 1PE12.3, /,
     &           '               - smoothing length: ', 1PE12.3, //)
         IF (encal.EQ.'r') THEN
            WRITE (iprint, 99113) tolerance
99113       FORMAT (' Radiative transfer tolerance    : ',1PE12.3, //)
            WRITE (iprint, 99613) opdenom
99613       FORMAT (' Opacity Denominator             : ',1PE12.3, //)
         ELSEIF (encal.EQ.'c') THEN
            WRITE (iprint, 99313) coolingrate
99313       FORMAT (' Cooling Rate                    : ',1PE12.3, //)
         ELSE
            WRITE (iprint, *)
         ENDIF
c
c--Simon's chemistry
c
         IF (iener.EQ.4) THEN
         WRITE (iprint,99401)dphot/pc,dlq/pc,dchem/pc,
     &           Rconst,sigma,rate_diss0,rate_cr,phrates
99401    FORMAT('Heating/cooling and H2 formation is ON',/,
     &           ' photoionisation distance (pc)       : ', 1PE12.3, /,
     &           ' dist. for heating/cooling (pc)      : ', 1PE12.3, /,
     &           ' dist. for HI col. density (pc)      : ', 1PE12.3, /,
     &           ' Rconst, rate of formation parameter : ', 1PE12.3, /,
     &           ' extinction (cm-2)                   : ', 1PE12.3, /,
     &           ' unshielded dissociation rate        : ', 1PE12.3, /,
     &           ' cosmic ray rate                     : ', 1PE12.3, /,
     &           ' dissociation parameter for CO       : ', 1PE12.3, //)
          ENDIF
c
c--Magnetic field details
c
         IF (imhd.EQ.idim) THEN
            IF (varmhd.EQ.'Brho') THEN
               WRITE(iprint,99213)
            ELSEIF (varmhd.EQ.'Bvol') THEN
               WRITE(iprint,99214)            
            ELSEIF (varmhd.EQ.'eulr') THEN
               WRITE(iprint,99215)
            ELSEIF (varmhd.EQ.'vecp') THEN
               WRITE(iprint,99216)
            ELSE
               WRITE(iprint,99217)
            ENDIF
99213       FORMAT(' Magnetic fields are ON, evolving B/rho',//)
99214       FORMAT(' Magnetic fields are ON, evolving B',//)
99215       FORMAT(' Magnetic fields are ON, using Euler potentials',//)
99216       FORMAT(' Magnetic fields are ON, using vector potential',//)
99217       FORMAT(' WARNING: Magnetic fields are ON',
     &             ', but unknown variable',//)
     
            IF (iresist.EQ.1) THEN
               WRITE(iprint,99221) etamhd
            ELSEIF (iresist.EQ.2) THEN
               WRITE(iprint,99222)
            ENDIF
99221       FORMAT(' Physical resistivity ON, using constant eta =',
     &             1PE12.3,//)
99222       FORMAT(' Physical resistivity ON, function of density',//)

         ENDIF
c
c--Print out massive point mass details
c
         IF (iptmass.NE.0) THEN
            WRITE(iprint,99014) iptmass, radcrit, ptmcrit, rhocrea,
     &                          rhocrea*udens
         ELSE
            WRITE(iprint,99015)
         ENDIF
99014    FORMAT (' Point mass creation ALLOWED, type ', I2, /,
     &           '  minimum creation radius      : ', 1PE12.3, /,
     &           '  creation density (in rhozero): ', 1PE12.3, /,
     &           '                (in code units): ', 1PE12.3, /,
     &           '            (in physical units): ', 1PE12.3, /)
99015    FORMAT (' Point mass creation NOT ALLOWED')

         IF (iptmass.NE.0.OR.nptmass.NE.0) THEN
            IF (iptintree.EQ.0) THEN
               WRITE (iprint,99146) 
99146          FORMAT (' Point masses done in GPTALL')
            ELSEIF (iptintree.EQ.1) THEN
               WRITE (iprint,99147) 
99147          FORMAT (' Point masses done in GFORSPT/TREE')
            ELSEIF (iptintree.EQ.2) THEN
               WRITE (iprint,99148) 
99148          FORMAT (' Point masses done entirely in TREE')
            ELSE
               WRITE (iprint,99149)
99149          FORMAT (' ERROR - Ptmasses in header ')
               CALL quit
            ENDIF

            WRITE(iprint,99016) nptmasstot
99016       FORMAT (' Number of point masses        : ', I4)

            IF (iaccevol.EQ.'v') THEN
               WRITE (iprint,99150) accfac
99150          FORMAT ('  Accretion radii: VARIABLE ROCHE ',1PE12.5)
            ELSEIF (iaccevol.EQ.'s') THEN
               WRITE (iprint,99151) accfac
99151       FORMAT ('  Accretion radii: VARIABLE SEPARATION ',1PE12.5)
            ELSE
               WRITE (iprint,99152)
99152          FORMAT ('  Accretion radii: FIXED ')
            ENDIF

            IF (hacc.GT.0.) THEN
               WRITE (iprint,99017) hacc
            ELSE
               CALL error(where2,1)
            END IF
99017       FORMAT ('  Outer accretion radius       : ', 1PE12.3)

            IF (haccall.GT.0.) THEN
               WRITE (iprint,99018) haccall
            ELSE
               CALL error(where2,2)
            END IF
99018       FORMAT ('  Inner accretion radius       : ', 1PE12.3)

            WRITE (iprint,99019) iptsoft
            IF (iptsoft.NE.0) WRITE (iprint,99020) ptsoft
99019       FORMAT ('  Gravity ptmass softening     : ', I2)
99020       FORMAT ('  Softening ptmass parameter   : ', 1PE12.3)

c            IF (nptmass.NE.0) THEN
c               DO i = 1, nptmass
c                  WRITE (iprint,99021) i,iphase(listpm(i)),
c     &                 xyzmh(5,listpm(i))
c99021             FORMAT (' Point mass: ',I4, ' type: ',I1,
c     &                 ' hacc: ',1PE12.3)
c               END DO
c            ENDIF
         ENDIF

      ENDIF
c
c--Smoothing length max min
c
      IF (hmin.NE.0.0) THEN
         WRITE (iprint, 99114) hmin
99114        FORMAT (/,' MINIMUM SMOOTHING LENGTH = ', 1PE12.3, //)
      ENDIF
      IF (hmaximum.NE.0.0) THEN
         WRITE (iprint, 99115) hmaximum
99115        FORMAT (/,' MAXIMUM SMOOTHING LENGTH = ', 1PE12.3, //)
      ENDIF
c
c--Write name of file used
c
      WRITE (iprint, 99025) file1
99025 FORMAT (//,' Name of input file : ', A7, //)

      CALL FLUSH (iprint)

#ifdef MPIALL
      ENDIF
#endif

      RETURN
      END
