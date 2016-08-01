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
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/accstat'
      INCLUDE 'COMMONS/rtstat'

      REAL*8 angmom

      CHARACTER*35 var
      CHARACTER*7 where, where2

      DATA where2/'header'/

      nparttot = npart
      nactivetot = nactive
      fmas1tot = fmas1
      fmas2tot = fmas2
c      n1tot = n1
c      n2tot = n2
      nptmasstot = nptmass
#ifdef MPI
      CALL MPI_REDUCE(npart,nparttot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
c      CALL MPI_REDUCE(n1,n1tot,1,MPI_INTEGER,MPI_SUM,0,
c     &     MPI_COMM_WORLD,ierr)
c      CALL MPI_REDUCE(n2,n2tot,1,MPI_INTEGER,MPI_SUM,0,
c     &     MPI_COMM_WORLD,ierr)
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
      WRITE (iprint, 99003) idim, nparttot, nactivetot, 
     &        (fmas1tot + fmas2tot), n1, fmas1tot, n2, fmas2tot
99003 FORMAT (/, ' idim                           : ', I10, /, 
     &        ' Total number of particles used : ', I10, /,
     &        ' Number of active particles     : ', I10, /,
     &        ' Total mass                     :', 1PE12.4, /,
     &        ' Distribution : object 1 :  number of particles : ', I10,
     &        /, '                            mass of object      : ',
     &        1PE12.4, /,
     &        '                object 2 :  number of particles : ', I10,
     &        /, '                            mass of object      : ',
     &        1PE12.4, /)

c
c--Calculate tff
c
      tcomp = SQRT((3 * pi) / (32 * rhozero))
      tff = tcomp * utime
      WRITE (iprint, 99004) tff, tcomp
99004 FORMAT(/,' The free fall time is  : ',1PE14.7,/,
     &            ' in computational units : ',1PE14.7)

      IF (where(1:6).NE.'newrun') THEN
c
c--Options of the code
c
         WRITE (iprint, 99005) encal, igrp, igphi, ifsvi, ifcor, 
     &           ibound, iexf, iener, ichoc, iexpan, damp
  
99005    FORMAT (/, ' The following options were set :', /,
     &           ' energy calculation     : ', A8, /,
     &           ' pressure gradients     : ', I8,
     &           '   self-gravity     : ', I8, /,
     &           ' artificial viscosity   : ', I8,
     &           '   coriolis         : ', I8, /,
     &           ' boundary type          : ', I8,
     &           '   external force   : ', I8, /,
     &           ' energy conservation    : ', I8,
     &           '   shock heating    : ', I8, /,
     &           ' homologous exp.        : ', I8,
     &           '   general damping  : ', F8.1, /)
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
99533       FORMAT (/,' Boundary type: planet/disc section ',/,
     &           '  radius variation: ',1PE12.4,/,
     &           '  phi boundary    : ',1PE12.4,/,
     &           '  flowrate        : ',1PE12.4,/,
     &           '  sigma normal.   : ',1PE12.4,/,
     &           '  disc H/R        : ',1PE12.4,/,
     &           '  planet mass     : ',1PE12.4,/,
     &           '  planet radius   : ',1PE12.4,/)
         ELSEIF (ibound.EQ.104) THEN 
            WRITE (iprint,99536) variation,phibound,
     &           signorm,hoverr
99536       FORMAT (/,' Boundary type: disc section/periodic phi  ',/,
     &           '  radius variation: ',1PE12.4,/,
     &           '  phi boundary    : ',1PE12.4,/,
     &           '  sigma normal.   : ',1PE12.4,/,
     &           '  disc H/R        : ',1PE12.4,/)
         ELSEIF (ibound.EQ.102 .OR. ibound.EQ.103) THEN
            IF (nptmass.GT.0) THEN
               WRITE (iprint,99535) rcyl,rmind,
     &              signorm,sdprof,tprof,hoverr

               DO i = 1, nptmass
                  WRITE (iprint, 99537) i, xyzmh(4,listpm(i)), i,
     &                 xyzmh(5,listpm(i))
               ENDDO
               IF (initialptm.EQ.5) WRITE (iprint, 99538) uniformtslim
            
99535          FORMAT (/,' Boundary type: whole disc  ',/,
     &              '  outer radius           : ',1PE12.4,/,
     &              '  inner radius           : ',1PE12.4,/,
     &              '  sigma normal           : ',1PE12.4,/,
     &              '  disc density profile   : ',1PE12.4,/,
     &              '  disc temp profile      : ',1PE12.4,/,
     &              '  disc H/R               : ',1PE12.4,/)
99537          FORMAT (/,'  planet (',I2,') mass     : ',1PE12.4,/,
     &              '  planet (',I2,') radius   : ',1PE12.4,/)
99538          FORMAT (/,' uniform ts rad  : ',1PE12.4,/)

            ELSEIF (initialptm.EQ.0) THEN
               WRITE (iprint,99534) rcyl,rmind,
     &              signorm,sdprof,tprof,hoverr,planetmass,rplanet
99534          FORMAT (/,' Boundary type: whole disc  ',/,
     &           '  outer radius           : ',1PE12.4,/,
     &           '  inner radius           : ',1PE12.4,/,
     &           '  sigma normal.          : ',1PE12.4,/,
     &           '  disc density profile   : ',1PE12.4,/,
     &           '  disc temp profile      : ',1PE12.4,/,
     &           '  disc H/R               : ',1PE12.4,/,
     &           '  planet mass            : ',1PE12.4,/,
     &           '  planet radius          : ',1PE12.4,/)
            ENDIF

            IF (imigrate.EQ.1) THEN
               WRITE (iprint, 11001) pmrate/(1.49598E13/3.1536E13)/
     &              (utime/udist), pmrate*rorbitmax/(1.49598E13/
     &              3.1536E13)/(utime/udist)*(utime/3.1536E13)
11001          FORMAT ('Imposed migration rate (AU/Myr)  :', 1PE12.4,/,
     &                 'Total migration distance (AU)    :', 1PE12.4,/)
            ENDIF
         ENDIF
c
c--Planetesimal output for Ben Ayliffe's explicit single-size 
c     implementation
c
         IF (iplanetesimals.GT.0) THEN
            WRITE (iprint, 11002) r_planetesimals(1)*udist/1.0E5,
     &           rho_planetesimals(1)*udens
11002       FORMAT (' Planetesimal radius (km) : ', 1PE12.4,/,
     &           ' Planetesimal density (g/cm^3) : ', 1PE12.4,/)
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
            WRITE (iprint, 99012) 1.0, gam, rhocrt, gam, gamdh, rhocrt2,
     &             gamdh, gamah, rhocrt3
99012       FORMAT (/,' Critical densities for changing gamma ', /,
     &      '           from ',G10.3,'to ',G10.3,': ', 1PE12.3, /,
     &      '           from ',G10.3,'to ',G10.3,': ', 1PE12.3, /,
     &      '           from ',G10.3,'to ',G10.3,': ', 1PE12.3, //)
         ELSEIF (encal.EQ.'x') THEN
            rp1 = rhophys1 * udens
            rc1 = rhochange1 * udens
            rc2 = rhochange2 * udens
            WRITE (iprint, 99040) gamphys1, rp1, gamphys1, 
     &             gamphys2, rc1, gamphys2, gamphys3, rc2
99040       FORMAT (/,' Physical equation of state ', /,
     &      '           from 1.00    to ',G10.3,': ', 1PE12.3, /,
     &      '           from ',G10.3,'to ',G10.3,': ', 1PE12.3, /,
     &      '           from ',G10.3,'to ',G10.3,': ', 1PE12.3, //)
         ENDIF
c
c--Print out constants used for integration
c
         WRITE (iprint, 99013) alpha, beta, acc, tol, tolptm, tolh
99013    FORMAT (/,' Numerical constants used in this run :', /,
     &           ' artificial viscosity  alpha      : ', 1PE12.3, /,
     &           '                       beta       : ', 1PE12.3, /,
     &           ' binary tree accuracy param.      : ', 1PE12.3, /,
     &           ' RK2 tolerance - gas              : ', 1PE12.3, /,
     &           '               - point masses     : ', 1PE12.3, /,
     &           '               - smoothing length : ', 1PE12.3, /)
         IF (encal.EQ.'r') THEN
            WRITE (iprint, 99113) tolerance_rt
99113       FORMAT (' Radiative transfer tolerance    : ',1PE12.3, /)
            WRITE (iprint, 99613) opdenom
99613       FORMAT (' Opacity Denominator             : ',1PE12.3)
            WRITE (iprint, 99614) iopmodel
99614       FORMAT (' Opacity model                   :   ',I3)
            WRITE (iprint, 99615) metallicity
99615       FORMAT (' Metallicity                     : ',1PE12.3, /)
            IF (idustRT.GT.0) THEN
               WRITE (iprint, 99620)
               WRITE (iprint, 99622) iDRT_cosmic_ray
               WRITE (iprint, 99624) iDRT_photoelectric
               WRITE (iprint, 99626) iH2formation_heating
               WRITE (iprint, 99628) iDRT_line
               WRITE (iprint, 99630) iCHEM_depletion
               WRITE (iprint, 99632) iDRT_gasdust
               WRITE (iprint, 99634) iISR_MIR
               WRITE (iprint, 99636) iISR_Draine
               IF (iatomic) THEN
                  WRITE (iprint, 99650)
               ELSE
                  WRITE (iprint, 99651)
               ENDIF
99620          FORMAT(' Diffuse ISM model used')
99622          FORMAT('    Heating: Cosmic ray                   : ',L1)
99624          FORMAT('             Photoelectric                : ',L1)
99626          FORMAT('             H2 formation on dust grains  : ',L1)
99628          FORMAT('    Cooling: Molecular line (Goldsmith01) : ',L1)
99630          FORMAT('             Molecular depletion allowed  : ',L1)
99632          FORMAT('    Gas-Dust Thermal Coupling             : ',I1)
99634          FORMAT('    Interstellar Radiation Field: Mid-IR  : ',L1)
99636          FORMAT('                                Draine UV : ',L1)
99650          FORMAT('    Initial conditions : Fully Atomic', //)
99651          FORMAT('    Initial conditions : Fully Molecular', //)
            ELSE
               WRITE (iprint, 99699)
99699          FORMAT(' Diffuse ISM model not used', //)
            ENDIF
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
99213       FORMAT(' Magnetic fields are ON, evolving B/rho',/)
99214       FORMAT(' Magnetic fields are ON, evolving B',/)
99215       FORMAT(' Magnetic fields are ON, using Euler potentials',/)
99216       FORMAT(' Magnetic fields are ON, using vector potential',/)
99217       FORMAT(' WARNING: Magnetic fields are ON',
     &             ', but unknown variable',/)

#ifdef NONIDEAL
            WRITE (iprint,99220)
99220       FORMAT(' Non-ideal MHD is being used',/)
#endif

            IF (imhdevol.EQ.4) THEN
               WRITE (iprint,99218)
               WRITE (iprint,99219) ndivBsubcycles, overclean_fac
99218          FORMAT(' Magnetic field has hyperbolic cleaning ON')
99219          FORMAT(' Cleaning done using ',I2,
     &              ' subcycles and overcleaning factor ',1F6.3,//)
            ENDIF

            IF (ABS(iresist).EQ.1) THEN
               WRITE(iprint,99221) etamhd
            ELSEIF (ABS(iresist).GE.2) THEN
               WRITE(iprint,99222) ABS(iresist)
            ENDIF
99221       FORMAT(' Physical resistivity ON, using constant eta =',
     &             1PE12.3,//)
99222       FORMAT(' Physical resistivity type ',I1,' ON',
     &           ' function of density and temperature',//)
            IF (iresist.LT.0) THEN
               WRITE(iprint,99223)
            ELSEIF (iresist.GT.0) THEN
               WRITE(iprint,99224)
            ENDIF
99223       FORMAT(' Physical resistivity evolved implicitly',//)
99224       FORMAT(' Physical resistivity evolved explicitly',//)

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
     &           '  minimum creation radius       : ', 1PE12.3, /,
     &           '  creation density (in rhozero) : ', 1PE12.3, /,
     &           '                (in code units) : ', 1PE12.3, /,
     &           '            (in physical units) : ', 1PE12.3, /)
99015    FORMAT (' Point mass creation NOT ALLOWED')

         IF (iptmass.NE.0.OR.nptmass.NE.0) THEN
            IF (iptintree.EQ.0) THEN
               WRITE (iprint,99146) 
99146          FORMAT (' Point masses done in GPTI')
            ELSEIF (iptintree.EQ.1) THEN
               WRITE (iprint,99147) 
99147          FORMAT (' Point masses done in GFORSPT/TREE')
            ELSEIF (iptintree.EQ.2) THEN
               WRITE (iprint,99148) 
99148          FORMAT (' Point masses done entirely in TREE')
            ELSE
               WRITE (iprint,99149)
99149          FORMAT (' ERROR - Ptmasses in header ')
               CALL quit(0)
            ENDIF

            WRITE(iprint,99016) nptmasstot
99016       FORMAT (' Number of point masses         : ', 8X, I4)

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
99017       FORMAT ('  Outer accretion radius        : ', 1PE12.3)

            IF (haccall.GT.0.) THEN
               WRITE (iprint,99018) haccall
            ELSE
               CALL error(where2,2)
            END IF
99018       FORMAT ('  Inner accretion radius        : ', 1PE12.3)

            WRITE (iprint,99019) iptsoft
            IF (iptsoft.NE.0) WRITE (iprint,99020) ptsoft
99019       FORMAT ('  Gravity ptmass softening      : ', 10X, I2)
99020       FORMAT ('  Softening ptmass parameter    : ', 1PE12.3)

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

c
c--Zero the RT statistics counters
c
#ifdef RT
      nintegfail = 0
      nintegsuc  = 0 
      nrtnc1     = 0 
      nrtnc2     = 0
      nrtnca     = 0
      nrtmaxit   = 0
      nintegtot  = 0
#endif
c
c--Zero the pointmass counters
c
      iaccf1  = 0
      iaccf2  = 0
      iaccf3  = 0
      iaccf4  = 0
      iaccf5  = 0
      iaccf6  = 0
      iaccph1 = 0
      iaccph2 = 0
      iaccsuc = 0
      iaccmg  = 0 
#ifdef MPIALL
      ENDIF
#endif

      RETURN
      END
