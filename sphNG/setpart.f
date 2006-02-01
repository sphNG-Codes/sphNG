       SUBROUTINE setpart
c************************************************************
c                                                           *
c  This subroutine handles the particle set up.             *
c                                                           *
c  The distribution can be   random  (idist=1)              *
c                            close-packed spheres (idist=2) *
c                            cubic lattice                  *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/new'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/latti'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/flag'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/setbin'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/crpart'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/Bzero'
      INCLUDE 'COMMONS/setlocal'

      CHARACTER*20 filevelx, filevely, filevelz
      CHARACTER*1 iok, iok2, iwhat, idens, ipres, icentral, irotatey
      CHARACTER*1 ien
c
c--Initialise time
c
      gt = 0.0

      CALL ktable
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine setpart')
99004 FORMAT (A1)

      CALL getime(ih, im, is, fhour)
      iseed = -1 * ((ih * 400) + (im * 17) + is)
      iseed = -6485
      WRITE(*,99500) iseed
99500 FORMAT(' The random seed is ',I5)
      rnd = ran1(iseed)
c
c--Point masses 
c
      npart = 0
      iaccevol = 'f'
      WRITE (*, 99085) 
99085 FORMAT (' Do you want point masses? (y/n)')
      READ (*, 99004) iok
      IF (iok.EQ.'y') THEN
         totptmass = 0.
         WRITE (*, 99086) 
99086    FORMAT (' Enter number of point masses')
         READ (*,*) numpt

         WRITE (*, 99091) 
99091    FORMAT (' Enter type of point mass (1,2,3,4)')
         READ (*,*) ipttype
         initialptm = ipttype

 
         WRITE (*, 99084) 
99084    FORMAT (' Do you want a binary? (y/n)')
         READ (*, 99004) iok
        
         IF (iok.EQ.'y' .AND. numpt.EQ.2) THEN
            WRITE (*, 99201) 
99201       FORMAT (' Enter binary mass')
            READ (*,*) bmass
            WRITE (*, 99202) 
99202       FORMAT (' Enter mass ratio')
            READ (*,*) qratio
            WRITE (*, 99203) 
99203       FORMAT (' Enter semi-major axis')
            READ (*,*) semiaxis
            WRITE (*, 99204) 
99204       FORMAT (' Enter eccentricity')
            READ (*,*) eccent
            qratio1 = 1.0 + qratio
            
            npart = npart + 2
            iphase(1) = ipttype
            iphase(2) = ipttype
            listpm(1) = 1
            listpm(2) = 2
            ra = semiaxis*(1.0+eccent)
            xyzmh(4,1) = bmass/qratio1
            xyzmh(4,2) = qratio*xyzmh(4,1)
            xyzmh(1,1) = -xyzmh(4,2)*ra/bmass
            xyzmh(2,1) = 0.
            xyzmh(3,1) = 0.
            xyzmh(1,2) = xyzmh(4,1)*ra/bmass
            xyzmh(2,2) = 0.
            xyzmh(3,2) = 0.
            vel = SQRT(semiaxis*(1.0-eccent*eccent)*bmass)/ra
            vxyzu(1,1) = 0.
            vxyzu(2,1) = -vel*qratio/qratio1
            vxyzu(3,1) = 0.
            vxyzu(1,2) = 0.
            vxyzu(2,2) = vel/qratio1
            vxyzu(3,2) = 0.
            WRITE (*,*) bmass, qratio, semiaxis, eccent
            WRITE (*,*) xyzmh(1,1), xyzmh(1,2), vxyzu(2,1), vxyzu(2,2)
            WRITE (*,*) xyzmh(4,1), xyzmh(4,2)
            totptmass = bmass
            spinx(1) = 0.
            spiny(1) = 0.
            spinz(1) = 0.               
            angaddx(1) = 0.
            angaddy(1) = 0.
            angaddz(1) = 0.               
            spinx(2) = 0.
            spiny(2) = 0.
            spinz(2) = 0.               
            angaddx(2) = 0.
            angaddy(2) = 0.
            angaddz(2) = 0.               
         ELSE
            DO i = 1, numpt
               npart = npart + 1
               iphase(npart) = ipttype
               listpm(i) = npart
               WRITE (*, 99087) 
99087          FORMAT (' Enter positions (x,y,z)')
               READ (*,*) tx,ty,tz
               xyzmh(1,npart) = tx
               xyzmh(2,npart) = ty
               xyzmh(3,npart) = tz
               WRITE (*, 99088) 
99088          FORMAT (' Enter velocities (vx,vy,vz)')
               READ (*,*) tvx,tvy,tvz
               vxyzu(1,npart) = tvx
               vxyzu(2,npart) = tvy
               vxyzu(3,npart) = tvz
               WRITE (*, 99089) 
99089          FORMAT (' Enter mass')
               READ (*,*) tmass
               xyzmh(4,npart) = tmass
               totptmass = totptmass + tmass
               spinx(npart) = 0.
               spiny(npart) = 0.
               spinz(npart) = 0.               
               angaddx(npart) = 0.
               angaddy(npart) = 0.
               angaddz(npart) = 0.
            END DO
         ENDIF

         IF (iok.EQ.'y' .AND. numpt.EQ.2) THEN
            WRITE (*, 99093) 
99093       FORMAT (' Do you want Roche-lobe accretion radii? (y/n)')
            READ (*, 99004) iok
            IF (iok.EQ.'y') THEN
               WRITE (*, 99094) 
99094          FORMAT (' Do you want variable accretion radii? (y/n)')
               READ (*, 99004) iok2
               IF (iok2.EQ.'y') THEN
                  iaccevol = 'v'
               ENDIF
               WRITE (*, 99095) 
99095          FORMAT (' Enter fraction of Roche-lobe size')
               READ (*,*) accfac
c
c--Sets accretion radii to be the roche lobe sizes 
c     (see Accretion Power in Astrophysics, Frank, King, & Raine)
c
               IF (qratio.GE.0.05) THEN
                  xyzmh(5,1) = accfac*semiaxis*(0.38-0.20*LOG10(qratio))
               ELSE
                  WRITE (*,*) 'qratio < 0.05'
                  CALL quit
               ENDIF
               IF (qratio.LT.0.5) THEN
                  xyzmh(5,2) = accfac*semiaxis*
     &                 (0.462*(qratio/qratio1)**(1.0/3.0))
               ELSE
                  xyzmh(5,2) = accfac*semiaxis*(0.38+0.20*LOG10(qratio))
               ENDIF
               WRITE (*,*) 'Accretion radii = ', xyzmh(5,1), xyzmh(5,2)
               hacc = xyzmh(5,2)
               haccall = xyzmh(5,2)
            ELSE
 555           WRITE (*, 99090) 
99090          FORMAT (' Enter Outer Accretion radius')
               READ (*,*) hacc
               WRITE (*, 99092) 
99092          FORMAT (' Enter Inner Accretion radius')
               READ (*,*) haccall
               IF (haccall.GT.hacc) THEN
                  WRITE(*,*) 'Inner <= Outer Radius'
                  GOTO 555
               ENDIF
               DO i = 1, numpt
                  xyzmh(5,i) = hacc
               END DO
            ENDIF
         ELSE
 556        WRITE (*, 99090) 
            READ (*,*) hacc
            WRITE (*, 99092) 
            READ (*,*) haccall
            IF (haccall.GT.hacc) THEN
               WRITE(*,*) 'Inner <= Outer Radius'
               GOTO 556
            ENDIF
            DO i = 1, numpt
               xyzmh(5,i) = hacc
            END DO
         ENDIF
         specang = SQRT(totptmass)
         ptmassin = totptmass
         nptmass = numpt
         iptmass = 0
      ELSE
         nptmass = 0
         iptmass = 0
      ENDIF

      DO i = 1, nptmass
         rho(i) = 0.
         dgrav(i) = 0.
      END DO
      
      WRITE (*,*) 'Number of point masses created: ',nptmass

      WRITE (*, 99061)
99061 FORMAT (' Do you want dynamic point mass creation (y/n)')
      READ (*, 99004) iok
      IF (iok.EQ.'y') THEN
         WRITE (*, 99091) 
         READ (*,*) ipttype
         iptmass = ipttype
 557     WRITE (*, 99090) 
         READ (*,*) hacc
         WRITE (*, 99092) 
         READ (*,*) haccall
         IF (haccall.GT.hacc) THEN
            WRITE(*,*) 'Inner <= Outer Radius'
            GOTO 557
         ENDIF
         WRITE (*, 99063)
99063    FORMAT (' Enter critical radius for point mass creation',/
     &           '    (rad. from another pt. mass) in code units')
         READ (*, *) radcrit
         WRITE (*, 99064)
99064    FORMAT (' Enter critical density for point mass creation',/
     &           '    in units of the initial density rhozero')
         READ (*, *) ptmcrit
      ENDIF
c
c--Specify Geometry
c
  50  WRITE (*, 88001)
88001 FORMAT (//, ' INITIAL CLOUD GEOMETRY ', //,
     &        '                       cube : 1 ',/,
     &        '                   cylinder : 2 ',/,
     &        '                     sphere : 3 ',/,
     &        '             accretion disk : 4 ',/,
     &        '                  ellipsoid : 5 ',/,
     &        '  constant angular momentum : 6 ',/,
     &        '           sphere with hole : 7 ')
      READ (*, *) igeom
      IF (igeom.LT.1 .OR. igeom.GT.7) GOTO 50
c
c--Specify Box Size
c
      rmax = 0.
      rcyl = 0.
      rmin = 0.
      rmind = 0.
      ichang = 1

 100  IF (igeom.EQ.1)  THEN 
         WRITE (*, 99002)udist
99002    FORMAT (//, ' PARTICLE SET UP ', //,
     &   ' Bounds of particle distribution in units of ',1pe14.5,
     &   ' (cm):  xmin, xmax, ymin, ymax, zmin, zmax')
         READ (*, *) xmin, xmax, ymin, ymax, zmin, zmax
         rmax = SQRT((xmax/2)**2 + (ymax/2)**2 +
     &               (zmax/2)**2)
         rcyl = SQRT((xmax/2)**2 + (ymax/2)**2) 

      ELSE IF (igeom.EQ.2) THEN
         WRITE (*, 88003) udist
88003    FORMAT (//, ' PARTICLE SET UP ', //,
     &   ' Bounds of particle distribution in units of ',1pe14.5,
     &   '(cm): rcyl, L/D')
         READ (*, *) rcyl, faclod
         zmax = rcyl * faclod
         zmin = - zmax
         xmax = rcyl
         ymax = rcyl
         xmin = - xmax
         ymin = - ymax
         rmax = SQRT(rcyl*rcyl + zmax*zmax)

      ELSE IF (igeom.EQ.3 .OR. igeom.EQ.7) THEN
         WRITE (*, 88002) udist
88002    FORMAT (//, ' PARTICLE SET UP ', //,
     &   ' Bounds of particle distribution in units of ',1pe14.5,
     &   '(cm):  rmax')
         READ (*, *) rmax

         IF (igeom.EQ.7) THEN
         WRITE (*, 88102)
88102    FORMAT (' Inner hole radius:  rmind')
         READ (*, *) rmind
         ENDIF

         xmax = rmax
         zmax = rmax
         ymax = rmax
         xmin = - xmax
         ymin = - ymax
         zmin = - zmax
         rcyl = rmax

      ELSE IF (igeom.EQ.4) THEN
         WRITE (*, 89003) udist
89003    FORMAT (//, ' PARTICLE SET UP ', //,
     &   ' Maximum and minimum radius of disk',/,' in units of ',
     &   1PE14.5,'(cm): rmaxd, rmind,',/,
     &   ' and the disk height to radius ratio, h/r')
         READ (*, *) rcyl, rmind, faclod
         zmax = rcyl * faclod
         zmin = - zmax
         xmax = rcyl
         ymax = rcyl
         xmin = - xmax
         ymin = - ymax
         rmax = SQRT(rcyl*rcyl + zmax*zmax)

      ELSE IF (igeom.EQ.5) THEN
         WRITE (*, 89004) udist
89004    FORMAT (//, ' PARTICLE SET UP ', //,
     &   ' Maximum and minimum radius of ellipsoid',/,
     &   ' in units of ',1pe14.5,'(cm): rmaxd, rmind,',/,
     &   ' and the ellipsoid height to radius ratio, h/r')
         READ (*, *) rcyl, rmind, faclod
         zmax = rcyl * faclod
         zmin = - zmax
         xmax = rcyl
         ymax = rcyl
         xmin = - xmax
         ymin = - ymax
         rmax = SQRT(rcyl*rcyl + zmax*zmax)

      ELSE IF (igeom.EQ.6) THEN
         WRITE (*, 89005) udist
89005    FORMAT (//, ' PARTICLE SET UP ', //,
     &   ' Maximum radius of cloud',/,
     &   ' in units of ',1pe14.5,'(cm): rmax,')
         READ (*, *) rmax
         angvel = SQRT(gg*totptmass*umass/udist**3.0)
         WRITE (*, 89006) angvel
89006    FORMAT(/, ' Angular velocity in rad/s at radius=1 is ',
     &        1PE14.6, ' for circular orbit')
         WRITE (*, 89007) 
89007    FORMAT(/, ' Enter angular velocity at radius=1 in fraction ',
     &        'of above units')
         READ (*, *) angvel
         angvel = angvel*SQRT(totptmass)
         WRITE (*, 89008) 
89008    FORMAT(/, ' Enter radial velocity in fraction of free-fall 
     &        units')
         READ (*, *) fracradial

         rmind = hacc
         xmax = rmax
         ymax = rmax
         zmax = rmax
         xmin = - xmax
         ymin = - ymax
         zmin = - zmax
         rcyl = rmax
      ENDIF

      deltax = xmax - xmin
      deltay = ymax - ymin
      deltaz = zmax - zmin
      rmax2 = rmax * rmax
      rmind2 = rmind * rmind
      rcyl2 = rcyl * rcyl
      zmax2 = zmax * zmax

 150  WRITE (*, 88008)
88008 FORMAT (//,' Do you want boundaries?',/,
     &  '                               0 : no boundaries',/,
     &  '   Reflective constant volume boundries with ghosts:',/,
     &  '                               1 : cartesian boundaries',/,
     &  '                               2 : cylindrical boundaries',/,
     &  '                               3 : spherical boundaries',/,
     &  '   Constant pressure boundaries - subtract pressure:',/,
     &  '                               7 : cons. pres. boundaries',/,
     &  '   Cartesian periodic boundaries with ghosts:',/,
     &  '                              11 : periodic boundaries',/,
     &  '   Dead particle boundaries:',/,
     &  '                               8 : constant N particles',/,
     &  '                              90 : constant N infall,',
     &     ' constant angular momentum',/,
     &  '                              91 : constant N infall,',
     &     ' constant omega')
      READ (*, *) ibound
      IF (ibound.LT.0 .OR. (ibound.GT.11 .AND. ibound.NE.90 .AND. 
     &     ibound.NE.91)) GOTO 150
      IF (ibound.EQ.7) THEN
         WRITE(*,88111)
88111    FORMAT(/,'   What is the external temperature (units K)?')
         READ (*,*) exttemp
         extu = 3.0/2.0*exttemp*Rg/gmw/uergg
         WRITE(*,88112)
88112    FORMAT(/,'   External density (gm/cc) (-ve=use rhozero)?')
         READ (*,*) extdens
         IF (extdens.GT.0.) THEN
            extdens = extdens/udens
            pext = 2.0/3.0*extu*extdens
         ENDIF
      ENDIF
      IF (ibound.EQ.8) THEN
         WRITE(*,88113)
88113    FORMAT(/,'   What is the dead particle radius?')
         READ (*,*) deadbound
         WRITE(*,88114)
88114    FORMAT(/,'   How many particles to accrete before stopping?')
         READ (*,*) nstop
         WRITE(*,88115)
88115    FORMAT(/,'   How many particles from end start fast dumping?')
         READ (*,*) nfastd
      ENDIF
      IF (ibound/10.EQ.9) THEN
         WRITE(*,88113)
         READ (*,*) deadbound
         WRITE(*,88116)
88116    FORMAT(/,'   How many particles in shell near outer boundary?')
         READ (*,*) nshell
         WRITE(*,88117)
88117    FORMAT(/,'   What is the minimum radius of the shell?')
         READ (*,*) rshell
         rmind = rshell
c         WRITE(*,88118)
c88118    FORMAT(/,'   What is the minimum radius of the STARTING shell?')
c         READ (*,*) rmind
      ENDIF
      nreassign = 0
      nkill = 0
      naccrete = 0
      anglostx = 0.
      anglosty = 0.
      anglostz = 0.
c
c--Determine Distribution
c
 200  IF (ichang.EQ.1 .OR. ichang.EQ.3) THEN
         WRITE (*, 99006)
99006    FORMAT (/,' What distribution :  cubic lattice = 1', /,
     &           '                      close packed  = 2', /,
     &           '                      random        = 3', /,
     &           '                      body centred  = 4', /,
     &           '                      custom        = 5', /,
     &           '                      rings         = 6')
         READ (*, *) idist
      ENDIF
 210  WRITE (*, 88009)
88009 FORMAT(' Do you want to enter the spacing between particles',
     &       ' or the number of particles? (s/n)')
      READ (*, 99004) iok

      IF (iok.eq.'n' .OR. iok.eq.'N') THEN
         WRITE (*, 88004)
88004    FORMAT(//,' Enter number of particles to use in ',
     &          ' the simulation ')
         READ (*, *) np
         IF (idist.EQ.4) THEN
            h1 = (deltax * deltay * deltaz * 2.0 / np) ** (1./3.)
            delmin = MIN(deltax,deltay)
            delmin = MIN(delmin,deltaz)
            h1 = delmin/FLOAT(INT(delmin/h1))
         ELSE
            h1 = (deltax*deltay*deltaz/np) ** (1./3.)
         ENDIF
         WRITE (*,*) 'Particle spacing, h1 = ', h1
         GOTO 320

      ELSE IF (iok.NE.'s' .AND. iok.NE.'S') THEN 
         GOTO 210
      ENDIF
c
c--Give spacing between particles
c
 300  IF (ichang.EQ.1 .OR. ichang.EQ.4) THEN
         WRITE (*, 99007) udist
99007    FORMAT (' Spacing between particles in ',1PE14.5,' (cm)?')
         READ (*, *) h1
      ENDIF
c
c--Compute number of particles and check if ok
c
 320  IF (idist.EQ.1 .OR. idist.EQ.3 .OR. idist.EQ.5 .OR. 
     &                                           idist.EQ.6) THEN
         facx = 1.
         facy = 1.
         facz = 1.
         dely = 0.
         delx = 0.
      ELSEIF (idist.EQ.4) THEN
         facx = 1.
         facy = 1.
         facz = 0.5
         dely = 0.
         delx = 0.
      ELSE
         facx = 1.
         facy = SQRT(3./4.)
         facz = SQRT(6.)/3.
         delx = 0.5*h1
         dely = h1*SQRT(3.)/6.
      ENDIF

      IF (idist.EQ.1 .OR. idist.EQ.2) THEN
         nx = INT(deltax/(facx*h1)) + 1
         ny = INT(deltay/(facy*h1)) + 2
         nz = INT(deltaz/(facz*h1)) + 1
         np = nx*ny*nz
      ELSEIF (idist.EQ.4 .OR. idist.EQ.5) THEN
         nx = INT(deltax/(facx*h1)) + 1
         ny = INT(deltay/(facy*h1)) + 1 
         nz = INT(deltaz/(facz*h1)) + 2
         np = nx*ny*nz
      ELSEIF (idist.EQ.6) THEN
         nx = INT(deltax/(2.0*facx*h1)) + 1
         ny = INT(pi*deltax/(facx*h1)) + 1
         nz = INT(deltaz/(facx*h1)) + 1
         np = nx*ny*nz
      ENDIF

  350 IF (np.LE.idim) THEN
         WRITE (*, 99008) np
99008    FORMAT (' Total number of particles needed :', I8, /,
     &           ' is that ok? (y/n)')
      ELSE
         WRITE (*, 99009) np, idim
99009    FORMAT (' Total number of particles needed :', I8, /,
     &           ' this number EXCEEDS the dimensions set to ', I8, /,
     &           ' is this still ok? (y/n)')
      ENDIF
      READ (*, 99004) iok
      IF (iok.EQ.'y') GOTO 450
 400  WRITE (*, 99010)
99010 FORMAT (' What do you want to change : all      = 1', /,
     &        '                              box size = 2', /,
     &        '                              distrib. = 3', /,
     &        '                              spacing  = 4', /,
     &        '                              npart    = 5')
      READ (*, *) ichang
      IF (ichang.LE.2) GOTO 100
      IF (ichang.EQ.3) GOTO 200
      IF (ichang.EQ.4) GOTO 300
      IF (ichang.EQ.5) GOTO 210
c
c--Set Density Distribution
c
 450  WRITE(*, 99014)
99014 FORMAT (' What density variations:', /,
     &   '  Uniform particle distribution, Non-uniform masses  (m)', /,
     &   '  Uniform particle masses, Non-uniform distribution  (d)', /,
     &   '  Uniform density                                    (u)')
      READ (*,99004) idens
      IF ((idens.NE.'m').AND.(idens.NE.'d').AND.(idens.NE.'u')) GOTO 450

      IF (idens.NE.'u') THEN
c
c--Set Non-Uniform Density Distribution
c
 550     WRITE(*, 99016)
99016    FORMAT (' Coords for density variations: Cartesian    (1)', /,
     &      '                                 Cylindrical (2)', /,
     &      '                                 Spherical   (3)')
         READ (*,*) icoord
         IF ((icoord.LT.1).OR.(icoord.GT.3)) GOTO 550

         IF (idens.EQ.'m') THEN
            CALL unifdis(igeom, idist, np, h1, facx, facy, facz, 
     &                      delx, dely, nx, ny, nz, ibound)
            IF (notdone.EQ.1) GOTO 400

            WRITE(*, 99105)
            READ (*, 99004) icentral
            IF (icentral.EQ.'y') CALL condense

            IF (icoord.EQ.1) THEN
               CALL cartmas
            ELSEIF (icoord.EQ.2) THEN
               CALL cylmas
            ELSE
               CALL sphmas
            ENDIF
         ELSE
            IF (icoord.EQ.1) THEN
               CALL cartdis(igeom, idist, np, h1)
            ELSE IF (icoord.EQ.2) THEN
               CALL cyldis(igeom,np)
            ELSE
               CALL sphdis(igeom, idist, np ,h1, facx, facy, facz, 
     &                   delx, dely, nx, ny, nz)
            ENDIF
         ENDIF
      ELSE
c
c--Set Uniform Density Distribution
c
         CALL unifdis(igeom, idist, np, h1, facx, facy, facz, 
     &                   delx, dely, nx, ny, nz, ibound)
         IF (notdone.EQ.1) GOTO 400
 575     WRITE(*, 99105)
99105    FORMAT (' Do you want to centrally condense the particles?')
         READ (*, 99004) icentral

         IF (icentral.EQ.'y') CALL condense
      ENDIF
c
c--Rotate particle distribution about y-axis
c
 500  WRITE (*, 99110)
99110 FORMAT (' Do you want to rotate particles about y-axis?')
      READ (*, *) irotatey
      IF (irotatey.EQ.'Y' .OR. irotatey.EQ.'y') THEN
         WRITE (*, 99115)
99115    FORMAT (' Enter angle (deg)?')
         READ (*,*) rangle
         rangle = rangle*pi/180.0
         DO i = 1, npart
            r=SQRT(xyzmh(1,i)**2 + xyzmh(3,i)**2)
            th=ATAN2(xyzmh(3,i), xyzmh(1,i)) + rangle
            xyzmh(1,i)=r*COS(th)
            xyzmh(3,i)=r*SIN(th)
         END DO
      ENDIF
c
c--Set Total Mass
c
      umassr = umass/solarm
      WRITE (*, 99011)
99011 FORMAT (' Do you want to enter: Total mass of the system (1)', /,
     &        '                       Particle mass            (2)')
      READ (*, *) ichmass
      IF (ichmass.EQ.1) THEN
         WRITE (*, 99012) umassr
99012    FORMAT (' Enter total mass of system in units of ',1pe14.5,
     &        ' solar masses')
         READ (*, *) totmas
         partm = totmas/(npart - nptmass)
      ELSEIF (ichmass.EQ.2) THEN
         WRITE (*, 99013) umassr
99013    FORMAT (' Enter particle masses in units of ',1pe14.5,
     &        ' solar masses')
         READ (*, *) partm
         totmas = partm*FLOAT(npart - nptmass)
      ELSE
         GOTO 500
      ENDIF

      IF (igeom.EQ.1) THEN
         rhozero = totmas / (deltax * deltay * deltaz)
      ELSEIF (igeom.EQ.2) THEN
         rhozero = totmas / (pi * deltaz * rcyl2)
      ELSEIF (igeom.EQ.3) THEN
         rhozero = totmas / (4.0 * pi * rmax2 * rmax / 3.0) 
      ELSEIF (igeom.EQ.4) THEN
         rhozero = totmas / (pi * deltaz * (rcyl2 - rmind2))
      ELSEIF (igeom.EQ.5) THEN
         rhozero = totmas / (4.0*pi/3.0*(rcyl2*zmax - rmind2*rmind))
         IF (rhozero.LE.0.0) THEN
            WRITE(*,*) 'Negative volume'
            CALL quit
         ENDIF
      ELSEIF (igeom.EQ.6) THEN
         rhozero = totmas / (4.0*pi*(rmax2*rmax - hacc*hacc*hacc)/3.0)
      ELSEIF (igeom.EQ.7) THEN
         rhozero = totmas / (4.0*pi/3.0 * (rmax2*rmax - rmind2*rmind)) 
      ELSE
         WRITE(*,*) 'ERROR igeom'
         CALL quit
      ENDIF
      WRITE(*,*) ' rhozero = ',rhozero
      IF (ibound.EQ.7) THEN
         IF (extdens.LE.0.) pext = 2./3.*extu*rhozero
         WRITE(*,*) ' setting external pressure = ',pext
      ENDIF
c
c--Set Particle Masses
c
      DO i = nptmass + 1, npart
         iphase(i) = 0
         rho(i) = 0.
c         IF (idens.EQ.'m') THEN
            xyzmh(4,i) = partm * disfrac(i)
c         ELSE   
c            xyzmh(4,i) = partm
c         ENDIF
         dgrav(i) = 0.
      END DO
c
c--Check if distribution is ok
c
      WRITE (*, 99020) npart - nptmass
99020 FORMAT (1X, I8, ' particles have been set, is this ok? (y/n)')
      READ (*, 99004) iok
      IF (iok.EQ.'n') THEN
         WRITE (*, 99010)
         READ (*, *) ichang
         IF (ichang.LE.2) GOTO 100
         IF (ichang.EQ.3) GOTO 200
         IF (ichang.EQ.4) GOTO 300
         IF (ichang.EQ.5) GOTO 210
      ENDIF
c
c--Magnetic field setup once we know rhozero
c
      IF (imhd.EQ.idim) call setBfield
c
c--Set e.o.s. related quantities
c
      WRITE (*, 99022)
99022 FORMAT (' What is the equation of state variable:', /,
     &        ' specific internal energy :  intener (i)', /,
     &        ' specific entropy         :  entropy (e)')
      READ (*, 99024) varsta
99024 FORMAT (A7)
      IF (varsta.EQ.'entropy' .OR. varsta.EQ.'e') THEN
         WRITE (*, 99026)
         varsta = 'entropy'
99026    FORMAT (' Value of specific entropy (code units)')
         READ (*, *) thermal
      ELSE
         varsta = 'intener'
99027    WRITE (*, 99028)
99028    FORMAT (' Do you want to enter: ',/,
     &     '    value of initial average temperature in kelvin: (t) ',/,
     &     '    initial average sound speed (code units)      : (s) ',/,
     &     '    initial average internal energy (code units)  : (u) ?')
         READ(*,99128) ien
99128    FORMAT(A1)	 
	 IF (ien.EQ.'t') THEN
	    WRITE(*,99129)
99129       FORMAT (' Value of initial average temperature in kelvin:')
	    READ (*,*) thermal
            thermal = 3.0/2.0*thermal*Rg/gmw/uergg
	 ELSE IF (ien.EQ.'u') THEN
	    WRITE(*,99130)
99130       FORMAT (' Value of internal energy (code units)') 
	    READ (*,*) thermal
	 ELSE IF (ien.EQ.'s') THEN
	    WRITE(*,99131) udist/utime
99131       FORMAT (' Value of sound speed in units of ',1pe10.4,'cm/s')
	    READ (*,*) thermal
	    vsound = ABS(thermal)
	    vsound2 = vsound*vsound
	 ELSE
            GOTO 99027
         ENDIF
      ENDIF
      rhocrt = rhocrit * udens
      rhocrt2 = rhocrit2 * udens
      rhocrt3 = rhocrit3 * udens

      IF (varsta.EQ.'intener') THEN
 616     WRITE (*, 99030) rhocrt, rhocrt2, rhocrt3
99030    FORMAT (' Equation of state/Energy calculation:',/,
     &        '   Isothermal,     p=2/3*u*rho   (i)',/,
     &        '   Adiabatic,      p=2/3*u*rho   (a)',/,
     &        '   Polytropic,     p=A*rho^gamma (p)',/,
     &        '   Variable gamma, p=A*rho^gamma (v)',/,
     &        '      critical rho (s) = ', 1PE14.5, 1PE14.5, 1PE14.5)
         READ (*, 99004) encal
         IF (encal.EQ.'p') THEN
            WRITE (*,99032)
99032       FORMAT (' Enter gamma')
            READ (*,*) gamma
            gm1 = gamma - 1.0
            WRITE (*,99033)
99033       FORMAT (' Do you want to specify the polytropic K? (y/n)',/,
     &           '(otherwise uses initial temp/internal energy/vsound)')
            READ (*,*) iok
            IF (iok.eq.'y'.or.iok.EQ.'Y') THEN
               WRITE(*,99133)
99133          FORMAT (' Enter polytropic K:')         
               READ (*,*) RK2
               RK2 = RK2*1.5
               WRITE(*,*) 'RK2 = ',RK2
            ELSE            
               RK2 = thermal/(rhozero**gm1)
            ENDIF            
         ELSE IF (encal.EQ.'a') THEN
            gamma = 5.0/3.0
            gm1 = gamma - 1.0
	    IF (ien.EQ.'s') thermal = vsound2/(gamma*gm1)
            RK2 = thermal/(rhozero**gm1)
         ELSE IF (encal.EQ.'i' .OR. encal.EQ.'t') THEN
            gamma = 1.0
            RK2 = thermal
	    IF (ien.EQ.'s') thermal = 1.5*vsound2
            RK2 = thermal
            tempiso = 2./3.*thermal/(Rg/gmw/uergg)
            WRITE(*,*) 'isothermal temperature = ',tempiso
         ELSE IF (encal.EQ.'v') THEN
c
c--Value of gamma is irrelevant for definition of variable e.o.s.
c
            gamma = 5.0/3.0
            gm1 = gamma - 1.0
	    IF (ien.EQ.'s') thermal = vsound2*1.5  !!!/(gamma*gm1)
            RK2 = thermal/(rhozero**gm1)
         ELSE IF (encal.EQ.'x') THEN
c
c--Value of gamma is irrelevant for definition of physical e.o.s.
c
            gamma = 5.0/3.0
            gm1 = gamma - 1.0
	    IF (ien.EQ.'s') thermal = vsound2/(gamma*gm1)
            RK2 = thermal/(rhozero**gm1)
         ELSE
            GOTO 616
         ENDIF
      ELSE         
         WRITE (*,99034)
99034    FORMAT (' Enter polytropic gamma')
         READ (*,*) gamma
      ENDIF
c
c--Set Pressure Distribution
c
 617  WRITE(*, 99036)
99036 FORMAT (' Do you want pressure variations (y/n)')
      READ (*,99004) ipres
      IF ((ipres.NE.'y').AND.(ipres.NE.'n')) GOTO 617

      IF (ipres.EQ.'y') THEN
c
c--Set Non-Uniform Pressure Distribution
c
 618     WRITE(*, 99038)
99038    FORMAT (' Coords for pressure variations: Cartesian (1)', /,
     &        '                                 Cylinderical (2)', /,
     &        '                                 Spherical    (3)')
         READ (*,*) icoord
         IF ((icoord.LT.1) .OR. (icoord.GT.3)) GOTO 618

         IF (icoord.EQ.1) THEN
            CALL cartpres
         ELSE IF (icoord.EQ.2) THEN
            CALL cylpres
         ELSE
            CALL sphpres
         ENDIF
         DO i = 1, npart
            vxyzu(4,i) = thermal*disfrac(i)
         END DO
      ELSE
c
c--Set Uniform Pressure Distribution
c
         DO i = 1, npart
            vxyzu(4,i) = thermal
         END DO
      ENDIF
c
c--Set in Perturbations in Position
c
      WRITE (*, 99040)
99040 FORMAT (' Random perturbations in particle positions? (y/n)')
      READ (*, 99004) iok
      IF (iok.EQ.'y') THEN
         WRITE (*, 99042)
99042    FORMAT (' Amplitude (in unit of spacing)')
         READ (*, *) dmax
         dmax05 = 0.5*dmax
         DO i = nptmass + 1, npart
 625        xi = xyzmh(1,i) + 2.0*h1*(ran1(1)*dmax - dmax05)
            IF (xi.LT.xmin .OR. xi.GT.xmax) GOTO 625
            xyzmh(1,i) = xi
 640        yi = xyzmh(2,i) + 2.0*h1*(ran1(1)*dmax - dmax05)
            IF (yi.LT.ymin .OR. yi.GT.ymax) GOTO 640
            xyzmh(2,i) = yi
 660        zi = xyzmh(3,i) + 2.0*h1*(ran1(1)*dmax - dmax05)
            IF (zi.LT.zmin .OR. zi.GT.zmax) GOTO 660
            xyzmh(3,i) = zi
         END DO
      ENDIF

c
c--Set Center of Mass at Zero
c
      cmx = 0. 
      cmy = 0.
      cmz = 0.
      DO i = 1,npart
         cmx = cmx + xyzmh(4,i)*xyzmh(1,i)
         cmy = cmy + xyzmh(4,i)*xyzmh(2,i)
         cmz = cmz + xyzmh(4,i)*xyzmh(3,i)
      END DO
      cmx = cmx/totmas
      cmy = cmy/totmas
      cmz = cmz/totmas
      WRITE (*,*) 'Centre of mass is at: ',cmx,cmy,cmz

      WRITE (*, 99045)
99045 FORMAT(' Do you want to set the centre of mass',/,
     & ' at zero ? (y/n) ')
      READ (*, 99004) iok
      IF (iok.EQ.'y' .OR. iok.EQ.'Y') THEN
         DO i = 1, npart
            xyzmh(1,i) = xyzmh(1,i) - cmx
            xyzmh(2,i) = xyzmh(2,i) - cmy
            xyzmh(3,i) = xyzmh(3,i) - cmz
         END DO
      END IF
c 
c--Find extrema of lattice
c
      xlmin = 1.E30
      xlmax = 0.
      ylmin = 1.E30
      ylmax = 0.
      zlmin = 1.E30
      zlmax = 0.
      DO i = 1, npart
         xlmax = MAX(xyzmh(1,i), xlmax)
         xlmin = MIN(xyzmh(1,i), xlmin)
         ylmax = MAX(xyzmh(2,i), ylmax)
         ylmin = MIN(xyzmh(2,i), ylmin)
         zlmax = MAX(xyzmh(3,i), zlmax)
         zlmin = MIN(xyzmh(3,i), zlmin)
      END DO
c
c--Self-Gravity?
c
      WRITE (*, 99046)
99046 FORMAT (' Is self-gravity included? (y/n)')
      READ (*, 99004) iok
      igphi = 0
      IF (iok.EQ.'y') igphi = 1
c
c--Get input file options
c
      CALL inopts
c
c--Set Velocities
c
      WRITE (*, 99047)
99047 FORMAT (' Do you want          rotation : (r)',/,
     &        '    or  a gaussian random P(k) : (p)',/,
     &        '    or a velocity distribution : (v)',/,
     &        '    or           no velocities : (n)')
      READ (*, 99004) iok
      IF (iok.EQ.'n') GOTO 777
      IF (iok.EQ.'r') GOTO 776
      IF (iok.EQ.'p') GOTO 779
c
c--Velocity Distribution
c
 765  WRITE(*, 99049)
99049 FORMAT (' Coords for velocity variations: Cartesian    (1)', /,
     &        '                                 Cylindrical  (2)', /,
     &        '                                 Spherical    (3)')
      READ (*,*) icoord
      IF ((icoord.LT.1) .OR. (icoord.GT.3)) GOTO 765

      IF (icoord.EQ.1) THEN
         CALL cartvel
      ELSE IF (icoord.EQ.2) THEN
         CALL cylvel
      ELSE
         CALL sphvel
      ENDIF

      GOTO 778
c
c--Set Rotation
c
 776  IF (igeom.EQ.6) THEN
         iok = 'd'
      ELSE
         WRITE (*, 99048)
         iok = 'n'
99048    FORMAT (' Do you want keplerian rotation :  (k)',/,
     &        '        or  solid body rotation :  (s)',/,
     &        '       or differential rotation :  (d)',/,
     &        '      or rotation perpendicular :  (p)',/,
     &        '        or no internal rotation :  (n)')
         READ (*, 99004) iok
      ENDIF

 777  IF (iok.EQ.'s' .OR. iok.EQ.'d' .OR. iok.EQ.'p') THEN
         IF (igeom.NE.6) THEN
            WRITE (*,99050)
99050       FORMAT (' Enter angular velocity in rad/sec at radius=1')
            READ (*,*) angvel
            angvel = angvel * utime
            IF (iok.EQ.'p' .OR. iok.EQ.'d') THEN
                 WRITE (*,99054) 
99054            FORMAT('     What rotation profile (omega)?',/,
     &                  '        omega ~ exponential (e) ?',/,
     &                  '        omega ~    1 over r (r) ?',/,
     &                  '        omega ~  1 over r^2 (s) ?')
                 READ (*, 99004) iwhat
            ENDIF
         ELSEIF (igeom.EQ.6) THEN
            fractan = angvel/specang
            iwhat = 's'
         ENDIF
      ENDIF

      IF (iexf.EQ.5 .OR. iexf.EQ.6 .OR. iok.EQ.'k') THEN
         WRITE (*,55504)
55504    FORMAT (' Enter mass for external forces')
         READ (*,*) xmass
      ELSE
         xmass = 0.
      ENDIF 

      IF (iok.EQ.'s') THEN
         WRITE (*,55506)
55506    FORMAT (' Enter fractional rotation gradient along z per',
     &        ' unit distance')
         READ (*,*) fracrotgrad
         WRITE (*,55507)
55507    FORMAT (' Enter offset for rotation gradient along z')
         READ (*,*) fracrotoffset
      ENDIF

      DO i = nptmass + 1, npart
         IF (iok.EQ.'k') THEN
            gg2 = 1.  
            radius = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
            xtild = xyzmh(1,i)/radius
            ytild = xyzmh(2,i)/radius
            alpha = SQRT(gg2 * xmass/radius)
            vxyzu(1,i) = -alpha * ytild 
            vxyzu(2,i) =  alpha * xtild
            vxyzu(3,i) = 0.
         ELSEIF (iok.EQ.'s') THEN
            vxyzu(1,i) = -angvel * xyzmh(2,i) * (fracrotoffset + 
     &           xyzmh(3,i)*fracrotgrad)
            vxyzu(2,i) =  angvel * xyzmh(1,i) * (fracrotoffset + 
     &           xyzmh(3,i)*fracrotgrad)
            vxyzu(3,i) = 0.
         ELSEIF (iok.EQ.'d' .AND. igeom.EQ.6) THEN 
            rxy2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
            rxy = SQRT(rxy2)
            r2 = rxy2 + xyzmh(3,i)**2
            r = SQRT(r2)
            IF (ibound.EQ.90) THEN
               fractanhere = fractan
            ELSEIF (ibound.EQ.91) THEN
               fractanhere = fractan*(rxy2/r2)
            ELSE
               WRITE (*,*) 'ERROR - ibound'
            ENDIF
            tanmag = fractanhere*specang/rxy
c            radmag = fracradial*SQRT(2.0*specang*specang/r - 
c     &           tanmag*tanmag)
c            vxyzu(1,i) = - radmag*xyzmh(1,i)/r - tanmag*xyzmh(2,i)/rxy
c            vxyzu(2,i) = - radmag*xyzmh(2,i)/r + tanmag*xyzmh(1,i)/rxy
c            vxyzu(3,i) = - radmag*xyzmh(3,i)/r
            vxyzu(1,i) = - tanmag*xyzmh(2,i)/rxy
            vxyzu(2,i) = tanmag*xyzmh(1,i)/rxy
            vxyzu(3,i) = 0.
         ELSEIF (iok.EQ.'d') THEN 
            rad2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
ccc            rad2 = rad2 / rcyl2
            IF (iwhat.eq.'e') THEN
               angvr = angvel / (exp(0.6667 * rad2))
            ELSE IF (iwhat.eq.'r') THEN
               angvr = angvel / (rad2 ** 0.5)
            ELSE IF (iwhat.eq.'s') THEN
               angvr = angvel / (rad2)
            ENDIF
            vxyzu(1,i) = -angvr * xyzmh(2,i) 
            vxyzu(2,i) =  angvr * xyzmh(1,i)
            vxyzu(3,i) = 0.
         ELSEIF (iok.EQ.'p') THEN 
            rad2 = xyzmh(3,i)**2 + xyzmh(2,i)**2
ccc            rad2 = rad2 / rmax2
            IF (iwhat.eq.'e') THEN
               angvz = angvel / (exp(0.6667 * rad2))
            ELSE IF (iwhat.eq.'r') THEN
               angvz = angvel / (rad2 ** 0.5)
            ELSE IF (iwhat.eq.'s') THEN
               angvz = angvel / (rad2)
            ENDIF
            vxyzu(3,i) = -angvz * xyzmh(2,i) 
            vxyzu(2,i) =  angvz * xyzmh(3,i)
            vxyzu(1,i) = 0.
         ELSE
            vxyzu(1,i) = 0.
            vxyzu(2,i) = 0.
            vxyzu(3,i) = 0.
         ENDIF
      END DO
      GOTO 778

 779  WRITE (*,*) 'Do you want random vels?'
      READ (*,99004) iok

      IF (iok.EQ.'Y' .OR. iok.EQ.'y') THEN
         WRITE (*,*) 'Enter amplitude'
         READ (*,*) amplitude

         ekinetic = 0.
         DO i = nptmass + 1, npart
66544       vxyzu(1,i) = ran1(1)
            vxyzu(2,i) = ran1(1)
            vxyzu(3,i) = ran1(1)
            vel2 = vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2
            vxyzu(1,i) = amplitude*vxyzu(1,i)
            vxyzu(2,i) = amplitude*vxyzu(2,i)
            vxyzu(3,i) = amplitude*vxyzu(3,i)
            IF (vel2.GT.1.0) GOTO 66544
            ekinetic = ekinetic + xyzmh(4,i)*amplitude**2*vel2
         END DO
         GOTO 66545
      ENDIF

      WRITE (*,*) 'Enter names of velocity files (x,y,z)'
      READ (*,77771) filevelx
      READ (*,77771) filevely
      READ (*,77771) filevelz
      
      WRITE (*,*) 'Enter size of velocity files (e.g.N=32^3)'
      READ (*,*) nspace
      IF (nspace.GT.nvelmax) THEN
         WRITE (*,*) 'Dimensions of velocity arrays too small'
         STOP
      ENDIF

      WRITE (*,*) 'Enter amplitude'
      READ (*,*) amplitude

      WRITE (*,*) 'Do you want fall off near edge?'
      READ (*,99004) iok

      OPEN (45,FILE=filevelx,FORM='unformatted')
      READ (45) (((velx(i,j,k), i=1,nspace),j=1,nspace), k=1,nspace)
      CLOSE (45)
      OPEN (45,FILE=filevely,FORM='unformatted')
      READ (45) (((vely(i,j,k), i=1,nspace),j=1,nspace), k=1,nspace)
      CLOSE (45)
      OPEN (45,FILE=filevelz,FORM='unformatted')
      READ (45) (((velz(i,j,k), i=1,nspace),j=1,nspace), k=1,nspace)
      CLOSE (45)

c      deli = 1.0/REAL(nspace/2)
      radnorm = rmax
      IF (igeom.EQ.1) radnorm = xmax
      deli = radnorm/REAL(nspace/2)
      DO i = nptmass + 1, npart
         iposx = INT(xyzmh(1,i)/radnorm*(nspace/2)+(nspace/2)+0.5)
         iposx = MIN(MAX(iposx, 1),nspace-1)
         iposy = INT(xyzmh(2,i)/radnorm*(nspace/2)+(nspace/2)+0.5)
         iposy = MIN(MAX(iposy, 1),nspace-1)
         iposz = INT(xyzmh(3,i)/radnorm*(nspace/2)+(nspace/2)+0.5)
         iposz = MIN(MAX(iposz, 1),nspace-1)
         delx = xyzmh(1,i) - (iposx-(nspace/2)-0.5)/REAL(nspace/2)*
     &        radnorm
         dely = xyzmh(2,i) - (iposy-(nspace/2)-0.5)/REAL(nspace/2)*
     &        radnorm
         delz = xyzmh(3,i) - (iposz-(nspace/2)-0.5)/REAL(nspace/2)*
     &        radnorm
         radius = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2)
c
c--Find interpolated x-velocity
c
         velx1 = velx(iposx,iposy,iposz) + delx/deli*
     &        (velx(iposx+1,iposy,iposz)-velx(iposx,iposy,iposz))
         velx2 = velx(iposx,iposy+1,iposz) + delx/deli*
     &        (velx(iposx+1,iposy+1,iposz)-velx(iposx,iposy+1,iposz))
         vely1 = velx1 + dely/deli*(velx2-velx1)

         velx1 = velx(iposx,iposy,iposz+1) + delx/deli*
     &        (velx(iposx+1,iposy,iposz+1)-velx(iposx,iposy,iposz+1))
         velx2 = velx(iposx,iposy+1,iposz+1) + delx/deli*
     &     (velx(iposx+1,iposy+1,iposz+1)-velx(iposx,iposy+1,iposz+1))
         vely2 = velx1 + dely/deli*(velx2-velx1)

         vxyzu(1,i) = amplitude*(vely1 + delz/deli*(vely2-vely1))
c
c--Find interpolated y-velocity
c
         velx1 = vely(iposx,iposy,iposz) + delx/deli*
     &        (vely(iposx+1,iposy,iposz)-vely(iposx,iposy,iposz))
         velx2 = vely(iposx,iposy+1,iposz) + delx/deli*
     &        (vely(iposx+1,iposy+1,iposz)-vely(iposx,iposy+1,iposz))
         vely1 = velx1 + dely/deli*(velx2-velx1)

         velx1 = vely(iposx,iposy,iposz+1) + delx/deli*
     &        (vely(iposx+1,iposy,iposz+1)-vely(iposx,iposy,iposz+1))
         velx2 = vely(iposx,iposy+1,iposz+1) + delx/deli*
     &     (vely(iposx+1,iposy+1,iposz+1)-vely(iposx,iposy+1,iposz+1))
         vely2 = velx1 + dely/deli*(velx2-velx1)

         vxyzu(2,i) = amplitude*(vely1 + delz/deli*(vely2-vely1))
c
c--Find interpolated z-velocity
c
         velx1 = velz(iposx,iposy,iposz) + delx/deli*
     &        (velz(iposx+1,iposy,iposz)-velz(iposx,iposy,iposz))
         velx2 = velz(iposx,iposy+1,iposz) + delx/deli*
     &        (velz(iposx+1,iposy+1,iposz)-velz(iposx,iposy+1,iposz))
         vely1 = velx1 + dely/deli*(velx2-velx1)

         velx1 = velz(iposx,iposy,iposz+1) + delx/deli*
     &        (velz(iposx+1,iposy,iposz+1)-velz(iposx,iposy,iposz+1))
         velx2 = velz(iposx,iposy+1,iposz+1) + delx/deli*
     &     (velz(iposx+1,iposy+1,iposz+1)-velz(iposx,iposy+1,iposz+1))
         vely2 = velx1 + dely/deli*(velx2-velx1)

         vxyzu(3,i) = amplitude*(vely1 + delz/deli*(vely2-vely1))

         IF (iok.EQ.'Y' .OR. iok.EQ.'y') THEN
            IF (radius/radnorm.GT.0.9) THEN
               factor = 0.0
            ELSEIF (radius/radnorm.LT.0.8) THEN
               factor = 1.0
            ELSE
               factor = (0.9-radius/radnorm)*10.0
            ENDIF
            vxyzu(1,i) = vxyzu(1,i)*factor
            vxyzu(2,i) = vxyzu(2,i)*factor
            vxyzu(3,i) = vxyzu(3,i)*factor
         ENDIF
      END DO
c
c--Set Velocity of Center of Mass to Zero
c
66545 cmx = 0. 
      cmy = 0.
      cmz = 0.
      DO i = 1,npart
         cmx = cmx + xyzmh(4,i)*vxyzu(1,i)
         cmy = cmy + xyzmh(4,i)*vxyzu(2,i)
         cmz = cmz + xyzmh(4,i)*vxyzu(3,i)
      END DO
      cmx = cmx/totmas
      cmy = cmy/totmas
      cmz = cmz/totmas
      WRITE (*,*) 'Vel of centre of mass is: ',cmx,cmy,cmz

      WRITE (*, 99886)
99886 FORMAT(' Do you want to set the centre of mass',/,
     & ' velocity to zero ? (y/n) ')
      READ (*, 99004) iok
      IF (iok.EQ.'y' .OR. iok.EQ.'Y') THEN
         DO i = 1, npart
            vxyzu(1,i) = vxyzu(1,i) - cmx
            vxyzu(2,i) = vxyzu(2,i) - cmy
            vxyzu(3,i) = vxyzu(3,i) - cmz
         END DO
      END IF
c
      ekinetic = 0.
      rootmeansquare = 0.
      DO i = nptmass + 1, npart
         vel2 = vxyzu(1,i)**2+vxyzu(2,i)**2+vxyzu(3,i)**2
         rootmeansquare = rootmeansquare + vel2
         ekinetic = ekinetic + xyzmh(4,i)*vel2
      END DO
      ekinetic = ekinetic/2.0
      WRITE (*,*) 'Root mean square velocity = ',
     &     SQRT(rootmeansquare/(npart-nptmass)),rootmeansquare
      WRITE (*,*) 'Sound speed = ',SQRT(2.0/3.0*gamma*thermal)
      WRITE (*,*) 'Kinetic energy = ',ekinetic
      WRITE (*,*) 'Is this okay?'
      READ (*,99004) iok
      IF (iok.EQ.'y' .OR. iok.EQ.'Y') GOTO 778

      WRITE (*,*) 'Enter multiplicative factor'
      READ (*,*) factor
      DO i = 1, npart
         vxyzu(1,i) = vxyzu(1,i)*factor
         vxyzu(2,i) = vxyzu(2,i)*factor
         vxyzu(3,i) = vxyzu(3,i)*factor
      END DO


77771 FORMAT(A20)
 778  n1 = npart
      n2 = 0
c
c--adjust smoothing lengths (also calculates initial density)
c  MUST be done (to get density) if evolving B/rho
c
      nactive = npart
      WRITE (*, 99044)
99044 FORMAT(' Do you want to adjust smoothing length',/,
     + ' to have similar number of neighbours ? (y/n) ')
      READ (*, 99004) iok
      IF (iok.EQ.'y' .OR. iok.EQ.'Y' .OR. imhd.EQ.idim) CALL hcalc

      IF (imhd.EQ.idim) THEN      
         Bextx = Bxzero
         Bexty = Byzero
         Bextz = Bzzero
         WRITE(*,99055) Bextx,Bexty,Bextz
99055    FORMAT(' Setting external B field = ',3(1PE12.4,2X)) 
      ELSE
         Bextx = 0.
         Bexty = 0.
         Bextz = 0.
      ENDIF
       
      WRITE (*, 99056)
99056 FORMAT (//, ' END OF SETUP')
c
c--Write options
c
      CALL wrinsph

      IF (idebug.EQ.'setpart') THEN
         WRITE (iprint, 99058) (xyzmh(1,i), i=1, npart)
         WRITE (iprint, 99058) (xyzmh(2,i), i=1, npart)
         WRITE (iprint, 99058) (xyzmh(3,i), i=1, npart)
         WRITE (iprint, 99058) (vxyzu(1,i), i=1, npart)
         WRITE (iprint, 99058) (vxyzu(2,i), i=1, npart)
         WRITE (iprint, 99058) (vxyzu(3,i), i=1, npart)
         WRITE (iprint, 99058) (vxyzu(4,i), i=1, npart)
         WRITE (iprint, 99058) (xyzmh(4,i), i=1, npart)
         WRITE (iprint, 99058) (rho(i), i=1, npart)
         IF (imhd.EQ.idim) THEN
            WRITE (iprint, 99058) (Bevolxyz(1,i), i=1, npart)
            WRITE (iprint, 99058) (Bevolxyz(2,i), i=1, npart)
            WRITE (iprint, 99058) (Bevolxyz(3,i), i=1, npart)
         ENDIF
99058    FORMAT (1X, 5(1PE12.5,1X))
      ENDIF

      RETURN
      END
