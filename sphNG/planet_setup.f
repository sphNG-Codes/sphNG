      SUBROUTINE newplanet

c************************************************************
c                                                           *
c  This subroutine handles the particle set up for          *
c  planets alone, to ease my lifes difficulties.            *
c                                                           *
c  Ben Ayliffe - 17/03/2011                                 *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/task'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/tming'

      CHARACTER*7 where

      DATA where/'newrun'/
c
c--General set up
c
      WRITE (*, 99001)
99001 FORMAT (' GENERAL SET UP', //)
      WRITE (*, 99002)
99002 FORMAT (' give name of the run (20 char.)')
      READ (*, 99003) namerun
      WRITE (namenextrun,99003) namerun
99003 FORMAT (A20)
c
c--Write label of run
c
      CALL labrun

      WRITE (*, 99004)
99004 FORMAT (' type of initialization?', /,
     &        ' start from scratch         : scratch (s)', /,
     &        ' using existing files       : exist (e)' )
      READ (*, 99005) what
99005 FORMAT (A7)
c
c--Start from scratch
c
      IF (what.EQ.'scratch' .OR. what.EQ.'s') THEN
         WRITE (*, 99006)
99006    FORMAT (' name of binary file (7 char. max)')
         READ (*, 99005) file1
         WRITE (*,*) 'MAXIMUM RECORD LENGTH = ',imaxrec,' ',file1
         OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         CALL planet_setup

         CALL inform(where)
         DO i = 1, npart
            isort(i) = i
            iorig(i) = i
         END DO
         nfullstep = 1
         CALL wdump_wrapper(idisk1)
      ENDIF
c
c--Use existing dumps to create new dump
c
      IF (what(1:5).EQ.'exist' .OR. what.EQ.'e') THEN
         WRITE (*, 99007)
99007    FORMAT (' name of the first binary file?')
         READ (*, 99005) file1
         WRITE (*, 99008)
99008    FORMAT (' name of the second binary file?')
         READ (*, 99005) file2
         WRITE (*, 99009)
99009    FORMAT (' name of the resulting binary file?')
         READ (*, 99005) file3
         WRITE(*,*) 'MAXIMUM RECORD LENGTH = ',imaxrec
         OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         OPEN (idisk2, FILE=file2, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         OPEN (idisk3, FILE=file3, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         CALL addump
         CALL inform(where)
         DO i = 1, npart
            isort(i) = i
            iorig(i) = i
         END DO
         nfullstep = 1
         CALL wdump_wrapper(idisk3)
         CLOSE (idisk3)
         CLOSE (idisk2)
         CLOSE (idisk1)
      ENDIF
c
c--Write output
c
      CALL header(where)
      CALL prout(where)
c
c--Terminate run
c
      CALL endrun

      RETURN
      END subroutine newplanet


c---------------------------------------------------------------


      SUBROUTINE planet_setup

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
      INCLUDE 'COMMONS/setlocal'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/setBfield'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/prdrag'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/dustfluid'

      REAL*8 inclination
      CHARACTER*1 iok, iwhat, idens
      CHARACTER*1 ien, irotref, uort, iplans, ians_rot
      INTEGER inum, inum2, nplanetesimals

      hzero = 0.
      irotpot = 0
      nlistinactive = 0

      sdprof = 0
      iplanetesimals = 0
      imigrate = 0
      eccentricity_p = 0.0
      inclination_p = 0.0
      gasdrag = .FALSE.

      iopmodel=0
      metallicity = 1.0
c
c--Initialise time
c
      gt = 0.0

      DO i = 1, idim
         iorig(i) = i
         isort(i) = i
         iunique(i) = i
      END DO

      CALL ktable

c
c--Get input file options
c
      CALL inopts
c
c--Check for consistency
c
      CALL chekopt

c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine planet_setup')
99004 FORMAT (A1)
11111 FORMAT (I1)

      CALL getime(ih, im, is, fhour)
      iseed = -1 * ((ih * 400) + (im * 17) + is)
      iseed = -6485
      WRITE(*,99500) iseed
99500 FORMAT(' The random seed is ',I5)
      rnd = ran1(iseed)

11001 FORMAT ('Model a planet (1), or an empty disc(2)?')
 1    WRITE (*, 11001)
      READ (*,11111) inum
      IF (inum.LT.1 .OR. inum.GT.2) GOTO 1

c
c--First deal with planet situations.
c

c--True for all but section case, in which case changed below.
      igeom = 10

      npart = 0
      iaccevol = 'f'

      inum2 = 0
      IF (inum.EQ.1) THEN
11002    FORMAT ('Model planet as potential (1) or sink (2)?')
 2       WRITE (*,11002)
         READ (*,11111) inum2
         IF (inum2.LT.1 .OR. inum2.GT.2) GOTO 2
         IF (inum2.EQ.1) irotpot = 1
      ENDIF
c
c--For a potential model choose between section and whole disc.
c
      IF (inum2.EQ.1 .OR. inum.EQ.2) THEN
 3       WRITE (*, 11003)
11003    FORMAT (//, ' INITIAL CLOUD GEOMETRY ', //,
     &        '      planet in disc section: 9 ',/,
     &        '      planets in whole disc : 10 ')
         READ (*, *) igeom
         IF (igeom.LT.9 .OR. igeom.GT.10) GOTO 3
      ENDIF
c
c--Route one: Planets modelled using potentials.
c
      irotref = 'n'
      planetmass(:) = 0.
      planetsemiaxis(:) = 0.
      planetradius(:) = 0.
      IF (inum2.EQ.1) THEN
         WRITE (*, 11201)
         READ (*,*) numplanet
         IF (numplanet.GT.nplanetmax) THEN
            WRITE (*,*) 'ERROR - numplanet.GT.nplanetmax'
            STOP
         ENDIF
         IF (numplanet.LT.1) THEN
            WRITE (*,*) 'ERROR - numplanet must be at least 1'
            STOP
         ENDIF

         DO i = 1, numplanet
            WRITE (*, 11203) i
            READ (*,*) planetsemiaxis(i)
            
c            WRITE (*, 11204)
c            READ (*,*) eccent
            
c            WRITE (*, 11205)
c            READ (*,*) inclination
c            inclination = inclination*pi/180.

            WRITE (*,11206)
            READ (*,*) planetmass(i)

            WRITE (*,11207) i
            READ (*,*) planetradius(i)

            print *,'R_planet ',i,' initial inflated = ',
     &           planetradius(i)*pradfac(i,0.)
         END DO
c
c--Currently only usable for 1 planet
c
         hmass = planetmass(1)
         coremass_orig = planetmass(1)
         coremass = planetmass(1)

      ENDIF

c--Route two: Planet(s) modelled using a sink particle.

      IF (inum.EQ.1 .AND. inum2.EQ.2) THEN
         iaccevol = 'f'
         totptmass = 0
         WRITE (*, 11201) 
11201    FORMAT (' Enter number of planets')
         READ (*,*) numpt
         
c--Print info on sink types.
      WRITE (*,*) '----------------------------------------------------'
      WRITE (*,*) '1: accretes everything regardless of tests.'
      WRITE (*,*) '2: accretes particles if they pass tests.'
      WRITE (*,*) '3: accretes part of a particle passing tests.'
      WRITE (*,*) '4: sink with accretion radius boundary corrections'
      WRITE (*,*) '   (a) smoothing length corrections'
      WRITE (*,*) '   (b) local density gradient correction of density'
      WRITE (*,*) '   (c) local pressure gradient correction to'
      WRITE (*,*) '       pressure force'
      WRITE (*,*) '   (d) local shear viscosity correction'
      WRITE (*,*) '5: sink with a surface that prevents sink accretion.'
      WRITE (*,*) '6: removes gas but loses its properties.'
      WRITE (*,*) '----------------------------------------------------'
c     hacc    = outer radius at which particle's accretion begins
c     haccall = radius at which all particles are accreted without test
         hacc = 1.0E10

         WRITE (*, 11202) 
11202    FORMAT (' Enter type of point mass (1,2,3,4,5,6)')
         READ (*,*) ipttype
         initialptm = ipttype
      
         DO i = 1, numpt
            npart = npart + 1
            iphase(npart) = ipttype
            listpm(i) = npart
            listrealpm(npart) = i
            
            WRITE (*, 11203) i
11203       FORMAT (' Enter semi-major axis of planet ', I2)
            READ (*,*) semiaxis
            
            WRITE (*, 11204)
11204       FORMAT (' Enter eccentricity')
            READ (*,*) eccent
            
            WRITE (*, 11205)
11205       FORMAT (' Enter inclination of planet in degrees')
            READ (*,*) inclination
            inclination = inclination*pi/180.
            
            WRITE (*, 11206) 
11206       FORMAT (' Enter planet mass')
            READ (*,*) tmass
            xyzmh(4,npart) = tmass
            totptmass = totptmass + tmass
            
            ra = semiaxis*(1.0+eccent)
            rorbit_orig = ra
            
            print *, 'Planet placed at ', ra
            xyzmh(1,npart) = ra*cos(inclination)
            xyzmh(2,npart) = 0.
            xyzmh(3,npart) = ra*sin(inclination)
            
            vel = (2./ra) - (1./semiaxis)
            vxyzu(1,npart) = 0.
            vxyzu(2,npart) = vel
            vxyzu(3,npart) = 0.
            
            spinx(npart) = 0.
            spiny(npart) = 0.
            spinz(npart) = 0.               
            angaddx(npart) = 0.
            angaddy(npart) = 0.
            angaddz(npart) = 0.
            
            WRITE (*,11207) i
11207       FORMAT (' Enter radius for planet ', I2)
            READ (*,*) xyzmh(5,npart)
            hacc = MIN(hacc, xyzmh(5,npart))
            haccall = hacc
         END DO
         specang = SQRT(totptmass)
         ptmassin = totptmass
         nptmass = numpt
         iptmass = 0

         IF (initialptm.EQ.5) THEN
11208       FORMAT (' Radius of planet ', I2, ' is ', 1PE12.5)
            DO i = 1, nptmass
               WRITE (*,11208) i, xyzmh(5,i)*pradfac(0,0.)
            END DO
         ENDIF
      ENDIF

      nptmasstot = nptmass

      DO i = 1, nptmass
         rho(i) = 0.
         dgrav(i) = 0.
      END DO
      
      WRITE (*,*) 'Number of point masses created: ',nptmass




c
c--Specify Box Size
c
      rmax = 0.
      rcyl = 0.
      rmin = 0.
      rmind = 0.
      ichang = 1
      igeomorig = igeom
      iplans = 'n'

      IF (igeom.EQ.9) THEN
         iexf = 7
         ifcor = 1
         omeg0 = sqrt(gg*solarm/(5.2*au)**3)*utime
         irotref = 'y'
         ibound = 100
c         iplans = 'g'

         WRITE (*, 11301) udist
11301    FORMAT (//, ' PARTICLE SET UP ', //,
     &        ' Variation of radius of disk',/,' in units of ',
     &        1PE14.5,'(cm): variation, and the disk height to',/,
     &        'radius ratio, h/r, and the phirange')
         READ (*, *) variation, hoverr, phibound
         radiusmax = 1.0+variation
         radiusmin = 1.0-variation
         zmax = radiusmax * 10.0 * hoverr
         zmin = - zmax
         xmax = radiusmax
         ymax = radiusmax
         xmin = - xmax
         ymin = - ymax
         rmax = SQRT(radiusmax*radiusmax + zmax*zmax)

 1139    WRITE (*, 11399)
11399    FORMAT ('Include heating by planetesimals (Fortier) (y/n)?')
         READ (*,*) iok
         IF (iok.EQ.'y' .OR. iok.EQ.'Y') iandrea = 1

         WRITE (*,11401)
         READ (*,*) iplans

         IF (iplans.NE.'g' .AND. iandrea.EQ.1) THEN
            print *, "Can't have planetesimals and artificial"
            print *, "planetesimal heating. Choose again"
            GOTO 1139
         ENDIF

      ELSE IF (igeom.EQ.10) THEN
         iexf = 5
         ibound = 102

         WRITE (*, 11302) udist
11302    FORMAT (//, ' PARTICLE SET UP ', //,
     &        ' Maximum and minimum radius of disk',/,' in units of ',
     &        1PE14.5,'(cm): rmaxd, rmind,',/,
     &        ' and the disk height to radius ratio, h/r')
         READ (*, *) rcyl, rmind, faclod

c--Choose between energy and temperature profile.
         use_tprof = .false.
 1133    WRITE (*, 11303)
11303    FORMAT (' Use standard u propto r^-0.5 initial disc',/
     &        ' (u) or set a specific temperature profile (t)?')
         READ (*,99004) uort

         IF (uort.EQ.'t') THEN
            use_tprof = .true.
            WRITE (*, 11304)
            READ (*,*) tprof
         ELSEIF (uort.NE.'u') THEN
            goto 1133
         ENDIF

11304    FORMAT ('Please choose a temperature profile for the disc',/
     &        'i.e. T propto r^tprof, where tprof = -1')

c--Specify density profile of disc.
         WRITE (*, 11305)
11305    FORMAT ('Please choose a surface density profile for the',/
     &        'disc i.e. sigma propto r^sdprof, where sdprof=-0.5')
         READ (*,*) sdprof

         hoverr = faclod
         zmax = rcyl * faclod
         zmin = - zmax
         xmax = rcyl
         ymax = rcyl
         xmin = - xmax
         ymin = - ymax
         rmax = SQRT(rcyl*rcyl + zmax*zmax)

c--Set planetesimal choice
         WRITE (*,11401)
11401    FORMAT ('Do you want: a gas disc (g)',/,
     &        '                an N-body planetesimal disc (p)',/,
     &        '                both a gas and planetesimal disc (b)')
         READ (*,*) iplans
      ENDIF

      deltax = xmax - xmin
      deltay = ymax - ymin
      deltaz = zmax - zmin
      rmax2 = rmax * rmax
      rmind2 = rmind * rmind
      rcyl2 = rcyl * rcyl
      zmax2 = zmax * zmax

      nreassign = 0
      nkill = 0
      naccrete = 0
      anglostx = 0.
      anglosty = 0.
      anglostz = 0.
c
c--Determine Distribution
c
      idist = 3 !random particle arrangement.

c--Always enter number of particles.
      IF (iplans.EQ.'b' .OR. iplans.EQ.'B') WRITE(*,11404)
11404 FORMAT ('You have opted to model gas and planetesimals.',/,
     &     ' This number will define the total no. of',/,
     &     ' particles for gas; idim must be large enough',/,
     &     ' for both gas and plantesimal particles.')

c
c--Specify number of particles to use.
c
 1145 WRITE (*, 11405)
11405 FORMAT(//,' Enter number of particles to use in ',
     &     ' the simulation ')
      READ (*,*) np

      h1 = (deltax*deltay*deltaz/np) ** (1./3.)
      WRITE (*,*) 'Particle spacing, h1 = ', h1

      facx = 1.
      facy = 1.
      facz = 1.
      dely = 0.
      delx = 0.

      IF (np.GT.idim) THEN
         WRITE (*, 11406) np, idim
11406    FORMAT (' Total number of particles needed :', I8, /,
     &           ' this number EXCEEDS the dimensions set to ', I8,/,
     &           ' try again.')
         GOTO 1145
      ENDIF

      WRITE (*, 11407) np
11407 FORMAT (' Total number of particles set :', I8)

      IF (igeom.EQ.10) THEN
         totvol = (pi * deltaz * (rcyl2 - rmind2))
      ELSEIF (igeom.EQ.9) THEN
         totvol = (pi*2.0*hoverr*(radiusmax**2-radiusmin**2))
      ELSE
         WRITE(*,*) 'ERROR igeom'
         CALL quit(0)
      ENDIF

c
c--make sure npart is initialised to zero
c
      npart = 0

c--Avoid NaNs by setting ekcle to non-zero value.
      DO i = 1, iradtrans
         ekcle(3,i) = 1.0
      END DO

c
c--Set particle distribution.
c      
      IF (iplans.EQ.'p') THEN
         IF (igeom.EQ.10) CALL planetesimals(np, nptmass)
         IF (igeom.EQ.9) THEN
            print *, 'Not implemented, not likely to be.'
            CALL quit(0) 
         ENDIF
      ELSE
         idens = 'd'
         icoord = 2
         CALL cyldis(igeom,np)
      ENDIF
      n1 = npart
      n2 = 0
c
c--Set Total Mass
c
      IF (igeom.EQ.9) THEN
         WRITE (*,99114)
99114    FORMAT(' Enter signorm at 5.2 au (units of 75 g/cm^2) ')
         READ (*, *) signorm
         rmindisc = 1.0-variation
         rmaxdisc = 1.0+variation
         totmas = 4.0/3.0*pi*75.0*signorm*(5.2*au)**2*
     &        (rmaxdisc**1.5-rmindisc**1.5)/(pi/phibound)/umass
         partm = totmas/(npart - nptmass)
         flowrate = 75.0*SQRT(gg*solarm*5.2*au)*((2./3.*rmindisc**
     &        (3./2.))+(2./3.*rmaxdisc**(3./2.))-LOG(rmindisc)-
     &        LOG(rmaxdisc)-(4./3.))/umass*utime*signorm/partm

         IF (iplans.EQ.'b') THEN

            iplanetesimal_info = 'p'
            n_planetesimal_types = 1

            DO i = nptmass + 1, npart
               iphase(i) = 0
               rho(i) = 0.
               xyzmh(4,i) = partm * disfrac(i)
               dgrav(i) = 0.
            END DO

            iplanetesimals = 2
c--Injection as I have it requires equal numbers of gas and dust.
            nplanetesimals = npart
c            WRITE (*,11094)
c 1184       READ (*,*) nplanetesimals
            IF (npart+nplanetesimals.GT.idim) THEN
               print *, 'idim too small for gas and planetesimals.'
               print *, 'idim, npart, nplanetesimals'
               print *, idim, npart, nplantesimals
               CALL quit(0) 
            ENDIF
            
            ntot = npart
            CALL planetesimals(nplanetesimals, ntot)
            
            WRITE (*, 11093)
            READ (*,*) psize
            r_planetesimals(1) = psize*1.E5/udist
            WRITE (*, 11092)
            READ (*,*) pdensity
            rho_planetesimals(1) = pdensity/udens
            
            WRITE (*, 11097)
            gasdrag = .true.
            WRITE (*, 11098)
            READ (*,*) idragscheme

            
            WRITE (*, 11996)
            READ (*,*) sgratio
            totgrainmass = sgratio*totmas
            
            partm2 = totgrainmass/nplanetesimals
            planetesimalmass = partm2
            print *, 'Total solids/gas mass ratio ', totgrainmass/totmas
            totmas = totmas + totgrainmass
            print*,'Planetesimal particle mass = ',
     &           totgrainmass/nplanetesimals               


            DO i = npart - nplanetesimals + 1, npart
               rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               rhotemp = ((75.*udist**2/umass)*signorm*rtemp**sdprof*
     &              exp(-(xyzmh(3,i)**2/(2.*(0.05*rtemp)**2))))/
     &              (sqrt(2.*pi)*0.05)
               xyzmh(5,i) = 1.2*(partm2/rhotemp)**(1./3.)
               iphase(i) = 11
               xyzmh(4,i) = partm2
            ENDDO
            
         ENDIF

      ELSEIF (igeom.EQ.10) THEN
         IF (iplans.EQ.'g' .OR. iplans.EQ.'G') THEN
            WRITE (*,99114)
            READ (*, *) signorm

            IF (abs(sdprof+2.0).LT.tiny) THEN
               totmas = (2.*pi*75.0*signorm*udist**2)*(log(rcyl)
     &              - log(rmind))/umass
            ELSE
c--udist in top line substituted for 5.2*au (assumes density set at r=1.0)
               totmas = (2.*pi*75.0*signorm/(udist**sdprof*(2.+
     &              sdprof)))*(rcyl**(2.+sdprof)-rmind**(2.+sdprof))
     &              *udist**(2.+sdprof)/umass
            ENDIF
            partm = totmas/(npart - nptmass)

            DO i = nptmass + 1, npart
               rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               rhotemp = ((75.*udist**2/umass)*signorm*rtemp**sdprof*
     &              exp(-(xyzmh(3,i)**2/(2.*(0.05*rtemp)**2))))/
     &              (sqrt(2.*pi)*0.05)
c               rhotemp = (75.*signorm*sqrt(1./rtemp)*udist**2/umass)/
c     &              (6.*0.05*rtemp)
               xyzmh(5,i) = 1.2*(partm/rhotemp)**(1./3.)
c
c--Alter radius for M6 Qunitic kernel (without this is gets too many
c     neighbours initially)
c
               IF (cnormk.LT.0.5/pi) xyzmh(5,i) = xyzmh(5,i)/3.0
            END DO
         ELSE
c
c--iplanetesimals = 0 - gas only, 1 - planetesimals only, 2 - both.
c
            iplanetesimals = 1
            IF (iplans.EQ.'b' .OR. iplans.EQ.'B') THEN

               iplanetesimal_info = 'p'
               n_planetesimal_types = 1

               iplanetesimals = 2
               WRITE (*,99114)
               READ (*, *) signorm
               
               IF (abs(sdprof+2.0).LT.tiny) THEN
                  totmas = (2.*pi*75.0*signorm*udist**2)*(log(rcyl)
     &                 - log(rmind))/umass
               ELSE
                  totmas = (2.*pi*75.0*signorm/(udist**sdprof*
     &                 (2.+sdprof)))*(rcyl**(2.+sdprof)-rmind**(2.+
     &                 sdprof))*udist**(2.+sdprof)/umass
               ENDIF
               partm = totmas/(npart - nptmass)
               
               DO i = nptmass + 1, npart
                  rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
                  rhotemp = ((75.*udist**2/umass)*signorm*rtemp**sdprof*
     &                 exp(-(xyzmh(3,i)**2/(2.*(0.05*rtemp)**2))))/
     &                 (sqrt(2.*pi)*0.05)
                  xyzmh(5,i) = 1.2*(partm/rhotemp)**(1./3.)

                  iphase(i) = 0
                  rho(i) = 0.
                  xyzmh(4,i) = partm * disfrac(i)
                  dgrav(i) = 0.
               END DO

               WRITE (*,11094)
11094          FORMAT ('Enter number of planetesimal particles ')
 1194          READ (*,*) nplanetesimals
               IF (npart+nplanetesimals.GT.idim) THEN
                  print *, 'idim too small for gas and planetesimals.'
                  print *, 'Try smaller number of planetesimals'
                  GOTO 1194
               ENDIF

               ntot = npart
               CALL planetesimals(nplanetesimals, ntot)
               
               WRITE (*, 11093)
11093          FORMAT ('Enter radius of plantesimals in km')
               READ (*,*) psize
               r_planetesimals(1) = psize*1.E5/udist
               WRITE (*, 11092)
11092          FORMAT ('Enter density of plantesimals (g/cm3)')
               READ (*,*) pdensity
               rho_planetesimals(1) = pdensity/udens

               WRITE (*, 11097)
               gasdrag = .true.

               WRITE (*, 11098)
               READ (*,*) idragscheme

               write (*,*) 'Drag scheme = ', idragscheme

11097          FORMAT ('It is assumed you want gas-planetesimal',/
     &              'interaction (i.e. gas drag). This can be',/
     &              'turned off in the ifile.')
11098          FORMAT ('Select a drag scheme:',/
     &              ' 0 - Perret & Murray-Cley, continuous scheme',/
     &              ' 1 - Brasser & Bains hybrid (Accumulation paper)',/
     &              ' 2 - Laibe & Price 2012 scheme')
            ENDIF

            IF (iplans.EQ.'p' .OR. iplans.EQ.'P') THEN
               WRITE (*, 11995)            
11995          FORMAT ('Enter total solids mass (Jupiter masses)')
               READ (*,*) totgrainmass
               totgrainmass = totgrainmass*1.8986E30/umass

               partm = totgrainmass/(npart-nptmass)
               totmas = totgrainmass
               print*,'Planetesimal particle mass = ',
     &              totgrainmass/(npart-nptmass)
            ENDIF

            IF (iplans.EQ.'b' .OR. iplans.EQ.'B') THEN
               WRITE (*,11996)
11996          FORMAT ('Enter solids/gas surface density ratio')
               READ (*,*) sgratio

               IF (abs(sdprof+2.0).LT.tiny) THEN
                  gasmass = (2.*pi*75.0*signorm*udist**2)*
     &                 (log(max_rplan) - log(min_rplan))/umass
                  print *, 'Using alpha = -2'
                  print *, 'Gassmass = ', gasmass
               ELSE
                  gasmass = (2.*pi*75.0*signorm/(udist**sdprof*
     &                 (2.+sdprof)))*(max_rplan**(2.+sdprof) - 
     &                 min_rplan**(2.+sdprof))*udist**(2.+sdprof)/umass

c                  gasmass = 2.*pi*(75.*udist**2/umass)*signorm*
c     &                 (max_rplan**(2.-sdprof) - min_rplan**(2.-sdprof))/
c     &                 (2.-sdprof)
                  print *, 'Gassmass = ', gasmass
               ENDIF
               totgrainmass = sgratio*gasmass

               partm2 = totgrainmass/nplanetesimals
               print *, 'Total solids/gas mass ratio ', totgrainmass/
     &              totmas
               totmas = totmas + totgrainmass
               print*,'Planetesimal particle mass = ',
     &              totgrainmass/nplanetesimals

            ENDIF

            IF (inum.EQ.1 .AND. inum2.EQ.1 .AND. numplanet.EQ.1) THEN
               WRITE (*,11095)
11095          FORMAT ('Do you want the planet to migrate to',/
     &              ' capture planetesimals in resonance?')
               READ (*,*) iok
               IF (iok.EQ.'y') THEN
                  IF (inum2.EQ.2 .AND.eccent.NE.0.0) THEN
                     print *, 'Migrating a sink with a non-circular'
                     print *, ' orbit is not currently implemented'
                     CALL quit(0) 
                  ENDIF
                  
                  imigrate = 1
                  WRITE (*,11096)
11096     FORMAT ('Enter migration rate in AU/Myr (+ve outwards)')
                  READ (*,*) pmrate
                  pmrate = pmrate*(1.49598E13/3.1536E13)*(utime/udist)
                  print *, 'Rate in code units = ', pmrate

                  iok = 'n'
                  WRITE (*,11104)
11104             FORMAT ('Should the migration stop at some radius?')
                  READ (*,*) iok
                  IF (iok.EQ.'y') THEN
                     WRITE (*,11105)
11105     FORMAT ('Enter radius at which to stop (code units)')
                     READ (*,*) rorbitmax
                     rorbitmax = (rorbitmax - planetsemiaxis(1))/pmrate
c--rorbitmax is stored as a time, at which point migration should stop.
                  ENDIF
               ENDIF
            ENDIF

            WRITE (*, 11994)
11994       FORMAT ('Radiation pressure and PR drag (y/n)?')
            READ (*,*) iok
            prswitch = .FALSE.
            IF (iok.EQ.'y') THEN
               prswitch = .TRUE.
               WRITE (*, 11998)
11998          FORMAT ('Enter a value for beta (Burns et al. 1979)')
               READ (*,*) pbeta
            ENDIF

            IF (iplans.EQ.'p' .OR. iplans.EQ.'P') THEN
               DO i = nptmass + 1, npart
                  rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2 +
     &                 xyzmh(3,i)**2)
                  rhotemp = ((totmas/totvol)*sqrt(1./rtemp))/
     &                 (6.*0.05*rtemp)
                  xyzmh(5,i) = 1.2*(partm/rhotemp)**(1./3.)
                  iphase(i) = 11
                  xyzmh(4,i) = partm
               ENDDO
            ELSEIF (iplans.EQ.'b' .OR. iplans.EQ.'B') THEN
               DO i = npart - nplanetesimals + 1, npart

                  rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
                  rhotemp = ((75.*udist**2/umass)*signorm*rtemp**sdprof*
     &                 exp(-(xyzmh(3,i)**2/(2.*(0.05*rtemp)**2))))/
     &                 (sqrt(2.*pi)*0.05)
                  xyzmh(5,i) = 1.2*(partm2/rhotemp)**(1./3.)
                  iphase(i) = 11
                  xyzmh(4,i) = partm2
               ENDDO
            ENDIF
         ENDIF
      ELSE
         WRITE (*,*) 'Error in mass setting, unknown igeom'
         CALL quit(0) 
      ENDIF

      rhozero = totmas/totvol
      WRITE(*,*) ' rhozero = ',rhozero,' rhozero(g/cc) = ',rhozero*udens
c
c--Set Particle Masses for gas cases
c
      IF (iplans.EQ.'g' .OR. iplans.EQ.'G') THEN
         DO i = nptmass + 1, npart
            iphase(i) = 0
            rho(i) = 0.
            xyzmh(4,i) = partm * disfrac(i)
            dgrav(i) = 0.
         END DO
      ENDIF
c     
c--Check if distribution is ok
c
      WRITE (*, 99020) npart - nptmass
      IF (npart-nptmass.LE.0) THEN 
         WRITE(*,99121)
         CALL quit(0)
      ENDIF
99020 FORMAT (1X, I8, ' particles have been set')
99121 FORMAT(' : ZERO PARTICLES SET')

c
c--Set e.o.s. related quantities
c
      IF (ibound.EQ.102 .AND. initialptm.EQ.5) THEN
99021    FORMAT ('For migration with a planet surface, momentum',/
     &        'conservation near the sink particle is key.',/
     &        'This is improved by setting a uniform timestep',/
     &        'for particles in the region.',/
     &        'Please choose the fraction of the Hill radius',/
     &        'out to which timesteps will all be the minimum:')
         WRITE (*, 99021)
         READ (*,*) uniformtslim
         uniformtslim2 = uniformtslim*uniformtslim
      ENDIF

      varsta = 'intener'
      ien = 'n'
      rhocrt = rhocrit * udens
      rhocrt2 = rhocrit2 * udens
      rhocrt3 = rhocrit3 * udens

 616  WRITE (*, 99030) rhocrt, rhocrt2, rhocrt3
99030 FORMAT (' Equation of state/Energy calculation:',/,
     &     '   Isothermal,     p=2/3*u*rho   (i)',/,
     &     '   Adiabatic,      p=2/3*u*rho   (a)',/,
     &     '   Polytropic,     p=A*rho^gamma (p)',/,
     &     '   Radiative,      p=R*rho*u/Cv  (r)',/,
     &     '   Variable gamma, p=A*rho^gamma (v)',/,
     &     '      critical rho (s) = ', 1PE14.5, 1PE14.5, 1PE14.5)
      READ (*, 99004) encal
      
      IF (encal.EQ.'i' .AND. (ibound.EQ.102.OR.ibound.EQ.100)) THEN
         gamma = 5.0/3.0  !Assuming temperature near planet < 100K.
      ELSE IF (encal.EQ.'p') THEN
         WRITE (*,99032)
99032    FORMAT (' Enter gamma')
         READ (*,*) gamma
         gm1 = gamma - 1.0
         WRITE (*,99033)
99033    FORMAT (' Do you want to specify the polytropic K? (y/n)',/,
     &        '(otherwise uses initial temp/internal energy/vsound)')
         READ (*,*) iok
         IF (iok.eq.'y'.or.iok.EQ.'Y') THEN
            WRITE(*,99133)
99133       FORMAT (' Enter polytropic K:')         
            READ (*,*) RK2
            RK2 = RK2*1.5
            WRITE(*,*) 'RK2 = ',RK2
         ELSE            
            RK2 = thermal/(rhozero**gm1)
         ENDIF            
      ELSE IF (encal.EQ.'a' .OR. encal.EQ.'c') THEN
         gamma = 5.0/3.0
         gm1 = gamma - 1.0
         IF (ien.EQ.'s') thermal = vsoundin2/(gamma*gm1)
         RK2 = thermal/(rhozero**gm1)
      ELSE IF (encal.EQ.'i' .OR. encal.EQ.'t') THEN
         gamma = 1.0
         RK2 = thermal
         IF (ien.EQ.'s') thermal = 1.5*vsoundin2
         RK2 = thermal
         tempiso = 2./3.*thermal/(Rg/gmw/uergg)
         WRITE(*,*) 'isothermal temperature = ',tempiso
c           spit out sound speed as used in eospg
         vsoundin2 = 2./3.*thermal
         WRITE(*,*) 'sound speed = ',SQRT(vsoundin2),
     &        ' (physical units=',SQRT(vsoundin2)*udist/utime,')'
      ELSE IF (encal.EQ.'v') THEN
c
c--Value of gamma is irrelevant for definition of variable e.o.s.
c
         gamma = 5.0/3.0
         gm1 = gamma - 1.0
         IF (ien.EQ.'s') thermal = vsoundin2*1.5 !!!/(gamma*gm1)
         RK2 = thermal/(rhozero**gm1)
         tempiso = 2./3.*thermal/(Rg/gmw/uergg)
         WRITE(*,*) 'isothermal temperature = ',tempiso
c           spit out sound speed as used in eospg
         vsoundin2 = (2./3.*RK2*rhozero**gm1)
         WRITE(*,*) 'sound speed = ',SQRT(vsoundin2),
     &        ' (physical units=',SQRT(vsoundin2)*udist/utime,')'
      ELSE IF (encal.EQ.'x') THEN
c
c--Value of gamma is irrelevant for definition of physical e.o.s.
c
         gamma = 5.0/3.0
         gm1 = gamma - 1.0
         IF (ien.EQ.'s') thermal = vsoundin2/(gamma*gm1)
         RK2 = thermal/(rhozero**gm1)
c
c--Radiative transfer
c
      ELSE IF (encal.EQ.'r') THEN
         gamma = 5.0/3.0
         gm1 = gamma - 1.0
         IF (ien.EQ.'s') thermal = vsoundin2/(gamma*gm1)
         RK2 = thermal/(rhozero**gm1)
         
         IF (ibound/10.EQ.10 .AND. ibound.NE.100) THEN
            WRITE (*,99053)
99053       FORMAT('Enter radiative transfer tolerance: ')
            READ (*,*) tolerance_rt
         ELSE
            WRITE (*,99035)
99035       FORMAT('Enter radiative transfer tolerance, boundary'//
     &          ' temperature, and scale height at which disk becomes'//
     &          ' optically thin:')
            READ (*,*) tolerance_rt, boundtemp, bounddens
         ENDIF
         
         WRITE (*, 89019)
89019    FORMAT ('Choose an opacity denominator')
         READ (*,*) opdenom
      ELSE
         GOTO 616
      ENDIF

      DO i = 1, npart
         IF (igeom.EQ.10 .AND. use_tprof) THEN
            radius = SQRT(xyzmh(1,i)**2+xyzmh(2,i)**2+xyzmh(3,i)**2)
            vxyzu(4,i)=(hoverr**2*radius**tprof)/
     &           (gamma-1.0)
         ELSEIF (igeom.EQ.9 .OR. igeom.EQ.10) THEN
            vxyzu(4,i)=hoverr**2/
     &           (SQRT(xyzmh(1,i)**2+xyzmh(2,i)**2+xyzmh(3,i)**2)*
     &           (gamma-1.0))
         ENDIF
      END DO
c
c--Set Center of Mass at Zero
c
      cmx = 0. 
      cmy = 0.
      cmz = 0.
      totmasn1 = 0.
      DO i = 1,n1
         cmx = cmx + xyzmh(4,i)*xyzmh(1,i)
         cmy = cmy + xyzmh(4,i)*xyzmh(2,i)
         cmz = cmz + xyzmh(4,i)*xyzmh(3,i)
         totmasn1 = totmasn1 + xyzmh(4,i)
      END DO
      cmx = cmx/totmasn1
      cmy = cmy/totmasn1
      cmz = cmz/totmasn1
      WRITE (*,*) 'Centre of mass (obj 1) is at: ',cmx,cmy,cmz

      WRITE (*, 99045)
99045 FORMAT(' Do you want to set the centre of mass',/,
     & ' at zero ? (y/n) ')
      READ (*, 99004) iok
      IF (iok.EQ.'y' .OR. iok.EQ.'Y') THEN
         DO i = 1, n1
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
c--Magnetic field setup once we know rhozero and thermal
c
      IF (imhd.EQ.idim) call setBfield
c
c--Set Velocities
c
      iok = 'n'
      IF (iplans.EQ.'g') THEN
         WRITE (*, 99048)
99048    FORMAT (' Do you want keplerian rotation :  (k)',/,
     &        '        or  solid body rotation :  (s)',/,
     &        '       or differential rotation :  (d)',/,
     &        '      or rotation perpendicular :  (p)',/,
     &        '      or Galactic rotation curve :  (v)',/,
     &        '        or no internal rotation :  (n)')
         READ (*, 99004) iok
      ENDIF

      IF (iok.EQ.'s' .OR. iok.EQ.'d' .OR. iok.EQ.'p') THEN
         WRITE (*,99050)
99050    FORMAT (' Enter angular velocity in rad/sec at radius=1')
         READ (*,*) angvel
         angvel = angvel * utime
         IF (iok.EQ.'p' .OR. iok.EQ.'d') THEN
            WRITE (*,99054) 
99054       FORMAT('     What rotation profile (omega)?',/,
     &           '        omega ~ exponential (e) ?',/,
     &           '        omega ~    1 over r (r) ?',/,
     &           '        omega ~  1 over r^2 (s) ?')
            READ (*, 99004) iwhat
         ENDIF
      ELSEIF (iexf.EQ.5 .OR. iexf.EQ.6 .OR. iexf.EQ.7. .OR.
     &        iok.EQ.'k') THEN
         WRITE (*,55504)
55504    FORMAT (' Enter mass for external forces')
         READ (*,*) xmass
         centralmass = xmass
c
c--Need to re-scale thermal energies if the central mass is not unity
c
         DO i = 1, npart
            vxyzu(4,i) = centralmass*vxyzu(4,i)
         END DO

         IF (inum2.NE.1) THEN
            DO i = 1, nptmass
               vxyzu(2,listpm(i)) = SQRT(xmass*vxyzu(2,listpm(i)))
               print *, 'Planet ', i, ' vel ', vxyzu(2,listpm(i))
            ENDDO
         ENDIF
         IF (prswitch) prcoefficient = prcoefficient*xmass
         IF (iplans.EQ.'p') THEN
            sqxmass = sqrt(xmass)
            print *, 'Sq mass ', sqxmass
            DO i = nptmass+1, npart
               vxyzu(1,i) = sqxmass*vxyzu(1,i)
               vxyzu(2,i) = sqxmass*vxyzu(2,i)
               vxyzu(3,i) = sqxmass*vxyzu(3,i)
            ENDDO
         ELSEIF (iplans.EQ.'b') THEN
            sqxmass = sqrt(xmass)
            print *, 'Sq mass ', sqxmass
            DO i = npart - nplanetesimals + 1, npart
               vxyzu(1,i) = sqxmass*vxyzu(1,i)
               vxyzu(2,i) = sqxmass*vxyzu(2,i)
               vxyzu(3,i) = sqxmass*vxyzu(3,i)
            ENDDO
            IF (iexf.EQ.7) THEN
               IF (irotref.NE.'y' .AND. irotref.NE.'Y') THEN
                  print *, 'iexf=7 but not rotating frame? You crazy?'
                  stop
               ENDIF
               
               print *, '##############################################'
               print *, 'WARNING: I am uncertain about changing the'
               print *, '         rotation frame for eccentric/inclined'
               print *, '         orbits of planetesimals.'
               print *, '##############################################'

c--below assumes r = 1.
               p_omega = sqxmass
               DO i = npart - nplanetesimals + 1, npart
                  rr = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
                  vphi = (-xyzmh(2,i)*vxyzu(1,i) + vxyzu(2,i)*
     &                 xyzmh(1,i))/rr
                  vr = (vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i))
     &                 /rr
                  omega = vphi/rr
                  omega = omega - p_omega
                  vphi = omega*rr
                  phi = ATAN(xyzmh(2,i)/xyzmh(1,i))

                  vxyzu(1,i) = vr*cos(phi) - vphi*sin(phi)
                  vxyzu(2,i) = vr*sin(phi) + vphi*cos(phi)
               ENDDO
            ENDIF
         ENDIF
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

      IF (iok.EQ.'v') THEN
      angvel=-2.16e+07
      angvel = angvel * utime/udist
c Set velocity dispersion parameter
        disp=5.
      ENDIF

      IF (iplans.EQ.'b') THEN
         DO i = nptmass + 1, npart - nplanetesimals
            gg2 = 1.
            radius = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
            xtild = xyzmh(1,i)/radius
            ytild = xyzmh(2,i)/radius
            velsub = 0.0
            IF (irotref.EQ.'y' .OR. irotref.EQ.'Y') THEN
               velsub = 1.0
            ENDIF
            alpha = SQRT(gg2 * xmass/radius)
            vxyzu(1,i) = -alpha * ytild + xyzmh(2,i)*velsub
            vxyzu(2,i) =  alpha * xtild - xyzmh(1,i)*velsub
            vxyzu(3,i) = 0.
         ENDDO
      ELSEIF (iplans.NE.'p') THEN
         DO i = nptmass + 1, n1
            IF (iok.EQ.'k') THEN
               gg2 = 1.  
               radius = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               xtild = xyzmh(1,i)/radius
               ytild = xyzmh(2,i)/radius
               velsub = 0.0
               IF (irotref.EQ.'y' .OR. irotref.EQ.'Y') THEN
                  velsub = 1.0
               ENDIF
               alpha = SQRT(gg2 * xmass/radius)
               vxyzu(1,i) = -alpha * ytild + xyzmh(2,i)*velsub
               vxyzu(2,i) =  alpha * xtild - xyzmh(1,i)*velsub
               vxyzu(3,i) = 0.
            ELSEIF (iok.EQ.'s') THEN
               vxyzu(1,i) = -angvel * xyzmh(2,i) * (fracrotoffset +
     &              xyzmh(3,i)*fracrotgrad)
               vxyzu(2,i) =  angvel * xyzmh(1,i) * (fracrotoffset + 
     &              xyzmh(3,i)*fracrotgrad)
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
               vxyzu(1,i) = - tanmag*xyzmh(2,i)/rxy
               vxyzu(2,i) = tanmag*xyzmh(1,i)/rxy
               vxyzu(3,i) = 0.
            ELSEIF (iok.EQ.'d') THEN 
               rad2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
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
            ELSEIF (iok.EQ.'v') THEN 
c circular velocities
               radius2 = xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i)   
               vcirc=SQRT(radius2/(1.0+radius2))
               phi=ATAN(xyzmh(2,i)/xyzmh(1,i))
               
               IF (xyzmh(1,i).GT.0 .and. xyzmh(2,i).GT.0) THEN
                  vxyzu(1,i) = -angvel*vcirc*SIN(phi)  
                  vxyzu(2,i) = angvel*vcirc*COS(phi)
                  vxyzu(3,i) = 0.
               ELSEIF (xyzmh(1,i).LE.0 .and. xyzmh(2,i).GT.0) THEN
                  vxyzu(1,i) = angvel*vcirc*SIN(phi)  
                  vxyzu(2,i) = -angvel*vcirc*COS(phi)
                  vxyzu(3,i) = 0.
               ELSEIF (xyzmh(1,i).LE.0 .and. xyzmh(2,i).LE.0) THEN
                  vxyzu(1,i) = angvel*vcirc*SIN(phi)  
                  vxyzu(2,i) = -angvel*vcirc*COS(phi)
                  vxyzu(3,i) = 0.
               ELSEIF (xyzmh(1,i).GT.0 .and. xyzmh(2,i).LE.0) THEN
                  vxyzu(1,i) = -angvel*vcirc*SIN(phi)  
                  vxyzu(2,i) = angvel*vcirc*COS(phi)
                  vxyzu(3,i) = 0.
               ENDIF            
               
c add Gaussian random velocity dispersion
c For 5 km/s dispersion, use Gaussian with sigma=disp=5
               
 250           vx1=40.*(ran1(1)-0.5)
               prob=EXP(-(vx1/disp)**2./2.)
               
               IF (prob.LT.ran1(1)) GOTO 250
               
               vxyzu(1,i)=vxyzu(1,i)+vx1*100000.*utime/udist
               
 260           vy1=40.*(ran1(1)-0.5)
               prob=EXP(-(vy1/disp)**2./2.)
               
               IF (prob.LT.ran1(1)) GOTO 260
               vxyzu(2,i)=vxyzu(2,i)+vy1*100000.*utime/udist
               
 270           vz1=40.*(ran1(1)-0.5)
               prob=EXP(-(vz1/disp)**2./2.)
               
               IF (prob.LT.ran1(1)) GOTO 270
               vxyzu(3,i)=vxyzu(3,i)+vz1*100000.*utime/udist
            ELSE
               vxyzu(1,i) = 0.
               vxyzu(2,i) = 0.
               vxyzu(3,i) = 0.
            ENDIF
         END DO
      ENDIF

c
c--option to set dtmax as fraction of free-fall time
c
      tcomp = SQRT((3 * pi) / (32 * rhozero))
      WRITE(*,*) 'The free-fall time in code units is ',tcomp
      WRITE (*, 78045)
78045 FORMAT(' Set dtmax as fraction of free-fall time? (y/n)')
      READ (*, 99004) iok
      IF (iok.EQ.'y' .OR. iok.EQ.'Y') THEN
         WRITE (*,78046)
78046    FORMAT(' Enter fraction:')
         READ (*,*) dtfrac
         dtmax = dtfrac*tcomp
         WRITE(*,*) 'setting dtmax = ',dtmax
      ENDIF
c
c--IF one-fluid dust, enter initial dust-to-gas ratio
c
      IF (idustFluid.EQ.1) THEN
         WRITE (*, 99043)
99043    FORMAT(' Enter initial dust to gas ratio (rho_d/rho_g)')
         READ (iread, *) dust_to_gas
         dust_epsilon = 1.0/(1.0+1.0/dust_to_gas)
      ENDIF
c
c--adjust smoothing lengths (also calculates initial density)
c  MUST be done (to get density) if evolving B/rho
c
      nactive = npart
      WRITE (*, 99044)
99044 FORMAT(' Do you want to adjust smoothing length',/,
     + ' to have similar number of neighbours ? (y/n) ')
      READ (*, 99004) iok

c Set abundances of main species
      IF (iener.EQ.4) THEN
      PRINT*,'entering in abundances'
      DO i=1,npart
        h2ratio(i)=0.
        abhpq(i)=0.
        abeq(i)=0.
        abHIq(i)=1.
        abco(i)=0.
      END DO
      END IF

      CALL preset(1)

      IF (ibound.EQ.100) THEN
         gapfac = 0.0
         CALL gapfactor(variation, hmass, gapfac)
      ENDIF
c
c--Set initial dust to gas ratio (to avoid divide by zero)
c
      IF (idustFluid.EQ.1) dustvar(1:npart) = 1.0

      IF (iok.EQ.'y' .OR. iok.EQ.'Y' .OR. imhd.EQ.idim .OR.
     &     idustFluid.EQ.1) CALL hcalc
c
c--Update initial dust to gas ratio
c
      IF (idustFluid.EQ.1) THEN
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               dustvar(i) = SQRT(dust_epsilon*rho(i))
            ENDIF
         END DO
      ENDIF

      IF (iok.NE.'y' .AND. iok.NE.'Y' .AND. imhd.NE.idim .AND. 
     &     idustFluid.NE.1)
     &     print *, 'Warning: uset being called without rho set'
      IF (igeom.EQ.10 .AND. use_tprof) CALL uset

      IF (imhd.NE.idim) THEN
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

      iuniquemax = npart

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
      END subroutine planet_setup
