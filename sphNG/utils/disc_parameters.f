      PROGRAM images
c***********************************************************
c                                                          *
c     This program reads sphNG dump files and extracts     *
c     the parameters of gaseous discs surrounding each     *
c     sink particle.                                       *
c                                                          *
c     Compile like: make mpi=no openmp=yes disc            *
c                                                          *
c***********************************************************

      INCLUDE 'idim'

      INTEGER*4 ilinex

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/varmhd'

      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      COMMON /aux  / hmin
      COMMON /gridv/ ilinex,iliney
      COMMON /origi/ cmx, cmy, cmz
      COMMON /files/ filein(10000)
      COMMON /phase2/ sinksize, sinkrho, denshigh, ihigh
      COMMON /numfile/ nfile

      PARAMETER (ndiscmax = 1000000)   ! max number of particles in a disc
      PARAMETER (nsinkmax =    1000)   ! max number of systems that can be analysed
      INTEGER*8 iunique_sink(nsinkmax),iunique_companion(3,nsinkmax)
      DIMENSION ilocalsink(nsinkmax),icompanion(3,nsinkmax)
      DIMENSION ilocalsinkindx(nsinkmax),ilocal2name(nsinkmax)
      DIMENSION indisc(nsinkmax)
      DIMENSION angmom(3,nsinkmax),discradius(13,nsinkmax)
      DIMENSION discmass(nsinkmax),discmasstot(nsinkmax)
      DIMENSION Bave(nsinkmax),Bmin(nsinkmax),Bmax(nsinkmax)
      DIMENSION Bup(nsinkmax),Blo(nsinkmax)
      DIMENSION Brpz_d(3,nsinkmax),Brpz_bg(3,nsinkmax),Brpz_di(3)
      DIMENSION radsearchmax(nsinkmax)
      DIMENSION nsink_mult(nsinkmax), node_mult(nsinkmax)
      DIMENSION listdone(nsinkmax)
      DIMENSION formtime(nsinkmax),formtime_unique(nsinkmax)
      DIMENSION nodelist(nsinkmax),pmassnode(nsinkmax),
     &     node_components(nsinkmax),node_comp_list(100,nsinkmax)
      DIMENSION nnode_mult_comp(nsinkmax)
      DIMENSION node_mult_comp(2,nsinkmax)
      DIMENSION xnode_axis(nsinkmax),xnode_ecc(nsinkmax)
      DIMENSION xnode_plane(3,nsinkmax)

      LOGICAL*1 ipartindisc(idim)

      DIMENSION iclist(6)
      DIMENSION cmtot(3,nsinkmax),cmtoti(3)
      DIMENSION vtot(3,nsinkmax)

      DIMENSION values(13),Bnax(6)

      DIMENSION listdisc(nsinkmax,ndiscmax),listnotdisc(ndiscmax)
      DIMENSION Bdisc(ndiscmax)
      SAVE listdisc,listnotdisc
C$OMP THREADPRIVATE(listdisc,listnotdisc)

      DIMENSION distance2(idim),indx(idim)
      SAVE distance2,indx
C$OMP THREADPRIVATE(distance2,indx)

      DIMENSION xyzm(4,2),vxyz(3,2)
      SAVE xyzm,vxyz
C$OMP THREADPRIVATE(xyzm,vxyz)

      CHARACTER*7  filein
      CHARACTER*3  prefix
      CHARACTER*21 contour1, contour2
      CHARACTER*30 contour3
      CHARACTER*10 conname
      CHARACTER*1  icent,icor,imove,isecrot
      CHARACTER*1  irepeat, ihigh, idt, idump, irot

      CHARACTER*4 ivalue
      CHARACTER*9 iuvalue

      LOGICAL     time_evol,print_pts,use_mhd,fexist
      LOGICAL     append2file,start_appending

      DATA pi/3.141592654/
c
c--set input parameters
c
      time_evol = .TRUE.     ! if true, must start from beginning of simulation or inlcude sink file
      print_pts = .FALSE.    ! if true, will print all the particles belonging to a disc to file
      append2file = .TRUE.   ! if true, will append to file, else will overwrite
      use_mhd   = .FALSE.    ! Will perform loops to calculate magnetic properties
      rhothresh_cgs = 1.d-14 ! Will only include gas particles in the disc with density greater than this [CLUSTER]
      rhothresh_cgs = 1.d-13 ! Will only include gas particles in the disc with density greater than this [ISOLATED STAR]
      rhothresh_cgs = 0.0    ! Will only include gas particles in the disc with density greater than this [DEFAULT]
      dthresh_au = 2000.0    ! Will only include gas particles in the disc with distance less than this  [DEFAULT]
      ethresh   =  0.3       ! Will only include gas particles in the disc with eccentricity less than this [DEFAULT]
      varmhd    = 'Brho'
      encal     = 'r'
c
c--read options
c
      DO i = 1, idim
         iorig(i) = i
         isort(i) = i
      END DO

      WRITE (*,*) 'Enter number of input files (or 1st in range*-1):'
      READ (*, *) numberfiles
      IF (numberfiles.LE.0) THEN
         WRITE (*,*) 'Enter the final file number in the ranges:'
         READ (*, *) numberfiles_f
         WRITE (*,*) 'Enter the simulation prefix:'
         READ (*, *) prefix
      ENDIF

      IF (numberfiles.EQ.1) THEN
         WRITE (*,*) 'DO YOU WANT TO DUMP DISCS AS A BINARY FILE?'
         READ (*,*) idump
         ifulldump = 0
         nfullstep = 1

         WRITE (*,*) 'DO YOU WANT TO ROTATE?'
         READ (*,*) irot
         IF (irot.EQ.'y') THEN
            WRITE (*,230)
 230        FORMAT('Enter start angles of viewing (0-360,0-180)?',/,
     &         '    First angle:  rotate in x-y plane, anticlockwise',/,
     &        '    Second angle: rotate from in y-z plane after 1st rot'
     &      ,/,'    Third angle: rotate from in x-z plane after 2nd rot'
     &      ,/,'  (i.e. 0,0,0 is looking down on the x-y plane')
            READ (*,*) angle1, angle2, angle3
            angle1rads = angle1*pi/180.
            angle2rads = angle2*pi/180.
            angle3rads = angle3*pi/180.
         ENDIF
      ENDIF

      IF (numberfiles.GT.0) THEN
         WRITE (*,*) 'Enter names of input files:'
         DO k = 1, numberfiles
            READ(*,90) filein(k)
 90         FORMAT(A7)
         END DO
      ELSE
         istart = abs(numberfiles)
         numberfiles = abs(numberfiles_f)-abs(numberfiles)+1
         DO k = 1, numberfiles
            WRITE(filein(k),'(a,I4.4)') prefix,istart+k-1
         END DO
      ENDIF

c     Resort the list to remove all nonexistent files
      k = 1
      DO WHILE (k.LE.numberfiles)
         inquire(file=trim(filein(k)),exist=fexist)
         IF (.not.fexist) THEN
            DO j = k,numberfiles-1
               filein(j) = filein(j+1)
            ENDDO
            numberfiles = numberfiles - 1
         ELSE
            k = k + 1
         ENDIF
      ENDDO

 444  ifile = 1
      izero = 0
      nunique_sink = 0
      start_appending = .FALSE.

      DO k = 1,numberfiles

         nfile=k
         ifile=1
         IF (start_appending) append2file = .TRUE.

         xeye = 0. ! JHW: This value was never defined.
         IF (xeye.NE.0.) THEN
            contour1 = 'CCL' // filein(nfile)
            contour2 = 'CCR' // filein(nfile)
         ELSE
            contour1 = 'CCD' // filein(nfile)
            contour2 = 'CCT' // filein(nfile)
         ENDIF

         OPEN (UNIT = 11, FILE = filein(nfile), FORM = 'unformatted')
c     
c--process file dumps
c
         IF (ifile.LE.9) THEN
            READ(contour1, 99108) conname
99108       FORMAT(A10)
            WRITE(contour1, 99002) conname, izero, izero, ifile
            READ(contour2, 99108) conname
            WRITE(contour2, 99002) conname, izero, izero, ifile
         ELSEIF (ifile.LE.99) THEN
            READ(contour1, 99108) conname
            WRITE(contour1, 99007) conname, izero, ifile
            READ(contour2, 99108) conname
            WRITE(contour2, 99007) conname, izero, ifile
         ELSE
            READ(contour1, 99108) conname
            WRITE(contour1, 99008) conname, ifile
            READ(contour2, 99108) conname
            WRITE(contour2, 99008) conname, ifile
         ENDIF
99002    FORMAT(A10, I1, I1, I1)
99007    FORMAT(A10, I1, I2)
99008    FORMAT(A10, I3)
c
c--skip files
c
 10      CONTINUE
         CALL rdump(11, ichkl, 0)
         CLOSE (11)

         print *,'Units are ',umassi, udisti, nptmass

         umass = umassi
         udist = udisti
         utime = utimei
         umagfd = umagfdi
         rhothresh = rhothresh_cgs/(umass/udist**3)
         dthresh   = dthresh_au*1.496d13 / udisti
c
c--Rotate if required
c
         IF (irot.EQ.'y') THEN
            IF (angle1rads.NE.0.0 .OR.
     &        angle2rads.NE.0.0 .OR.
     &        angle3rads.NE.0.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,vxyzu,angle1rads,angle1rade,npart)
C$OMP& shared(angle2rads,angle2rade,angle3rads,angle3rade,frac)
C$OMP& private(i,r,th1,th2,th3)
               DO i = 1, npart
c--First rotation
                  r = SQRT(xyzmh(1,i)*xyzmh(1,i)+xyzmh(2,i)*xyzmh(2,i))

                  th1 = ATAN2(xyzmh(2,i),xyzmh(1,i))

                  th1 = th1 + angle1rads
                  xyzmh(1,i) = r*COS(th1)
                  xyzmh(2,i) = r*SIN(th1)

                  r = SQRT(vxyzu(1,i)*vxyzu(1,i)+vxyzu(2,i)*vxyzu(2,i))

                  th1 = ATAN2(vxyzu(2,i),vxyzu(1,i))

                  th1 = th1 + angle1rads
                  vxyzu(1,i) = r*COS(th1)
                  vxyzu(2,i) = r*SIN(th1)
c--Second rotation
                  r = SQRT(xyzmh(2,i)*xyzmh(2,i)+xyzmh(3,i)*xyzmh(3,i))

                  th2 = ATAN2(xyzmh(3,i),xyzmh(2,i))

                  th2 = th2 + angle2rads
                  xyzmh(2,i) = r*COS(th2)
                  xyzmh(3,i) = r*SIN(th2)

                  r = SQRT(vxyzu(2,i)*vxyzu(2,i)+vxyzu(3,i)*vxyzu(3,i))

                  th2 = ATAN2(vxyzu(3,i),vxyzu(2,i))

                  th2 = th2 + angle2rads
                  vxyzu(2,i) = r*COS(th2)
                  vxyzu(3,i) = r*SIN(th2)
c--Third rotation
                  r = SQRT(xyzmh(1,i)*xyzmh(1,i)+xyzmh(3,i)*xyzmh(3,i))

                  th3 = ATAN2(xyzmh(3,i),xyzmh(1,i))

                  th3 = th3 + angle3rads
                  xyzmh(1,i) = r*COS(th3)
                  xyzmh(3,i) = r*SIN(th3)

                  r = SQRT(vxyzu(1,i)*vxyzu(1,i)+vxyzu(3,i)*vxyzu(3,i))

                  th3 = ATAN2(vxyzu(3,i),vxyzu(1,i))
                  th3 = th3 + angle3rads
                  vxyzu(1,i) = r*COS(th3)
                  vxyzu(3,i) = r*SIN(th3)
               END DO
C$OMP END PARALLEL DO
               
               OPEN (UNIT = 11, FILE = 'ROTATE', FORM = 'unformatted')
               CALL wdump(11)
               CLOSE (11)
            ENDIF
         ENDIF
c     
c--Convert gas to some other particle type so that as discs are identified
c     they can be turned back to gas to allow dump file of discs to be
c     written.
c
         DO i = 1, npart
            IF (iphase(i).EQ.0) iphase(i)=10
         END DO

         IF (k.EQ.1 .and. time_evol) THEN
            IF (nptmass.GT.1) THEN
               IF (nptmass.LT.6) THEN
                  PRINT *,'Found ',nptmass,' sinks'
                  DO i = 1, nptmass
                     PRINT *,'Found sink ',iunique(listpm(i))
                  END DO
               ENDIF
               
               OPEN (19,FILE='Names')
               READ (19,*,END=677) nptmass_in,nunique_sink
               READ (19,*,END=677) (iunique_sink(iii),
     &              iii=1,nunique_sink)
               READ (19,*,END=677) (formtime_unique(iii),
     &              iii=1,nunique_sink)
               CLOSE(19)
               GOTO 679

 677           PRINT *,'ERROR - need to start with <=1 sinks'
               STOP
            ELSE
               nunique_sink = nptmass
               IF (nunique_sink.EQ.1) THEN
                  iunique_sink(1) = iunique(listpm(1))
                  formtime_unique(1) = gt
               ENDIF
            ENDIF
         ELSE
            IF (nptmass.GT.nunique_sink) THEN
               DO iii = 1, nptmass
c
c--Check not already known
c
                  DO jjj = 1, nunique_sink
                     IF (iunique_sink(jjj).EQ.iunique(listpm(iii)))
     &                    GOTO 678
                  END DO
c
c--Add unknown sinks to list
c
                  nunique_sink = nunique_sink + 1
                  iunique_sink(nunique_sink) = iunique(listpm(iii))
                  formtime_unique(nunique_sink) = gt
 678              CONTINUE
               END DO
            ENDIF
         ENDIF
 679     CONTINUE
c         IF (nptmass.NE.nunique_sink) THEN
         IF (nptmass.GT.nunique_sink) THEN
            PRINT *,'ERROR - nptmass.GT.nunique_sink ',nptmass,
     &           nunique_sink
            STOP
         ENDIF

         OPEN (19,FILE='Names',STATUS='replace')
         WRITE(19,*) nptmass,nunique_sink
         WRITE(19,*) (iunique_sink(iii),
     &        iii=1,nunique_sink)
         WRITE(19,*) (formtime_unique(iii),
     &        iii=1,nunique_sink)
         CLOSE(19)

         PRINT *,'Number of sinks ',nptmass,nunique_sink
         PRINT *,'Names of sinks ',(iunique_sink(iii),
     &        iii=1,nunique_sink)
c
c--Find local index
c
         DO jjj = 1, nptmass
            DO iii = 1, nunique_sink
               IF (iunique(listpm(jjj)).EQ.iunique_sink(iii)) THEN
                  ilocalsink(iii) = listpm(jjj)
                  ilocalsinkindx(iii) = jjj
                  ilocal2name(jjj) = iii
                  formtime(jjj) = formtime_unique(iii)
                  GOTO 680
               ENDIF
            END DO
            PRINT *,'ERROR - not found ',iii,iunique_sink(iii)
            STOP
 680        CONTINUE
            WRITE (*,*) 'Identified ',jjj,' as ',iii,iunique_sink(iii)
         END DO
c
c--Values for disc radii
c
         values(1) = 0.02
         values(2) = 0.05
         values(3) = 0.10
         values(4) = 0.20
         values(5) = 0.30
         values(6) = 0.40
         values(7) = 0.50
         values(8) = 0.632
         values(9) = 0.70
         values(10) = 0.80
         values(11) = 0.90
         values(12) = 0.95
         values(13) = 1.00
c
c--Process sinks
c
c        Initialise arrays
         ipartindisc = .FALSE.
         nsink_mult = 1
         discradius = 0.
         indisc     = 0
         Bmin       = 1.0d30
         Brpz_d     = 0.
         Brpz_bg    = 0.
         Bave       = 0.
         Bmax       = 0.
         print*, 'arrays initialised.  nptmass = ',nptmass
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nunique_sink,ilocalsink,xyzmh,vxyzu,Bxyz,node_components)
C$OMP& shared(node_comp_list,cmtot,vtot,discmass,angmom,npart,indisc)
C$OMP& shared(icompanion,iunique_companion,iphase,nsink_mult,iunique)
C$OMP& shared(iunique_sink,radsearchmax,udisti,ipartindisc,ilocal2name)
C$OMP& shared(discmasstot,discradius,values,nptmass,listpm,use_mhd)
C$OMP& shared(print_pts,Bmin,Bave,Bmax,Bup,Blo,Brpz_d,Brpz_bg,umagfdi)
C$OMP& shared(rhothresh,dthresh,ethresh,rho,umass,udist)
C$OMP& private(isink,i,j,l,iii,jjj,iptcur,xsink,ysink,zsink,sinkmass)
C$OMP& private(dx,dy,dz,ipart,xipart,yipart,zipart,dvx,dvy,dvz,Bi)
C$OMP& private(Bdisc,Brpz_di,kn,kx,rmass)
C$OMP& private(ij,idisc_id,ns,n0,cmtoti)
C$OMP& private(ndisc,nnotdisc,nsinks,distcomp2,jpart,r2,jjj_mindist)
C$OMP& private(totalmass,dist,etot,semimajor,eccentricity)
C$OMP& private(radapastron,ipos)
C$OMP& private(jval,jjval)
         DO isink = 1, nptmass
            iptcur = listpm(isink)
            xsink = xyzmh(1,iptcur)
            ysink = xyzmh(2,iptcur)
            zsink = xyzmh(3,iptcur)
            sinkmass = xyzmh(4,iptcur)
            node_components(isink) = 1
            node_comp_list(1,isink) = ilocal2name(isink)
c
c--Initialise centre of mass to that of sink, and for multiple
c
            cmtot(1:3,isink) = sinkmass*xyzmh(1:3,iptcur)
            vtot(1:3,isink)  = sinkmass*vxyzu(1:3,iptcur)

            print *,'Sink location, mass ',xsink,ysink,zsink,sinkmass

            discmass(isink) = 0.
            DO l = 1, 3
               angmom(l,isink) = 0.
            END DO
            ndisc = 0
            nnotdisc = 0
            listdisc(isink,:) = 0
c
c--Compute distance^2 of all particle from sink
c
            DO i = 1, npart
               dx = xyzmh(1,i)-xsink
               dy = xyzmh(2,i)-ysink
               dz = xyzmh(3,i)-zsink
               distance2(i) = dx**2 + dy**2 + dz**2
            END DO
c
c--Sort by distance
c
            print *,'sorting ',isink
            CALL indexx2(npart,distance2,indx)

c            print *,'sorted ',indx(1),indx(2),distance2(indx(1)),
c     &           distance2(indx(2)),iptcur
c
c--Set multiplicity value of this sink, and zero companion info
c
            nsinks = 1
            DO j = 1, 3
               icompanion(j,isink) = 0
               iunique_companion(j,isink) = 0
            END DO
c
c--Begin loop over all particles to find discs and companions
c
            DO j = 1, npart
               i = indx(j)
c
c--Skip self
c
c               print *,'  ',isink,i,nsinks
               IF (i.EQ.iptcur) CYCLE
c
c--Keep track of other sinks.  Increase centre of mass and centre of mass
c     velocity to include other sinks.
c
               IF (iphase(i).GT.0 .AND. iphase(i).LE.9) THEN
c
c--Signify that isink has another sink particle within cut-off radius
c     (nominally 2000 AU).  This allows us to specify what we mean by
c     'single star' -- one that does not have another star within 
c     2000.  The other star does NOT have to be bound.  It is a
c     definition similar to an observer's definition.
c
                  nsinks = nsinks + 1
                  nsink_mult(isink) = nsinks
c
c--Find identification of sink
c
                  DO iii = 1, nunique_sink
                     IF (iunique(i).EQ.iunique_sink(iii)) THEN
c
c--Determine if these two sinks are mutual nearest neighbours.  This
c     is the nearest neighbour of isink, but isink may not be the 
c     nearest neighbour of iii.  If they are, they may be a binary 
c     system (if they are bound).
c
c
c--Find distance^2 of closest sink to iii
c
                        ipart = ilocalsink(iii)
                        xipart = xyzmh(1,ipart)
                        yipart = xyzmh(2,ipart)
                        zipart = xyzmh(3,ipart)
                        distcomp2 = 1.0E+30
                        DO jval = 1, nptmass
                           DO jjval = 1, nunique_sink
                  IF (iunique(listpm(jval)).EQ.iunique_sink(jjval)) THEN
                                 jjj = jjval
                                 GOTO 789
                              ENDIF
                           END DO
                        WRITE (*,*) 'iunique(jval).NE.iunique_sink(jjj)'
                           STOP
 789                       CONTINUE
                           IF (jjj.NE.iii) THEN
                              jpart = ilocalsink(jjj)
                              r2 = (xyzmh(1,jpart)-xipart)**2 + 
     &                             (xyzmh(2,jpart)-yipart)**2 +
     &                             (xyzmh(3,jpart)-zipart)**2
                              IF (r2.LT.distcomp2) THEN
                                 distcomp2 = r2
                                 jjj_mindist = jjj
                              ENDIF
                           ENDIF
                        END DO
c
c--This is the test to see if they are mutual nearest neighbours or not
c     Cannot do the test for boundedness here if we want to include the
c     circustellar disc mass because the other sink particle may not
c     have been done yet.
c     NOTE: If the companion number is NEGATIVE, then the companion
c     is not the mutual nearest neighbour, but we record its information
c     anyway (with a minus sign).
c
                        IF (jjj_mindist.EQ.ilocal2name(isink)) THEN
                           icompanion(nsinks-1,isink) = iii
                           iunique_companion(nsinks-1,isink) =
     &                          iunique(i)
                        ELSE
                           icompanion(nsinks-1,isink) = - iii
                           iunique_companion(nsinks-1,isink) =
     &                          - iunique(i)
                        ENDIF
                        GOTO 880
                     ENDIF
                  END DO
                  PRINT *,'ERROR - sink miss ',i,iunique(i),nunique_sink
                  STOP
 880              CONTINUE
c
c--The exit ends the computation of the disc for this sink particle.
c     The assumption is that a circumstellar disc cannot exist beyond 
c     the distance of another sink particle.
c
                  radsearchmax(isink) = SQRT(distance2(i))
                  EXIT
               ENDIF
c
c--Set total mass to consider whether new particle is bound or not
c
               totalmass = sinkmass + discmass(isink)

c               print *,'  ',j,xsink,cmtot(1,isink)/totalmass,
c     &              vxyzu(1,iptcur),vtot(1,isink)/totalmass,
c     &              nsinks,iphase(i)
c
c--Find parameters of circumstellar disc (nsinks = 1 only)
c
               IF (nsinks.EQ.1 .AND. iphase(i).EQ.10) THEN
                  dx = xyzmh(1,i) - cmtot(1,isink)/totalmass
                  dy = xyzmh(2,i) - cmtot(2,isink)/totalmass
                  dz = xyzmh(3,i) - cmtot(3,isink)/totalmass
                  dist = SQRT(dx**2 + dy**2 + dz**2)
c
c--Make distance cut: dthresh = 2000 AU
c
c                  print *,'dist = ',dist*udisti/1.496E+13
                  IF (dist.GT.dthresh) THEN
                     print *,' exiting A ',isink
                     radsearchmax(isink) = SQRT(distance2(i))
                     EXIT
                  ELSE
c
c--Calculate relative energy and orbit
c
                     xyzm(4,1) = xyzmh(4,i)
                     xyzm(4,2) = totalmass
                     xyzm(1:3,1) = xyzmh(1:3,i)
                     xyzm(1:3,2) = cmtot(1:3,isink)/totalmass
                     vxyz(1:3,1) = vxyzu(1:3,i)
                     vxyz(1:3,2) = vtot(1:3,isink)/totalmass
                     CALL find_orbit(xyzm,vxyz,etot,semimajor,
     &                    eccentricity)
                     radapastron = semimajor*(1.0+eccentricity)
c
c--Add to disc if passes tests
c
c                     print *,i,dist,eccentricity,
c    &                    radapastron*udisti/1.496E+13,
c    &                    etot,semimajoraxis

                     IF (eccentricity.LT.ethresh     .AND. 
     &                    radapastron.LT.dthresh     .AND.
     &                    rho(i).     GE.rhothresh ) THEN
                        IF (i.GT.idim) THEN
                           PRINT *,'i.GT.idim'
                           STOP
                        ELSE
                           IF (ipartindisc(i)) THEN
                              PRINT *,'C-S Disc particle conflict ',
     &                             isink,i
                              CYCLE
                           ENDIF
                           ipartindisc(i) = .TRUE.
                        ENDIF
                        
                        iphase(i) = 0

                        discmass(isink) = discmass(isink) + xyzmh(4,i)
                        cmtot(1:3,isink) = cmtot(1:3,isink) + 
     &                       xyzmh(1:3,i)*xyzmh(4,i)
                        vtot(1:3,isink) = vtot(1:3,isink) + 
     &                       vxyzu(1:3,i)*xyzmh(4,i)

                        dvx = vxyzu(1,i) - vtot(1,isink)/totalmass
                        dvy = vxyzu(2,i) - vtot(2,isink)/totalmass
                        dvz = vxyzu(3,i) - vtot(3,isink)/totalmass

                        angmom(1,isink) = angmom(1,isink) +
     &                       xyzmh(4,i)*(dvy*dz - dvz*dy)
                        angmom(2,isink) = angmom(2,isink) +
     &                       xyzmh(4,i)*(dvz*dx - dvx*dz)
                        angmom(3,isink) = angmom(3,isink) +
     &                       xyzmh(4,i)*(dvx*dy - dvy*dx)
                        ndisc = ndisc + 1
                        IF (ndisc.GT.ndiscmax) THEN
                           WRITE (*,*) 'ndisc.GT.ndiscmax'
                           STOP
                        ENDIF

                        listdisc(isink,ndisc)  = i

                     ELSE IF (radapastron*udisti.LT.1.496E+13*200.) THEN
                        nnotdisc = nnotdisc + 1
                        IF (nnotdisc.GT.ndiscmax) THEN
                           WRITE (*,*) 'nnotdisc.GT.ndiscmax'
                           STOP
                        ENDIF
                        listnotdisc(nnotdisc) = i
                     ENDIF
                  ENDIF
               ENDIF
            END DO  ! End of the loop over SPH particles

            discmasstot(isink) = discmass(isink)
c
c--Set values of discradius that contain certain fractions of the total
c     disc mass
c
            IF (ndisc.GT.0) THEN
               DO iii = 1, 13
                  ipos = MAX(1,INT(values(iii)*ndisc))
                  IF (ipos.NE.0) THEN
                     discradius(iii,isink) =
     &                    SQRT(distance2(listdisc(isink,ipos)))
                  ELSE
                     discradius(iii,isink) = 0.
                  ENDIF
               END DO
c              Determine which particles are in the 0.632M disc
               DO iii = 1,ndisc
                  rmass = xyzmh(4,jjj)
                  jjj = listdisc(isink,iii)
                  dx = xyzmh(1,jjj) - cmtot(1,isink)/totalmass
                  dy = xyzmh(2,jjj) - cmtot(2,isink)/totalmass
                  dz = xyzmh(3,jjj) - cmtot(3,isink)/totalmass
                  dist = SQRT(dx**2 + dy**2 + dz**2)
                  if (dist.LE.discradius(8,isink)) then
                     indisc(isink) = indisc(isink) + 1
                     listdisc(isink,iii) = -abs(listdisc(isink,iii))
                  endif
               enddo
c              Calculate the magnetic field of the 0.632M disc
               ij = 0
               DO iii = 1,ndisc
                  jjj = abs(listdisc(isink,iii))
                  if (listdisc(isink,iii) < 0) then
                     Bi  = dot_product(Bxyz(1:3,jjj),Bxyz(1:3,jjj))
                     Bi  = sqrt(Bi)
                     ij  = ij + 1
                     if (ij.le.ndiscmax) Bdisc(ij) = Bi
                     Bave(isink) = Bave(isink) + log10(Bi)
                     Bmin(isink) = min(Bi,Bmin(isink))
                     Bmax(isink) = max(Bi,Bmax(isink))
                     idisc_id = 1
                  else
                     idisc_id = 0
                  endif
                  IF (print_pts) THEN
                     WRITE (44+isink,'(7(1PE12.5,1X),I12)')
     &                 xyzmh(1:3,jjj),Bxyz(1:3,jjj)*umagfdi,
     &                 rho(i)*(umass/udist**3),idisc_id
                  ENDIF
               enddo
               IF (indisc(isink).gt.0 .and. use_mhd) THEN
c                 Calculate the upper and lower percentiles
                  ns = indisc(isink)
                  CALL indexx2(ns,Bdisc(:),indx)
                  if (ns.le.2) then
                     kn = 1
                     kx = ns
                  else
                     kn = indx(max(2,   int(0.025*ns)  ))
                     kx = indx(min(ns-1,int(0.975*ns)+1))
                  endif
                  Blo(isink) = Bdisc(kn)
                  Bup(isink) = Bdisc(kx)
c                 Calculate the B-components in the disc
                  print*, 'ndisc, n63.2 = ',ndisc,indisc(isink)
                  n0 = ndisc
                  if (totalmass.gt.0.0) then
                     cmtoti = cmtot(1:3,isink)/totalmass
                  else
                     cmtoti = 0.
                  endif
                  call update_B(idim,imhd2,n0,listdisc(isink,1:n0),
     &                          xyzmh,Bxyz,cmtoti,
     &                          angmom(1:3,isink),Brpz_d(:,isink))
               ENDIF
            ELSE
               print*, "There are no particles in the disc: ",
     &                 ndisc,indisc(isink)
            ENDIF
            print*, 'nnotdisc = ',nnotdisc
            IF (nnotdisc.GT.0 .and. use_mhd) THEN
c              Calculate Bz component of the ambient B-field
               if (totalmass.gt.0.0) then
                  cmtoti = cmtot(1:3,isink)/totalmass
               else
                  cmtoti = 0.
               endif
               n0 = nnotdisc
               call update_B(idim,imhd2,n0,-listnotdisc(1:n0),
     &                       xyzmh,Bxyz,cmtoti,
     &                       angmom(1:3,isink),Brpz_bg(:,isink))
            ENDIF

            PRINT*,'Discmass ',isink,' : ',discmass(isink),ndisc,nsinks
            PRINT*,'Disc radius ',isink,' : ',discradius(7:10,isink)
            PRINT*,'Angular momomentum ',isink,' : ',angmom(1:3,isink)
            PRINT*,' '
         END DO ! End of the loop over sink particles for c-s discs
C$OMP END PARALLEL DO
c=====================================
c
c--Now look for discs around multiples
c
c--Need to find closest bound mutual nearest neighbours.  For this, 
c     loop over all nodes (maybe single, binary, triple) recursively.
c
c--Set up initial list of nodes (sink particles and their c-s discs).
c     The initial positions and velocities come from cmtot and vtot.
c
         nnode_mult_comp = 0  ! array
         nbinary = 0
         ntriple = 0
         nquadpair = 0
         nquadhier = 0
         nnodes = nptmass
         nnodes_list = nptmass
         DO inode = 1, nnodes_list
            nodelist(inode) = inode
            pmassnode(inode) = xyzmh(4,listpm(inode))
            node_mult(inode) = 1
         END DO
c
c--Now loop over all nodes to find closest bound pair of mutual cloest
c     neighbours
c
 1500    r2min = 1.0E+30
         DO in = 1, nnodes_list
            n = nodelist(in)
            DO in2 = in+1, nnodes_list
               n2 = nodelist(in2)
               xyzm(4,1) = pmassnode(n) + discmasstot(n)
               xyzm(4,2) = pmassnode(n2) + discmasstot(n2)
               xyzm(1:3,1) = cmtot(1:3,n)/xyzm(4,1)
               xyzm(1:3,2) = cmtot(1:3,n2)/xyzm(4,2)
               vxyz(1:3,1) = vtot(1:3,n)/xyzm(4,1)
               vxyz(1:3,2) = vtot(1:3,n2)/xyzm(4,2)
               rx = xyzm(1,1) - xyzm(1,2)
               ry = xyzm(2,1) - xyzm(2,2)
               rz = xyzm(3,1) - xyzm(3,2)
               vxdiff = vxyz(1,1) - vxyz(1,2)
               vydiff = vxyz(2,1) - vxyz(2,2)
               vzdiff = vxyz(3,1) - vxyz(3,2)
               r2 = (rx**2 + ry**2 + rz**2)
               CALL find_orbit(xyzm,vxyz,etot,semimajor,ecc)

               IF (etot.LT.0.0 .AND. r2.LT.r2min .AND. 
     &              node_mult(n)+node_mult(n2).LT.5) THEN
c
c--Check for *mutual* nearest neighbours
c
                  radius_test2_keep = 1.0E+30
                  DO in3 = 1, nnodes_list
                     n3 = nodelist(in3)
                     IF (n3.NE.n2) THEN
                        totmass = pmassnode(n3) + discmasstot(n3)
                        xn3 = cmtot(1,n3)/totmass
                        yn3 = cmtot(2,n3)/totmass
                        zn3 = cmtot(3,n3)/totmass
                        radius_test2 = (xyzm(1,2)-xn3)**2 + 
     &                       (xyzm(2,2)-yn3)**2 + (xyzm(3,2)-zn3)**2
                        IF (radius_test2.LT.radius_test2_keep) THEN
                           radius_test2_keep = radius_test2
                           n23_keep = n3
                        ENDIF
                     ENDIF
                  END DO
 
                  radius_test2_keep = 1.0E+30
                  DO in3 = 1, nnodes_list
                     n3 = nodelist(in3)
                     IF (n3.NE.n) THEN
                        totmass = pmassnode(n3) + discmasstot(n3)
                        xn3 = cmtot(1,n3)/totmass
                        yn3 = cmtot(2,n3)/totmass
                        zn3 = cmtot(3,n3)/totmass
                        radius_test2 = (xyzm(1,1)-xn3)**2 + 
     &                       (xyzm(2,1)-yn3)**2 + (xyzm(3,1)-zn3)**2
                        IF (radius_test2.LT.radius_test2_keep) THEN
                           radius_test2_keep = radius_test2
                           n13_keep = n3
                        ENDIF
                     ENDIF
                  END DO

                  IF (n23_keep.EQ.n .AND. n13_keep.EQ.n2) THEN
                     r2min = MIN(r2min,r2)

                     IF (pmassnode(n).GE.pmassnode(n2)) THEN
                        nkeep1 = n
                        nkeep2 = n2
                     ELSE
                        nkeep1 = n2
                        nkeep2 = n
                     ENDIF

                     semimajorkeep = semimajor
                     eccentricitykeep = ecc
c
c--Orbital plane
c
                     planex = vydiff*rz - vzdiff*ry
                     planey = vzdiff*rx - vxdiff*rz
                     planez = vxdiff*ry - vydiff*rx
                     planelength = SQRT(planex**2+planey**2+planez**2)
                     planex = planex/planelength
                     planey = planey/planelength
                     planez = planez/planelength
                  ENDIF
               ENDIF
            END DO
         END DO
         IF (r2min.EQ.1.0E+30) GOTO 1600

         PRINT *,'Found multiple ',nkeep1, nkeep2,
     &        SQRT(r2min)*udisti/(1.496E+13),
     &        semimajorkeep*udisti/(1.496E+13), eccentricitykeep
c
c--Store information for multiple system
c
         node_mult_combined = node_mult(nkeep1) + node_mult(nkeep2)
         IF (node_mult_combined.EQ.2) THEN
            nbinary = nbinary + 1
         ELSEIF (node_mult_combined.EQ.3) THEN
            ntriple = ntriple + 1
         ELSEIF (node_mult_combined.EQ.4) THEN
            IF (node_mult(nkeep1).EQ.node_mult(nkeep2)) THEN
               nquadpair = nquadpair + 1
            ELSE
               nquadhier = nquadhier + 1
            ENDIF
         ELSE
            PRINT *,'ERROR - node_mult_combined > 4 ',
     &           node_mult_combined
            STOP
         ENDIF
c
c--Store information for new multiple system in new node
c
         nnodes = nnodes + 1
         IF (nnodes.GT.nsinkmax) THEN
            PRINT *,'nnodes>nsinkmax'
            STOP
         ENDIF
         formtime(nnodes) = MIN(formtime(nkeep1),formtime(nkeep2))
         node_mult(nnodes) = node_mult_combined
         pmassnode(nnodes) = pmassnode(nkeep1) + pmassnode(nkeep2)
         discmass(nnodes) = 0.
         discmasstot(nnodes) = discmasstot(nkeep1) + discmasstot(nkeep2)
         cmtot(1:3,nnodes) = cmtot(1:3,nkeep1) + cmtot(1:3,nkeep2)
         vtot(1:3,nnodes) = vtot(1:3,nkeep1) + vtot(1:3,nkeep2)
         xnode_axis(nnodes) = semimajorkeep
         xnode_ecc(nnodes) = eccentricitykeep
         xnode_plane(1,nnodes) = planex
         xnode_plane(2,nnodes) = planey
         xnode_plane(3,nnodes) = planez
c
c--Keep list of component sink particles of node
c
         ncomponents = 0
         DO i = 1, node_components(nkeep1)
            ncomponents = ncomponents + 1
            IF (ncomponents.GT.100) THEN
               PRINT *,'ERROR - ncomponents.GT.100 a'
            ELSE
c
c--Only keep list of components up to 100 components
c     (Only multiples up to quads are considered later anyway)
c
               node_comp_list(ncomponents,nnodes) = 
     &              node_comp_list(i,nkeep1)
            ENDIF
         END DO
         DO i = 1, node_components(nkeep2)
            ncomponents = ncomponents + 1
            IF (ncomponents.GT.100) THEN
               PRINT *,'ERROR - ncomponents.GT.100 b'
            ELSE
               node_comp_list(ncomponents,nnodes) = 
     &              node_comp_list(i,nkeep2)
            ENDIF
         END DO
         node_components(nnodes) = ncomponents
c
c--IF number of components is greater than 4, stop computing more
c     information (don't need to know the sub-components or whether
c     it has circum-multiple disc).
c
         IF (ncomponents.GT.4) GOTO 8500
c
c--If node consists of sub-nodes, keep separate list of these, but only
c     for triples or quads (2+1) or (2+2) or (2+1)+1
c
         IF (ncomponents.GT.2 .AND. ncomponents.LT.5) THEN
            ncomponents = 0
            IF (node_components(nkeep1).EQ.2 .OR. 
     &           node_components(nkeep1).EQ.3) THEN
               IF (node_components(nkeep1).EQ.3) THEN
                  ncomponents = ncomponents + 1
                  IF (ncomponents.GT.2) THEN
                     PRINT *,'ERROR - nnode_mult_comp.GT.2 a'
                  ELSE
                     node_mult_comp(ncomponents,nnodes) =
     &                    node_mult_comp(1,nkeep1)
                  ENDIF
               ENDIF
               ncomponents = ncomponents + 1
               IF (ncomponents.GT.2) THEN
                  PRINT *,'ERROR - nnode_mult_comp.GT.2 b'
               ELSE
                  node_mult_comp(ncomponents,nnodes) =
     &                 nkeep1
               ENDIF
            ENDIF
            IF (node_components(nkeep2).EQ.2 .OR. 
     &           node_components(nkeep2).EQ.3) THEN
               IF (node_components(nkeep2).EQ.3) THEN
                  ncomponents = ncomponents + 1
                  IF (ncomponents.GT.2) THEN
                     PRINT *,'ERROR - nnode_mult_comp.GT.2 c'
                  ELSE
                     node_mult_comp(ncomponents,nnodes) =
     &                    node_mult_comp(1,nkeep2)
                  ENDIF
               ENDIF
               ncomponents = ncomponents + 1
               IF (ncomponents.GT.2) THEN
                  PRINT *,'ERROR - nnode_mult_comp.GT.2 d'
               ELSE
                  node_mult_comp(ncomponents,nnodes) =
     &                 nkeep2
               ENDIF
            ENDIF
            nnode_mult_comp(nnodes) = ncomponents
         ENDIF
c
c--Now compute circum-multiple disc masses.  Involves looping over all
c     particles again, but excluding those that are closer to one of
c     the sink particles than the distance between the two sinks.
c     Also exclude the component sink particles.
c
c--Initialise centre of mass to that of binary.
c     NOTE: cmtot and vtot are postion and velocity TIMES total mass.
c
         totmass = pmassnode(nnodes) + discmasstot(nnodes)

         xnode = cmtot(1,nnodes)/totmass
         ynode = cmtot(2,nnodes)/totmass
         znode = cmtot(3,nnodes)/totmass

         discmass(nnodes) = 0.
         DO l = 1, 3
            angmom(l,nnodes) = 0.
         END DO
         ndisc = 0
         nnotdisc = 0
         listdisc(nnodes,:) = 0
c
c--Compute distance^2 of all particle from centre of mass
c
         DO i = 1, npart
            dx = xyzmh(1,i) - xnode
            dy = xyzmh(2,i) - ynode
            dz = xyzmh(3,i) - znode
            distance2(i) = dx**2 + dy**2 + dz**2
         END DO
c
c--Sort by distance
c
         CALL indexx2(npart,distance2,indx)
c
c--Set multiplicity value of this sink, and zero companion info
c
         DO j = 1, 3
            icompanion(j,nnodes) = 0
            iunique_companion(j,nnodes) = 0
         END DO
c
c--Begin loop over all particles to find discs and companions
c
         nsinks = 0
         DO j = 1, npart
            i = indx(j)
c
c--Skip sink particle components of this node, and particles that are
c     already in the disc of a component
c
            IF (ipartindisc(i)) CYCLE
            DO ii = 1, node_components(nnodes)
               icomp = node_comp_list(ii,nnodes)
               IF (ilocalsink(icomp).EQ.i) GOTO 8000
            END DO
c
c--Keep track of other sinks.  Increase centre of mass and centre of 
c     mass velocity to include other sinks.
c
            IF (iphase(i).GT.0 .AND. iphase(i).LE.9) THEN
c
c--Signify that node has another sink particle within cut-off radius
c     (nominally 2000 AU).  This allows us to specify what we mean by
c     'binary star' -- one that does not have another star within 
c     2000.  The other star does NOT have to be bound.  It is a
c     definition similar to an observer's definition.
c
               nsinks = nsinks + 1
               nsink_mult(nnodes) = nsinks
c
c--Find identification of sink
c
               DO iii = 1, nunique_sink
                  IF (iunique(i).EQ.iunique_sink(iii)) THEN
                     icompanion(nsinks,nnodes) = iii
                     iunique_companion(nsinks,nnodes) = iunique(i)
                     GOTO 890
                  ENDIF
               END DO
               PRINT *,'ERROR - sink id ',i,iunique(i),nunique_sink
               STOP
 890           CONTINUE
c
c--The exit ends the computation of the disc for this node.
c     The assumption is that THIS disc cannot exist beyond 
c     the distance of another sink particle.
c
               radsearchmax(nnodes) = SQRT(distance2(i))
               EXIT
            ENDIF
c
c--Set total mass to consider whether new particle is bound or not
c
            totalmass = pmassnode(nnodes) + discmasstot(nnodes)
c
c--Find parameters of circum-multiple disc
c
            IF (iphase(i).EQ.10) THEN
               dx = xyzmh(1,i) - cmtot(1,nnodes)/totalmass
               dy = xyzmh(2,i) - cmtot(2,nnodes)/totalmass
               dz = xyzmh(3,i) - cmtot(3,nnodes)/totalmass
               dist = SQRT(dx**2 + dy**2 + dz**2)
c
c--Make distance cut - 2000 AU
c
               IF (dist.GT.dthresh) THEN
                  print *,' exiting A'
                  radsearchmax(nnodes) = SQRT(distance2(i))
                  EXIT
               ELSE
                  dvx = vxyzu(1,i) - vtot(1,nnodes)/totalmass
                  dvy = vxyzu(2,i) - vtot(2,nnodes)/totalmass
                  dvz = vxyzu(3,i) - vtot(3,nnodes)/totalmass
c
c--Calculate relative energy and angular momentum and orbit
c                  
                  divvr2 = dvx*dvx + dvy*dvy + dvz*dvz
                  vrad = (dvx*dx + dvy*dy + dvz*dz)/dist
                  vrad2 = vrad*vrad
                  vkep2 = totalmass/dist
                  vpotdif = -vkep2 + divvr2/2.
                  vtan2 = divvr2 - vrad2
                  vtan = SQRT(vtan2)
                  specangmom2 = vtan2*dist**2

                  semimajoraxis = - totalmass/(2.0*vpotdif)

                  IF (vpotdif.LT.0) THEN
                     tempnum=1.0 + 2.0*vpotdif*specangmom2/
     &                    (totalmass**2)
                     IF (tempnum.LT.0.0) THEN
                        WRITE(iprint,*) 'ERROR - tempnum'
                        WRITE(iprint,*) vpotdif,specangmom2
                        WRITE(iprint,*) vkep2,divvr2,r2,vtan2
                        STOP
                     ELSE
                        eccentricity = SQRT(tempnum)
                        radapastron = - totalmass*
     &                       (1.0+eccentricity)/(2.0*vpotdif)
                     ENDIF
                  ELSE
                     eccentricity = 1.0E+10
                     radapastron = 1.0E+10
                  ENDIF
c
c--Add to disc if passes tests
c
                  IF (eccentricity.LT.ethresh     .AND. 
     &                 radapastron.LT.dthresh     .AND.
     &                 rho(i)     .GE.rhothresh ) THEN
                     IF (i.GT.idim) THEN
                        PRINT *,'i.GT.idim'
                        STOP
                     ELSE
                        ipartindisc(i) = .TRUE.
                     ENDIF

                     iphase(i) = 11

                     discmass(nnodes) = discmass(nnodes) + xyzmh(4,i)
                     discmasstot(nnodes) = discmasstot(nnodes) + 
     &                    xyzmh(4,i)
                     cmtot(1:3,nnodes) = cmtot(1:3,nnodes) + 
     &                    xyzmh(1:3,i)*xyzmh(4,i)
                     vtot(1:3,nnodes) = vtot(1:3,nnodes) + 
     &                    vxyzu(1:3,i)*xyzmh(4,i)
                     angmom(1,nnodes) = angmom(1,nnodes) +
     &                    xyzmh(4,i)*(dvy*dz - dvz*dy)
                     angmom(2,nnodes) = angmom(2,nnodes) +
     &                    xyzmh(4,i)*(dvz*dx - dvx*dz)
                     angmom(3,nnodes) = angmom(3,nnodes) +
     &                    xyzmh(4,i)*(dvx*dy - dvy*dx)
                     ndisc = ndisc + 1
                     IF (ndisc.GT.ndiscmax) THEN
                        WRITE (*,*) 'ndisc.GT.ndiscmax'
                        STOP
                     ENDIF
                     listdisc(nnodes,ndisc)  = i
                  ELSE IF (radapastron*udisti.LT.1.496E+13*200.) THEN
                     nnotdisc = nnotdisc + 1
                     IF (nnotdisc.GT.ndiscmax) THEN
                        WRITE (*,*) 'nnotdisc.GT.ndiscmax'
                        STOP
                     ENDIF
                     listnotdisc(nnotdisc)  = i
                  ENDIF
               ENDIF
            ENDIF
 8000       CONTINUE
         END DO                 ! End of the loop over SPH particles
c
c--Set values of discradius that contain certain fractions of the total
c     disc mass
c
         IF (ndisc.GT.0) THEN
            DO iii = 1, 13
               ipos = MAX(1,INT(values(iii)*ndisc))
               IF (ipos.NE.0) THEN
                  discradius(iii,nnodes) = 
     &            SQRT(distance2(listdisc(nnodes,ipos)))
               ELSE
                  discradius(iii,nnodes) = 0.
               ENDIF
            END DO

c           calculate which particles are in the 0.632M disc
            DO iii = 1,ndisc
               jjj = listdisc(nnodes,iii)
               dx = xyzmh(1,jjj) - cmtot(1,nnodes)/totalmass
               dy = xyzmh(2,jjj) - cmtot(2,nnodes)/totalmass
               dz = xyzmh(3,jjj) - cmtot(3,nnodes)/totalmass
               dist = SQRT(dx**2 + dy**2 + dz**2)
               if (dist.LE.discradius(8,nnodes)) then
                  indisc(nnodes) = indisc(nnodes) + 1
                  listdisc(nnodes,iii) = -abs(listdisc(nnodes,iii))
               endif
            enddo
c           Calculate the magnetic field of the 0.632M disc
            ij = 0
            DO iii = 1,ndisc
               jjj = abs(listdisc(nnodes,iii))
               if (listdisc(nnodes,iii) < 0) then
                  Bi = dot_product(Bxyz(1:3,jjj),Bxyz(1:3,jjj))
                  Bi = sqrt(Bi)
                  ij = ij + 1
                  if (ij.le.ndiscmax) Bdisc(ij) = Bi
                  Bave(nnodes) = Bave(nnodes) + log10(Bi)
                  Bmin(nnodes) = min(Bi,Bmin(nnodes))
                  Bmax(nnodes) = max(Bi,Bmax(nnodes))
                  idisc_id = 2
               else
                  idisc_id = 0
               endif
               IF (print_pts) THEN
                  WRITE (44+nnodes,'(7(1PE12.5,1X),I12)')
     &               xyzmh(1:3,jjj),Bxyz(1:3,jjj)*umagfdi,
     &               rho(i)*(umass/udist**3),idisc_id
               ENDIF
            enddo

            IF (use_mhd) THEN
c              Calculate the B-components in the disc
               print*, 'ndisc, n63.2 = ',ndisc,indisc(nnodes)
               if (totalmass.gt.0.0) then
                  cmtoti = cmtot(1:3,nnodes)/totalmass
               else
                  cmtoti = 0.
               endif
               n0 = ndisc
               call update_B(idim,imhd2,n0,listdisc(nnodes,1:n0),
     &                       xyzmh,Bxyz,cmtoti,
     &                       angmom(1:3,nnodes),Brpz_d(:,nnodes))
            ENDIF
         ENDIF

         IF (nnotdisc.GT.0 .and. use_mhd) THEN
            print*, 'nnotdisc = ',nnotdisc
c           Calculate Bz component of the ambient B-field
            if (totalmass.gt.0.0) then
               cmtoti = cmtot(1:3,nnodes)/totalmass
            else
               cmtoti = 0.
            endif
            n0 = nnotdisc
            call update_B(idim,imhd2,n0,-listnotdisc(1:n0),xyzmh,Bxyz,
     &                    cmtoti,angmom(1:3,nnodes),Brpz_bg(:,nnodes))
         ENDIF

         PRINT *,'Discmass ',nnodes,discmass(nnodes),ndisc
         PRINT *,'Disc radius ',discradius(7:10,nnodes)
         PRINT *,'Angular momomentum ',angmom(1:3,nnodes)
c
c--Remove the two component nodes that have been merged into a multiple
c     system from the list of active nodes
c
 8500    nnodes_list = nnodes_list + 1
         nodelist(nnodes_list) = nnodes

         nnodes_list_new = 0
         DO j = 1, nnodes_list
            i = nodelist(j)
            IF (i.NE.nkeep1 .AND. i.NE.nkeep2) THEN
               nnodes_list_new = nnodes_list_new + 1
               nodelist(nnodes_list_new) = i
            ENDIF
         END DO
         IF (nnodes_list_new.NE.nnodes_list-2) THEN
            PRINT *,'ERROR - nnodes_nlist_ew.NE.nnodes_list-2'
            STOP
         ENDIF
         nnodes_list = nnodes_list_new
c
c--Go back and look for next multiple system
c
         GOTO 1500
c
c--Output results for this file
c
 1600    PRINT *,'File ',nfile,': Found no more multiple systems'
         PRINT *,''
         PRINT *,'File ',nfile,' has ',nbinary,' binaries (pairs)'
         PRINT *,'File ',nfile,' has ',ntriple,' triples'
         PRINT *,'File ',nfile,' has ',nquadpair,' quads (pairs)'
         PRINT *,'File ',nfile,' has ',nquadhier,' quads (hierarchical)'
         PRINT *,''
c
c--Write out information for all individual sink particles
c
         DO i = 1, nptmass
            DO jval = 1, nunique_sink
               IF (iunique(listpm(i)).EQ.iunique_sink(jval)) THEN
                  ival = jval
                  GOTO 790
               ENDIF
            END DO
            WRITE (*,*) 'Not found X ',nptmass,nunique_list,listpm(1),
     &           iunique_sink(1)
            STOP
 790        CONTINUE
            IF (ival.GT.99) THEN
               WRITE (ivalue,"(I3)") ival
            ELSEIF (ival.GT.9) THEN
               WRITE (ivalue,"('0',I2)") ival
            ELSE
               WRITE (ivalue,"('00',I1)") ival
            ENDIF
            IF (iunique_sink(ival).GT.9999999) THEN
               WRITE (iuvalue,"(I8)") iunique_sink(ival)
            ELSEIF (iunique_sink(ival).GT.999999) THEN
               WRITE (iuvalue,"('0',I7)") iunique_sink(ival)
            ELSEIF (iunique_sink(ival).GT.99999) THEN
               WRITE (iuvalue,"('00',I6)") iunique_sink(ival)
            ELSEIF (iunique_sink(ival).GT.9999) THEN
               WRITE (iuvalue,"('000',I5)") iunique_sink(ival)
            ELSEIF (iunique_sink(ival).GT.999) THEN
               WRITE (iuvalue,"('0000',I4)") iunique_sink(ival)
            ELSEIF (iunique_sink(ival).GT.99) THEN
               WRITE (iuvalue,"('00000',I3)") iunique_sink(ival)
            ELSEIF (iunique_sink(ival).GT.9) THEN
               WRITE (iuvalue,"('000000',I2)") iunique_sink(ival)
            ELSE
               WRITE (iuvalue,"('0000000',I1)") iunique_sink(ival)
            ENDIF
            contour1 = 'Disc_' // ivalue(1:3) // '_' // iuvalue(1:8) 
     &           // '.txt'

            print *,contour1
c
c--NOTE: Sink particle spins are defined as the negative of the way
c     disc and orbital angular momentum are defined in this code,
c     hence the need to multiply the spins by -1.
c
         Bnax = 0.
         IF (use_mhd) THEN
            IF (indisc(i) .GT. 0) THEN
               Bnax(1) = Bmin(i)
               Bnax(2) = 10**(Bave(i)/indisc(i))
               Bnax(3) = Bmax(i)
            ENDIF
         ENDIF

         IF (append2file) THEN
            OPEN (16,file=contour1,ACCESS='append')
         ELSE
            OPEN (16,file=contour1)
            start_appending = .TRUE.
         ENDIF
         WRITE (16,"(17(1PE12.5,1x),I4,1x,17(1PE12.5,1x),1x,2(I9,1x))") 
     &           gt,formtime(i),xyzmh(4,listpm(i)),
     &           discmass(i),
     &           (discradius(j,i)*udisti/1.496E+13,j=1,13),
     &           nsink_mult(i),
     &           (angmom(j,i),j=1,3),
     &           -spinx(i),-spiny(i),-spinz(i),
     &           Bnax(1:3)*umagfdi,Blo(i)*umagfdi,Bup(i)*umagfdi,        !25-27,28,29
     &           Brpz_d(1:3,i)*umagfdi,Brpz_bg(1:3,i)*umagfdi,           !30-32,33-35
     &           (icompanion(j,i),j=1,1),
     &           (iunique_companion(j,i),j=1,1)
            CLOSE (16)
         END DO

c         WRITE (95,'(1465(1PE12.5,1X))') gt,(discmass(i),i=1,183),
c     &        ((discradius(j,i),j=1,4),i=1,183)
c     &        ((angmom(j,i),j=1,3),i=1,183)

c
c--Write out multiple system information
c     Start with highest index and work backward, keeping a list of those
c     sinks that have already been included in an output.
c     Do not consider any systems that have more than 4 components.
c
         nlistdone = 0
         DO i = nnodes, 1, -1
            IF (node_components(i).LE.4) THEN
c
c--Output data for all pairs, regardless of whether they are in higher
c     order multiples or not
c
               IF (node_components(i).EQ.2) THEN
c
c--Write out information on the binary (actually a pair, which maybe
c     a sub-component of a higher-order system).
c     File name is Pair_N_[component_].txt
c     Where N is the number of stars in the system, and [] is repeated
c     for the N sink particle numbers that make up the system.
c
                  contour3 = 'Pair'

                  DO j = 1, node_components(i)
                     icompnumber = node_comp_list(j,i)
                     IF (icompnumber.GT.99) THEN
                        WRITE (ivalue,"(I3)") icompnumber
                     ELSEIF (icompnumber.GT.9) THEN
                        WRITE (ivalue,"('0',I2)") icompnumber
                     ELSE
                        WRITE (ivalue,"('00',I1)") icompnumber
                     ENDIF
                     iend = 4 + (j-1)*4
                     contour3 = contour3(1:iend) // '_' // ivalue(1:3)
                  END DO
                  iend = 4 + 4*node_components(i)
                  contour3 = contour3(1:iend)  // '.txt'
                  print *,contour3
c
c--Need to compute total disc quantities for the system as a whole
c
c
c--Determine primary mass
c
                  primarymass = 0.
c
c--Add in discs of individual protostars
c
                  DO j = 1, node_components(i)
                     icomp = ilocalsinkindx(node_comp_list(j,i))
                     IF (pmassnode(icomp).GT.primarymass) THEN
                        primarymass = pmassnode(icomp)
                        iprimary = icomp
                     ENDIF
                  END DO
                  IF (iprimary.EQ.
     &                 ilocalsinkindx(node_comp_list(1,i))) THEN
                     isecondary = ilocalsinkindx(node_comp_list(2,i))
                  ELSE
                     isecondary = ilocalsinkindx(node_comp_list(1,i))
                  ENDIF
c
c--Check if component of higher order multiple up to 4 components
c
                  iordermax = 2
                  DO j = i+1, nnodes
                     IF (node_components(j).LE.4) THEN
                        DO l = 1, node_components(j)
                           IF (iprimary.EQ.
     &                        ilocalsinkindx(node_comp_list(l,j))) THEN
                              iordermax = node_components(j)
                              EXIT
                           ENDIF
                        END DO
                     ENDIF
                  END DO
c
c--Determine disc radius that contains certain percentage of mass
c     This would be possible if circumbinary disc mass contains more
c     than ~1/2 of the mass of the total disc mass, because then the
c     radius would lie within the radius of the c-b disc.  However, if
c     the c-b mass is less, it doesn't make sense because the radius
c     containing 60% of the total disc mass will lie between the
c     radii of the c-s discs and the inner radius of the c-b disc.
c     In fact, the c-s discs would then be "resolved" into two separate
c     discs.
c     So better just to output the distributions of all three discs.
c

c
c--If pair, then want to write separate disc masses and angular 
c     momentum of two c-s discs and the c-b disc
c
                 IF (append2file) THEN
                    OPEN (16,file=contour3,ACCESS='append')
                 ELSE
                    OPEN (16,file=contour3)
                    start_appending = .TRUE.
                 ENDIF
                 Bnax   = 0.
                 itotal = 0
                 IF (use_mhd) THEN
                    IF (indisc(i).GT.0) THEN
                       Bnax(1) = Bmin(i)
                       Bnax(2) = 10**(Bave(i)/indisc(i))
                       Bnax(3) = Bmax(i)
                       Bnax(4) = Bmin(i)
                       Bnax(5) = Bave(i)
                       Bnax(6) = Bmax(i)
                       itotal  = indisc(i)
                    ENDIF
                    IF (indisc(iprimary).GT.0) THEN
                       Bnax(4) = min(Bnax(4),Bmin(iprimary))
                       Bnax(5) = Bnax(5)+Bave(iprimary)
                       Bnax(6) = max(Bnax(6),Bmax(iprimary))
                       itotal  = itotal + indisc(iprimary)
                    ENDIF
                    IF (indisc(isecondary).GT.0) THEN
                       Bnax(4) = min(Bnax(4),Bmin(isecondary))
                       Bnax(5) = Bnax(5)+Bave(isecondary)
                       Bnax(6) = max(Bnax(6),Bmax(isecondary))
                       itotal  = itotal + indisc(isecondary)
                    ENDIF
                    Bnax(5) = 10**(Bnax(5)/itotal)
c                   Calculate the upper and lower percentiles
                    Bdisc = 0.
                    jb    = 0
                    do bloop = 1,3
                       if (bloop==1) ibval = i
                       if (bloop==2) ibval = iprimary
                       if (bloop==3) ibval = isecondary
                       j = 1
                       do while (listdisc(ibval,j).ne.0)
                          jj = listdisc(ibval,j)
                          if (jj < 0) then
                             jb = jb + 1
                             if (jb.le.ndiscmax) then
                                Bdisc(jb) = dot_product(Bxyz(1:3,-jj),
     &                                                  Bxyz(1:3,-jj))
                             else
                                print*, 'exceeding ndiscmax (pair)',jb
                             endif
                          endif
                          j = j + 1
                       enddo
                    enddo
                    CALL indexx2(jb,Bdisc(:),indx)
                    if (jb.le.2) then
                       kn = 1
                       kx = jb
                    else
                       kn = indx(max(2,   int(0.025*jb)  ))
                       kx = indx(min(jb-1,int(0.975*jb)+1))
                    endif
                    Blo(i) = sqrt(Bdisc(kn))
                    Bup(i) = sqrt(Bdisc(kx))
                  ENDIF

      WRITE (16,"(I3,1x,61(1PE12.5,1x),4(I4,1x),14(1PE12.5,1x))")
     &              iordermax,gt,formtime(i),pmassnode(i),
     &              primarymass,
     &              discmasstot(i),discmass(i),discmass(iprimary),
     &              discmass(isecondary),
     &              (discradius(j,i)*udisti/1.496E+13,j=1,13),
     &              (discradius(j,iprimary)*udisti/1.496E+13,j=1,13),
     &              (discradius(j,isecondary)*udisti/1.496E+13,j=1,13),
     &              (xnode_axis(i)*udisti/1.496E+13),
     &              xnode_ecc(i), (xnode_plane(j,i),j=1,3),
     &              (angmom(j,i),j=1,3),
     &              (angmom(j,iprimary),j=1,3),
     &              (angmom(j,isecondary),j=1,3),
     &              (node_comp_list(j,i),j=1,4),
     &              -spinx(iprimary),
     &              -spiny(iprimary),
     &              -spinz(iprimary),
     &              -spinx(isecondary),
     &              -spiny(isecondary),
     &              -spinz(isecondary),
     &              Bnax*umagfdi,Blo(i)*umagfdi,Bup(i)*umagfdi
                  CLOSE (16)
               END IF
c--END OF PAIRS-----------------------
c
c--For systems, check that components of node do not appear in list 
c     of components already done.  In other words, if a pair is part of
c     a higher-order multiple system, it will not be written out as a
c     System because it has already been included in a higher-order
c     System.
c
               DO j = 1, node_components(i)
                  DO l = 1, nlistdone
                     IF (node_comp_list(j,i).EQ.listdone(l)) GOTO 70 
                  END DO
               END DO
c
c--Add to list of 
c
               DO j = 1, node_components(i)
                  nlistdone = nlistdone + 1
                  listdone(nlistdone) = node_comp_list(j,i)
               END DO
c
c--Write out information on the multiple system
c     File name is System_N_[component_].txt
c     Where N is the number of stars in the system, and [] is repeated
c     for the N sink particle numbers that make up the system.
c
               WRITE (ivalue,"(I1)") node_components(i)
               contour3 = 'System_' // ivalue(1:1)

               DO j = 1, node_components(i)
                  icompnumber = node_comp_list(j,i)
                  IF (icompnumber.GT.99) THEN
                     WRITE (ivalue,"(I3)") icompnumber
                  ELSEIF (icompnumber.GT.9) THEN
                     WRITE (ivalue,"('0',I2)") icompnumber
                  ELSE
                     WRITE (ivalue,"('00',I1)") icompnumber
                  ENDIF
                  iend = 8+(j-1)*4
                  contour3 = contour3(1:iend) // '_' // ivalue(1:3)
               END DO
               iend = 8+4*node_components(i)
               contour3 = contour3(1:iend)  // '.txt'
               print *,contour3
c
c--Need to compute total disc quantities for the system as a whole
c
c
c--Determine primary mass, and total disc mass and angular momentum
c
               primarymass = 0.
c
c--Add in discs of individual protostars (including itself if the
c     system is a single star).
c
               DO j = 1, node_components(i)
                  icomp = ilocalsinkindx(node_comp_list(j,i))
                  iclist(j) = icomp
                  primarymass = MAX(primarymass,pmassnode(icomp))
               END DO
c
c--Add in discs of components that are multiple (i.e. binary or triple)
c
               DO j = 1, nnode_mult_comp(i)
                  icomp = node_mult_comp(j,i)
                  iclist(node_components(i)+j) = icomp
               END DO
c
c--Pad the rest of the array so that index does not give invalid value
c
               DO j = node_components(i)+nnode_mult_comp(i)+1, 6
                  iclist(j) = nsinkmax
               END DO
c
c--Determine disc radius that contains certain percentage of mass
c             NOTE: This assumes a maximimum of a quad system, which yeilds a maximum of 7 discs
              IF (append2file) THEN
                 OPEN (16,file=contour3,ACCESS='append')
               ELSE
                  OPEN (16,file=contour3)
                  start_appending = .TRUE.
               ENDIF
               Bnax   = 0.
               itotal = 0
               IF (use_mhd) THEN
                  IF (indisc(i).GT.0) THEN
                     Bnax(1) = Bmin(i)
                     Bnax(2) = 10**(Bave(i)/indisc(i))
                     Bnax(3) = Bmax(i)
                     Bnax(4) = Bmin(i)
                     Bnax(5) = Bave(i)
                     Bnax(6) = Bmax(i)
                     Brpz_di = Brpz_d(:,i)*indisc(i)
                     itotal  = indisc(i)
                  ENDIF
c                 If this is system > 1, then discs of higher order do not include the 
c                 particles from the lower order disc; else disc(i)==disc(iclist(1))
                  IF (node_components(i).GT.1) THEN
                     DO j = 1,6
                        IF (indisc(iclist(j)).GT.0) THEN
                           itotal  = itotal + indisc(iclist(j))
                           Bnax(4) = min(Bnax(4),Bmin(iclist(j)))
                           Bnax(5) = Bnax(5)+    Bave(iclist(j))
                           Bnax(6) = max(Bnax(6),Bmax(iclist(j)))
                           Brpz_di = Brpz_di
     &                   + Brpz_d(:,iclist(j))*indisc(iclist(j))
                        ENDIF
                     ENDDO
                     IF (itotal.GT.0.0) THEN
                        Bnax(5) = 10**(Bnax(5)/itotal)
                        Brpz_di = Brpz_di/itotal
                     ENDIF
                  ELSE
                     Bnax(5) = Bnax(2)
                  ENDIF
c                 Calculate the upper and lower percentiles
                  Bdisc = 0.
                  jb    = 0
                  do bloop = 0,6
                     if (bloop==0) then
                        ibval = i
                     else
                        ibval = iclist(bloop)
                     endif
                     j = 1
                     do while (listdisc(ibval,j).ne.0)
                        jj = listdisc(ibval,j)
                        if (jj < 0) then
                           jb = jb + 1
                           if (jb.le.ndiscmax) then
                              Bdisc(jb) = dot_product(Bxyz(1:3,-jj),
     &                                                Bxyz(1:3,-jj))
                           else
                              print*, 'exceeding ndiscmax (system)',jb
                           endif
                        endif
                        j = j + 1
                     enddo
                  enddo
                  CALL indexx2(jb,Bdisc(:),indx)
                  kn = indx(max(2,   int(0.025*jb)  ))
                  kx = indx(min(jb-1,int(0.975*jb)+1))
                  Blo(i) = sqrt(Bdisc(kn))
                  Bup(i) = sqrt(Bdisc(kx))
               ENDIF

               WRITE (16,"(I2,1x,125(1PE12.5,1x),4(I4,1x))") 
     &              node_components(i),gt,formtime(i),
     &              pmassnode(i),primarymass,   ! Total mass of all the stars, mass of the most massive star.
     &              discmasstot(i),discmass(i), ! Summed mass of all 7 possible discs, the mass of only the gas surrounding the total system.
     &              discmass(iclist(1)),        ! Masses of the individual discs followd by masses of the pairs;
     &              discmass(iclist(2)),        ! arrays are filled from 1 to 6, so end values will always be zero.
     &              discmass(iclist(3)),
     &              discmass(iclist(4)),
     &              discmass(iclist(5)),
     &              discmass(iclist(6)),
     &              (discradius(j,i)*udisti/1.496E+13,j=1,13),         ! disc radii of either the total or circumsystem disc
     &              (discradius(j,iclist(1))*udisti/1.496E+13,j=1,13), ! disc radii in the same order as the masses
     &              (discradius(j,iclist(2))*udisti/1.496E+13,j=1,13),
     &              (discradius(j,iclist(3))*udisti/1.496E+13,j=1,13),
     &              (discradius(j,iclist(4))*udisti/1.496E+13,j=1,13),
     &              (discradius(j,iclist(5))*udisti/1.496E+13,j=1,13),
     &              (discradius(j,iclist(6))*udisti/1.496E+13,j=1,13),
     &              (xnode_axis(i)*udisti/1.496E+13),
     &              xnode_ecc(i), (xnode_plane(j,i),j=1,3),
     &              (angmom(j,i),j=1,3),      ! Angular momentum of the circumsystem disc
     &              Bnax*umagfdi,             ! 1:3=Bfields of the circumsystem disc;4:6=Bfields of the total disc
     &              Blo(i)*umagfdi,Bup(i)*umagfdi,
     &              Brpz_di*umagfdi,Brpz_bg(1:3,i)*umagfdi,
     &              (node_comp_list(j,i),j=1,4)
               CLOSE (16)
            ENDIF
 70         CONTINUE
         END DO

         IF (idump.EQ.'Y' .OR. idump.EQ.'y') THEN
            DO i = 1, npart
               IF (iphase(i).EQ.10) iphase(i)=-1
            END DO
            
            OPEN (UNIT = 11, FILE = 'DISCDMP', FORM = 'unformatted')
            CALL wdump(11)
            CLOSE (11)
         ENDIF
         
      END DO

      STOP
      END

      SUBROUTINE quit
      STOP
      END


      SUBROUTINE endrun
      CALL quit
      END


      SUBROUTINE find_orbit(xyzm,vxyz,etot,semimajor,eccentricity)

      DIMENSION xyzm(4,2),vxyz(3,2)
c
c--Set total mass to consider whether two objects are bound or not
c
      totalmass = xyzm(4,1) + xyzm(4,2)
      IF (totalmass.LE.0.) THEN
         PRINT *,'ERROR - find_orbit'
         STOP
      ENDIF
c
c--Find relative positions and velocities
c
      rx = xyzm(1,1) - xyzm(1,2)
      ry = xyzm(2,1) - xyzm(2,2)
      rz = xyzm(3,1) - xyzm(3,2)
      vxdiff = vxyz(1,1) - vxyz(1,2)
      vydiff = vxyz(2,1) - vxyz(2,2)
      vzdiff = vxyz(3,1) - vxyz(3,2)
      r2 = (rx**2 + ry**2 + rz**2)
c
c--Energies and angular momenta
c
      ekin = 0.5*(vxdiff**2 + vydiff**2 + vzdiff**2)
      epot = -totalmass/SQRT(r2)
      etot = ekin + epot
      dangx = vydiff*rz - vzdiff*ry
      dangy = vzdiff*rx - vxdiff*rz
      dangz = vxdiff*ry - vydiff*rx
      dang2 = (dangx**2 + dangy**2 + dangz**2)
c
c--Orbital parameters
c
cc      IF (etot.NE.0.) THEN
      IF (etot.LE.0.) THEN
         semimajor = - totalmass/(2.0*etot)
         value = 1.0 - dang2/totalmass/semimajor
         IF (value.GE.0.0) THEN
            eccentricity = SQRT(value)
         ELSE
cc            eccentricity = -1.0
            eccentricity = 1.0
         ENDIF
      ELSE
         semimajor = 1.0E+30
         eccentricity = 1.0
      ENDIF

      RETURN
      END
c
c--Update magnetic field components in the new orientation
c
      SUBROUTINE update_B(idim,imhd2,ndisc,listdisc,xyzmh,Bxyz,cmtot,
     &                    angmom,Brpz)
      INTEGER   idim,ndisc
      DIMENSION xyzmh(5,idim),cmtot(3),Bxyz(3,imhd2),angmom(3)
      DIMENSION listdisc(ndisc)
      DIMENSION uhat(3),vhat(3),what(3)
      DIMENSION nbdisc(6),Brpz(3),Brpz0(6)
      LOGICAL   has_ang

c
      has_ang = .true.
      if (abs(angmom(1)).le.epsilon(angmom(1)) .and. 
     &    abs(angmom(2)).le.epsilon(angmom(2)) .and. 
     &    abs(angmom(2)).le.epsilon(angmom(3))) has_ang = .false.

c --Calculate the relevant vectors for angular momentum
      if (has_ang) then
         what(1) = angmom(1)
         what(2) = angmom(2)
         what(3) = angmom(3)
         what    = what/sqrt(what(1)**2+what(2)**2+what(3)**2)
         uhat(1) = 0.
         uhat(2) = what(2)
         uhat(3) = -what(2)*what(2)/what(3)
         uhat    = uhat/sqrt(uhat(1)**2+uhat(2)**2+uhat(3)**2)
         vhat(1) = what(1)
         vhat(2) = -what(1)**2*what(2)/(what(3)**2+what(2)**2)
         vhat(3) = -what(1)**2*what(3)/(what(3)**2+what(2)**2)
         vhat    = vhat/sqrt(vhat(1)**2+vhat(2)**2+vhat(3)**2)
      else
         what(1) = 0.
         what(2) = 0.
         what(3) = 1.
         uhat(1) = 1.
         uhat(2) = 0.
         uhat(3) = 0.
         vhat(1) = 0.
         vhat(2) = 1.
         vhat(3) = 0.
      endif
      nbdisc = 0
      Brpz0  = 0.
      Brpz   = 0.
      do j = 1,ndisc
         if (listdisc(j) < 0) then
            i   = -listdisc(j)
            dx  = xyzmh(1,i) - cmtot(1)
            dy  = xyzmh(2,i) - cmtot(2)
            dz  = xyzmh(3,i) - cmtot(3)
            Bx  = Bxyz(1,i)
            By  = Bxyz(2,i)
            Bz  = Bxyz(3,i)
            xp  = uhat(1)*dx + uhat(2)*dy + uhat(3)*dz
            yp  = vhat(1)*dx + vhat(2)*dy + vhat(3)*dz 
            zp  = what(1)*dx + what(2)*dy + what(3)*dz
            Bxp = uhat(1)*Bx + uhat(2)*By + uhat(3)*Bz
            Byp = vhat(1)*Bx + vhat(2)*By + vhat(3)*Bz 
            Bzp = what(1)*Bx + what(2)*By + what(3)*Bz
            r2  = xp*xp + yp*yp
            if (r2.gt.0.0) then
               Br   = (xp*Bxp + yp*Byp)/sqrt(r2)
               Bphi = (yp*Bxp - xp*Byp)/sqrt(r2)
            else
               Br   = 0.
               Bphi = 0.
            endif
            if (Br > 0.) then
               nbdisc(1) = nbdisc(1) + 1
               Brpz0(1)  = Brpz0(1)  + log10(Br)
            else
               nbdisc(2) = nbdisc(2) + 1
               Brpz0(2)  = Brpz0(2)  + log10(-Br)
            endif
            if (Bphi > 0.) then
               nbdisc(3) = nbdisc(3) + 1
               Brpz0(3)  = Brpz0(3)  + log10(Bphi)
            else
               nbdisc(4) = nbdisc(4) + 1
               Brpz0(4)  = Brpz0(4)  + log10(-Bphi)
            endif
            if (Bzp > 0.) then
               nbdisc(5) = nbdisc(5) + 1
               Brpz0(5)  = Brpz0(5)  + log10(Bzp)
            else
               nbdisc(6) = nbdisc(6) + 1
               Brpz0(6)  = Brpz0(6)  + log10(-Bzp)
            endif
            B2 = Bx**2 + By**2 + Bz**2
            Bp2 = Bxp**2 + Byp**2 + Bzp**2
            if (abs(B2-Bp2) > 100.*Bp2*epsilon(Bp2)) then
               print*, "B & Bp are not the same magnitude!"
               print*,B2,Bp2,abs(B2-Bp2),100.*Bp2*epsilon(Bp2)
               stop
            endif
            Bp2 = Br**2 + Bphi**2 + Bzp**2
            if (abs(B2-Bp2) > 100.*Bp2*epsilon(Bp2)) then
               print*, "B & Brpz are not the same magnitude!"
               print*, B2,Bp2,abs(B2-Bp2),100.*Bp2*epsilon(Bp2)
               stop
            endif
         endif
      enddo

      do i = 1,6
         if (nbdisc(i)>0) Brpz0(i) = 10**(Brpz0(i)/nbdisc(i))*nbdisc(i)
      enddo
      do i = 1,3
         if (nbdisc(2*i-1) > 0 .or. nbdisc(2*i) > 0) then 
            Brpz(i) = (Brpz0(2*i-1) - Brpz0(2*i))
     &              / (nbdisc(2*i-1)+nbdisc(2*i))
         endif
      enddo

      RETURN
      END
