      PROGRAM images
c***********************************************************
c                                                          *
c  this program generates images suitable for imagetool    *
c                                                          *
c***********************************************************

      INCLUDE 'idim'

      REAL*8 udist, umass, utime
      INTEGER*4 iline, ntallpix, nwidepix
      REAL*4 timeout, xmindensout, xmaxdensout

      PARAMETER (nwidemax = 5000)
      PARAMETER (ntallmax = 3300)
      REAL*4 global
      COMMON /gimage/ global(nwidemax,ntallmax)
      COMMON /gvals / nwide, ntall
      COMMON /gpix  / ntallpix, nwidepix

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/gtime'
      
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      
      COMMON /phase2/ sinksize, sinkrho, denshigh, ihigh
      COMMON /gridv/ iline
      COMMON /files/ filein(10000)
      COMMON /origi/ cmx, cmy, cmz

      REAL*4 rhonew(idim)
      INTEGER*1 iphasenew(idim)
      CHARACTER*100 fileident
      INTEGER nums(8),numssink(8),numsrt(8),numsmhd(8)

      INTEGER*8 iuniquelist(2000)
      DIMENSION formtime(2000)
      DIMENSION list(idim)

      CHARACTER*7 filein
      CHARACTER*13 contour2
      CHARACTER*12 contour
      CHARACTER*11 imageout
      CHARACTER*9  conname
      CHARACTER*1 icent, icor, isink, irt, irt2, inewdump, idt, iselect
      CHARACTER*1 ihigh, isingle, iglobal

      DATA pi/3.141592654/
c
c--read options
c
      WRITE (*,*) 'Enter number of input files:'
      READ (*, *) number
c      number =1
      
      WRITE (*,*) 'Enter names of input files:'
      DO k=1, number
         READ(*,90) filein(k)
 90      FORMAT(A7)
      END DO

      WRITE (*, *) 'Give plotting window (only 4 numbers needed)'
      READ (*, *) xmin, xmax, ymin, ymax
 
 110  FORMAT(A1)
c
c--loop over files present
c
      WRITE (*,184)
 184  FORMAT('New dump format (y/n)?')
      READ (*,110) inewdump

      IF (inewdump.NE.'Y' .AND. inewdump.NE.'y') THEN
         WRITE (*,185)
 185     FORMAT('Sink particle code (y/n)?')
         READ (*,110) isink
         
         IF (isink.EQ.'y' .OR. isink.EQ.'Y') THEN
            WRITE (*,187)
 187        FORMAT('Radiative transfer code (y/n)?')
            READ (*,110) irt
            IF (irt.EQ.'y' .OR. irt.EQ.'Y') THEN
 188           FORMAT('With cv (y/n)?')
               READ (*,110) irt2
            ENDIF
         ENDIF
      ENDIF

      denshigh = 1.0E+30
      IF (isink.EQ.'y' .OR. isink.EQ.'Y' .OR. inewdump.EQ.'Y' .OR. 
     &     inewdump.EQ.'y') THEN
         WRITE (*,186)
 186     FORMAT('Enter size and density of sink on image')
         READ (*,*) sinksize,sinkrho

         WRITE (*,189)
 189     FORMAT('Do you want very high density plotted as sink(y/n)?')
         READ (*,*) ihigh

         IF (ihigh.EQ.'y' .OR. ihigh.EQ.'Y') THEN
            WRITE (*,191)
 191        FORMAT('Enter critical density (g cm^-3)')
            READ (*,*) denshigh
         ENDIF
      ENDIF

      WRITE (*,190)
 190  FORMAT('What dumps to process - every one, every 2nd, etc?')
      READ (*,*) ifreq

      WRITE (*,200)
 200  FORMAT('What dump to start from?')
      READ (*,*) istart

      WRITE (*,210)
 210    FORMAT('Do you want to centre on the centre of mass (m) or',
     &     ' on the maximum density (d) or file (f) or not at all?')
      READ (*,110) icent

      iglobal = 'n'
      IF (icent.EQ.'f') THEN
         WRITE (*,216)
 216     FORMAT('Do you want to do only one sink? (y/n)')
         READ (*,110) isingle
         IF (isingle.EQ.'y' .OR. isingle.EQ.'Y') THEN
            WRITE (*,217)
 217        FORMAT('Enter number of sink')
            READ (*,*) isinkdo
         ELSE
            isingle = 'n'
            WRITE (*,218)
 218        FORMAT('Do you want all sinks in one global file? (y/n)')
            READ (*,110) iglobal
            IF (iglobal.EQ.'y' .OR. iglobal.EQ.'Y') THEN
               iglobal = 'y'
               global = 0.
            ENDIF
         ENDIF

         OPEN (11,FILE='fort.88')
         i = 1
 212     READ (11,*,END=214) formtime(i),dummy,ival,iuniquelist(i)
         print *,i,formtime(i),iuniquelist(i)
         i = i+1
         GOTO 212
 214     CLOSE(11)
         nform = i - 1
      ELSE
         nform = 1
      ENDIF

      WRITE (*,*) 'Undo rotating reference frame (y/n/s(pecial)) ?'
      READ (*,110) icor

      IF (icor.EQ.'y' .OR. icor.EQ.'s') THEN
         WRITE (*,*) 'Enter rotation frequency (sec) ?'
         READ (*, *) rotfreq
      ENDIF

      WRITE (*,220)
 220  FORMAT('How many grid points in the plots?')
      READ (*,*) iline
c
c--Determine number of images wide and tall for global images
c
      IF (iglobal.EQ.'y') THEN
         ntall = INT(SQRT(nform/1.6))+1
         nwide = INT(ntall*1.6)
 215     IF (ntall*nwide.LT.nform+5) THEN
            nwide = nwide+1
            GOTO 215
         ENDIF
         ntallpix = ntall*(iline+1)
         nwidepix = nwide*(iline+1)
         print *,'nwidepix, ntallpix: ',nwidepix, ntallpix
         IF (ntallpix.GT.ntallmax .OR. 
     &        nwidepix.GT.nwidemax) THEN
            PRINT *,'ERROR - global image array too small ',
     &           ntallmax,nwidemax,ntallpix,nwidepix
            STOP
         ENDIF
      ENDIF
            
      WRITE (*,230)
 230  FORMAT('Enter angles of viewing (0-360,0-180)?',/,
     &     '    First angle:  rotate in x-y plane, anticlockwise',/,
     &     '    Second angle: rotate from in y-z plane after 1st rot'
     &    ,/,'    Third angle: rotate from in x-z plane after 2nd rot'
     &     ,/,'  (i.e. 0,0,0 is looking down on the x-y plane')
      READ (*,*) angle1, angle2, angle3
      angle1rad = angle1*pi/180.
      angle2rad = angle2*pi/180.
      angle3rad = angle3*pi/180.

      WRITE (*,240)
 240  FORMAT('Do you want (d)ensity or (t)emperature or (r)adiation',
     &     ' temperature or mass (w)eighted radiation temp',
     &     ' or (o)ptical depth?')
      READ (*,110) idt

      WRITE (*,250)
 250  FORMAT('Do you want to select files by time (y/n)?')
      READ (*,110) iselect
      IF (iselect.EQ.'y' .OR. iselect.EQ.'Y') THEN
         inumber = 0
         WRITE (*,260)
 260          FORMAT('Enter time step')
         READ (*,*) timestep
      ENDIF
c
c--Process data files
c
      CALL buildcol

      DO k=1,number
         contour = 'CC' // filein(k)
         imageout = 'I' // filein(k)
         OPEN(15,file=imageout)

         OPEN (UNIT = 11, FILE = filein(k), FORM = 'unformatted')
c     
c--process file dumps
c
         ifile = 1
         zero = 0
         WRITE(contour, 99002) contour, zero, zero, ifile
99002    FORMAT(A9, I1, I1, I1)
c
c--skip files
c
 10      CONTINUE

         CALL rdump(11, ichkl, 0)

         print *,'Units are ',umassi,udisti
         print *,'Time is ',gt
         time = gt
         
         print *,'Writing file fort.22 which links i and iunique'

         IF (icent.EQ.'f') THEN
            DO i = 1, npart
               IF (iphase(i).GT.0) THEN
                  WRITE (22,*) i,iphase(i),iunique(i)
                  print *,'link ',i,iphase(i),iunique(i),xyzmh(1,i),
     &                 xyzmh(2,i),xyzmh(3,i),xyzmh(4,i)
               ENDIF
c--MUST avoid dead particles (iphase=-1) because they might have copies
c     of iunique values associated with them
c
               IF (iphase(i).GE.0) THEN
                  IF (iunique(i).LE.idim) THEN
                     ivalue = iunique(i)
                     list(ivalue) = i

c                     print *,'Setting ',ivalue,list(ivalue),i,iunique(i)
                  ELSE
                     print *,'iunique(i).GT.idim'
                     STOP
                  ENDIF
               ENDIF

            END DO
         ENDIF

         IF (iselect.EQ.'y' .OR. iselect.EQ.'Y') THEN
            inumber = inumber + 1
            IF (inumber.EQ.1) starttime = time - timestep
            fraction = (time-starttime)/timestep 
            IF (fraction.LT.0.999999) THEN
               WRITE (*,*) 'Ignoring ',contour,time,starttime,timestep,
     &              fraction,(time-starttime)/timestep
               GOTO 430
            ELSE
               starttime = time
            ENDIF
         ENDIF
c
c--rotate for rotating reference frame
c
         IF (icor.EQ.'y' .OR. icor.EQ.'s') THEN
            romega = rotfreq*utime
            WRITE(*,*) romega
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,vxyzu,romega,time,icor,npart)
C$OMP& private(r,th,i)
            DO i=1,npart
               r=sqrt(xyzmh(1,i)*xyzmh(1,i)+xyzmh(2,i)*xyzmh(2,i))
               th=ATAN2(xyzmh(2,i),xyzmh(1,i))
               IF (icor.EQ.'s') THEN
                  th=th+8.26553*romega
               ELSE
                  th=th+time*romega
               ENDIF
               xyzmh(1,i)=r*COS(th)
               xyzmh(2,i)=r*SIN(th)

               r=sqrt(vxyzu(1,i)*vxyzu(1,i)+vxyzu(2,i)*vxyzu(2,i))
               th=ATAN2(vxyzu(2,i),vxyzu(1,i))
               IF (icor.EQ.'s') THEN
                  th=th+8.26553*romega
               ELSE
                  th=th+time*romega
               ENDIF
               vxyzu(1,i)=r*COS(th)
               vxyzu(2,i)=r*SIN(th)
            END DO
C$OMP END PARALLEL DO
         ENDIF
c
c--reset coordinates to center of mass
c
         frac = 1.0
         IF (icent.EQ.'m') THEN
            CALL origin (n1, n2, frac)
         ELSEIF (icent.EQ.'d') THEN
            irhomax = 1
            DO i = 1, npart
               IF (rho(i).GT.rho(irhomax)) irhomax = i
            END DO
            WRITE (*,*) 'irhomax ',irhomax, rho(irhomax), 
     &           xyzmh(1,irhomax), xyzmh(2,irhomax), xyzmh(3,irhomax)
            ncen = 0
            xcen = 0.
            ycen = 0.
            zcen = 0.
            DO i = 1, npart
               dx = xyzmh(1,i) - xyzmh(1,irhomax)
               dy = xyzmh(2,i) - xyzmh(2,irhomax)
               dz = xyzmh(3,i) - xyzmh(3,irhomax)
               r2 = dx*dx + dy*dy + dz*dz
               IF (r2.LT.16.0*xyzmh(5,irhomax)*xyzmh(5,irhomax)) THEN
                  ncen = ncen + 1
                  xcen = xcen + xyzmh(1,i)
                  ycen = ycen + xyzmh(2,i)
                  zcen = zcen + xyzmh(3,i)
               ENDIF
            END DO
            xcen = xcen/ncen
            ycen = ycen/ncen
            zcen = zcen/ncen
            WRITE (*,*) 'centre: ',ncen, xcen, ycen, zcen
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,xcen,ycen,zcen,frac,npart)
C$OMP& private(i)
            DO i = 1, npart
               xyzmh(1,i) = xyzmh(1,i) - xcen*frac
               xyzmh(2,i) = xyzmh(2,i) - ycen*frac
               xyzmh(3,i) = xyzmh(3,i) - zcen*frac
            END DO
C$OMP END PARALLEL DO
         ENDIF
c
c--If centering on file, process each file nform times
c
         IF (isingle.EQ.'y' .OR. isingle.EQ.'Y') THEN
            istartsink = isinkdo
            iendsink = isinkdo
         ELSE
            istartsink = 1
            iendsink = nform
         ENDIF
         DO nformloop = istartsink, iendsink

            IF (icent.EQ.'f') THEN
               print *,'time ',time,formtime(nformloop)-0.0106,nformloop

               IF (time.LT.formtime(nformloop)-0.0106) GOTO 777
               print *,'Doing ',nformloop
c
c--Centre on each object
c
               ivalue = iuniquelist(nformloop)
               icentre = list(ivalue)
               IF (icentre.LE.0) THEN
                  WRITE (*,*) 'ERROR - icentre.LE.0'
                  STOP
               ENDIF
               IF(iphase(icentre).LT.0) THEN
                  WRITE (*,*) 'ERROR - iphase(icentre).LE.0 ',
     &                 iphase(icentre),icentre,iuniquelist(nformloop),
     &                 nformloop,istartsink,iendsink
                  STOP
               ENDIF
               xcentre = xyzmh(1,icentre)
               ycentre = xyzmh(2,icentre)
               zcentre = xyzmh(3,icentre)
               DO i = 1, npart
                  xyzmh(1,i) = xyzmh(1,i) - xcentre
                  xyzmh(2,i) = xyzmh(2,i) - ycentre
                  xyzmh(3,i) = xyzmh(3,i) - zcentre
               END DO
            ENDIF
c
c--Rotate by angles
c
            IF (angle1rad.NE.0.0 .OR. 
     &           angle2rad.NE.0.0 .OR. 
     &           angle3rad.NE.0.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,vxyzu,angle1rad,npart)
C$OMP& shared(angle2rad,angle3rad)
C$OMP& private(i,r,th1,th2,th3)
               DO i = 1, npart
                  r = SQRT(xyzmh(1,i)*xyzmh(1,i) +xyzmh(2,i)*xyzmh(2,i))

                  th1 = ATAN2(xyzmh(2,i),xyzmh(1,i))

                  th1 = th1 + angle1rad
                  xyzmh(1,i) = r*COS(th1)
                  xyzmh(2,i) = r*SIN(th1)

                  r = SQRT(xyzmh(2,i)*xyzmh(2,i) +xyzmh(3,i)*xyzmh(3,i))

                  th2 = ATAN2(xyzmh(3,i),xyzmh(2,i))

                  th2 = th2 + angle2rad
                  xyzmh(2,i) = r*COS(th2)
                  xyzmh(3,i) = r*SIN(th2)

                  r = SQRT(xyzmh(1,i)*xyzmh(1,i) +xyzmh(3,i)*xyzmh(3,i))

                  th3 = ATAN2(xyzmh(3,i),xyzmh(1,i))

                  th3 = th3 + angle3rad
                  xyzmh(1,i) = r*COS(th3)
                  xyzmh(3,i) = r*SIN(th3)
                  
               END DO
            ENDIF
c
c--Centre on frame
c
         xcen = (xmax+xmin)/2.0
         ycen = (ymax+ymin)/2.0
         xmaxp = xmax - xcen
         xminp = xmin - xcen
         ymaxp = ymax - ycen
         yminp = ymin - ycen

         DO i = 1, npart
            xyzmh(1,i) = xyzmh(1,i) - xcen
            xyzmh(2,i) = xyzmh(2,i) - ycen
         END DO
         
         WRITE(15,*) time
         WRITE(15,*) contour
c
c--create images
c
         IF (iglobal.EQ.'n') THEN
            IF (icent.EQ.'f') THEN
99101          FORMAT('Z0000000',I1,A4)
99102          FORMAT('Z000000',I2,A4)
99103          FORMAT('Z00000',I3,A4)
99104          FORMAT('Z0000',I4,A4)
99105          FORMAT('Z000',I5,A4)
99106          FORMAT('Z00',I6,A4)
99107          FORMAT('Z0',I7,A4)
99108          FORMAT('Z',I8,A4)
               ival = iuniquelist(nformloop)
               IF (ival.LT.10) THEN
                  WRITE(contour2, 99101) ival,contour(6:9)
               ELSEIF (ival.LT.100) THEN
                  WRITE(contour2, 99102) ival,contour(6:9)
               ELSEIF (ival.LT.1000) THEN
                  WRITE(contour2, 99103) ival,contour(6:9)
               ELSEIF (ival.LT.10000) THEN
                  WRITE(contour2, 99104) ival,contour(6:9)
               ELSEIF (ival.LT.100000) THEN
                  WRITE(contour2, 99105) ival,contour(6:9)
               ELSEIF (ival.LT.1000000) THEN
                  WRITE(contour2, 99106) ival,contour(6:9)
               ELSEIF (ival.LT.10000000) THEN
                  WRITE(contour2, 99107) ival,contour(6:9)
               ELSEIF (ival.LT.100000000) THEN
                  WRITE(contour2, 99108) ival,contour(6:9)
               ELSE
                  WRITE(contour2, 99108)MOD(ival,100000000),contour(6:9)
               ENDIF
               OPEN(16,file=contour2, FORM = 'unformatted')
            ELSE
               OPEN(16,file=contour, FORM = 'unformatted')
            ENDIF
            WRITE(*,*) contour,' ',contour(6:9),' ',contour2
            WRITE(*,*) time
         ENDIF

         CALL grid(xminp, xmaxp, yminp, ymaxp, nformloop, idt, iglobal)

         IF (iglobal.EQ.'n') CLOSE(16)

 777     CONTINUE

c         STOP

         END DO ! nformloop

         IF (iglobal.EQ.'y') THEN
            contour(2:2) = 'G'

            timeout = time
            DO i = 1, nwidepix
               DO j = 1, ntallpix
                  xmaxdensout = MAX(xmaxdensout,global(i,j))
                  IF (global(i,j).NE.0.) xmindensout = 
     &                 MIN(xmindensout,global(i,j))
               END DO
            END DO

            PRINT *,contour,': ',timeout,xmaxdensout,xmindensout

            OPEN(16,file=contour, FORM = 'unformatted')
      WRITE(16,IOSTAT=io) nwidepix,ntallpix,timeout,xmaxdensout,
     &           xmindensout
      WRITE(16,IOSTAT=io) (global(i,1:ntallpix), i=1,nwidepix)
            CLOSE(16)

            global = 0.
         ENDIF

         PRINT *,'Updating filename'

 430     READ(contour, 99008) conname, ifile
99007    FORMAT(A9,I1,I2)
99008    FORMAT(A9,I3)
         ifile = ifile + 1
         IF (ifile.LE.9) THEN
            WRITE(contour, 99002) conname, zero, zero, ifile
         ELSEIF (ifile.LE.99) THEN
            WRITE(contour, 99007) conname, zero, ifile
         ELSE
            WRITE(contour, 99008) conname, ifile
         ENDIF

         PRINT *,'Updated filename'

         IF (inewdump.NE.'Y' .AND. inewdump.NE.'y') GOTO 10

   20    CONTINUE

         CLOSE(11)

         PRINT *,'Closing 15'
         
         CLOSE(15)

      END DO

      STOP
      END

      SUBROUTINE quit(i)
      STOP
      END

      SUBROUTINE endrun
      CALL quit(0)
      END
      
