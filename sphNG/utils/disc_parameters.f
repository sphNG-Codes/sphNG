      PROGRAM images
c***********************************************************
c                                                          *
c  This program generates 2D binary image files of column  *
c     density and density-weighted temperature, but does   *
c     so from two slightly different viewing positions so  *
c     as to allow 3D movies to be generated.               *
c     It uses rdump.F from the main SPH code to read the   *
c     sphNG dump files.  So the appropriate parameters     *
c     need to be set in the sphNG files for the dump files *
c     that are being read.  The value of 'idim' in 'idim'  *
c     needs to be large enough.  If the file contains      *
c     radiative transfer or MHD information, then          *
c     iradtrans or imhd need to be set too.                *
c                                                          *
c     Compile like: make mpi=no openmp=yes movie3DRG       *
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

      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      COMMON /aux  / hmin
      COMMON /gridv/ ilinex,iliney
      COMMON /origi/ cmx, cmy, cmz
      COMMON /files/ filein(10000)
      COMMON /phase2/ sinksize, sinkrho, denshigh, ihigh
      COMMON /numfile/ nfile

      INTEGER*8 iunique_sink(1000)
      DIMENSION ilocalsink(1000),indx(idim)
      DIMENSION discmass(1000),angmom(3,1000),discradius(4,1000)
      DIMENSION distance2(idim)
      DIMENSION listdisc(1000000)

      CHARACTER*7 filein
      CHARACTER*21 contour1, contour2
      CHARACTER*11 imageout
      CHARACTER*10  conname
      CHARACTER*1 icent,icor,imove,isecrot
      CHARACTER*1 irepeat, ihigh, idt

      CHARACTER*4 ivalue
      CHARACTER*9 iuvalue

      DATA pi/3.141592654/
c
c--read options
c
      DO i = 1, idim
         iorig(i) = i
         isort(i) = i
      END DO
      numberstart = 1

      WRITE (*,*) 'Enter number of input files:'
      READ (*, *) numberfiles
      
      WRITE (*,*) 'Enter names of input files:'
      DO k=1, numberfiles
         READ(*,90) filein(k)
 90      FORMAT(A7)
      END DO

      GOTO 444

      WRITE (*, *) 'Give location of observer from centre'
      WRITE (*, *) '     (distance, x, y, half-separation of eyes)'
      WRITE (*, *) '     (If distance between eyes is zero will'
      WRITE (*, *) '     only do one view point, but both density'
      WRITE (*, *) '     and temperature images)'
      READ (*, *) zobs, xobs, yobs, xeye
 
      WRITE (*, *) 'Give plotting window (opening angle)'
      READ (*, *) wangle
 
 110  FORMAT(A1)
c
c--loop over files present
c
      denshigh = 1.0E+30
      WRITE (*,186)
 186  FORMAT('Enter size (code units) and density of sink on image')
      READ (*,*) sinksizestart,sinkrho
      sinksize = sinksizestart
         
      WRITE (*,189)
 189  FORMAT('Do you want very high density plotted as sink(y/n)?')
      READ (*,*) ihigh

      IF (ihigh.EQ.'y' .OR. ihigh.EQ.'Y') THEN
         WRITE (*,191)
 191     FORMAT('Enter critical density (g cm^-3)')
         READ (*,*) denshigh
      ENDIF

      WRITE (*,210)
 210    FORMAT('Do you want to centre on the centre of mass (m) or',
     &     ' on the maximum density (d) or not at all?')
      READ (*,110) icent

      WRITE (*,*) 'Undo rotating reference frame (y/n/s(pecial)) ?'
      READ (*,110) icor

      IF (icor.EQ.'y' .OR. icor.EQ.'s') THEN
         WRITE (*,*) 'Enter rotation frequency (sec) ?'
         READ (*, *) rotfreq
      ENDIF

      WRITE (*,220)
 220  FORMAT('How many grid points in the plots (horizontal:1.6)?')
      READ (*,*) ilinex
      iliney=ilinex/1.6

      WRITE (*,230)
 230    FORMAT('Enter start angles of viewing (0-360,0-180)?',/,
     &     '    First angle:  rotate in x-y plane, anticlockwise',/,
     &     '    Second angle: rotate from in y-z plane after 1st rot'
     &     ,/,'    Third angle: rotate from in x-z plane after 2nd rot'
     &     ,/,'  (i.e. 0,0,0 is looking down on the x-y plane')
      READ (*,*) angle1, angle2, angle3
      angle1rads = angle1*pi/180.
      angle2rads = angle2*pi/180.
      angle3rads = angle3*pi/180.

      WRITE (*,233)
 233  FORMAT('Do you want to repeat frames?')
      READ (*,110) irepeat

      WRITE (*,234)
 234  FORMAT('Enter minimum h?')
      READ (*,*) hminstart
      hmin = hminstart

      IF (xeye.NE.0.) THEN
         WRITE (*,286)
 286     FORMAT('Density or temperature (d/t)?')
         READ (*,110) idt
      ENDIF

      IF (irepeat.EQ.'y') THEN
         WRITE (*,235)
 235     FORMAT('Enter end angles of viewing')
         READ (*,*) angle1, angle2, angle3
         angle1rade = angle1*pi/180.
         angle2rade = angle2*pi/180.
         angle3rade = angle3*pi/180.
         
         WRITE (*,237)
 237     FORMAT('Second rotation?')
         READ (*,110) isecrot
         zcentre = 0.
         IF (isecrot.EQ.'y') THEN
            WRITE (*,236)
 236        FORMAT('Enter 6 angles for second rotation')
            READ (*,*) angle1, angle2, angle3, angle4, angle5, angle6
            angle1rad2s = angle1*pi/180.
            angle2rad2s = angle2*pi/180.
            angle3rad2s = angle3*pi/180.
            angle1rad2e = angle4*pi/180.
            angle2rad2e = angle5*pi/180.
            angle3rad2e = angle6*pi/180.
            WRITE (*,238)
 238        FORMAT('Enter z for centring')
            READ (*,*) zcentre
         ENDIF
         
         WRITE (*,240)
 240     FORMAT('Enter number of output files?')
         READ (*,*) numberfiles
         
         WRITE (*,245)
 245     FORMAT('Enter number of output file to start from?')
         READ (*,*) numberstart
         
         WRITE (*,250)
 250     FORMAT('Enter ratio for subsequent frame?',/,
     &        '    >1.0 for zoom out',/,
     &        '    <1.0 for zoom in ')
         READ (*,*) ratio
         
         WRITE (*,251)
 251     FORMAT('Enter end wangle?')
         READ (*,*) wangleend
         
         WRITE (*,260)
 260     FORMAT('Do you wish to move to center the of mass(m)',
     &        ' or density maximum(d) or position (p)?')
         READ (*,110) imove
         IF (imove.EQ.'p') THEN
            WRITE (*,270)
 270        FORMAT('Enter coordinates of position')
            READ (*,*) xpos, ypos
         ENDIF
         WRITE (*,274)
 274     FORMAT('Enter minimum h to end with?')
         READ (*,*) hminend

      ELSE
         angle1rade = angle1rads
         angle2rade = angle2rads
         angle3rade = angle3rads
      ENDIF
c     
c--Process data files
c
c      CALL buildcol

c      ifile = 1
 444  ifile = numberstart
      izero = 0

      DO k=numberstart,numberfiles

         IF (irepeat.EQ.'y') THEN
            nfile = 1
         ELSE
            nfile=k
            ifile=1
         ENDIF

         IF (xeye.NE.0.) THEN
            contour1 = 'CCL' // filein(nfile)
            contour2 = 'CCR' // filein(nfile)
         ELSE
            contour1 = 'CCD' // filein(nfile)
            contour2 = 'CCT' // filein(nfile)
         ENDIF
         imageout = 'I' // filein(nfile)
         OPEN(15,file=imageout)

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
c
c--skip files
c
 10      CONTINUE

         CALL rdump(11, ichkl, 0)

         CLOSE (11)

         print *,'Units are ',umassi,udisti

         IF (k.EQ.1) THEN
            IF (nptmass.GT.1) THEN
               PRINT *,'ERROR - need to start with <=1 sinks'
               STOP
            ENDIF
            nunique_sink = nptmass
            IF (nunique_sink.EQ.1) THEN
               iunique_sink(1) = iunique(listpm(1))
            ENDIF
         ELSE
            IF (nptmass.GT.nunique_sink) THEN
               DO iii = 1, nptmass
c
c--Check not already known
c
                  DO jjj = 1, nunique_sink
                     IF (iunique_sink(jjj).EQ.iunique(listpm(isink)))
     &                    GOTO 678
                  END DO
c
c--Add unknown sinks to list
c
                  nunique_sink = nunique_sink + 1
                  iunique_sink(nunique_sink) = iunique(listpm(iii))
 678              CONTINUE
               END DO
            ENDIF
         ENDIF
         IF (nptmass.NE.nunique_sink) THEN
            PRINT *,'ERROR - nptmass.NE.nunique_sink ',nptmass,
     &           nunique_sink
            STOP
         ENDIF

         PRINT *,'Number of sinks ',nptmass,nunique_sink
         PRINT *,'Names of sinks ',(iunique_sink(iii),
     &        iii=1,nunique_sink)
c
c--Find local index
c
         DO iii = 1, nunique_sink
            DO jjj = 1, nptmass
               IF (iunique(listpm(jjj)).EQ.iunique_sink(iii)) THEN
                  ilocalsink(iii) = listpm(jjj)
                  GOTO 679
               ENDIF
            END DO
            PRINT *,'ERROR - not found ',iii,iunique_sink(iii)
            STOP
 679        CONTINUE
         END DO
c
c--Process sinks
c
         DO isink = 1, nunique_sink
            iptcur = ilocalsink(isink)
            xsink = xyzmh(1,iptcur)
            ysink = xyzmh(2,iptcur)
            zsink = xyzmh(3,iptcur)
            sinkmass = xyzmh(4,iptcur)

            print *,'Sink location ',xsink,ysink,zsink,sinkmass

            discmass(isink) = 0.
            DO l = 1, 3
               angmom(l,isink) = 0.
            END DO
            ndisc = 0

            DO i = 1, npart
               dx = xyzmh(1,i)-xsink
               dy = xyzmh(2,i)-ysink
               dz = xyzmh(3,i)-zsink
               distance2(i) = dx**2 + dy**2 + dz**2
            END DO

            print *,'sorting'

            CALL indexx2(npart,distance2,indx)

            print *,'sorted ',indx(1),indx(2),distance2(indx(1)),
     &           distance2(indx(2))

            DO j = 1, npart
               i = indx(j)

               totalmass = sinkmass + discmass(isink)

               IF (iphase(i).EQ.0 .AND. i.NE.iptcur) THEN
                  dx = xyzmh(1,i)-xsink
                  dy = xyzmh(2,i)-ysink
                  dz = xyzmh(3,i)-zsink
                  dist = SQRT(dx**2 + dy**2 + dz**2)
c
c--Make distance cut - 1000 AU or 1/2 distance to nearest sink
c
                  IF (dist*udisti.LT.1.496E+13*1000.) THEN
                     dvx = vxyzu(1,i)-vxyzu(1,iptcur)
                     dvy = vxyzu(2,i)-vxyzu(2,iptcur)
                     dvz = vxyzu(3,i)-vxyzu(3,iptcur)

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
     &                       (totalmass**2)
                        IF (tempnum.LT.0.0) THEN
                           WRITE(iprint,*) 'ERROR - tempnum'
                           WRITE(iprint,*) vpotdif,specangmom2,pmassi
                           WRITE(iprint,*) vkep2,divvr2,r2,vtan2
                           STOP
                        ELSE
                           eccentricity = SQRT(tempnum)
                           radapastron = - totalmass*
     &                          (1.0+eccentricity)/(2.0*vpotdif)
                        ENDIF
                     ELSE
                        eccentricity = 1.0E+10
                        radapastron = 1.0E+10
                     ENDIF
c
c--Add to disc if passes tests
c
c                     print *,i,dist,eccentricity,vkep2,divvr2,
c     &                    semimajoraxis

                     IF (eccentricity.LT.0.5 .AND. 
     &                    radapastron*udisti.LT.1.496E+13*1000.0) THEN
                        discmass(isink) = discmass(isink) + xyzmh(4,i)
                        angmom(1,isink) = angmom(1,isink) +
     &                       xyzmh(4,i)*(dvy*dz - dvz*dy)
                        angmom(2,isink) = angmom(2,isink) +
     &                       xyzmh(4,i)*(dvz*dx - dvx*dz)
                        angmom(3,isink) = angmom(3,isink) +
     &                       xyzmh(4,i)*(dvx*dy - dvy*dx)
                        ndisc = ndisc + 1
                        listdisc(ndisc) = i

                        WRITE (44,'(3(1PE12.5,1X))') xyzmh(1,i),
     &                       xyzmh(2,i),xyzmh(3,i)
                     ENDIF
                  ENDIF
               ENDIF
            END DO
            discradius(1,isink) = 
     &           SQRT(distance2(listdisc(INT(0.5*ndisc))))
            discradius(2,isink) = 
     &           SQRT(distance2(listdisc(INT(0.9*ndisc))))
            discradius(3,isink) = 
     &           SQRT(distance2(listdisc(INT(0.95*ndisc))))
            discradius(4,isink) = 
     &           SQRT(distance2(listdisc(INT(0.99*ndisc))))
            PRINT *,'Discmass ',isink,discmass(isink),ndisc
            PRINT *,'Disc radius ',discradius(1,isink),
     &           discradius(2,isink),discradius(3,isink),
     &           discradius(4,isink)
            PRINT *,'Angular momomentum ',angmom(1,isink),
     &           angmom(2,isink),angmom(3,isink)
         END DO

         DO i = 1, nunique_sink
            IF (i.GT.99) THEN
               WRITE (ivalue,"(I3)") i
            ELSEIF (i.GT.9) THEN
               WRITE (ivalue,"('0',I2)") i
            ELSE
               WRITE (ivalue,"('00',I1)") i
            ENDIF
            IF (iunique_sink(i).GT.9999999) THEN
               WRITE (iuvalue,"(I8)") iunique_sink(i)
            ELSEIF (iunique_sink(i).GT.999999) THEN
               WRITE (iuvalue,"('0',I7)") iunique_sink(i)
            ELSEIF (iunique_sink(i).GT.99999) THEN
               WRITE (iuvalue,"('00',I6)") iunique_sink(i)
            ELSEIF (iunique_sink(i).GT.9999) THEN
               WRITE (iuvalue,"('000',I5)") iunique_sink(i)
            ELSEIF (iunique_sink(i).GT.999) THEN
               WRITE (iuvalue,"('0000',I4)") iunique_sink(i)
            ELSEIF (iunique_sink(i).GT.99) THEN
               WRITE (iuvalue,"('00000',I3)") iunique_sink(i)
            ELSEIF (iunique_sink(i).GT.9) THEN
               WRITE (iuvalue,"('000000',I2)") iunique_sink(i)
            ELSE
               WRITE (iuvalue,"('0000000',I1)") iunique_sink(i)
            ENDIF
            contour1 = 'Disc_' // ivalue(1:3) // '_' // iuvalue(1:8) 
     &           // '.dat'

            print *,contour1

            OPEN (16,file=contour1,ACCESS='append')
            WRITE (16,"(9(1PE12.5,1x))") gt,discmass(i),
     &           (discradius(j,i),j=1,4),(angmom(j,i),j=1,3)
            CLOSE (16)

         END DO

c         WRITE (95,'(1465(1PE12.5,1X))') gt,(discmass(i),i=1,183),
c     &        ((discradius(j,i),j=1,4),i=1,183)
c     &        ((angmom(j,i),j=1,3),i=1,183)

         GOTO 20



         OPEN(16,file=contour1, FORM = 'unformatted')
         OPEN(17,file=contour2, FORM = 'unformatted')
         WRITE(*,99005) contour1,contour2
         WRITE(*,*) time
99005    FORMAT(A13,1X,A13)

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
         IF (irepeat.EQ.'y') frac = REAL(k)/numberfiles
         write (*,*) 'frac ',frac

         IF (icent.EQ.'m' .OR. imove.EQ.'m') THEN
c            CALL origin (n1, n2, frac)
         ENDIF
         IF (icent.EQ.'d' .OR. imove.EQ.'d') THEN
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
c--Rotate by angles
c
         IF (angle1rads.NE.0.0 .OR. angle1rade.NE.0.0 .OR.
     &        angle2rads.NE.0.0 .OR. angle2rade.NE.0.0 .OR.
     &        angle3rads.NE.0.0 .OR. angle3rade.NE.0.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,vxyzu,angle1rads,angle1rade,npart)
C$OMP& shared(angle2rads,angle2rade,angle3rads,angle3rade,frac)
C$OMP& private(i,r,th1,th2,th3)
            DO i = 1, npart
            r = SQRT(xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i))

            th1 = ATAN2(xyzmh(2,i),xyzmh(1,i))

            th1 = th1 + angle1rads+frac*(angle1rade-angle1rads)
            xyzmh(1,i) = r*COS(th1)
            xyzmh(2,i) = r*SIN(th1)

c            r = SQRT(vx(i)*vx(i) + vy(i)*vy(i))

c            th1 = ATAN2(vy(i),vx(i))

c            th1 = th1 + angle1rads+frac*(angle1rade-angle1rads)
c            vx(i) = r*COS(th1)
c            vy(i) = r*SIN(th1)

            r = SQRT(xyzmh(2,i)*xyzmh(2,i) + xyzmh(3,i)*xyzmh(3,i))

            th2 = ATAN2(xyzmh(3,i),xyzmh(2,i))

            th2 = th2 + angle2rads+frac*(angle2rade-angle2rads)
            xyzmh(2,i) = r*COS(th2)
            xyzmh(3,i) = r*SIN(th2)

c            r = SQRT(vy(i)*vy(i) + vz(i)*vz(i))

c            th2 = ATAN2(vz(i),vy(i))

c            th2 = th2 + angle2rads+frac*(angle2rade-angle2rads)
c            vy(i) = r*COS(th2)
c            vz(i) = r*SIN(th2)

            r = SQRT(xyzmh(1,i)*xyzmh(1,i) + xyzmh(3,i)*xyzmh(3,i))

            th3 = ATAN2(xyzmh(3,i),xyzmh(1,i))

            th3 = th3 + angle3rads+frac*(angle3rade-angle3rads)
            xyzmh(1,i) = r*COS(th3)
            xyzmh(3,i) = r*SIN(th3)

c            r = SQRT(vx(i)*vx(i) + vz(i)*vz(i))

c            th3 = ATAN2(vz(i),vx(i))

c            th3 = th3 + angle3rads+frac*(angle3rade-angle3rads)
c            vx(i) = r*COS(th3)
c            vz(i) = r*SIN(th3)
            END DO
C$OMP END PARALLEL DO
         ENDIF
c
c--Move to specified position
c
         IF (imove.EQ.'p') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,xpos,ypos,frac,npart)
C$OMP& private(i)
            DO i = 1, npart
               xyzmh(1,i) = xyzmh(1,i) - xpos*SQRT(frac)
               xyzmh(2,i) = xyzmh(2,i) - ypos*SQRT(frac)
            END DO
C$OMP END PARALLEL DO
         ENDIF
c
c--Rotate by angles
c
         IF (isecrot.EQ.'y') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,vxyzu,angle1rad2s,angle1rad2e,npart)
C$OMP& shared(angle2rad2s,angle2rad2e,angle3rad2s,angle3rad2e,frac)
C$OMP& private(i,r,th1,th2,th3)
            DO i = 1, npart
            r = SQRT(xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i))

            th1 = ATAN2(xyzmh(2,i),xyzmh(1,i))

            th1 = th1 + angle1rad2s+frac*(angle1rad2e-angle1rad2s)
            xyzmh(1,i) = r*COS(th1)
            xyzmh(2,i) = r*SIN(th1)

            r = SQRT(vxyzu(1,i)*vxyzu(1,i) + vxyzu(2,i)*vxyzu(2,i))

            th1 = ATAN2(vxyzu(2,i),vxyzu(1,i))

            th1 = th1 + angle1rad2s+frac*(angle1rad2e-angle1rad2s)
            vxyzu(1,i) = r*COS(th1)
            vxyzu(2,i) = r*SIN(th1)

            r = SQRT(xyzmh(2,i)*xyzmh(2,i) + xyzmh(3,i)*xyzmh(3,i))

            th2 = ATAN2(xyzmh(3,i),xyzmh(2,i))

            th2 = th2 + angle2rad2s+frac*(angle2rad2e-angle2rad2s)
            xyzmh(2,i) = r*COS(th2)
            xyzmh(3,i) = r*SIN(th2)

            r = SQRT(vxyzu(2,i)*vxyzu(2,i) + vxyzu(3,i)*vxyzu(3,i))

            th2 = ATAN2(vxyzu(3,i),vxyzu(2,i))

            th2 = th2 + angle2rad2s+frac*(angle2rad2e-angle2rad2s)
            vxyzu(2,i) = r*COS(th2)
            vxyzu(3,i) = r*SIN(th2)

            r = SQRT(xyzmh(1,i)*xyzmh(1,i) + xyzmh(3,i)*xyzmh(3,i))

            th3 = ATAN2(xyzmh(3,i),xyzmh(1,i))

            th3 = th3 + angle3rad2s+frac*(angle3rad2e-angle3rad2s)
            xyzmh(1,i) = r*COS(th3)
            xyzmh(3,i) = r*SIN(th3)

            r = SQRT(vxyzu(1,i)*vxyzu(1,i) + vxyzu(3,i)*vxyzu(3,i))

            th3 = ATAN2(vxyzu(3,i),vxyzu(1,i))

            th3 = th3 + angle3rad2s+frac*(angle3rad2e-angle3rad2s)
            vxyzu(1,i) = r*COS(th3)
            vxyzu(3,i) = r*SIN(th3)
            END DO
C$OMP END PARALLEL DO
         ENDIF
c
c--Centre on frame
c
         IF (xobs.NE.0.0 .OR. yobs.NE.0.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(xyzmh,xobs,yobs,npart)
C$OMP& private(i)
            DO i = 1, npart
               xyzmh(1,i) = xyzmh(1,i) - xobs
               xyzmh(2,i) = xyzmh(2,i) - yobs
            END DO
C$OMP END PARALLEL DO
         ENDIF
         
         WRITE(15,*) time
         WRITE(15,*) contour1, contour2

         DO i = 1, npart
            IF (iphase(i).EQ.11) THEN


               rxy = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)

               GOTO 332 

               IF (rxy.GT.1.99 .AND. rxy.LT.2.01 .AND. 
     &              xyzmh(3,i).GT.0.15) THEN
                  print *,i,rxy,xyzmh(3,i),vxyzu(3,i),iunique(i)
               ENDIF

c               GOTO 333

c 332           IF (i.EQ.470405) THEN
 332           IF (iunique(i).EQ.2900916) THEN
                  DO j = 1, npart
                     IF (iphase(j).EQ.0) THEN
                        dist2 = 1.e+10
                        r2 = (xyzmh(1,i)-xyzmh(1,j))**2 + 
     &                       (xyzmh(2,i)-xyzmh(2,j))**2 +
     &                       (xyzmh(3,i)-xyzmh(3,j))**2
                        IF (r2.LT.dist2) THEN
                           dist2 = r2
                           rhogas = rho(j)
                        ENDIF
                     ENDIF
                  END DO
                  print *,'Track ',gt,rxy,xyzmh(3,i),vxyzu(3,i)
                  WRITE (92,3400) gt,rxy,xyzmh(3,i),vxyzu(3,i),rhogas,
     &                 rho(i)
 3400             FORMAT(6(1PE12.5,1X))
 333              CONTINUE
               ENDIF
            ENDIF
         END DO

         GOTO 20
c
c--create images
c
         zobsdo = zobs
         xeyedo = xeye
         wangledo = wangle
         IF (irepeat.EQ.'y') THEN
            zobsdo = zobs*ratio**k
            xeyedo = xeye*ratio**k
            sinksize = sinksizestart*ratio**k
            wangledo = wangle + frac*(wangleend - wangle)
            hmin = hminstart + frac*(hminend - hminstart)
         ENDIF
         WRITE (*,*) 'Doing ',zobsdo,wangledo,xeyedo,ratio,k
         IF (xeyedo.EQ.0.0) THEN
c            CALL grid (zobsdo,wangledo,16,17,xeyedo,idt)
         ELSE
c            CALL grid (zobsdo,wangledo,16,17,xeyedo,idt)
c            CALL grid (zobsdo,wangledo,16,17,-xeyedo,idt)
         ENDIF

         CLOSE(16)
         CLOSE(17)
         READ(contour1, 99008) conname, ifile
99007    FORMAT(A10,I1,I2)
99008    FORMAT(A10,I3)
         ifile = ifile + 1
c         IF (ifile.LE.9) THEN
c            WRITE(contour1, 99002) conname, zero, zero, ifile
c         ELSEIF (ifile.LE.99) THEN
c            WRITE(contour1, 99007) conname, zero, ifile
c         ELSE
c            WRITE(contour1, 99008) conname, ifile
c         ENDIF

c         READ(contour2, 99008) conname, ifile
c         ifile = ifile + 1
c         IF (ifile.LE.9) THEN
c            WRITE(contour2, 99002) conname, zero, zero, ifile
c         ELSEIF (ifile.LE.99) THEN
c            WRITE(contour2, 99007) conname, zero, ifile
c         ELSE
c            WRITE(contour2, 99008) conname, ifile
c         ENDIF

   20    CONTINUE
         
         CLOSE(15)

      END DO

      STOP
      END

      SUBROUTINE quit
      STOP
      END


      SUBROUTINE endrun
      CALL quit(0)
      END
