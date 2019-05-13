      PROGRAM images
c***********************************************************
c                                                          *
c  This programme generates a .pix file from an sphNG      *
c     dump file.  It can be viewed using ssplash as        *
c     ssplash -readpix ascii <dump_file_name>              *
c                                                          *
c     For some reason, ssplash assumes the axes are        *
c     x and z, so the labels and numbers may need to be    *
c     changed.                                             *
c                                                          *
c***********************************************************

      INCLUDE 'idim'

      INTEGER*4 iline, ntallpix, nwidepix

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/units'

      COMMON /gridv/ iline, isec
      COMMON /origi/ cmx, cmy, cmz
      COMMON /files/ fileout, filein(10000)
      COMMON /optio/ what, idump
      COMMON /eost/ temp(idim)

      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      PARAMETER (ndensity = 600)
      PARAMETER (ntemp = 2.7*200)
      REAL pixmap(ndensity,ntemp)
      INTEGER npixmap(ndensity,ntemp)

      CHARACTER*7 filein
      CHARACTER*11 fileout
      CHARACTER*1  icent, icor
c
c--read options
c
      WRITE (*,*) 'Enter number of input files:'
      READ (*, *) number
      
      WRITE (*,*) 'Enter names of input files:'
      DO k=1, number
         READ(*,90) filein(k)
 90      FORMAT(A7)
      END DO

 110  FORMAT(A1)
c
c--loop over files present
c
c      WRITE (*,210)
c 210  FORMAT('Do you want to centre on the centre of mass (m) or',
c     &     ' on the maximum density (d) or not at all (n)?')
c      READ (*,110) icent

      WRITE (*,*) 'Undo rotating reference frame (y/n) ?'
      READ (*,110) icor

      IF (icor.EQ.'y') THEN
         WRITE (*,*) 'Enter rotation frequency (sec) ?'
         READ (*, *) rotfreq
      ENDIF
c
c--Set units
c
      CALL unit

      WRITE (*,*) 'Density unit ',udens
c
c--Process data files
c
      DO k=1,number

         fileout = filein(k) // '.pix'
         OPEN(21,file=fileout)

         OPEN (UNIT = 11, FILE = filein(k), FORM = 'unformatted',
     &        RECL=imaxrec)
c
c--process file dumps
c
         ifile = 1
         izero = 0

 10      CONTINUE

         CALL rdump(11, ichkl, 0)

         print *,'Units are ',umassi,udisti
         IF (umass.NE.umassi) THEN
            WRITE (*,*) 'Mass units inconsistent ',umass, umassi
         ENDIF
         IF (udist.NE.udisti) THEN
            WRITE (*,*) 'Distance units inconsistent ',udens, udensi
         ENDIF
            
         print *,'Time is ',gt
         time = gt

         WRITE (*,*) udist, umass, utime
         uergg = udist**2/utime**2
         uergcc = (umass)/((udist)*(utime)**2)
         uradconst = 7.5646e-15/uergcc
         WRITE (*,*) npart, n1, n2
         WRITE (*,*) time,gamma,rhozero,RK2

         WRITE(*,*) npart, rhozero, nptmass, rho(3), iphase(1)
         WRITE(*,*) 'RK2=',RK2,' gamma=',gamma,' rhozero=',rhozero
         WRITE(*,*) 'Units=',udist, umass, utime

c
c--rotate for rotating reference frame
c
         IF (icor.EQ.'y') THEN
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

         tcomp = sqrt((3*pi)/(32*rhozero))
         tff = tcomp * utime
         timeff = time/tcomp
c
c--Make pixel map
c
         pixmap = 0.
         npixmap = 0

         ikount = 0
         totmass = 0.
         distance = 1.0e+20
         denscen = 0.
         densmean = 0.
         nactive = 0
         DO 500 i = 1, npart
c         IF(MOD(i,100).EQ.0) write (*,*) npart,i

            IF (iphase(i).NE.-1) THEN
               nactive = nactive + 1            

               IF (ekcle(3,i).NE.0.0) THEN
                  temperature = vxyzu(4,i)/ekcle(3,i)
                  temprad = (ekcle(1,i)*rho(i)/uradconst)**0.25
               ELSEIF (gamma.EQ.1.0) THEN
                  temperature = 2.0*uergg/Rg*2.0/3.0*vxyzu(4,i)
                  temprad = 0.0
               ELSE
                  temperature = 2.0*uergg/Rg*(gamma-1.0)*vxyzu(4,i)
                  temprad = 0.0
               ENDIF

c               WRITE(21,89771) rho(i),temperature,
c     &              chemistry(1,i),
c     &              chemistry(3,i)

               idensity = LOG10(rho(i)*udens/4.0E-24*2.0)*100 -100
c               idensity = LOG10(1000.)*100 -100
               IF (idensity.LT.1) THEN
                  idensity = 1
               ELSEIF (idensity.GT.ndensity) THEN
                  idensity = ndensity
               ENDIF
               itemp = LOG10(temperature)*200
c               itemp = LOG10(50.)*200
               IF (itemp.LT.1) THEN
                  itemp = 1
               ELSEIF (itemp.GT.ntemp) THEN
                  itemp = ntemp
               ENDIF
               pixmap(idensity, itemp) = pixmap(idensity, itemp) +
     &              chemistry(3,i)
c               write (*,*) idensity,itemp,pixmap(idensity,itemp),
c     &              temperature
               npixmap(idensity, itemp) = npixmap(idensity, itemp) + 1

            ENDIF

89005       FORMAT(5(1PE12.3,1X))
89771       FORMAT(17(1PE16.9,1X),I10)

 500     CONTINUE

         WRITE (21,88001) ndensity, ntemp
88001    FORMAT ('# ',I5,I5)

         DO j = 1, ntemp
            DO i = 1, ndensity
               IF (npixmap(i,j).GT.0) THEN
                  pixmap(i,j) = LOG10(pixmap(i,j)/npixmap(i,j) * 1.4E-4)
               ELSE
                  pixmap(i,j) = -10.
               ENDIF
c            IF (i.GT.100 .AND. i.LT.120 .AND. j.EQ.400) pixmap(i,j) = 
c     &           -6.0
c            IF (i.EQ.100 .AND. j.EQ.100) pixmap(i,j) = -6.0
            END DO
            WRITE (21,"(256(es10.3,1x))") (pixmap(i,j),i=1,ndensity)
         END DO

         CLOSE(21)
      END DO

      STOP
      END

      SUBROUTINE quit(i)
      STOP
      END

      SUBROUTINE endrun
      CALL quit(0)
      END
