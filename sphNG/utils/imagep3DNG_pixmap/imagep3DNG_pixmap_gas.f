      PROGRAM images
c***********************************************************
c                                                          *
c  this program generates images suitable for imagetool    *
c                                                          *
c***********************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/units'
      
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      PARAMETER (ndensity = 700)
      PARAMETER (ntemp = 2.7*200)
      REAL pixmap(ndensity,ntemp)
      INTEGER npixmap(ndensity,ntemp)

      CHARACTER*100 fileident
      INTEGER*8 number8
      DIMENSION nums(8),numssink(8),numsRT(8),index(idim),list(idim),
     &     numsMHD(8)

      INTEGER*2 nneigh(idim)

      CHARACTER*7 filein
      DIMENSION filein(1000)
      CHARACTER*11 inname, imageout
      CHARACTER*11 fileout
      CHARACTER*9  name
      CHARACTER*1  ien, icont, icent, isink, itorq, ivarbnd

      DATA pi/3.141592654/
c
c--read options
c
      WRITE (*,*) 'Enter number of input files:'
      READ (*, *) numberfiles

      WRITE (*,*) 'Enter names of input files:'
      DO k=1, numberfiles
         READ(*,90) filein(k)
 90      FORMAT(A7)
      END DO

 110  FORMAT(A1)

      WRITE (*,210)
 210  FORMAT('Do you want to centre on the centre of mass (m) or',
     &     ' on the maximum density (d) or ',
     &     ' on a position (p) or ',
     &     ' on the most massive or xth sink (s,x) or not at all?')
      READ (*,110) icent
      IF (icent.EQ.'x') READ (*,*) ixthsink
      IF (icent.EQ.'p') THEN
         READ (*,*) xpos, ypos, zpos
      ENDIF

      WRITE (*,*) 'Undo rotating reference frame (0,1,2,3) ?'
      WRITE (*,*) '    (3=1 but change velocities to rotating frame)'
      READ (*,*) icor

      IF (icor.NE.0) THEN
         WRITE (*,*) 'Enter rotation frequency (sec) ?'
         READ (*, *) rotfreq
      ENDIF

      DO k=1,numberfiles
c
c--Process data files
c
         fileout = filein(k) // '.pix'
         imageout = 'IP' // filein(k)
         maxrec = 300*idim
         OPEN(15,file=imageout)

         OPEN (UNIT = 11, FILE = filein(k), FORM = 'unformatted')
c
c--process file dumps
c
 10      OPEN(21,file=fileout)
         WRITE(*,99004) fileout
99004    FORMAT(A11)

         CALL rdump(11, ichkl, 0)

         print *,'Units are ',umassi,udisti
         print *,'Time is ',gt
         time = gt

         CALL unit
         
         WRITE (*,*) umass, udist, utime
         uergg = udist**2/utime**2
         uergcc = (umass)/((udist)*(utime)**2)
         uradconst = 7.5646e-15/uergcc
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

         tcomp = sqrt((3*pi)/(32*rhozero))
         tff = tcomp * utime
         timeff = time/tcomp
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
c--Make pixel map
c
         DO i = 1, ndensity
            DO j = 1, ntemp
               pixmap(i,j) = 0.
            END DO
         END DO
         DO i = 1, ndensity
            DO j = 1, ntemp
               npixmap(i,j) = 0
            END DO
         END DO

         ikount = 0
         totmass = 0.
         distance = 1.0e+20
         denscen = 0.
         densmean = 0.
         nactive = 0
         DO i=1, npart
            IF(MOD(i,10000).EQ.0) write (*,*) i

            IF (iphase(i).NE.-1) THEN
               nactive = nactive + 1

               IF (ekcle(3,i).NE.0.0) THEN
                  temperature = vxyzu(4,i)/ekcle(3,i)
                  temprad = (ekcle(1,i)*rho(i)/uradconst)**0.25
               ELSEIF (gamma.EQ.1.0) THEN
                  temperature = 2.0*uergg/8.314E+7*2.0/3.0*vxyzu(4,i)
                  temprad = 0.0
               ELSE
                  temperature =2.0*uergg/8.314E+7*(gamma-1.0)*vxyzu(4,i)
                  temprad = 0.0
               ENDIF
               tmax = MAX(tmax,temperature)
               IF (temperature.GT.5000.0) THEN
                  dmass = dmass + xyzmh(4,i)
               ELSEIF (temperature.GT.4000.0) THEN
                  radiusstar = r1
               ENDIF
               pressurei = 8.314E+7*rho(i)*temperature*1.61/uergg

c               WRITE(21,89771) rho(i),temperature,dust_tk(1,i),
c     &              chemistry(1,i),
c     &              chemistry(3,i)
89771          FORMAT(17(1PE16.9,1X),I10)

               idensity = LOG10(rho(i)*6.7746E-20/4.0E-24*2.0)*50
c               idensity = LOG10(rho(i)*6.7746E-20/4.0E-24*2.0)*100-100
c               idensity = LOG10(1000.)*100 -100
               IF (idensity.LT.1) THEN
                  idensity = 1
               ELSEIF (idensity.GT.ndensity) THEN
                  idensity = ndensity
                  GOTO 99
               ENDIF
               itemp = LOG10(temperature)*200
c               itemp = LOG10(50.)*200
               IF (itemp.LT.1) THEN
                  itemp = 1
                  GOTO 99
               ELSEIF (itemp.GT.ntemp) THEN
                  itemp = ntemp
                  GOTO 99
               ENDIF
               pixmap(idensity, itemp) = pixmap(idensity, itemp) +
     &              chemistry(3,i)
c               write (*,*) idensity,itemp,pixmap(idensity,itemp),
c     &              temperature
               npixmap(idensity, itemp) = npixmap(idensity, itemp) + 1

            ENDIF

 99         CONTINUE
            
         END DO
         
         WRITE(15,99080) time, timeff, fileout
99080    FORMAT(1PE12.5,1X,1PE12.5,1X,A11)
c
c--Write out pixel map
c
         WRITE (21,88001) ndensity, ntemp
88001    FORMAT ('# ',I5,I5)

         DO j = 1, ntemp
            DO i = 1, ndensity
               IF (npixmap(i,j).GT.0) THEN
                  pixmap(i,j) = LOG10(pixmap(i,j)/npixmap(i,j) * 1.4E-4)
               ELSE
                  pixmap(i,j) = -12.
               ENDIF
c            IF (i.GT.100 .AND. i.LT.120 .AND. j.EQ.400) pixmap(i,j) = 
c     &           -6.0
c            IF (i.EQ.100 .AND. j.EQ.100) pixmap(i,j) = -6.0
            END DO
c
c--Put one pixel at max of scale
c
            IF (j.EQ.ntemp) pixmap(ndensity,ntemp)=LOG10(3.0*1.4E-4)

            WRITE (21,"(256(es10.3,1x))") (pixmap(i,j),i=1,ndensity)
         END DO
         
         CLOSE(21)

      END DO


      CLOSE(15)

      STOP
      END

      SUBROUTINE quit(i)
      STOP
      END

      SUBROUTINE endrun
      CALL quit(0)
      END
