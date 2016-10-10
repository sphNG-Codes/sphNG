      SUBROUTINE makeplot
c************************************************************
c                                                           *
c  This subroutine reads in an existing dump file and       *
c     brings dead particles to life as dust particles.      *
c     Originally used to add dust to a settled disc.        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/makeplt'

      CHARACTER*11 ifile(1000), ofile

      CHARACTER*4 dummy4
      CHARACTER*3 nsuffix
      CHARACTER*1 idesire

      iplotdump = 1

1000  FORMAT (A40)
1001  FORMAT (A1)
c
c--assume that MHD variable is B
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'Bvol'
      ENDIF

      PRINT *, 'Enter number of files'
      READ (*, *) nfilenames

      PRINT *, 'Do you want to subtract time (0.0=no)?'
      READ (*, *) offset

      PRINT *, 'Names of binary dump files ?'
      DO k = 1, nfilenames
         READ (*, 1000) ifile(k)
      END DO

      PRINT *, 'Name of output binary dump file ?'
      READ (*, 1000) ofile

      PRINT *, 'Do you want:'
      PRINT *, '   (a) azimuthal averaging, with dust density'
      PRINT *, '   (r) azimuthal averaging, with dust-to-gas ratio'
      PRINT *, '   (h) double h'
      PRINT *, '   (v) compute velocity dispersion: v_z(1)'
      READ (*,1001) idesire

      READ (ofile,88) dummy4, ioffset
 88   FORMAT(A4,I3)

      print *,'ioffset ',ioffset

      DO k = 1, nfilenames

         ifilenum = k + ioffset - 1

         print *,'ifilenum ',ifilenum,k,ioffset

         WRITE (nsuffix,('(I3)')) ifilenum

         print *,'nsuffix ',nsuffix

         IF (ifilenum.LE.9) THEN
            ofile = ofile(1:4) // '00' // nsuffix(3:3)
         ELSEIF (ifilenum.LE.99) THEN
            ofile = ofile(1:4) // '0' // nsuffix(2:3)
         ELSEIF (ifilenum.LE.999) THEN
            ofile = ofile(1:4) // nsuffix(1:3)
         ENDIF

         print *,'Opening ',ofile

         IF (idesire.NE.'v')
     &        OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

         imax = 1073741824
         imaxstep = imax/2

         nactive = 0
         DO i = 1, idim
            isort(i) = i
            iorig(i) = i
         END DO

         OPEN (UNIT = 8, FILE = ifile(k), FORM = 'unformatted')
         PRINT *, 'reading file ', ifile(k)
c
c--Read dump file
c
         CALL options

         file1 = ofile
         print *, 'Name0=', file1

         CALL rdump_wrapper(8,ichkl,1)
         IF (ichkl.EQ.1) THEN
            PRINT*, 'ERROR READING DUMP FILE'
            CALL quit(0)
         ENDIF

         gt = gt - offset
         print *, 'New time ',gt
c
c--End reading of dump file
c--------------------------
c
c--Do the azimutal average, change of smoothing length, or velocity
c     dispersion calculation
c
         veldisp = 0.
         nveldisp = 0
         DO i = 1, npart
            r = DSQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
            vr = (vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i))/r
c            IF (xyzmh(2,i).LT.0.) iphase(i) = -1
            IF (idesire.EQ.'a' .OR. idesire.EQ.'r') THEN
               xyzmh(1,i) = r
               xyzmh(2,i) = 0.0
               xyzmh(4,i) = xyzmh(4,i)/(2*3.14159265*r) * xyzmh(5,i)*1.4
               vxyzu(1,i) = vr
               vxyzu(2,i) = 0.0
            ELSEIF (idesire.EQ.'h') THEN
               xyzmh(5,i) = xyzmh(5,i)*2.0
            ELSEIF (idesire.EQ.'v') THEN
               IF (r.GT.0.99 .AND. r.LT.1.01) THEN
                  veldisp = veldisp + vxyzu(3,i)**2
                  nveldisp = nveldisp + 1
               ENDIF
            ENDIF
            IF (iphase(i).EQ.0) THEN
               xtype(1,i) = 1.
               xtype(2,i) = 0.
            ELSEIF (iphase(i).GE.11) THEN
               xtype(1,i) = 0.
               xtype(2,i) = 1.
            ENDIF
         END DO
c
c--Normalise dust masses by gas masses
c
         IF (idesire.EQ.'r') THEN

            udens = umass/(udist**3)
C$OMP PARALLEL DO default(none)
C$OMP& shared(npart,iphase,xyzmh,rho,udens,xtype)
C$OMP& private(i,j,ngasclose,gasmeandensity,r2)
            DO i = 1, npart
               IF (MOD(i,10000).EQ.0) print *,i

               IF (iphase(i).GE.11) THEN
                  ngasclose = 0
                  gasmeandensity = 0.
c
c--Loop over all gas particles to find those that are close (d(r_xy)<h)
c
                  DO j = 1, npart
                     r2 = (xyzmh(1,i)-xyzmh(1,j))**2 + 
     &                    (xyzmh(3,i)-xyzmh(3,j))**2
                     IF (r2.LT.xyzmh(5,i)**2) THEN
                        IF (iphase(j).EQ.0) THEN
                           ngasclose = ngasclose + 1
                           gasmeandensity = gasmeandensity +rho(j)*udens
                        ENDIF
                     ENDIF
                  END DO
                  IF (ngasclose.GT.0) THEN
c
c--NOTE: Because the dust masses are used for the dust-to-gas ratio, this
c     will dominate the velocities if velocity vectors are plotted in
c     SPLASH.  To avoid this, can divide by a large number (or don't plot
c     velocity vectors with dust-to-gas ratio).
c
c                  xtype(2,i) = xtype(2,i)/
c     &                 (gasmeandensity/ngasclose)*1.0E-10
                     xyzmh(4,i) = xyzmh(4,i)/(gasmeandensity/ngasclose)
                  ENDIF
               ENDIF
            END DO
C$OMP END PARALLEL DO
         ENDIF
c
c--Finish computing velocity dispersion, and exit without dumping
c
         IF (idesire.EQ.'v') THEN
            veldisp = SQRT(veldisp/nveldisp)
            PRINT *,'Veldisp = ',veldisp,nveldisp
            GOTO 200
         ENDIF
c
c--Write output dump file
c
         DO i = 1, npart
            IF (iphase(i).GE.0) nactive = nactive + 1
         END DO

         ifulldump = 0
         nfullstep = 1
         PRINT *,'writing dump file'
         CALL wdump_wrapper(7)
c
c--End writing of full dump file
c-------------------------------
c
         CLOSE (7)

         PRINT 2000, ofile
 2000    FORMAT ('file ', A10, 'has been created')

 200     CLOSE (8)

      END DO
 
      RETURN
      END
