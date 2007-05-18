      SUBROUTINE extract
c************************************************************
c                                                           *
c  This subroutine extracts one dump from a binary output   *
c     file                                                  *
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
c      INCLUDE 'COMMONS/torq'
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

      CHARACTER*40 ifile(10), ofile
      CHARACTER*1 iok, iok2, iokm, iaddmhd
      
1000  FORMAT (A40)
1001  FORMAT (A1)

      nfile = 1
c
c--assume that MHD variable is B
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'Bvol'
      ENDIF
      PRINT *, 'name of file to be extracted ?'

      DO i = 1, nfile
         READ (*, 1000) ifile(i)
      END DO

      PRINT *, 'name of reduced file ?'
      READ (*, 1000) ofile
      OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

      PRINT *, 'dump to be extracted'
      READ *, ireduct
      PRINT *, 'do you want time reset to zero ?'
      READ (*,1001) iok
      PRINT *, 'do you want to reset centre of mass ?'
      READ (*,1001) iokm
      IF (imhd.EQ.idim) THEN
         PRINT *, 'do you want to add/reset magnetic fields ?'
         READ (*,1001) iaddmhd
      ENDIF

      image = 0
      imo = 0

      DO 15 k = 1, nfile

         OPEN (UNIT = 8, FILE = ifile(k), FORM = 'unformatted')
         PRINT *, 'reading file ', ifile(k)

         DO 10 image = 1, 1

            PRINT *, 'reading image number ', image
c
c--Read dump file
c
            CALL rdump(8,ichkl,1)
            IF (ichkl.EQ.1) PRINT*, 'ERROR READING DUMP FILE'

c
c--End reading of dump file
c--------------------------
c
      IF (image.EQ.ireduct) THEN

         IF (iokm.EQ.'y' .OR. iokm.EQ.'Y') THEN
            cmx = 0. 
            cmy = 0.
            cmz = 0.
            vcmx = 0. 
            vcmy = 0.
            vcmz = 0.
            pmtot = 0.
            DO i =1, npart
               pmassi = xyzmh(4,i)
               cmx = cmx + pmassi*xyzmh(1,i)
               cmy = cmy + pmassi*xyzmh(2,i)
               cmz = cmz + pmassi*xyzmh(3,i)
               vcmx = vcmx + pmassi*vxyzu(1,i)
               vcmy = vcmy + pmassi*vxyzu(2,i)
               vcmz = vcmz + pmassi*vxyzu(3,i)
               pmtot = pmtot + pmassi
            END DO
            cmx = cmx/pmtot
            cmy = cmy/pmtot
            cmz = cmz/pmtot
            vcmx = vcmx/pmtot
            vcmy = vcmy/pmtot
            vcmz = vcmz/pmtot
            DO i = 1, npart
               xyzmh(1,i) = xyzmh(1,i) - cmx
               xyzmh(2,i) = xyzmh(2,i) - cmy
               xyzmh(3,i) = xyzmh(3,i) - cmz
               vxyzu(1,i) = vxyzu(1,i) - vcmx
               vxyzu(2,i) = vxyzu(2,i) - vcmy
               vxyzu(3,i) = vxyzu(3,i) - vcmz
            END DO
         END IF

         PRINT *, 'writing image just read on output file'
         imo = imo + 1
         IF (iok.EQ.'y') gt = 0.0
         
         PRINT *, ' do you want to add point masses ? '
         READ (*,1001) iok2

         IF (iok2.EQ.'y') THEN
            PRINT *, ' reading from file externcluster '
            OPEN (UNIT=13, FILE='externcluster')
            READ (13,*) nptmassnew
            DO in = 1, nptmassnew
               inew = in + npart
               READ(13,*) xyzmh(1,inew), xyzmh(2,inew), 
     &                    xyzmh(3,inew), xyzmh(4,inew),
     &                    vxyzu(1,inew), vxyzu(2,inew), vxyzu(3,inew), 
     &                    xyzmh(5,inew)
               iphase(inew) = 1
               listpm(nptmass+in) = inew
               spinx(nptmass+in) = 0.
               spiny(nptmass+in) = 0.
               spinz(nptmass+in) = 0.
               vxyzu(4,inew) = tiny
               rho(inew) = tiny
               dgrav(inew) = 0.
            END DO
            n2 = nptmassnew
            hacc = xyzmh(5,inew)
            haccall = 0.2*hacc
            npart = npart + nptmassnew
            nptmass = nptmass + nptmassnew
         ENDIF
         PRINT *, ' do you want to change masses and temps ? '
         READ (*,1001) iok2
         IF (iok2.EQ.'y') THEN
            PRINT *, ' by how much ?'
            READ *, howmuch
            DO ii = 1, npart
               IF (iphase(ii).EQ.0) THEN
                  xyzmh(4,ii) = xyzmh(4,ii) * howmuch
c                        vxyzu(4,ii) = vxyzu(4,ii) * howmuch
               ENDIF
            END DO
         ENDIF

         IF (imhd.EQ.idim) THEN
            IF (iaddmhd.EQ.'y') THEN
               CALL setBfield
c               CALL hcalc ! to get B/rho from B
            ENDIF
         ENDIF               
c
c--Dump new file
c
      ifulldump = 0
      DO i = 1, npart
         isort(i) = i
         iorig(i) = i
      END DO
      nfullstep = 1
      PRINT *,'writing dump file'
      CALL wdump(7)
c
c--End writing of full dump file
c-------------------------------
c
      ENDIF

 10   CONTINUE

      CLOSE (8)

 15   CONTINUE      

      CLOSE (7)

      PRINT 2000, ofile, imo
 2000 FORMAT ('file ', A10, 'has been created and contains ', I3, 
     &     ' images')
 
      RETURN
      END
