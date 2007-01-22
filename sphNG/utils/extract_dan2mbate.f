      SUBROUTINE extract
c************************************************************
c                                                           *
c  This subroutine extracts one dump from a binary output   *
c     file                                                  *
c  This version converts from Dan's file format to sphNG    *
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
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/presb'
      
      CHARACTER*40 ifile(10), ofile
      CHARACTER*1 iok, iok2, iokm
      REAL dummy(idim)
      INTEGER nprint

1000  FORMAT (A40)
1001  FORMAT (A1)

      nfile = 1

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

      image = 1
      imo = 0
      CALL unit

      DO 15 k = 1, nfile

         OPEN (UNIT = 8, FILE = ifile(k), FORM = 'unformatted')
         PRINT *, 'reading file ', ifile(k)

         READ (8, END=100) gt, npart, nprint, gamma
         WRITE(*,*) 't = ',gt,' nprint = ',nprint,' gamma = ',gamma
c         npart = nprint
         n1 = npart
         n2 = 0
         nreassign = 0
         naccrete = 0
         nkill = 0
         rhozero = 1.0
         dtmax = 1.0
         tkin = 0.
         tgrav = 0.
         tterm = 0.
         anglostx = 0.
         anglosty = 0.
         anglostz = 0.
         specang = 0.
         ptmassin = 0.
         
         DO j=1,3
            READ (8, END=100) (xyzmh(j,i),i=1,npart)
         ENDDO
         DO j=1,3
            READ (8, END=100) (vxyzu(j,i),i=1,npart)
         ENDDO
         READ (8, END=100) (xyzmh(5,i),i=1,npart)
         READ (8, END=100) (dummy(i),i=1,npart)
         DO i=1,npart
            IF (i.le.5) WRITE(*,*) ' density ',i,' = ',dummy(i)
            rho(i) = dummy(i) ! type conversion
         ENDDO
         READ (8, END=100) (vxyzu(4,i),i=1,npart)
         READ (8, END=100) (xyzmh(4,i),i=1,npart)
         IF (imhd.EQ.idim) THEN
            WRITE(*,*) 'attempting to read MHD quantities'
c--skip alpha (art visc. parameters)
            DO j=1,3
               READ(8, END=100)
            ENDDO
            DO j=1,3
               READ(8, END=100) (Bxyz(j,i),i=1,npart)
            ENDDO
            varmhd = 'Bvol'
c
c--dump B directly
c
            DO i=1,npart
               DO j = 1,3
                  Bevolxyz(j,i) = Bxyz(j,i)
               ENDDO
            ENDDO
         ELSE
            WRITE(*,*) 'finished reading (hydro) file'
         ENDIF

         DO i=1,5
            WRITE(*,*) 'particle ',i,': RK2 = ',
     &                  vxyzu(4,i)/(rho(i)**(gamma-1.))
         ENDDO
         WRITE(*,*) 'Enter RK2 '
         READ(*,*) RK2
       
         DO i=1,npart
            iphase(i) = 0
            isteps(i) = 0
         ENDDO
c
c--set ptmass quantities to zero
c
         nptmass = 0
         DO i=1,nptmass
            listpm(i) = 0
            spinx(i) = 0.
            spiny(i) = 0.
            spinz(i) = 0.
            angaddx(i) = 0.
            angaddy(i) = 0.
            angaddz(i) = 0.
            spinadx(i) = 0.
            spinady(i) = 0.
            spinadz(i) = 0.
         ENDDO

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
c
c--Dump new file
c
      ifulldump = 0
      DO i = 1, npart
         isort(i) = i
         iorig(i) = i
      END DO
      nfullstep = 1
      CALL wdump(7)
c
c--End writing of full dump file
c-------------------------------
c
            ENDIF

 100     CLOSE (8)

 15   CONTINUE      

      CLOSE (7)

      PRINT 2000, ofile, imo
 2000 FORMAT ('file ', A10, 'has been created and contains ', I3, 
     &     ' images')
 
      RETURN
      END
