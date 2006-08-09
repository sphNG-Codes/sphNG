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

      CHARACTER*40 ifile(10), ofile
      CHARACTER*1 iok, iok2, iokm, iaddmhd
      CHARACTER*100 fileident
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      DIMENSION nums1(8),nums2(8),nums3(8),nums4(8)
      DIMENSION Bxyz(3,imhd)
      
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
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Write ouput file
c
      READ (8, END=100) int1i,r1i,int2i,i1i,int3i
      IF (int1i.NE.int1) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         CALL quit
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default integer size wrong'
         CALL quit
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default real size wrong'
         CALL quit
      ENDIF
      READ (8, END=100) fileident
c
c--Single values
c
c--Default int
      READ (8, END=100) number
      IF (number.LT.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump: not enough default ints'
         CALL quit
      ENDIF
      READ (8, END=100) npart,n1,n2,nreassign,naccrete,nkill
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (8, END=100) number
      END DO
c--Default real
      READ (8, END=100) number
      IF (number.LT.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump: not enough default reals'
         CALL quit
      ENDIF
      IF (imhd.EQ.idim .AND. number.GE.18) THEN
         READ (8, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin, tmag, Bextx, Bexty, Bextz
      ELSE
         IF (imhd.EQ.idim) THEN
            WRITE(*,*) 'WARNING: dump does not contain external field'
            WRITE(*,*) '         (setting to zero)'
            Bextx = 0.
            Bexty = 0.
            Bextz = 0.
         ENDIF
         READ (8, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin
      ENDIF
c--real*4
      READ (8, END=100) number
c--real*8
      READ (8, END=100) number
      IF (number.LT.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump: nreal8 too small header'
         CALL quit
      ENDIF
      IF (imhd.EQ.idim) THEN
         IF (number.GT.3) THEN
            READ (8, END=100) udisti, umassi, utimei, umagfdi
         ELSE
            WRITE (*,*) 'WARNING: no mag field units in rdump'
            READ (8, END=100) udisti, umassi, utimei
            umagfdi = umagfd         
         ENDIF
      ELSE
         READ (8, END=100) udisti, umassi, utimei
      ENDIF
c
c--Arrays
c
c--Number of array lengths
c
      READ (8, END=100) number
      IF (number.LT.2 .OR. number.GT.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump'
         CALL quit
      ENDIF
c
c--Read array type 1 header
c
      READ (8, END=100) number8, (nums1(i), i=1,8)
      IF (number8.NE.npart) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong'
         CALL quit
      ENDIF
      npart = number8
      PRINT*,' npart = ',npart
c
c--Read array type 2 header
c
      READ (8, END=100) number8, (nums2(i), i=1,8)
      nptmass = number8
      PRINT*,' nptmasses = ',nptmass
c
c--Read array type 3 header
c
      IF (number.GE.3) THEN
         READ (8, END=100) number8, (nums3(i), i=1,8)
         IF (number8.GT.iradtrans .OR. number8.NE.0 .AND. 
     &        number8.NE.npart) THEN
            WRITE (*,*) 'ERROR 9 in rdump: iradtrans wrong ',number8,
     &           iradtrans,npart
            CALL quit
         ENDIF
         nradtrans = number8
      ENDIF
c
c--Read array type 4 header
c
      IF (number.GE.4) THEN
         READ (8, END=100) number8, (nums4(i), i=1,8)
         IF (number8.GT.imhd .OR. number8.NE.1 .AND. 
     &        number8.NE.npart) THEN
            WRITE (*,*) 'ERROR 10 in rdump: imhd wrong ',number8,
     &           imhd,npart
            CALL quit
         ENDIF
         nmhd = number8
      ENDIF
c
c--Read array type 1 arrays
c
c--Default int
      READ (8, END=100) (isteps(i), i=1, npart)
c--int*1
      READ (8, END=100) (iphase(i), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      DO j = 1, 5
         READ (8, END=100) (xyzmh(j,i), i=1, npart)
      END DO
      DO j = 1, 4
         READ (8, END=100) (vxyzu(j,i), i=1, npart)
      END DO      
c--skip unnecessary reals
      IF (nums1(6).GT.9) THEN
         DO j=1,nums1(6)-9
            READ (8, END=100)
         ENDDO
      ENDIF    
c--real*4
      READ (8, END=100) (rho(i), i=1, npart)
      IF (nlmax.EQ.1) THEN
         iread = 1
      ELSE
         iread = 2
         READ (8, END=100) (dgrav(i), i=1, npart)
         IF (gt.EQ.0.0) THEN
            DO j = 1, npart
               dgrav(j) = 0.
            ENDDO
         ENDIF
      ENDIF
c--skip unnecessary real*4's
      IF (nums1(7).GT.iread) THEN
         DO j=1,nums1(7)-iread
            READ (8, END=100)
         ENDDO
      ENDIF
c     READ (8, END=100) (alphaMM(1,i), i=1, npart)
c--real*8

c
c--Read array type 2 arrays
c
c--Default int
      READ (8, END=100) (listpm(i), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      READ (8, END=100) (spinx(i),i=1,nptmass)
      READ (8, END=100) (spiny(i),i=1,nptmass)
      READ (8, END=100) (spinz(i),i=1,nptmass)
      READ (8, END=100) (angaddx(i),i=1,nptmass)
      READ (8, END=100) (angaddy(i),i=1,nptmass)
      READ (8, END=100) (angaddz(i),i=1,nptmass)
      READ (8, END=100) (spinadx(i),i=1,nptmass)
      READ (8, END=100) (spinady(i),i=1,nptmass)
      READ (8, END=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

      IF (number.GE.3 .AND. nradtrans.EQ.npart 
     &    .AND. encal.EQ.'r') THEN
c
c--Array length 3 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 5
            READ (8, END=100) (ekcle(j,i), i=1, npart)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (number.GE.4 .AND. nmhd.EQ.npart .AND. imhd.EQ.idim) THEN
          PRINT *,' dump file contains magnetic fields...' 
c
c--Array length 4 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
c
c--read B field from dump
c
         IF (nums4(6).LT.3) THEN
            WRITE(*,*) 'ERROR: no MHD variables in dump file'
            CALL quit
         ENDIF
         DO j = 1, 3
            READ (8, END=100) (Bxyz(j,i), i=1, npart)
         END DO
c
c--read Euler potentials from dump if necessary
c
         IF (varmhd.EQ.'eulr') THEN
            IF (nums4(6).GE.5) THEN
               DO j = 1, 2
                  READ (8, END=100) (Bevolxyz(j,i), i=1, npart)
               END DO
            ELSE
               WRITE(*,*) 'ERROR: Cannot start Euler potentials run '
               WRITE(*,*) '       from non-Euler potentials dump'
               CALL quit
            ENDIF
         ELSEIF (varmhd.EQ.'Brho') THEN
c
c--convert from B to B/rho for evolution
c
            DO i=1,npart
               IF (rho(i).LE.0.) THEN
                  WRITE(*,*) 'ERROR: rho <= 0 in rdump, evolving B/rho'
                  CALL quit
               ENDIF
               DO j = 1,3
                  Bevolxyz(j,i) = Bxyz(j,i)/rho(i)
               ENDDO
            ENDDO
         ELSEIF (varmhd.EQ.'Bvol') THEN
c
c--Bevol = B if evolving B
c
            DO i=1,npart
               DO j = 1,3
                  Bevolxyz(j,i) = Bxyz(j,i)
               ENDDO
            ENDDO
         ELSE
            STOP 'unknown MHD variable in rdump'
         ENDIF
c--real*4

c--real*8

      ELSEIF (imhd.EQ.idim) THEN
          PRINT *,' no magnetic fields detected in dump file' 
      ENDIF

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
               
               IF (iaddmhd.EQ.'y') THEN
                  CALL setBfield
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

 10      CONTINUE

 100     CLOSE (8)

 15   CONTINUE      

      CLOSE (7)

      PRINT 2000, ofile, imo
 2000 FORMAT ('file ', A10, 'has been created and contains ', I3, 
     &     ' images')
 
      RETURN
      END
