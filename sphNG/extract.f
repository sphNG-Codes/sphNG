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
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/physcon'

      CHARACTER*11 ifile(10), ofile
      CHARACTER*1 iok, iok2, iokm, iaddmhd, iexs
      
1001  FORMAT (A1)
1002  FORMAT (A11)

      nfile = 1
c
c--assume that MHD variable is B
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'Bvol'
      ENDIF
      PRINT *, 'name of file to be extracted ?'

      DO i = 1, nfile
         READ (*, 1002) ifile(i)
      END DO

      PRINT *, 'name of reduced file ?'
      READ (*, 1002) ofile
      print *,'xxx',ofile,'xxx',ifile(1)
      OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

      PRINT *, 'name of corresponding ifile ?'
      READ (*, 1002) inname


      PRINT *, 'dump to be extracted'
      READ *, ireduct
      PRINT *, 'do you want time reset to zero ?'
      READ (*,1001) iok
      PRINT *, 'do you want to reset (m) centre of mass or (s) sink(1)?'
      READ (*,1001) iokm
      IF (imhd.EQ.idim) THEN
         PRINT *, 'do you want to add/reset magnetic fields ?'
         READ (*,1001) iaddmhd
      ENDIF

      imo = 0

      DO 15 k = 1, nfile

         OPEN (UNIT = 8, FILE = ifile(k), FORM = 'unformatted')
         PRINT *, 'reading file ', ifile(k)
c
c--Read dump file
c
         CALL options

         file1 = ofile
         print *, 'Name0=', file1

         CALL rdump_wrapper(8,ichkl,1)
         IF (ichkl.EQ.1) PRINT*, 'ERROR READING DUMP FILE'
c
c--End reading of dump file
c--------------------------
c
         PRINT *, ' do you want to reduce the number of particles? '
         READ (*,1001) iok2
         IF (iok2.EQ.'y') THEN
            iseed = -6485
            CALL ran1(iseed)
            PRINT *, ' by how much (e.g. 0.1 means remove 90%) ?'
            READ *, howmuch
            nlistinactive = 0
            icount = 0
            DO ii = 1, npart
               IF (iphase(ii).EQ.0) THEN
                  IF (ran1(1).LT.howmuch ) THEN
                     icount = icount + 1
c
c--Increase h to compensate for loss of particles, and increase mass
c     to compensate for loss of particles
c
                     xyzmh(4,ii) = xyzmh(4,ii)/howmuch
                     xyzmh(5,ii) = xyzmh(5,ii)*howmuch**(-1.0/3.0)
                     IF (icount.NE.ii) CALL move_particle(ii,icount)
                  ENDIF
               ENDIF
            END DO
            npart = icount
            n1 = npart
         ENDIF

         IF (iokm.EQ.'m' .OR. iokm.EQ.'s') THEN
            IF (iokm.EQ.'m') THEN
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
            ELSE
               cmx = xyzmh(1,listpm(1))
               cmy = xyzmh(2,listpm(1))
               cmz = xyzmh(3,listpm(1))
               vcmx = 0.
               vcmy = 0.
               vcmz = 0.
            ENDIF
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
            PRINT *, ' from file ? '
            READ (*,1001) iok2
            IF (iok2.EQ.'y') THEN
               PRINT *, ' reading from file externcluster '
               OPEN (UNIT=13, FILE='externcluster')
               READ (13,*) nptmassnew
               DO in = 1, nptmassnew
                  inew = in + npart
                  READ(13,*) xyzmh(1,inew), xyzmh(2,inew), 
     &                 xyzmh(3,inew), xyzmh(4,inew),
     &                 vxyzu(1,inew), vxyzu(2,inew), vxyzu(3,inew), 
     &                 xyzmh(5,inew)
                  iphase(inew) = 1
                  listpm(nptmass+in) = inew
                  spinx(nptmass+in) = 0.
                  spiny(nptmass+in) = 0.
                  spinz(nptmass+in) = 0.
                  vxyzu(4,inew) = tiny
                  rho(inew) = tiny
                  dgrav(inew) = 0.
                  iunique(inew) = inew
               END DO
               n2 = nptmassnew
               hacc = xyzmh(5,inew)
               haccall = 0.2*hacc
               npart = npart + nptmassnew
               nptmass = nptmass + nptmassnew
            ELSE
               WRITE (*, 99765)
99765          FORMAT (' Are there any existing sink particles? (y/n)')
               READ (*, 1001) iexs
               totptmass = 0.
               IF (iexs.EQ.'y') THEN
                  DO i = 1, nptmass
                     totptmass = totptmass + xyzmh(4,listpm(i))
                  END DO
               ELSE
                  nptmass = 0
               ENDIF
               npttemp = nptmass
               WRITE (*, 99086)
99086          FORMAT (' Enter number of point masses')
               READ (*,*) numpt
               WRITE (*, 99091)
99091          FORMAT (' Enter type of point mass (1,2,3,4,5)')
               READ (*,*) ipttype
               initialptm = ipttype
              
               DO i = 1, numpt
                  npart = npart + 1
                  iphase(npart) = ipttype
                  listpm(nptmass+i) = npart
                  listrealpm(npart) = i+nptmass
                  WRITE (*, 99087)
99087             FORMAT (' Enter positions (x,y,z)')
                  READ (*,*) tx,ty,tz
                  xyzmh(1,npart) = tx
                  xyzmh(2,npart) = ty
                  xyzmh(3,npart) = tz
                  WRITE (*, 99088)
99088             FORMAT (' Enter velocities (vx,vy,vz)')
                  READ (*,*) tvx,tvy,tvz
                  vxyzu(1,npart) = tvx
                  vxyzu(2,npart) = tvy
                  vxyzu(3,npart) = tvz
                  WRITE (*, 99089)
99089             FORMAT (' Enter mass')
                  READ (*,*) tmass
                  xyzmh(4,npart) = tmass
                  totptmass = totptmass + tmass
                  spinx(nptmass+i) = 0.
                  spiny(nptmass+i) = 0.
                  spinz(nptmass+i) = 0.
                  angaddx(nptmass+i) = 0.
                  angaddy(nptmass+i) = 0.
                  angaddz(nptmass+i) = 0.
                  IF (iexs.EQ.'n' .AND. ipttype.NE.5) THEN
99090                FORMAT (' Enter Outer Accretion radius')
 556                 WRITE (*, 99090)
                     READ (*,*) hacc
99092                FORMAT (' Enter Inner Accretion radius')
                     WRITE (*, 99092)
                     READ (*,*) haccall
                     IF (haccall.GT.hacc) THEN
                        WRITE(*,*) 'Inner <= Outer Radius'
                        GOTO 556
                     ENDIF
                     xyzmh(5,npart) = hacc
                     print *, 'h = ',xyzmh(5, npart)
                  ENDIF
               END DO
               specang = SQRT(totptmass)
               ptmassin = totptmass
               nptmass = nptmass + numpt
c--Moves particles near new and existing point masses to avoid the
c--disc struggling to adapt to the sudden addition.
               IF ((ibound.EQ.102 .OR. ibound.EQ.103) 
     &              .AND. ipttype.EQ.5) THEN
                  print *, 'Enter a radius for the planet surface'
                  READ (*,*) rplanet
                  hacc = 1.0E-8
                  haccall = 1.0E-8
                  xyzmh(5,npart) = hacc
                  print *, 'h = ',xyzmh(5, npart)
                  DO j = npttemp, nptmass
                     DO i = 1, npart
                        IF (i.NE.listpm(j)) THEN
 789                    rloc = sqrt((xyzmh(1,i)-xyzmh(1,listpm(j)))**2 +
     &                       (xyzmh(2,i)-xyzmh(2,listpm(j)))**2 +
     &                       (xyzmh(3,i)-xyzmh(3,listpm(j)))**2)
                        IF (rloc.LT.2.1*(rplanet+0.01)) THEN
                  print *, i, rloc, xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)
                           xyzmh(3,i) = xyzmh(3,i)*1.5
                           GOTO 789
                        ENDIF
                        ENDIF
                     END DO
                  END DO
                  print *, 'Enter a uniform timestep limit (in RHill)'
                  READ (*,*) uniformtslim
                  uniformtslim2 = uniformtslim*uniformtslim
               ENDIF
            ENDIF
         ENDIF

         PRINT *, ' do you want to rotate gas ?'
         READ (*,1001) iok2
         IF (iok2.EQ.'y') THEN
            PRINT *, ' set inner radius (>=0)?'
            READ (*,*)  rinnercut
            DO i = 1, npart
               IF (iphase(i).EQ.0) THEN
                  xtmp = xyzmh(1,i)
                  ytmp = xyzmh(2,i)
                  ztmp = xyzmh(3,i)
                  radius = SQRT(xtmp**2+ytmp**2+ztmp**2)
                  IF (radius.LT.rinnercut) THEN
                     iphase(i) = -1
                     nlistinactive = nlistinactive + 1
                     listinactive(nlistinactive) = i
                  ELSE
                  vxtmp = vxyzu(1,i)
                  vytmp = vxyzu(2,i)
                  vztmp = vxyzu(3,i)
c Alison - added rotation matrix
                  posang = (-176d0)*pi/180d0
c Rotation if LHS is nearest to us (assumed)
                  ang_inc = -142d0 * pi/ 180d0
c Rotation if LHS is FATHEREST from us
c                  ang_inc = +142d0 * pi/ 180d0
c rotfac is -1 for clockwise rotation on the sky
                  rotfac = 1.0
                  vxyzu(1,i) = rotfac*(vxtmp*COS(posang) + 
     &                 vytmp*SIN(posang)*COS(ang_inc) +
     &                 vztmp*SIN(posang)*SIN(ang_inc))
                  vxyzu(2,i) = rotfac*(-vxtmp*SIN(posang) +
     &                 vytmp*COS(ang_inc)*COS(posang) + 
     &                 vztmp*COS(posang)*SIN(ang_inc))
                  vxyzu(3,i) = rotfac*(-vytmp*SIN(ang_inc) +
     &                 vztmp*COS(ang_inc))

                  xyzmh(1,i) = -(xtmp*COS(posang) +
     &                 ytmp*SIN(posang)*COS(ang_inc) +
     &                 ztmp*SIN(posang)*SIN(ang_inc))
                  xyzmh(2,i) = -(-xtmp*SIN(posang) +
     &                 ytmp*COS(ang_inc)*COS(posang) +
     &                 ztmp*COS(posang)*SIN(ang_inc))
                  xyzmh(3,i) =-(-ytmp*SIN(ang_inc) +ztmp*COS(ang_inc))
                  ENDIF
               ENDIF
            END DO
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

         CALL wrinsph

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
         CALL wdump_wrapper(7)
c
c--End writing of full dump file
c-------------------------------
c
 10   CONTINUE

      CLOSE (8)

 15   CONTINUE      

      CLOSE (7)

      PRINT 2000, ofile, imo
 2000 FORMAT ('file ', A10, 'has been created and contains ', I3, 
     &     ' images')
 
      RETURN
      END
