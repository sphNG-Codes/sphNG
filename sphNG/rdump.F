      SUBROUTINE rdump(idisk1, ichkl, iskip)
c************************************************************
c                                                           *
c  This routine reads a dump into memory                    *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      REAL*8 umassi, udisti, utimei, umagfdi

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/gtdble'
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
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/xforce'      

      DIMENSION itempsort(idim)
      EQUIVALENCE (itempsort, next1)

      CHARACTER*7 where
      CHARACTER*100 fileident
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      DIMENSION nums1(8),nums2(8),nums3(8),nums4(8)
      DIMENSION Bxyz(3,imhd)

      DATA icall/2/
      DATA where/'rdump'/
c     
c--Read
c
      IF (itrace.EQ.'all') WRITE (*, 99001)
99001 FORMAT (' entry subroutine rdump')

      ichkl = 0
      icount = 0
c
c--Dump file
c-------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Read ouput file
c

      print *,'reading header'

      READ (idisk1, END=100) int1i,r1i,int2i,i1i,int3i

      print *,'have read header'

      IF (int1i.NE.int1) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         CALL quit
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default real size wrong'
         CALL quit
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default integer size wrong'
         CALL quit
      ENDIF
      READ (idisk1, END=100) fileident
c
c--Single values
c
c--Default int
      READ (idisk1, END=100) number
      IF (number.LT.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump: not enough default ints'
         CALL quit
      ENDIF
      IF (number.EQ.7) THEN
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill,
     &        nblocks
      ELSE
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill
         nblocks = 1
      ENDIF
      IF (npart.GT.idim .AND. nblocks.EQ.1) THEN
         WRITE (*,*) 'ERROR in rdump: npart>idim'
         CALL quit
      ENDIF
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (idisk1, END=100) number
      END DO
c--Default real
      READ (idisk1, END=100) number
      IF (number.LT.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump: not enough default reals'
         CALL quit
      ENDIF
c--need to read hzero for problems which utilise 
c  the monaghan-style external boundary force
c  D.Price 14.5.07

      IF (number.GE.19) THEN
         READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin, tmag, Bextx, Bexty, Bextz, hzero

         WRITE(*,*) ' read hzero = ',hzero
         IF (imhd.EQ.idim) THEN
            WRITE(*,*) 'External field found, Bext = ',Bextx,Bexty,Bextz
         ENDIF
      ELSE
         IF (iexf.EQ.10) THEN
            WRITE(*,*) 'ERROR: hzero not found in rdump'
            CALL quit
         ENDIF
         IF (imhd.EQ.idim .AND. number.GE.18) THEN
            READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin, tmag, Bextx, Bexty, Bextz
            WRITE(*,*) 'External field found, Bext = ',Bextx,Bexty,Bextz
         ELSE
            IF (imhd.EQ.idim) THEN
              WRITE(*,*) 'WARNING: dump does not contain external field'
              WRITE(*,*) '         (setting to zero)'
               Bextx = 0.
               Bexty = 0.
               Bextz = 0.
            ENDIF
            READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin
         ENDIF
      ENDIF

c--real*4
      READ (idisk1, END=100) number
c--real*8
      READ (idisk1, END=100) number
      IF (number.LT.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump: nreal8 too small header'
         CALL quit
      ENDIF
      IF (imhd.EQ.idim) THEN
         IF (number.GT.3) THEN
            READ (idisk1, END=100) udisti, umassi, utimei, umagfdi
         ELSE
            WRITE (*,*) 'WARNING: no mag field units in rdump'
            READ (idisk1, END=100) udisti, umassi, utimei
            umagfdi = umagfd         
         ENDIF
      ELSE
         READ (idisk1, END=100) udisti, umassi, utimei
      ENDIF
c
c--Arrays
c
c--Number of array lengths
c     This has up to 4 blocks (hydro variables, sink particles, radiative
c     transfer, and MHD) repeated nblocks times for each of nblocks MPI
c     processes
c
      READ (idisk1, END=100) number
#ifdef MPI
      print *,iproc,': Number of array lengths & blocks',number,nblocks
#else
      print *,'Number of array lengths & blocks ',number,nblocks
#endif
c
c--Expected
c
      numberexp = 2
      IF (encal.EQ.'r') numberexp = 3
      IF (imhd.EQ.idim) numberexp = 4

      IF (number.LT.2 .OR. number/nblocks.LT.2 .OR. 
     &     numberexp.EQ.4 .AND.number/nblocks.NE.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump ', number, numberexp
         CALL quit
#ifdef MPI
      ELSEIF (nblocks.NE.1 .AND. nblocks.NE.numproc) THEN
         WRITE (*,*) 'ERROR in rdump: nblocks.NE.numproc ',
     &        nblocks,numberexp,numproc
         WRITE (*,*) 'Cannot start up job with different number of MPI',
     &        ' processes'
         CALL quit
#endif
      ENDIF
      numberarray = number/nblocks
c
c--Loop over blocks
c
      DO iblock = 1, nblocks
#ifdef MPI
         print *,iproc,': Doing block ',iblock
#else
         print *,'Doing block ',iblock
#endif
c
c--Read array type 1 header
c
      READ (idisk1, END=100) number8, (nums1(i), i=1,8)
      IF (number8.LT.npart .AND. nblocks.EQ.1 .OR. number8.GT.idim) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong',number8,
     &           npart, idim
         CALL quit
      ENDIF
#ifdef MPI
      IF (nblocks.EQ.1 .OR. nblocks.NE.1 .AND. iproc.EQ.iblock-1) THEN
#endif
      nhydro = number8
#ifdef MPI
      ENDIF
#endif

#ifdef MPIDEBUG
      print *,iproc,': nhydro ',number8
#endif
c
c--Read array type 2 header
c
      READ (idisk1, END=100) number8, (nums2(i), i=1,8)
      IF (number8.GT.iptdim) THEN
         WRITE(*,*) 'ERROR 8.5 in rdump: nptmass exceeds dimensions!!',
     &            number8,iptdim
         CALL quit
      ENDIF
#ifdef MPI
      IF (nblocks.EQ.1 .OR. nblocks.NE.1 .AND. iproc.EQ.iblock-1) THEN
#endif
      nptmass = number8
#ifdef MPI
      ENDIF
#endif
#ifdef MPIDEBUG
      print *,iproc,': nptmass ',number8
#endif
c
c--Read array type 3 header
c
      IF (numberarray.GE.3) THEN
         READ (idisk1, END=100) number8, (nums3(i), i=1,8)
         IF (iradtrans.EQ.idim .AND. number8.GT.iradtrans) THEN
            WRITE (*,*) 'ERROR 9 in rdump: iradtrans wrong ',number8,
     &           iradtrans,npart
            CALL quit
         ENDIF
         nradtrans = number8
      ELSE
         nradtrans = 0
      ENDIF
c
c--Read array type 4 header
c
      IF (numberarray.GE.4) THEN
         READ (idisk1, END=100) number8, (nums4(i), i=1,8)
         IF (imhd.EQ.idim.AND.number8.GT.imhd) THEN
            WRITE (*,*) 'ERROR 10 in rdump: imhd wrong ',number8,
     &           imhd,npart
            CALL quit
         ENDIF
         nmhd = number8
      ELSE
         nmhd = 0
      ENDIF

c      print *,'nums4(7) ',nums4(7)
c      nums4(7) = 4

#ifdef MPI
c
c--IF MPI code, then need to skip blocks that are not for this process
c
      IF (nblocks.NE.1 .AND. iproc.NE.iblock-1) THEN
#ifdef MPIDEBUG
         print *,iproc,': Skipping block ',iblock
#endif
         DO j = 1, 8
            DO k = 1, nums1(j)
               READ (idisk1, END=100)
            END DO
         END DO
         DO j = 1, 8
            DO k = 1, nums2(j)
               READ (idisk1, END=100)
            END DO
         END DO
         IF (numberarray.GE.3) THEN
            DO j = 1, 8
               DO k = 1, nums3(j)
                  READ (idisk1, END=100)
               END DO
            END DO
            IF (numberarray.GE.4) THEN
               DO j = 1, 8
                  DO k = 1, nums4(j)
                     READ (idisk1, END=100)
                  END DO
               END DO
            ENDIF
         ENDIF
      ELSE
#endif
c
c--Read array type 1 arrays
c
#ifdef MPIDEBUG
         print *,iproc,': Reading data block'
#endif
c--Default int
         READ (idisk1, END=100) (isteps(i), i=icount+1, icount+nhydro)
#ifdef MPIDEBUG
         print *,iproc,': read isteps ',isteps(1)
#endif
c--int*1
         READ (idisk1, END=100) (iphase(i), i=icount+1, icount+nhydro)
c--int*2

c--int*4

c--int*8
         IF (nums1(5).GE.1) THEN
            READ (idisk1, END=100) (iunique(i), i=icount+1,
     &           icount+nhydro)
         ENDIF
         IF (nums1(5).LT.1 .OR. iunique(icount+1).EQ.0) THEN
            print *,'Read iunique ',(iunique(i), i=icount+1,icount+10)
            DO i=icount+1,icount+nhydro
               iunique(i) = i
            END DO
         ENDIF
c--Default real
         DO j = 1, 5
            READ (idisk1, END=100) (xyzmh(j,i), i=icount+1, 
     &           icount+nhydro)
         END DO
         DO j = 1, 4
            READ (idisk1, END=100) (vxyzu(j,i), i=icount+1, 
     &           icount+nhydro)
         END DO
#ifdef MPIDEBUG
         print *,iproc,': Read xyzmh, vxyzu'
#endif
c--skip unnecessary reals
         IF (nums1(6).GT.9) THEN
            DO j=1,nums1(6)-9
               READ (idisk1, END=100)
            ENDDO
         ENDIF
c--real*4
         READ (idisk1, END=100) (rho(i), i=icount+1, icount+nhydro)
         iread = 1
         IF (nlmax.EQ.1) THEN
            IF (nums1(7).GE.2) THEN
               READ (idisk1, END=100) (gradhs(1,i), i=icount+1, 
     &              icount+nhydro)
               iread = iread + 1
            ELSEIF (nums1(7).GE.3) THEN
               READ (idisk1, END=100) (gradhs(2,i), i=icount+1, 
     &              icount+nhydro)
               iread = iread + 1
            ENDIF
         ELSE
            iread = 2
            READ (idisk1, END=100) (dgrav(i), i=icount+1, icount+nhydro)
            IF (gt.EQ.0.0) THEN
               DO j = icount+1, icount+nhydro
                  dgrav(j) = 0.
               ENDDO
            ENDIF
         ENDIF
         IF (nums1(7).GT.iread) THEN
            iread = iread + 1
            READ (idisk1, END=100) (alphaMM(1,i), i=icount+1, 
     &           icount+nhydro)
         ENDIF
c--skip unnecessary real*4's
         IF (nums1(7).GT.iread) THEN
            DO j=1,nums1(7)-iread
               READ (idisk1, END=100)
            ENDDO
         ENDIF
c--real*8

#ifdef MPIDEBUG
         print *,iproc,': Read hydro'
#endif
c
c--Read array type 2 arrays
c
c--Default int
         READ (idisk1, END=100) (listpm(i), i=1,nptmass)

         DO i = 1, nptmass
            write (*,*) 'Setting ',i,listpm(i),xyzmh(5,listpm(i)),hacc,
     &           iphase(listpm(i))
            xyzmh(5,listpm(i)) = hacc
         END DO


c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         READ (idisk1, END=100) (spinx(i),i=1,nptmass)
         READ (idisk1, END=100) (spiny(i),i=1,nptmass)
         READ (idisk1, END=100) (spinz(i),i=1,nptmass)
         READ (idisk1, END=100) (angaddx(i),i=1,nptmass)
         READ (idisk1, END=100) (angaddy(i),i=1,nptmass)
         READ (idisk1, END=100) (angaddz(i),i=1,nptmass)
         READ (idisk1, END=100) (spinadx(i),i=1,nptmass)
         READ (idisk1, END=100) (spinady(i),i=1,nptmass)
         READ (idisk1, END=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

         IF (numberarray.GE.3 .AND. nradtrans.EQ.nhydro 
     &        .AND. encal.EQ.'r') THEN
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
               READ (idisk1, END=100) (ekcle(j,i), i=icount+1, 
     &              icount+nradtrans)
            END DO
c--real*4

c--real*8

         ENDIF
         IF (numberarray.GE.4 .AND. nmhd.EQ.nhydro .AND. 
     &        imhd.EQ.idim) THEN
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
               READ (idisk1, END=100) (Bxyz(j,i), i=icount+1, 
     &              icount+nmhd)
            END DO
#ifdef MPIDEBUG
         print *,iproc,': Read Bxyz'
#endif
c
c--read Euler potentials from dump if necessary
c
            IF (varmhd.EQ.'eulr') THEN
               IF (nums4(6).GE.5) THEN
                  DO j = 1, 2
                     READ (idisk1, END=100) (Bevolxyz(j,i), i=icount+1, 
     &                    icount+nmhd)
                  END DO
#ifdef MPIDEBUG
         print *,iproc,': Read Bevolxyz'
#endif
               ELSE
                  WRITE(*,*) 'ERROR: Cannot start Euler potentials run'
                  WRITE(*,*) '       from non-Euler potentials dump'
                  CALL quit
               ENDIF
            ELSEIF (varmhd.EQ.'Brho') THEN
c
c--convert from B to B/rho for evolution
c
               DO i=icount+1,icount+nmhd
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
               DO i=icount+1,icount+nmhd
                  DO j = 1,3
                     Bevolxyz(j,i) = Bxyz(j,i)
                  ENDDO
               ENDDO
            ELSE
               STOP 'unknown MHD variable in rdump'
            ENDIF
c--real*4
            iread = 0
            IF (nums4(7).GT.4.AND.ifsvi.EQ.6) THEN
c--skip current/div B
               DO j = 1, 4
                  iread = iread + 1
                  READ (idisk1, END=100)
               END DO
c--read artificial resistivity parameter
               READ (idisk1, END=100) (alphaMM(isizealphaMM,i), 
     &              i=icount+1, icount+nmhd)
               iread = read + 1
            ELSE
               IF (ifsvi.EQ.6) THEN
           WRITE(*,*) 'WARNING: resistivity parameter not found in dump'
               ENDIF
               WRITE(*,*) ' => setting alphaB = alphamin everywhere'
               DO i=icount+1,icount+nmhd
                  alphaMM(2,i) = alphamin(2)
               ENDDO
            ENDIF
            DO i = iread + 1, nums4(7)
               READ (idisk1, END=100)
            END DO
c--real*8
         ELSEIF (imhd.EQ.idim) THEN
            WRITE(*,*) ' no magnetic fields detected in dump file' 
         ENDIF

#ifdef MPIDEBUG
         print *,iproc,': Read MHD'
#endif

#ifdef MPI
c
c--Split particles between different MPI processes
c
c      GOTO 10000

c      IF (nblocks.NE.1) THEN
c         IF (iproc.EQ.iblock-1) THEN
            IF (nblocks.EQ.1) THEN
               CALL split(npart,n1,icount,nhydro,isteps,iphase,xyzmh,
     &              vxyzu,rho,dgrav,Bevolxyz,Bxyz,gradhs,alphaMM,
     &              iunique)
            ELSE
               print *,iproc,': MAINTAIN PARTICLE SPLIT ',npart,nhydro
               npart = nhydro
               n1 = nhydro
            ENDIF
c         ENDIF
c      ENDIF

10000 CONTINUE

      IF (iblock.EQ.nblocks) print *,iproc,': KEEPING ',npart
c
c--End if to allow skipping of blocks of file belonging to other processes
c
      ENDIF
#endif
c
c--End loop over blocks
c
      END DO
#ifdef MPI
c
c--Set global number of sink particles
c
      CALL MPI_ALLREDUCE(nptmass,nptmasstot,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
#endif

      IF (iskip.EQ.1) RETURN
c
c--End reading of dump file
c--------------------------
c
      gtdouble = DBLE(gt)
c
c--Sort particles to ensure most efficient running.  Note that this 
c     should not be visible to the outside observer.  In other words,
c     an array must be kept of original list of particles and this
c     must be used to index *ANY* value from an array which is written 
c     to the outside.  This requires modification to almost every output
c     line in the code.  Done 21 Nov 2000.
c
      xminimum = 1.0E+30
      xmaximum = -1.0E+30
      DO i = 1, npart
         xmaximum = MAX(xmaximum, xyzmh(1,i))
         xminimum = MIN(xminimum, xyzmh(1,i))
      END DO
      xrange = xmaximum-xminimum
      istepmin = imax
      istepmax = 0
      DO i = 1, npart
         llist(i) = i
         IF (iphase(i).EQ.-1) THEN
            tempsort(i) = LOG(REAL(imax))/LOG(2.0)
            istepmax = imax
         ELSE
            IF (gt.EQ.0.0 .OR. isteps(i).EQ.0) THEN
               tempsort(i) = (xyzmh(1,i)-xminimum)/xrange
            ELSE
               tempsort(i) = LOG(REAL(isteps(i)))/LOG(2.0) +
     &           (xyzmh(1,i)-xminimum)/xrange
            ENDIF
            istepmin = MIN(istepmin, isteps(i))
            istepmax = MAX(istepmax, isteps(i))
         ENDIF
      END DO
c
c--Initialise timesteps to be consistent with integration method
c
      IF (individualtimesteps.EQ.0) THEN
         DO i = 1, npart
            isteps(i) = istepmin
         END DO
      ELSEIF (individualtimesteps.EQ.1) THEN
         DO i = 1, npart
            IF (iphase(i).GE.1 .AND. isteps(i).GT.istepmin) 
     &           isteps(i) = istepmin
         END DO
      ENDIF
c
c--Sort particles based on their individual timesteps and x
c

      CALL indexx(npart, llist, tempsort, iorig)

      DO i = 1, npart
         isort(iorig(i)) = i
      ENDDO
      DO k = 1, 5
         DO i = 1, npart
            tempsort(i) = xyzmh(k,iorig(i))
         END DO
         DO i = 1, npart
            xyzmh(k,i) = tempsort(i)
         END DO
      ENDDO
      DO k = 1, 4
         DO i = 1, npart
            tempsort(i) = vxyzu(k,iorig(i))
         END DO
         DO i = 1, npart
            vxyzu(k,i) = tempsort(i)
         END DO
      ENDDO
      DO i = 1, npart
         tempsort(i) = rho(iorig(i))
      END DO
      DO i = 1, npart
         rho(i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dgrav(iorig(i))
      END DO
      DO i = 1, npart
         dgrav(i) = tempsort(i)
      END DO
      DO k = 1, isizealphaMM
         DO i = 1, npart
            tempsort(i) = alphaMM(k,iorig(i))
         END DO
         DO i = 1, npart
            alphaMM(k,i) = tempsort(i)
         END DO
      ENDDO
      
      IF (encal.EQ.'r') THEN
         DO k = 1, 5
            DO i = 1, npart
               tempsort(i) = ekcle(k,iorig(i))
            END DO
            DO i = 1, npart
               ekcle(k,i) = tempsort(i)
            END DO
         ENDDO
      ENDIF
      IF (imhd.EQ.idim) THEN
         DO k = 1, 3
            DO i = 1, npart
               tempsort(i) = Bevolxyz(k,iorig(i))
            END DO
            DO i = 1, npart
               Bevolxyz(k,i) = tempsort(i)
            END DO
         ENDDO
      ENDIF

      DO i = 1, npart
         itempsort(i) = isteps(iorig(i))
      END DO
      DO i = 1, npart
         isteps(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = iphase(iorig(i))
      END DO
      DO i = 1, npart
         iphase(i) = itempsort(i)
      END DO
      DO i = 1, nptmass
         listpm(i) = isort(listpm(i))
      END DO
c
c--Zero torques
c
c      DO i = 1, idim
c         torqt(i) = 0.
c         torqg(i) = 0.
c         torqp(i) = 0.
c         torqv(i) = 0.
c         torqc(i) = 0.
c      END DO
c
c--Check units in file the same as in the code!
c
      IF (udisti.LT.0.99999*udist .OR. udisti.GT.1.00001*udist) THEN
         CALL error(where,1)
      ELSEIF (umassi.LT.0.99999*umass .OR.umassi.GT.1.00001*umass) THEN
         CALL error(where,2)
      ELSEIF (imhd.EQ.idim) THEN
         IF (umagfdi.LT.0.9999*umagfd 
     &      .OR.umagfdi.GT.1.00001*umagfd) THEN
            CALL error(where,4)
         ENDIF
      ENDIF
      IF (npart.GT.idim) THEN
         CALL error(where,3)
      ENDIF
c
c--Check that dtmax times are the same.  If not, modify isteps(i) as in mesop.f
c
ccc      GOTO 50

      IF (gt.NE.0.0 .AND. 
     &     (dtmaxdp.LT.0.9999*dtmax .OR. dtmaxdp.GT.1.0001*dtmax)) THEN
         ipower = INT(LOG10(dtmax/dtmaxdp)/LOG10(2.0))

         ifactor = 2**ABS(ipower)
         imaxnew = imaxstep/ifactor
         iminnew = 2*ifactor

         IF (ipower.LT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  isteps(j) = MIN(isteps(j), imaxnew)
                  isteps(j) = isteps(j)*ifactor
               ENDIF
            END DO
         ELSEIF (ipower.GT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  IF (isteps(j)/ifactor .LT. 2) CALL error(where, 4)
                  isteps(j) = isteps(j)/ifactor
               ENDIF
            END DO
         ENDIF
      ENDIF
c
c--Change reference frame
c
c      IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2.AND.gt.NE.0.0)) THEN
      IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2)) THEN
         CALL chanref(icall)
      ELSEIF (ifcor.GT.2) THEN
         ifcor = ifcor - 2
      ENDIF

#ifdef MPIDEBUG
      xxmin = 1.0E+30
      DO i = 1, npart
         xxnew = xyzmh(1,i)**2 + xyzmh(2,i)**2 + (xyzmh(3,i)-0.0)**2
         IF (xxnew.LT.xxmin) THEN
            xxmin = xxnew
            ipos1 = i
         ENDIF
      END DO
      print *,iproc,': CLOSEST z=+0.0 ',ipos1,xxmin
      xxmin = 1.0E+30
      DO i = 1, npart
         xxnew = xyzmh(1,i)**2 + xyzmh(2,i)**2 + (xyzmh(3,i)+0.0)**2
         IF (xxnew.LT.xxmin) THEN
            xxmin = xxnew
            ipos2 = i
         ENDIF
      END DO
      print *,iproc,': CLOSEST z=-0.0 ',ipos2,xxmin
#endif

      IF (itrace.EQ.'all') WRITE (*, 99002)
99002 FORMAT (' exit subroutine rdump')
      RETURN

 100  CONTINUE
      WRITE(*,*) 'END OF FILE REACHED IN RDUMP'
      ichkl = 1

      RETURN
      END


c============================================================


      SUBROUTINE split(npart,n1,icount,nhydro,isteps,iphase,xyzmh,
     &     vxyzu,rho,dgrav,Bevolxyz,Bxyz,gradhs,alphaMM,iunique)
c************************************************************
c                                                           *
c--Splits based on spatial coordinates (x, y, z) for 2, 4   *
c     or 8 processes. For more processes (10, 12, 14, 16,   *
c     32, 38, 48, 64, 74, 122, 128) it does a spatial split *
c     that tries to give equal numbers of particles from a  *
c     sphere to each process in wedges that all extend into *
c     the middle of the sphere.                             *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      INTEGER*1 iphase(idim)
      INTEGER*8 iunique(idim)
      REAL*4 rho(idim), dgrav(idim), gradhs(2,idim),
     &     alphaMM(isizealphaMM,idim)
      DIMENSION isteps(idim), xyzmh(5,idim), vxyzu(4,idim), 
     &     Bevolxyz(3,imhd), Bxyz(3,imhd)

      INTEGER classify_disc,classify_sphere,classify_cartesian

      DIMENSION numdiv(10)

      icountnew = icount
c
c--Type of domain decomposition
c
c--Disc-type domain decomposition
c
      IF (mpitype.EQ.'d') THEN
         DO i = icount+1, nhydro
            ibelong = classify_disc(i,xyzmh)
            IF (iproc.EQ.ibelong) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)
            ENDIF
         END DO
c
c--Sphere-type domain decomposition
c
      ELSEIF (mpitype.EQ.'s') THEN
         IF (numproc.EQ.2) THEN
            DO i = icount+1, nhydro
               IF (iproc.EQ.0 .AND. 
     &              xyzmh(1,i).GE.0.0) THEN
                  icountnew = icountnew + 1
                  isteps(icountnew) = isteps(i)
                  iphase(icountnew) = iphase(i)
                  DO j = 1, 5
                     xyzmh(j,icountnew) = xyzmh(j,i)
                  END DO
                  DO j = 1, 4
                     vxyzu(j,icountnew) = vxyzu(j,i)
                  END DO
                  rho(icountnew) = rho(i)
                  dgrav(icountnew) = dgrav(i)
                  IF (nlmax.EQ.1) THEN
                     DO j = 1, 2
                        gradhs(j,icountnew) = gradhs(j,i)
                     END DO
                  ENDIF
                  IF (imhd.EQ.idim) THEN
                     DO j = 1, 3
                        Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                     END DO
                     DO j = 1, 3
                        Bxyz(j,icountnew) = Bxyz(j,i)
                     END DO
                  ENDIF
                  DO j = 1, isizealphaMM
                     alphaMM(j,icountnew) = alphaMM(j,i)
                  END DO
                  iunique(icountnew) = iunique(i)

               ELSEIF (iproc.EQ.1 .AND. 
     &                 xyzmh(1,i).LT.0.0) THEN
                  icountnew = icountnew + 1
                  isteps(icountnew) = isteps(i)
                  iphase(icountnew) = iphase(i)
                  DO j = 1, 5
                     xyzmh(j,icountnew) = xyzmh(j,i)
                  END DO
                  DO j = 1, 4
                     vxyzu(j,icountnew) = vxyzu(j,i)
                  END DO
                  rho(icountnew) = rho(i)
                  dgrav(icountnew) = dgrav(i)
                  IF (nlmax.EQ.1) THEN
                     DO j = 1, 2
                        gradhs(j,icountnew) = gradhs(j,i)
                     END DO
                  ENDIF
                  IF (imhd.EQ.idim) THEN
                     DO j = 1, 3
                        Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                     END DO
                     DO j = 1, 3
                        Bxyz(j,icountnew) = Bxyz(j,i)
                     END DO
                  ENDIF
                  DO j = 1, isizealphaMM
                     alphaMM(j,icountnew) = alphaMM(j,i)
                  END DO
                  iunique(icountnew) = iunique(i)
                  
               ENDIF
            END DO
         ELSEIF (numproc.EQ.4) THEN
            DO i = icount+1, icount+nhydro
               IF (iproc.EQ.0 .AND. xyzmh(3,i).GE.0.0 .AND.
     &              xyzmh(1,i).GE.0.0) THEN
                  icountnew = icountnew + 1
                  isteps(icountnew) = isteps(i)
                  iphase(icountnew) = iphase(i)
                  DO j = 1, 5
                     xyzmh(j,icountnew) = xyzmh(j,i)
                  END DO
                  DO j = 1, 4
                     vxyzu(j,icountnew) = vxyzu(j,i)
                  END DO
                  rho(icountnew) = rho(i)
                  dgrav(icountnew) = dgrav(i)
                  IF (nlmax.EQ.1) THEN
                     DO j = 1, 2
                        gradhs(j,icountnew) = gradhs(j,i)
                     END DO
                  ENDIF
                  IF (imhd.EQ.idim) THEN
                     DO j = 1, 3
                        Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                     END DO
                     DO j = 1, 3
                        Bxyz(j,icountnew) = Bxyz(j,i)
                     END DO
                  ENDIF
                  DO j = 1, isizealphaMM
                     alphaMM(j,icountnew) = alphaMM(j,i)
                  END DO
                  iunique(icountnew) = iunique(i)
                  
               ELSEIF (iproc.EQ.1 .AND. xyzmh(3,i).GE.0.0 .AND.
     &                 xyzmh(1,i).LT.0.0) THEN
                  icountnew = icountnew + 1
                  isteps(icountnew) = isteps(i)
                  iphase(icountnew) = iphase(i)
                  DO j = 1, 5
                     xyzmh(j,icountnew) = xyzmh(j,i)
                  END DO
                  DO j = 1, 4
                     vxyzu(j,icountnew) = vxyzu(j,i)
                  END DO
                  rho(icountnew) = rho(i)
                  dgrav(icountnew) = dgrav(i)
                  IF (nlmax.EQ.1) THEN
                     DO j = 1, 2
                        gradhs(j,icountnew) = gradhs(j,i)
                     END DO
                  ENDIF
                  IF (imhd.EQ.idim) THEN
                     DO j = 1, 3
                        Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                     END DO
                     DO j = 1, 3
                        Bxyz(j,icountnew) = Bxyz(j,i)
                     END DO
                  ENDIF
                  DO j = 1, isizealphaMM
                     alphaMM(j,icountnew) = alphaMM(j,i)
                  END DO
                  iunique(icountnew) = iunique(i)
                  
               ELSEIF (iproc.EQ.2 .AND. xyzmh(3,i).LT.0.0 .AND.
     &              xyzmh(1,i).GE.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.3 .AND. xyzmh(3,i).LT.0.0 .AND.
     &              xyzmh(1,i).LT.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ENDIF
         END DO
      ELSEIF (numproc.EQ.8) THEN
         DO i = icount+1, icount+nhydro
            IF (iproc.EQ.0 .AND. xyzmh(3,i).GE.0.0 .AND.
     &           xyzmh(1,i).GE.0.0 .AND. xyzmh(2,i).GE.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.1 .AND. xyzmh(3,i).GE.0.0 .AND.
     &              xyzmh(1,i).LT.0.0 .AND. xyzmh(2,i).GE.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.2 .AND. xyzmh(3,i).LT.0.0 .AND.
     &              xyzmh(1,i).GE.0.0 .AND. xyzmh(2,i).GE.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.3 .AND. xyzmh(3,i).LT.0.0 .AND.
     &              xyzmh(1,i).LT.0.0 .AND. xyzmh(2,i).GE.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.4 .AND. xyzmh(3,i).GE.0.0 .AND.
     &           xyzmh(1,i).GE.0.0 .AND. xyzmh(2,i).LT.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.5 .AND. xyzmh(3,i).GE.0.0 .AND.
     &              xyzmh(1,i).LT.0.0 .AND. xyzmh(2,i).LT.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.6 .AND. xyzmh(3,i).LT.0.0 .AND.
     &              xyzmh(1,i).GE.0.0 .AND. xyzmh(2,i).LT.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ELSEIF (iproc.EQ.7 .AND. xyzmh(3,i).LT.0.0 .AND.
     &              xyzmh(1,i).LT.0.0 .AND. xyzmh(2,i).LT.0.0) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ENDIF
         END DO
      ELSEIF (numproc.EQ.10 .OR. 
     &        numproc.EQ.12 .OR. numproc.EQ.14 .OR. numproc.EQ.16 
     &        .OR. numproc.EQ.32 .OR. numproc.EQ.38 .OR. numproc.EQ.48
     &        .OR. numproc.EQ.64 .OR. numproc.EQ.74 
     &        .OR. numproc.EQ.122.OR. numproc.EQ.128) THEN
         numdiv(1) = 0
         IF (numproc.EQ.10) THEN
            xpower = 1.0
            ndivisor = 2
            numdiv(2) = 4
         ELSEIF (numproc.EQ.12) THEN
            xpower = 1.0
            ndivisor = 2
            numdiv(2) = 5
         ELSEIF (numproc.EQ.14) THEN
            xpower = 0.65
            ndivisor = 2
            numdiv(2) = 6
c         ELSEIF (numproc.EQ.8) THEN
c            xpower = 1.0
c            ndivisor = 2
         ELSEIF (numproc.EQ.16) THEN
            xpower = 0.61
            ndivisor = 2
            numdiv(2) = 7
         ELSEIF (numproc.EQ.32) THEN
            xpower = 0.74
            ndivisor = 3
            numdiv(2) = 5
            numdiv(3) = 10
         ELSEIF (numproc.EQ.38) THEN
            xpower = 0.70
            ndivisor = 3
            numdiv(2) = 6
            numdiv(3) = 12
         ELSEIF (numproc.EQ.48) THEN
            xpower = 0.65
            ndivisor = 3
            numdiv(2) = 8
            numdiv(3) = 15
         ELSEIF (numproc.EQ.64) THEN
            xpower = 0.76
            ndivisor = 4
            numdiv(2) = 5
            numdiv(3) = 10
            numdiv(4) = 16
         ELSEIF (numproc.EQ.74) THEN
            xpower = 0.72
            ndivisor = 4
            numdiv(2) = 6
            numdiv(3) = 12
            numdiv(4) = 18
         ELSEIF (numproc.EQ.122) THEN
            xpower = 0.75
            ndivisor = 5
c--Given by 6*(j-1)
            numdiv(2) = 6
            numdiv(3) = 12
            numdiv(4) = 18
            numdiv(5) = 24
         ELSEIF (numproc.EQ.128) THEN
            xpower = 0.74
            ndivisor = 5
            numdiv(2) = 6
            numdiv(3) = 12
            numdiv(4) = 19
            numdiv(5) = 26
         ENDIF
         iprocdiv2 = numproc/2
         pilocal = ATAN2(0.0,-1.0)
c         print *,pilocal,ATAN(0.2/1.0)

         DO i = icount+1, icount+nhydro
            xlatitude = ABS(ABS(ATAN(xyzmh(3,i)/SQRT(xyzmh(1,i)**2 +
     &           xyzmh(2,i)**2))) - pilocal/2.0)
            IF (xlatitude.LT.0.0) xlatitude = 0.0
            xlongitude = ATAN2(xyzmh(2,i),xyzmh(1,i)) + pilocal

            jvalue = INT(ndivisor*(xlatitude/(pilocal/2.0))**xpower)+1
            IF (jvalue.GT.ndivisor) jvalue = ndivisor
            
            j2value = numdiv(jvalue)
            IF (jvalue.GT.1) THEN
               lvalue = INT(numdiv(jvalue)*xlongitude/(2.0*pilocal)) + 1
               IF (lvalue.GT.numdiv(jvalue)) lvalue = numdiv(jvalue)
            ELSE
               lvalue = 0
            ENDIF

            iprocvalue = 0
            DO ii = 1, jvalue-1
               iprocvalue = iprocvalue + numdiv(ii)
            END DO
            iprocvalue = iprocvalue + lvalue

c            IF (iproc.EQ.0 .AND. MOD(i,1000).EQ.0) 
c     &       print *,iproc,jvalue,numdiv(jvalue),lvalue,iprocvalue

            IF (xyzmh(3,i).GE.0.0) iprocvalue = iprocvalue + iprocdiv2

            IF (iproc.EQ.iprocvalue) THEN
               icountnew = icountnew + 1
               isteps(icountnew) = isteps(i)
               iphase(icountnew) = iphase(i)
               DO j = 1, 5
                  xyzmh(j,icountnew) = xyzmh(j,i)
               END DO
               DO j = 1, 4
                  vxyzu(j,icountnew) = vxyzu(j,i)
               END DO
               rho(icountnew) = rho(i)
               dgrav(icountnew) = dgrav(i)
               IF (nlmax.EQ.1) THEN
                  DO j = 1, 2
                     gradhs(j,icountnew) = gradhs(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN               
                  DO j = 1, 3
                     Bevolxyz(j,icountnew) = Bevolxyz(j,i)
                  END DO
                  DO j = 1, 3
                     Bxyz(j,icountnew) = Bxyz(j,i)
                  END DO
               ENDIF
               DO j = 1, isizealphaMM
                  alphaMM(j,icountnew) = alphaMM(j,i)
               END DO
               iunique(icountnew) = iunique(i)

            ENDIF
         END DO
      ELSE
         PRINT *,' Currently invalid number of processes ',numproc
         CALL quit
      ENDIF

      ELSE
         PRINT *,' Currently invalid MPI domain decompostion ',numproc,
     &        mpitype
         CALL quit
      ENDIF
      npart = icountnew
      n1 = icountnew
      icount = icountnew

      RETURN
      END
