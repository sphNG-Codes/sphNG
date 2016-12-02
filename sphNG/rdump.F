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
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

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
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/xforce'      
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/savernd'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'

      CHARACTER*100 fileident
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      CHARACTER*16 tagsreal(128),tagexp(10),tagi
      DIMENSION nums1(8),nums2(8),nums3(8),nums4(8)
      PARAMETER (idimhead = 30)
      DIMENSION rheader(idimhead)
      INTEGER iversion, maxversion
      LOGICAL tagged, isequential
c
c--Is an invalid cleaning wave array a fatal error?
c    Do not change this value unless you know
c    *exactly what you are doing*.
      LOGICAL psierror
      PARAMETER (psierror=.TRUE.)
c     
c--Read
c
      IF (itrace.EQ.'all') WRITE (*, 99001)
99001 FORMAT (' entry subroutine rdump')

      ichkl = 0
      icount = 0
      icountinactive = 0
      icountsink = 0
      nhydrocumulative = 0

      iplanetesimals = 0
      irotpot = 0
      idragscheme = 100
c
c--Dump file
c-------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
      int1alt = 060769
      int2alt = 060878
c
c--format version readable by this routine
c
      maxversion = 1
      tagged = .FALSE.
c
c--Read ouput file
c
c--This section (and another section below) force MPI processes to read
c     the dump file sequentially (for file systems that do not cope well
c     with multiple simultaneous reads).
c
      isequential = .FALSE.
      IF (isequential) THEN
#ifdef MPI
         nowgo = 0
#ifdef MPIDEBUG
         print *,iproc,': About to read file'
#endif
         IF (iproc.GT.0) THEN
            CALL MPI_RECV(nowgo,1,MPI_INTEGER,iproc-1,99,MPI_COMM_WORLD,
     &           istatus,ierr)
#ifdef MPIDEBUG
            print *,iproc,': Reading file'
#endif
         ENDIF
#endif
      ENDIF
c
c--End of sequential reads section
c
      print *,'reading header'

      READ (idisk1, END=100) int1i,r1i,int2i,iversion,int3i

      print *,'have read header'

      IF (int1i.NE.int1 .and. int1i.NE.int1alt) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         CALL quit(0)
      ENDIF
      IF (int2i.NE.int2 .and. int2i.NE.int2alt) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default real size wrong'
         CALL quit(0)
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default integer size wrong'
         CALL quit(0)
      ENDIF
c
c--handle version numbers
c  version 0 had iversion=690706
c
      IF (iversion.EQ.int1) THEN
         iversion = 0
      ENDIF
      IF (iversion.GT.maxversion) THEN
         WRITE (*,*) 'ERROR 4 in rdump: format version is ',iversion,
     &        ' but this version of sphNG can only read ',maxversion
      ENDIF

      READ (idisk1, END=100) fileident
      IF (fileident(2:2).EQ.'T') THEN
         tagged = .true.
      ENDIF
c
c--Single values
c
c--Default int
      READ (idisk1, END=100) number
      print *,'number ints ',number
      IF (number.LT.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump: not enough default ints'
         CALL quit(0)
      ENDIF
      IF (tagged) READ (idisk1, END=100) ! skip tags
      IF (number.EQ.9+NTAB) THEN
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill,
     &        nblocks,iyr,idum,(iv(i), i=1,NTAB)
         print *,'USING 9+NTAB'
      ELSEIF (number.EQ.10+NTAB) THEN
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill,
     &        nblocks,iyr,idum,(iv(i), i=1,NTAB), iplanetesimals
         print *,'USING 10+NTAB'
      ELSEIF (number.EQ.11+NTAB) THEN
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill,
     &        nblocks,iyr,idum,(iv(i), i=1,NTAB), iplanetesimals,
     &        irotpot
         print *,'USING 11+NTAB', iplanetesimals, irotpot
      ELSEIF (number.EQ.12+NTAB) THEN
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill,
     &        nblocks,iyr,idum,(iv(i), i=1,NTAB), iplanetesimals,
     &        irotpot, idragscheme
         print *,'USING 12+NTAB', irotpot, idragscheme
      ELSEIF (number.EQ.7) THEN
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill,
     &        nblocks
      ELSE
         READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill
         nblocks = 1
      ENDIF
      IF (npart.GT.idim .AND. 
     &     (nblocks.EQ.1 
#ifndef MPI
     &     .OR. nblocks.GT.1)
#else
     &     )
#endif
     &     ) THEN
         WRITE (*,*) 'ERROR in rdump: npart>idim'
         CALL quit(0)
      ENDIF
c--int*1, int*2, int*4, int*8
      DO i = 1, 3
         READ (idisk1, END=100) number
         print *,'number int1,2,4 ',number
      END DO
      READ (idisk1, END=100) number
      print *,'number int8 ',number
      IF (number.EQ.1) THEN
         IF (tagged) READ (idisk1, END=100) ! skip tags
         READ (idisk1, END=100) iuniquemax
      ELSE
         iuniquemax = npart
      ENDIF
c--Default real
      READ (idisk1, END=100) number
      print *,'number real ',number
      IF (number.LT.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump: not enough default reals'
         CALL quit(0)
      ENDIF

      IF (tagged) THEN
c--read the actual tags from tagged dump files
         READ (idisk1, END=100) tagsreal(1:number)
      ELSE
c--fake the tags for older formats
         tagsreal(1:30) = 
     &  (/'gt              ','dtmax           ','gamma           ',
     &    'rhozero         ','RK2             ','escap           ',
     &    'tkin            ','tgrav           ','tterm           ',
     &    'anglostx        ','anglosty        ','anglostz        ',
     &    'specang         ','ptmassin        ','tmag            ',
     &    'Bextx           ','Bexty           ','Bextz           ',
     &    'hzero           ','uzero_n2        ','hmass           ',
     &    'gapfac          ','                ','sdprof          ',
     &    'rorbit_orig     ','min_rplan       ','max_rplan       ',
     &    'planetesimalmass','coremass_orig   ','coremass        '/)
      ENDIF
      READ (idisk1, END=100) (rheader(i),i=1,min(number,idimhead))
c
c--extract variables from the header
c     
      CALL extract('gt',gt,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) CALL quit(0)
      CALL extract('dtmax',dtmaxdp,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) CALL quit(0)
      CALL extract('gamma',gamma,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) CALL quit(0)
      CALL extract('rhozero',rhozero,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) CALL quit(0)
      CALL extract('RK2',RK2,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) CALL quit(0)
      CALL extract('escap',escap,rheader,tagsreal,number,ierr)
      CALL extract('tkin',tkin,rheader,tagsreal,number,ierr)
      CALL extract('tgrav',tgrav,rheader,tagsreal,number,ierr)
      CALL extract('tterm',tterm,rheader,tagsreal,number,ierr)
      CALL extract('anglostx',anglostx,rheader,tagsreal,number,ierr)
      CALL extract('anglosty',anglosty,rheader,tagsreal,number,ierr)
      CALL extract('anglostz',anglostz,rheader,tagsreal,number,ierr)
      CALL extract('specang',specang,rheader,tagsreal,number,ierr)
      CALL extract('ptmassin',ptmassin,rheader,tagsreal,number,ierr)
      IF (imhd.EQ.idim) THEN
         CALL extract('tmag',tmag,rheader,tagsreal,number,ierr)
         CALL extract('Bextx',Bextx,rheader,tagsreal,number,ierr)
         CALL extract('Bexty',Bexty,rheader,tagsreal,number,ierr)
         CALL extract('Bextz',Bextz,rheader,tagsreal,number,ierr)
         WRITE(*,*) 'External field found, Bext = ',Bextx,Bexty,Bextz
      ENDIF
      CALL extract('hzero',hzero,rheader,tagsreal,number,ierr)
      IF (iexf.EQ.10 .and. ierr.NE.0) CALL quit(0)      
      CALL extract('uzero_n2',uzero_n2,rheader,tagsreal,number,ierr)
      IF (n2.GT.0) THEN
         WRITE(*,*) ' read u for surrounding medium = ',uzero_n2
      ENDIF
      CALL extract('hmass',hmass,rheader,tagsreal,number,ierr)
      CALL extract('gapfac',gapfac,rheader,tagsreal,number,ierr)
      CALL extract('sdprof',sdprof,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) THEN
         sdprof = -0.5
         print *, 'Surface density pre vary (goes as r^-0.5)'
      ENDIF
      CALL extract('rorbit_orig',rorbit_orig,
     &             rheader,tagsreal,number,ierr)
      IF (imigrate.GT.0) rorbitmax = (rorbitmax - rorbit_orig)/pmrate
      CALL extract('min_rplan',min_rplan,rheader,tagsreal,number,ierr)
      CALL extract('max_rplan',max_rplan,rheader,tagsreal,number,ierr)
      IF (ierr.NE.0) max_rplan = 1.0E6
      print *, 'Setting planetesimal boundaries; check values.'
      print *, 'r_min = ', min_rplan
      print *, 'r_max = ', max_rplan
      IF (.NOT.tagged) THEN
         print *, 'DO NOT USE CODE WITH DUMPS MADE FROM APRIL 2011'
         print *, ' AND BEFORE 16/05/2011. For these, reals 26, 27 are'
         print *, ' planetesimal radius and density respectively.'
      ENDIF
      CALL extract('planetesimalmass',planetesimalmass,
     &             rheader,tagsreal,number,ierr)
      print *, 'Planetesimal mass ', planetesimalmass
      CALL extract('coremass_orig',coremass_orig,
     &             rheader,tagsreal,number,ierr)
      print *, 'Core mass orig ', coremass_orig
      CALL extract('coremass',coremass,rheader,tagsreal,number,ierr)
      print *, 'Core mass running record ', coremass

c--real*4
      READ (idisk1, END=100) number
c--real*8
      READ (idisk1, END=100) number
      IF (number.LT.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump: nreal8 too small header'
         CALL quit(0)
      ENDIF
      IF (tagged) READ (idisk1, END=100) ! skip tags
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
      IF (iradtrans.EQ.idim) numberexp = 3
      IF (imhd.EQ.idim) numberexp = 4

      IF (number.LT.2 .OR. number/nblocks.LT.2 .OR. 
     &     numberexp.EQ.4 .AND.number/nblocks.NE.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump ', number, numberexp
         CALL quit(0)
#ifdef MPI
      ELSEIF (nblocks.NE.1 .AND. nblocks.NE.numproc) THEN
         WRITE (*,*) 'ERROR in rdump: nblocks.NE.numproc ',
     &        nblocks,numberexp,numproc
         WRITE (*,*) 'Cannot start up job with different number of MPI',
     &        ' processes'
         CALL quit(0)
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
         CALL quit(0)
      ENDIF
#ifdef MPI
      IF (nblocks.EQ.1 .OR. nblocks.NE.1 .AND. iproc.EQ.iblock-1) THEN
#endif
      nhydro = min(number8,npart)
      IF (iskip.EQ.1 .AND. icount + nhydro.GT.idim) THEN
         WRITE (*,*) 'ERROR - idim exceeded'
         CALL quit(0)
      ENDIF
#ifdef MPI
      ELSE
         nhydrocumulative = nhydrocumulative + number8
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
         CALL quit(0)
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
         print *,'nums3 ',(nums3(i),i=1,8)
         IF (iradtrans.EQ.idim .AND.
     &      (number8.GT.iradtrans)) THEN
            WRITE (*,*) 'ERROR 9 in rdump: iradtrans wrong ',number8,
     &           iradtrans,npart
            CALL quit(0)
         ENDIF
         nradtrans = min(number8,npart)
      ELSE
         nradtrans = 0
      ENDIF
c
c--Read array type 4 header
c
      IF (numberarray.GE.4) THEN
         READ (idisk1, END=100) number8, (nums4(i), i=1,8)
         print *,'nums4 ',(nums4(i),i=1,8)
         IF (imhd.EQ.idim.AND. 
     &      (number8.GT.imhd)) THEN
            WRITE (*,*) 'ERROR 10 in rdump: imhd wrong ',number8,
     &           imhd,npart
            CALL quit(0)
         ENDIF
         nmhd = min(number8,npart)
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
               IF (tagged) READ (idisk1, END=100) ! skip tags
               READ (idisk1, END=100)
            END DO
         END DO
         DO j = 1, 8
            DO k = 1, nums2(j)
               IF (tagged) READ (idisk1, END=100) ! skip tags
               READ (idisk1, END=100)
            END DO
         END DO
         IF (numberarray.GE.3) THEN
            DO j = 1, 8
               DO k = 1, nums3(j)
                  IF (tagged) READ (idisk1, END=100) ! skip tags
                  READ (idisk1, END=100)
               END DO
            END DO
            IF (numberarray.GE.4) THEN
               DO j = 1, 8
                  DO k = 1, nums4(j)
                     IF (tagged) READ (idisk1, END=100) ! skip tags
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
         IF (tagged) THEN
           READ (idisk1, END=100) tagi
           CALL checktag('isteps', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) (isteps(i), i=icount+1, icount+nhydro)
#ifdef MPIDEBUG
         print *,iproc,': read isteps ',isteps(1)
#endif
         IF (nums1(1).GE.2) THEN
            IF (tagged) THEN 
              READ (idisk1, END=100) tagi
              CALL checktag('listinactive', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) nlistinactive,
     &           listinactive(icountinactive+1:icountinactive+
     &           nhydro-1)
c
c--If reading multiple blocks onto single process, need to re-index the
c     inactive particles
c
            IF (icount.GT.0) THEN
               DO i = icountinactive+1,icountinactive+nlistinactive
                  listinactive(i) = listinactive(i) + icount
               END DO
            ENDIF
            print *,'READ nlistinactive ',nlistinactive
         ENDIF

c--int*1
         IF (tagged) THEN
            READ (idisk1, END=100) tagi
            CALL checktag('iphase', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) (iphase(i), i=icount+1, icount+nhydro)
         print *,'iphases ',iphase(1),iphase(nhydro),nums1(5),icount

         DO i = icountinactive+1,icountinactive+nlistinactive
            IF (iphase(listinactive(i)).NE.-1) THEN
               print *,'ERROR - Read active particle as inactive ',
     &              i,icountinactive,listinactive(i),icount,nhydro,
     &              iphase(listinactive(i))
               CALL quit(0)
            ENDIF
         END DO

c--int*2

c--int*4

c--int*8
         IF (nums1(5).GE.1) THEN
            IF (tagged) THEN
              READ (idisk1, END=100) tagi
              CALL checktag('iunique', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) (iunique(i), i=icount+1,
     &           icount+nhydro)
         ENDIF
#ifdef MPIDEBUG
         print *,iproc,': read iunique ',iunique(icount+1)
#endif

c         DO ijk = icount+1, icount+nhydro
c            WRITE (45,*) iunique(ijk), iphase(ijk)
c         END DO

         IF (nums1(5).LT.1 .OR. iunique(icount+1).EQ.0) THEN
            print *,'Read iunique ',(iunique(i), i=icount+1,icount+10)
            DO i=icount+1,icount+nhydro
               iunique(i) = nhydrocumulative + i
            END DO
         ENDIF
c--Default real
         tagexp(1:5) = (/'x','y','z','m','h'/)
         DO j = 1, 5
            IF (tagged) THEN
               READ (idisk1, END=100) tagi
               CALL checktag(tagexp(j), tagi, ierr)          
            ENDIF
            READ (idisk1, END=100) xyzmh(j,icount+1:icount+nhydro)
         END DO

         tagexp(1:4) = (/'vx','vy','vz','u '/)
         DO j = 1, 4
            IF (tagged) THEN
              READ (idisk1, END=100) tagi
              CALL checktag(tagexp(j), tagi, ierr)
            ENDIF
            READ (idisk1, END=100) vxyzu(j,icount+1:icount+nhydro)
         END DO

#ifdef MPIDEBUG
         print *,iproc,': Read xyzmh, vxyzu'
#endif
         iwanted = 9
         IF (iener.EQ.4) THEN
            iwanted = iwanted + 5
 
            IF (nums1(6).LT.iwanted) THEN
               WRITE (*,*) 'ERROR - nums1(6).LT.iwanted ',
     &              nums1(6),iwanted
               CALL quit(0)
            ENDIF

            IF (tagged) THEN
              READ (idisk1, END=100) tagi
              CALL checktag('h2ratio', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) h2ratio(icount+1:icount+nhydro)

            IF (tagged) THEN
              READ (idisk1, END=100) tagi
              CALL checktag('abHIq', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) abHIq(icount+1:icount+nhydro)

            IF (tagged) THEN
              READ (idisk1, END=100) tagi
              CALL checktag('abhpq', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) abhpq(icount+1:icount+nhydro)

            IF (tagged) THEN 
              READ (idisk1, END=100) tagi
              CALL checktag('abeq', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) abeq(icount+1:icount+nhydro)

            IF (tagged) THEN 
              READ (idisk1, END=100) tagi
              CALL checktag('abco', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) abco(icount+1:icount+nhydro)
          ENDIF
c--skip unnecessary reals
         IF (nums1(6).GT.iwanted) THEN
            DO j=1,nums1(6)-iwanted
               IF (tagged) READ (idisk1, END=100)
               READ (idisk1, END=100)
            END DO
         ENDIF
c--real*4
         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('rho', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) rho(icount+1:icount+nhydro)
         iread = 1
         IF (nlmax.EQ.1) THEN
            IF (nums1(7).GE.2) THEN
               IF (tagged) THEN
                 READ (idisk1, END=100) tagi
                 CALL checktag('gradh', tagi, ierr)
               ENDIF
               READ (idisk1, END=100) gradhs(1,icount+1:icount+nhydro)
               iread = iread + 1
            ENDIF
            IF (nums1(7).GE.3) THEN
               IF (tagged) THEN 
                 READ (idisk1, END=100) tagi
                 CALL checktag('gradsoft', tagi, ierr)
               ENDIF
               READ (idisk1, END=100) gradhs(2,icount+1:icount+nhydro)
               iread = iread + 1
            ENDIF
            DO j = 1, idim
               dgrav(j) = 0.
            END DO
         ELSE
            iread = 2
            IF (tagged) THEN 
              READ (idisk1, END=100) tagi
              CALL checktag('dgrav', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) dgrav(icount+1:icount+nhydro)
            IF (gt.EQ.0.0) THEN
               DO j = icount+1, icount+nhydro
                  dgrav(j) = 0.
               END DO
            ENDIF
         ENDIF
         IF (nums1(7).GT.iread) THEN
            iread = iread + 1
            IF (tagged) THEN 
              READ (idisk1, END=100) tagi
              CALL checktag('alpha', tagi, ierr)
            ENDIF
            READ (idisk1, END=100) alphaMM(1,icount+1:icount+nhydro)
            IF (gt.NE.0.0 .AND. ifsvi.GE.6) THEN
               DO i = icount+1, icount+nhydro
                  IF (alphaMM(1,i).LT.alphamin(1) .OR. 
     &                 alphaMM(1,i).GT.alphamax(1)) THEN
                     alphaMM(1,i) = alphamin(1)
c                     WRITE (*,*) 'ERROR in rdump - alphaMM(1,i).LE.0 ',
c     &                    i,(alphaMM(1,j),j=1,10)
c                     CALL quit(0)
                  ENDIF
               END DO
            ENDIF
         ENDIF
c--skip unnecessary real*4's
         IF (nums1(7).GT.iread) THEN
            DO j=1,nums1(7)-iread
               IF (tagged) READ (idisk1, END=100) ! skip tags
               READ (idisk1, END=100)
            ENDDO
         ENDIF
c--real*8

#ifdef MPIDEBUG
         print *,iproc,': Read hydro ',nptmass
#endif
c
c--Read array type 2 arrays
c
c--Default int
         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('listpm', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) listpm(icountsink+1:icountsink+nptmass)

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('spinx', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) spinx(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('spiny', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) spiny(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('spinz', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) spinz(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN
           READ (idisk1, END=100) tagi
           CALL checktag('angaddx', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) angaddx(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('angaddy', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) angaddy(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('angaddz', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) angaddz(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('spinadx', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) spinadx(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN
           READ (idisk1, END=100) tagi
           CALL checktag('spinady', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) spinady(icountsink+1:icountsink+nptmass)

         IF (tagged) THEN 
           READ (idisk1, END=100) tagi
           CALL checktag('spinadz', tagi, ierr)
         ENDIF
         READ (idisk1, END=100) spinadz(icountsink+1:icountsink+nptmass)
c--real*4

c--real*8
         DO i = 1, nums2(8)
            IF (tagged) READ (idisk1, END=100) ! skip tags
            READ (idisk1, END=100)
         END DO

         IF (numberarray.GE.3 .AND. nradtrans.EQ.nhydro 
     &        .AND. iradtrans.GE.nradtrans) THEN

c
c--Array length 3 arrays
c      
c--Default int

c--int*1

c--int*2
            DO i = 1, nums3(3)
               IF (tagged) READ (idisk1, END=100) ! skip tags
               READ (idisk1, END=100)
            END DO
c--int*4

c--int*8

c--Default real
            tagexp(1:5) = (/'Erad           ','opacity        ',
     &                      'Cv             ','lambda         ',
     &                      'EddingtonFactor'/)
            DO j = 1, 5
               IF (tagged) THEN 
                 READ (idisk1, END=100) tagi
                 CALL checktag(tagexp(j), tagi, ierr)
               ENDIF
               READ (idisk1, END=100) ekcle(j,icount+1:icount+nradtrans)
            END DO
            ireadekcle = 5

            IF (idustRT.GT.0 .AND. nums3(6).GE.7) THEN
               tagexp(1:2)=(/'DustTemperature','DustOpacity    '/)
               DO j = 1, 2
                  IF (tagged) THEN
                    READ (idisk1, END=100) tagi
                    CALL checktag(tagexp(j), tagi, ierr)
                  ENDIF
                  READ (idisk1, END=100) 
     &                 dust_tk(j,icount+1:icount+nradtrans)
               END DO
               ireadekcle = ireadekcle + 2
            ENDIF

            DO j = 1, nums3(6) - ireadekcle
               IF (tagged) READ (idisk1, END=100) ! skip tags
               READ (idisk1, END=100)
            END DO
c--real*4
            IF (idustRT.GT.0 .AND. nums3(7).GE.
     &           nchemistry + nheatingISR + 2) THEN
                  tagexp(1:3) =(/'Cplus  ','AtomicC','CO     '/)
                  DO j = 1, nchemistry
                     IF (tagged) THEN
                       READ (idisk1, END=100) tagi
                       CALL checktag(tagexp(j), tagi, ierr)
                     ENDIF
                     READ (idisk1, END=100)
     &                    chemistry(j,icount+1:icount+nradtrans)
                  END DO

                  tagexp(1:4)=(/'HeatingISR     ','Gphotoelectric ',
     &                          'Extinction     ','H2SelfShielding'/)
                  DO j = 1, nheatingISR
                     IF (tagged) THEN 
                       READ (idisk1, END=100) tagi
                       CALL checktag(tagexp(j), tagi, ierr)
                     ENDIF
                     READ (idisk1, END=100) 
     &                    heatingISR(j,icount+1:icount+nradtrans)
                  END DO

                  IF (tagged) THEN 
                    READ (idisk1, END=100) tagi
                    CALL checktag('H2Frac', tagi, ierr)
                  ENDIF
                  READ (idisk1, END=100) 
     &                 h2frac(icount+1:icount+nradtrans)

                  IF (tagged) THEN 
                    READ (idisk1, END=100) tagi
                    CALL checktag('dH2dt', tagi, ierr)
                  ENDIF
                  READ (idisk1, END=100) 
     &                 dh2dt(1,icount+1:icount+nradtrans)

                  DO j = 1, nums3(7) - nchemistry - nheatingISR - 2
                     READ (idisk1, END=100)
                  END DO
            ELSE
               DO j = 1, nums3(7)
                  IF (tagged) READ (idisk1, END=100) ! skip tags
                  READ (idisk1, END=100)
               END DO
            ENDIF
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
               CALL quit(0)
            ENDIF
            tagexp(1:3) = (/'Bx','By','Bz'/)
            DO j = 1, 3
               IF (tagged) THEN 
                 READ (idisk1, END=100) tagi
                 CALL checktag(tagexp(j), tagi, ierr)
               ENDIF
               READ (idisk1, END=100) Bxyz(j,icount+1:icount+nmhd)
            END DO

c
c--Uncomment to modify input B field
c
c               DO i = icount+1,
c     &              icount+nmhd
c                  Bxyz(1,i) = 0.
c                  Bxyz(2,i) = 0.
c                  Bxyz(3,i) = 10.0 + xyzmh(3,i)
c               END DO

            IF (imhdevol.EQ.4
     &          .AND.varmhd.NE.'eulr'.AND.varmhd.NE.'vecp') THEN
               IF (nums4(6).LT.4) THEN
                  WRITE(*,*) 'ERROR: imhdevol.EQ.4 but nums4(6) < 4'
                  CALL quit(0)
               ELSE
                  IF (tagged) THEN 
                    READ (idisk1, END=100) tagi
                    CALL checktag('psi/c_h', tagi, ierr)
c Dump files from before r1799 cannot be used.
                    IF(ierr.NE.0.AND.psierror) THEN
                      WRITE(*,99981)
                      CALL quit(0)
                    ENDIF
                  ENDIF
99981             FORMAT(' Dump files from before r1799, when psi was',
     &                   ' changed to psi/c_h are incompatible.') 

                  READ (idisk1,END=100) Bevolxyz(4,icount+1:icount+nmhd)
               ENDIF
            ENDIF
#ifdef MPIDEBUG
         print *,iproc,': Read Bxyz'
#endif
c
c--read Euler potentials from dump if necessary
c
            IF (varmhd.EQ.'eulr') THEN
               IF (nums4(6).GE.5) THEN
                  tagexp(1:2)=(/'EulerAlpha','EulerBeta '/)
                  DO j = 1, 2
                     IF (tagged) THEN 
                       READ (idisk1, END=100) tagi
                       CALL checktag(tagexp(j), tagi, ierr)
                       IF (ierr.NE.0) THEN
c Fatal error on tag mismatch here 
                         WRITE(*,99982)
                         WRITE(*,99983)
                         CALL quit(0)
                       ENDIF
                     ENDIF
                     READ (idisk1, END=100) 
     &                    Bevolxyz(j,icount+1:icount+nmhd)
                  END DO
#ifdef MPIDEBUG
         print *,iproc,': Read Euler Bevolxyz'
#endif
               ELSE
                  WRITE(*,99982)
                  WRITE(*,99983)
                  CALL quit(0)
               ENDIF
99982          FORMAT(' ERROR: Cannot start Euler potentials run')
99983          FORMAT('        from non-Euler potentials dump')

            ELSEIF (varmhd.EQ.'vecp') THEN
               IF (nums4(6).GE.6) THEN
                  tagexp(1:3) = (/'Ax','Ay','Az'/)
                  DO j = 1, 3
                     IF (tagged) THEN 
                       READ (idisk1, END=100) tagi
                       CALL checktag(tagexp(j), tagi, ierr) 
                       IF (ierr.NE.0) THEN
c Fatal error on tag mismatch here
                         WRITE(*,99984)
                         WRITE(*,99985)
                         CALL quit(0)
                       ENDIF
                     ENDIF
                     READ (idisk1, END=100) 
     &                    Bevolxyz(j,icount+1:icount+nmhd)
                  END DO
#ifdef MPIDEBUG
         print *,iproc,': Read vecp Bevolxyz'
#endif
               ELSE
                  WRITE(*,99984)
                  WRITE(*,99985)
                  CALL quit(0)
               ENDIF
99984          FORMAT(' ERROR: Cannot start vec potentials run')
99985          FORMAT('        from non-vec potentials dump')

            ELSEIF (varmhd.EQ.'Brho') THEN
c
c--convert from B to B/rho for evolution
c
               DO i=icount+1,icount+nmhd
                  IF (iphase(i).EQ.0) THEN
                     IF (rho(i).LE.0.) THEN
                  WRITE(*,*) 'ERROR: rho <= 0 in rdump, evolving B/rho'
                        CALL quit(0)
                     ENDIF
                     DO j = 1,3
                        Bevolxyz(j,i) = Bxyz(j,i)/rho(i)
                     ENDDO
                  ENDIF
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
               WRITE (*,*) 'unknown MHD variable in rdump'
               DO j = imhdevol+1, MIN(nums4(6),6)
                  IF (tagged) READ (idisk1, END=100) ! skip tags
                  READ (idisk1,END=100) Bevolxyz(j,icount+1:icount+nmhd)
               END DO
               IF (iskip.EQ.0) CALL quit(0)
            ENDIF
c--real*4
            iread = 0
            IF (nums4(7).GT.4.AND.ifsvi.GE.6) THEN

c--skip current/div B
               DO j = 1, 4
                  iread = iread + 1
                  IF (tagged) READ (idisk1, END=100) ! skip tags
                  READ (idisk1, END=100)
               END DO

c--read artificial resistivity parameter
               IF (tagged) THEN 
                 READ (idisk1, END=100) tagi
                 CALL checktag('alphaB', tagi, ierr)
               ENDIF
               READ (idisk1, END=100) alphaMM(isizealphaMM, 
     &              icount+1:icount+nmhd)
               iread = iread + 1
            ELSE
               IF (ifsvi.GE.6) THEN
           WRITE(*,*) 'WARNING: resistivity parameter not found in dump'
               ENDIF
               WRITE(*,*) '=> setting alphaB=',alphamin(2),' everywhere'
               DO i=icount+1,icount+nmhd
                  alphaMM(2,i) = alphamin(2)
               ENDDO
            ENDIF

            DO i = iread + 1, nums4(7)
               IF (tagged) READ (idisk1, END=100) ! skip tags
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
               CALL split(npart,n1,nptmass,icount,nhydro,
     &              nlistinactive,isteps,iphase,
     &              xyzmh,vxyzu,rho,dgrav,Bevolxyz,Bxyz,gradhs,alphaMM,
     &              iunique,ekcle,listpm,spinx,spiny,spinz,angaddx,
     &              angaddy,angaddz,spinadx,spinady,spinadz,
     &              listinactive)
            ELSE
               print *,iproc,': MAINTAIN PARTICLE SPLIT ',npart,nhydro
               npart = nhydro
c               n1 = nhydro
            ENDIF
c         ENDIF
c      ENDIF

10000 CONTINUE
c
c--End if to allow skipping of blocks of file belonging to other processes
c
      ENDIF
#else
c
c--If not MPI calculation, then increment icount so that next block can be read
c
      IF (nblocks.NE.1) THEN
         DO i = icountsink+1, icountsink+nptmass
            listpm(i) = listpm(i) + icount
         END DO
      ENDIF
      icount = icount + nhydro
      icountinactive = icountinactive + nlistinactive
      icountsink = icountsink + nptmass
      print *,'Totals at this point npart, nptmass ',icount,icountsink
      IF (nblocks.NE.1) THEN
         IF (iblock.EQ.nblocks) THEN
            IF (icount.NE.npart) THEN
               WRITE (*,*) 'ERROR - rdump: icount.NE.npart ',
     &              icount,npart,nblocks,iblock
               CALL quit(0)
            ENDIF
            nptmass = icountsink
            nlistinactive = icountinactive
         ENDIF
      ENDIF
#endif
c
c--If found block for this process, can skip rest of file
c
#ifdef MPI
      IF (nblocks.NE.1 .AND. iproc.EQ.iblock-1) THEN
#ifdef MPIDEBUG
         print *,iproc,': Skipping remainder of dump file ',iblock
#endif
         GOTO 11000
      ENDIF
#endif
c
c--End loop over blocks
c
      END DO
#ifdef MPI
11000 CONTINUE
#endif
#ifdef MPI
c
c--This is the second section that is required to force sequential MPI
c     reads of the dump file by different MPI processors.
c
      IF (isequential) THEN
         IF (iproc.LT.numproc-1) THEN
#ifdef MPIDEBUG
            print *,iproc,' Finished reading'
#endif
            nowgo = 1
            CALL MPI_SEND(nowgo,1,MPI_INTEGER,iproc+1,99,MPI_COMM_WORLD,
     &           ierr)
         ELSE
#ifdef MPIDEBUG
            print *,iproc,' Finished reading'
#endif
         ENDIF
      ENDIF
c
c--Set global number of sink particles
c
      CALL MPI_ALLREDUCE(nptmass,nptmasstot,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
      print *,iproc,' FINISHED READ ',npart,nptmass,nptmasstot
#else
      nptmasstot = nptmass
      print *,' FINISHED READ ',npart,nptmass
#endif

      IF (iskip.EQ.1) RETURN
c
c--End reading of dump file
c--------------------------
c
      IF (itrace.EQ.'all') WRITE (*, 99002)
99002 FORMAT (' exit subroutine rdump')
      RETURN

c
c--File truncated
c
 100  CONTINUE
      WRITE(*,*) 'END OF FILE REACHED IN RDUMP'
      ichkl = 1

      RETURN
      
      CONTAINS

      !------------------------------------------
      ! extraction of single real*8 variables
      !------------------------------------------
      SUBROUTINE extract(tag,rval,rarr,tags,ntags,ierr)
      IMPLICIT NONE
      CHARACTER(LEN=*),      INTENT(IN)  :: TAG
      REAL,                  INTENT(OUT) :: RVAL
      REAL,                  INTENT(IN)  :: RARR(:)
      CHARACTER(LEN=16),     INTENT(IN)  :: TAGS(:)
      INTEGER,               INTENT(IN)  :: NTAGS
      INTEGER,               INTENT(OUT) :: IERR
      LOGICAL :: MATCHED
      INTEGER :: I

      ierr = 1
      matched = .FALSE.
      rval = 0.0 ! default if not found
      over_tags: do i=1,min(ntags,size(tags))
         if (trim(tags(i))==trim(adjustl(tag))) then
            if (size(rarr) >= i) then
               rval = rarr(i)
               matched = .true.
            endif
            exit over_tags  ! only match first occurrence
         endif
      enddo over_tags
      if (matched) ierr = 0
      if (ierr.NE.0) print "(a)",' ERROR: could not find '//
     &                            trim(adjustl(tag))//' in header'

      END SUBROUTINE extract

c     --------------------------------------------
c       Verify that tag matches expectation
c     --------------------------------------------
      SUBROUTINE checktag(tag,tagi,ierr)
      IMPLICIT NONE
      CHARACTER(LEN=*),      INTENT(IN)  :: tag
      CHARACTER(LEN=*),      INTENT(IN)  :: tagi
      INTEGER,               INTENT(OUT) :: ierr

      ierr = 0
      IF (TRIM(tag).NE.TRIM(tagi)) THEN
        ierr = 1
        WRITE(*,99980)tagi,tag
      ENDIF
99980 FORMAT(' WARNING: read tag ',A16,' but expected ',A16)

      END SUBROUTINE checktag

      END SUBROUTINE rdump

c============================================================

#ifdef MPI
      SUBROUTINE split(npart,n1,nptmass,icount,nhydro,nlistinactive,
     &     isteps,iphase,
     &     xyzmh,vxyzu,rho,dgrav,Bevolxyz,Bxyz,gradhs,alphaMM,iunique,
     &     ekcle,listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &     spinadx,spinady,spinadz,listinactive)
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

      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'

      INTEGER*1 iphase(idim2)
      INTEGER*8 iunique(idim)
      REAL*4 rho(idim), dgrav(idim), gradhs(2,idim),
     &     alphaMM(isizealphaMM,idim)
      DIMENSION isteps(idim), xyzmh(5,idim), vxyzu(4,idim), 
     &     Bevolxyz(imhdevol,imhd), Bxyz(3,imhd), ekcle(5,iradtrans)
      DIMENSION listpm(iptdim),spinx(iptdim),spiny(iptdim),
     &     spinz(iptdim),angaddx(iptdim),angaddy(iptdim),
     &     angaddz(iptdim),spinadx(iptdim),spinady(iptdim),
     &     spinadz(iptdim)
      DIMENSION listinactive(idim)

      INTEGER classify

      DIMENSION numdiv(10)

      DIMENSION imap(idim)

      DO i = 1, idim
         imap(i) = 0
      END DO

      icountnew = icount
c
c--Type of domain decomposition
c
c--Disc-type domain decomposition
c
      IF (mpitype.EQ.'d' .OR. mpitype.EQ.'r' .OR. mpitype.EQ.'p') THEN
         DO i = icount+1, nhydro
            ibelong = classify(mpitype,xyzmh(1,i))
            IF (iproc.EQ.ibelong) THEN
               icountnew = icountnew + 1
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
                  imap(i) = icountnew
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
                  IF (iradtrans.EQ.idim) THEN
                     DO j = 1, 5
                        ekcle(j,icountnew) = ekcle(j,i)
                     END DO
                  ENDIF
                  IF (imhd.EQ.idim) THEN
                     DO j = 1, imhdevol
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
                  imap(i) = icountnew
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
                  IF (iradtrans.EQ.idim) THEN
                     DO j = 1, 5
                        ekcle(j,icountnew) = ekcle(j,i)
                     END DO
                  ENDIF
                   IF (imhd.EQ.idim) THEN
                     DO j = 1, imhdevol
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
                  imap(i) = icountnew
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
                  IF (iradtrans.EQ.idim) THEN
                     DO j = 1, 5
                        ekcle(j,icountnew) = ekcle(j,i)
                     END DO
                  ENDIF
                   IF (imhd.EQ.idim) THEN
                     DO j = 1, imhdevol
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
                  imap(i) = icountnew
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
                  IF (iradtrans.EQ.idim) THEN
                     DO j = 1, 5
                        ekcle(j,icountnew) = ekcle(j,i)
                     END DO
                  ENDIF
                   IF (imhd.EQ.idim) THEN
                     DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j = 1, imhdevol
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
               imap(i) = icountnew
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
               IF (iradtrans.EQ.idim) THEN
                  DO j = 1, 5
                     ekcle(j,icountnew) = ekcle(j,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN               
                  DO j = 1, imhdevol
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
         CALL quit(0)
      ENDIF

      ELSE
         PRINT *,' Currently invalid MPI domain decompostion ',numproc,
     &        mpitype
         CALL quit(0)
      ENDIF
      
      nptmassnew = 0
      DO i = 1, nptmass
         IF (imap(listpm(i)).NE.0) THEN
            nptmassnew = nptmassnew + 1
            listpm(nptmassnew) = imap(listpm(i))
            IF (iphase(listpm(nptmassnew)).LE.0 .OR. 
     &           iphase(listpm(nptmassnew)).GE.10) THEN
               WRITE (*,*) 'ERROR - iphase(listpm(nptmassnew)).LE.0 ',
     &              nptmassnew,iphase(listpm(nptmassnew))
               CALL quit(0)
            ENDIF
            spinx(nptmassnew) = spinx(i)
            spiny(nptmassnew) = spiny(i)
            spinz(nptmassnew) = spinz(i)
            angaddx(nptmassnew) = angaddx(i)
            angaddy(nptmassnew) = angaddy(i)
            angaddz(nptmassnew) = angaddz(i)
            spinadx(nptmassnew) = spinadx(i)
            spinady(nptmassnew) = spinady(i)
            spinadz(nptmassnew) = spinadz(i)
         ENDIF
      END DO
      nptmass = nptmassnew

      nlistinactivenew = 0
      DO i = 1, nlistinactive
         IF (imap(listinactive(i)).NE.0) THEN
            nlistinactivenew = nlistinactivenew + 1
            listinactive(nlistinactivenew) = imap(listinactive(i))
            IF (iphase(listinactive(nlistinactivenew)).NE.-1) THEN
               WRITE (*,*) 'ERROR - iphase(listinactive(',
     &              'nlistinactivenew)).NE.-1 ',
     &              nlistinactivenew,listinactive(nlistinactivenew),
     &              iphase(listinactive(nlistinactivenew))
               CALL quit(0)
            ENDIF
         ENDIF
      END DO
      nlistinactive = nlistinactivenew

      npart = icountnew
c      n1 = icountnew
      icount = icountnew

      RETURN
      END
#endif

