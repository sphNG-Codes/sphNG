      SUBROUTINE wdump(idisk1)
c************************************************************
c                                                           *
c  This routine writes a dump on disk                       *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
c      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/savernd'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/f2'
      INCLUDE 'COMMONS/h2'
      INCLUDE 'COMMONS/perform'

      DIMENSION itempsort(idim)
      EQUIVALENCE (itempsort, next1)

      CHARACTER*7 where
      CHARACTER*100 fileident
      INTEGER*4 int1, int2
      INTEGER*8 number8
      DIMENSION nums(8)
c
c--Standard numbers
c
      PARAMETER (int1=690706,int2=780806)

      DATA where/'wdump'/
c
c=-set nprint = npart by default
c
      nprint = npart
c
c--uncomment these lines to dump ghosts as well
c  (useful for debugging)
c
c#ifdef DUMPGHOSTS
c      nprint = npart + nghost
c      DO i=npart+1,npart+nghost
c         isort(i) = i
c      ENDDO
c#endif
c
c--Write
c
      irec = irec + 1
      iresort = iresort + 1
      ifulldump = ifulldump + 1
      IF (ifulldump.EQ.nfullstep) THEN
         ifulldump = 0
c
c--Write full dump file
c----------------------
      i1 = int1
      r1 = real(int2)
c
c--Write output file
c
c#ifdef MPI
c      IF (iproc.EQ.0) THEN
c#endif
      WRITE (idisk1, ERR=100) int1,r1,int2,i1,int1
      WRITE (idisk1, ERR=100) fileident('F')
c#ifdef MPI
c      ENDIF
c#endif
c
c--Single values
c
#ifdef MPI
#ifdef MPIDEBUG
      print *,iproc,': Reduction '
#endif
      nblocks = numproc
      CALL MPI_REDUCE(npart,nparttot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(n1,n1tot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(n2,n2tot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nreassign,nreassigntot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(naccrete,naccretetot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nkill,nkilltot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(escap,escaptot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tkin,tkintot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tgrav,tgravtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tterm,ttermtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(anglostx,anglostxtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(anglosty,anglostytot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(anglostz,anglostztot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(specang,specangtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(ptmassin,ptmassintot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmag,tmagtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
#ifdef MPIDEBUG
      print *,iproc,': Reduction finished '
#endif
#else
      nblocks = 1

      nparttot = npart
      n1tot = n1
      n2tot = n2
      nreassigntot = nreassign
      naccretetot = naccrete
      nkilltot = nkill
      escaptot = escap
      tkintot = tkin
      tgravtot = tgrav
      ttermtot = tterm
      anglostxtot = anglostx
      anglostytot = anglosty
      anglostztot = anglostz
      specangtot = specang
      ptmassintot = ptmassin
      tmagtot = tmag
#endif

c#ifdef MPI
c      IF (iproc.EQ.0) THEN
c#endif
c
c--Default int
      number = 10 + NTAB
      WRITE (idisk1, ERR=100) number
      WRITE (idisk1, ERR=100) nparttot,n1tot,n2tot,nreassigntot,
     &     naccretetot,nkilltot,nblocks,iyr,idum,
     &     nlistinactive,(iv(i),i=1,NTAB)

c--int*1, int*2, int*4
      number = 0
      DO i = 1, 3
         WRITE (idisk1, ERR=100) number
      END DO
c--int*8
      number = 1
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100) iuniquemax
c--Default real
      number = 21
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100) gt, dtmax, gamma, rhozero, RK2,
     &     escaptot, tkintot, tgravtot, ttermtot, 
     &     anglostxtot, anglostytot, anglostztot,
     &     specangtot, ptmassintot, tmagtot, Bextx, Bexty, Bextz,
     &     hzero, uzero_n2, hmass
c--real*4
      number = 0
      WRITE (idisk1, ERR=100) number
c--real*8
      number = 4
      WRITE (idisk1, ERR=100) number
      WRITE (idisk1, ERR=100) udist, umass, utime, umagfd
c
c--Arrays
c
c
c--Number of array lengths
c     This has up to 4 blocks (hydro variables, sink particles, radiative
c     transfer, and MHD) repeated nblocks times for each of nblocks MPI
c     processes
c
      number = 2
      IF (encal.EQ.'r') number = 3
      IF (imhd.EQ.idim) number = 4
      nblockarrays = number*nblocks
      WRITE (idisk1, ERR=100) nblockarrays
#ifdef MPI
c      ENDIF
      nowgo = 0
#ifdef MPIDEBUG
      print *,iproc,': About to write'
#endif
      IF (iproc.GT.0) THEN
         CALL MPI_RECV(nowgo,1,MPI_INTEGER,iproc-1,99,MPI_COMM_WORLD,
     &        istatus,ierr)
#ifdef MPIDEBUG
         print *,iproc,': Opening file ',nprint
#endif
         CALL file
#ifdef MPIDEBUG
         print *,iproc,': Opened file'
#endif
      ENDIF
#endif
c
c--Array length 1 header
c
      number8 = nprint
      nums(1) = 2
      nums(2) = 1
      nums(3) = 0
      nums(4) = 0
      nums(5) = 1
      nums(6) = 9
      IF (nlmax.EQ.1) THEN
         nums(7) = 4
      ELSE
         nums(7) = 3
      ENDIF
      nums(8) = 0
      WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 2 header
c
      number8 = nptmass
      nums(1) = 1
      nums(2) = 0
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 9
      nums(7) = 0
      nums(8) = 0
      WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 3 header
c
      IF (number.GE.3) THEN
         IF (encal.EQ.'r') THEN
            number8 = nprint
         ELSE
            number8 = 0
         ENDIF
         nums(1) = 0
         nums(2) = 0
         nums(3) = 0
         nums(4) = 0
         nums(5) = 0
         IF (encal.EQ.'r') THEN
            nums(6) = 5
         ELSE
            nums(6) = 0
         ENDIF
         nums(7) = 0
         nums(8) = 0
         WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
      ENDIF
c
c--Array length 4 header
c
      IF (number.GE.4) THEN
         IF (imhd.EQ.idim) THEN
            number8 = nprint
         ELSE
            number8 = 0
         ENDIF
         nums(1) = 0
         nums(2) = 0
         nums(3) = 0
         nums(4) = 0
         nums(5) = 0
         IF (imhd.EQ.idim) THEN
            IF (varmhd.EQ.'eulr') THEN
               nums(6) = 5
            ELSEIF (varmhd.EQ.'vecp') THEN
               nums(6) = 6
            ELSE
               nums(6) = 3
            ENDIF
         ENDIF
         nums(7) = 5
         nums(8) = 0
         WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
      ENDIF
#ifdef MPIDEBUG
         print *,iproc,': Written array headers'
#endif
c
c--Array length 1 arrays
c      
c--Default int
      WRITE (idisk1, ERR=100) (isteps(isort(i)), i=1, nprint)
      WRITE (idisk1, ERR=100) (iorig(listinactive(i)),  
     &     i=1, nprint)
c--int*1
      WRITE (idisk1, ERR=100) (iphase(isort(i)), i=1, nprint)
c--int*2

c--int*4

c--int*8
c
c--NOTE: iunique is not sorted
c
      WRITE (idisk1, ERR=100) (iunique(i), i=1, nprint)
c--Default real
      DO j = 1, 5
         WRITE (idisk1, ERR=100) (xyzmh(j,isort(i)), i=1, nprint)
      END DO
      DO j = 1, 4
         WRITE (idisk1, ERR=100) (vxyzu(j,isort(i)), i=1, nprint)
      END DO      
c--real*4
      WRITE (idisk1, ERR=100) (rho(isort(i)), i=1, nprint)
      IF (nlmax.EQ.1) THEN
         WRITE (idisk1, ERR=100) (gradhs(1,isort(i)), i=1, nprint)      
         WRITE (idisk1, ERR=100) (gradhs(2,isort(i)), i=1, nprint)
      ELSE
         WRITE (idisk1, ERR=100) (dgrav(isort(i)), i=1, nprint)      
      ENDIF
      WRITE (idisk1, ERR=100) (alphaMM(1,isort(i)), i=1, nprint)
c--real*8

#ifdef MPIDEBUG
         print *,iproc,': Written hydro arrays'
#endif
c
c--Array length 2 arrays
c

c--Default int
      WRITE (idisk1, ERR=100) (iorig(listpm(i)), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      WRITE (idisk1, ERR=100) (spinx(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (spiny(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (spinz(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (angaddx(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (angaddy(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (angaddz(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (spinadx(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (spinady(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

#ifdef MPIDEBUG
         print *,iproc,': Written sink particle arrays'
#endif
      IF (encal.EQ.'r' .AND. iradtrans.EQ.idim) THEN
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
            WRITE (idisk1, ERR=100) (ekcle(j,isort(i)), i=1, nprint)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (imhd.EQ.idim) THEN
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
c--Dump B (not whatever the evolved MHD variable is)
c
         IF (varmhd.EQ.'Brho') THEN
            DO j = 1, 3
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i))*rho(isort(i)), i=1, nprint)
            END DO
         ELSEIF (varmhd.EQ.'Bvol') THEN
            DO j = 1, 3
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO         
         ELSE
            DO j = 1, 3
               WRITE (idisk1, ERR=100) 
     &            (Bxyz(j,isort(i)), i=1, nprint)
            END DO         
         ENDIF
c
c--Dump Euler/vector potentials if necessary
c
         IF (varmhd.EQ.'eulr') THEN
            DO j = 1, 2
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO
         ELSEIF (varmhd.EQ.'vecp') THEN
            DO j = 1, 3
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO
         ENDIF
c--real*4
         DO j = 1, 4
            WRITE (idisk1, ERR=100) (divcurlB(j,isort(i)), i=1, nprint)
         ENDDO
         WRITE(idisk1, ERR=100) (alphaMM(isizealphaMM,isort(i)),
     &        i=1,nprint)
c--real*8

      ENDIF
      CALL FLUSH (idisk1)
      CLOSE (idisk1)
#ifdef MPI
      IF (iproc.LT.numproc-1) THEN
#ifdef MPIDEBUG
         print *,iproc,' Finished writing'
#endif
         nowgo = 1
         CALL MPI_SEND(nowgo,1,MPI_INTEGER,iproc+1,99,MPI_COMM_WORLD,
     &        ierr)
      ELSE
#ifdef MPIDEBUG
         print *,iproc,' Finished writing'
#endif
      ENDIF
#endif
c
c--End writing of full dump file
c-------------------------------
c
      ELSE
c
c--Write small dump file
c----------------------------

      i1 = int1
      r1 = real(int2)
c
c--Write output file
c
      WRITE (idisk1, ERR=100) int1,r1,int2,i1,int1
      WRITE (idisk1, ERR=100) fileident('S')
c
c--Single values
c
c--Default int
      number = 6
      WRITE (idisk1, ERR=100) number
      WRITE (idisk1, ERR=100) npart,n1,n2,nreassign,naccrete,nkill
c--int*1, int*2, int*4, int*8
      number = 0
      DO i = 1, 4
         WRITE (idisk1, ERR=100) number
      END DO
c--Default real
      number = 15
c
c     dump mass array ONLY if the particle masses are not equal
c     for equal mass particles, pmassinitial = particle mass
c     for non-equal mass particles, pmassinitial = 0 and the mass array is dumped
c
      DO i = 1,npart
         IF (iphase(i).EQ.0) THEN
            pmassinitial = xyzmh(4,i)
            GOTO 35
         ENDIF
      ENDDO
 35   jlen = 3
      DO i = 1, npart
         IF (iphase(i).EQ.0) THEN
            IF (xyzmh(4,i).NE.pmassinitial) THEN
               pmassinitial = 0.
               jlen = 4
               GOTO 40
            ENDIF
         ENDIF
      END DO
 40   WRITE (idisk1, ERR=100) number
      WRITE (idisk1, ERR=100) gt, dtmax, gamma, rhozero, RK2,
     &     escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &     specang, ptmassin, pmassinitial
c--real*4
      number = 0
      WRITE (idisk1, ERR=100) number
c--real*8
      number = 4
      WRITE (idisk1, ERR=100) number
      WRITE (idisk1, ERR=100) udist, umass, utime, umagfd
c
c--Arrays
c
c--Number of array lengths
c
      number = 2
      IF (encal.EQ.'r') number = 3
      IF (imhd.EQ.idim) number = 4
      WRITE (idisk1, ERR=100) number
c
c--Array length 1 header
c
      number8 = npart
      nums(1) = 0
      nums(2) = 1
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = jlen
      nums(7) = 2
      nums(8) = 0
      WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 2 header
c
      number8 = nptmass
      nums(1) = 1
      nums(2) = 0
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 1
      nums(7) = 0
      nums(8) = 0
      WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 3 header
c
      IF (number.GE.3) THEN
         IF (encal.EQ.'r') THEN
            number8 = npart
         ELSE
            number8 = 0
         ENDIF
         nums(1) = 0
         nums(2) = 0
         nums(3) = 0
         nums(4) = 0
         nums(5) = 0
         IF (encal.EQ.'r') THEN
            nums(6) = 1
         ELSE
            nums(6) = 0
         ENDIF
         nums(7) = 0
         nums(8) = 0
         WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
      ENDIF
c
c--Array length 4 header
c
      IF (number.GE.4) THEN
         IF (imhd.EQ.idim) THEN
            number8 = npart
         ELSE
            number8 = 0
         ENDIF
         nums(1) = 0
         nums(2) = 0
         nums(3) = 0
         nums(4) = 0
         nums(5) = 0
         IF (imhd.EQ.idim) THEN
            nums(6) = 3
         ELSE
            nums(6) = 0
         ENDIF
         nums(7) = 0
         nums(8) = 0
         WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
      ENDIF
c
c--Array length 1 arrays
c
c--Default int
c      WRITE (idisk1, ERR=100) (isteps(isort(i)), i=1, npart)
c--int*1
      WRITE (idisk1, ERR=100) (iphase(isort(i)), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      DO j = 1, jlen
         WRITE (idisk1, ERR=100) (xyzmh(j,isort(i)), i=1, npart)
      END DO
c      DO j = 1, 4
c         WRITE (idisk1, ERR=100) (vxyzu(j,isort(i)), i=1, npart)
c      END DO
c--real*4
      WRITE (idisk1, ERR=100) (rho(isort(i)), i=1, npart)
c     dump smoothing length as a real*4 to save space
      DO i = 1, npart
         dq(i) = xyzmh(5,i)
      END DO
      WRITE (idisk1, ERR=100) (dq(isort(i)), i=1, npart)      
c      WRITE (idisk1, ERR=100) (dgrav(isort(i)), i=1, npart)      
c     WRITE (idisk1, ERR=100) (alphaMM(1,isort(i)), i=1, npart)
c--real*8

c
c--Array length 2 arrays
c
c--Default int
      WRITE (idisk1, ERR=100) (iorig(listpm(i)), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      WRITE (idisk1, ERR=100) (xyzmh(4,listpm(i)), i=1,nptmass)
c      WRITE (idisk1, ERR=100) (spinx(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (spiny(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (spinz(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (angaddx(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (angaddy(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (angaddz(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (spinadx(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (spinady(i),i=1,nptmass)
c      WRITE (idisk1, ERR=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

      IF (encal.EQ.'r' .AND. iradtrans.EQ.idim) THEN
c
c--Array length 3 arrays
c
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 1
            WRITE (idisk1, ERR=100) (ekcle(j,isort(i)), i=1, npart)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (imhd.EQ.idim) THEN
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
c  Dump B not B/rho
c
         DO j = 1, 3
            WRITE (idisk1, ERR=100) (Bxyz(j,isort(i)), i=1, npart)
         END DO
c--real*4

c--real*8

      ENDIF
      CALL FLUSH (idisk1)
      CLOSE (idisk1)
c
c--End writing of small dump file
c--------------------------------
c
      ENDIF
c
c--Sort particles to ensure most efficient running.  Note that this 
c     should not be visible to the outside observer.  In other words,
c     an array must be kept of original list of particles and this
c     must be used to index *ANY* value from an array which is written 
c     to the outside.  This requires modification to almost every output
c     line in the code.  Done 21 Nov 2000.
c
      IF (itiming) CALL getused(tsort1)

      IF (gt.EQ.0.) GOTO 777
c      GOTO 777

      iresort = 0

      IF (nptmass.GT.0) THEN
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
               tempsort(i) = LOG(REAL(isteps(i)))/LOG(2.0) + 
     &              (xyzmh(1,i)-xminimum)/xrange
               istepmin = MIN(istepmin, isteps(i))
               istepmax = MAX(istepmax, isteps(i))
            ENDIF
         END DO
c
c--Sort particles based on their individual timesteps and x
c
         CALL indexx(npart, llist, tempsort, itempsort)
c
c--Move into llist.  This now stores 'iorignew'
c
         DO i = 1, npart
            llist(i) = itempsort(i)
         END DO
      ELSE
c
c--Sort particles based on a walk through the tree
c
         CALL treelist(npart, llist, dumxyzmh)
      ENDIF
c
c--Now re-order arrays
c
      DO i = 1, npart
         tempsort(i) = xyzmh(1,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(2,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(3,llist(i))
       END DO
      DO i = 1, npart
         xyzmh(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(4,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(5,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(5,i) = tempsort(i)
      END DO
c
c--vxyzu
c
      DO i = 1, npart
        tempsort(i) = vxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(4,i) = tempsort(i)
      END DO
c
c--dumxyzmh
c
      DO i = 1, npart
         tempsort(i) = dumxyzmh(1,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(2,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(3,llist(i))
       END DO
      DO i = 1, npart
         dumxyzmh(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(4,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(5,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(5,i) = tempsort(i)
      END DO
c
c--dumvxyzu
c
      DO i = 1, npart
        tempsort(i) = dumvxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumvxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumvxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumvxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(4,i) = tempsort(i)
      END DO
c
c--rho
c
      DO i = 1, npart
         tempsort(i) = rho(llist(i))
      END DO
      DO i = 1, npart
         rho(i) = tempsort(i)
      END DO
c
c--divv
c
      DO i = 1, npart
         tempsort(i) = divv(llist(i))
      END DO
      DO i = 1, npart
         divv(i) = tempsort(i)
      END DO
c
c--curlv
c
      DO i = 1, npart
         tempsort(i) = curlv(llist(i))
      END DO
      DO i = 1, npart
         curlv(i) = tempsort(i)
      END DO
c
c--pr
c
      DO i = 1, npart
         tempsort(i) = pr(llist(i))
      END DO
      DO i = 1, npart
         pr(i) = tempsort(i)
      END DO
c
c--vsound
c
      DO i = 1, npart
         tempsort(i) = vsound(llist(i))
      END DO
      DO i = 1, npart
         vsound(i) = tempsort(i)
      END DO
c
c--poten
c
      DO i = 1, npart
         tempsort(i) = poten(llist(i))
      END DO
      DO i = 1, npart
         poten(i) = tempsort(i)
      END DO
c
c--dgrav
c
      DO i = 1, npart
         tempsort(i) = dgrav(llist(i))
      END DO
      DO i = 1, npart
         dgrav(i) = tempsort(i)
      END DO
c
c--gradhs
c
      DO k = 1, 2
         DO i = 1, npart
            tempsort(i) = gradhs(k,llist(i))
         END DO
         DO i = 1, npart
            gradhs(k,i) = tempsort(i)
         END DO
      END DO
c
c--isteps
c
      DO i = 1, npart
         itempsort(i) = isteps(llist(i))
      END DO
      DO i = 1, npart
         isteps(i) = itempsort(i)
      END DO
c
c--iphase
c
      DO i = 1, npart
         itempsort(i) = iphase(llist(i))
      END DO
      DO i = 1, npart
         iphase(i) = itempsort(i)
      END DO
c
c--it0
c
      DO i = 1, npart
         itempsort(i) = it0(llist(i))
      END DO
      DO i = 1, npart
         it0(i) = itempsort(i)
      END DO
c
c--it1
c
      DO i = 1, npart
         itempsort(i) = it1(llist(i))
      END DO
      DO i = 1, npart
         it1(i) = itempsort(i)
      END DO
c
c--it2
c
      DO i = 1, npart
         itempsort(i) = it2(llist(i))
      END DO
      DO i = 1, npart
         it2(i) = itempsort(i)
      END DO
c
c--f1vxyzu
c
      DO i = 1, npart
         tempsort(i) = f1vxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(4,i) = tempsort(i)
      END DO
c
c--f2vxyzu
c
      DO i = 1, npart
         tempsort(i) = f2vxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f2vxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f2vxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f2vxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(4,i) = tempsort(i)
      END DO
c
c--f1ha
c
      DO k = 1, 1+isizealphaMM
         DO i = 1, npart
            tempsort(i) = f1ha(k,llist(i))
         END DO
         DO i = 1, npart
            f1ha(k,i) = tempsort(i)
         END DO
c
c--f2ha
c
         DO i = 1, npart
            tempsort(i) = f2ha(k,llist(i))
         END DO
         DO i = 1, npart
            f2ha(k,i) = tempsort(i)
         END DO
      END DO
c
c--dumrho
c
      DO i = 1, npart
         tempsort(i) = dumrho(llist(i))
      END DO
      DO i = 1, npart
         dumrho(i) = tempsort(i)
      END DO
c
c--alphaMM
c
      DO k = 1, isizealphaMM
         DO i = 1, npart
            tempsort(i) = alphaMM(k,llist(i))
         END DO
         DO i = 1, npart
            alphaMM(k,i) = tempsort(i)
         END DO
c
c--dumalpha
c
         DO i = 1, npart
            tempsort(i) = dumalpha(k,llist(i))
         END DO
         DO i = 1, npart
            dumalpha(k,i) = tempsort(i)
         END DO
      END DO
c
c--ddv
c
      DO i = 1, npart
         tempsort(i) = ddv(llist(i))
      END DO
      DO i = 1, npart
         ddv(i) = tempsort(i)
      END DO
c
c--h2mol
c
      DO i = 1, npart
         tempsort(i) = h2mol(llist(i))
      END DO
      DO i = 1, npart
         h2mol(i) = tempsort(i)
      END DO
c
c--h2ratio
c
      DO i = 1, npart
         tempsort(i) = h2ratio(llist(i))
      END DO
      DO i = 1, npart
         h2ratio(i) = tempsort(i)
      END DO
c
c--ekcle, dumekcle
c
      IF (encal.EQ.'r') THEN
         DO k = 1, 5
            DO i = 1, npart
               tempsort(i) = ekcle(k,llist(i))
            END DO
            DO i = 1, npart
               ekcle(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = dumekcle(k,llist(i))
            END DO
            DO i = 1, npart
               dumekcle(k,i) = tempsort(i)
            END DO
         END DO
      ENDIF
c
c--Bevolxyz, dumBevolxyz, Bxyz, f1Bxyz, f2Bxyz
c
      IF (imhd.EQ.idim) THEN
         DO k = 1, 3
            DO i = 1, npart
               tempsort(i) = Bevolxyz(k,llist(i))
            END DO
            DO i = 1, npart
               Bevolxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = dumBevolxyz(k,llist(i))
            END DO
            DO i = 1, npart
               dumBevolxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = Bxyz(k,llist(i))
            END DO
            DO i = 1, npart
               Bxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = f1Bxyz(k,llist(i))
            END DO
            DO i = 1, npart
               f1Bxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = f2Bxyz(k,llist(i))
            END DO
            DO i = 1, npart
               f2Bxyz(k,i) = tempsort(i)
            END DO
         END DO
      ENDIF
c
c--iorig
c
      DO i = 1, npart
         itempsort(i) = iorig(llist(i))
      END DO
      DO i = 1, npart
         iorig(i) = itempsort(i)
      END DO
c
c--isortnew
c
      DO i = 1, npart
         itempsort(llist(i)) = i
      END DO
      DO i = 1, nptmass
         listpm(i) = itempsort(listpm(i))
         listrealpm(listpm(i)) = i
      END DO
      DO i = 1, nlstacc
         listacc(i) = itempsort(listacc(i))
      END DO
      DO i = 1, nlistinactive
         listinactive(i) = itempsort(listinactive(i))
      END DO
      DO i = 1, nbinmax
         DO j = 1, nlstbins(i)
            listbins(j,i) = itempsort(listbins(j,i))
         END DO
      END DO
      DO i = 1, npart
         isort(i) = itempsort(isort(i))
      END DO
c
c--Must rebuild the ghosts and make the tree again after sorting
c
      CALL ghostp(ntot, npart, xyzmh, vxyzu, ekcle, Bevolxyz)

      DO i = 1, ntot
         DO j = 1, 5
            dumxyzmh(j,i) = xyzmh(j,i)
         END DO
         DO j = 1, 4
            dumvxyzu(j,i) = vxyzu(j,i)
         END DO
      END DO

      IF (itiming) CALL getused(tsortmtree1)

      CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)

      IF (itiming) THEN
         CALL getused(tsortmtree2)
         tsortmtree = tsortmtree + (tsortmtree2 - tsortmtree1)
      ENDIF
 777  CONTINUE

      IF (itiming) THEN
         CALL getused(tsort2)
         tsort = tsort + (tsort2 - tsort1)
      ENDIF
c
c--Zero torques
c
c      DO i = 1, npart
c         torqt(i) = 0.
c         torqg(i) = 0.
c         torqp(i) = 0.
c         torqv(i) = 0.
c         torqc(i) = 0.
c      END DO

      ENDFILE idisk1
      BACKSPACE idisk1

      ipos = irec

      RETURN
c
c--An error as occured while writing
c
 100  CALL error(where, 1)

      RETURN
      END

      FUNCTION fileident(firstchar)
c*********************************************************************
c  contruct header string based on compile-time options
c  these are for information only (ie. not important for restarting)
c
c*********************************************************************
      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/varmhd'
      CHARACTER*100 fileident
      CHARACTER*1 firstchar

      fileident = firstchar//'sphNG:'
      IF (nlmax.EQ.1) THEN
         fileident = fileident(1:7) // 'gradh=on, '
      ELSE
         fileident = fileident(1:7) // 'gradh=off,'
      ENDIF
      IF (imhd.EQ.idim) THEN
         fileident = fileident(1:17) // 'MHD='//varmhd(1:4)//','
      ELSE
         fileident = fileident(1:17) // 'MHD=none,'
      ENDIF
      IF (iradtrans.EQ.idim) THEN
         fileident = fileident(1:26) // 'RT=on '
      ELSE
         fileident = fileident(1:26) // 'RT=off'
      ENDIF

      END FUNCTION fileident
