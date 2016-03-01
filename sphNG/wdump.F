      SUBROUTINE wdump(idisk1)
c************************************************************
c                                                           *
c  This routine ONLY writes a dump on disk                  *
c                                                           *
c************************************************************
      !IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
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
c#ifdef DUMPGHOSTS
c      INCLUDE 'COMMONS/ghost'
c#endif
      INCLUDE 'COMMONS/numpa'
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
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'

      COMMON /plotinfo/ xtype(2,idim), iplotdump

      DIMENSION itempsort(idim)
c      EQUIVALENCE (itempsort, next1)

      CHARACTER*7 where
      CHARACTER*100 fileident
      CHARACTER*16 tagarr(10)
      INTEGER*4 int1, int2, int1o
      INTEGER*8 number8
      DIMENSION nums(8)
      INTEGER iversion
      CHARACTER*1 labelx(3)
      PARAMETER (labelx = (/'x','y','z'/))
c
c--Standard numbers
c
      PARAMETER (int1=060769,int2=060878)
      PARAMETER (int1o=690706)

      DATA where/'wdump'/
    
 
      tagarr = ''
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
      ifulldump = ifulldump + 1
      IF (ifulldump.EQ.nfullstep) THEN
         ifulldump = 0
      ENDIF
c
c--Write dump file
c----------------------
      r1 = real(int2)
c
c--Format version. Change this ONLY if the dump format
c  is incompatible with previous formats
c
      iversion = 1
c
c--Write output file
c
#ifdef MPICOPY
      IF (iproc.EQ.0) THEN
#endif
      WRITE (idisk1, ERR=100) int1,r1,int2,iversion,int1o
      IF (ifulldump.EQ.0) THEN
         WRITE (idisk1, ERR=100) fileident('FT')
      ELSE
         WRITE (idisk1, ERR=100) fileident('ST')
      ENDIF
#ifdef MPICOPY
      ENDIF
#endif
c
c--Single values
c
      nparttot = npart
#ifdef MPI
#ifdef MPIDEBUG
      print *,iproc,': Reduction '
#endif
      nblocks = numproc
      CALL MPI_REDUCE(npart,nparttot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
c      CALL MPI_REDUCE(n1,n1tot,1,MPI_INTEGER,MPI_SUM,0,
c     &     MPI_COMM_WORLD,ierr)
c      CALL MPI_REDUCE(n2,n2tot,1,MPI_INTEGER,MPI_SUM,0,
c     &     MPI_COMM_WORLD,ierr)
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

c      n1tot = n1
c      n2tot = n2
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

#ifdef MPICOPY
      IF (iproc.EQ.0) THEN
#endif
c
c--Default int
      number = 12 + NTAB
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100) tag('nparttot'),tag('n1'),tag('n2'),
     &      tag('nreassigntot'),tag('naccretetot'),tag('nkilltot'),
     &      tag('nblocks'),tag('iyr'),tag('idum'),(tag('iv'),i=1,NTAB),
     &      tag('iplanetesimals'),tag('irotpot'),tag('idragscheme')

      WRITE (idisk1, ERR=100) nparttot,n1,n2,nreassigntot,
     &     naccretetot,nkilltot,nblocks,iyr,idum,
     &     (iv(i),i=1,NTAB), iplanetesimals, irotpot,
     &     idragscheme
c--int*1, int*2, int*4
      number = 0
      DO i = 1, 3
         WRITE (idisk1, ERR=100) number
      END DO
c--int*8
      number = 1
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100) tag('iuniquemax')
      WRITE (idisk1, ERR=100) iuniquemax
c--Default real
c
c--Full dump
c
      number = 37

      pmassinitial = 0.
      IF (ifulldump.NE.0) THEN
c
c--Small dump
c
c--Dump mass array ONLY if the particle masses are not equal
c     For equal mass particles, pmassinitial = particle mass
c     For non-equal mass particles, pmassinitial = 0 and mass array is dumped
c
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               pmassinitial = xyzmh(4,i)
               GOTO 35
            ENDIF
         ENDDO
 35      jlen = 3
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               IF (xyzmh(4,i).NE.pmassinitial) THEN
                  pmassinitial = 0.
                  jlen = 4
                  GOTO 40
               ENDIF
            ENDIF
         END DO
 40      CONTINUE
#ifdef MPI
         CALL MPI_REDUCE(jlen,jlentot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (jlentot.EQ.numproc*3) THEN
            jlen = 3
         ELSE
            pmassinitial = 0.
            jlen = 4
         ENDIF
#endif
      ENDIF

      WRITE (idisk1, ERR=100) number

c      print *,'stuff ',gt, dtmax, gamma, rhozero, RK2
      
      WRITE (idisk1, ERR=100) 
     &  tag('gt'),tag('dtmax'),tag('gamma'),tag('rhozero'),tag('RK2'),
     &  tag((/'escap','tkin ','tgrav','tterm'/)),
     &  tag((/'anglostx','anglosty','anglostz'/)),
     &  tag((/'specang ','ptmassin','tmag    '/)),
     &  tag((/'Bextx','Bexty','Bextz'/)),
     &  tag((/'hzero   ','uzero_n2','hmass   ','gapfac  '/)), 
     &  tag('pmassinitial'),tag('sdprof'),
     &  tag('rorbit_orig'),tag((/'min_rplan','max_rplan'/)),
     &  tag('planetesimalmass'),
     &  tag('coremass_orig'),tag('coremass'),
     &  tag('hfact'),tag('xmin'),tag('xmax'),tag('ymin'),
     &  tag('ymax'),tag('zmin'),tag('zmax')

      WRITE (idisk1, ERR=100) gt, dtmax, gamma, rhozero, RK2,
     &     escaptot, tkintot, tgravtot, ttermtot, 
     &     anglostxtot, anglostytot, anglostztot,
     &     specangtot, ptmassintot, tmagtot, Bextx, Bexty, Bextz,
     &     hzero, uzero_n2, hmass, gapfac, pmassinitial, sdprof,
     &     rorbit_orig, min_rplan, max_rplan, planetesimalmass,
     &     coremass_orig, coremass, hfact, xmin, xmax, ymin, ymax,
     &     zmin, zmax

c--real*4
      number = 0
      WRITE (idisk1, ERR=100) number
c--real*8
      number = 4
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100)
     &      tag((/'udist ','umass ','utime ','umagfd'/))
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
      IF (iradtrans.EQ.idim) number = 3
      IF (imhd.EQ.idim) number = 4
      nblockarrays = number*nblocks
      WRITE (idisk1, ERR=100) nblockarrays

#ifdef MPI
c      ENDIF
      nowgo = 0
#ifdef MPIDEBUG
      print *,iproc,': About to write ',nprint,ifulldump
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
c--Full dump file quantities
c
      IF (ifulldump.EQ.0) THEN
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
c
c--Extra information for plot dump
c
      IF (iplotdump.EQ.1) nums(6) = nums(6) + 2
c
c--For Clare's chemistry
c
      IF (iener.EQ.4) nums(6) = nums(6) + 5

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
         IF (iradtrans.EQ.idim) THEN
            number8 = nprint
            nums(1) = 0
            nums(2) = 0
            nums(3) = 0
            nums(4) = 0
            nums(5) = 0
            nums(6) = 5
            IF (idustRT.GT.0) nums(6) = nums(6) + 2
            nums(7) = 0
            IF (idustRT.GT.0) nums(7) = nums(7) + nchemistry + 
     &           nheatingISR + 2
            nums(8) = 0
         ELSE
            number8 = 0
            nums(:) = 0
         ENDIF
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
               nums(6) = imhdevol
            ENDIF
         ENDIF
         nums(7) = 7
c         nums(7) = 13
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
      WRITE (idisk1, ERR=100) tag('isteps')
      WRITE (idisk1, ERR=100) (isteps(isort(i)), i=1, nprint)
      DO i = nlistinactive + 1, nprint
         listinactive(i) = 1
      END DO
      WRITE (idisk1, ERR=100) tag('listinactive')
      WRITE (idisk1, ERR=100) nlistinactive,(iorig(listinactive(i)),  
     &     i=1, nprint-1)
c--int*1
      WRITE (idisk1, ERR=100) tag('iphase')
      WRITE (idisk1, ERR=100) (iphase(isort(i)), i=1, nprint)
c--int*2

c--int*4

c--int*8
c
c--NOTE: iunique is not sorted
c
      WRITE (idisk1, ERR=100) tag('iunique')
      WRITE (idisk1, ERR=100) (iunique(i), i=1, nprint)
c--Default real
      tagarr(1:5) = (/'x','y','z','m','h'/)
      DO j = 1, 5
         WRITE (idisk1, ERR=100) tag(tagarr(j))
         WRITE (idisk1, ERR=100) (xyzmh(j,isort(i)), i=1, nprint)
      END DO
      tagarr(1:4) = (/'vx','vy','vz','u '/)
      DO j = 1, 4
         WRITE (idisk1, ERR=100) tag(tagarr(j))
         WRITE (idisk1, ERR=100) (vxyzu(j,isort(i)), i=1, nprint)
      END DO     

      IF (iener.EQ.4) THEN
         WRITE (idisk1, ERR=100) tag('h2ratio')
         WRITE (idisk1, ERR=100) (h2ratio(isort(i)), i=1, nprint)
         WRITE (idisk1, ERR=100) tag('abHIq')
         WRITE (idisk1, ERR=100) (abHIq(isort(i)), i=1, nprint)
         WRITE (idisk1, ERR=100) tag('abhpq')
         WRITE (idisk1, ERR=100) (abhpq(isort(i)), i=1, nprint)
         WRITE (idisk1, ERR=100) tag('abeq')
         WRITE (idisk1, ERR=100) (abeq(isort(i)), i=1, nprint)
         WRITE (idisk1, ERR=100) tag('abco')
         WRITE (idisk1, ERR=100) (abco(isort(i)), i=1, nprint)         
      ENDIF       

      IF (iplotdump.EQ.1) THEN
         tagarr(1:2) = (/'gas ','dust'/)
         DO j = 1, 2
            WRITE (idisk1, ERR=100) tag(tagarr(j))
            WRITE (idisk1, ERR=100) (xtype(j,isort(i)), i=1, nprint)
         END DO     
      ENDIF
c--real*4
      WRITE (idisk1, ERR=100) tag('rho')
      WRITE (idisk1, ERR=100) (rho(isort(i)), i=1, nprint)
      IF (nlmax.EQ.1) THEN
         WRITE (idisk1, ERR=100) tag('gradh')
         WRITE (idisk1, ERR=100) (gradhs(1,isort(i)), i=1, nprint)      
         WRITE (idisk1, ERR=100) tag('gradsoft')
         WRITE (idisk1, ERR=100) (gradhs(2,isort(i)), i=1, nprint)
      ELSE
         WRITE (idisk1, ERR=100) tag('dgrav')
         WRITE (idisk1, ERR=100) (dgrav(isort(i)), i=1, nprint)      
      ENDIF
      WRITE (idisk1, ERR=100) tag('alpha')
      WRITE (idisk1, ERR=100) (alphaMM(1,isort(i)), i=1, nprint)
c--real*8

#ifdef MPIDEBUG
         print *,iproc,': Written hydro arrays'
#endif
c
c--Array length 2 arrays
c

c--Default int
      WRITE (idisk1, ERR=100) tag('listpm')
      WRITE (idisk1, ERR=100) (iorig(listpm(i)), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      WRITE (idisk1, ERR=100) tag('spinx')
      WRITE (idisk1, ERR=100) (spinx(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('spiny')
      WRITE (idisk1, ERR=100) (spiny(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('spinz')
      WRITE (idisk1, ERR=100) (spinz(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('angaddx')
      WRITE (idisk1, ERR=100) (angaddx(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('angaddy')
      WRITE (idisk1, ERR=100) (angaddy(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('angaddz')
      WRITE (idisk1, ERR=100) (angaddz(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('spinadx')
      WRITE (idisk1, ERR=100) (spinadx(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('spinady')
      WRITE (idisk1, ERR=100) (spinady(i),i=1,nptmass)
      WRITE (idisk1, ERR=100) tag('spinadz')
      WRITE (idisk1, ERR=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

#ifdef MPIDEBUG
         print *,iproc,': Written sink particle arrays'
#endif
      IF (iradtrans.EQ.idim) THEN
c
c--Array length 3 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         tagarr(1:5) = (/tag('Erad'),tag('opacity'),tag('Cv'),
     &                   tag('lambda'),tag('EddingtonFactor')/)
         DO j = 1, 5
            WRITE (idisk1, ERR=100) tag(tagarr(j))
            WRITE (idisk1, ERR=100) (ekcle(j,isort(i)), i=1, nprint)
         END DO
         IF (idustRT.GT.0) THEN
            tagarr(1:2) = (/tag('DustTemperature'),
     &           tag('DustOpacity')/)
            DO j = 1, 2
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) (dust_tk(j,isort(i)),i=1,nprint)
            END DO
         ENDIF
c--real*4
         IF (idustRT.GT.0) THEN
            tagarr(1:3) = (/tag('Cplus'),tag('AtomicC'),tag('CO')/)
            DO j = 1, nchemistry
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) (chemistry(j,isort(i)), i=1, 
     &              nprint)
            END DO

            tagarr(1:4) = (/tag('HeatingISR'),tag('Gphotoelectric'),
     &           tag('Extinction'),tag('H2SelfShielding')/)
            DO j = 1, nheatingISR
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) (heatingISR(j,isort(i)), i=1, 
     &              nprint)
            END DO

            WRITE (idisk1, ERR=100) tag('H2Frac')
            WRITE (idisk1, ERR=100) (h2frac(isort(i)), i=1, 
     &           nprint)
            WRITE (idisk1, ERR=100) tag('dH2dt')
            WRITE (idisk1, ERR=100) (dh2dt(1,isort(i)), i=1, 
     &           nprint)
         ENDIF
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
               WRITE (idisk1, ERR=100) tag('B'//labelx(j))
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i))*rho(isort(i)), i=1, nprint)
            END DO
            DO j = 4, imhdevol
               WRITE (idisk1, ERR=100) tag('psi/c_h')
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO
         ELSEIF (varmhd.EQ.'Bvol') THEN
            DO j = 1,imhdevol
               IF (j.LE.3) THEN
                  WRITE (idisk1, ERR=100) tag('B'//labelx(j))
               ELSE
                  WRITE (idisk1, ERR=100) tag('psi/c_h')
               ENDIF
               WRITE (idisk1, ERR=100)
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO         
         ELSE
            DO j = 1, 3
               WRITE (idisk1, ERR=100) tag('B'//labelx(j))
               WRITE (idisk1, ERR=100) 
     &            (Bxyz(j,isort(i)), i=1, nprint)
            END DO         
         ENDIF
c
c--Dump Euler/vector potentials if necessary
c
         tagarr(1:2)=(/'EulerAlpha','EulerBeta '/)
         IF (varmhd.EQ.'eulr') THEN
            DO j = 1, 2
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO
         ELSEIF (varmhd.EQ.'vecp') THEN
            DO j = 1, 3
               WRITE (idisk1, ERR=100) tag('A'//labelx(j))
               WRITE (idisk1, ERR=100) 
     &            (Bevolxyz(j,isort(i)), i=1, nprint)
            END DO
         ENDIF
c--real*4
         DO j = 1, 4
            IF (j.EQ.1) THEN
               WRITE (idisk1, ERR=100) tag('divB')
            ELSE
               WRITE (idisk1, ERR=100) tag('curlB'//labelx(j-1))
            ENDIF
            WRITE (idisk1, ERR=100) (divcurlB(j,isort(i)), i=1, nprint)
         ENDDO
         WRITE(idisk1, ERR=100) tag('alphaB')
         WRITE(idisk1, ERR=100) (alphaMM(isizealphaMM,isort(i)),
     &        i=1,nprint)
         WRITE(idisk1, ERR=100) tag('EtaReal')
         WRITE(idisk1, ERR=100) (etareal(isort(i)),
     &        i=1,nprint)
         WRITE(idisk1, ERR=100) tag('EtaArtificial')
         WRITE(idisk1, ERR=100) (etaartificial(isort(i)),
     &        i=1,nprint)
c--real*8

      ENDIF
c
c--End writing of full dump quantities
c-------------------------------------
c
      ELSE
c
c--Write small dump quantities
c-----------------------------
c
c--Array length 1 header
c
      number8 = npart
      nums(1) = 0
      nums(2) = 1
      nums(3) = 0
      nums(4) = 1
      nums(5) = 0
      nums(6) = jlen
      nums(7) = 3
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
         number8 = npart
         nums(1) = 0
         nums(2) = 0
         nums(3) = 0
         nums(4) = 0
         nums(5) = 0
         nums(6) = 0
         nums(7) = 2
         nums(8) = 0
         WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
      ENDIF
c
c--Array length 4 header
c
      IF (number.GE.4) THEN
         number8 = npart
         nums(1) = 0
         nums(2) = 0
         nums(3) = 0
         nums(4) = 0
         nums(5) = 0
         nums(6) = 0
         nums(7) = 3
         nums(8) = 0
         WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
      ENDIF
c
c--Array length 1 arrays
c
c--Default int
c      WRITE (idisk1, ERR=100) tag('isteps')
c      WRITE (idisk1, ERR=100) (isteps(isort(i)), i=1, npart)
c--int*1
      WRITE (idisk1, ERR=100) tag('iphase')
      WRITE (idisk1, ERR=100) (iphase(isort(i)), i=1, npart)
c--int*2

c--int*4
c
c--NOTE: iunique is not sorted
c     NOTE: For small dumps this is ASSUMED to fit within INT*4
c
      DO i = 1, npart
         IF (iunique(i).LT.huge(i)) THEN
            itempsort(i) = iunique(i)
         ELSE
            WRITE (*,*) 'ERROR - iunique(i).GE.',huge(i)
            CALL quit(0)
         ENDIF
      END DO
      WRITE (idisk1, ERR=100) tag('iunique')
      WRITE (idisk1, ERR=100) (itempsort(i), i=1, npart)

c--int*8

c--Default real
      tagarr(1:5) = (/'x','y','z','m','h'/)
      DO j = 1, jlen
         WRITE (idisk1, ERR=100) tag(tagarr(j))
         WRITE (idisk1, ERR=100) (xyzmh(j,isort(i)), i=1, npart)
      END DO
c      tagarr(1:4) = (/'vx','vy','vz','u '/)
c      DO j = 1, 4
c         WRITE (idisk1, ERR=100) tag(tagarr(j))
c         WRITE (idisk1, ERR=100) (vxyzu(j,isort(i)), i=1, npart)
c      END DO
c--real*4
      WRITE (idisk1, ERR=100) tag('rho')
      WRITE (idisk1, ERR=100) (rho(isort(i)), i=1, npart)
c--Dump smoothing length as a real*4 to save space
      DO i = 1, npart
         dq(i) = real(xyzmh(5,i),kind=4)
      END DO
      WRITE (idisk1, ERR=100) tag('h')
      WRITE (idisk1, ERR=100) (dq(isort(i)), i=1, npart)
c--Dump internal energy as a real*4 to save space
      DO i = 1, npart
         dq(i) = real(vxyzu(4,i),kind=4)
      END DO
      WRITE (idisk1, ERR=100) tag('u')
      WRITE (idisk1, ERR=100) (dq(isort(i)), i=1, npart)

c      WRITE (idisk1, ERR=100) tag('dgrav')
c      WRITE (idisk1, ERR=100) (dgrav(isort(i)), i=1, npart)
c      WRITE (idisk1, ERR=100) tag('alpha')
c     WRITE (idisk1, ERR=100) (alphaMM(1,isort(i)), i=1, npart)
c--real*8

c
c--Array length 2 arrays
c
c--Default int
      WRITE (idisk1, ERR=100) tag('listpm')
      WRITE (idisk1, ERR=100) (iorig(listpm(i)), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      WRITE (idisk1, ERR=100) tag('m')
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

      IF (number.GE.3) THEN
c
c--Array length 3 arrays
c
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real

c--real*4
         DO i = 1, npart
            dq(i) = real(ekcle(1,i),kind=4)
         END DO
         WRITE (idisk1, ERR=100) tag('Erad')
         WRITE (idisk1, ERR=100) (dq(isort(i)), i=1, npart)
         DO i = 1, npart
            dq(i) = real(ekcle(3,i),kind=4)
         END DO
         WRITE (idisk1, ERR=100) tag('Cv')
         WRITE (idisk1, ERR=100) (dq(isort(i)), i=1, npart)

c--real*8

      ENDIF
      IF (number.GE.4) THEN
c
c--Array length 4 arrays
c
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real

c--real*4
c
c  Dump B not B/rho
c
         DO j = 1, 3
            DO i = 1, npart
               dq(i) = real(Bxyz(j,i),kind=4)
            END DO
            WRITE (idisk1, ERR=100) tag('B'//labelx(j))
            WRITE (idisk1, ERR=100) (dq(isort(i)), i=1, npart)
         END DO

c--real*8

      ENDIF
c
c--End of writing small quantities
c
      ENDIF
c
c--End of writing block quantities
c
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
#ifdef MPICOPY
      ENDIF
#endif
c
c--End writing of dump file
c--------------------------
c

      RETURN
c
c--An error as occured while writing
c
 100  CALL error(where, 1)

      RETURN
      
      CONTAINS

      ELEMENTAL FUNCTION tag(label)
c*********************************************************************
c  contruct tag, padded or truncated to 16-characters based on input
c  string
c*********************************************************************
      IMPLICIT NONE
      CHARACTER*16 tag
      CHARACTER(LEN=*), INTENT(IN) :: label

      tag = adjustl(label)

      END FUNCTION tag

      END SUBROUTINE wdump

c---------------------------------------------------------------------

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
      CHARACTER*2 firstchar

      fileident = firstchar//'sphNG:'
      IF (nlmax.EQ.1) THEN
         fileident = fileident(1:8) // 'gradh=on, '
      ELSE
         fileident = fileident(1:8) // 'gradh=off,'
      ENDIF
      IF (imhd.EQ.idim) THEN
         fileident = fileident(1:18) // 'MHD='//varmhd(1:4)//','
      ELSE
         fileident = fileident(1:18) // 'MHD=none,'
      ENDIF
      IF (iradtrans.EQ.idim) THEN
         fileident = fileident(1:27) // 'RT=on '
      ELSE
         fileident = fileident(1:27) // 'RT=off'
      ENDIF

      END FUNCTION fileident
