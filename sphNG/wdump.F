      SUBROUTINE wdump(idisk1)
c************************************************************
c                                                           *
c  This routine ONLY writes a dump on disk                  *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif
      !IMPLICIT NONE
#ifdef USEKROME
      INCLUDE 'COMMONS/krome_mods'
#endif
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif

      INCLUDE 'idim'
      INCLUDE 'igrape'
#ifdef USEKROME
      INCLUDE 'COMMONS/kromevar'
#endif

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif
#ifdef MPIMOD
      INCLUDE 'COMMONS/mpimod'
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
      INCLUDE 'COMMONS/makeplt'
      INCLUDE 'COMMONS/radsink'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/dustfluid'
      INCLUDE 'COMMONS/dustfluidvelu'
      INCLUDE 'COMMONS/dustfluidgrains'
      INCLUDE 'COMMONS/tstop'
      INCLUDE 'COMMONS/tstopvar'
      INCLUDE 'COMMONS/dustimplicit'
      INCLUDE 'COMMONS/grain'
      INCLUDE 'COMMONS/compv'
      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09rho'
      INCLUDE 'COMMONS/f2HY09'
      INCLUDE 'COMMONS/HY09check'
      INCLUDE 'COMMONS/HY09accel'
      INCLUDE 'COMMONS/eosq'

      DIMENSION itempsort(idim)
c      EQUIVALENCE (itempsort, next1)

      DIMENSION dustfrac(ndusttypes,idim_dustFluid)
      DIMENSION dustfracsum1(idim_dustFluid)
      DIMENSION deltavsum(3,idim_dustFluid)
      DIMENSION dustproperties(4)

      CHARACTER*7 where
      CHARACTER*100 fileident
      CHARACTER*16 tagarr(10), tagwrite
      INTEGER*4 int1, int2, int1o
      INTEGER*8 number8
      DIMENSION nums(8)
      INTEGER iversion
      CHARACTER*1 labelx(3)
      PARAMETER (labelx = (/'x','y','z'/))
#ifdef USEKROME
      PARAMETER (maxmols=150)
      CHARACTER*16 mollist(maxmols)
      INTEGER dumpmols, index
      REAL*8 masses(krome_nmols), abundj(idim)
#endif

c
c--Standard numbers
c
      PARAMETER (int1=060769,int2=060878)
      PARAMETER (int1o=690706)

      DATA where/'wdump'/
#ifdef USEKROME
      dumpmols = 0
      if (usekrome .EQ. 1) then
         CALL readmolfile(maxmols,mollist,dumpmols)
      end if
#endif
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
#ifdef MPIMOD
      ifulldump = 0
#endif

#if defined MPIALL || MPIMOD
c  On MPI only iproc=0 writes the file ident. 
      IF (iproc.EQ.0) THEN
#endif
      WRITE (idisk1, ERR=100) int1,r1,int2,iversion,int1o
      IF (ifulldump.EQ.0) THEN
         WRITE (idisk1, ERR=100) fileident('FT')
      ELSE
         WRITE (idisk1, ERR=100) fileident('ST')
      ENDIF
#if defined MPIALL || MPIMOD
      ENDIF
#endif
c
c--Single values
c
#ifndef MPIMOD
      nparttot = npart
#endif
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
#ifdef MPIMOD
      nblocks = numproc
#endif

#if defined MPIALL || MPIMOD
c  On MPI only iproc=0 writes the headers. 
c  Other MPI processes skip forward to the lock, then write their data arrays
c  On MPICOPY only iproc=0 writes anything, the corresponding ENDIF is at the
c  end of the subroutine.
      IF (iproc.EQ.0) THEN
#endif
c
c--Default int
      number = 15 + NTAB
c      number = 14 + NTAB
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100) tag('nparttot'),tag('n1'),tag('n2'),
     &      tag('nreassigntot'),tag('naccretetot'),tag('nkilltot'),
     &      tag('nblocks'),tag('iyr'),tag('idum'),(tag('iv'),i=1,NTAB),
     &      tag('iplanetesimals'),tag('irotpot'),tag('idragscheme'),
     &      tag('HY09_dust_bins'),tag('idustFluid'),tag('ndusttypes')

      WRITE (idisk1, ERR=100) nparttot,n1,n2,nreassigntot,
     &     naccretetot,nkilltot,nblocks,iyr,idum,
     &     (iv(i),i=1,NTAB), iplanetesimals, irotpot,
     &     idragscheme, HY09_ndust_bins, idustFluid, ndusttypes
      print *,'Wrote to dump file ',HY09_ndust_bins
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
      number = 46

      IF (idustFluid.NE.0) THEN
         IF (ndusttypes.GT.1) THEN
            ndustheaders = 4
            tagarr(1:ndustheaders) =
     &           (/'smincgs     ','smaxcgs     ',
     &           'sindex      ','graindenscgs'/)
            dustproperties(1) = smincgs
            dustproperties(2) = smaxcgs
            dustproperties(3) = sindex
            dustproperties(4) = densgrain(1)*udens
         ELSE
            ndustheaders = 2
            tagarr(1:ndustheaders) =
     &         (/'grainsizecgs','graindenscgs'/)
            dustproperties(1) = sgrain(1)*udist
            dustproperties(2) = densgrain(1)*udens
         ENDIF
      ELSE
         ndustheaders = 0
      ENDIF
      number = number + ndustheaders
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
     &  tag('ymax'),tag('zmin'),tag('zmax'),tag('rmax'),tag('rcyl'),
     &  tag('HY09_size_min'),tag('HY09_size_max'),tag('HY09_cutoff'),
     &  tag('HY09_slope'),tag('HY09_dustdensity'),
     &  tag('HY09dustgasratio'),tag('HY09_vcoag_coeff'),
     &  (tag(tagarr(i)),i=1,ndustheaders)

      WRITE (idisk1, ERR=100) gt, dtmax, gamma, rhozero, RK2,
     &     escaptot, tkintot, tgravtot, ttermtot, 
     &     anglostxtot, anglostytot, anglostztot,
     &     specangtot, ptmassintot, tmagtot, Bextx, Bexty, Bextz,
     &     hzero, uzero_n2, hmass, gapfac, pmassinitial, sdprof,
     &     rorbit_orig, min_rplan, max_rplan, planetesimalmass,
     &     coremass_orig, coremass, hfact, xmin, xmax, ymin, ymax,
     &     zmin, zmax, rmax, rcyl,
     &     HY09_size_min, HY09_size_max, HY09_cutoff, HY09_slope, 
     &     HY09_dust_density, HY09_dustgas_ratio, HY09_vcoag_coeff,
     &     (dustproperties(i),i=1,ndustheaders)

c--real*4
      number = 0
      WRITE (idisk1, ERR=100) number
c--real*8
      number = 4
      WRITE (idisk1, ERR=100) number

      WRITE (idisk1, ERR=100)
     &      tag((/'udist ','umass ','utime ','umagfd'/))
      WRITE (idisk1, ERR=100) udist, umass, utime, umagfd
#if defined MPI || MPIMOD
c Note: this is MPI *not* MPIALL because for MPICOPY builds only iproc=0
c       writes to the output file. A correspending #ifdef MPICOPY ENDIF
c       is located at the end of this subroutine.
      ENDIF
#endif
c
c--Arrays
c
c
c--Number of array lengths
c     This has up to 4 blocks (hydro variables, sink particles, radiative
c     transfer, and MHD) repeated nblocks times for each of nblocks MPI
c     processes
c
c     This code runs on all MPI processes
      number = 2
      IF (iradtrans.EQ.idim) number = 3
      IF (imhd.EQ.idim) number = 4
      nblockarrays = number*nblocks

#if defined MPI || MPIMOD
c     Only write the total number of arrays on first process.
      IF (iproc.EQ.0) THEN
#endif
         WRITE (idisk1, ERR=100) nblockarrays
#if defined MPI || MPIMOD
      ENDIF
#endif

#ifdef MPI
      nowgo = 0
#ifdef MPIDEBUG
      print *,iproc,': About to write ',nprint,ifulldump
#endif
c     iproc=0 doesn't wait for the lock. All other processes wait until
c     the previous processes have written their arrays
      IF (iproc.GT.0) THEN
         CALL MPI_RECV(nowgo,1,MPI_INTEGER,iproc-1,99,MPI_COMM_WORLD,
     &        istatus,ierr)
#ifdef MPIDEBUG
         print *,iproc,': Opening file ',nprint
#endif
c     MPI processes call CLOSE after writing the arrays. 
c     file will then open the binary file with STATUS='append'
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

      IF (istellar_ionisation) nums(6) = nums(6) + 6

      IF (idustFluid.NE.0) THEN
         !--dustfrac, deltav, and tstop
         IF (iDumpDeltavTstop) THEN
            nums(6) = nums(6) + 5*ndusttypes
         ELSE
            nums(6) = nums(6) + 1*ndusttypes
         ENDIF

         !--dustfracsum and deltavsum
         IF (ndusttypes.GT.1) nums(6) = nums(6) + 4
      ENDIF
c
c--Pablo: Change to include dust grain evol
c
      IF (idim_grow.EQ.idim) nums(6) = nums(6) + 3

      IF (nlmax.EQ.1) THEN
         nums(7) = 4
      ELSE
         nums(7) = 3
      ENDIF
      IF (idimHY09.EQ.idim) THEN
         nums(6) = nums(6) + HY09_ndust_bins + 5
         nums(7) = nums(7) + 2
      ENDIF
      IF (idustIMPL.EQ.1 .AND. imakedust) THEN
         nums(6) = nums(6) + 3
         nums(7) = nums(7) + 1
      ENDIF
      nums(8) = 0
      WRITE (idisk1, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 2 header
c
      number8 = nptmass
      nums(1) = 1
      IF (istellarfeedback.EQ.6) nums(1) = nums(1) + 1
      nums(2) = 0
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 9
      IF (istellar_ionisation) nums(6) = nums(6) + 1
      IF (istellarfeedback.EQ.5 .OR. istellarfeedback.EQ.6)
     &     nums(6) = nums(6) + 6
      IF (istellarfeedback.EQ.6) nums(6) = nums(6) + 4
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
#ifdef USEKROME
         number8 = nprint
c         nums(4) = nums(4) + 1
         IF (usekrome .EQ. 1) THEN
            nums(6) = nums(6) + dumpmols + 2
         ELSE IF (usekrome .EQ. 2 .OR. usekrome.EQ.3) THEN
            nums(6) = nums(6) + krome_nmols + 2
         END IF
#endif
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
c        nums(7) = 13
#ifdef NONIDEAL
         nums(7) = nums(7) + 4
#endif
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

      IF (istellar_ionisation) THEN
         WRITE (idisk1, ERR=100) tag('HIIFrac')
         WRITE (idisk1, ERR=100) (HIIion(isort(i)), i=1, 
     &        nprint)
         WRITE (idisk1, ERR=100) tag('dHIIdt')
         WRITE (idisk1, ERR=100) (dHIIdt(1,isort(i)), i=1, 
     &           nprint)
         tagarr(1:4) = tag((/'StellarExt  ','StellarCol  ',
     &        'StellarH2Col','StellarHII  '/))
         DO j = 1, 4
            WRITE (idisk1, ERR=100) tagarr(j)
            WRITE (idisk1, ERR=100) 
     &           (stellarrad(j,1,isort(i)),i=1,nprint)
         END DO
      ENDIF

      IF (idustFluid.NE.0) THEN
         dustfrac     = 0.
         DO k = 1, ndusttypes
            IF (idustFluid.EQ.1) THEN
               !--sqrt(rho*epsilon) method
               dustfrac(k,1:nprint) = dustvar(k,isort(1:nprint))**2
     &              /rho(isort(1:nprint))
            ELSEIF (ABS(idustFluid).EQ.2) THEN
               !--sqrt(epsilon/(1-epsilon)) method
               dustfrac(k,1:nprint) = dustvar(k,isort(1:nprint))**2/
     &              (1.0+dustvar(k,isort(1:nprint))**2)
            ELSEIF (idustFluid.EQ.3) THEN
               !--asin(sqrt(epsilon)) method
               dustfrac(k,1:nprint) = SIN(dustvar(k,isort(1:nprint)))**2
            ELSEIF (idustFluid.EQ.4) THEN
               !--epsilon method
               dustfrac(k,1:nprint) = dustvar(k,isort(1:nprint))
            ELSE
               WRITE (*,*) 'ERROR - idustFluid not valid ',idustFluid
               CALL quit(0)
            ENDIF
         END DO
c
c--If BOTH idustFluid and idimHY09, then use the dust set up from HY09
c
         IF (idimHY09.EQ.idim .AND. gt.EQ.0.) THEN
            dustfrac(:,1:nprint) = HY09bin_rho(:,isort(1:nprint))
         ENDIF
c
c--Note that dustfracsum is dumped to allow the dust and gas velocities
c     to be computed, but it is not used in the integration at all.
c     NB:  Not using isort() as precomputed above.
c
         IF (ndusttypes.GT.1) THEN
            WRITE (idisk1, ERR=100) tag('dustfractot')
            WRITE (idisk1, ERR=100) (SUM(dustfrac(:,i)),i=1,nprint)
         ENDIF
c
c--Note dustfrac (epsilon) is dumped, but dustvar(eg. SQRT(epsilon*rho))
c     is the variable that is integrated.  NB: Not using isort() as 
c     precomputed above.
c
         DO k = 1, ndusttypes
            WRITE (idisk1, ERR=100) tag('dustfrac')
            WRITE (idisk1, ERR=100) (dustfrac(k,i),i=1,nprint)
         END DO
c
c--Note that deltavsum and dustfluidvxyz
c     are dumped to allow the dust and gas velocities to be computed,
c     but they are not used in the integration at all.
c
         IF (ndusttypes.GT.1) THEN
            dustfracsum1   = 0.
            deltavsum      = 0.
            tagarr(1:3)    = (/'deltavsumx','deltavsumy','deltavsumz'/)
            WHERE (SUM(dustfrac(:,1:nprint),1).GT.0.)
               dustfracsum1(1:nprint) = 1./SUM(dustfrac(:,1:nprint),1)
            END WHERE
            DO i = 1,nprint
               deltavsum(1,i) = dustfracsum1(i)*SUM(dustfrac(:,i)
     &                          *dustfluidvxyz(1,:,isort(i)))
               deltavsum(2,i) = dustfracsum1(i)*SUM(dustfrac(:,i)
     &                          *dustfluidvxyz(2,:,isort(i)))
               deltavsum(3,i) = dustfracsum1(i)*SUM(dustfrac(:,i)
     &                          *dustfluidvxyz(3,:,isort(i)))
            ENDDO
            DO j = 1, 3
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) (deltavsum(j,i),i=1,nprint)
            ENDDO
         ENDIF
c
         IF (iDumpDeltavTstop) THEN
            tagarr(1:3) = (/'deltavx','deltavy','deltavz'/)
            DO k = 1, ndusttypes
               DO j = 1, 3
                  WRITE (idisk1, ERR=100) tag(tagarr(j))
                  WRITE (idisk1, ERR=100) (dustfluidvxyz(j,k,isort(i)),
     &                 i=1,nprint)
               END DO
            END DO
c
c--Note that stoppingtime
c     is dumped for analysis, debugging, and/or comparision
c
            DO k = 1, ndusttypes
               WRITE (idisk1, ERR=100) tag('tstop')
               WRITE (idisk1, ERR=100)
     &              (stoppingtime(k,isort(i)),i=1,nprint)
            END DO
         ENDIF
      ENDIF
c
c--Outputs gas pressure gradient so that it can be used to set dust
c     at terminal velocity initially for two-fluid dust.
c
      IF (idustIMPL.EQ.1 .AND. imakedust) THEN
         tagarr(1:3) = (/'pgradx','pgrady','pgradz'/)
         DO j = 1, 3
            WRITE (idisk1, ERR=100) tag(tagarr(j))
            WRITE (idisk1, ERR=100) (dustfluidvxyz(j,1,isort(i)),
     &           i=1,nprint)
         END DO         
      ENDIF
c
c--Pablo: Change to include dust grain evol
c
      IF (idim_grow.EQ.idim) THEN
         WRITE (idisk1, ERR=100) tag('rgrain')
         WRITE (idisk1, ERR=100) (rgrain(isort(i)),i=1,nprint)

         WRITE (idisk1, ERR=100) tag('rhoG')
         WRITE (idisk1, ERR=100) (rhoG(isort(i)),i=1,nprint)

         WRITE (idisk1, ERR=100) tag('vcomp')
         WRITE (idisk1, ERR=100) (vcomp(isort(i)),i=1,nprint)
      ENDIF

      IF (iplotdump.EQ.1) THEN
         tagarr(1:2) = (/'gas ','dust'/)
         DO j = 1, 2
            WRITE (idisk1, ERR=100) tag(tagarr(j))
            WRITE (idisk1, ERR=100) (xtype(j,isort(i)), i=1, nprint)
         END DO     
      ENDIF

      IF (idimHY09.EQ.idim) THEN
         DO j = 1, HY09_ndust_bins
c
c--Name dust bin sizes in microns
c
            size_mu = HY09binsizes(j)*1.0E+04
            IF (size_mu.LT.0.01) THEN
               WRITE (tagwrite,99001) size_mu
99001          FORMAT ('Dust:',F6.4)
            ELSEIF (size_mu.LT.0.1) THEN
               WRITE (tagwrite,99002) size_mu
99002          FORMAT ('Dust:',F5.3)
            ELSEIF (size_mu.LT.1.0) THEN
               WRITE (tagwrite,99003) size_mu
99003          FORMAT ('Dust:',F4.2)
            ELSEIF (size_mu.LT.10.0) THEN
               WRITE (tagwrite,99004) size_mu
99004          FORMAT ('Dust:',F3.1)
            ELSEIF (size_mu.LT.100.0) THEN
               WRITE (tagwrite,99005) INT(size_mu)
99005          FORMAT ('Dust:',I2) 
            ELSEIF (size_mu.LT.1000.0) THEN
               WRITE (tagwrite,99006) INT(size_mu)
99006          FORMAT ('Dust:',I3) 
            ELSEIF (size_mu.LT.10000.0) THEN
               WRITE (tagwrite,99007) INT(size_mu)
99007          FORMAT ('Dust:',I4) 
            ELSEIF (size_mu.LT.1.0E+05) THEN
               WRITE (tagwrite,99008) INT(size_mu)
99008          FORMAT ('Dust:',I5) 
            ELSE
               WRITE (tagwrite,99009) INT(size_mu)
99009          FORMAT ('Dust:',I6) 
            ENDIF
            WRITE (idisk1, ERR=100) tag(tagwrite)
            WRITE (idisk1, ERR=100) (HY09bin_rho(j,isort(i)),i=1,nprint)
         END DO

         WRITE (idisk1, ERR=100) tag('hfrac')
         WRITE (idisk1, ERR=100) (HY09frac_h_moved(isort(i)),i=1,nprint)
c
c--Output the gravitational acceleration in code units.  Only required
c     for diagnostic/testing purposes.
c
         tagarr(1:3) = (/'gravaccelx','gravaccely','gravaccelz'/)
         DO j = 1, 3
            WRITE (idisk1, ERR=100) tag(tagarr(j))
            WRITE (idisk1, ERR=100) (grav_accel(j,isort(i)),
     &           i=1,nprint)
         END DO
c         tagarr(1:3) = (/'gasaccelx','gasaccely','gasaccelz'/)
c         DO j = 1, 3
c            WRITE (idisk1, ERR=100) tag(tagarr(j))
c            WRITE (idisk1, ERR=100) (gas_accel(j,isort(i)),
c     &           i=1,nprint)
c         END DO
c
c--Output the estimate of the relative grain speed relative to the gas
c     for dust grains at the peak of the mass distribution.
c     This is stored in f1HY09vrel_h as the speed/h in code units
c     so to get the speed in cgs units need to multiply by h and
c     convert to cgs units.
c
         WRITE (idisk1, ERR=100) tag('dustvrel')
         WRITE (idisk1, ERR=100) (f1HY09vrel_h(isort(i))*
     &        xyzmh(5,isort(i))*udist/utime,i=1,nprint)
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

      IF (idustIMPL.EQ.1 .AND. imakedust) THEN
         WRITE (idisk1, ERR=100) tag('vsound')
         WRITE (idisk1, ERR=100) (vsound(isort(i)), i=1, nprint)
      ENDIF

      IF (idimHY09.EQ.idim) THEN
c
c--Determine the grain size associated with the peak of the grain mass
c     distribution (in microns)
c
         DO i = 1, nprint
            ilocation = MAXLOC(HY09bin_rho(:,isort(i)), DIM=1 )
            HY09peak(i) = HY09binsizes(ilocation)*1.0E+04
         END DO
         WRITE (idisk1, ERR=100) tag('dustmax')
         WRITE (idisk1, ERR=100) (HY09peak(i), i=1, nprint)
c
c--Output flag for whether particle is disc (1) or envelope (0)
c
         WRITE (idisk1, ERR=100) tag('discflag')
         WRITE (idisk1, ERR=100) (HY09discflag(isort(i)), i=1, nprint)
      ENDIF
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

      IF (istellarfeedback.EQ.6) THEN
         WRITE (idisk1, ERR=100) tag('listlyon')
         WRITE (idisk1, ERR=100) (listlyon(i), i=1,nptmass)
      ENDIF
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

      IF (istellar_ionisation) THEN
         WRITE (idisk1, ERR=100) tag('ionflux')
         WRITE (idisk1, ERR=100) (stellar_radtempion(3,i), i=1,nptmass)
         print *,(stellar_radtempion(3,i), i=1,nptmass)
      ENDIF

      IF (istellarfeedback.EQ.5 .OR. istellarfeedback.EQ.6) THEN
         WRITE (idisk1, ERR=100) tag('ptmassform')
         WRITE (idisk1, ERR=100) (ptmassform(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('diskinit')
         WRITE (idisk1, ERR=100) (diskinit(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('stellartime')
         WRITE (idisk1, ERR=100) (stellartime(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('stellarmass')
         WRITE (idisk1, ERR=100) (stellarmass(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('stellarradius')
         WRITE (idisk1, ERR=100) (stellarradius(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('stellarlum')
         WRITE (idisk1, ERR=100) (stellarluminosity(i),i=1,nptmass)
      ENDIF

      IF (istellarfeedback.EQ.6) THEN
         WRITE (idisk1, ERR=100) tag('lyontime')
         WRITE (idisk1, ERR=100) (lyontime(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('lyonmass')
         WRITE (idisk1, ERR=100) (lyonmass(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('lyonradius')
         WRITE (idisk1, ERR=100) (lyonradius(i),i=1,nptmass)
         WRITE (idisk1, ERR=100) tag('lyonlum')
         WRITE (idisk1, ERR=100) (lyonluminosity(i),i=1,nptmass) 
      ENDIF

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
         tagarr(1:4) = tag((/'Erad   ','opacity','Cv     ','lambda '/))
         tagarr(5)   = tag('EddingtonFactor')
         DO j = 1, 5
            WRITE (idisk1, ERR=100) tag(tagarr(j))
            WRITE (idisk1, ERR=100) (ekcle(j,isort(i)), i=1, nprint)
         END DO
         IF (idustRT.GT.0) THEN
            tagarr(1:2) = tag((/'DustTemperature','DustOpacity    '/))
            DO j = 1, 2
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) (dust_tk(j,isort(i)),i=1,nprint)
            END DO
         ENDIF
#ifdef USEKROME
         masses = krome_get_mass()
         IF(usekrome.EQ.1) THEN
c  -     dump selected abundances /H2 abundance             
            DO j=1, dumpmols
               index = krome_get_index(adjustl(mollist(j)))
c               DO i=1,nprint
c                  IF (mfracs(krome_idx_H2,isort(i)) .GT. 0d0) THEN
c                     abundj(isort(i)) = mfracs(index,isort(i))*
c     &                  masses(krome_idx_H2)/
c     &                  (masses(index) * mfracs(krome_idx_H2,isort(i)))
c                 ELSE
c                     abundj(isort(i)) = 0d0
c                 END IF
c               END DO
c     - I'm now going to dump mass fractions not abundances /H2
                  WRITE (idisk1,ERR=100) tag(molnames(index))
                  WRITE (idisk1,ERR=100) 
     &                 (mfracs(index,isort(i)),i=1,nprint)
c               WRITE (idisk1,ERR=100) (abundj(isort(i)),i=1,nprint)
               END DO
            ELSE IF (usekrome.EQ.2 .OR. usekrome.EQ.3) THEN
c -      dump mass fractions of all species
            PRINT *, "Dumping all mass fractions"
            DO j=1, krome_nmols
               WRITE (idisk1,ERR=100) tag(molnames(j))
               WRITE (idisk1,ERR=100) (mfracs(j,isort(i)),i=1,nprint)
            END DO
         ENDIF
         WRITE (idisk1,ERR=100) tag('tlastconv')
         WRITE (idisk1,ERR=100) (tlastconv(isort(i)),i=1,nprint)
#ifdef USEKROME
         print *, "writing numfails", usekrome
         IF (usekrome .GT. 0) THEN  
            WRITE (idisk1, ERR=100) tag('numfails')
            WRITE (idisk1, ERR=100) 
     &           (real(numfails(isort(i))),i=1,nprint)
         END IF
#endif

         PRINT *, "Written chemistry"
#endif

c--real*4
         IF (idustRT.GT.0) THEN
            tagarr(1:3) = tag((/'Cplus  ','AtomicC','CO     '/))
            DO j = 1, nchemistry
               WRITE (idisk1, ERR=100) tag(tagarr(j))
               WRITE (idisk1, ERR=100) (chemistry(j,isort(i)), i=1, 
     &              nprint)
            END DO

            tagarr(1:4) = tag((/'HeatingISR     ','Gphotoelectric ',
     &           'Extinction     ','H2SelfShielding'/))
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
           WRITE (idisk1, ERR=100) (divcurlB(j,isort(i)), i=1, nprint)
          ELSE
           WRITE (idisk1, ERR=100) tag('curlB'//labelx(j-1))
#ifdef NONIDEAL
           WRITE (idisk1, ERR=100) (real(jcurrent(j-1,isort(i)),kind=4),
     &            i=1, nprint)
#else
           WRITE (idisk1, ERR=100) (divcurlB(j,isort(i)), i=1, nprint)
#endif
          ENDIF
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
#ifdef NONIDEAL
         tagarr(1:4)=(/'eta_{OR}','eta_{HE}','eta_{AD}','ne/n    '/)
         DO j = 1,4
           WRITE(idisk1, ERR=100) tag(tagarr(j))
           WRITE(idisk1, ERR=100) (real(eta_nimhd(j,isort(i)),kind=4),
     &          i=1, nprint)
         END DO
#endif
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
c
c--Pablo: Change to include dust grain evol
c
      IF (idim_grow.EQ.idim) THEN
        nums(6) = jlen+3
      ELSE
        nums(6) = jlen
      ENDIF
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
c
c--Pablo: Change to include dust grain evol
c
      IF (idim_grow.EQ.idim) THEN
         WRITE (idisk1, ERR=100) tag('rgrain')
         WRITE (idisk1, ERR=100) (rgrain(isort(i)), i=1, npart)

         WRITE (idisk1, ERR=100) tag('rhoG')
         WRITE (idisk1, ERR=100) (rhoG(isort(i)), i=1, npart)

         WRITE (idisk1, ERR=100) tag('vcomp')
         WRITE (idisk1, ERR=100) (vcomp(isort(i)),i=1,npart)
      ENDIF

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
#ifdef MPIMOD
      IF (iproc+1 .EQ. numproc) THEN
#endif
         CALL FLUSH (idisk1)
         CLOSE (idisk1)
#ifdef MPIMOD
      ENDIF
#endif
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
      IF (idustFluid.NE.0) THEN
         fileident = fileident(1:34) // '1dust'
      ENDIF

      END FUNCTION fileident
