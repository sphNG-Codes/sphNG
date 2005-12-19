      SUBROUTINE preset
c************************************************************
c                                                           *
c  This routine makes sure that everything not entered as   *
c     option is defined before integration begins.          *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/stop'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/binfile'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/optbl'
      INCLUDE 'COMMONS/tgtbl'
      INCLUDE 'COMMONS/mutbl'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/cgas'

      CHARACTER*7 where
      CHARACTER*21 ptfile, accfile, killfile, reassfile

      DATA where/'preset'/
c
c--Open file1
c
      OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
c
c--Find correct dump to start
c
      CALL place(idisk1, ipos, irec, 1)
c
c--Read dump
c
      CALL rdump(idisk1, ichkl)
      IF (ichkl.EQ.1) CALL error(where, ichkl)
c
c--Open point mass data output file
c
      WRITE (ptfile,99990) namerun
99990 FORMAT ('P', A20)
      WRITE (accfile,99991) namerun
99991 FORMAT ('A', A20)
      IF (nptmass.NE.0 .OR. iptmass.NE.0) THEN
         OPEN (iptprint, FILE=ptfile,STATUS='unknown',
     &        FORM='unformatted')
         OPEN (iaccpr,FILE=accfile,STATUS='unknown',
     &        FORM='unformatted')
      ENDIF
c
c--Open files for killing and reassignment of particles
c
      WRITE (killfile,99992) namerun
99992 FORMAT ('K', A20)
      WRITE (reassfile,99993) namerun
99993 FORMAT ('R', A20)
      IF (ibound.EQ.8 .OR. ibound/10.EQ.9 .OR. ibound.EQ.100) THEN
         OPEN (ikillpr, FILE=killfile, STATUS='unknown',
     &        FORM='unformatted')
         OPEN (ireasspr, FILE=reassfile, STATUS='unknown',
     &        FORM='unformatted')
      ENDIF
c
c--Open notify file
c
      OPEN (inotify,FILE='notify')
c
c--Initialize viscosity switch
c
      DO i = 1, npart
         IF (ifsvi.EQ.6 .AND. gt.EQ.0.0) alphaMM(i) = alphamin
      END DO
c
c--Set constant for artificial viscosity
c
      IF (alpha.EQ.0. .AND. beta.EQ.0.) THEN
         alpha = 1.0
         beta = 2.0
      ENDIF
c
c--Set accuracy parameter for tree force calculation
c     theoretical limit for 3D tree is 0.57 = 1/SQRT(3)
c
c      acc = 0.7
      acc = 0.5
c      acc = 0.3
c      acc = 0.0
c
c--Set stop flag
c
      istop = 0
c
c--Set min and max limit of neighbours the code tries to inforce
c
      neimin = 30
      neimax = 70
      nrange = 12
c
c--Total mass
c
      pmassmin = 1.0E+10
      fmas1 = 0.
      DO i = 1, n1
         IF (iphase(i).GE.0) THEN
            fmas1 = fmas1 + xyzmh(4,i)
            IF (pmassmin.GT.xyzmh(4,i)) pmassmin = xyzmh(4,i)
         ENDIF
      END DO
      fmas2 = 0.
      DO i = n1 + 1, npart
         IF (iphase(i).GE.0) THEN
            fmas2 = fmas2 + xyzmh(4,i)
            IF (pmassmin.GT.xyzmh(4,i)) pmassmin = xyzmh(4,i)
         ENDIF
      END DO
c
c--Set mass for particle being partially accreted at which it is 
c     completely accreted
c
      pmassleast = pmassmin/100.
c
c--Point Mass Presets
c
c--Set critical density for point mass creation
c
      rhocrea = rhozero*ptmcrit
c
c--Standardise point mass types and accretion radii
c
      DO i = 1, nptmass
         iptcur = listpm(i)
         IF (initialptm.LT.1 .OR. initialptm.GT.4) CALL error(where, 2)
c         iphase(iptcur) = initialptm
cccc         xyzmh(5,iptcur) = hacc
      END DO
c
c--For Massive Accretion, only accrete mass and change ptmass properties
c     for certain particles - exclude these:
c
      DO i = 1, idim
         notacc(i) = .FALSE.
      END DO
      nnotacc = 0
c      WRITE (iprint,*) '********* IDELAYACC = 1 *********'
c      OPEN (22, FILE='WANTRUN')
c 50   READ (22, *, END=100) inum
c      nnotacc = nnotacc + 1
c      notacc(inum) = .TRUE.
c      GOTO 50
c 100  CLOSE(22)
c
c--Gravitational softening for ptmass-ptmass interactions
c
      iptsoft = 1
      ptsoft = 1.000E-05
c
c--Set psoft for softening the gravitational potential when the 1/(r+psoft)
c     potential law is used (i.e. when igrape=1, or isoft=1).
c
      psoft = 0.01
c
c--Value of minimum h in order to save computing time (if hmin is different
c     from 0 then program does not follow high density regions accurately)
c
      hmin = 0.
c      hmin = 0.01

      IF (isoft.EQ.1 .AND. hmin.LT.psoft) hmin = psoft
c
c--Compute tables for kernel quantities
c
      CALL ktable
c
c--Load in tables of opacity and cv
c
      IF (encal.EQ.'r') THEN
         OPEN(UNIT=8,FILE='/home/mbate/tables/opacitytbl')
         DO i=1, opmxtg
            READ(8,*) (optable(i,j), j=1, opmxrh)
         ENDDO
         CLOSE(8)
      
         OPEN(UNIT=8,FILE='/home/mbate/tables/gasttbl',
     &        FORM='unformatted')
         DO i=1, tgmxu
            READ(8) (tgtable(i,j), j=1, tgmxrh)
         ENDDO
         CLOSE(8)
 
         OPEN(UNIT=8,FILE='/home/mbate/tables/molmasstbl',
     &        FORM='unformatted')
         DO i=1, mumxu
            READ(8) (mutable(i,j), j=1, mumxrh)
         ENDDO
         CLOSE(8)
      ENDIF
c
c--Build table for choosing inflow into planet/disc simulation
c 
      IF (ibound.EQ.100) CALL cumradtable
c
c--Build table for choosing inflow into planet/disc simulation from ZEUS
c     output
c 
      IF (ibound.EQ.100) CALL zeustable
 
      IF (itrace.EQ.'all') WRITE(iprint,250)
 250  FORMAT(' exit subroutine preset')

      RETURN
      END
