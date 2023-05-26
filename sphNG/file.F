      SUBROUTINE file
c************************************************************
c                                                           *
c  This subroutine handles the files. A new file is         *
c     created if the size of the current file gets too      *
c     large.                                                *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif

      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/planetesimal'

      LOGICAL ifirst
      DATA ifirst/.TRUE./

c      DATA maxsize/10000000/

      IF (ifirst) THEN
         ifirst = .FALSE.
c         GOTO 100
      ENDIF
c
c--Check for file status
c
c      CALL statfile(file1, ifsize)
c
c--If smaller than maxsize return
c
c      IF (ifsize.LT.maxsize) RETURN
c
c--If not create new file and close old
c
c      CLOSE (idisk1)

c  Do not increment the name on initialization. Needed for MPI set ups to match
c     non-MPI.
      IF (job(1:7).NE.'initial' .AND. job(1:8).NE.'makedust') THEN
         CALL increment_dump_filename(file1)
      ENDIF

#ifdef MPI
c      CALL appendiproc(7,file1,iproc)
      IF (iproc.EQ.0) THEN
         OPEN (idisk1, FILE=file1, FORM='unformatted', RECL=imaxrec)
      ELSE
         OPEN (idisk1, FILE=file1, FORM='unformatted',
     &        POSITION='APPEND', RECL=imaxrec)
      ENDIF
#else
#ifdef MPICOPY
      IF (iproc.EQ.0) THEN
#endif
      OPEN (idisk1, FILE=file1, FORM='unformatted', RECL=imaxrec)
#ifdef DGRID
c
c--Open dragfile for writing drag energy table.
c
      IF (gasdrag) THEN
11403    FORMAT (A4, I4.4)
         WRITE (dragfile, 11403) 'drag', inum
         OPEN (idragfile, FILE=dragfile, FORM='unformatted',
     &        RECL=imaxrec)
      ENDIF
#endif
#ifdef MPICOPY
      ENDIF
#endif
#endif

c
c--Write witness on listing
c
#ifdef MPIALL
      IF (iproc.EQ.0) THEN      
#endif
      WRITE (iprint, 99003)
99003 FORMAT (/)
      WRITE (iprint, 99004)
99004 FORMAT (' ------------------------------------------------')
      WRITE (iprint, 99005) file1
99005 FORMAT (' ------------ file ', A7, ' created --------------')
      WRITE (iprint, 99004)
      WRITE (iprint, 99003)
#ifdef MPIALL
      ENDIF
#endif

c 100  CONTINUE
      RETURN
      END
c
c====================================================================
c
      SUBROUTINE increment_dump_filename(file1)

      CHARACTER*11 file1
      CHARACTER*5 name

      IF ((file1(4:4).GE.'0') .AND. (file1(4:4).LE.'9') .AND. 
     &     (file1(5:5).GE.'0') .AND. (file1(5:5).LE.'9')) THEN
         READ (file1, 99303) name, inum
99303    FORMAT (A3, I4)
         inum = inum + 1
         iadd = 0
         IF (inum.LE.9) THEN
            WRITE (file1, 99311) name, iadd, iadd, iadd, inum
99311       FORMAT (A3, I1, I1, I1, I1)
         ELSEIF (inum.LE.99) THEN
            WRITE (file1, 99312) name, iadd, iadd, inum
99312       FORMAT (A3, I1, I1, I2)
         ELSEIF (inum.LE.999) THEN
            WRITE (file1, 99313) name, iadd, inum
99313       FORMAT (A3, I1, I3)
         ELSE
            WRITE (file1, 99303) name, inum
         ENDIF
      ELSEIF ((file1(5:5).GE.'0') .AND. (file1(5:5).LE.'9')) THEN
         READ (file1, 99403) name, inum
99403    FORMAT (A4, I3)
         inum = inum + 1
         iadd = 0
         IF (inum.LE.9) THEN
            WRITE (file1, 99411) name, iadd, iadd, inum
99411       FORMAT (A4, I1, I1, I1)
         ELSEIF (inum.LE.99) THEN
            WRITE (file1, 99412) name, iadd, inum
99412       FORMAT (A4, I1, I2)
         ELSE
            WRITE (file1, 99403) name, inum
         ENDIF
      ELSE
         READ (file1, 99502) name, inum
99502    FORMAT (A5, I2)
         inum = inum + 1
         IF (inum.LE.9) THEN
            iadd = 0
            WRITE (file1, 99511) name, iadd, inum
99511       FORMAT (A5, I1, I1)
         ELSE
            WRITE (file1, 99502) name, inum
         ENDIF
      ENDIF

      END
