      SUBROUTINE file
c************************************************************
c                                                           *
c  This subroutine handles the files. A new file is         *
c     created if the size of the current file gets too      *
c     large.                                                *
c                                                           *
c************************************************************
#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/actio'

      LOGICAL ifirst
      DATA ifirst/.TRUE./
      CHARACTER*5 name

      DATA maxsize/10000000/

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
      CLOSE (idisk1)
      IF ((file1(5:5).GE.'0') .AND. (file1(5:5).LE.'9')) THEN
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

#ifdef MPI
      CALL appendiproc(7,file1,iproc)
#endif

      OPEN (idisk1, FILE=file1, FORM='unformatted', RECL=imaxrec)
c
c--Write witness on listing
c
#ifdef MPI
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
#ifdef MPI
      ENDIF
#endif
c
c--Reset dump counter
c
      irec = 0

 100  RETURN
      END
