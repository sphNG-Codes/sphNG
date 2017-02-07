      SUBROUTINE newrun
c************************************************************
c                                                           *
c  This subroutine is the driver for the set up of a new    *
c     simulation.                                           *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif
      INCLUDE 'idim'

#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/task'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/tming'

      CHARACTER*7 where

      DATA where/'newrun'/
#ifdef MPI
c
c--IF MPI setup, then read set up input from file
c
c      iread = 10
      OPEN (iread, FILE='setup.txt', STATUS='old',
     &      ACTION='read')
#endif
c
c--General set up
c
      WRITE (*, 99001)
99001 FORMAT (/, ' GENERAL SET UP', /)
      WRITE (*, 99002)
99002 FORMAT (' give name of the run (20 char.)')
      READ (iread, 99003) namerun
      WRITE (namenextrun,99003) namerun
99003 FORMAT (A20)
      WRITE (*,*)
c
c--Write label of run
c
      CALL labrun

      WRITE (*, 99004)
99004 FORMAT (/, ' type of initialization?', /,
     &        ' start from scratch         : scratch (s)', /,
     &        ' using existing files       : exist (e)' )
      READ (iread, 99005) what
99005 FORMAT (A7)
c
c--Start from scratch
c
      IF (what.EQ.'scratch' .OR. what.EQ.'s') THEN
         WRITE (*, 99006)
99006    FORMAT (' name of binary file (7 char. max)')
         READ (iread, 99005) file1
         WRITE (*,*) 'MAXIMUM RECORD LENGTH = ',imaxrec,' ',file1
#ifdef MPI
         IF (iproc.EQ.0) CALL file
#else
         OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
#endif
         CALL setpart
         CALL inform(where)
         nfullstep = 1
         CALL wdump_wrapper(idisk1)
c         CLOSE (idisk1)
      ENDIF
c
c--Use existing dumps to create new dump
c
      IF (what(1:5).EQ.'exist' .OR. what.EQ.'e') THEN
         WRITE (*, 99007)
99007    FORMAT (' name of the first binary file?')
         READ (iread, 99005) file1
         WRITE (*, 99008)
99008    FORMAT (' name of the second binary file?')
         READ (iread, 99005) file2
         WRITE (*, 99009)
99009    FORMAT (' name of the resulting binary file?')
         READ (iread, 99005) file3
         WRITE(*,*) 'MAXIMUM RECORD LENGTH = ',imaxrec
         OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         OPEN (idisk2, FILE=file2, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         OPEN (idisk3, FILE=file3, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         CALL addump
         CALL inform(where)
         DO i = 1, npart
            isort(i) = i
            iorig(i) = i
         END DO
         nfullstep = 1
         CALL wdump_wrapper(idisk3)
         CLOSE (idisk3)
         CLOSE (idisk2)
         CLOSE (idisk1)
      ENDIF
c
c--Write output
c
      CALL header(where)
      CALL prout(where)
#ifdef MPI
      CLOSE(iread)
#endif
c
c--Terminate run
c
      CALL endrun

      RETURN
      END
