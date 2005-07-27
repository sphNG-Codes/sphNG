      SUBROUTINE newrun
c************************************************************
c                                                           *
c  This subroutine is the driver for the set up of a new    *
c     simulation.                                           *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/task'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/part'

      CHARACTER*7 where

      DATA where/'newrun'/
c
c--General set up
c
      WRITE (*, 99001)
99001 FORMAT (' GENERAL SET UP', //)
      WRITE (*, 99002)
99002 FORMAT (' give name of the run (20 char.)')
      READ (*, 99003) namerun
      WRITE (namenextrun,99003) namerun
99003 FORMAT (A20)
c
c--Write label of run
c
      CALL labrun

      WRITE (*, 99004)
99004 FORMAT (' type of initialization?', /,
     &        ' start from scratch         : scratch (s)', /,
     &        ' using existing files       : exist (e)' )
      READ (*, 99005) what
99005 FORMAT (A7)
c
c--Start from scratch
c
      IF (what.EQ.'scratch' .OR. what.EQ.'s') THEN
         WRITE (*, 99006)
99006    FORMAT (' name of binary file (7 char. max)')
         READ (*, 99005) file1
         WRITE (*,*) 'MAXIMUM RECORD LENGTH = ',imaxrec,' ',file1
         OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
         CALL setpart
         CALL inform(where)
         DO i = 1, npart
            isort(i) = i
            iorig(i) = i
         END DO
         CALL wdump(idisk1)
         CLOSE (idisk1)
      ENDIF
c
c--Use existing dumps to create new dump
c
      IF (what(1:5).EQ.'exist' .OR. what.EQ.'e') THEN
         WRITE (*, 99007)
99007    FORMAT (' name of the first binary file?')
         READ (*, 99005) file1
         WRITE (*, 99008)
99008    FORMAT (' name of the second binary file?')
         READ (*, 99005) file2
         WRITE (*, 99009)
99009    FORMAT (' name of the resulting binary file?')
         READ (*, 99005) file3
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
         CALL wdump(idisk3)
         CLOSE (idisk3)
         CLOSE (idisk2)
         CLOSE (idisk1)
      ENDIF
c
c--Write output
c
      CALL header(where)
      CALL prout(where)
c
c--Terminate run
c
      CALL endrun

      RETURN
      END
