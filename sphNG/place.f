      SUBROUTINE place(idisk1, ipos, irec, iflag)
c************************************************************
c                                                           *
c  This routine puts the pointer at beginning of dump ipos  *
c                                                           *
c************************************************************

      CHARACTER*7 where

      DATA where/'place'/

      irec = 0
      IF (ipos.NE.9999) THEN
c
c--Go to beginning of dump ipos
c
         DO i = 1, ipos - 1
            READ (idisk1, END=300)
            irec = i
         END DO
         irec = irec + 1
         RETURN
      ELSE
         DO i = 1, 10000
            READ (idisk1, END=200)
            irec = i
         END DO
      ENDIF
c
c--Go back one reccord
c
 200  BACKSPACE idisk1

      IF (iflag.EQ.1) BACKSPACE idisk1
      GOTO 400

 300  CALL error(where, 1)

 400  RETURN
      END
