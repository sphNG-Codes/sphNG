      SUBROUTINE secmes
c************************************************************
c                                                           *
c  This subroutine handles the messages received from the   *
c     operator.                                             *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/stop'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/secret'

      CHARACTER*20 string(10)
      CHARACTER*5 dummy5
      CHARACTER*3 dummy3
c
c--Open message file first
c
      iwait = 0
      OPEN (idisk2, FILE='secret', FORM='formatted')
c
c--Read message file (one mes. per line , max=10)
c
      nmes = 0
      DO i = 1, 10
         READ (idisk2, 99000, END=200) string(i)
99000    FORMAT (A20)
         nmes = i
      END DO
c
c--Read messages one ofter another
c
 200  DO 300 i = 1, nmes
c
c--Identify each message
c
c  0) wait before processing messages
c
         IF (string(i)(1:2).EQ.'at') THEN
            READ (string(i), *) dummy3, attime
            IF (gt.LT.attime) THEN
               iwait = 1
               WRITE (*, 99001) i
99001          FORMAT (' message ', I2, ' read. waiting to process.')
               GOTO 400
            ELSE
               GOTO 300
            ENDIF
         ENDIF
c
c  4) change minimum time for keeping GRAPE
c
         IF (string(i)(1:5).EQ.'tkeep') THEN
            READ (string(i), *) dummy5, tkeep
            IF (tkeep.LT.3.0) tkeep = 3.0
            WRITE (*, 99010) i, tkeep
99010    FORMAT (' message ', I2, ' read. tkeep changed to : ',1PE12.5)
            GOTO 300
         ENDIF
c
c--Unexpected messages
c
         WRITE (*, 99090) string(i)
99090    FORMAT (' unexpected message received : ', /, A20)

 300  CONTINUE
c
c--Erase message file
c
 400  CONTINUE
c      IF (iwait.EQ.0) THEN
c         CLOSE (idisk2, STATUS='delete')
c      ELSE
         CLOSE (idisk2)
c      ENDIF

      RETURN
      END
