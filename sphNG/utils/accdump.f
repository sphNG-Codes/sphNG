      PROGRAM accdump
      IMPLICIT NONE
c This programme will produce an ASCII file from an sphNG A-file   
c Build with something like
c     ifort -o accdump accdump.f 

      INTEGER idim
      PARAMETER (idim = 2000000)

      INTEGER j, lastline

      INTEGER*8 index, index2, dummy, dummy2 
      REAL*8 realtime, xyzmh, vxyzu
      DIMENSION index(idim), index2(idim), realtime(idim)
      DIMENSION xyzmh(3,idim), vxyzu(3,idim)

      INTEGER iunitin, iunitout
      PARAMETER (iunitin = 20)
      PARAMETER (iunitout = 30)

      CHARACTER*21 infile, outfile


      WRITE (*,*) 'Enter input filename:'
      READ (*,99001) infile
      WRITE (*,*) 'Enter output filename:'
      READ (*,99001) outfile
99001 FORMAT(A21)
 
      realtime(:) = 99
      WRITE(*,*)realtime(2),realtime(50)

      OPEN (iunitin, FILE=infile, STATUS='OLD', FORM='UNFORMATTED',
     >               CONVERT="BIG_ENDIAN")
      OPEN (iunitout, FILE=outfile, STATUS='NEW', FORM='FORMATTED')

      j = 1
 100  READ(iunitin, END=200, ERR=199) 
     >     index(j), realtime(j),
     >     xyzmh(1,j), xyzmh(2,j), xyzmh(3,j),
     >     vxyzu(1,j), vxyzu(2,j), vxyzu(3,j),
     >     index2(j)
      
      lastline = j
      j = j + 1
      IF (j.GT.idim) THEN
        WRITE(*,99002)
99002   FORMAT('Error: idim exceeded')
        GOTO 200
      END IF
      GOTO 100
 199  WRITE(*,99003)
99003 FORMAT('Error: file read error..?')
 200  WRITE(*,99004)
99004 FORMAT('Read completed. Beginning write...')

      WRITE(*,*)realtime(2),realtime(50)

      CLOSE(iunitin)

      DO 300 j = 1, lastline
        WRITE(iunitout,99100)
     >  index(j), realtime(j),
     >  xyzmh(1,j), xyzmh(2,j), xyzmh(3,j),
     >  vxyzu(1,j), vxyzu(2,j), vxyzu(3,j),
     >  index2(j),dummy,dummy2
 300  CONTINUE
99100 FORMAT(' ',I7,7E15.7,I9,2I9)

      CLOSE(iunitout) 
      
      WRITE(*,99005)
99005 FORMAT('Done.')
      END PROGRAM
      
