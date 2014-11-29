      PROGRAM ptdump
      IMPLICIT NONE
c This programme will produce an ASCII dump from an sphNG P-file
c It should be safe for runs with multiple sinks, but is untested 
c nptmasstot > 1.
c Build with something like
c     ifort -o ptdump ptdump.f
      
      INTEGER idim1, idim2
      PARAMETER (idim1 = 20000)
      PARAMETER (idim2 = 2000000)

      INTEGER j, k, p, q
      
      INTEGER nptmasstot
      REAL*8  rt_tcomp, realtime
      DIMENSION nptmasstot(idim1), rt_tcomp(idim1), realtime(idim1)
 
      INTEGER*8 iunique
      INTEGER*2 nactotal
      REAL*4 rho
      REAL*8 xyzm, vxyz
      DIMENSION iunique(idim2), nactotal(idim2), rho(idim2)
      DIMENSION xyzm(4,idim2), vxyz(3,idim2)

      INTEGER iunitin, iunitout
      PARAMETER (iunitin = 20)
      PARAMETER (iunitout = 30)

      CHARACTER*21 infile, outfile


      WRITE (*,*) 'Enter input filename:'
      READ (*,99001) infile
      WRITE (*,*) 'Enter output filename:'
      READ (*,99001) outfile
99001 FORMAT(A21)

      OPEN (iunitin, FILE=infile, STATUS='OLD', FORM='UNFORMATTED',
     >               CONVERT="BIG_ENDIAN")
      OPEN (iunitout, FILE=outfile, STATUS='NEW', FORM='FORMATTED')

      j = 1
      k = 1
 100  READ(iunitin, END=200, ERR=199)
     > rt_tcomp(j), realtime(j), nptmasstot(j)
      
      IF (nptmasstot(j).EQ.0) GOTO 100
      DO 130 p = 1, nptmasstot(j)
        READ(iunitin, END=200, ERR=199)
     >  iunique(k),xyzm(1,k),xyzm(2,k),xyzm(3,k),
     >  vxyz(1,k),vxyz(2,k),vxyz(3,k),xyzm(4,k),
     >  rho(k),nactotal(k)
        k = k + 1
 130  CONTINUE
      IF (j.GT.idim1) THEN
        WRITE(*,99002)
99002   FORMAT('Error: idim1 exceeded')
        GOTO 200
      END IF
      IF (k.GT.idim2) THEN
        WRITE(*,99005)
99005   FORMAT('Error: idim2 exceeded')
        GOTO 200
      END IF
      j = j + 1
      GOTO 100
 199  WRITE(*,99003)
99003 FORMAT('Error: file read error..?')
 200  WRITE(*,99004)
99004 FORMAT('Read completed. Beginning write...')
      
      WRITE(*,99006)j,k
99006 FORMAT(2I6)

      p = 1
      DO 500 q = 1, j
        DO 520 p = p, p + nptmasstot(j)
          WRITE(iunitout,99007)
     >    realtime(q),nptmasstot(q),iunique(p),nactotal(p),
     >    xyzm(1,p),xyzm(2,p),xyzm(3,p),xyzm(4,p),
     >    vxyz(1,p),vxyz(2,p),vxyz(3,p),rho(p)
 520    CONTINUE  
 500  CONTINUE
99007 FORMAT(1E15.7,3I8,8E15.7)

      CLOSE(iunitin)
      CLOSE(iunitout)
      STOP
      ENDPROGRAM
