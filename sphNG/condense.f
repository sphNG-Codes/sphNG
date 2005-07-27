      SUBROUTINE condense
c************************************************************
c                                                           *
c  This subroutine re-positions uniform density particles   *
c     in a centrally condensed distribution                 *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/flag'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'

      DIMENSION gauss(1001)

      CHARACTER*1 ians
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine condense ')

      third = 1./3.
      rmax2 = rmax*rmax

      WRITE (*,*) 'Do you want exp[-3 x^2], r^-1 or r^-2 (0, 1 or 2)?'
      READ (*,*) islope

      IF (islope.EQ.0) THEN
         OPEN (19, FILE='/home/mbate/Important/Gaussian_Data')
         DO i = 1, 1001
            READ (19,*,END=50) dummy, dummy, gauss(i)
         ENDDO
         CLOSE (19)
         GOTO 100
 50      WRITE (*,*) 'ERROR - Gaussian_Data'
         CALL quit
 100  ENDIF

      ians = 'n'
      IF (islope.EQ.1 .OR. islope.EQ.2) THEN
         WRITE (*,*) 'Change particle masses to be uniform density?'
         READ (*,88001) ians
      ENDIF
88001 FORMAT (A1)

      fractot = 0.
      DO i = nptmass + 1, npart
         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         r2 = xi*xi + yi*yi + zi*zi
         r1 = SQRT(r2)
         r1rm = r1/rmax

         IF (islope.EQ.0) THEN
            xmassfraccontained = r1rm**3
            DO ipos = 1, 1001
               IF (xmassfraccontained.LT.gauss(ipos)) GOTO 150
            END DO
 150        value1 = gauss(ipos)
            value2 = gauss(ipos-1)
            diff = (value1 - xmassfraccontained)/(value1 - value2)
            rnew = rmax*(ipos - 1 - diff)/1000.0
         ENDIF

         IF (islope.EQ.1) rnew = rmax*(r1rm**(1.50))
         IF (islope.EQ.2) rnew = rmax*(r1rm**(3.00))

         xyzmh(1,i) = rnew*xi/r1
         xyzmh(2,i) = rnew*yi/r1
         xyzmh(3,i) = rnew*zi/r1

         IF (rnew/r1.LT.0.001) THEN
            WRITE(*,*)'ERROR in condense - particle too near the origin'
         ENDIF
         IF (rnew/r1.LT.0.00001) THEN
            WRITE(*,*)'ERROR in condense - particle at r=0 '
         ENDIF
         IF (ians.EQ.'y') THEN
            IF (islope.EQ.1) disfrac(i) = (rnew/r1)**1.5
            IF (islope.EQ.2) disfrac(i) = (rnew/r1)**3.0
         ELSE
            disfrac(i) = 1.0
         ENDIF
         fractot = fractot + disfrac(i)
         xyzmh(5,i) = xyzmh(5,i)*(rnew/r1)
      END DO

      fractot = fractot / FLOAT(npart - nptmass)
      WRITE (*,*) 'fractot=',fractot
      DO i = nptmass + 1, npart
         disfrac(i) = disfrac(i)/fractot
      END DO

      RETURN
      END
