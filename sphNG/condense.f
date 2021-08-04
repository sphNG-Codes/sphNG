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
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bonnortbl'
      INCLUDE 'COMMONS/phase'

      DIMENSION gauss(1001)

      CHARACTER*1 ians
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine condense ')

      third = 1./3.
      rmax2 = rmax*rmax

      WRITE (*,*) 'Do you want exp[-3 x^2], r^-1, or r^-2, or ',
     &     'Bonnor-Ebert sphere, or polytrope (0, 1, 2, 3, 4)?'
      READ (iread,*) islope

      ians = 'n'
88001 FORMAT (A1)

      IF (islope.EQ.0) THEN
         OPEN (itablerd, FILE='/home/mbate/Important/Gaussian_Data',
     &         STATUS='old', ACTION='read')
         DO i = 1, 1001
            READ (itablerd,*,END=50) dummy, dummy, gauss(i)
         ENDDO
         CLOSE (itablerd)
         GOTO 100
 50      WRITE (*,*) 'ERROR - Gaussian_Data'
         CALL quit(0)
      ELSEIF (islope.EQ.1 .OR. islope.EQ.2) THEN
         WRITE (*,*) 'Change particle masses to be uniform density?'
         READ (iread,88001) ians
      ELSEIF (islope.EQ.3) THEN
         WRITE (*,*) 'Enter concentration parameter, xi'
         READ (iread,*) xi
         CALL bonnorebert(xi)
      ELSEIF (islope.EQ.4) THEN
         WRITE (*,*) 'Enter polytropic index (3 for solar-type star)'
         READ (iread,*) poly_index
         WRITE (*,*) 'Enter concentration parameter, xi'
         READ (iread,*) xi
         CALL polytrope(poly_index,xi)
      ELSE
         WRITE (*,*) 'ERROR - Invalid choice'
         CALL quit(0)
      ENDIF

 100  fractot = 0.
      rnew = 1.
      DO i = nptmass + 1, npart
         IF (MOD(i,100000).EQ.0) print *,'Done ',i,' particles'
         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         r2 = xi*xi + yi*yi + zi*zi
         r1 = SQRT(r2)
         r1rm = r1/rmax
c
c--Gaussian density profile
c
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
c
c--Power-law density profiles
c
         IF (islope.EQ.1) rnew = rmax*(r1rm**(1.50))
         IF (islope.EQ.2) rnew = rmax*(r1rm**(3.00))
c
c--Bonnor-Ebert or polytrope density profile
c
         IF (islope.EQ.3 .OR. islope.EQ.4) THEN
            xmassfraccontained = r1rm**3
            ipos1 = 1
            ipos2 = ibelast
 240        iposnew = (ipos1+ipos2)/2

            IF (xmassfraccontained.LT.bonnor_radmass(2,iposnew)) THEN
               ipos2=iposnew
            ELSE
               ipos1=iposnew
            ENDIF
            IF (ipos2-ipos1.LE.1) THEN
               ipos = ipos1
               GOTO 250
            ENDIF
            GOTO 240

 250        IF (ipos.LT.2) THEN
               ipos = 2
            ELSEIF (ipos.GT.ibelast) THEN
               ipos = ibelast-1
            ENDIF
            value1 = bonnor_radmass(2,ipos)
            value2 = bonnor_radmass(2,ipos-1)
            diff = (value1 - xmassfraccontained)/(value1 - value2)
            rnew = rmax*(bonnor_radmass(1,ipos) - 
     &           diff*(bonnor_radmass(1,ipos) - 
     &           bonnor_radmass(1,ipos-1)))
         ENDIF
c
c--Set actual particle quantities
c
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
c
c--Leave hole?
c
      WRITE (*,*) 'Do you want a central hole (y/n)?'
      READ (iread,88001) ians
      IF (ians.EQ.'y' .OR. ians.EQ.'Y') THEN
         WRITE (*,*) 'Enter radius of hole (code units)'
         READ (iread,*) radhole
         
         DO i = nptmass + 1, npart
            xi = xyzmh(1,i)
            yi = xyzmh(2,i)
            zi = xyzmh(3,i)
            r2 = xi*xi + yi*yi + zi*zi
            IF (r2.LT.radhole**2) iphase(i) = -1
         END DO
c
c--Completely remove particles within the hole
c
         fractot = 0.
         inew = nptmass + 1
         DO i = nptmass + 1, npart
            IF (iphase(i).EQ.0) THEN
               IF (i.NE.inew) THEN
                  xyzmh(:,inew) = xyzmh(:,i)
                  disfrac(inew) = disfrac(i)
               ENDIF
               fractot = fractot + disfrac(inew)
               inew = inew + 1
            ENDIF
         END DO
         npart = inew - 1
      ENDIF

      fractot = fractot / FLOAT(npart - nptmass)
      WRITE (*,*) 'fractot=',fractot
      DO i = nptmass + 1, npart
         disfrac(i) = disfrac(i)/fractot
      END DO

      RETURN
      END
