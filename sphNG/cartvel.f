      SUBROUTINE cartvel
c************************************************************
c                                                           *
c     This subroutine gives particles velocities in a       *
c              cartesian coordinate distribution            *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/regionslocal'

      CHARACTER*1 iok
c
c--Allow for tracing flow
c
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine cartvel ')
c
c--Set Velocity Areas In Cartesian Coordinates
c
      rcyl2 = rcyl * rcyl 
      rmind2 = rmind * rmind
      rmax2 = rmax * rmax

 40   WRITE (*, 88002)
88002 FORMAT (' Enter number of different regions (max 10)')   
      READ (*,*) ireg
      IF (ireg.GT.10) GOTO 40

      icount = 1
      DO 120 i = 1, ireg
 60      WRITE (*,*) '   Enter xmin of region ',icount
         READ (*,*) rxmin(icount)
         WRITE (*,*) '   Enter xmax of region ',icount
         READ (*,*) rxmax(icount)
         WRITE (*,*) '   Enter ymin of region ',icount
         READ (*,*) rymin(icount)
         WRITE (*,*) '   Enter ymax of region ',icount
         READ (*,*) rymax(icount)
         WRITE (*,*) '   Enter zmin of region ',icount
         READ (*,*) rzmin(icount)
         WRITE (*,*) '   Enter zmax of region ',icount
         READ (*,*) rzmax(icount)

         WRITE (*,*) '   Enter x velocity of region ',icount,
     &         ' in units of av. sound speed'
         READ (*,*) xvel(icount)
         WRITE (*,*) '   Enter y velocity of region ',icount,
     &         ' in units of av. sound speed'
         READ (*,*) yvel(icount)
         WRITE (*,*) '   Enter z velocity of region ',icount,
     &         ' in units of av. sound speed'
         READ (*,*) zvel(icount)

 100     WRITE (*,*) ' Is region ',icount,' correct (y/n)? '
         READ (*, 99004) iok
         IF ((iok.NE.'y') .AND. (iok.NE.'n')) GOTO 100
         IF (iok.NE.'y') GOTO 60
         icount = icount + 1

 120  CONTINUE

      IF (encal.EQ.'i') THEN
         total = 0.
         DO i = 1, npart
            total = total + sqrt((2./3.)*vxyzu(4,i))
         END DO
         vsound = total/FLOAT(npart)
      ELSE IF (encal.EQ.'a') THEN
         total = 0.
         DO i = 1, npart
            total = total + sqrt(gamma*(gamma-1.0)*vxyzu(4,i))
         END DO
         vsound = total/FLOAT(npart)
      ELSE IF (encal.EQ.'p') THEN
         total = 0.
         gama1 = gamma - 1.
         DO i = 1, npart
            total = total + sqrt(gamma*(2./3.)*RK2*rhozero**gama1)
         END DO
         vsound = total/FLOAT(npart)
      ELSE
         WRITE (*,*) ' NOT IMPLEMENTED'
         CALL quit
      END IF

      WRITE (*,*) 'Average sound speed is ',vsound

      DO i = 1, npart
         vxyzu(1,i) = 1.0E10
      END DO

      DO 240 i = 1, npart
         DO 220 j = 1, ireg
            IF ((xyzmh(1,i).GE.rxmin(j)) .AND. 
     &           (xyzmh(1,i).LE.rxmax(j)) .AND.
     &           (xyzmh(2,i).GE.rymin(j)) .AND. 
     &           (xyzmh(2,i).LE.rymax(j)) .AND.
     &           (xyzmh(3,i).GE.rzmin(j)) .AND. 
     &           (xyzmh(3,i).LE.rzmax(j))) THEN
               IF (vxyzu(1,i).EQ.1.0E10) THEN
                  vxyzu(1,i) = xvel(j)*vsound
                  vxyzu(2,i) = yvel(j)*vsound
                  vxyzu(3,i) = zvel(j)*vsound
               ELSE
                  vxyzu(1,i) = (vxyzu(1,i) + xvel(j)*vsound)/2.0
                  vxyzu(2,i) = (vxyzu(2,i) + yvel(j)*vsound)/2.0
                  vxyzu(3,i) = (vxyzu(3,i) + zvel(j)*vsound)/2.0
               ENDIF
            ENDIF
 220     CONTINUE
 240  CONTINUE 

      RETURN
      END
