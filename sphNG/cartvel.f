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
      INCLUDE 'COMMONS/nonidealsetup'

      CHARACTER*1 iok, igr, cvf
      INTEGER ivf(10)
      REAL kpiLx, kpiLy, kpiLz
      REAL vampx(10),vampy(10),vampz(10),kx(10),ky(10),kz(10)
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
      READ (iread,*) ireg
      IF (ireg.GT.10) GOTO 40

      icount = 1
      DO 120 i = 1, ireg
c
c--Ask for global region
c
         igr = 'n'
 60      IF (ireg.EQ.1) THEN
            WRITE (*,*) '   Is region global (y/n)?'
            READ (iread,*) igr
         ENDIF
         IF (igr.EQ.'y') THEN
            rxmin(icount) = xmin
            rxmax(icount) = xmax
            rymin(icount) = ymin
            rymax(icount) = ymax
            rzmin(icount) = zmin
            rzmax(icount) = zmax
         ELSE
            WRITE (*,*) '   Enter xmin of region ',icount
            READ (iread,*) rxmin(icount)
            WRITE (*,*) '   Enter xmax of region ',icount
            READ (iread,*) rxmax(icount)
            WRITE (*,*) '   Enter ymin of region ',icount
            READ (iread,*) rymin(icount)
            WRITE (*,*) '   Enter ymax of region ',icount
            READ (iread,*) rymax(icount)
            WRITE (*,*) '   Enter zmin of region ',icount
            READ (iread,*) rzmin(icount)
            WRITE (*,*) '   Enter zmax of region ',icount
            READ (iread,*) rzmax(icount)
         ENDIF
c
c--Is velocity constant or sin wave?
c
         cvf = "a"
         DO WHILE (cvf.NE.'c' .and. cvf.NE.'w')
            WRITE (*,*) '   Constant velocity in region', icount, '(c)'
            WRITE (*,*) '   Sin-wave velocity in region', icount, '(w)'
            READ (iread,*) cvf
         END DO
         IF (cvf.eq.'c') THEN
c
c--Constant velocity
c
            ivf(icount) = 0

            WRITE (*,*) '   Enter x velocity of region ',icount,
     &           ' in units of av. sound speed'
            READ (iread,*) xvel(icount)
            WRITE (*,*) '   Enter y velocity of region ',icount,
     &           ' in units of av. sound speed'
            READ (iread,*) yvel(icount)
            WRITE (*,*) '   Enter z velocity of region ',icount,
     &           ' in units of av. sound speed'
            READ (iread,*) zvel(icount)
         ELSE
c
c--Sine wave velocity
c
           ivf(icount) = 1

           WRITE (*,*) '   Enter x-velocity amplitude of region ',icount
     &           ,' in units of av. sound speed'
           READ (iread,*) vampx(icount)
           WRITE (*,*) '   Enter y-velocity amplitude of region ',icount
     &           ,' in units of av. sound speed'
           READ (iread,*) vampy(icount)
           WRITE (*,*) '   Enter z-velocity amplitude of region ',icount
     &           ,' in units of av. sound speed'
           READ (iread,*) vampz(icount)
           WRITE (*,*) '   Enter x wavenumber of region ',icount
     &           ,' in units 2pi/L'
           READ (iread,*) kx(icount)
           WRITE (*,*) '   Enter y wavenumber of region ',icount
     &           ,' in units 2pi/L'
           READ (iread,*) ky(icount)
           WRITE (*,*) '   Enter z wavenumber of region ',icount
     &           ,' in units 2pi/L'
           READ (iread,*) kz(icount)
         ENDIF

 100     WRITE (*,*) ' Is region ',icount,' correct (y/n)? '
         READ (iread, 99004) iok
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
      ELSE IF (encal.EQ.'a' .OR. encal.EQ.'c') THEN
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
         CALL quit(0)
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
               IF (ivf(j).EQ.0) THEN
c
c--Constant speed
c
                  IF (vxyzu(1,i).EQ.1.0E10) THEN
                     vxyzu(1,i) = xvel(j)*vsound
                     vxyzu(2,i) = yvel(j)*vsound
                     vxyzu(3,i) = zvel(j)*vsound
                  ELSE
                     vxyzu(1,i) = (vxyzu(1,i) + xvel(j)*vsound)/2.0
                     vxyzu(2,i) = (vxyzu(2,i) + yvel(j)*vsound)/2.0
                     vxyzu(3,i) = (vxyzu(3,i) + zvel(j)*vsound)/2.0
                  ENDIF
               ELSE
c
c--Sine wave speed                                                
c
                  kpiLx = 2.0*pi/ABS(rxmax(j)-rxmin(j))*kx(j)
                  kpiLy = 2.0*pi/ABS(rymax(j)-rymin(j))*ky(j)
                  kpiLz = 2.0*pi/ABS(rzmax(j)-rzmin(j))*kz(j)
                  
                  print *,'kpiLz',i,kpiLz,pi,kz(j),kx(j),valfven

                  IF (vxyzu(1,i).EQ.1.0E10) THEN
                     vxyzu(1,i) = vampx(j)*valfven*SIN( kpiLx*xyzmh(1,i)
     1                    + kpiLy*xyzmh(2,i)
     2                    + kpiLz*xyzmh(3,i) )
                     vxyzu(2,i) = vampy(j)*valfven*SIN( kpiLx*xyzmh(1,i)
     1                    + kpiLy*xyzmh(2,i)
     2                    + kpiLz*xyzmh(3,i) )
                     vxyzu(3,i) = vampz(j)*valfven*SIN( kpiLx*xyzmh(1,i)
     1                    + kpiLy*xyzmh(2,i)
     2                    + kpiLz*xyzmh(3,i) )
                  ELSE
                     print*, 'This is called.  Find out what this is'
                     vxyzu(1,i) = 0.5*(vxyzu(1,i)
     1                    + vampx(j)*valfven*SIN(kpiLx*xyzmh(1,i)
     2                    + kpiLy*xyzmh(2,i)
     3                    + kpiLz*xyzmh(3,i) ))
                     vxyzu(2,i) = 0.5*(vxyzu(2,i)
     1                    + vampy(j)*valfven*SIN(kpiLx*xyzmh(1,i)
     2                    + kpiLy*xyzmh(2,i)
     3                    + kpiLz*xyzmh(3,i) ))
                     vxyzu(3,i) = 0.5*(vxyzu(3,i)
     1                    + vampz(j)*valfven*sin(kpiLx*xyzmh(1,i)
     2                    + kpiLy*xyzmh(2,i)
     3                    + kpiLz*xyzmh(3,i) ))
                  ENDIF
               ENDIF
            ENDIF
 220     CONTINUE
 240  CONTINUE 

      RETURN
      END
