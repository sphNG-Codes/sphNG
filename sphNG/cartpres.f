      SUBROUTINE cartpres
c************************************************************
c                                                           *
c     This subroutine gives particles masses in a cartesian *
c              coordinate distribution                      *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/regionslocal'

      CHARACTER*1 iok
c
c--Allow for tracing flow
c
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine cartpres ')
c
c--Set Variable Pressure Areas In Cartesian Coordinates
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
 60      WRITE (*, 88004) icount
88004    FORMAT ('   Enter xmin of region ', I2)
         READ (*,*) rxmin(icount)
         WRITE (*, 88006) icount
88006    FORMAT ('   Enter xmax of region ', I2)
         READ (*,*) rxmax(icount)
         WRITE (*, 88008) icount
88008    FORMAT ('   Enter ymin of region ', I2)
         READ (*,*) rymin(icount)
         WRITE (*, 88010) icount
88010    FORMAT ('   Enter ymax of region ', I2)
         READ (*,*) rymax(icount)
         WRITE (*, 88012) icount
88012    FORMAT ('   Enter zmin of region ', I2)
         READ (*,*) rzmin(icount)
         WRITE (*, 88014) icount
88014    FORMAT ('   Enter zmax of region ', I2)
         READ (*,*) rzmax(icount)

 80      WRITE (*, 88016) icount
88016    FORMAT ('   Enter relative pressure of region ', I2, 
     &          ' (0.0 to 1.0)')
         READ (*,*) pres(icount)
         IF ((pres(icount).LT.0.).OR.(pres(icount).GT.1.)) GOTO 80

 100     WRITE (*, 88018) icount
88018    FORMAT (' Is region ', I2,' correct (y/n)? ')
         READ (*, 99004) iok
         IF ((iok.NE.'y').AND.(iok.NE.'n')) GOTO 100
         IF (iok.NE.'y') GOTO 60

         icount = icount + 1
 120  CONTINUE

      xmax5 = xmax * 0.5
      ymax5 = ymax * 0.5
      zmax5 = zmax * 0.5
      fractot = 0.

      DO 240 i = 1, npart
         DO 220 j = 1, ireg
            IF ((xyzmh(1,i).GE.rxmin(j)) .AND. 
     &           (xyzmh(1,i).LE.rxmax(j)) .AND.
     &           (xyzmh(2,i).GE.rymin(j)) .AND. 
     &           (xyzmh(2,i).LE.rymax(j)) .AND.
     &           (xyzmh(3,i).GE.rzmin(j)) .AND. 
     &           (xyzmh(3,i).LE.rzmax(j))) THEN
               disfrac(i) = pres(j)
               fractot = fractot + disfrac(i)
               GOTO 240
            ENDIF
 220     CONTINUE

 240  CONTINUE 

      fractot = fractot/FLOAT(npart)
      WRITE (*,*) fractot
      DO 300 i = 1, npart
         disfrac(i) = disfrac(i)/fractot
 300  CONTINUE

      RETURN
      END
