      SUBROUTINE cartmas
c************************************************************
c                                                           *
c     This subroutine gives particles masses in a Cartesian *
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
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/regionslocal'

      CHARACTER*1 iok, iform
c
c--Allow for tracing flow
c
      third = 1./3.
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine cartmas ')

      WRITE (*, 88011)
88011 FORMAT (' Do you want a density gradient, or stepped density',
     &     '(g/s)?')   
      READ (*,99004) iform

      fractot = 0.

      IF (iform.EQ.'g') THEN
         WRITE (*, 88032)
88032    FORMAT (' Enter fractional gradient over 1 unit length')   
         WRITE (*, 88033)
88033    FORMAT ('   In x-direction:')   
         READ (*,*) fracgradx
         WRITE (*, 88034)
88034    FORMAT ('   In y-direction:')   
         READ (*,*) fracgrady
         WRITE (*, 88036)
88036    FORMAT ('   In y-direction:')   
         READ (*,*) fracgradz

         DO i = nptmass + 1, npart
            disfrac(i) = (1.0 + xyzmh(1,i)*fracgradx)*
     &           (1.0 + xyzmh(2,i)*fracgrady)*(1.0 + xyzmh(3,i)*
     &           fracgradz)*disfrac(i)
            fractot = fractot + disfrac(i)
         END DO

      ELSE
c
c--Set Condensed Areas In Cartesian Coordinates
c
         rcyl2 = rcyl * rcyl 
         rmind2 = rmind * rmind
         rmax2 = rmax * rmax

 40      WRITE (*, 88002)
88002    FORMAT (' Enter number of different regions (max 10)')   
         READ (*,*) ireg
         IF (ireg.GT.10) GOTO 40

         icount = 1
         DO 120 i = 1, ireg
 60         WRITE (*, 88004) icount
88004       FORMAT ('   Enter xmin of region ', I5)
            READ (*,*) rxmin(icount)
            WRITE (*, 88006) icount
88006       FORMAT ('   Enter xmax of region ', I5)
            READ (*,*) rxmax(icount)
            WRITE (*, 88008) icount
88008       FORMAT ('   Enter ymin of region ', I5)
            READ (*,*) rymin(icount)
            WRITE (*, 88010) icount
88010       FORMAT ('   Enter ymax of region ', I5)
            READ (*,*) rymax(icount)
            WRITE (*, 88012) icount
88012       FORMAT ('   Enter zmin of region ', I5)
            READ (*,*) rzmin(icount)
            WRITE (*, 88014) icount
88014       FORMAT ('   Enter zmax of region ', I5)
            READ (*,*) rzmax(icount)

 80         WRITE (*, 88016) icount
88016       FORMAT ('   Enter relative density of region ', I5)
            READ (*,*) dens(icount)
            IF ((dens(icount).LT.0.)) GOTO 80

 100        WRITE (*, 88018) icount
88018       FORMAT (' Is region ', I2,' correct (y/n)? ')
            READ (*, 99004) iok
            IF ((iok.NE.'y').AND.(iok.NE.'n')) GOTO 100
            IF (iok.NE.'y') GOTO 60

            icount = icount + 1
 120     CONTINUE

         xmax5 = xmax * 0.5
         ymax5 = ymax * 0.5
         zmax5 = zmax * 0.5

         DO 240 i = nptmass + 1, npart
            DO 220 j = 1, ireg
               IF ((xyzmh(1,i).GE.rxmin(j)) .AND. 
     &              (xyzmh(1,i).LE.rxmax(j)) .AND.
     &              (xyzmh(2,i).GE.rymin(j)) .AND. 
     &              (xyzmh(2,i).LE.rymax(j)) .AND.
     &              (xyzmh(3,i).GE.rzmin(j)) .AND. 
     &              (xyzmh(3,i).LE.rzmax(j))) THEN
                  disfrac(i) = dens(j)*disfrac(i)
                  fractot = fractot + disfrac(i)
                  GOTO 240
               ENDIF
 220        CONTINUE
            WRITE (*,*) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)
 240     CONTINUE
      ENDIF

      fractot = fractot/FLOAT(npart-nptmass)
      WRITE (*,*) fractot
      DO 300 i = nptmass + 1, npart
         disfrac(i) = disfrac(i)/fractot
 300  CONTINUE

      RETURN
      END
