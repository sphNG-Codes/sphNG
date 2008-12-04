      SUBROUTINE cyldis(igeom,np)
c************************************************************
c                                                           *
c  This subroutine positions particles in a cylindrical     *
c     distribution                                          *
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
      INCLUDE 'COMMONS/typef'

      npart = np + nptmass
      IF (igeom.EQ.9) THEN
         DO i = nptmass + 1, npart
 100        radius=(((1.0+variation)**1.5-(1.0-variation)**1.5)*
     &        ran1(1) + (1.0-variation)**1.5)**(2.0/3.0)
            phi = 2.0*(ran1(1)-0.5)*phibound
            xyzmh(1,i) = radius*COS(phi)
            xyzmh(2,i) = radius*SIN(phi)
            xyzmh(3,i) = radius*hoverr*gasdev(1)

            IF (iexf.EQ.7 .AND. ibound.EQ.100) THEN
               rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2 +
     &              xyzmh(3,i)**2)
               IF ((rtemp.GT.(1.0+variation)) .OR. 
     &              ((rtemp.LT.(1.0-variation)))) GOTO 100
               IF (sqrt((xyzmh(1,i)-1.0)**2 + xyzmh(2,i)**2
     &              + xyzmh(3,i)**2) .LE. (rplanet*pradfac*2.1))
     &              GOTO 100
            ENDIF

            xyzmh(5,i) = (pi*2.0*hoverr*((1.0+variation)**2-
     &           (1.0-variation)**2)/(pi/phibound)/
     &           ((npart-nptmass)/50.0))**(1.0/3.0)
            xyzmh(5,i) = xyzmh(5,i)/2.
         END DO


      ELSE IF (igeom.EQ.10) THEN
         DO i = nptmass + 1, npart
 101        radius=(((rcyl)**1.5-(rmind)**1.5)*
     &        ran1(1) + (rmind)**1.5)**(2.0/3.0)
            phi = 2.0*(ran1(1)-0.5)*pi
            xyzmh(1,i) = radius*COS(phi)
            xyzmh(2,i) = radius*SIN(phi)
            xyzmh(3,i) = radius*hoverr*gasdev(1)

            rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2 +
     &           xyzmh(3,i)**2)
            IF (nptmass.GE.1) THEN
               DO j = 1, nptmass
                  IF (sqrt((xyzmh(1,i)-xyzmh(1,listpm(j)))**2 +
     &                 (xyzmh(2,i)-xyzmh(2,listpm(j)))**2
     &                 + (xyzmh(3,i)-xyzmh(3,listpm(j)))**2) .LE.
     &                 (rplanet*pradfac*2.1))
     &                 GOTO 101
               ENDDO
            ENDIF
            xyzmh(5,i) = 0.1
         END DO
      ELSE
         
         WRITE (*,99001)
99001        FORMAT ('NOT IMPLEMENTED')
         CALL quit
      ENDIF

      DO i = nptmass + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END
