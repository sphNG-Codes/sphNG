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
      INCLUDE 'COMMONS/pxpy'

      npart = np + nptmass
      IF (igeom.EQ.9) THEN
         DO i = nptmass + 1, npart
 100        CONTINUE
            radius=(((1.0+variation)**1.5-(1.0-variation)**1.5)*
     &        ran1(1) + (1.0-variation)**1.5)**(2.0/3.0)
            phi = 2.0*(ran1(1)-0.5)*phibound
            xyzmh(1,i) = radius*COS(phi)
            xyzmh(2,i) = radius*SIN(phi)

            xyzmh(3,i) = radius*hoverr*gasdev(1)
c
c--No radial height profile
c
c            xyzmh(3,i) = hoverr*gasdev(1)
c
c--Uniform density
c
            IF (.FALSE.) THEN
               xyzmh(1,i) = (1.0-variation)+2.5*variation*ran1(1)
               xyzmh(2,i) = (1.0+variation)*
     &              (-phibound+2.0*phibound*ran1(1))
               phi = ATAN2(xyzmh(2,i),xyzmh(1,i))

               IF (phi.GT.phibound .OR. phi.LT.-phibound) GOTO 100
c
c--Uniform in z
c
               xyzmh(3,i) = -hoverr+2.0*hoverr*ran1(1)
            ENDIF

            IF (iexf.EQ.7 .AND. (ibound.EQ.100 .OR. ibound.EQ.104)) THEN
               IF (ibound.EQ.104) THEN
                  rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               ELSE
                  rtemp = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2 +
     &                 xyzmh(3,i)**2)
               ENDIF
               IF ((rtemp.GT.(1.0+variation)) .OR. 
     &              ((rtemp.LT.(1.0-variation)))) GOTO 100
               IF (ibound.EQ.100) THEN
                  IF (sqrt((xyzmh(1,i)-1.0)**2 + xyzmh(2,i)**2
     &                 + xyzmh(3,i)**2) .LE. (rplanet*pradfac(1)*2.1))
     &                 GOTO 100
               ENDIF
            ENDIF

            xyzmh(5,i) = (pi*2.0*hoverr*((1.0+variation)**2-
     &           (1.0-variation)**2)/(pi/phibound)/
     &           ((npart-nptmass)/50.0))**(1.0/3.0)
            xyzmh(5,i) = xyzmh(5,i)/2.
         END DO


      ELSE IF (igeom.EQ.10) THEN
         DO i = nptmass + 1, npart
 101        IF (abs(sdprof+2.0).LT.tiny) THEN
               radius = exp((log(rcyl)-log(rmind))*ran1(1) + log(rmind))
            ELSE
               radius=(((rcyl)**(2.+sdprof)-(rmind)**(2.+sdprof))*
     &              ran1(1) + (rmind)**(2.+sdprof))**(1./(2.+sdprof))
            ENDIF
            phi = 2.0*(ran1(1)-0.5)*pi
            xyzmh(1,i) = radius*COS(phi)
            xyzmh(2,i) = radius*SIN(phi)
            IF (use_tprof) THEN
               xyzmh(3,i) = radius**(0.5*(tprof+1)+1)*hoverr*gasdev(1)
            ELSE
               xyzmh(3,i) = radius*hoverr*gasdev(1)
            ENDIF

            IF (nptmass.GE.1) THEN
               DO j = 1, nptmass
                  IF (sqrt((xyzmh(1,i)-xyzmh(1,listpm(j)))**2 +
     &                 (xyzmh(2,i)-xyzmh(2,listpm(j)))**2
     &                 + (xyzmh(3,i)-xyzmh(3,listpm(j)))**2) .LE.
     &                 (xyzmh(5,listpm(j))*pradfac(j)*2.1))
     &                 GOTO 101
               ENDDO

            ELSEIF ((ibound.EQ.102 .OR. ibound.EQ.103)
     &              .AND. irotpot.EQ.1) THEN
               rtemp = sqrt((xyzmh(1,i)-rorbit_orig)**2 +
     &              xyzmh(2,i)**2 + xyzmh(3,i)**2)
               IF (rtemp.LE.(rplanet*pradfac(1)*2.1)) GOTO 101
            ENDIF
            xyzmh(5,i) = 0.1
         END DO
      ELSE
         
         WRITE (*,99001)
99001        FORMAT ('NOT IMPLEMENTED')
         CALL quit(0)
      ENDIF

      DO i = nptmass + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END
