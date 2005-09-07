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
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'

      npart = np + nptmass
      IF (igeom.EQ.8) THEN
         DO i = nptmass + 1, npart
            radius=(((1.0+variation)**1.5-(1.0-variation)**1.5)*
     &        ran1(1) + (1.0-variation)**1.5)**(2.0/3.0)
            phi = 2.0*(ran1(1)-0.5)*phibound
            xyzmh(1,i) = radius*COS(phi)
            xyzmh(2,i) = radius*SIN(phi)
            xyzmh(3,i) = radius*hoverr*gasdev(1)
            xyzmh(5,i) = (pi*2.0*hoverr*((1.0+variation)**2-
     &           (1.0-variation)**2)/(pi/phibound)/
     &           ((npart-nptmass)/50.0))**(1.0/3.0)
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
