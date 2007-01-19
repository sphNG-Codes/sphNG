      SUBROUTINE boundry(npart,list,nlst0,xyzmh,vxyzu,idonebound)
c************************************************************
c                                                           *
c   This subroutine checks to see if particles will surpass *
c   boundary and resets them to the boundary if they do     *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim), vxyzu(4,idim)
      DIMENSION list(idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/mhd'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine boundary')

      iinner = 0
      iouter = 0
      idonebound = 0
      IF ( ibound.EQ.1) THEN
         DO 250 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 250
            delta = (1.3-ran1(1))*0.5*xyzmh(5,i)
            ichan = 0
            IF (xyzmh(1,i).LE.xmin) THEN
               ichan = ichan + 1
               vxyzu(1,i) = 0.
               xyzmh(1,i) = xmin + delta
            ENDIF
            IF (xyzmh(1,i).GE.xmax) THEN
               ichan = ichan + 1
               vxyzu(1,i) = 0.
               xyzmh(1,i) = xmax - delta
            ENDIF
            IF (xyzmh(2,i).LE.ymin) THEN
               ichan = ichan + 1
               vxyzu(2,i) = 0.
               xyzmh(2,i) = ymin + delta
            ENDIF
            IF (xyzmh(2,i).GE.ymax) THEN
               ichan = ichan + 1
               vxyzu(2,i) = 0.
               xyzmh(2,i) = ymax - delta
            ENDIF
            IF (xyzmh(3,i).LE.zmin) THEN
               ichan = ichan + 1
               vxyzu(3,i) = 0.
               xyzmh(3,i) = zmin + delta
            ENDIF
            IF (xyzmh(3,i).GE.zmax) THEN
               ichan = ichan + 1
               vxyzu(3,i) = 0.
               xyzmh(3,i) = zmax - delta
            ENDIF
            IF (ichan.NE.0) THEN
               iouter = iouter + 1
               idonebound = 1
            ENDIF
 250     CONTINUE

      ELSE IF ( ibound.EQ.2 ) THEN
         rcyl2 = rcyl * rcyl
         rmind2 = rmind * rmind
         DO 260 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 260
            delta = (1.3-ran1(1))*0.5*xyzmh(5,i)
            rcyld = rcyl - delta
            rmindd = rmind + delta
            zmind = zmin + delta
            zmaxd = zmax - delta
            ichan = 0
            r2 = xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i) 
            IF ( r2.GT.rcyl2) THEN
               ichan = ichan + 1
               r = SQRT(r2)
               zi = xyzmh(3,i)
               yi = xyzmh(2,i)*rcyld/r
               xi = xyzmh(1,i)*rcyld/r
c
c--Remove velocities perpendicular to boundary
c
               vr = (vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i))/r
               IF (vr.GT.0) THEN
                  vxyzu(1,i) = vxyzu(1,i) - vr*xyzmh(1,i)/r
                  vxyzu(2,i) = vxyzu(2,i) - vr*xyzmh(2,i)/r
               ENDIF
c
c--Adjust to conserve angular momentum
c
               IF (ifcor.EQ.1) THEN
                  vxyzu(1,i) = vxyzu(1,i) * xyzmh(2,i) / yi
                  vxyzu(2,i) = vxyzu(2,i) * xyzmh(1,i) / xi
               ELSEIF (ifcor.EQ.2) THEN
                  vxyzu(2,i) = vxyzu(2,i) * xyzmh(3,i) / zi
                  vxyzu(3,i) = vxyzu(3,i) * xyzmh(2,i) / yi        
               ENDIF
               xyzmh(2,i) = yi
               xyzmh(1,i) = xi
            ENDIF

            IF ( r2.LT.rmind2) THEN
               iinner = iinner + 1
               r = sqrt(r2)
               zi = xyzmh(3,i)
               yi = xyzmh(2,i)*rmindd/r
               xi = xyzmh(1,i)*rmindd/r
c
c--Remove velocities perpendicular to boundary
c
               vr = (vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i))/r
               IF (vr.LT.0) THEN
                  vxyzu(1,i) = vxyzu(1,i) - vr*xyzmh(1,i)/r
                  vxyzu(2,i) = vxyzu(2,i) - vr*xyzmh(2,i)/r
               ENDIF
c
c--Adjust to conserve angular momentum
c
               IF (ifcor.EQ.1) THEN
                  vxyzu(1,i) = vxyzu(1,i) * xyzmh(2,i) / yi
                  vxyzu(2,i) = vxyzu(2,i) * xyzmh(1,i) / xi
               ELSEIF (ifcor.EQ.2) THEN
                  vxyzu(2,i) = vxyzu(2,i) * xyzmh(3,i) / zi
                  vxyzu(3,i) = vxyzu(3,i) * xyzmh(2,i) / yi
               ENDIF
               xyzmh(2,i) = yi
               xyzmh(1,i) = xi
            ENDIF

            IF (xyzmh(3,i).LE.zmin) THEN
               ichan = ichan + 1
               vxyzu(3,i) = 0.
               xyzmh(3,i) = zmind
            ENDIF
            IF (xyzmh(3,i).GE.zmax) THEN
               ichan = ichan + 1
               vxyzu(3,i) = 0.
               xyzmh(3,i) = zmaxd
            ENDIF
            IF (ichan.NE.0) THEN
               iouter = iouter + 1
               idonebound = 1
            ENDIF
 260     CONTINUE

      ELSE IF ( ibound.EQ.3 .OR. ibound.EQ.6) THEN
         rmax2 = rmax * rmax
         DO 300 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 300
            r2 = (xyzmh(1,i)*xyzmh(1,i)+xyzmh(2,i)*xyzmh(2,i)+
     &           xyzmh(3,i)*xyzmh(3,i)) 
            IF ( r2.GT.rmax2 ) THEN
               iouter = iouter + 1
               idonebound = 1
               r=sqrt(r2)
c
c--Put particle inside boundary
c
               delta = (1.3-ran1(1))*0.5*xyzmh(5,i)
               rmaxd = rmax - delta
               zi = xyzmh(3,i)*rmaxd/r
               yi = xyzmh(2,i)*rmaxd/r
               xi = xyzmh(1,i)*rmaxd/r
c
c--Remove velocities perpendicular to boundary
c
               vr = ( vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i)
     &                + vxyzu(3,i)*xyzmh(3,i) ) / r
               IF ( vr.GT.0 ) THEN
                   vxyzu(1,i) = vxyzu(1,i) - vr * xyzmh(1,i) / r
                   vxyzu(2,i) = vxyzu(2,i) - vr * xyzmh(2,i) / r
                   vxyzu(3,i) = vxyzu(3,i) - vr * xyzmh(3,i) / r
               ENDIF
c
c--Adjust to conserve angular momentum
c
               IF (ifcor.EQ.1) THEN
                  vxyzu(1,i) = vxyzu(1,i) * xyzmh(2,i) / yi
                  vxyzu(2,i) = vxyzu(2,i) * xyzmh(1,i) / xi
               ELSEIF (ifcor.EQ.2) THEN
                  vxyzu(2,i) = vxyzu(2,i) * xyzmh(3,i) / zi
                  vxyzu(3,i) = vxyzu(3,i) * xyzmh(2,i) / yi
               ENDIF
               xyzmh(3,i) = zi
               xyzmh(2,i) = yi
               xyzmh(1,i) = xi
            ENDIF
 300     CONTINUE
      ELSE IF ( ibound.EQ.11) THEN
         DO 350 j = 1, nlst0
            i = list(j)
            IF (iphase(i).LT.0) GOTO 350
            ichan = 0
            xi = xyzmh(1,i)
            yi = xyzmh(2,i)
            zi = xyzmh(3,i)
            IF (xyzmh(1,i).LT.xmin) THEN
               ichan = ichan + 1
               xyzmh(1,i) = xyzmh(1,i) + (xmax - xmin)
            ENDIF
            IF (xyzmh(1,i).GE.xmax) THEN
               ichan = ichan + 1
               xyzmh(1,i) = xyzmh(1,i) - (xmax - xmin)
            ENDIF
            IF (xyzmh(2,i).LT.ymin) THEN
               ichan = ichan + 1
               xyzmh(2,i) = xyzmh(2,i) + (ymax - ymin)
            ENDIF
            IF (xyzmh(2,i).GE.ymax) THEN
               ichan = ichan + 1
               xyzmh(2,i) = xyzmh(2,i) - (ymax - ymin)
            ENDIF
            IF (xyzmh(3,i).LT.zmin) THEN
               ichan = ichan + 1
               xyzmh(3,i) = xyzmh(3,i) + (zmax - zmin)
            ENDIF
            IF (xyzmh(3,i).GE.zmax) THEN
               ichan = ichan + 1
               xyzmh(3,i) = xyzmh(3,i) - (zmax - zmin)
            ENDIF
            IF (imhd.EQ.idim) CALL copyBevol(i,i,Bevolxyz,
     &         xyzmh(1,i)-xi,xyzmh(2,i)-yi,xyzmh(3,i)-zi)

            IF (ichan.NE.0) THEN
               iouter = iouter + 1
               idonebound = 1
            ENDIF
 350     CONTINUE
      ENDIF

      IF (ibound.EQ.2 .AND. iinner.NE.0) WRITE (iprint,99009) iinner
      IF (iouter.NE.0) WRITE (iprint,99010) iouter

99009 FORMAT(' number of corrections inner boundary:',I6)
99010 FORMAT(' number of corrections outer boundary:',I6)

      RETURN
      END

c************************************************************
c                                                           *
c   This subroutine is for periodic boundaries without      *
c   ghost particles                                         *
c                                                           *
c   Given the separation of a particle pair, returns the    *
c   less of the separation and the separation across the    *
c   boundary.                                               *
c                                                           *
c   DJP 18.01.07                                            *
c                                                           *
c************************************************************
      SUBROUTINE modbound(dx,dy,dz)
      INCLUDE 'idim'
      INCLUDE 'COMMONS/rbnd'
      
c      print*,'xmin,max = ',xmin,xmax
c      print*,dx,dx - (xmax - xmin)*dx/abs(dx),
c     &       dy,dy - (ymax - ymin)*dy/abs(dy),
c     &       dz,dz - (zmax - zmin)*dz/abs(dz)
      IF (abs(dx).GT.tiny) THEN
         term = dx - (xmax - xmin)*dx/abs(dx)
         IF (abs(term).lt.abs(dx)) dx = term
      ENDIF
      IF (abs(dy).GT.tiny) THEN
         term = dy - (ymax - ymin)*dy/abs(dy)
         IF (abs(term).lt.abs(dy)) dy = term
      ENDIF
      IF (abs(dz).GT.tiny) THEN
         term = dz - (zmax - zmin)*dz/abs(dz)
         IF (abs(term).lt.abs(dz)) dz = term
      ENDIF
c      print*,'using ',dx,dy,dz,abs(term).lt.abs(dx)
      RETURN
      END SUBROUTINE modbound
