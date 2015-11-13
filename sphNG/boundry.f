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
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rotat'
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
      ELSE IF ( ibound.EQ.12 ) THEN  
         rcyl2 = rcyl * rcyl
         DO 450 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 450
            rcyld = rcyl - delta
            ichan = 0
            r2 = xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i) 
            IF ( r2.GT.rcyl2) THEN
               iinner = iinner + 1
               r = sqrt(r2)
               delta = r - rcyl            
               ichan = ichan + 1
               zi = xyzmh(3,i)
               yi = xyzmh(2,i)*(rcyl-delta)/r
               xi = xyzmh(1,i)*(rcyl-delta)/r
c
c--Remove velocities perpendicular to boundary
c
               vr = (vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i))/r
               IF (vr.GT.0) THEN
                  vxyzu(1,i) = vxyzu(1,i) - 2.*vr*xyzmh(1,i)/r
                  vxyzu(2,i) = vxyzu(2,i) - 2.*vr*xyzmh(2,i)/r
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

            IF (xyzmh(3,i).LT.zmin) THEN
               ichan = ichan + 1
               iouter = iouter + 1               
               xyzmh(3,i) = xyzmh(3,i) + (zmax - zmin)
            ENDIF
            IF (xyzmh(3,i).GE.zmax) THEN
               ichan = ichan + 1
               iouter = iouter + 1               
               xyzmh(3,i) = xyzmh(3,i) - (zmax - zmin)
            ENDIF
            IF (iouter.NE.0) THEN
               idonebound = 1
            ENDIF
 450     CONTINUE
      ELSE IF ( ibound.EQ.13 ) THEN
         DO 500 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 500      
            xi = xyzmh(1,i)
            yi = xyzmh(2,i)
            zi = xyzmh(3,i)
c get cylindrical r
            rxy = sqrt(xi**2 + yi**2)
            IF (rxy.GT.tiny) THEN
               drxy = 1./rxy
            ELSE
               drxy = 0.
            ENDIF
c rintorus is radius from centre of torus
            rintorus2 = (rxy - Rtorus)**2 + zi**2
            IF (rintorus2.GE.atorus**2) THEN
               vxyzu(1,i) = 0.
               vxyzu(2,i) = 0.
               vxyzu(3,i) = 0.
               vxyzu(4,i) = 0.  !internal energy=0 => zero pressure grad
            ENDIF
c
c-- Reinject particles which escape in z direction, around the torus axis
c   
c            IF (ABS(zi).GE.1.5*atorus) THEN
c               ichan = ichan +1
c               iouter = iouter + 1
c               print *, iouter, i
c               delta = (0.5-ran1(1))*0.01
c               xyzmh(3,i) = delta
c               xyzmh(1,i) = (Rtorus+delta)*cos(2.*pi*ran1(1))
c               xyzmh(2,i) = (Rtorus+delta)*sin(2.*pi*ran1(1))
c               vxyzu(1,i) = 0.
c               vxyzu(2,i) = 0.
c               vxyzu(3,i) = 0.
c            ENDIF   
c            IF (ichan.NE.0) THEN
c               idonebound = 1
c            ENDIF         
 500     CONTINUE         
      ELSE IF ( ibound.EQ.103 ) THEN
c
c--Global disc boundaries, but with enforced sub-Keplerian motions near
c     inner and outer boundaries.
c
         DO 550 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 550
            delta = (1.3-ran1(1))*0.5*xyzmh(5,i)
            rcyld = rcyl - 0.05
            rmindd = rmind + 0.05
            zmind = zmin + delta
            zmaxd = zmax - delta
            ichan = 0
            r2 = xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i)

            IF (r2.GT.rcyld*rcyld) THEN
               iouter = iouter + 1
               r = SQRT(r2)
               yi = xyzmh(2,i)
               xi = xyzmh(1,i)
c
c--Set particles in the boundary to move with sub-Keplerian velocity
c
               vk = xmass/SQRT(r)
               vg = vk*SQRT((1.0-3.0*(vsound(i)/vk)**2))
               ang1 = ATAN2(yi,xi)
               vxyzu(1,i) = - vg*SIN(ang1)
               vxyzu(2,i) = vg*COS(ang1)
c               vxyzu(3,i) = 0.0
            ENDIF
            IF (r2.LT.rmindd*rmindd) THEN
               iinner = iinner + 1
               r = sqrt(r2)
               yi = xyzmh(2,i)
               xi = xyzmh(1,i)
c
c--Set particles in the boundary to move with sub-Keplerian velocity
c
               vk = xmass/SQRT(r)
               vg = vk*SQRT((1.0-3.0*(vsound(i)/vk)**2))
               ang1 = ATAN2(yi,xi)
               vxyzu(1,i) = - vg*SIN(ang1)
               vxyzu(2,i) = vg*COS(ang1)
c               vxyzu(3,i) = 0.0
            ENDIF

c            IF (xyzmh(3,i).LE.zmin) THEN
c               ichan = ichan + 1
c               vxyzu(3,i) = 0.
c               xyzmh(3,i) = zmind
c            ENDIF
c            IF (xyzmh(3,i).GE.zmax) THEN
c               ichan = ichan + 1
c               vxyzu(3,i) = 0.
c               xyzmh(3,i) = zmaxd
c            ENDIF
c
c--idonebound should be set if a particle has been moved as it 
c     affects the tree structure (changing its velocity doesn't).
c
            IF (ichan.NE.0) THEN
               iouter = iouter + 1
               idonebound = 1
            ENDIF
 550     CONTINUE

      ELSE IF ( ibound.EQ.104 ) THEN
c
c--Disc section boundaries.  Periodic in phi (moves particles from
c     one edge to the other.  But has enforced sub-Keplerian motions 
c     near inner and outer radial boundaries (like ibound=103)
c
         DO 600 j = 1, nlst0
            i = list(j)
c
c--Only do for gas and dust
c
            IF (iphase(i).LT.0 .OR.
     &           iphase(i).GT.0 .AND. iphase(i).LT.11) GOTO 600

            delta = (1.3-ran1(1))*0.5*xyzmh(5,i)

            rcyld = rcyl - 0.02
            rmindd = rmind + 0.02

            ichan = 0
            r2 = xyzmh(1,i)*xyzmh(1,i) + xyzmh(2,i)*xyzmh(2,i)
c
c--To avoid radial boundary
c
c            GOTO 777
c
c--Gas only for radial boundaries
c
            IF (iphase(i).EQ.0) THEN
               IF (r2.GT.rcyld*rcyld) THEN
                  iouter = iouter + 1
                  r = SQRT(r2)
                  yi = xyzmh(2,i)
                  xi = xyzmh(1,i)
c
c--Set particles in the boundary to move with sub-Keplerian velocity
c     NOTE: This will be wrong for large dust particles.
c
                  vk = xmass/SQRT(r) ! Keplerian
c                  vk = xmass*r  ! Omega = const
c                  vk = xmass*SQRT(r)
c                  vk = xmass

c                  vg = vk  !  *SQRT((1.0-3.0*(vsound(i)/vk)**2))
                  vg = vk*SQRT((1.0-3.0*(vsound(i)/vk)**2))

c                  vg = 0.
c                  phi = ATAN2(yi,xi)
c                  xyzmh(1,i) = 0.995*rcyl*COS(phi)
c                  xyzmh(2,i) = 0.995*rcyl*SIN(phi)
c                  GOTO 222

c
c--If calculation done in rotating reference frame, need to calculate
c     velocity in rotating frame
c
                  IF (ifcor.EQ.1) vg = vg - omeg0*r
                  ang1 = ATAN2(yi,xi)
                  vxyzu(1,i) = - vg*SIN(ang1)
                  vxyzu(2,i) = vg*COS(ang1)
c                  vxyzu(3,i) = 0.0
 222              CONTINUE
               ENDIF
               IF (r2.LT.rmindd*rmindd) THEN
                  iinner = iinner + 1
                  r = sqrt(r2)
                  yi = xyzmh(2,i)
                  xi = xyzmh(1,i)
c
c--Set particles in the boundary to move with sub-Keplerian velocity
c     NOTE: This will be wrong for large dust particles.
c
                  vk = xmass/SQRT(r)
c                  vk = xmass*r
c                  vk = xmass*SQRT(r)
c                  vk = xmass

c                  vg = vk    ! *SQRT((1.0-3.0*(vsound(i)/vk)**2))
                  vg = vk*SQRT((1.0-3.0*(vsound(i)/vk)**2))

c                  vg = 0.
c                  phi = ATAN2(yi,xi)
c                  xyzmh(1,i) = 1.005*rmind*COS(phi)
c                  xyzmh(2,i) = 1.005*rmind*SIN(phi)
c                  GOTO 223


c
c--If calculation done in rotating reference frame, need to calculate
c     velocity in rotating frame
c
                  IF (ifcor.EQ.1) vg = vg - omeg0*r
                  ang1 = ATAN2(yi,xi)
                  vxyzu(1,i) = - vg*SIN(ang1)
                  vxyzu(2,i) = vg*COS(ang1)
c                  vxyzu(3,i) = 0.0
 223              CONTINUE
               ENDIF
            ENDIF
c
c--Now do periodic boundaries in phi
c
c            GOTO 888

 777        phiparticle = ATAN2(xyzmh(2,i),xyzmh(1,i))
            IF (phiparticle.LT.-phibound) THEN
               rcylpart = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               vxpart = vxyzu(1,i)
               vypart = vxyzu(2,i)

               ichan = ichan + 1
               phinew = phiparticle + 2.0*phibound
               xyzmh(1,i) = rcylpart*COS(phinew)
               xyzmh(2,i) = rcylpart*SIN(phinew)
               cos2phi = COS(-2.0*phibound)
               sin2phi = SIN(-2.0*phibound)
               vxyzu(1,i) = vxpart*cos2phi + vypart*sin2phi
               vxyzu(2,i) = - vxpart*sin2phi + vypart*cos2phi
            ENDIF
            IF (phiparticle.GT.+phibound) THEN
               rcylpart = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               vxpart = vxyzu(1,i)
               vypart = vxyzu(2,i)

               ichan = ichan + 1
               phinew = phiparticle - 2.0*phibound
               xyzmh(1,i) = rcylpart*COS(phinew)
               xyzmh(2,i) = rcylpart*SIN(phinew)
               cos2phi = COS(2.0*phibound)
               sin2phi = SIN(2.0*phibound)
               vxyzu(1,i) = vxpart*cos2phi + vypart*sin2phi
               vxyzu(2,i) = - vxpart*sin2phi + vypart*cos2phi
            ENDIF
c
c--Periodic z boundaries if required (usually not)
c
            IF (.FALSE.) THEN
               IF (xyzmh(3,i).LT.zmin) THEN
                  ichan = ichan + 1
                  xyzmh(3,i) = xyzmh(3,i) + (zmax - zmin)
               ENDIF
               IF (xyzmh(3,i).GE.zmax) THEN
                  ichan = ichan + 1
                  xyzmh(3,i) = xyzmh(3,i) - (zmax - zmin)
               ENDIF
            ENDIF
c
c--idonebound should be set if a particle has been moved as it 
c     affects the tree structure (changing its velocity doesn't).
c
 888        IF (ichan.NE.0) THEN
               iouter = iouter + 1
               idonebound = 1
            ENDIF
 600     CONTINUE

      ENDIF

      IF (iinner.NE.0) WRITE (iprint,99009) iinner
      IF (iouter.NE.0) WRITE (iprint,99010) iouter

99009 FORMAT(' number of corrections inner boundary:',I6)
99010 FORMAT(' number of corrections outer boundary:',I6)

      IF (itrace.EQ.'all') WRITE (iprint, 99002)
99002 FORMAT (' exit subroutine boundary')

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
      
      IF (abs(dx).GT.0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
      IF (abs(dy).GT.0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
      IF (abs(dz).GT.0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)

      RETURN
      END SUBROUTINE modbound

c************************************************************
c                                                           *
c  same as modbound but also adjusts the difference in      *
c  euler potentials across the boundary                     *
c                                                           *
c************************************************************
      SUBROUTINE modboundeulr(dx,dy,dz,dalpha,dbeta)
      INCLUDE 'idim'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/presb'
      
      deltax = 0.
      deltay = 0.
      deltaz = 0.
      IF (abs(dx).GT.0.5*dxbound) THEN
         deltax = -dxbound*SIGN(1.0,dx)
         dx = dx + deltax
      ENDIF
      IF (abs(dy).GT.0.5*dybound) THEN
         deltay = -dybound*SIGN(1.0,dy)
         dy = dy + deltay
      ENDIF
      IF (abs(dz).GT.0.5*dzbound) THEN
         deltaz = -dzbound*SIGN(1.0,dz)
         dz = dz + deltaz
      ENDIF

c--   WARNING! DOES NOT DO MIXED CARTESIAN FIELDS!
      IF (abs(Bextz).GT.0.) THEN
         dalpha = dalpha - Bextz*deltay
         dbeta = dbeta + deltax
      ELSEIF (abs(Bexty).GT.0) THEN
         dalpha = dalpha - Bexty*deltax
         dbeta = dbeta + deltaz         
      ELSEIF (abs(Bextx).GT.0) THEN
         dalpha = dalpha - Bextx*deltaz
         dbeta = dbeta + deltay         
      ENDIF

      RETURN
      END SUBROUTINE modboundeulr


c************************************************************
c                                                           *
c  Avoid crossing boundaries, check on timestep             *
c                                                           *
c************************************************************

       FUNCTION dtcrossbound(xi,yi,zi,vxi,vyi,vzi,fxi,fyi,fzi,dt)
       
       IMPLICIT NONE
       
       INCLUDE 'idim'
       
       INCLUDE 'COMMONS/rbnd'
       INCLUDE 'COMMONS/typef'
       INCLUDE 'COMMONS/tokamak'
     
       REAL xi, yi, zi, vxi, vyi, vzi, dt, dtcrossbound
       REAL rr, drr, vr, deltar     
       REAL fxi, fyi, fzi, fr
       REAL costheta, sintheta, cosphi, sinphi, rintorus, rintorus2
       REAL drintorus, rrcyl, drcyl
       
       dtcrossbound = dt
       
       IF (iexf.EQ.9. OR. iBext.NE.0) THEN
          CALL get_torus_factors(xi,yi,zi,costheta,sintheta,cosphi,
     &              sinphi,rrcyl,drcyl,rintorus,rintorus2,drintorus)
          vr =  vxi*costheta*cosphi+vyi*costheta*sinphi+vzi*sintheta
          fr =  fxi*costheta*cosphi+fyi*costheta*sinphi+fzi*sintheta
          deltar = atorus - rintorus
          IF (fr.GT.0.) THEN
             dtcrossbound = 0.5*(sqrt(vr*vr+2.*fr*deltar)-vr)/(fr+tiny)
          ELSEIF (vr.GT.0.) THEN
             dtcrossbound = 0.5*deltar/(vr+tiny)
          ENDIF
       ELSEIF (iexf.EQ.10) THEN
          rr = sqrt(xi*xi+yi*yi)
          IF (rr.GT.tiny) THEN
             drr = 1./rr
          ELSE
             drr = 0.
          ENDIF    
          vr = (vxi*xi + vyi*yi)*drr
          fr = (fxi*xi+fyi*yi)*drr
          deltar = rcyl - rr
          IF (fr.GT.0.) THEN
             dtcrossbound = 0.5*(sqrt(vr*vr+2.*fr*deltar)-vr)/(fr+tiny)
          ELSEIF (vr.GT.0.) THEN
             dtcrossbound = 0.5*deltar/(vr+tiny)
          ENDIF  
       ENDIF       
       
       IF (dtcrossbound.LT.dt) print *, dtcrossbound, dt, vr, deltar
       
       RETURN
       END 
