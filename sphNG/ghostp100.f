      SUBROUTINE ghostp100(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries - for a part of a disc        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevolxyz(imhdevol,imhd)

      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/eosq'

      REAL*4 rhoreal4

      CHARACTER*7 where

      DATA where/'ghos100'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp100')

      nghost = 0
      iinner = 0
      uradconst = radconst/uergcc

      rmind = (1.0-variation)
      rmaxd = (1.0+variation)
      rmind2 = rmind*rmind
      rmaxd2 = rmaxd*rmaxd

      icentre = 0
      xcentre = 0.0
      DO 300 i = 1, npart
         nghostold = nghost
         nptot = MIN0(npart + nghost, idim)
         hasghost(i) = .FALSE.
         IF (iphase(i).NE.0) GOTO 300
         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         pmassi = xyzmh(4,i)
         hi = xyzmh(5,i)

         vxi = vxyzu(1,i)
         vyi = vxyzu(2,i)
         vzi = vxyzu(3,i)
         ui = vxyzu(4,i)
         rhoi = rho(i)
         vsoundi = vsound(i)
         presi = pr(i)

         delta = 0.001*hi
         hi2 = hi*hi
         r2 = xi**2 + yi**2 + zi**2
         r = SQRT(r2)
         r2centre = (xi-1.0)**2 + yi**2 + zi**2
         IF (r2centre.LT.xcentre) THEN
            xcentre = r2centre
            icentre = i
         ENDIF
c
c--Radial boundaries
c
         delta2 = delta*delta
         drmin = rmaxd - radkernel*hi
         drmin2 = drmin*drmin
         IF (r2.GT.drmin2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            rgr = 2*rmaxd - r
            xyzmh(1,nptot) = xi*rgr / r
            xyzmh(2,nptot) = yi*rgr / r
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vg = xmass/SQRT(rgr)
            ang1 = ATAN2(yi,xi)
            vxyzu(1,nptot) = - vg*SIN(ang1) + xyzmh(2,nptot)
            vxyzu(2,nptot) = vg*COS(ang1) - xyzmh(1,nptot)
            vxyzu(3,nptot) = 0.0
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
         ENDIF
         IF (rmind.GT.0) THEN
            drmax = rmind + radkernel*hi
            drmax2 = drmax*drmax
            IF (r2.LT.drmax2) THEN
               rgr = 2*rmind - r
               IF (rgr.GT.0.) THEN
                  hasghost(i) = .TRUE.
                  nghost = nghost + 1
                  iinner = iinner + 1
                  nptot = MIN0(npart + nghost, idim)
                  ireal(nptot) = i
                  xyzmh(1,nptot) = xi*rgr / r
                  xyzmh(2,nptot) = yi*rgr / r
                  xyzmh(3,nptot) = zi
                  xyzmh(4,nptot) = pmassi
                  xyzmh(5,nptot) = hi
                  vg = xmass/SQRT(rgr)
                  ang1 = ATAN2(yi,xi)
                  vxyzu(1,nptot) = - vg*SIN(ang1) + xyzmh(2,nptot)
                  vxyzu(2,nptot) = vg*COS(ang1) - xyzmh(1,nptot)
                  vxyzu(3,nptot) = 0.0
                  vxyzu(4,nptot) = ui
                  rho(nptot) = rhoi
                  vsound(nptot) = vsoundi
                  pr(nptot) = presi
                  iphase(nptot) = 0
               ENDIF
            ENDIF
         ENDIF
c
c--Phi boundaries
c
         phimax = phibound - radkernel*hi/r
         phimin = - phibound + radkernel*hi/r
         ang1 = ATAN2(yi,xi)
         IF (ang1.GT.phimax) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            iinner = iinner + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            ang2 = 2*phibound - ang1
            xyzmh(1,nptot) = r*COS(ang2)
            xyzmh(2,nptot) = r*SIN(ang2)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vg = xmass/SQRT(r)
            vxyzu(1,nptot) = - vg*SIN(ang2) + xyzmh(2,nptot)
            vxyzu(2,nptot) = vg*COS(ang2) - xyzmh(1,nptot)
            vxyzu(3,nptot) = 0.0
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
         ENDIF
         IF (ang1.LT.phimin) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            iinner = iinner + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            ang2 = - 2*phibound - ang1
            xyzmh(1,nptot) = r*COS(ang2)
            xyzmh(2,nptot) = r*SIN(ang2)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vg = xmass/SQRT(r)
            vxyzu(1,nptot) = - vg*SIN(ang2) + xyzmh(2,nptot)
            vxyzu(2,nptot) = vg*COS(ang2) - xyzmh(1,nptot)
            vxyzu(3,nptot) = 0.0
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
         ENDIF
c
c--Z axis does not have a boundary.  Do not worry about corners - assume
c     flow near corners is supersonic so pressure boundaries shouldn't
c     matter much.
c
         IF (nghostold.NE.nghost) THEN
            DO k=nptot-(nghost-nghostold)+1,nptot
               IF (encal.EQ.'r') THEN
                  DO j=1,5
                     ekcle(j,k) = ekcle(j,i)
                  END DO
                  rad = SQRT(xyzmh(1,k)**2 + xyzmh(2,k)**2 +
     &                 xyzmh(3,k)**2)
                  vxyzu(4,k) = hoverr**2/(rad*(gamma-1.0))
c
c--Assumes 75 g/cm^2 surface density at planet's radius, multiplied by signorm
c
                  profile = 0.5
      rhoreal4 = signorm*75.0/(SQRT(2.0*pi)*hoverr*rad*udist)/
     &        (rad**profile)*EXP(-xyzmh(3,i)**2/
     &        (2.0*(hoverr*rad)**2))/umass*udist**3
c         ekcle(3,i) = getcv(rhoreal4,vxyzu(4,i))
c         ekcle(1,i) = uradconst*(vxyzu(4,i)/ekcle(3,i))**4/rhoreal4
c         ekcle(2,i) = getkappa(vxyzu(4,i),ekcle(3,i),rhoreal4)
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j=1,imhdevol
                     Bevolxyz(j,k) = Bevolxyz(j,i)
                  END DO
               ENDIF
            ENDDO
         ENDIF
         
 300  CONTINUE
c
c--Add boundary for planet's surface -closest particle to planet is icentre
c
c--But might be better to do it with repulsive boundary condition otherwise
c     particles may drift inside planet
c
      GOTO 400
      rplanet = 2.4E-05
      IF (xyzmh(5,icentre).GT.rplanet) THEN
         nghost = nghost + 1
         nptot = MIN0(npart + nghost, idim)
         ireal(nptot) = icentre
         xyzmh(1,nptot) = 1.0
         xyzmh(2,nptot) = 0.0
         xyzmh(3,nptot) = 0.0
         xyzmh(4,nptot) = xyzmh(4,icentre)
         xyzmh(5,nptot) = rplanet/2.0
         vxyzu(1,nptot) = 0.0
         vxyzu(2,nptot) = 0.0
         vxyzu(3,nptot) = 0.0
         vxyzu(4,nptot) = vxyzu(4,icentre)
         rho(nptot) = rho(icentre)
         iphase(nptot) = 0

         IF (encal.EQ.'r') THEN
            DO j=1,5
               ekcle(j,nptot) = ekcle(j,icentre)
            END DO
         ENDIF
      ENDIF

 400  WRITE (iprint, *) 'npart, nghost', npart, nghost
      ntot = npart + nghost
      IF (iinner.NE.0) WRITE (iprint,99090) iinner
99090 FORMAT(' adding ',I6,' inner ghosts')

      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
