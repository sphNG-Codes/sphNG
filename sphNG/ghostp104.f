      SUBROUTINE ghostp104(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries - for a part of a disc        *
c     Similar to ibound=100, except that 104 uses periodic  *
c        boundaries in phi whereas 100 uses reflected       *
c        ghosts with enforced Keplerian rotation.           *
c     Assumes that the disc section is centred on (1,0,0)   *
c     Unlike most ghost particle routines, this creates     *
c        ghosts of dust particles as well.
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
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/rotat'

      REAL*4 rhoreal4

      CHARACTER*7 where

      DATA where/'ghos104'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp104')

      IF (ifcor.GE.2) THEN
         WRITE (*,*) 'ERROR - ghostp104 with ifcor = ',ifcor
         CALL quit(0)
      ENDIF

      nghost = 0
      iinner = 0
      uradconst = radconst/uergcc

      rmind = (1.0-variation)
      rmaxd = (1.0+variation)
      rmind2 = rmind*rmind
      rmaxd2 = rmaxd*rmaxd

      DO 300 i = 1, npart
         nghostold = nghost
         nptot = MIN0(npart + nghost, idim)
         hasghost(i) = .FALSE.
c
c--Only create ghosts of gas and dust
c
         IF (iphase(i).LT.0 .OR. 
     &        iphase(i).GT.0 .AND. iphase(i).LT.11) GOTO 300
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
         iphasei = iphase(i)

         r2 = xi**2 + yi**2
         r = SQRT(r2)
c
c--Radial boundaries
c
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
c
c--Assume Keplerian velocity structure at boundary
c
            vg = xmass/SQRT(rgr)
c            vg = xmass*(rgr)
c            vg = xmass*SQRT(rgr)
c            vg = xmass

            ang1 = ATAN2(yi,xi)
            vxyzu(1,nptot) = - vg*SIN(ang1) + ifcor*xyzmh(2,nptot)
            vxyzu(2,nptot) = vg*COS(ang1) - ifcor*xyzmh(1,nptot)
c            vxyzu(1,nptot) = 0.0
c            vxyzu(2,nptot) = 0.0

            vxyzu(3,nptot) = 0.0
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = iphasei
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
c
c--Assume Keplerian velocity structure at boundary
c
                  vg = xmass/SQRT(rgr)
c                  vg = xmass*(rgr)
c                  vg = xmass*SQRT(rgr)
c                  vg = xmass

                  ang1 = ATAN2(yi,xi)
                  vxyzu(1,nptot) = - vg*SIN(ang1) + ifcor*xyzmh(2,nptot)
                  vxyzu(2,nptot) = vg*COS(ang1) - ifcor*xyzmh(1,nptot)
c                  vxyzu(1,nptot) = 0.0
c                  vxyzu(2,nptot) = 0.0

                  vxyzu(3,nptot) = 0.0
                  vxyzu(4,nptot) = ui
                  rho(nptot) = rhoi
                  vsound(nptot) = vsoundi
                  pr(nptot) = presi
                  iphase(nptot) = iphasei
               ENDIF
            ENDIF
         ENDIF
c
c--Phi boundaries (periodic in phi)
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
            ang2 = ang1 - 2.0*phibound
            xyzmh(1,nptot) = r*COS(ang2)
            xyzmh(2,nptot) = r*SIN(ang2)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi

            sin2phi = SIN(2.0*phibound)
            cos2phi = COS(2.0*phibound)
            vxyzu(1,nptot) = vxi*cos2phi + vyi*sin2phi
            vxyzu(2,nptot) = -vxi*sin2phi + vyi*cos2phi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = iphasei
         ENDIF
         IF (ang1.LT.phimin) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            iinner = iinner + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            ang2 = ang1 + 2.0*phibound
            xyzmh(1,nptot) = r*COS(ang2)
            xyzmh(2,nptot) = r*SIN(ang2)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi

            sin2phi = SIN(-2.0*phibound)
            cos2phi = COS(-2.0*phibound)
            vxyzu(1,nptot) = vxi*cos2phi + vyi*sin2phi
            vxyzu(2,nptot) = -vxi*sin2phi + vyi*cos2phi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = iphasei
         ENDIF
c
c--Z axis does not have a boundary.  Do not worry about corners - assume
c     flow near corners is supersonic so pressure boundaries shouldn't
c     matter much.
c
         GOTO 888
c
c--If want to run a section with no vertical gravity, make reflective
c     vertical boundaries by uncommenting the above GOTO
c
         dzmin = (zi - zmin)/hi
         IF (dzmin.GT.0. .AND. dzmin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF

         dzmax = (zmax - zi)/hi
         IF (dzmax.GT.0. .AND. dzmax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
c
c--Skip the rest if not gas
c
 888     IF (iphase(nptot).NE.0) GOTO 300

         IF (nghostold.NE.nghost) THEN
            DO k=nptot-(nghost-nghostold)+1,nptot
               IF (encal.EQ.'r') THEN
                  DO j=1,5
                     ekcle(j,k) = ekcle(j,i)
                  END DO
                  rad = SQRT(xyzmh(1,k)**2 + xyzmh(2,k)**2 +
     &                 xyzmh(3,k)**2)
                  vxyzu(4,k) = hoverr**2/(rad*gamma*(gamma-1.0))
c
c--Assumes 75 g/cm^2 surface density at planet's radius, multiplied by signorm
c
                  profile = 0.5
c         rhoreal4 = signorm*75.0/(SQRT(2.0*pi)*hoverr*rad*udist)/
c     &        (rad**profile)*EXP(-xyzmh(3,i)**2/
c     &        (2.0*(hoverr*rad)**2))/umass*udist**3
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

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      ntot = npart + nghost
      IF (iinner.NE.0) WRITE (iprint,99090) iinner
99090 FORMAT(' adding ',I6,' inner ghosts')

      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
