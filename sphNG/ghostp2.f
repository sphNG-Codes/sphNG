      SUBROUTINE ghostp2(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(4,iradtrans)
      DIMENSION Bevolxyz(3,imhd)

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

      CHARACTER*7 where

      DATA where/'ghostp2'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp2')

      nghost = 0
      iinner = 0
c
c--If ibound = 2 then we have cylindrical boundaries
c--now find ghost for cylindrical boundaries
c
      rmind2 = rmind*rmind
      rcyl2 = rcyl*rcyl

      DO 300 i = 1, npart
         nghostold = nghost
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
         delta = 0.001*hi
         hi2 = hi*hi
         r2 = xi**2 + yi**2
         delta2 = delta*delta
         drmin = rcyl - radkernel*hi
         drmin2 = drmin*drmin
         drmax2 = rcyl2 + delta2 - 2*rcyl*delta
         IF (r2.GT.drmin2 .AND. r2.LT.drmax2 .OR. drmin.LT.0.) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            r = SQRT(r2)
            rgr = 2*rcyl - r
            xyzmh(1,nptot) = xi*rgr / r
            xyzmh(2,nptot) = yi*rgr / r
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            v2 = vxi*vxi + vyi*vyi
            vg = SQRT(v2)
            ang1 = ATAN2(yi,xi)
            ang2 = ATAN2(vyi,vxi)
            angr = pi - ang2 + ang1
            vxyzu(1,nptot) = vg*COS(angr)
            vxyzu(2,nptot) = vg*SIN(angr)
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         IF (rmind.NE.0) THEN
            drmax = rmind + radkernel*hi
            drref = rmind - radkernel*hi
            drmax2 = drmax*drmax
            drmin2 = rmind2 + delta2 + 2*rmind*delta
            IF (r2.GT.drmin2 .AND. r2.LT.drmax2 .AND. drref.GT.0.) THEN
               hasghost(i) = .TRUE.
               nghost = nghost + 1
               iinner = iinner + 1
               nptot = MIN0(npart + nghost, idim)
               ireal(nptot) = i
               r = SQRT(r2)
               rgr = 2*rmind - r
               xyzmh(1,nptot) = xi*rgr / r
               xyzmh(2,nptot) = yi*rgr / r
               xyzmh(3,nptot) = zi
               xyzmh(4,nptot) = pmassi
               xyzmh(5,nptot) = hi
               v2 = vxi*vxi + vyi*vyi
               vg = SQRT(v2)
               ang1 = ATAN2(yi,xi)
               ang2 = ATAN2(vyi,vxi)
               angr = pi - ang2 + ang1
               vxyzu(1,nptot) = vg*COS(angr)
               vxyzu(2,nptot) = vg*SIN(angr)
               vxyzu(3,nptot) = vzi
               vxyzu(4,nptot) = ui
               rho(nptot) = rhoi
               iphase(nptot) = 0
            ENDIF
         ENDIF
c
c--Z axis
c
         dzmin = (zi - zmin)/hi
         IF (dzmin.GT.delta .AND. dzmin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zmin - dzmin*hi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF

         dzmax = (zmax - zi)/hi
         IF (dzmax.GT.delta .AND. dzmax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zmax + dzmax*hi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF

         IF (nghostold.NE.nghost) THEN
            IF (encal.EQ.'r') THEN
               DO j=1,5
                  ekcle(j,nptot) = ekcle(j,i)
               END DO
            ENDIF
            IF (imhd.EQ.idim) THEN
               DO j=1,3
                  Bevolxyz(j,nptot) = Bevolxyz(j,i)
               END DO
            ENDIF
         ENDIF

 300  CONTINUE

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      ntot = npart + nghost
      IF (iinner.NE.0) WRITE (iprint,99090) iinner
99090 FORMAT(' adding ',I6,' inner ghosts')

      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
