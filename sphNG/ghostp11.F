      SUBROUTINE ghostp11(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
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
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/raddust'

      CHARACTER*7 where

      DATA where/'ghost11'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp11')

      nghost = 0
#ifdef PERIODIC_NO_GHOSTS
      RETURN
#endif
      uradconst = radconst/uergcc
c
c--Find ghost particles (for all particles within radkernel*h of boundary)
c
c
c  For periodic boundary conditions, safest is to use hmax
c
      hmax = 0.
      DO i = 1, npart
         hmax = max(2.1*xyzmh(5,i),hmax)
      ENDDO

      DO 200 i = 1, npart
         nghostold = nghost
         nptot = MIN0(npart + nghost, idim)
         hasghost(i) = .FALSE.
         IF (iphase(i).GE.1 .AND. iphase(i).LE.10) GOTO 200
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

         delta = 0.0
c
c--X axis
c
         dxmin = (xi - xmin)/hmax
         IF (dxmin.GT.delta .AND. dxmin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF

         dxmax = (xmax - xi)/hmax
         IF (dxmax.GT.delta .AND. dxmax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
c
c--Y axis
c
         dymin = (yi - ymin)/hmax
         IF (dymin.GT.delta .AND. dymin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF

         dymax = (ymax - yi)/hmax
         IF (dymax.GT.delta .AND. dymax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
c
c--Z axis
c
         dzmin = (zi - zmin)/hmax
         IF (dzmin.GT.delta .AND. dzmin.LT.radkernel) THEN
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
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF

         dzmax = (zmax - zi)/hmax
         IF (dzmax.GT.delta .AND. dzmax.LT.radkernel) THEN
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
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
c
c--Edges
c
         dxmin2 = dxmin*dxmin
         dymin2 = dymin*dymin
         dzmin2 = dzmin*dzmin
         dxmax2 = dxmax*dxmax
         dymax2 = dymax*dymax
         dzmax2 = dzmax*dzmax
         delta2 = delta*delta

         radius2 = dxmin2 + dymin2
         radk2 = radkernel*radkernel
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmin2 + dymax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dymin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dymax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF

         radius2 = dxmin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dzmin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dzmax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dzmin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dzmax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
 
         radius2 = dzmin2 + dymin2
         IF (radius2.GT.delta2 .AND. dzmin.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dzmin2 + dymax2
         IF (radius2.GT.delta2 .AND. dzmin.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dzmax2 + dymin2
         IF (radius2.GT.delta2 .AND. dzmax.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dzmax2 + dymax2
         IF (radius2.GT.delta2 .AND. dzmax.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
c
c--Corners
c
         radius2 = dxmin2 + dymin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmin2 + dymin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmin2 + dymax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmin2 + dymax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi + (xmax-xmin)
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dymin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dymin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi + (ymax-ymin)
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dymax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF
         radius2 = dxmax2 + dymax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi - (xmax-xmin)
            xyzmh(2,nptot) = yi - (ymax-ymin)
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = iphase(i)
            IF (imhd.EQ.idim) CALL copyBevol(nptot,i,Bevolxyz,
     &         xyzmh(1,nptot)-xi,xyzmh(2,nptot)-yi,xyzmh(3,nptot)-zi)
         ENDIF

         IF (nghostold.NE.nghost) THEN
             DO k=nptot-(nghost-nghostold)+1,nptot
               IF (encal.EQ.'r') THEN
                  DO j=1,5
                     ekcle(j,k) = ekcle(j,i)
                  END DO
                  vxyzu(4,k) = 0.704097133431896
                  ekcle(1,k) = uradconst*(vxyzu(4,k)/
     &              ekcle(3,k))**4/50.226017
                  IF (idustRT.GT.0) THEN
                     dust_tk(1,k) = dust_tk(1,i)
                     dust_tk(2,k) = dust_tk(2,i)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

 200  CONTINUE

      ntot = npart + nghost

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      IF (itrace.EQ.'all') WRITE (iprint, 99002)
99002 FORMAT (' exit subroutine ghostp11')

      RETURN
      END
      
      SUBROUTINE copyBevol(inew,iold,Bevolxyz,deltax,deltay,deltaz)
      INCLUDE 'idim'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'

      DIMENSION Bevolxyz(imhdevol,imhd)
c
c--for Euler potentials need gradient to be continuous across boundary
c      
      IF (varmhd.EQ.'eulr') THEN
c--   WARNING! DOES NOT DO MIXED CARTESIAN FIELDS YET
         IF (abs(Bextz).GT.0.) THEN
            Bevolxyz(1,inew) = Bevolxyz(1,iold) - Bextz*deltay
            Bevolxyz(2,inew) = Bevolxyz(2,iold) + deltax
         ELSEIF (abs(Bexty).GT.0) THEN
            Bevolxyz(1,inew) = Bevolxyz(1,iold) - Bexty*deltax
            Bevolxyz(2,inew) = Bevolxyz(2,iold) + deltaz         
         ELSEIF (abs(Bextx).GT.0) THEN
            Bevolxyz(1,inew) = Bevolxyz(1,iold) - Bextx*deltaz
            Bevolxyz(2,inew) = Bevolxyz(2,iold) + deltay         
         ENDIF
      ELSE
         DO j=1,imhdevol
            Bevolxyz(j,inew) = Bevolxyz(j,iold)
         END DO
      ENDIF
      
      RETURN
      END
