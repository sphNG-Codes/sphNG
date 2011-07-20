      SUBROUTINE ghostp1(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
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

      CHARACTER*7 where

      DATA where/'ghostp1'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp1')

      nghost = 0
      uradconst = radconst/uergcc
c
c--Find ghost particles (for all particles within radkernel*h of boundary)
c
      DO 200 i = 1, npart
         nghostold = nghost
         nptot = MIN0(npart + nghost, idim)
         hasghost(i) = .FALSE.
         IF (iphase(i).NE.0) GOTO 200
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

         delta = 0.1*hi
c
c--X axis
c
         dxmin = (xi - xmin)/hi
         IF (dxmin.GT.delta .AND. dxmin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF

         dxmax = (xmax - xi)/hi
         IF (dxmax.GT.delta .AND. dxmax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
c
c--Y axis
c
         dymin = (yi - ymin)/hi
         IF (dymin.GT.delta .AND. dymin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF

         dymax = (ymax - yi)/hi
         IF (dymax.GT.delta .AND. dymax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
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
            xyzmh(3,nptot) = 2.0*zmin - zi
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
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
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
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF

         radius2 = dxmin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dzmin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &                    dzmax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dzmin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &                    dzmax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
 
         radius2 = dzmin2 + dymin2
         IF (radius2.GT.delta2 .AND. dzmin.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dzmin2 + dymax2
         IF (radius2.GT.delta2 .AND. dzmin.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dzmax2 + dymin2
         IF (radius2.GT.delta2 .AND. dzmax.GT.delta .AND. 
     &                    dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dzmax2 + dymax2
         IF (radius2.GT.delta2 .AND. dzmax.GT.delta .AND. 
     &                    dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
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
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmin - xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = 2.0*ymin - yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = 2.0*zmin - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2.LT.radk2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = 2.0*xmax - xi
            xyzmh(2,nptot) = 2.0*ymax - yi
            xyzmh(3,nptot) = 2.0*zmax - zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = -vxi
            vxyzu(2,nptot) = -vyi
            vxyzu(3,nptot) = -vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
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
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j=1,imhdevol
                     Bevolxyz(j,k) = Bevolxyz(j,i)
                  END DO
               ENDIF
            ENDDO
         ENDIF

 200  CONTINUE

      ntot = npart + nghost

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
