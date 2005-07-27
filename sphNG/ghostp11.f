      SUBROUTINE ghostp11(npart, xyzmh, vxyzu)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)

      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/kerne'

      CHARACTER*7 where

      DATA where/'ghosp11'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp11')

      nghost = 0
c
c--Find ghost particles (for all particles within radkernel*h of boundary)
c
      DO 200 i = 1, npart
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
         delta = 0.0
c
c--Corners
c
         radk2 = radkernel*radkernel

         dxmin = (xi - xmin)/hi
         dxmax = (xmax - xi)/hi
         dymin = (yi - ymin)/hi
         dymax = (ymax - yi)/hi
         dzmin = (zi - zmin)/hi
         dzmax = (zmax - zi)/hi

         dxmin2 = dxmin*dxmin
         dymin2 = dymin*dymin
         dzmin2 = dzmin*dzmin
         dxmax2 = dxmax*dxmax
         dymax2 = dymax*dymax
         dzmax2 = dzmax*dzmax
         delta2 = delta*delta

         radius2iii = dxmin2 + dymin2 + dzmin2
         radius2iia = dxmin2 + dymin2 + dzmax2
         radius2iai = dxmin2 + dymax2 + dzmin2
         radius2iaa = dxmin2 + dymax2 + dzmax2
         radius2aii = dxmax2 + dymin2 + dzmin2
         radius2aia = dxmax2 + dymin2 + dzmax2
         radius2aai = dxmax2 + dymax2 + dzmin2
         radius2aaa = dxmax2 + dymax2 + dzmax2

         IF ((radius2iii.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2iii.LT.radk2) .OR.
     &      (radius2iia.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2iia.LT.radk2) .OR.
     &      (radius2iai.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2iai.LT.radk2) .OR.
     &      (radius2iaa.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2iaa.LT.radk2) .OR.
     &      (radius2aii.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2aii.LT.radk2) .OR.
     &      (radius2aia.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymin.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2aia.LT.radk2) .OR.
     &      (radius2aai.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmin.GT.delta .AND.
     &          radius2aai.LT.radk2) .OR.
     &      (radius2aaa.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &          dymax.GT.delta .AND. dzmax.GT.delta .AND.
     &          radius2aaa.LT.radk2)) THEN
            hasghost(i) = .TRUE.
            DO kk=1,7
               nghost = nghost + 1
               nptot = MIN0(npart + nghost, idim)
               ireal(nptot) = i
               xyzmh(1,nptot) = xi
               IF (kk.GT.3) THEN
                  IF (ABS(xi-xmax).LT.ABS(xi-xmin)) THEN
                     xyzmh(1,nptot) = xi-(xmax-xmin)
                  ELSE
                     xyzmh(1,nptot) = xi+(xmax-xmin)
                  ENDIF
               ENDIF
               xyzmh(2,nptot) = yi
               IF (MOD(kk/2,2).EQ.1) THEN
                  IF (ABS(yi-ymax).LT.ABS(yi-ymin)) THEN
                     xyzmh(2,nptot) = yi-(ymax-ymin)
                  ELSE
                     xyzmh(2,nptot) = yi+(ymax-ymin)
                  ENDIF
               ENDIF
               xyzmh(3,nptot) = zi
               IF (MOD(kk,2).EQ.1) THEN
                  IF (ABS(zi-zmax).LT.ABS(zi-zmin)) THEN
                     xyzmh(3,nptot) = zi-(zmax-zmin)
                  ELSE
                     xyzmh(3,nptot) = zi+(zmax-zmin)
                  ENDIF
               ENDIF

c         write (*,*) i, kk, nptot,xi,yi,zi,xyzmh(1,nptot),xyzmh(2,nptot),xyzmh(3,nptot)

               xyzmh(4,nptot) = pmassi
               xyzmh(5,nptot) = hi
               vxyzu(1,nptot) = vxi
               vxyzu(2,nptot) = vyi
               vxyzu(3,nptot) = vzi
               vxyzu(4,nptot) = ui
               rho(nptot) = rhoi
               iphase(nptot) = 0
            END DO
         ENDIF
         IF (hasghost(i)) GOTO 200
c
c--Edges
c
         radius2iio = dxmin2 + dymin2
         radius2iao = dxmin2 + dymax2
         radius2aio = dxmax2 + dymin2
         radius2aao = dxmax2 + dymax2
         radius2ioi = dxmin2 + dzmin2
         radius2ioa = dxmin2 + dzmax2
         radius2aoi = dxmax2 + dzmin2
         radius2aoa = dxmax2 + dzmax2
         radius2oii = dzmin2 + dymin2
         radius2oai = dzmin2 + dymax2
         radius2oia = dzmax2 + dymin2
         radius2oaa = dzmax2 + dymax2
         IF ((radius2iio.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &           dymin.GT.delta .AND. radius2iio.LT.radk2) .OR.
     &     (radius2iao.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &           dymax.GT.delta .AND. radius2iao.LT.radk2) .OR.
     &     (radius2aio.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &           dymin.GT.delta .AND. radius2aio.LT.radk2) .OR.
     &     (radius2aao.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &           dymax.GT.delta .AND. radius2aao.LT.radk2) .OR.
     &     (radius2ioi.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &           dzmin.GT.delta .AND. radius2ioi.LT.radk2) .OR.
     &     (radius2ioa.GT.delta2 .AND. dxmin.GT.delta .AND. 
     &           dzmax.GT.delta .AND. radius2ioa.LT.radk2) .OR.
     &     (radius2aoi.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &           dzmin.GT.delta .AND. radius2aoi.LT.radk2) .OR.
     &     (radius2aoa.GT.delta2 .AND. dxmax.GT.delta .AND. 
     &           dzmax.GT.delta .AND. radius2aoa.LT.radk2) .OR.
     &     (radius2oii.GT.delta2 .AND. dzmin.GT.delta .AND. 
     &           dymin.GT.delta .AND. radius2oii.LT.radk2) .OR.
     &     (radius2oia.GT.delta2 .AND. dzmin.GT.delta .AND. 
     &           dymax.GT.delta .AND. radius2oia.LT.radk2) .OR.
     &     (radius2oai.GT.delta2 .AND. dzmax.GT.delta .AND. 
     &           dymin.GT.delta .AND. radius2oai.LT.radk2) .OR.
     &     (radius2oaa.GT.delta2 .AND. dzmax.GT.delta .AND. 
     &           dymax.GT.delta .AND. radius2oaa.LT.radk2)) THEN
            hasghost(i) = .TRUE.
            DO kk=1,3
               nghost = nghost + 1
               nptot = MIN0(npart + nghost, idim)
               ireal(nptot) = i
               xyzmh(1,nptot) = xi
               IF (kk.GT.1) THEN
                  IF (radius2iio.LT.radk2 .OR. radius2iao.LT.radk2 .OR.
     &              radius2ioi.LT.radk2 .OR. radius2ioa.LT.radk2) THEN
                     xyzmh(1,nptot) = xi + (xmax-xmin)
               ELSEIF (radius2aio.LT.radk2 .OR. radius2aao.LT.radk2 .OR.
     &              radius2aoi.LT.radk2 .OR. radius2aoa.LT.radk2) THEN
                     xyzmh(1,nptot) = xi - (xmax-xmin)
                  ENDIF
               ENDIF
               xyzmh(2,nptot) = yi
               IF (MOD(kk,2).EQ.1) THEN
                  IF (radius2iio.LT.radk2 .OR. radius2aio.LT.radk2 .OR.
     &              radius2oii.LT.radk2 .OR. radius2oia.LT.radk2) THEN
                     xyzmh(2,nptot) = yi + (ymax-ymin)
               ELSEIF (radius2iao.LT.radk2 .OR. radius2aao.LT.radk2 .OR.
     &              radius2oai.LT.radk2 .OR. radius2oaa.LT.radk2) THEN
                     xyzmh(2,nptot) = yi - (ymax-ymin)
                  ENDIF
               ENDIF
               xyzmh(3,nptot) = zi
               IF (radius2ioi.LT.radk2 .OR. radius2aoi.LT.radk2 .OR.
     &              radius2ioa.LT.radk2 .OR. radius2aoa.LT.radk2) THEN
                  IF (MOD(kk,2).EQ.1) THEN
                  IF (radius2oii.LT.radk2 .OR. radius2oai.LT.radk2 .OR.
     &              radius2ioi.LT.radk2 .OR. radius2aoi.LT.radk2) THEN
                     xyzmh(3,nptot) = zi + (zmax-zmin)
               ELSEIF (radius2oia.LT.radk2 .OR. radius2oaa.LT.radk2 .OR.
     &              radius2ioa.LT.radk2 .OR. radius2aoa.LT.radk2) THEN
                     xyzmh(3,nptot) = zi - (zmax-zmin)
                  ENDIF
                  ENDIF
               ELSE
                  IF (kk.GT.1) THEN
                  IF (radius2oii.LT.radk2 .OR. radius2oai.LT.radk2 .OR.
     &              radius2ioi.LT.radk2 .OR. radius2aoi.LT.radk2) THEN
                     xyzmh(3,nptot) = zi + (zmax-zmin)
               ELSEIF (radius2oia.LT.radk2 .OR. radius2oaa.LT.radk2 .OR.
     &              radius2ioa.LT.radk2 .OR. radius2aoa.LT.radk2) THEN
                     xyzmh(3,nptot) = zi - (zmax-zmin)
                  ENDIF
                  ENDIF
               ENDIF

c         write (*,*) i, kk, nptot,xi,yi,zi,xyzmh(1,nptot),xyzmh(2,nptot),xyzmh(3,nptot)

               xyzmh(4,nptot) = pmassi
               xyzmh(5,nptot) = hi
               vxyzu(1,nptot) = vxi
               vxyzu(2,nptot) = vyi
               vxyzu(3,nptot) = vzi
               vxyzu(4,nptot) = ui
               rho(nptot) = rhoi
               iphase(nptot) = 0
            END DO
         ENDIF
         IF (hasghost(i)) GOTO 200
c
c--X axis
c
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
            iphase(nptot) = 0
         ENDIF

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
            iphase(nptot) = 0
         ENDIF
c
c--Y axis
c
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
            iphase(nptot) = 0
         ENDIF

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
            iphase(nptot) = 0
         ENDIF
c
c--Z axis
c
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
            iphase(nptot) = 0
         ENDIF

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
            iphase(nptot) = 0
         ENDIF
 200  CONTINUE

      ntot = npart + nghost

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
