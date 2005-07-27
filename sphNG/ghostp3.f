      SUBROUTINE ghostp3(npart, xyzmh, vxyzu)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/crpart'
      INCLUDE 'COMMONS/kerne'

      CHARACTER*7 where

      DATA where/'ghostp3'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp3')

      nghost = 0
      inshell = 0
c
c--If ibound = 3, 8, or 9 then we have spherical boundaries
c--now find ghost for spherical boundaries
c
      rmax2 = rmax*rmax
      DO 300 i = 1, npart
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
         r2 = xi*xi + yi*yi + zi*zi + tiny
         delta2 = delta*delta
         drmin = rmax - radkernel*hi
         drmin2 = drmin*drmin
         drmax2 = rmax2 + delta2 - 2*rmax*delta

         IF (r2.GT.rshell*rshell) inshell = inshell + 1

         IF (r2.GT.drmin2 .AND. r2.LT.drmax2) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            r = SQRT(r2)
            rgr = 2*rmax - r
            xir = xi/r
            yir = yi/r
            zir = zi/r
            xyzmh(1,nptot) = xir*rgr
            xyzmh(2,nptot) = yir*rgr
            xyzmh(3,nptot) = zir*rgr
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vr = vxi*xir + vyi*yir + vzi*zir
            vxr = vr*xir
            vyr = vr*yir
            vzr = vr*zir
            vxyzu(1,nptot) = vxi - 2*vxr
            vxyzu(2,nptot) = vyi - 2*vyr
            vxyzu(3,nptot) = vzi - 2*vzr
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            iphase(nptot) = 0
         ENDIF

 300  CONTINUE

      WRITE (iprint, *) 'nghost ', nghost

      ntot = npart + nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
