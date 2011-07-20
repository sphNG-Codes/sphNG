      SUBROUTINE ghostp3(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
c                                                           *
c************************************************************

      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif

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
      INCLUDE 'COMMONS/crpart'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/cgas'

      CHARACTER*7 where
      REAL*4 rhoreal4

      DATA where/'ghostp3'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp3')

      nghost = 0
      inshell = 0
c
c--Set boundary temperature for radiative transfer
c
c      boundtemp=11.9
c      boundtemp=12.037
c      boundtemp=10.115

      uradconst = radconst/uergcc
c
c--If ibound = 3, 8, or 9 then we have spherical boundaries
c--now find ghost for spherical boundaries
c
      rmax2 = rmax*rmax
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
         r2 = xi*xi + yi*yi + zi*zi + tiny
         delta2 = delta*delta
         drmin = rmax - radkernel*hi
         drmin2 = drmin*drmin
         drmax2 = rmax2 + delta2 - 2*rmax*delta

         IF (r2.GT.rshell*rshell) inshell = inshell + 1

         IF (r2.GT.drmin2 .AND. r2.LT.drmax2) THEN
c            print *,i,r2,hi
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

         IF (nghostold.NE.nghost) THEN
            IF (encal.EQ.'r') THEN
               IF (rhoi.EQ.0.0) THEN
                  ekcle(3,nptot) = 1.5*Rg/(gmw*uergg)
                  vxyzu(4,nptot) = boundtemp*ekcle(3,nptot)
                  ekcle(1,nptot) = uradconst*(vxyzu(4,nptot)/
     &                 ekcle(3,nptot))**4/(1.0/(4.0*pi/3.0*3.2**3))
                  rhoreal4 = (1.0/(4.0*pi/3.0*3.2**3))
                  ekcle(2,nptot) = getkappa(vxyzu(4,nptot),
     &                 ekcle(3,nptot),rhoreal4)
               ELSE
                  ekcle(3,nptot) = 1.5*Rg/(gmw*uergg)
c                  vxyzu(4,nptot) = boundtemp*ekcle(3,nptot)
                  ekcle(1,nptot) = uradconst*(boundtemp)**4/rhoi
                  ekcle(2,nptot) = getkappa(vxyzu(4,nptot),
     &                 ekcle(3,nptot),rho(i))
               ENDIF
               ekcle(4,nptot) = ekcle(4,i)
               ekcle(5,nptot) = ekcle(5,i)
            ENDIF
            IF (imhd.EQ.idim) THEN
               DO j=1,imhdevol
                  Bevolxyz(j,nptot) = Bevolxyz(j,i)
               END DO
            ENDIF
         ENDIF

 300  CONTINUE

#ifdef MPIALL
      IF (iproc.EQ.0) WRITE (iprint, *) 'nghost ', nghost
#else
      WRITE (iprint, *) 'nghost ', nghost
#endif

      ntot = npart + nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
