      SUBROUTINE ghostp102(npart,xyzmh,vxyzu,ekcle,Bevolxyz,dustvar)
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
      DIMENSION dustvar(idim_dustFluid)

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

c      REAL*4 rhoreal4

      CHARACTER*7 where

      DATA where/'ghos102'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp102')

      nghost = 0
      iinner = 0
      uradconst = radconst/uergcc


      rmaxd = rcyl
      rmind2 = rmind*rmind
      rmaxd2 = rmaxd*rmaxd

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
         vsoundi = vsound(i)
         presi = pr(i)
         IF (idustFluid.NE.0) dustvari = dustvar(i)

         r2 = xi**2 + yi**2 + zi**2
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
            vg = xmass/SQRT(rgr)
            ang1 = ATAN2(yi,xi)
            vxyzu(1,nptot) = - vg*SIN(ang1)
            vxyzu(2,nptot) = vg*COS(ang1)
            vxyzu(3,nptot) = 0.0
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            IF (idustFluid.NE.0) dustvar(nptot) = dustvari
            iphase(nptot) = 0
         ENDIF
         goto 234
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
                  vxyzu(1,nptot) = - vg*SIN(ang1)
                  vxyzu(2,nptot) = vg*COS(ang1)
                  vxyzu(3,nptot) = 0.0
                  vxyzu(4,nptot) = ui
                  rho(nptot) = rhoi
                  vsound(nptot) = vsoundi
                  pr(nptot) = presi
                  IF (idustFluid.NE.0) dustvar(nptot) = dustvari
                  iphase(nptot) = 0
               ENDIF
            ENDIF
         ENDIF
 234     continue
c
c--Z axis does not have a boundary.  Do not worry about corners,
c  assume flow near corners is supersonic so pressure boundaries
c  shouldn't matter much.

         IF (nghostold.NE.nghost) THEN
            DO k=nptot-(nghost-nghostold)+1,nptot
               IF (encal.EQ.'r') THEN
                  DO j=1,5
                     ekcle(j,k) = ekcle(j,i)
                  END DO
                  rad = SQRT(xyzmh(1,k)**2 + xyzmh(2,k)**2 +
     &                 xyzmh(3,k)**2)
                  IF (use_tprof) THEN
                     vxyzu(4,k) = centralmass*
     &                    (hoverr**2*rad**tprof)/(gamma-1.0)
                  ELSE
                     vxyzu(4,k) = centralmass*
     &                    hoverr**2/(rad*(gamma-1.0))
                  ENDIF
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO j=1,imhdevol
                     Bevolxyz(j,k) = Bevolxyz(j,i)
                  END DO
               ENDIF
            ENDDO
         ENDIF
         
 300  CONTINUE

c      OPEN (unit=36, file='ghostout.txt', access='append')
c      DO i = 1, npart+nghost
c         rlocal = sqrt(xyzmh(1,i)**2+xyzmh(2,i)**2+xyzmh(3,i)**2)
c         write(36,999) rlocal, vxyzu(1,i)*(-xyzmh(2,i)/rlocal) + 
c     &        vxyzu(2,i)*(xyzmh(1,i)/rlocal), xyzmh(1,i),
c     &        xyzmh(2,i),xyzmh(3,i), vxyzu(4,i)
c      END DO
c      CLOSE (36)
c 999  FORMAT (6(1PE12.5,2X))

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      ntot = npart + nghost
      IF (iinner.NE.0) WRITE (iprint,99090) iinner
99090 FORMAT(' adding ',I6,' inner ghosts')

      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
