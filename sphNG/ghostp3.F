      SUBROUTINE ghostp3(npart, xyzmh, vxyzu, ekcle, Bevolxyz, dustvar)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevolxyz(imhdevol,imhd)
      DIMENSION dustvar(ndusttypes,idim_dustFluid)

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
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/stellarradiation'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
      REAL    jcurrenti(3),nden_nimhdi(n_nden),eta_nimhdi(4)
#endif

      DIMENSION dustvari(ndusttypes)

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
         rhoreal4 = rho(i)
         vsoundi = vsound(i)
         presi = pr(i)
         IF (idustFluid.NE.0) dustvari(:) = dustvar(:,i)

         IF (istellar_ionisation) THEN
            HIIioni = HIIion(i)
         ENDIF

         delta = 0.001*hi
         hi2 = hi*hi
         r2 = xi*xi + yi*yi + zi*zi + tiny
         delta2 = delta*delta
         drmin = rmax - radkernel*hi
         drmin2 = drmin*drmin
         drmax2 = rmax2 + delta2 - 2*rmax*delta
#ifdef NONIDEAL
         jcurrenti(1:3)        = jcurrent(1:3,i)
         nden_nimhdi(1:n_nden) = nden_nimhd(1:n_nden,i)
         eta_nimhdi(1:4)       = eta_nimhd(1:4,i)
#endif

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
            rho(nptot) = rhoreal4
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
            IF (idustFluid.NE.0) dustvar(:,nptot) = dustvari(:)

            IF (istellar_ionisation) THEN
               HIIion(nptot) = HIIioni
            ENDIF
#ifdef NONIDEAL
            jcurrent(1:3,nptot)  = jcurrenti(1:3)
            nden_nimhd(1:n_nden,nptot) = nden_nimhdi(1:n_nden)
            eta_nimhd(1:4,nptot) = eta_nimhdi(1:4)
#endif
         ENDIF

         IF (nghostold.NE.nghost) THEN
            IF (encal.EQ.'r') THEN
               IF (rhoreal4.EQ.0.0) THEN
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
                  ekcle(1,nptot) = uradconst*(boundtemp)**4/rhoreal4
                  ekcle(2,nptot) = getkappa(vxyzu(4,nptot),
     &                 ekcle(3,nptot),rhoreal4)
               ENDIF
               ekcle(4,nptot) = ekcle(4,i)
               ekcle(5,nptot) = ekcle(5,i)
               IF (idustRT.GT.0) THEN
                  dust_tk(1,nptot) = dust_tk(1,i)
                  dust_tk(2,nptot) = dust_tk(2,i)
               ENDIF
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

