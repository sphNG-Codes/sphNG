      SUBROUTINE inform(where)
c************************************************************
c                                                           *
c  This routine computes relevant quantities for print out  *
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

#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/out1'
      INCLUDE 'COMMONS/out2'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/nonideal'
      INCLUDE 'COMMONS/eosq'
#endif

      CHARACTER*7 where, where2
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine inform')
c
c--Compute : center of mass, velocity of cm, mean density and dispersion
c     for first object
c
      uradconst = radconst/uergcc

      cmx1 = 0.
      cmy1 = 0.
      cmz1 = 0.
      vcmx1 = 0.
      vcmy1 = 0.
      vcmz1 = 0.
      cmx2 = 0.
      cmy2 = 0.
      cmz2 = 0.
      vcmx2 = 0.
      vcmy2 = 0.
      vcmz2 = 0.

      romean1 = 0.
      romax1 = 0.
      rocen1 = 0.
      romean2 = 0.
      romax2 = 0.
      rocen2 = 0.

      tgmean1 = 0.
      tgmax1 = 0.
      tgcen1 = 0.
      tgmean2 = 0.
      tgmax2 = 0.
      tgcen2 = 0.

      trmean1 = 0.
      trmax1 = 0.
      trcen1 = 0.
      trmean2 = 0.
      trmax2 = 0.
      trcen2 = 0.
      valphamax1 = -1.0
      valphamin1 = 10.0
      valphamax2 = -1.0
      valphamin2 = 10.0
      hmi1 = 1.E30
      hma1 = 0.
      hmi2 = 1.E30
      hma2 = 0.

      h2mass1 = 0.
      tdmax1 = 0.
      tdmin1 = 1.0E+10
      tdmean1 = 0.

      h2mass2 = 0.
      tdmax2 = 0.
      tdmin2 = 1.0E+10
      tdmean2 = 0.

      fmas1 = 0.
      fmas2 = 0.

      n1new = 0
      n2new = 0.
      DO i = 1, npart
         IF (iunique(iorig(i)).LE.n1 .AND. 
     &                    (iphase(i).EQ.0 .OR. iphase(i).GE.10)) THEN
            n1new = n1new + 1
            pmassi = xyzmh(4,i)
            fmas1 = fmas1 + pmassi
            cmx1 = cmx1 + pmassi*xyzmh(1,i)
            cmy1 = cmy1 + pmassi*xyzmh(2,i)
            cmz1 = cmz1 + pmassi*xyzmh(3,i)

            vcmx1 = vcmx1 + pmassi*vxyzu(1,i)
            vcmy1 = vcmy1 + pmassi*vxyzu(2,i)
            vcmz1 = vcmz1 + pmassi*vxyzu(3,i)

            hmi1 = MIN(hmi1, xyzmh(5,i))
            hma1 = MAX(hma1, xyzmh(5,i))

            romean1 = romean1 + rho(i)
            romax1 = MAX(romax1, rho(i))

            valphamax1 = MAX(valphamax1, alphaMM(1,i))
            valphamin1 = MIN(valphamin1, alphaMM(1,i))

            IF (encal.EQ.'r') THEN
               tempgas = vxyzu(4,i)/ekcle(3,i)
               tgmean1 = tgmean1 + tempgas
               tgmax1 = MAX(tgmax1, tempgas)
               temprad = (ekcle(1,i)*rho(i)/uradconst)**0.25
               trmean1 = trmean1 + temprad
               trmax1 = MAX(trmax1, temprad)

               IF (idustRT.GT.0) THEN
                  tdmean1 = tdmean1 + dust_tk(1,i)
                  tdmax1 = MAX(tdmax1, dust_tk(1,i))
                  tdmin1 = MIN(tdmin1, dust_tk(1,i))
c
c--The factor of two is because h2frac=0.5 when fully molecular
c
                  h2mass1 = h2mass1 + pmassi*2.0*h2frac(i)
               ENDIF

            ELSEIF (encal.EQ.'m') THEN
               tempgas = vxyzu(4,i)/getcv(rho(i),vxyzu(4,i))
               tgmean1 = tgmean1 + tempgas
               tgmax1 = MAX(tgmax1, tempgas)
            ELSEIF (gamma.EQ.1.0) THEN
               tempgas = gmw*uergg/Rg*2.0/3.0*vxyzu(4,i)
               tgmean1 = tgmean1 + tempgas
               tgmax1 = MAX(tgmax1, tempgas)
            ELSE
               tempgas = gmw*uergg/Rg*(gamma-1.0)*vxyzu(4,i)
               tgmean1 = tgmean1 + tempgas
               tgmax1 = MAX(tgmax1, tempgas)
            ENDIF
         ENDIF
      END DO
      DO j = 1, nptmass
         i = listpm(j)
         pmassi = xyzmh(4,i)
         fmas1 = fmas1 + pmassi
         cmx1 = cmx1 + pmassi*xyzmh(1,i)
         cmy1 = cmy1 + pmassi*xyzmh(2,i)
         cmz1 = cmz1 + pmassi*xyzmh(3,i)

         vcmx1 = vcmx1 + pmassi*vxyzu(1,i)
         vcmy1 = vcmy1 + pmassi*vxyzu(2,i)
         vcmz1 = vcmz1 + pmassi*vxyzu(3,i)
      END DO

      cmx1 = cmx1/fmas1
      cmy1 = cmy1/fmas1
      cmz1 = cmz1/fmas1

      vcmx1 = vcmx1/fmas1
      vcmy1 = vcmy1/fmas1
      vcmz1 = vcmz1/fmas1

      IF (n1new.GT.0) THEN
         romean1 = romean1/n1new
         tgmean1 = tgmean1/n1new
         trmean1 = trmean1/n1new
      ENDIF
      IF (idustRT.GT.0) THEN
         tdmean1 = tdmean1/n1new
      ENDIF
c
c--If doing accretion on to binary, then recentre the centre of mass and
c     set the centre of mass velocity to zero
c
      IF (ibound.EQ.8 .OR. ibound/10.EQ.9) THEN
         WRITE (iprint,*) 'ZERO CENTRE OF MASS'
         DO i=1, npart
            IF (iphase(i).GE.0) THEN
               xyzmh(1,i) = xyzmh(1,i) - cmx1
               xyzmh(2,i) = xyzmh(2,i) - cmy1
               xyzmh(3,i) = xyzmh(3,i) - cmz1
               vxyzu(1,i) = vxyzu(1,i) - vcmx1
               vxyzu(2,i) = vxyzu(2,i) - vcmy1
               vxyzu(3,i) = vxyzu(3,i) - vcmz1
            ENDIF
         END DO
      ENDIF
c
c--Compute : center of mass, velocity of cm, mean density and dispersion
c     for second object (if existing)
c
      IF (n2.NE.0) THEN
         DO i = 1, npart
            IF (iunique(iorig(i)).GT.n1 .AND. 
     &           iunique(iorig(i)).LE.n1+n2 .AND. iphase(i).GE.0) THEN
               pmassi = xyzmh(4,i)
               n2new = n2new + 1
               fmas2 = fmas2 + xyzmh(4,i)
               cmx2 = cmx2 + pmassi*xyzmh(1,i)
               cmy2 = cmy2 + pmassi*xyzmh(2,i)
               cmz2 = cmz2 + pmassi*xyzmh(3,i)

               vcmx2 = vcmx2 + pmassi*vxyzu(1,i)
               vcmy2 = vcmy2 + pmassi*vxyzu(2,i)
               vcmz2 = vcmz2 + pmassi*vxyzu(3,i)

               hmi2 = MIN(hmi2, xyzmh(5,i))
               hma2 = MAX(hma2, xyzmh(5,i))

               romean2 = romean2 + rho(i)
               romax2 = MAX(romax2, rho(i))

               valphamax2 = MAX(valphamax2, alphaMM(1,i))
               valphamin2 = MIN(valphamin2, alphaMM(1,i))
               IF (encal.EQ.'r') THEN
                  tempgas = vxyzu(4,i)/ekcle(3,i)
                  tgmean2 = tgmean2 + tempgas
                  tgmax2 = MAX(tgmax2, tempgas)
                  temprad = (ekcle(1,i)*rho(i)/uradconst)**0.25
                  trmean2 = trmean2 + temprad
                  trmax2 = MAX(trmax2, temprad)

                  IF (idustRT.GT.0) THEN
                     tdmean2 = tdmean2 + dust_tk(1,i)
                     tdmax2 = MAX(tdmax2, dust_tk(1,i))
                     tdmin2 = MIN(tdmin2, dust_tk(1,i))

                     h2mass2 = h2mass2 + pmassi*2.0*h2frac(i)
                  ENDIF

               ELSEIF (encal.EQ.'m') THEN
                  tempgas = vxyzu(4,i)/getcv(rho(i),vxyzu(4,i))
                  tgmean2 = tgmean2 + tempgas
                  tgmax2 = MAX(tgmax2, tempgas)
               ELSEIF (gamma.EQ.1.0) THEN
                  tempgas = gmw*uergg/Rg*2.0/3.0*vxyzu(4,i)
                  tgmean2 = tgmean2 + tempgas
                  tgmax2 = MAX(tgmax2, tempgas)
               ELSE
                  tempgas = gmw*uergg/Rg*(gamma-1.0)*vxyzu(4,i)
                  tgmean2 = tgmean2 + tempgas
                  tgmax2 = MAX(tgmax2, tempgas)
               ENDIF
            ENDIF
         END DO

         cmx2 = cmx2/fmas2
         cmy2 = cmy2/fmas2
         cmz2 = cmz2/fmas2

         vcmx2 = vcmx2/fmas2
         vcmy2 = vcmy2/fmas2
         vcmz2 = vcmz2/fmas2

         romean2 = romean2/n2new
         tgmean2 = tgmean2/n2new
         trmean2 = trmean2/n2new
         IF (idustRT.GT.0) THEN
            tdmean2 = tdmean2/n1new
         ENDIF
      ENDIF
c
c--Compute maximum distance
c
      dmax1 = 0.
      zmax1 = 0.
      rmin2 = 1.0E+30
      imin = 0
      rocen1 = 0.
      DO i = 1, npart
         IF (iunique(iorig(i)).LE.n1 .AND. iphase(i).GE.0) THEN
#ifdef MPI
            dz = xyzmh(3,i) - zlocationmpi
            d2 = (xyzmh(1,i) - xlocationmpi)**2 + 
     &           (xyzmh(2,i) - ylocationmpi)**2 + dz*dz
#else
            dz = xyzmh(3,i) - cmz1
            d2 = (xyzmh(1,i) - cmx1)**2 + (xyzmh(2,i) - cmy1)**2 +dz*dz
#endif
            dmax1 = MAX(dmax1, d2)
            dmin1 = MIN(dmin1, d2)
            zmax1 = MAX(zmax1, ABS(dz))
            IF (d2.LT.rmin2 .AND. 
     &                      (iphase(i).EQ.0 .OR. iphase(i).GE.10)) THEN
               imin = i
               rmin2 = d2
            ENDIF
         ENDIF
      END DO
      dmax1 = SQRT(dmax1)
      rmin1 = SQRT(rmin2)
      IF (imin.NE.0) THEN
         rocen1 = rho(imin)
         IF (encal.EQ.'r') THEN
            tgcen1 = vxyzu(4,imin)/ekcle(3,imin)
            trcen1 = (ekcle(1,imin)*rocen1/uradconst)**0.25
         ELSEIF (encal.EQ.'m') THEN
            tgcen1 = vxyzu(4,imin)/getcv(rho(imin),vxyzu(4,imin))
         ENDIF
      ENDIF

      dmax2 = 0.
      zmax2 = 0.
      rmin2 = 1.0E+30
      imin = 0
      rocen2 = 0.
      DO i = 1, npart
         IF (iunique(iorig(i)).GT.n1 .AND. iphase(i).GE.0) THEN
#ifdef MPI
            dz = xyzmh(3,i) - zlocationmpi
            d2 = (xyzmh(1,i) - xlocationmpi)**2 + 
     &           (xyzmh(2,i) - ylocationmpi)**2 + dz*dz
#else
            dz = xyzmh(3,i) - cmz2
            d2 = (xyzmh(1,i) - cmx2)**2 + (xyzmh(2,i) - cmy2)**2 +dz*dz
#endif
            dmax2 = MAX(dmax2, d2)
            zmax2 = MAX(zmax2, ABS(dz))
            IF (d2.LT.rmin2 .AND. 
     &                      (iphase(i).EQ.0 .OR. iphase(i).GE.10)) THEN
               imin = i
               rmin2 = d2
            ENDIF
         ENDIF
      END DO
      dmax2 = SQRT(dmax2)
      IF (imin.NE.0) rocen2 = rho(imin)
c
c--Compute energies
c
      CALL toten
c
c--Compute total angular momentum
c
      CALL angmom
c
c--Compute various parameters to do with the magnetic field
c
      IF (imhd.EQ.idim) CALL mhdparams
c
c--Compute various properties to do with non-ideal MHD
#ifdef NONIDEAL
      CALL nimhd_get_stats(imhd2,npart,xyzmh,vxyzu,ekcle,Bxyz,rho,
     &               vsound,
     &               alphaMM,alphamin,n_R,n_electronT,iphase,gamma,
     &               ifsvi,encal,np_nimhd,et_nimhd,evdata,
     &               ionfrac_eta,n1,n2,iunique,iorig)
#endif
c
      IF (where(1:6).NE.'newrun') THEN
c
c--Write dump on disk
c
         IF (nstep.LT.1) nstep = 1
c         IF (nbuild.EQ.1 .OR. MOD(ncount, nstep).EQ.0 .OR. 
         IF (MOD(ncount, nstep).EQ.0 .OR. 
     &                       iptcreat.EQ.1) THEN
#ifdef MPI
            IF (iproc.EQ.0) CALL file
#else
            CALL file
#endif
            CALL wdump_wrapper(idisk1)
         ENDIF

#ifdef DGRID
c
c--Write out 3D grid of drag deposited energy.
c
         IF (gasdrag) THEN
#ifdef MPICOPY
            IF (iproc.EQ.0) THEN
#endif
            CALL write_dragtable(idragfile)
#ifdef MPICOPY
            ENDIF
#endif
         ENDIF
#endif
c
c--Write global results on listing
c
         where2 = 'inform'
         CALL prout(where2)
         iptcreat = 0
c
c--Update input file
c
         IF (ifulldump.EQ.0) CALL wrinsph

      ELSE
c
c--Write global results on listing
c
         where2 = 'inform'
         CALL prout(where2)         
      ENDIF

      IF (idebug(1:6).EQ.'inform') THEN
         WRITE (iprint, 99002) cmx1, cmy1, cmz1, vcmx1, vcmy1, vcmz1,
     &                         romean1, zmax1
99002    FORMAT (1X, 5(1PE12.5,1X))
      ENDIF

      RETURN
      END
