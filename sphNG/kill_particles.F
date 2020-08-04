       SUBROUTINE kill_particles (i,xyzmh,vxyzu,time,ikilled)
c************************************************************
c                                                           *
c  This subroutine kills a particle that has exceeded       *
c     boundaries specified by the boundary type, ibound     *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif

      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/binfile'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/planetesimal'

      INTEGER i, ikilled
      REAL time, xyzmh(5,idim), vxyzu(4,idim)

      REAL r2, pmassi, phiparticle, rlim
      INTEGER ipkill, ig

      REAL rp, pvel, theta, pvx, pvy, dvx, dvy, dvz, vmagnitude2
      REAL dragr, dragp, dragt, eccentricity, smaxis
      INTEGER idragr, idragp, idragt

      LOGICAL killme
      CHARACTER*7 where

      DATA where/'killp'/
c
c--Dead particle boundaries - constant N boundaries
c     (infall onto accreting binary; Bate 1995; Bate & Bonnell 1997)
c
      IF (ibound.EQ.8 .OR. ibound/10.EQ.9) THEN
         r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
         IF (r2.GT.(1.05*deadbound)**2) THEN
            ipkill = iphase(i)
            iphase(i) = -1
            pmassi = xyzmh(4,i)
#ifndef MPICOPY
C$OMP CRITICAL (killparticle)
#endif
            ikilled = 1

            nlstacc = nlstacc + 1
            IF (nlstacc.GT.nlstaccmax) THEN
               WRITE (iprint,*) 'ERROR step nlstacc'
               CALL quit(1)
            ENDIF
            listacc(nlstacc) = i

            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

            nactive = nactive - 1
            nkill = nkill + 1
            WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &           xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           poten(i),dgrav(i),ipkill,0.0,0.0,0.0,0.0,0.0,0.0
            CALL FLUSH (ikillpr)
            anglostx = anglostx + pmassi*(xyzmh(2,i)*vxyzu(3,i) - 
     &           vxyzu(2,i)*xyzmh(3,i))
            anglosty = anglosty + pmassi*(vxyzu(1,i)*xyzmh(3,i) - 
     &           xyzmh(1,i)*vxyzu(3,i))
            anglostz = anglostz + pmassi*(xyzmh(1,i)*vxyzu(2,i) - 
     &           vxyzu(1,i)*xyzmh(2,i))
#ifndef MPICOPY
C$OMP END CRITICAL (killparticle)
#endif
            RETURN
         ENDIF
      ENDIF

c
c--Dead particle boundaries - planet embedded in disc
c
      IF (ibound.EQ.100 .AND. iexf.EQ.7) THEN
         phiparticle = ATAN2(xyzmh(2,i),xyzmh(1,i))
         r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
         IF (r2.LT.(1.0-variation)**2 .OR. r2.GT.(1.0+variation)**2
     &        .OR. ABS(phiparticle).GT.phibound) THEN
            ipkill = iphase(i)
            iphase(i) = -1
#ifndef MPICOPY
C$OMP CRITICAL (killparticle)
#endif
            ikilled = 1

            nlstacc = nlstacc + 1
            IF (nlstacc.GT.nlstaccmax) THEN
               WRITE (iprint,*) 'ERROR step nlstacc'
               CALL quit(1)
            ENDIF
            listacc(nlstacc) = i

            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

            nactive = nactive - 1
            nkill = nkill + 1
            WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &           xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           poten(i),dgrav(i),ipkill,0.0,0.0,0.0,0.0,0.0,0.0
            CALL FLUSH(ikillpr)
#ifndef MPICOPY
C$OMP END CRITICAL (killparticle)
#endif
            RETURN
         ENDIF
      ENDIF

c
c--Dead particle boundaries - disc section, no planet, periodic phi
c
      IF (ibound.EQ.104) THEN
         r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
         IF (r2.LT.(1.0-variation)**2 .OR. r2.GT.(1.0+variation)**2
     &        ) THEN
            ipkill = iphase(i)
            iphase(i) = -1
#ifndef MPICOPY
C$OMP CRITICAL (killparticle)
#endif
            ikilled = 1

            nlstacc = nlstacc + 1
            IF (nlstacc.GT.nlstaccmax) THEN
               WRITE (iprint,*) 'ERROR step nlstacc'
               CALL quit(1)
            ENDIF
            listacc(nlstacc) = i

            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

            nactive = nactive - 1
            nkill = nkill + 1
            WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &           xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           poten(i),dgrav(i),ipkill,0.0,0.0,0.0,0.0,0.0,0.0
            CALL FLUSH(ikillpr)
#ifndef MPICOPY
C$OMP END CRITICAL (killparticle)
#endif
            RETURN
         ENDIF
      ENDIF

c
c--Dead particle boundaries - planet embedded in whole disc
c
      IF (ibound.EQ.102 .OR. ibound.EQ.103 .OR. ibound.EQ.105) THEN
         killme = .false.
c
c--Note: Ben Ayliffe used spherical radius, whereas Pablo Loren-Aguilar
c     uses cylindrical radius.  Ben also did not kill gas inside rmind,
c     and he did not kill planetesimals at r>rcyl or r<min_rplan
c
         IF (ibound.EQ.102 .OR. ibound.EQ.105) THEN
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
         ELSE
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
         ENDIF
         IF (iphase(i).EQ.0) THEN
            IF (ibound.EQ.102) THEN
               IF (r2.GT.rcyl**2) killme = .true.
            ELSEIF (ibound.EQ.103 .OR. ibound.EQ.105) THEN
               IF (r2.GT.rcyl**2 .OR. r2.LT.rmind**2) killme = .true.
            ELSE
               IF (r2.LT.rmind**2) killme = .true.
            ENDIF
         ELSEIF (iphase(i).GE.11) THEN
            CALL semimajoraxis(i, eccentricity, smaxis)
            IF (smaxis.LT.min_rplan .OR. smaxis.GT.max_rplan)
     &           killme = .true.
            IF (ibound.EQ.103) THEN
               IF (r2.GT.rcyl**2.OR.r2.LT.min_rplan**2)
     &              killme = .true.
            ENDIF
         ENDIF
         IF (killme) THEN
            ipkill = iphase(i)
            iphase(i) = -1
#ifndef MPICOPY
C$OMP CRITICAL (killparticle)
#endif
            ikilled = 1

            nlstacc = nlstacc + 1
            IF (nlstacc.GT.nlstaccmax) THEN
               WRITE (iprint,*) 'ERROR step nlstacc'
               CALL quit(1)
            ENDIF
            listacc(nlstacc) = i

            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

            nactive = nactive - 1
            nkill = nkill + 1
            WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &           xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           poten(i),dgrav(i),ipkill,0.0,0.0,0.0,0.0,0.0,0.0
            CALL FLUSH(ikillpr)
#ifndef MPICOPY
C$OMP END CRITICAL (killparticle)
#endif
            RETURN
         ENDIF
      ENDIF
c
c--Kill planetesimals approaching planet modelled as potential.
c
      IF ((irotpot.EQ.1 .OR. iexf.EQ.7) .AND. iplanetesimals.GT.0
     &     .AND. iphase(i).GE.11) THEN
         IF (irotpot.EQ.1) r2 = (xyzmh(1,i)-px)**2 +
     &        (xyzmh(2,i)-py)**2 + xyzmh(3,i)**2
         IF (iexf.EQ.7) r2 = (xyzmh(1,i)-1.0)**2 +
     &        xyzmh(2,i)**2 + xyzmh(3,i)**2

         rlim = rplanet*pradfac(1)
         IF (r2.LT.rlim**2) THEN
            ipkill = iphase(i)
            iphase(i) = -1
#ifndef MPICOPY
C$OMP CRITICAL (killparticle)
#endif
            ikilled = 1
c
c--Add solids mass into core as it is accreted.
c
            planetmass = planetmass + xyzmh(4,i)
            coremass = planetmass

            IF (encal.EQ.'r') THEN
c
c--Nearest gas particle to site of impact to which energy can be given
c  is stored in nearestgas().
                  
               ig = nearestgas(i)
c
c--Calculate planetesimal KE relative to planetary core
c
               rp = SQRT(px*px + py*py)
               pvel = SQRT(xmass/rp)
               theta = ATAN2(py, px)
               pvx = -pvel*sin(theta)
               pvy = pvel*cos(theta)

               dvx = vxyzu(1,i) - pvx
               dvy = vxyzu(2,i) - pvy
               dvz = vxyzu(3,i)
               vmagnitude2 = dvx*dvx + dvy*dvy + dvz*dvz
               pcolen(ig) = 0.5*xyzmh(4,i)*vmagnitude2/xyzmh(4,ig)

#ifdef DGRID
               IF (gasdrag) THEN
                  CALL sphericalpos(xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),
     &                 dragr, dragp, dragt)
                  CALL dragcoordinates(dragr, dragp, dragt, idragr,
     &                 idragp, idragt)
                  IF (idragr.GT.1 .AND.
     &                 (pradfac(1)-1.0).LT.1.0E-8) THEN
                     print *, 'Solids killed at surface not going in'
                     print *, 'the innermost bin. Weird.'
                     print *, 'idragr = ', idragr, dragr
                     print *, 'pradfac = ', pradfac(1), SQRT(r2)
                     print *, 'rplanet = ', rplanet, rplanet*
     &                    pradfac(1)
                     CALL quit(1) 
                  ENDIF

c---DragEnergyTest output. Checking coordinate transform and indexing.
c                     write (44, 146) xyzmh(1,i), xyzmh(2,i), xyzmh(3,i),
c     &                    px, py, pz, dragr, dragp, dragt, idragr,
c     &                    idragp, idragt
c 146                 FORMAT (9(1PE15.8,1X),3(I4,1X))

c
c--Actually always want the colision energy in inner most bin, even
c  when the planet is shrinking, therefore:
                  idragr = 1
c--Makes MAX statements below redundant
                  dragenergy(MAX(idragr,1), idragp, idragt) = 
     &                 dragenergy(MAX(idragr,1), idragp, idragt) +
     &                 pcolen(ig)
               ENDIF
#endif
            ENDIF

            nlstacc = nlstacc + 1
            IF (nlstacc.GT.nlstaccmax) THEN
               WRITE (iprint,*) 'ERROR step nlstacc'
               CALL quit(1)
            ENDIF
            listacc(nlstacc) = i

            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

            nactive = nactive - 1
            nkill = nkill + 1
            WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &           xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           poten(i),dgrav(i),ipkill,px,py,pcolen(ig),
     &           xyzmh(1,ig),xyzmh(2,ig),xyzmh(3,ig)
            CALL FLUSH(ikillpr)
#ifndef MPICOPY
C$OMP END CRITICAL (killparticle)
#endif
            RETURN
         ENDIF
      ENDIF

 999  RETURN
      END
