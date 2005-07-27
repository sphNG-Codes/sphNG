      SUBROUTINE revtree (nnatom, npart, xyzmh)
c************************************************************
c                                                           *
c  Subroutine by W. Press.  Revises the tree structure.     *
c     Assumes that only the positions have been altered     *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,idim)

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/logun'

      natom = nnatom

c
c--REVISE ENTIRE TREE (STANDARD REVTREE)
c
      IF (nlst.GT.nnatom/5000 .OR. (.NOT. ipartialrevtree)) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(natom,npart,imfac)
C$OMP& shared(nactatom,listmap,iphase)
C$OMP& private(j,ipart)

C$OMP DO SCHEDULE(static)
      DO j = 1, nactatom
         ipart = listmap(j)
         IF (ipart.GT.npart .OR. iphase(ipart).EQ.-1 .OR.
     &               (iphase(ipart).GE.1 .AND. iptintree.EQ.1)) THEN
            imfac(ipart) = 0
         ELSE
            imfac(ipart) = 1
         ENDIF
      END DO
C$OMP END DO
C$OMP END PARALLEL

      DO ilevel = 1, nlevel

         IF ((level(ilevel + 1)-level(ilevel)).GT.16) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(ilevel,level,natom,isibdaupar)
C$OMP& shared(qrad,xyzmh,imfac)
C$OMP& private(new,l,ll,fl,fll,emred,difx,dify,difz,rr,pmassl,pmassll)

C$OMP DO SCHEDULE(static)
         DO new = level(ilevel), level(ilevel + 1) - 1

            l = isibdaupar(2,new)
            ll = isibdaupar(1,l)

            IF (l.LE.natom) THEN
               pmassl = imfac(l)*xyzmh(4,l)
            ELSE
               pmassl = xyzmh(4,l)
            ENDIF
            IF (ll.LE.natom) THEN
               pmassll = imfac(ll)*xyzmh(4,ll)
            ELSE
               pmassll = xyzmh(4,ll)
            ENDIF
            xyzmh(4,new) = pmassl + pmassll

            IF (xyzmh(4,new).NE.0) THEN
               fl = pmassl/xyzmh(4,new)
               fll = pmassll/xyzmh(4,new)
            ELSE
               fl = 0.5
               fll = 0.5
            ENDIF
            emred = fl*fll*xyzmh(4,new)
            difx = xyzmh(1,ll) - xyzmh(1,l)
            dify = xyzmh(2,ll) - xyzmh(2,l)
            difz = xyzmh(3,ll) - xyzmh(3,l)
            IF (fl.GT.fll) THEN
               xyzmh(1,new) = xyzmh(1,l) + fll*difx
               xyzmh(2,new) = xyzmh(2,l) + fll*dify
               xyzmh(3,new) = xyzmh(3,l) + fll*difz
            ELSE
               xyzmh(1,new) = xyzmh(1,ll) - fl*difx
               xyzmh(2,new) = xyzmh(2,ll) - fl*dify
               xyzmh(3,new) = xyzmh(3,ll) - fl*difz
            ENDIF
c
c--Find radius
c
            rr = SQRT(difx**2 + dify**2 + difz**2) + tiny
            qrad(1,new) = MAX(fll*rr + qrad(1,l), fl*rr + qrad(1,ll))
            xyzmh(5,new) = MAX(xyzmh(5,l), xyzmh(5,ll))
c
c--Find quadrupole moments
c
            qrad(2,new) = (emred*difx)*difx
            qrad(5,new) = (emred*difx)*dify
            qrad(7,new) = (emred*difx)*difz
            qrad(3,new) = (emred*dify)*dify
            qrad(6,new) = (emred*dify)*difz
            qrad(4,new) = (emred*difz)*difz
            IF (l.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,l)
               qrad(3,new) = qrad(3,new) + qrad(3,l)
               qrad(4,new) = qrad(4,new) + qrad(4,l)
               qrad(5,new) = qrad(5,new) + qrad(5,l)
               qrad(6,new) = qrad(6,new) + qrad(6,l)
               qrad(7,new) = qrad(7,new) + qrad(7,l)
            ENDIF
            IF (ll.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,ll)
               qrad(3,new) = qrad(3,new) + qrad(3,ll)
               qrad(4,new) = qrad(4,new) + qrad(4,ll)
               qrad(5,new) = qrad(5,new) + qrad(5,ll)
               qrad(6,new) = qrad(6,new) + qrad(6,ll)
               qrad(7,new) = qrad(7,new) + qrad(7,ll)
            ENDIF
         END DO
C$OMP END DO
C$OMP END PARALLEL

         ELSE

         DO new = level(ilevel), level(ilevel + 1) - 1

            l = isibdaupar(2,new)
            ll = isibdaupar(1,l)

            IF (l.LE.natom) THEN
               pmassl = imfac(l)*xyzmh(4,l)
            ELSE
               pmassl = xyzmh(4,l)
            ENDIF
            IF (ll.LE.natom) THEN
               pmassll = imfac(ll)*xyzmh(4,ll)
            ELSE
               pmassll = xyzmh(4,ll)
            ENDIF
            xyzmh(4,new) = pmassl + pmassll

            IF (xyzmh(4,new).NE.0) THEN
               fl = pmassl/xyzmh(4,new)
               fll = pmassll/xyzmh(4,new)
            ELSE
               fl = 0.5
               fll = 0.5
            ENDIF
            emred = fl*fll*xyzmh(4,new)
            difx = xyzmh(1,ll) - xyzmh(1,l)
            dify = xyzmh(2,ll) - xyzmh(2,l)
            difz = xyzmh(3,ll) - xyzmh(3,l)
            IF (fl.GT.fll) THEN
               xyzmh(1,new) = xyzmh(1,l) + fll*difx
               xyzmh(2,new) = xyzmh(2,l) + fll*dify
               xyzmh(3,new) = xyzmh(3,l) + fll*difz
            ELSE
               xyzmh(1,new) = xyzmh(1,ll) - fl*difx
               xyzmh(2,new) = xyzmh(2,ll) - fl*dify
               xyzmh(3,new) = xyzmh(3,ll) - fl*difz
            ENDIF
c
c--Find radius
c
            rr = SQRT(difx**2 + dify**2 + difz**2) + tiny
            qrad(1,new) = MAX(fll*rr + qrad(1,l), fl*rr + qrad(1,ll))
            xyzmh(5,new) = MAX(xyzmh(5,l), xyzmh(5,ll))
c
c--Find quadrupole moments
c
            qrad(2,new) = (emred*difx)*difx
            qrad(5,new) = (emred*difx)*dify
            qrad(7,new) = (emred*difx)*difz
            qrad(3,new) = (emred*dify)*dify
            qrad(6,new) = (emred*dify)*difz
            qrad(4,new) = (emred*difz)*difz
            IF (l.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,l)
               qrad(3,new) = qrad(3,new) + qrad(3,l)
               qrad(4,new) = qrad(4,new) + qrad(4,l)
               qrad(5,new) = qrad(5,new) + qrad(5,l)
               qrad(6,new) = qrad(6,new) + qrad(6,l)
               qrad(7,new) = qrad(7,new) + qrad(7,l)
            ENDIF
            IF (ll.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,ll)
               qrad(3,new) = qrad(3,new) + qrad(3,ll)
               qrad(4,new) = qrad(4,new) + qrad(4,ll)
               qrad(5,new) = qrad(5,new) + qrad(5,ll)
               qrad(6,new) = qrad(6,new) + qrad(6,ll)
               qrad(7,new) = qrad(7,new) + qrad(7,ll)
            ENDIF
         END DO

         ENDIF
      END DO
c
c--ONLY REVISE PART OF THE TREE - CURRENT PARTICLES AND THEIR NEIGHBOURS
c     AND ANY ACCRETED OR KILLED PARTICLES FROM THE LAST STEP
c
      ELSE

      IF (nlst.GT.50) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(nlst,llist,iflagtree,nneigh,neighb,neighover)
C$OMP& private(i,j,k,ipart,jpart,kpart)
C$OMP DO SCHEDULE(static)
         DO i = 1, nlst
            ipart = llist(i)
C$OMP CRITICAL (flagtree)
            iflagtree(ipart) = .TRUE.
C$OMP END CRITICAL (flagtree)
            DO j = 1, nneigh(ipart)
               IF (j.GE.nlmax) THEN
                  jpart = neighover(j-nlmax+1,ABS(neighb(nlmax,ipart)))
               ELSE
                  jpart = neighb(j,ipart)
               ENDIF
C$OMP CRITICAL (flagtree)
               iflagtree(jpart) = .TRUE.
C$OMP END CRITICAL (flagtree)
               DO k = 1, nneigh(jpart)
                  IF (k.GE.nlmax) THEN
                  kpart = neighover(k-nlmax+1,ABS(neighb(nlmax,jpart)))
                  ELSE
                     kpart = neighb(k,jpart)
                  ENDIF
C$OMP CRITICAL (flagtree)
                  iflagtree(kpart) = .TRUE.
C$OMP END CRITICAL (flagtree)
               END DO
            END DO
         END DO
C$OMP END DO
C$OMP END PARALLEL

      ELSE

         DO i = 1, nlst
            ipart = llist(i)
            iflagtree(ipart) = .TRUE.
            DO j = 1, nneigh(ipart)
               IF (j.GE.nlmax) THEN
                  jpart = neighover(j-nlmax+1,ABS(neighb(nlmax,ipart)))
               ELSE
                  jpart = neighb(j,ipart)
               ENDIF
               iflagtree(jpart) = .TRUE.
C$OMP PARALLEL DO SCHEDULE(static) default(none)
C$OMP& shared(jpart,nneigh,neighb,neighover,iflagtree)
C$OMP& private(k,kpart)
               DO k = 1, nneigh(jpart)
                  IF (k.GE.nlmax) THEN
                  kpart = neighover(k-nlmax+1,ABS(neighb(nlmax,jpart)))
                  ELSE
                     kpart = neighb(k,jpart)
                  ENDIF
                  iflagtree(kpart) = .TRUE.
               END DO
C$OMP END PARALLEL DO
            END DO
         END DO
      ENDIF

      IF (nptmass.GT.10) THEN
         
C$OMP PARALLEL default(none)
C$OMP& shared(nptmass,nptlist,nearpt,iflagtree)
C$OMP& private(i,j,jpart)
C$OMP DO SCHEDULE(static)
         DO i = 1, nptmass
            DO j = 1, nptlist(i)
               jpart = nearpt(j,i)
C$OMP CRITICAL (flagtree)
               iflagtree(jpart) = .TRUE.
C$OMP END CRITICAL (flagtree)
            END DO
         END DO
C$OMP END DO
C$OMP END PARALLEL

      ELSE
         DO i = 1, nptmass
            nneighipart = nptlist(i)
C$OMP PARALLEL default(none)
C$OMP& shared(i,nneighipart,nearpt,iflagtree)
C$OMP& private(j,jpart)
C$OMP DO SCHEDULE(static)
            DO j = 1, nneighipart
               jpart = nearpt(j,i)
               iflagtree(jpart) = .TRUE.
            END DO
C$OMP END DO
C$OMP END PARALLEL
         END DO
      ENDIF

      DO i = 1, nlstacc
         ipart = listacc(i)
         IF (ipart.LT.1 .OR. ipart.GT.npart) THEN
            WRITE (iprint, *) 'ERROR - revtree listacc',ipart,nlstacc,i
            CALL quit
         ENDIF
         iflagtree(ipart) = .TRUE.
      END DO

C$OMP PARALLEL default(none)
C$OMP& shared(natom,npart,xyzmh,imfac,isibdaupar)
C$OMP& shared(nactatom,listmap,iflagtree,ipar,iphase)
C$OMP& private(j,ipart,iparent)

C$OMP DO SCHEDULE(static)
      DO j = 1, nactatom
         ipart = listmap(j)
         IF (iflagtree(ipart)) THEN
            iflagtree(ipart) = .FALSE.
            iparent = isibdaupar(3,ipart)
C$OMP CRITICAL (flagtree)
            iflagtree(iparent) = .TRUE.
C$OMP END CRITICAL (flagtree)
            IF (ipart.GT.npart .OR. iphase(ipart).EQ.-1 .OR.
     &               (iphase(ipart).GE.1 .AND. iptintree.EQ.1)) THEN
               imfac(ipart) = 0
            ELSE
               imfac(ipart) = 1
            ENDIF
         ENDIF
      END DO
C$OMP END DO
C$OMP END PARALLEL

      DO ilevel = 1, nlevel

         IF ((level(ilevel + 1)-level(ilevel)).GT.16) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(ilevel,level,natom,isibdaupar)
C$OMP& shared(qrad,xyzmh,imfac)
C$OMP& shared(iflagtree,ipar,nlevel)
C$OMP& private(new,l,ll,fl,fll,emred,difx,dify,difz,rr,iparent)
C$OMP& private(pmassl,pmassll)

C$OMP DO SCHEDULE(static)
         DO new = level(ilevel), level(ilevel + 1) - 1

            IF (iflagtree(new)) THEN

            iflagtree(new) = .FALSE.
            IF (ilevel.LT.nlevel) THEN
               iparent = isibdaupar(3,new)
C$OMP CRITICAL (flagtree)
               iflagtree(iparent) = .TRUE.
C$OMP END CRITICAL (flagtree)
            ENDIF

            l = isibdaupar(2,new)
            ll = isibdaupar(1,l)

            IF (l.LE.natom) THEN
               pmassl = imfac(l)*xyzmh(4,l)
            ELSE
               pmassl = xyzmh(4,l)
            ENDIF
            IF (ll.LE.natom) THEN
               pmassll = imfac(ll)*xyzmh(4,ll)
            ELSE
               pmassll = xyzmh(4,ll)
            ENDIF
            xyzmh(4,new) = pmassl + pmassll

            IF (xyzmh(4,new).NE.0) THEN
               fl = pmassl/xyzmh(4,new)
               fll = pmassll/xyzmh(4,new)
            ELSE
               fl = 0.5
               fll = 0.5
            ENDIF
            emred = fl*fll*xyzmh(4,new)
            difx = xyzmh(1,ll) - xyzmh(1,l)
            dify = xyzmh(2,ll) - xyzmh(2,l)
            difz = xyzmh(3,ll) - xyzmh(3,l)
            IF (fl.GT.fll) THEN
               xyzmh(1,new) = xyzmh(1,l) + fll*difx
               xyzmh(2,new) = xyzmh(2,l) + fll*dify
               xyzmh(3,new) = xyzmh(3,l) + fll*difz
            ELSE
               xyzmh(1,new) = xyzmh(1,ll) - fl*difx
               xyzmh(2,new) = xyzmh(2,ll) - fl*dify
               xyzmh(3,new) = xyzmh(3,ll) - fl*difz
            ENDIF
c
c--Find radius
c
            rr = SQRT(difx**2 + dify**2 + difz**2) + tiny
            qrad(1,new) = MAX(fll*rr + qrad(1,l), fl*rr + qrad(1,ll))
            xyzmh(5,new) = MAX(xyzmh(5,l), xyzmh(5,ll))
c
c--Find quadrupole moments
c
            qrad(2,new) = (emred*difx)*difx
            qrad(5,new) = (emred*difx)*dify
            qrad(7,new) = (emred*difx)*difz
            qrad(3,new) = (emred*dify)*dify
            qrad(6,new) = (emred*dify)*difz
            qrad(4,new) = (emred*difz)*difz
            IF (l.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,l)
               qrad(3,new) = qrad(3,new) + qrad(3,l)
               qrad(4,new) = qrad(4,new) + qrad(4,l)
               qrad(5,new) = qrad(5,new) + qrad(5,l)
               qrad(6,new) = qrad(6,new) + qrad(6,l)
               qrad(7,new) = qrad(7,new) + qrad(7,l)
            ENDIF
            IF (ll.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,ll)
               qrad(3,new) = qrad(3,new) + qrad(3,ll)
               qrad(4,new) = qrad(4,new) + qrad(4,ll)
               qrad(5,new) = qrad(5,new) + qrad(5,ll)
               qrad(6,new) = qrad(6,new) + qrad(6,ll)
               qrad(7,new) = qrad(7,new) + qrad(7,ll)
            ENDIF

            ENDIF
         END DO
C$OMP END DO
C$OMP END PARALLEL

         ELSE

         DO new = level(ilevel), level(ilevel + 1) - 1

            IF (iflagtree(new)) THEN

            iflagtree(new) = .FALSE.
            IF (ilevel.LT.nlevel) THEN
               iparent = isibdaupar(3,new)
               iflagtree(iparent) = .TRUE.
            ENDIF

            l = isibdaupar(2,new)
            ll = isibdaupar(1,l)

            IF (l.LE.natom) THEN
               pmassl = imfac(l)*xyzmh(4,l)
            ELSE
               pmassl = xyzmh(4,l)
            ENDIF
            IF (ll.LE.natom) THEN
               pmassll = imfac(ll)*xyzmh(4,ll)
            ELSE
               pmassll = xyzmh(4,ll)
            ENDIF
            xyzmh(4,new) = pmassl + pmassll

            IF (xyzmh(4,new).NE.0) THEN
               fl = pmassl/xyzmh(4,new)
               fll = pmassll/xyzmh(4,new)
            ELSE
               fl = 0.5
               fll = 0.5
            ENDIF
            emred = fl*fll*xyzmh(4,new)
            difx = xyzmh(1,ll) - xyzmh(1,l)
            dify = xyzmh(2,ll) - xyzmh(2,l)
            difz = xyzmh(3,ll) - xyzmh(3,l)
            IF (fl.GT.fll) THEN
               xyzmh(1,new) = xyzmh(1,l) + fll*difx
               xyzmh(2,new) = xyzmh(2,l) + fll*dify
               xyzmh(3,new) = xyzmh(3,l) + fll*difz
            ELSE
               xyzmh(1,new) = xyzmh(1,ll) - fl*difx
               xyzmh(2,new) = xyzmh(2,ll) - fl*dify
               xyzmh(3,new) = xyzmh(3,ll) - fl*difz
            ENDIF
c
c--Find radius
c
            rr = SQRT(difx**2 + dify**2 + difz**2) + tiny
            qrad(1,new) = MAX(fll*rr + qrad(1,l), fl*rr + qrad(1,ll))
            xyzmh(5,new) = MAX(xyzmh(5,l), xyzmh(5,ll))
c
c--Find quadrupole moments
c
            qrad(2,new) = (emred*difx)*difx
            qrad(5,new) = (emred*difx)*dify
            qrad(7,new) = (emred*difx)*difz
            qrad(3,new) = (emred*dify)*dify
            qrad(6,new) = (emred*dify)*difz
            qrad(4,new) = (emred*difz)*difz
            IF (l.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,l)
               qrad(3,new) = qrad(3,new) + qrad(3,l)
               qrad(4,new) = qrad(4,new) + qrad(4,l)
               qrad(5,new) = qrad(5,new) + qrad(5,l)
               qrad(6,new) = qrad(6,new) + qrad(6,l)
               qrad(7,new) = qrad(7,new) + qrad(7,l)
            ENDIF
            IF (ll.GT.natom) THEN
               qrad(2,new) = qrad(2,new) + qrad(2,ll)
               qrad(3,new) = qrad(3,new) + qrad(3,ll)
               qrad(4,new) = qrad(4,new) + qrad(4,ll)
               qrad(5,new) = qrad(5,new) + qrad(5,ll)
               qrad(6,new) = qrad(6,new) + qrad(6,ll)
               qrad(7,new) = qrad(7,new) + qrad(7,ll)
            ENDIF

            ENDIF
         END DO

         ENDIF

      END DO

      ENDIF

      RETURN
      END
