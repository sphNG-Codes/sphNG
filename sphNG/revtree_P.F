      SUBROUTINE revtree (nnatom, npart, xyzmh)
c************************************************************
c                                                           *
c  Subroutine by W. Press.  Revises the tree structure.     *
c     Assumes that only the positions have been altered     *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax)

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/perform'

      DIMENSION listpar(idim),listparnew(idim)

#ifdef _OPENMP
      INTEGER*8 forcei_lock,revtree_lock
      COMMON /locks / forcei_lock(idim),revtree_lock
#endif

      natom = nnatom

c
c--REVISE ENTIRE TREE (STANDARD REVTREE)
c
c      IF (nlst.GT.nnatom/5000 .OR. ipartialrevtree.EQ.0) THEN
      IF (nlst.GT.nnatom/10000 .OR. ipartialrevtree.EQ.0) THEN
c      IF (nlst.GT.100 .OR. ipartialrevtree.EQ.0) THEN

      IF (itiming) CALL getused(revtreeptemp1)

C$OMP PARALLEL default(none)
C$OMP& shared(natom,npart,imfac)
C$OMP& shared(nactatom,listmap,iphase,iflagtree)
C$OMP& private(j,ipart)

C$OMP DO SCHEDULE(static)
      DO j = 1, nactatom
         ipart = listmap(j)
         iflagtree(ipart) = .FALSE.
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

      IF (itiming) THEN
         CALL getused(revtreeptemp2)
         revtreep1 = revtreep1 + (revtreeptemp2 - revtreeptemp1)
      ENDIF
c
c--ONLY REVISE PART OF THE TREE - PARTICLES WHOSE POSITIONS HAVE BEEN
c     UPDATED IN STEP AND ANY ACCRETED OR KILLED PARTICLES FROM THE 
c     LAST STEP
c
      ELSE

      IF (itiming) CALL getused(revtreeptemp1)

      iupdatenode = 0

      IF (nlst.GT.nptmass) iupdatenode = 1




      GOTO 888




      IF (nptmass.GT.10) THEN
         
C$OMP PARALLEL default(none)
C$OMP& shared(nptmass,nptlist,nearpt,iflagtree)
C$OMP& private(i,j,jpart)
C$OMP& reduction(+:iupdatenode)
C$OMP DO SCHEDULE(static)
         DO i = 1, nptmass
            DO j = 1, nptlist(i)
               jpart = nearpt(j,i)
C$OMP CRITICAL (flagtree)
               iflagtree(jpart) = .TRUE.
C$OMP END CRITICAL (flagtree)
               iupdatenode = iupdatenode + 1
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
C$OMP& reduction(+:iupdatenode)
C$OMP DO SCHEDULE(static)
            DO j = 1, nneighipart
               jpart = nearpt(j,i)
               iflagtree(jpart) = .TRUE.
               iupdatenode = iupdatenode + 1
            END DO
C$OMP END DO
C$OMP END PARALLEL
         END DO
      ENDIF

      IF (itiming) THEN
         CALL getused(revtreeptemp2)
         revtreep3 = revtreep3 + (revtreeptemp2 - revtreeptemp1)
      ENDIF

 888  CONTINUE

      iupdatenode = iupdatenode + nlstacc

      nlistpar = 0

c      print *,'Revtree nlistpar ',nlistpar

      IF (nlstacc.GT.0) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(nlstacc,listacc,npart,iflagtree,listpar,isibdaupar)
C$OMP& shared(iphase,imfac)
C$OMP& private(i,ipart)

C$OMP DO SCHEDULE(static)
         DO i = 1, nlstacc
            ipart = listacc(i)
            IF (ipart.LT.1 .OR. ipart.GT.npart) THEN
               WRITE (*, *) 'ERROR - revtree listacc',ipart,nlstacc,i
               CALL quit
            ENDIF
            listpar(i) = isibdaupar(3,ipart)
            IF (iphase(ipart).EQ.-1) imfac(ipart) = 0
         END DO
C$OMP END DO
C$OMP END PARALLEL

         nlistpar = nlstacc

c         print *,'Revtree nlistpar ',nlistpar,nlstacc

         IF (itiming) THEN
            CALL getused(revtreeptemp2)
            revtreep4 = revtreep4 + (revtreeptemp2 - revtreeptemp1)
         ENDIF
      ENDIF

      IF (iupdatenode.GT.0) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(natom,npart,xyzmh,imfac,isibdaupar,nlst,llist)
C$OMP& shared(revtree_lock)
C$OMP& shared(nactatom,listmap,iflagtree,ipar,iphase,nlistpar)
C$OMP& private(j,ipart,iparent)

C$OMP DO SCHEDULE(static)
         DO j = 1, nlst
            ipart = llist(j)
            iflagtree(ipart) = .FALSE.
            iparent = isibdaupar(3,ipart)

            IF (iparent.NE.0) THEN
#ifdef _OPENMP
               CALL OMP_SET_LOCK(revtree_lock)
#endif
               IF (.NOT.iflagtree(iparent)) THEN
                  iflagtree(iparent) = .TRUE.
                  nlistpar = nlistpar + 1
                  listpar(nlistpar) = iparent
               ENDIF
#ifdef _OPENMP
               CALL OMP_UNSET_LOCK(revtree_lock)
#endif
            ENDIF

            IF (ipart.GT.npart .OR. iphase(ipart).EQ.-1 .OR.
     &           (iphase(ipart).GE.1 .AND. iptintree.EQ.1)) THEN
               imfac(ipart) = 0
            ELSE
               imfac(ipart) = 1
            ENDIF
         END DO
C$OMP END DO
C$OMP END PARALLEL

c         nlistpar = nlistpar + nlst

c         print *,'Revtree nlistpar pre ',nlistpar, nlst

         IF (itiming) THEN
            CALL getused(revtreeptemp2)
            revtreep5 = revtreep5 + (revtreeptemp2 - revtreeptemp1)
         ENDIF

 456     CONTINUE

c         print *,'Revtree nlistpar doing ',nlistpar

         nlistparnew = 0

C$OMP PARALLEL default(none)
C$OMP& shared(ilevel,level,natom,isibdaupar,revtree_lock)
C$OMP& shared(qrad,xyzmh,imfac,listpar,listparnew)
C$OMP& shared(iflagtree,ipar,nlevel,nlistpar,nlistparnew)
C$OMP& private(i,new,l,ll,fl,fll,emred,difx,dify,difz,rr,iparent)
C$OMP& private(pmassl,pmassll,oldmass,oldh,oldqrad)

C$OMP DO SCHEDULE(static)
         DO i = 1, nlistpar

c            print *,i

            new = listpar(i)
            iflagtree(new) = .FALSE.

            iparent = isibdaupar(3,new)
            IF (iparent.NE.0) THEN
#ifdef _OPENMP
               CALL OMP_SET_LOCK(revtree_lock)
#endif
               IF (.NOT.iflagtree(iparent)) THEN
                  iflagtree(iparent) = .TRUE.
                  nlistparnew = nlistparnew + 1
                  listparnew(nlistparnew) = iparent
               ENDIF
#ifdef _OPENMP
               CALL OMP_UNSET_LOCK(revtree_lock)
#endif
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

            oldmass = xyzmh(4,new)
            oldh = xyzmh(5,new)

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
            oldqrad = qrad(1,new)
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

c            IF (ilevel.LT.nlevel .OR. 
c     &           ABS(oldmass-xyzmh(4,new))/oldmass.GT.1.0E-4 .OR.
c     &           ABS(oldh-xyzmh(5,new))/oldh.GT.1.0E-0 .OR.
c     &           ABS(oldqrad-qrad(1,new))/oldqrad.GT.1.0E-0) THEN
c               iflagtree(iparent) = .TRUE.
c            ENDIF

c            ENDIF
         END DO
C$OMP END DO
C$OMP DO SCHEDULE(static)
         DO i = 1, nlistparnew
            listpar(i) = listparnew(i)
         END DO
C$OMP END DO
C$OMP END PARALLEL
         nlistpar = nlistparnew

         IF (nlistpar.GT.0) GOTO 456

         IF (itiming) THEN
            CALL getused(revtreeptemp2)
            revtreep2 = revtreep2 + (revtreeptemp2 - revtreeptemp1)
         ENDIF

      ENDIF

      ENDIF

      RETURN
      END
