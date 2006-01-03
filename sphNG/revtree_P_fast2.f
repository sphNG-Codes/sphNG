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
      INCLUDE 'COMMONS/timeextra'

      natom = nnatom

c
c--REVISE ENTIRE TREE (STANDARD REVTREE)
c
c      IF (nlst.GT.nnatom/5000 .OR. (.NOT. ipartialrevtree)) THEN
      IF (nlst-nptmass.GT.10000 .OR. itbinupdate.GE.nbinmax-1 .OR. 
     &     (.NOT. ipartialrevtree)) THEN

      IF (itiming) CALL getused(revtreeptemp1)

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

         IF ((level(ilevel + 1)-level(ilevel)).GT.2) THEN

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
c--ONLY REVISE PART OF THE TREE - CURRENT PARTICLES AND THEIR NEIGHBOURS
c     AND ANY ACCRETED OR KILLED PARTICLES FROM THE LAST STEP
c
      ELSE

      IF (itiming) CALL getused(revtreeptemp1)

c      IF (itiming) THEN
c         CALL getused(revtreeptemp2)
c         revtreep3 = revtreep3 + (revtreeptemp2 - revtreeptemp1)
c      ENDIF

c      IF (itiming) THEN
c         CALL getused(revtreeptemp2)
c         revtreep4 = revtreep4 + (revtreeptemp2 - revtreeptemp1)
c      ENDIF

      IF (nlst.GT.nptmass .OR. iptintree.EQ.2) THEN

         DO i = 1, nlstacc
            ipart = listacc(i)
            IF (ipart.LT.1 .OR. ipart.GT.npart) THEN
               WRITE (iprint, *) 'ERROR - revtree listacc',
     &                           ipart,nlstacc,i
               CALL quit
            ENDIF
            imfac(ipart) = 0
            iparent = isibdaupar(3,ipart)
            iflagtree(iparent) = .TRUE.
            numberparents = numberparents + 1
            listparents(numberparents) = iparent
         END DO
c
c--Note: could sort list of parents here!
c
         IF (itiming) THEN
            CALL getused(revtreeptemp2)
            revtreep5 = revtreep5 + (revtreeptemp2 - revtreeptemp1)
         ENDIF

         numberstart = 1
 400     numnextlevel = numberparents

         IF (numberparents-numberstart.GT.2) THEN

C$OMP PARALLEL default(none)
C$OMP& shared(natom,isibdaupar)
C$OMP& shared(qrad,xyzmh,imfac)
C$OMP& shared(iflagtree,ipar,numberstart,nroot,iprint)
C$OMP& shared(numberparents,listparents,numnextlevel)
C$OMP& shared(nlst,nptmass,itbinupdate,nlstacc)
C$OMP& private(new,l,ll,fl,fll,emred,difx,dify,difz,rr,iparent,i)
C$OMP& private(pmassl,pmassll)
C$OMP& private(qrad1old,xold,yold,zold)

C$OMP DO SCHEDULE(static)
            DO i = numberstart, numberparents
               new = listparents(i)

               IF (.NOT.iflagtree(new)) THEN
                  WRITE (iprint, *) 'ERROR - revtreeX1',new,i
                  WRITE (iprint, *) nlst,nptmass,iptintree,itbinupdate,
     &                 nbinmax,nlstacc,numberparents,numberstart
                  CALL quit
               ENDIF

               iflagtree(new) = .FALSE.

               qrad1old = qrad(1,new)
               xold = xyzmh(1,new)
               yold = xyzmh(2,new)
               zold = xyzmh(3,new)

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

               IF (new.NE.nroot) THEN
                  IF ((xyzmh(1,new)-xold)**2 + (xyzmh(2,new)-yold)**2 + 
     &                 (xyzmh(3,new)-zold)**2.GT.1.0E-6*qrad1old) THEN
                     iparent = isibdaupar(3,new)
C$OMP CRITICAL(parentlist5)
                     IF (.NOT.iflagtree(iparent)) THEN
                        iflagtree(iparent) = .TRUE.
                        numnextlevel = numnextlevel + 1
                        listparents(numnextlevel) = iparent
                     ENDIF
C$OMP END CRITICAL(parentlist5)                     
                  ENDIF
               ENDIF
               
            END DO
C$OMP END DO
C$OMP END PARALLEL
c
c--Else don't bother to start up parallel threads
c
         ELSE
            DO i = numberstart, numberparents
               new = listparents(i)

               IF (.NOT.iflagtree(new)) THEN
                  WRITE (iprint, *) 'ERROR - revtreeX2',new,i
                  CALL quit
               ENDIF

               iflagtree(new) = .FALSE.

               qrad1old = qrad(1,new)
               xold = xyzmh(1,new)
               yold = xyzmh(2,new)
               zold = xyzmh(3,new)

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

               IF (new.NE.nroot) THEN
                  IF ((xyzmh(1,new)-xold)**2 + (xyzmh(2,new)-yold)**2 + 
     &                 (xyzmh(3,new)-zold)**2.GT.1.0E-6*qrad1old) THEN
                     iparent = isibdaupar(3,new)
                     IF (.NOT.iflagtree(iparent)) THEN
                        iflagtree(iparent) = .TRUE.
                        numnextlevel = numnextlevel + 1
                        listparents(numnextlevel) = iparent
                     ENDIF
                  ENDIF
               ENDIF
               
            END DO
         ENDIF

         numberstart = numberparents + 1
         numberparents = numnextlevel
         IF (numberparents.GT.numberstart) GOTO 400

      ENDIF

      IF (itiming) THEN
         CALL getused(revtreeptemp2)
         revtreep2 = revtreep2 + (revtreeptemp2 - revtreeptemp1)
      ENDIF

      ENDIF

      RETURN
      END
