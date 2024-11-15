      SUBROUTINE revtree (nnatom, npart, xyzmh)
c************************************************************
c                                                           *
c  Subroutine by W. Press.  Revises the tree structure.     *
c     Assumes that only the positions have been altered     *
c                                                           *
c     Partial tree revision completed by TB December 2024   *
c     This was started by MRB pre-2006                      *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax2)

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
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/timeextra'

      INTEGER iflag

      natom = nnatom

c
c--REVISE ENTIRE TREE (STANDARD REVTREE)
c
      IF (ifullrevision) THEN

      IF (itiming) CALL getused(revtreeptemp1)
      nfullrev = nfullrev + 1

C$OMP PARALLEL default(none)
C$OMP& shared(natom,npart,imfac,h2frac)
C$OMP& shared(nactatom,listmap,iphase)
C$OMP& private(j,ipart)

C$OMP DO SCHEDULE(static)
      DO j = 1, nactatom
         ipart = listmap(j)
         IF (ipart.GT.npart .OR. iphase(ipart).EQ.-1 .OR.
     &        (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10 .AND. 
     &        iptintree.EQ.1)) THEN
            imfac(ipart) = 0
         ELSE
            imfac(ipart) = 1
         ENDIF
         IF (idustRT.GT.0) THEN
            IF (ipart.GT.npart .OR. iphase(ipart).NE.0) THEN
               h2frac(ipart) = 0.
            ENDIF
         ENDIF
      END DO
C$OMP END DO
C$OMP END PARALLEL

      DO ilevel = 1, nlevel
         nnode_on_level = level(ilevel + 1)-level(ilevel)

C$OMP PARALLEL IF(nnode_on_level > 2) default(none)
C$OMP& shared(ilevel,level,natom,isibdaupar)
C$OMP& shared(qrad,xyzmh,imfac,h2frac)
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
            IF (ileveltreeacc .EQ. 1) THEN
               rr = SQRT(difx**2 + dify**2 + difz**2) + tiny
               qrad(1,new) = MAX(fll*rr+qrad(1,l), fl*rr+qrad(1,ll))
            ELSE
               qrad(1,new) = tree_node_radius(isibdaupar,xyzmh,qrad,
     &              new,l,ll,natom)
            ENDIF
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
            IF (idustRT.GT.0) THEN
               h2frac(new) = h2frac(l)*fl+h2frac(ll)*fll
            ENDIF
         END DO
C$OMP END DO
C$OMP END PARALLEL
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

      IF (nlst.LE.nptmass .AND. iptintree.NE.2 .AND. nlstacc.LE.0) THEN
         RETURN
      ENDIF
      IF (itiming) CALL getused(revtreeptemp1)
      npartialrev =  npartialrev + 1

      DO i = 1, nlstacc
         ipart = listacc(i)
         IF (ipart.LT.1 .OR. ipart.GT.npart) THEN
            WRITE (iprint, *) 'ERROR - revtree listacc',
     &           ipart,nlstacc,i
            CALL quit(1)
         ENDIF
         imfac(ipart) = 0
         iparent = isibdaupar(3,ipart)
         IF (.NOT.iflagtree(iparent)) THEN
            iflagtree(iparent) = .TRUE.
            numberparents = numberparents + 1
            listparents(numberparents) = iparent
         ENDIF
      END DO
      IF (itiming) THEN
         CALL getused(revtreeptemp2)
         revtreep5 = revtreep5 + (revtreeptemp2 - revtreeptemp1)
      ENDIF
c
c--Level zero are leaves, not nodes, so start from nodes (level 1)
c
      numberremaining = numberparents
      nlevelupdate = 0
         
      DO WHILE (numberremaining .GT. 0)
         nlevelupdate = nlevelupdate + 1
         numlistpar = numberremaining
         numberremaining = 0
         numthislevel = 0
c
c--Find nodes to revise on this level
c
         DO i=1,numlistpar
            ipart = listparents(i)
            IF (levelnum(ipart) .EQ. nlevelupdate) THEN               
               numthislevel = numthislevel + 1
               list(numthislevel) = ipart
            ELSE
               numberremaining = numberremaining + 1
               listparents(numberremaining) = ipart
            ENDIF
         ENDDO

C$OMP PARALLEL IF (numthislevel > 2) default(none)
C$OMP& shared(natom,isibdaupar)
C$OMP& shared(qrad,xyzmh,imfac,h2frac)
C$OMP& shared(iflagtree,numberparents,nroot,iprint)
C$OMP& shared(numberremaining,listparents,numthislevel,list)
C$OMP& shared(nlst,nptmass,itbinupdate,nlstacc)
C$OMP& shared(levelnum,nlevelupdate)
C$OMP& private(new,l,ll,fl,fll,emred,difx,dify,difz,rr,iparent,i,j)
C$OMP& private(pmassl,pmassll,iflag,numberadd)
C$OMP& private(qrad1old,xold,yold,zold)

C$OMP DO SCHEDULE(static)
         DO i = 1, numthislevel
            new = list(i)

            IF (.NOT.iflagtree(new)) THEN
               WRITE (iprint, *) 'ERROR - revtreeX1',new,i
               WRITE (iprint, *) nlst,nptmass,iptintree,itbinupdate,
     &              nbinmax,nlstacc,numberremaining,numberparents
               CALL quit(1)
            ENDIF

            iflagtree(new) = .FALSE.

            l = isibdaupar(2,new)
            ll = isibdaupar(1,l)

            qrad1old = qrad(1,new)

            xold = xyzmh(1,new)
            yold = xyzmh(2,new)
            zold = xyzmh(3,new)

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
            IF (ileveltreeacc .EQ. 1) THEN
               rr = SQRT(difx**2 + dify**2 + difz**2)
               qrad(1,new) = MAX(fll*rr+qrad(1,l), fl*rr+qrad(1,ll))
            ELSE
               qrad(1,new) = tree_node_radius(isibdaupar,xyzmh,qrad,
     &              new,l,ll,natom)
            ENDIF
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
            IF (idustRT.GT.0) THEN
               h2frac(new) = h2frac(l)*fl+h2frac(ll)*fll
            ENDIF

            IF (new.NE.nroot) THEN
               iparent = isibdaupar(3,new)

               IF (nlevelupdate.GE.levelnum(iparent)) THEN
                  WRITE (*,*) 'ERROR revtree level ',nlevelupdate,
     &                 levelnum(iparent),new,iparent,levelnum(new)
               ENDIF
C$OMP ATOMIC CAPTURE
               iflag = iflagtree(iparent)
               iflagtree(iparent) = .TRUE.
C$OMP END ATOMIC
               IF (.NOT.iflag) THEN
C$OMP ATOMIC CAPTURE
                  numberremaining = numberremaining + 1
                  numberadd = numberremaining
C$OMP END ATOMIC
                  IF (numberadd.GT.idim) THEN
                     WRITE (*,*) 'parentlist5 ',
     &                    numberremaining,idim
                     CALL quit(1)
                  ENDIF
                  listparents(numberadd) = iparent
               ENDIF
            ENDIF
         END DO
C$OMP END DO
C$OMP END PARALLEL
      ENDDO

      IF (itiming) THEN
         CALL getused(revtreeptemp2)
         revtreep2 = revtreep2 + (revtreeptemp2 - revtreeptemp1)
      ENDIF

      ENDIF

      RETURN
      END
