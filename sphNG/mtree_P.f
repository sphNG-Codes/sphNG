      SUBROUTINE mtree(nnatom, npart, xyzmh)
c************************************************************
c                                                           *
c  Subroutine by W. Press (11/21/86).  Given the positions  *
c     and masses of NATOM points, this subroutine          *
c     constructs a nearest neighbor tree (including         *
c     quadrupole moments) for use by TREEFORCE, filling EM, *
c     R, Q, ISIB, IPAR, IDAU, and NROOT, the pointer to the *
c     root of the tree (which is also the largest value     *
c     used in the range 1..MMAX)                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      PARAMETER (bigno=1.E30, ncbrt=330, hashfac=1.0)

      DIMENSION xyzmh(5,mmax2)

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/interstellar'

      DIMENSION xyzmap(3,ncbrt)
      LOGICAL*1 iused(idim)
      INTEGER*8 ibin
      CHARACTER*7 where

      DATA where/'mtree'/

      INTEGER levelcount(nmaxlevel),lc(nmaxlevel)
      INTEGER inodeloc(idim)
c
c--Allow for tracing flow
c
c      IF (itrace.EQ.'all') WRITE (iprint, 99001)
c99001 FORMAT (' entry subroutine mtree')

      third = 1./3.

      IF (idustRT.GT.0) THEN
         IF (nh2frac.NE.mmax) CALL error(where,5)
      ENDIF

      natom = nnatom
      nactatom = 0
      DO j = 1, natom
         IF (iphase(j).EQ.0 .OR. iphase(j).GE.10 .OR. 
     &        (iphase(j).GE.1 .AND. iphase(j).LT.10 .AND. 
     &        iptintree.GT.0)) THEN
            nactatom = nactatom + 1
            listmap(nactatom) = j
         ENDIF
      END DO

      nactive = nactatom
      nactold = 0
      newend = natom
      niteration = 0
      nlevel = 0
      levelnum(1:nnatom) = 0
      levelcount(:) = 0

C$OMP PARALLEL default(none)
C$OMP& shared(natom,npart,list,isibdaupar)
C$OMP& shared(xyzmh,qrad,listmap,iflagtree)
C$OMP& private(j,ipart)
C$OMP& shared(nactatom,nactive,newend,where,iused)
C$OMP& shared(nay,iphase)
C$OMP& shared(ihash,nhash,key,nglob)
C$OMP& private(l,new,n,ll,fl,fll,emred,difx,dify,difz,rr)
C$OMP& shared(icbrt,icbrt1,xyzmap,level,imfac,levelnum)
C$OMP& shared(next1,next2,next3,third,nactold,iprint,newold,niteration)
C$OMP& shared(h2frac)
C$OMP& private(m,np,mm,iz,iy,ix,nc,ddd,xnp,ynp,znp)
C$OMP& private(lllx,llux,llly,lluy,lllz,lluz,llx,lux,lly,luy,llz,luz)
C$OMP& private(jx,icjx,jy,icjy,jz,k,d,dmin,nwal)
C$OMP& private(ibin,i,itemp,nj,ibin1,pmassl,pmassll)
C$OMP& reduction(+:levelcount)
C$OMP& reduction(MAX:nlevel)

C$OMP DO SCHEDULE(runtime)
      DO j = 1, nactatom
         ipart = listmap(j)
         list(j) = ipart
         IF (ipart.GT.npart .OR. 
     &        (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10 
     &        .AND. iptintree.EQ.1)) THEN
            imfac(ipart) = 0
         ELSE
            imfac(ipart) = 1
         ENDIF
         IF (idustRT.GT.0) THEN
            IF (ipart.GT.npart .OR. iphase(ipart).NE.0) THEN
               h2frac(ipart) = 0.
            ENDIF
         ENDIF
         isibdaupar(1,ipart) = 0
         isibdaupar(2,ipart) = 0
c 1=mono, 2=xx,3=yy,4=zz,5=xy,6=yz,7=zx
         qrad(1,ipart) = 0.
         qrad(2,ipart) = 0.
         qrad(3,ipart) = 0.
         qrad(4,ipart) = 0.
         qrad(5,ipart) = 0.
         qrad(6,ipart) = 0.
         qrad(7,ipart) = 0.
         iflagtree(ipart) = .FALSE.
      END DO
C$OMP END DO
c
c--Return here to fill in each new level of hierarchy:
c     find nearest neighbors:
c
 50   CONTINUE

c************************************************************
c                                                           *
c  Used to be a subroutine called naybor(nactive)           *
c                                                           *
c************************************************************

C$OMP SINGLE
      icbrt = hashfac*float(nactive)**third
      icbrt1 = icbrt - 1
      nhash = icbrt**3
      IF (nhash.GT.nactive) CALL error(where, 1)
      IF (icbrt.GT.ncbrt) CALL error(where, 2)
C$OMP END SINGLE

C$OMP SECTIONS
      CALL indexx3(nactive, list, xyzmh, 1, next1)
C$OMP SECTION
      CALL indexx3(nactive, list, xyzmh, 2, next2)
C$OMP SECTION
      CALL indexx3(nactive, list, xyzmh, 3, next3)
C$OMP END SECTIONS
c
c--Bin the points by X and fill XMAP
c
C$OMP DO SCHEDULE(runtime)
      DO ibin = 1, icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) xyzmap(1,ibin)=0.5*(xyzmh(1,list(next1(i)))
     &                              + xyzmh(1,list(next1(i + 1))))
         itemp = icbrt*(ibin - 1)
         DO j = ((ibin-1)*nactive)/icbrt+1, i
            key(next1(j)) = itemp
         END DO
      END DO
C$OMP END DO
c
c--Bin the points by Y and fill YMAP
c
C$OMP DO SCHEDULE(runtime)
      DO ibin = 1, icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) xyzmap(2,ibin)=0.5*(xyzmh(2,list(next2(i))) 
     &                             + xyzmh(2,list(next2(i + 1))))
         ibin1 = ibin - 1
         DO j = ((ibin-1)*nactive)/icbrt+1, i
            nj = next2(j)
            key(nj) = icbrt*(key(nj) + ibin1)
         END DO
      END DO
C$OMP END DO
c
c--Bin the points by Z and fill ZMAP
c
C$OMP DO SCHEDULE(runtime)
      DO ibin = 1, icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) xyzmap(3,ibin)=0.5*(xyzmh(3,list(next3(i))) 
     &                             + xyzmh(3,list(next3(i + 1))))
         DO j = ((ibin-1)*nactive)/icbrt+1, i
            nj = next3(j)
            key(nj) = key(nj) + ibin
         END DO
      END DO
C$OMP END DO
c
c--Now fill the head-of-list table and the linked list
c
C$OMP DO SCHEDULE(runtime)
      DO m = 1, nhash
         ihash(m) = 0
      END DO
C$OMP END DO

C$OMP SINGLE
      DO j = 1, nactive
         m = key(j)
c--next1 is used for something different here to save memory !
         next1(j) = ihash(m)
         ihash(m) = j
      END DO
C$OMP END SINGLE
c
c--Now loop over the particles to find nearest neighbors
c
C$OMP DO SCHEDULE(runtime)
      DO np = 1, nactive
c
c--Find its cell
c
         mm = key(np) - 1
         iz = MOD(mm, icbrt)
         mm = mm/icbrt
         iy = MOD(mm, icbrt)
         ix = mm/icbrt
c
c--Set closest distance so far and coordinates
c
         nc = 0
         ddd = bigno
         xnp = xyzmh(1,list(np))
         ynp = xyzmh(2,list(np))
         znp = xyzmh(3,list(np))
c
c--Initialize the limits of the search loop
c
         lllx = ix
         llux = ix
         llly = iy
         lluy = iy
         lllz = iz
         lluz = iz
         llx = ix
         lux = ix
         lly = iy
         luy = iy
         llz = iz
         luz = iz
c
c--Loop for finding closest particle in the volume of new cells
c
 100     DO jx = llx, lux
            icjx = icbrt*jx
            DO jy = lly, luy
               icjy = icbrt*(jy + icjx) + 1
               DO jz = llz, luz
                  m = icjy + jz
                  k = ihash(m)
                  IF (k.NE.0) THEN
 200                 IF (k.NE.np) THEN
                        d = (xyzmh(1,list(k)) -xnp)**2 +
     &                       (xyzmh(2,list(k)) -ynp)
     &                      **2 + (xyzmh(3,list(k)) - znp)**2

                        IF (d.LT.ddd) THEN
                           ddd = d
                           nc = k
                        ENDIF
                     ENDIF
                     k = next1(k)
                     IF (k.NE.0) GOTO 200
                  ENDIF
               END DO
            END DO
         END DO
c
c--We must now find the closest wall, if any
c
         dmin = bigno
         nwal = 0
         IF (lllx.GT.0) THEN
            d = xnp - xyzmap(1,lllx)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 1
               ENDIF
            ENDIF
         ENDIF
         IF (llux.LT.icbrt1) THEN
            d = xyzmap(1,llux + 1) - xnp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 2
               ENDIF
            ENDIF
         ENDIF
         IF (llly.GT.0) THEN
            d = ynp - xyzmap(2,llly)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 3
               ENDIF
            ENDIF
         ENDIF
         IF (lluy.LT.icbrt1) THEN
            d = xyzmap(2,lluy + 1) - ynp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 4
               ENDIF
            ENDIF
         ENDIF
         IF (lllz.GT.0) THEN
            d = znp - xyzmap(3,lllz)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 5
               ENDIF
            ENDIF
         ENDIF
         IF (lluz.LT.icbrt1) THEN
            d = xyzmap(3,lluz + 1) - znp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 6
               ENDIF
            ENDIF
         ENDIF
c
c--Reset the search volume to its augmented value
c
         llx = lllx
         lux = llux
         lly = llly
         luy = lluy
         llz = lllz
         luz = lluz
c
c--Augment it and set the next search according to which wall is closest
c
         IF (nwal.NE.0) THEN
            IF (nwal.EQ.1) THEN
               lllx = lllx - 1
               llx = lllx
               lux = lllx
            ELSEIF (nwal.EQ.2) THEN
               llux = llux + 1
               llx = llux
               lux = llux
            ELSEIF (nwal.EQ.3) THEN
               llly = llly - 1
               lly = llly
               luy = llly
            ELSEIF (nwal.EQ.4) THEN
               lluy = lluy + 1
               lly = lluy
               luy = lluy
            ELSEIF (nwal.EQ.5) THEN
               lllz = lllz - 1
               llz = lllz
               luz = lllz
            ELSEIF (nwal.EQ.6) THEN
               lluz = lluz + 1
               llz = lluz
               luz = lluz
            ENDIF
            GOTO 100
         ENDIF
         nay(np) = nc
      END DO
C$OMP END DO

c************************************************************
c                                                           *
c  End of old subroutine called naybor(nactive)             *
c                                                           *
c************************************************************

C$OMP SINGLE
      IF (nactive.EQ.nactold) THEN
         WRITE (iprint,99100) nactive
99100    FORMAT (' MTREE IS IN AN INFINITE LOOP ',I8)
         DO k = 1, nactive
            l = list(k)
            WRITE (iprint,99101) l, xyzmh(1,l), xyzmh(2,l), xyzmh(3,l)
99101       FORMAT (I8, 1X, 1F14.7, 1X, 1F14.7, 1X, 1F14.7)
         END DO   
         CALL quit(1)
      ENDIF
C$OMP END SINGLE
c
c--Find new nodes:
c
C$OMP DO SCHEDULE(runtime)
      DO j = 1, nactive
         l = list(j)
         new = newend + j
         IF (new.GT.mmax-1) CALL error(where, 3)
         iused(j) = .FALSE.
         n = nay(j)

         IF (isibdaupar(1,l).EQ.0 .AND. n.GT.j) THEN
            IF (nay(n).EQ.j) THEN
               iused(j) = .TRUE.
               ll = list(n)
               isibdaupar(1,ll) = l
               isibdaupar(1,l) = ll
               isibdaupar(3,ll) = new
               isibdaupar(3,l) = new
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
c               levelnum(new) = niteration + 1
               nodelevel = MAX(levelnum(l),levelnum(ll))+1
               levelnum(new) = nodelevel
               levelcount(nodelevel) = levelcount(nodelevel) + 1
               nlevel = MAX(nlevel,nodelevel)
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
               xyzmh(1,new) = xyzmh(1,l) + fll*difx
               xyzmh(2,new) = xyzmh(2,l) + fll*dify
               xyzmh(3,new) = xyzmh(3,l) + fll*difz
c
c--Find radius and maximum h
c
               IF (ileveltreeacc .EQ. 1) THEN
                  rr = SQRT(difx**2 + dify**2 + difz**2)
                  qrad(1,new) = MAX(fll*rr+qrad(1,l), fl*rr+qrad(1,ll))
               ELSE
                  qrad(1,new) = tree_node_radius(isibdaupar,xyzmh,qrad,
     &                 new,l,ll,natom)
               ENDIF
               IF (qrad(1,new).EQ.0.) THEN
                  print *,'ERROR - mtree node has zero size: ',new,l,ll,
     &                 natom,pmassl,pmassll,fl,fll,rr,
     &                 qrad(1,l),qrad(1,ll),xyzmh(:,l),xyzmh(:,ll)
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
               qrad(4,new) = emred*(difz*difz)
               IF (l.GT.natom) THEN
                  qrad(2,new) = qrad(2,new) + qrad(2,l)
                  qrad(5,new) = qrad(5,new) + qrad(5,l)
                  qrad(7,new) = qrad(7,new) + qrad(7,l)
                  qrad(3,new) = qrad(3,new) + qrad(3,l)
                  qrad(6,new) = qrad(6,new) + qrad(6,l)
                  qrad(4,new) = qrad(4,new) + qrad(4,l)
               ENDIF
               IF (ll.GT.natom) THEN
                  qrad(2,new) = qrad(2,new) + qrad(2,ll)
                  qrad(5,new) = qrad(5,new) + qrad(5,ll)
                  qrad(7,new) = qrad(7,new) + qrad(7,ll)
                  qrad(3,new) = qrad(3,new) + qrad(3,ll)
                  qrad(6,new) = qrad(6,new) + qrad(6,ll)
                  qrad(4,new) = qrad(4,new) + qrad(4,ll)
               ENDIF
               IF (idustRT.GT.0) THEN
                  h2frac(new) = h2frac(l)*fl+h2frac(ll)*fll
               ENDIF
               isibdaupar(1,new) = 0
               isibdaupar(2,new) = l
               iflagtree(new) = .FALSE.
            ENDIF
         ENDIF
      END DO
C$OMP END DO
c
c--Compactify list:
c
C$OMP SINGLE
      newold = newend
      nglob = 1
      DO j = 1, nactive
         IF (iused(j)) THEN
            new = newold + j
            newend = newend + 1
            CALL move_node(new,newend,xyzmh)
         ELSEIF (isibdaupar(1,list(j)).EQ.0) THEN
            list(nglob) = list(j)
            nglob = nglob + 1
         ENDIF`
      END DO
C$OMP END SINGLE

C$OMP DO SCHEDULE(runtime)
      DO j = newold + 1, newend
         list(nglob+j-(newold + 1)) = j
      END DO
C$OMP END DO

C$OMP SINGLE
      nglob = nglob + (newend - newold)

      niteration = niteration + 1
      IF (nlevel+1.GT.nmaxlevel) CALL error(where,4)
c      level(niteration) = newold + 1
      nactold = nactive
      IF (nglob - 1.LE.idim) THEN
         nactive = nglob - 1
      ELSE
         WRITE (iprint,*) 'ERROR - nactive>idim'
         CALL quit(1)
      ENDIF
C$OMP END SINGLE

      IF (nactive.GT.1) GOTO 50

C$OMP END PARALLEL

c      niteration = niteration + 1
c      IF (nlevel.GT.nmaxlevel) CALL error(where,4)
c      level(niteration) = newend + 1
c      niteration = niteration - 1

      nroot = newend
c
c--Find the number of tree nodes on each level.
c
      ncpg = levelcount(nlevel)
      lc(1) = 0
      DO i=1,nlevel-1
         lc(i+1) = lc(i) + levelcount(i)
         level(i) = lc(i) + 1
         ncpg = ncpg + levelcount(i)
      ENDDO
      level(nlevel) = lc(nlevel) + 1
      level(nlevel+1) = level(nlevel) + levelcount(nlevel)
c
c--Re-order tree so that nodes on the same level are contiguous.
c
      istart = nnatom + 1
      itemp = mmax
      icount = 0
      DO ilevel=1,nlevel
         ilevelcount = 0
         n_in_level = levelcount(ilevel)
         igap = 0
         iwrong = 0
         iend = istart + n_in_level - 1
c
c--Find indices of nodes not in ilevel that are in assigned slots.
c
         DO inode=istart,iend
            IF (levelnum(inode) .NE. ilevel) THEN
               inodeloc(idim-igap) = inode
               igap = igap + 1
            ENDIF
            IF (ilevelcount .EQ. n_in_level) EXIT
         ENDDO
c
c--Find indices of nodes on ilevel that are beyond assigned slots.
c
         ilevelcount = n_in_level - igap
         DO inode=iend+1,nroot
            IF (levelnum(inode) .EQ. ilevel) THEN
               ilevelcount = ilevelcount + 1
               iwrong = iwrong + 1
               inodeloc(iwrong) = inode
            ENDIF
            IF (ilevelcount .EQ. n_in_level) EXIT
         ENDDO
         IF (igap .NE. iwrong) THEN
            WRITE(iprint,*)"igap.NE.iwrong",igap,iwrong
            CALL quit(1)
         ENDIF
c
c--Switch nodes into their assigned slots.
c
         DO i=1,igap
            ngap = inodeloc(idim-i+1)
            nwrong = inodeloc(i)
            CALL move_node(ngap,itemp)
            CALL move_node(nwrong,ngap)
            CALL move_node(itemp,nwrong)
            iwrong = iwrong - 1
            icount = icount + 1
         ENDDO
         istart = istart + n_in_level
      ENDDO
      WRITE(iprint,*)"number of switched nodes",icount

      RETURN
      END


c
c--This routine moves a tree node from location ifrom to ito.
c
      SUBROUTINE move_node(ifrom,ito,xyzmh)

      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/interstellar'

      INTEGER, INTENT(IN) :: ifrom,ito
      LOGICAL, INTENT(IN) :: iswitch
      REAL, INTENT(IN) :: xyzmh(5,mmax2)

      levelnum(ito) = levelnum(ifrom)
      xyzmh(1,ito) = xyzmh(1,ifrom)
      xyzmh(2,ito) = xyzmh(2,ifrom)
      xyzmh(3,ito) = xyzmh(3,ifrom)
      xyzmh(4,ito) = xyzmh(4,ifrom)
      xyzmh(5,ito) = xyzmh(5,ifrom)
      qrad(1,ito) = qrad(1,ifrom)
      qrad(2,ito) = qrad(2,ifrom)
      qrad(3,ito) = qrad(3,ifrom)
      qrad(4,ito) = qrad(4,ifrom)
      qrad(5,ito) = qrad(5,ifrom)
      qrad(6,ito) = qrad(6,ifrom)
      qrad(7,ito) = qrad(7,ifrom)
      IF (idustRT.GT.0) h2frac(ito) =  h2frac(ifrom)
      isibdaupar(1,ito) = isibdaupar(1,ifrom)
      isibdaupar(2,ito) = isibdaupar(2,ifrom)
      isibdaupar(3,isibdaupar(2,ito)) = ito
      isibdaupar(3,isibdaupar(1,isibdaupar(2,ito))) = ito
      iflagtree(ito) = iflagtree(ifrom)

      END SUBROUTINE move_node


c
c--Walks the tree to a defined depth (ileveltreeacc). This allows the
c     node radius to be refined to be less of an overestimate.
c
      FUNCTION tree_node_radius(isibdaupar,xyzmh,qrad,new,l,ll,natom)

      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INTEGER, DIMENSION(3,mmax), INTENT(IN) :: isibdaupar
      REAL*4, DIMENSION(7,mmax), INTENT(IN) :: qrad
      REAL, DIMENSION(5,mmax), INTENT(IN) :: xyzmh
      INTEGER, INTENT(IN) :: new,l,ll,natom

      INTEGER :: node,nd,istack,new_dau
      REAL :: distmax,dist2,tree_node_radius

      INTEGER, DIMENSION(istacksizemrev), SAVE :: nstack,ndepth
      REAL, DIMENSION(3), SAVE :: xyznew
C$OMP THREADPRIVATE(nstack,ndepth,xyznew)

      distmax = 0.
      xyznew = xyzmh(1:3,new)
      nstack(1) = l
      nstack(2) = ll
      ndepth(1) = 1
      ndepth(2) = 1
      istack = 2
      DO WHILE (istack.NE.0)
         node = nstack(istack)
         nd = ndepth(istack)
         istack = istack - 1
         IF (node.GT.natom .AND. nd.LT.ileveltreeacc) THEN
            new_dau = isibdaupar(2,node)
            istack = istack + 1
            nstack(istack) = new_dau
            ndepth(istack) = nd + 1
            istack = istack + 1
            nstack(istack) = isibdaupar(1,new_dau)
            ndepth(istack) = nd + 1
         ELSE
            dist2 = SUM((xyznew - xyzmh(1:3,node))**2)
            IF (node.GT.natom) THEN
               distmax = MAX(distmax,SQRT(dist2)+qrad(1,node))
            ELSE
               distmax = MAX(distmax,SQRT(dist2))
            ENDIF
         ENDIF
      ENDDO
      tree_node_radius = distmax

      END FUNCTION tree_node_radius
