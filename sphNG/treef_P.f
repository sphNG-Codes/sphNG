      SUBROUTINE treef(m,npart,xyzmh,acc,igrav,fsx,fsy,fsz,epot)
c
c************************************************************
c                                                           *
c  Subroutine by W. Press (11/21/86). Given particle number *
c     M and cutoff radius R, returns a list NEARL(1..NLIST) *
c     of particles within a distance RCUT, and returns      *
c     three components of gravitational force on M due go   *
c     all particles outside RCUT.  This routine presumes    *
c     that subroutine MTREE has previously been called.     *
c  Updated 15/3/93 for individual timesteps as in Julio's   *
c     code.                                                 *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax)

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/sort'

      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/call'

      INTEGER istack, nlistga, nlistgn

      DIMENSION nstack(100), listga(idim), listgn(idim)
      SAVE nstack, listga, listgn
C$OMP THREADPRIVATE(nstack,listga,listgn)

      CHARACTER*7 where

      DATA where/'treef'/
c
c--Allow for tracing flow
c
c      IF (itrace.EQ.'all') WRITE(iprint,250)
  250 FORMAT(' entry subroutine treef')
c
c--Accreted particles not in tree
c
      IF (iphase(m).EQ.-1) CALL error(where, 1)
c
c--Point masses no longer in tree - gravity and neighbours done in gforspt.f
c
      IF (iphase(m).GE.1 .AND. iptintree.EQ.0) CALL error(where, 2)

      IF (iphase(m).GE.1) THEN
         DO i = 1, nptmass
            IF (listpm(i).EQ.m) THEN
               iptn = i
               nptlist(iptn) = 0
               GOTO 50
            ENDIF
         END DO
      ENDIF

 50   nneigh(m) = 0 
      gravmag = SQRT(gravxyzstore(1,m)**2 + gravxyzstore(2,m)**2 + 
     &     gravxyzstore(3,m)**2)
      accpar = acc**2
      fsx = 0.
      fsy = 0.
      fsz = 0.
      epot = 0.
      rrx = xyzmh(1,m)
      rry = xyzmh(2,m)
      rrz = xyzmh(3,m)
      nlistga = 0
      nlistgn = 0
      istack = 0
      mpar = m
      hm = xyzmh(5,m)

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , M1',m
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , M1',m
c            CALL quit
c         ENDIF
c      END DO

 100  node = isibdaupar(1,mpar)
      IF (node.NE.0) THEN
         istack = istack + 1
         nstack(istack) = node
 150     IF (istack.NE.0) THEN
            n = nstack(istack)
            istack = istack - 1
            dxi = xyzmh(1,n) - rrx
            dyi = xyzmh(2,n) - rry
            dzi = xyzmh(3,n) - rrz
            rr = dxi*dxi + dyi*dyi + dzi*dzi
c
c--Decide whether to open up the node:
c     new if structure
c
            IF (n.GT.natom) THEN
               IF (nlmax.EQ.1) THEN
                  rcutmn = MAX(hm,xyzmh(5,n))*radkernel
               ELSE
                  rcutmn = (hm + xyzmh(5,n))*radkernel/2.0
               ENDIF
               qradn = qrad(1,n)
               qcut2 = (qradn + rcutmn)**2
               qrr = qradn**2
               IF (n.LE.natom) THEN
                  pmassn = imfac(n)*xyzmh(4,n)
               ELSE
                  pmassn = xyzmh(4,n)
               ENDIF
               IF (rr.LT.qcut2) THEN
                  j = isibdaupar(2,n)
                  istack = istack + 1
                  nstack(istack) = j
                  istack = istack + 1
                  nstack(istack) = isibdaupar(1,j)
               ELSEIF (
     &                 icall.EQ.1 .AND. 
     &         (igrav.NE.0 .OR. (iphase(m).GE.1 .AND. iptintree.EQ.1))
     &                 .AND. accpar*rr.LT.qrr) THEN
                  j = isibdaupar(2,n)
                  istack = istack + 1
                  nstack(istack) = j
                  istack = istack + 1
                  nstack(istack) = isibdaupar(1,j)
               ELSEIF (icall.NE.1 .AND. 
     &         (igrav.NE.0 .OR. (iphase(m).GE.1 .AND. iptintree.EQ.1))
     &                 .AND. 
c--Quadrupole
c     &                 0.005*gravmag*rr**3.LT.
     &                 0.0035*gravmag*rr**3.LT.
     &                 pmassn*qrr**2) THEN
c--Monopole
c     &                 0.0001*gravmag*rr**2.LT.
c     &                 pmassn*qrr) THEN
                  j = isibdaupar(2,n)
                  istack = istack + 1
                  nstack(istack) = j
                  istack = istack + 1
                  nstack(istack) = isibdaupar(1,j)
               ELSE
                  nlistgn = nlistgn + 1
                  IF (nlistgn.GT.idim) THEN
                     WRITE (iprint,*) 'ERROR - nlistgn.GT.idim'
                     CALL quit
                  ENDIF
                  listgn(nlistgn) = n
               ENDIF
            ELSE
c
c--This line was altered for individual timesteps
c
               hn = xyzmh(5,n)
               IF (nlmax.EQ.1) THEN
                  rcut2 = (MAX(hm,hn)*radkernel)**2
               ELSE
                  rcut2 = ((hm + hn)*radkernel/2.0)**2
               ENDIF
               IF (rr.LT.rcut2 .AND. iphase(n).EQ.0 .AND.
     &                                         iphase(m).EQ.0) THEN
                  nneigh(m) = nneigh(m) + 1

                  IF (nlmax.EQ.1) THEN
                     IF (nneigh(m).LE.nneighmax) THEN
                        neighlist(nneigh(m)) = n
                     ELSE
                        WRITE (iprint,*) 'nneigh exceeds ',nneighmax
                        CALL quit
                     ENDIF
                  ELSE
                     IF (nneigh(m).LT.nlmax) THEN
                        neighb(nneigh(m),m) = n
                     ELSEIF (nneigh(m).GT.nneighmax) THEN
                        WRITE (iprint,*) 'nneigh exceeds INT*2'
                        CALL quit
                     ELSE
                        IF (nneigh(m).EQ.nlmax) THEN
c                        write (*,*) 'over ',m,n,neighb(nlmax,m),
c     &                       numoverflow
                           IF (neighb(nlmax,m).LT.0) THEN
                              numoverflowlocal = -neighb(nlmax,m)
                        neighover(nneigh(m)-nlmax+1,numoverflowlocal)=n
                           ELSE
C$OMP CRITICAL (neighlistoverflow)
                              IF (numoverflow.LT.noverflow) THEN
                                 numoverflow = numoverflow + 1
                                 numoverflowlocal = numoverflow
                                 neighb(nneigh(m),m)= -numoverflowlocal
                        neighover(nneigh(m)-nlmax+1,numoverflowlocal)=n
                              ELSE
                                 WRITE (iprint,*) 'noverflow EXCEEDED',
     &                                nneigh(m), nlmax, m, iorig(m),
     &                                numoverflow
                                 CALL quit
                              ENDIF
C$OMP END CRITICAL (neighlistoverflow)
                           ENDIF
                        ELSEIF (nneigh(m)-nlmax+1.LE.nlovermax) THEN
                        neighover(nneigh(m)-nlmax+1,numoverflowlocal)=n
                        ELSE
                           WRITE (iprint,*) 'nlovermax EXCEEDED',
     &                          nneigh(m), nlmax, m, iorig(m), 
     &                          numoverflowlocal
                           CALL quit
                        ENDIF
                     ENDIF
                  ENDIF

               ELSE
                  nlistga = nlistga + 1
                  IF (nlistga.GT.idim) THEN
                     WRITE (iprint,*) 'ERROR - nlistga.GT.idim'
                     CALL quit
                  ENDIF
                  listga(nlistga) = n
               ENDIF

               IF (iptintree.GT.0 .AND. iphase(m).GE.1 .AND.
     &                               iphase(n).EQ.0) THEN
                  hmean2 = hm*hm
ccc                  hmean2 = rcut2
                  IF (rr.LT.hmean2) THEN
                     nptlist(iptn) = nptlist(iptn) + 1
                     IF (nptlist(iptn).GT.iptneigh) THEN
                        WRITE (iprint,*) 'ERROR - iptneigh '
                        CALL quit
                     ELSE
                        nearpt(nptlist(iptn),iptn) = n
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
            GOTO 150
         ENDIF
         mpar = isibdaupar(3,mpar)
         GOTO 100
      ENDIF

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , M2',m
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , M2',m
c         ENDIF
c      END DO
c
c--Calculate gravity forces if needed
c
      IF (igrav.NE.0 .OR. (iphase(m).GE.1 .AND. iptintree.EQ.1)) THEN
c
c--Get force from the atomic list
c
         CALL gforsa(m, nlistga, listga, xyzmh, fsx, fsy, fsz, epot)
c
c--Get force from the node list
c
         CALL gforsn(m, nlistgn, listgn, xyzmh, fsx, fsy, fsz, epot)
      ENDIF
      epot = -epot

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , M3',m
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , M3',m
c         ENDIF
c      END DO

      RETURN
      END


      SUBROUTINE getneighi(m,rrx,rry,rrz,rcut,numneigh,neighb,xyzmh)
c
c*********************************************************************
c                                                                    *
c     Return list of neighbours for a single particle within         *
c     fixed search radius rcut                                       *
c     DJP:10.11.2005 and MRB:15.12.2005                              *
c                                                                    *
c*********************************************************************
c
      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax),neighb(idim)

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/sort'

      DIMENSION nstack(100)
      SAVE nstack
C$OMP THREADPRIVATE(nstack)

c*************************************
c--find neighbours within rcut       *
c*************************************
      numneigh = 0

      IF (iphase(m).EQ.-1 .OR. iphase(m).GE.1 .AND. iptintree.EQ.0) THEN
         WRITE (iprint,*) 'ERROR: getneighi with iphase ',iphase(m)
         CALL quit
      ENDIF

      istack = 0
      mpar = m
      rcut2 = rcut*rcut
 100  node = isibdaupar(1,mpar)
      IF (node.NE.0) THEN
         istack = istack + 1
         nstack(istack) = node
 150     IF (istack.NE.0) THEN
            n = nstack(istack)
            istack = istack - 1

            dxi = xyzmh(1,n) - rrx
            dyi = xyzmh(2,n) - rry
            dzi = xyzmh(3,n) - rrz
            rr = dxi*dxi + dyi*dyi + dzi*dzi
c
c--Decide whether to open up the node:
c     new if structure
c
            IF (n.GT.natom) THEN
               qcut2 = (qrad(1,n) + rcut)**2
               IF (rr.LT.qcut2) THEN
                  j = isibdaupar(2,n)
                  istack = istack + 1
                  nstack(istack) = j
                  istack = istack + 1
                  nstack(istack) = isibdaupar(1,j)
               ENDIF
            ELSE
               IF (rr.LT.rcut2) THEN
                  numneigh = numneigh + 1
                  IF (numneigh.LE.nneighmax) THEN
                     neighb(numneigh) = n
                  ELSE
                     WRITE (iprint,*) 'ERROR: nneighmax exceeded'
                     STOP
                  ENDIF
               ENDIF
            ENDIF
            GOTO 150
         ENDIF
         mpar = isibdaupar(3,mpar)
         GOTO 100
      ENDIF

      RETURN
      END

