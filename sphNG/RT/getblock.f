      SUBROUTINE getblock(m,rrx,rry,rrz,numneigh,neighb,xyzmh)
c
c*********************************************************************
c                                                                    *
c     Return list of GAS particles whose smoothing sphere (2h)       *
c     covers the line joining particle m and position (rrx,rry,rrz)  *
c     MRB: 4.1.2007                                                  *
c                                                                    *
c*********************************************************************
c
      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax2),neighb(idim)

      INCLUDE '../COMMONS/treecom_P'
      INCLUDE '../COMMONS/nlim'
      INCLUDE '../COMMONS/logun'
      INCLUDE '../COMMONS/phase'
      INCLUDE '../COMMONS/sort'
      INCLUDE '../COMMONS/kerne'

      DIMENSION itest(100)

      PARAMETER (istacksize=200)
      DIMENSION nstack(istacksize)
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

      dxm = xyzmh(1,m) - rrx
      dym = xyzmh(2,m) - rry
      dzm = xyzmh(3,m) - rrz
      rrm = dxm**2 + dym**2 + dzm**2

      if (m.EQ.9999999) print *,'Examine ',dxm,dym,dzm,numneigh

      if (m.EQ.9999999) then
         itest(1) = isibdaupar(3,25920)
         dxi = xyzmh(1,itest(1)) - rrx
         dyi = xyzmh(2,itest(1)) - rry
         dzi = xyzmh(3,itest(1)) - rrz
         rr = dxi*dxi + dyi*dyi + dzi*dzi

         print *,'    ',itest(1),SQRT(rr),qrad(1,itest(1))
         DO i = 2, 20
            itest(i) = isibdaupar(3,itest(i-1))
            dxi = xyzmh(1,itest(i)) - rrx
            dyi = xyzmh(2,itest(i)) - rry
            dzi = xyzmh(3,itest(i)) - rrz
            rr = dxi*dxi + dyi*dyi + dzi*dzi

            print *,'    ',itest(i),SQRT(rr),qrad(1,itest(i)),
     &           isibdaupar(2,itest(i)),isibdaupar(1,isibdaupar
     &           (2,itest(i)))

         END DO
      end if
      if (m.EQ.9999999) print *,'Want ',25920,isibdaupar(3,25920),
     &     isibdaupar(3,isibdaupar(3,25920))

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
c
c--Should give the closest possible particle in the node to (rrx,rry,rrz)
c     If it is <=0 then this means there could be particles on top of the
c     location and/or on both sides of the position
c
               rr_closest = (MAX(0.0,SQRT(rr) - qrad(1,n)))**2

               if (m.EQ.9999999)  then
          print *,'node ',n,rr_closest,rrm,rr,dxi,dyi,dzi,numneigh,
     &              qrad(1,n),xyzmh(5,n)

               idau = isibdaupar(2,n)
               isibdau = isibdaupar(1,idau)

               DO itestit = 1, 20
                  IF (itest(itestit).EQ.idau) print *,' HERE****idau ',
     &                 idau
               IF (itest(itestit).EQ.isibdau) print *,' HERE****idau ',
     &                 isibdau
               END DO
               print *,'components ',idau,isibdau
               print *,'component1 ',xyzmh(1,idau),
     & xyzmh(2,idau),xyzmh(3,idau)
               print *,'component2 ',xyzmh(1,isibdau),
     & xyzmh(2,isibdau),xyzmh(3,isibdau)

               endif

c               qcut2 = (qrad(1,n) + rcut)**2
               if (m.EQ.9999999) print *,'TEST ',rr_closest,rrm,
     &              rr_closest.LT.rrm, SQRT(rr),qrad(1,n)

               IF (istack+2.GT.istacksize) THEN
                  WRITE(iprint,*)'ERROR: stack too small in getblock!',
     &                 istack+2,istacksize
                  CALL quit
               ENDIF

               IF (rr_closest.LT.rrm) THEN
                  if (m.EQ.9999999) print *,'INN'
                  j = isibdaupar(2,n)
                  istack = istack + 1
                  nstack(istack) = j
                  istack = istack + 1
                  nstack(istack) = isibdaupar(1,j)
                  if (m.EQ.9999999) print *,'ADD ',j,nstack(istack)
               ENDIF
            ELSE
               IF (iphase(n).EQ.0) THEN
c
c--Must be closer in order to block
c
                  IF (rr.LT.rrm) THEN

                     if (m.EQ.9999999) print *,dxi,dyi,dzi,numneigh
                     costheta_mn = dxm*dxi + dym*dyi + dzm*dzi
                  if (m.EQ.9999999) print *,'costheta ',costheta_mn
                     IF (costheta_mn.GE.0.0) THEN
c
c     Gives costheta of angle between i-(rrx,rry,rrz)-j
c
                        costheta_mn = costheta_mn/SQRT(rrm*rr)
c
c     Dot product of unit vector to particle i and potential blocker j gives
c     cos(theta).  Angle of blocking is cos(theta)=A/H where H is |r(j)| and
c     A is the distance to the tangent with 2h circle which is also given
c     by H^2 = A^2 + O^2 since the triangle is a right angle triangle.
c     Then, since theta1<theta2 if cos(theta1)<cos(theta2) don't need to
c     use arccos to find actual angles to find whether particle i is blocked.
c
                        hn2 = (xyzmh(5,n)*radkernel)**2
                        tempval = (rr - hn2)/rr
                        IF (tempval.LT.0.0) tempval = 0.0
                        costheta_block = SQRT(tempval)

                        if (m.EQ.9999999)
     &   print *,'costheta B MN ',costheta_block,costheta_mn,SQRT(hn2)
                        IF (costheta_block.LT.costheta_mn) THEN
c
c--If potentially blocking particle overlaps with
c     particle m, then it is *not* considered to block particle m
c
                           IF ((dxm-dxi)**2+(dym-dyi)**2+
     &                          (dzm-dzi)**2.GE.hn2) THEN
c                  IF (rr.LT.rcut2) THEN
                           numneigh = numneigh + 1
                           RETURN
c                        IF (numneigh.LE.nneighmax) THEN
c                           neighb(numneigh) = n
c                        ELSE
c                           WRITE (iprint,*) 'ERROR: nneighmax exceeded'
c                           STOP
c                        ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ELSE
                     if (m.EQ.9999999) print *,'Not closer'
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
