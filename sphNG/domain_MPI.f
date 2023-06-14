c--set domain extent coordinates
      SUBROUTINE reset_domain_extent

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/mpidomains'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/kerne'

      INTEGER i,k,l,iphasei
      REAL rcs,xyzki
      
      domain_extent =  1E30
      de_cs_gas     =  1E30

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,iphase,xyzmh,radkernel)
C$OMP& reduction(MIN:domain_extent,de_cs_gas)
C$OMP& private(i,k,l,rcs,xyzki,iphasei)
      DO i=1,npart
         iphasei = iphase(i)
         IF (iphasei .GE. 0) THEN
            rcs = radkernel * xyzmh(5,i)
            DO k=1,3
               l = k + 3
               xyzki = xyzmh(k,i)
               domain_extent(k) = MIN(domain_extent(k),xyzki)
               domain_extent(l) = MIN(domain_extent(l),-xyzki)
               IF (iphasei .EQ. 0 .OR. iphasei .GE. 10) THEN
                  de_cs_gas(k) = MIN(de_cs_gas(k),xyzki-rcs)
                  de_cs_gas(l) = MIN(de_cs_gas(l),rcs-xyzki)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C$OMP END PARALLEL DO
      
      END SUBROUTINE reset_domain_extent


c--sorting routine for sorting by index high to low
c     O(N^2) with very low overhead,
c     fine for short arrays < 150 values
c     very slow for long arrays
      SUBROUTINE insert_sort_index(N,ARR,indx)
      INTEGER ARR(N)
      INTEGER indx(N),N

      indx(1) = 1
      DO j=2, N
         l = ARR(j)
         indx(j) = j
         DO i=j-1,1,-1
            IF (ARR(indx(i)) .GE. l) goto 10
            indx(i+1) = indx(i)
         ENDDO
         i = 0
 10      indx(i+1) = j
      ENDDO
      RETURN
      END


c--as above but sorts array itself (not by index,) low to high    
      SUBROUTINE insert_sort(N,ARR)
      INTEGER ARR(N)
      INTEGER N

      DO j=2, N
         l = ARR(j)
         i = j - 1
         DO i=j-1,1,-1
            IF (ARR(i) .LE. l) goto 10
            ARR(i+1) = ARR(i)
         ENDDO
         i = 0
 10      ARR(i+1) = l
      ENDDO

      RETURN
      END SUBROUTINE insert_sort

      
c--find minimum distance squared between two axis aligned cuboids a & b
c     ext_a/b(1:3): minimum xyz coords of cuboid
c     ext_a/b(4:6): maximum xyz coords of cuboid (inverted signs)
      PURE FUNCTION cuboid_dist2(ext_a,ext_b)

      IMPLICIT NONE

      REAL, INTENT(IN) :: ext_a(6),ext_b(6)
      REAL :: xdiff,ydiff,zdiff,cuboid_dist2

      xdiff = MAX(ext_a(1)+ext_b(4), 0.0, ext_b(1)+ext_a(4))
      ydiff = MAX(ext_a(2)+ext_b(5), 0.0, ext_b(2)+ext_a(5))
      zdiff = MAX(ext_a(3)+ext_b(6), 0.0, ext_b(3)+ext_a(6))
      cuboid_dist2 = xdiff**2 + ydiff**2 + zdiff**2

      END FUNCTION cuboid_dist2


c--as cuboid_dist2 but one cuboid has 0 extent i.e. particle
      PURE FUNCTION point_cuboid_d2(xyzp,ext)

      IMPLICIT NONE

      REAL, INTENT(IN) :: xyzp(3),ext(6)
      REAL :: xdiff,ydiff,zdiff,point_cuboid_d2

      xdiff = MAX(xyzp(1)+ext(4), 0.0, ext(1)-xyzp(1))
      ydiff = MAX(xyzp(2)+ext(5), 0.0, ext(2)-xyzp(2))
      zdiff = MAX(xyzp(3)+ext(6), 0.0, ext(3)-xyzp(3))
      point_cuboid_d2 = xdiff**2 + ydiff**2 + zdiff**2

      END FUNCTION point_cuboid_d2


c--walks the tree using breadth-first traversal
c     finds all tree nodes/leaves to send to a neighbouring MPI rank
c     this walk leads to partially sorted array of nodes/leaves
c     probably a way to do this traversal with istacksize < 2*idim
      SUBROUTINE grav_only_treef(de_ip,xyzmh,acc,nplist,nblist,listg)

      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

      REAL xyzmh(5,mmax2)

      INCLUDE 'COMMONS/mpidomains'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/logun'

      INTEGER n,j,jj,k,itemp,istack,icurrent,nplist,nblist
      REAL de_ip(6)
      REAL sep2,sep,accpar,acc,rr,point_cuboid_d2,qrr
      LOGICAL addnode

      INTEGER nstack(2*idim), listg(ngrav_nodes)
      SAVE nstack
C$OMP THREADPRIVATE(nstack)

      accpar = acc**2
c--testing gravity accuarcy parameter at lower tolerance
c     acc = 0.5 for local gravity calculations
      accpar = 0.7**2

      nplist = 0
      nblist = 0
      icurrent = 0
      istack = 0
      n = nroot

      DO WHILE (icurrent .LE. istack)
         rr = point_cuboid_d2(xyzmh(1:3,n),de_ip)
         IF (n.GT.natom) THEN
            qrr = qrad(1,n)**2
            IF (accpar*rr.LT.qrr) THEN
               IF (istack+2.GT.idim*2) THEN
                  WRITE(iprint,*)'ERROR: stack too small in grav_only!',
     &                 istack+2,idim*2
                  CALL quit(1)
               ENDIF
               j = isibdaupar(2,n)
               nstack(istack+2) = j
               nstack(istack+1) = isibdaupar(1,j)
               istack = istack + 2
            ELSE
               nblist = nblist + 1
               listg(nblist) = n - 1
c               WRITE(iprint,*) "nn,lev",n,levelnum(n)
            ENDIF
         ELSE
            nplist = nplist + 1
            nblist = nblist + 1
c            WRITE(iprint,*) "np,lev",n,levelnum(n)
c--if any particles need to be sent give and send all particles on rank
            RETURN
            listg(nblist) = n - 1
         ENDIF
         IF (nblist+1 .GT. ngrav_nodes) THEN
            WRITE (iprint,*)'ERROR - nlistg.GT.ngrav_nodes part',
     &           nblist,ngrav_nodes
            CALL quit(0)
         ENDIF
         icurrent = icurrent + 1
         n = nstack(icurrent)
c         WRITE(iprint,*) "end",icurrent,n
      ENDDO

c--reverse order of original array so that partial sort is increasing
      DO j=1,nblist/2
         k = nblist-j+1
         itemp = listg(j)
         listg(j) = listg(k)
         listg(k) = itemp
      ENDDO

      END SUBROUTINE grav_only_treef
