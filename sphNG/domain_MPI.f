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
                  de_cs_gas(l) = MIN(de_cs_gas(l),-xyzki-rcs)
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

      
c--find minimum distance squared between two axis aligned cuboids a & b
c     ext_a/b(1:3): minimum xyz coords of cuboid
c     ext_a/b(4:6): maximum xyz coords of cuboid
      PURE FUNCTION cuboid_dist2(ext_a,ext_b)

      IMPLICIT NONE

      REAL, INTENT(IN) :: ext_a(6),ext_b(6)
      REAL :: xdiff,ydiff,zdiff,cuboid_dist2

      xdiff = MAX(ext_a(1)-ext_b(4), 0.0, ext_b(1)-ext_a(4))
      ydiff = MAX(ext_a(2)-ext_b(5), 0.0, ext_b(2)-ext_a(5))
      zdiff = MAX(ext_a(3)-ext_b(6), 0.0, ext_b(3)-ext_a(6))
      cuboid_dist2 = xdiff**2 + ydiff**2 + zdiff**2

      END FUNCTION cuboid_dist2
