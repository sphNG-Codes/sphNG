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
      de_cs_star    =  1E30

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,iphase,xyzmh,radkernel)
C$OMP& reduction(MIN:domain_extent,de_cs_gas,de_cs_star)
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
               IF (iphasei .EQ. 0) THEN
                  de_cs_gas(k) = MIN(de_cs_gas(k),xyzki-rcs)
                  de_cs_gas(l) = MIN(de_cs_gas(l),rcs-xyzki)
               ELSEIF (iphasei .GE. 10) THEN
                  de_cs_star(k) = MIN(de_cs_star(k),xyzki-rcs)
                  de_cs_star(l) = MIN(de_cs_star(l),rcs-xyzki)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C$OMP END PARALLEL DO
      
      END SUBROUTINE reset_domain_extent
