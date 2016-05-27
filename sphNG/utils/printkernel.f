      PROGRAM printkernel
      IMPLICIT NONE
      
      INCLUDE '../COMMONS/kerne'
      INCLUDE '../COMMONS/table'
      INTEGER i
      REAL vi, phi
      
      CALL ktable
      
      WRITE(*,*) '# sphNG kernel table, radkern = ',radkernel, 
     &           ' cnormk = ',cnormk
      DO i=0,itable
         vi = sqrt(i*dvtable)
         phi = fpoten(i)
         IF (vi.GT.part2kernel) THEN
            phi = phi + part2potenkernel/vi
         ELSEIF (vi.GT.part1kernel) THEN
            phi = phi + part1potenkernel/vi
         ENDIF
         WRITE(*,"(6(es12.4,1x))") vi,
     &        cnormk*wij(i),cnormk*grwij(i),fmass(i),phi,
     &        dphidh(i)
      ENDDO
      
      END
