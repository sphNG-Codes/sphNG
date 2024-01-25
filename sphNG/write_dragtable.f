      SUBROUTINE write_dragtable(idragfile)

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/xforce'
      
      INTEGER i, j, k
      INTEGER idragfile
      REAL hmass, hillmass

      hmass = hillmass(planetmass(1), 0, xmass, 0)

      WRITE (idragfile, ERR=100) idragresr, idragresp, idragrest
      WRITE (idragfile, ERR=100) dragtrmin, dragtrmax, dragtdr, 
     &     dragtpmin, dragtpmax, dragtdp, dragttmin, dragttmax, dragtdt
      WRITE (idragfile, ERR=100) planetmass(1), hmass

      DO i = 1, idragresr
         DO j = 1, idragresp
            DO k = 1, idragrest
               WRITE (idragfile, ERR=100) dragenergy(i,j,k)
            END DO
         END DO
      END DO

      CLOSE (idragfile)

      RETURN
      
 100  PRINT *, 'An error has occurred writing drag file'
      CALL quit(0) 

      END SUBROUTINE write_dragtable

