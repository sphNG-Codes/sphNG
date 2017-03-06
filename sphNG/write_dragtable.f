      subroutine write_dragtable(idragfile)

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/xforce'
      
      INTEGER i, j, k
      INTEGER idragfile
      REAL hmass, hillmass

      hmass = hillmass(planetmass, 0, xmass, 0)

      WRITE (idragfile, ERR=100) idragresr, idragresp, idragrest
      WRITE (idragfile, ERR=100) dragtrmin, dragtrmax, dragtdr, 
     &     dragtpmin, dragtpmax, dragtdp, dragttmin, dragttmax, dragtdt
      WRITE (idragfile, ERR=100) planetmass, hmass

      do i = 1, idragresr
         do j = 1, idragresp
            do k = 1, idragrest
               WRITE (idragfile, ERR=100) dragenergy(i,j,k)
            enddo
         enddo
      enddo

      CLOSE (idragfile)

      RETURN
      
 100  print *, 'An error has occurred writing drag file'
      call quit(0) 

      end subroutine write_dragtable

