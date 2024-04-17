c***************************************************************
c
c  This subroutine computes the dust fraction from the 
c     dust variable for one-fluid multi-grain dust
c     for a single SPH particle.
c
c  Written by MRB (2 Apr 2024)
c
c***************************************************************
      SUBROUTINE dustvar2dustfrac(rho1i,dustvari,dustfraci,
     &     dustfracisum)

      IMPLICIT NONE

      INCLUDE 'idim'
      
      REAL,    INTENT(IN)  :: rho1i
      REAL,    INTENT(IN)  :: dustvari(ndusttypes)
      REAL,    INTENT(OUT) :: dustfraci(ndusttypes), dustfracisum

      IF (idustFluid.EQ.1) THEN
         !--sqrt(rho*epsilon) method
         dustfraci(:) = MIN(dustvari(:)**2*rho1i,1.)
      ELSEIF (ABS(idustFluid).EQ.2) THEN
         !--sqrt(epsilon/(1-epsilon)) method (Ballabio et al. 2018)
         dustfraci(:) = dustvari(:)**2/(1.0+dustvari(:)**2)
      ELSEIF (idustFluid.EQ.3) THEN
         !--asin(sqrt(epsilon)) method (Hutchison)
         dustfraci(:) = SIN(dustvari(:))**2
      ELSEIF (idustFluid.EQ.4) THEN
         !--epsilon method                                                
         dustfraci(:) = dustvari(:)
      ELSE
         WRITE (*,*) 'ERROR - idustFluid ',idustFluid
         CALL quit(0)
      ENDIF

      dustfracisum = SUM(dustfraci)
 
      RETURN
      END SUBROUTINE dustvar2dustfrac
