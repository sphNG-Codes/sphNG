c***************************************************************
c
c  This subroutine initialises the dust-to-gas ratio and
c  grain size for the one-fluid dust.
c  Written by Mark Hutchison in 2018.
c
c***************************************************************
      SUBROUTINE initialise_dust(dust_to_gas,dustfrac,ndust)

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/tstop'
      INCLUDE 'COMMONS/dustfluidgrains'
      INCLUDE 'COMMONS/planetesimal'

      INTERFACE set_dustfrac_single
         SUBROUTINE set_dustfrac_single(dust_to_gas,dustfrac)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: dust_to_gas
         REAL, INTENT(OUT) :: dustfrac(:)
         END SUBROUTINE set_dustfrac_single
      END INTERFACE set_dustfrac_single

      INTERFACE set_dustfrac_power_law
         SUBROUTINE set_dustfrac_power_law(dust_to_gas_tot,dustfrac,
     &                                     smin,smax,sind)
         IMPLICIT NONE
         REAL, INTENT(IN)  :: dust_to_gas_tot,smin,smax,sind
         REAL, INTENT(OUT) :: dustfrac(:)
         END SUBROUTINE set_dustfrac_power_law
      END INTERFACE set_dustfrac_power_law
      
      INTEGER, INTENT(IN),  OPTIONAL :: ndust
      REAL,    INTENT(OUT), OPTIONAL :: dustfrac(:)
      REAL,    INTENT(OUT)           :: dust_to_gas
      REAL              :: psize,pdens
      INTEGER           :: i,ngas
      CHARACTER(len=1)  :: iok

99004 FORMAT (A1)
      WRITE (*, 99342)
99342 FORMAT(' Enter total initial dust to gas ratio',
     &       ' (rho_d/rho_g)')
      READ (iread, *) dust_to_gas
      !
      !--Set the numerical drag scheme
      !
      WRITE (*, 99343)
99343 FORMAT (' Enter the desired drag formulation : ',/,
     &     '           Epstein+Stokes (Laibe/Price/2012b) :  (1)',/,
     &     '             Stokes (Perets/Murray-Cley/2011) :  (2)',/,
     &     '  Epstein+Stokes (Baines/1968 + Brasser/2007) :  (3)',/,
     &     '                    Constant drag coefficient :  (4)',/,
     &     '                       Constant stopping time :  (5)')
      READ (iread,*) idragscheme
      PRINT *,'Selected drag scheme = ',idragscheme
      !
      !--Set the grain size/density and/or drag coefficient
      !
      SELECT CASE(idragscheme)
       CASE(1,2,3)
          !
          !--Set grain size and/or dust fraction
          !
          IF (ndusttypes.GT.1) THEN
99334        FORMAT(' Enter minimum grain size for',
     &              ' dust distribution (cm)')
             WRITE (*, 99334)
             READ (iread, *) smincgs
99335        FORMAT(' Enter maximum grain size for',
     &              ' dust distribution (cm)')
             WRITE (*, 99335)
             READ (iread, *) smaxcgs
99336        FORMAT(' Enter power-law index for dust (e.g. 3.5)')
             WRITE (*, 99336)
             READ (iread, *) sindex
          ELSE
99344        FORMAT(' Enter grain size for dust (cm)')
             WRITE (*, 99344)
             READ (iread, *) psize
             IF (idustFluid.NE.0) THEN
                sgrain(1) = psize
                sgrain    = sgrain/udist
             ELSE
                r_planetesimals(1) = psize
                r_planetesimals    = r_planetesimals/udist
             ENDIF
          ENDIF

99345     FORMAT(' Enter intrinsic grain density (g/cm^3)')
          WRITE (*, 99345)
          IF (ndusttypes.GT.1) THEN
             PRINT*, ' (currently the same for all dust types)'
          ENDIF
          READ (iread, *) pdens
          IF (idustFluid.NE.0) THEN
             densgrain(1) = pdens
             densgrain    = pdens/udens
          ELSE
             rho_planetesimals(1) = pdens
             rho_planetesimals    = rho_planetesimals/udens
          ENDIF
 
       CASE(4)
          WRITE (*, 99346)
99346     FORMAT(' Enter drag constant K (code units)')
          READ (iread, *) K_code
       CASE(5)
          WRITE (*, 99347)
99347     FORMAT(' Enter (constant) stopping time (code units)')
          READ (iread, *) K_code
       CASE DEFAULT
      END SELECT
      IF (idragscheme == 4 .OR. idragscheme == 5) THEN
         IF (ndusttypes.GT.1) THEN
             !--IMPORTANT: the values here don't matter, but they
             !  have to be defined and smin must equal smax
             sindex  = 1.
             smaxcgs = 1. 
             smincgs = smaxcgs
         ENDIF
 
         IF (K_code < 0.) WRITE (iprint,*) 'ERROR preset: K_code < 0' 
      ENDIF

      IF (idustFluid.NE.0) THEN
         !
         !--Set dust fraction
         !
         IF (ndusttypes.EQ.1) THEN
            CALL set_dustfrac_single(dust_to_gas,dustfrac)
         ELSE
            CALL set_dustfrac_power_law(dust_to_gas,dustfrac,
     &                                  smincgs,smaxcgs,sindex)
         ENDIF
         !
         !--Notify which dustvar parameterisation is being used
         !
         WRITE (*, 99349) idustFluid
99349    FORMAT (' Parameterisation for the dust fraction is: ',I2,/,
     &        '                           sqrt(epsilon*rho)  :   (1)',/,
     &        '                   sqrt(epsilon/(1-epsilon))  :(+/-2)',/,
     &        '                         asin(sqrt(epsilon))  :   (3)',/,
     &        '                                     epsilon  :   (4)')
      ELSE
         IF (present(ndust)) THEN
            IF (npart+ndust.GT.idim) THEN
               print *, 'idim too small for gas and planetesimals'
               print *, '...recompile with large idim'
               stop
            ENDIF
            ngas = npart
            npart = ngas + ndust

            DO i = ngas + 1, npart
               xyzmh(:,i) = xyzmh(:,i-ngas)
               iphase(i) = 11
               xyzmh(4,i) = dust_to_gas*xyzmh(4,i) 
            ENDDO
         ENDIF
      ENDIF
 
      END SUBROUTINE initialise_dust

c***************************************************************
c
c  Subroutine to set the dust fraction given the
c  dust-to-gas ratio. Equation (57) in Price & Laibe (2015)
c
c***************************************************************
      SUBROUTINE set_dustfrac_single(dust_to_gas,dustfrac)
      IMPLICIT NONE

      REAL, INTENT(IN)  :: dust_to_gas
      REAL, INTENT(OUT) :: dustfrac(:)

      dustfrac(:) = dust_to_gas/(1. + dust_to_gas)

      END SUBROUTINE set_dustfrac_single

c***************************************************************
c
c  Subroutine to set the dust fraction given the
c  dust-to-gas ratio.
c
c***************************************************************
      SUBROUTINE set_dustfrac_power_law(dust_to_gas_tot,dustfrac,
     &                                  smin,smax,sind)
      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/dustfluidgrains'

      INTERFACE set_grainsize
         SUBROUTINE set_grainsize(smin,smax,grid)
         REAL, INTENT(IN)  :: smin,smax
         REAL, OPTIONAL, INTENT(OUT) :: grid(:)
         END SUBROUTINE set_grainsize
      END INTERFACE set_grainsize

      REAL, INTENT(IN)  :: dust_to_gas_tot,smin,smax,sind
      REAL, INTENT(OUT) :: dustfrac(:)
      INTEGER :: i
      REAL :: dustfrac_tot
      REAL :: norm
      REAL :: rhodtot
      REAL :: grid(ndusttypes+1) = 0.
      REAL :: rhodusti(ndusttypes)
      REAL :: exact
      REAL :: power = 0.
      REAL, PARAMETER :: tol = 1.e-10

      !--reset global power-law index
      sindex = sind

      IF (smax==smin .or. ndusttypes==1) THEN
         !--If all the same grain size, just scale the dust fraction
         dustfrac(:) = dust_to_gas_tot/(1. + dust_to_gas_tot)
     &              *1./REAL(ndusttypes)
      ELSE
         CALL set_grainsize(smin,smax,grid)

         !--Dust density is computed from drhodust \propto dn*mdust
         !  where dn \propto s**(-p)*ds and mdust \propto s**(3). This
         !  is then integrated across each cell to account for
         !  mass contributions from unrepresented grain sizes
         DO i = 1, ndusttypes
            IF (sindex == 4.) THEN
               rhodusti(i) = LOG(grid(i+1)/grid(i))
            ELSE
               power = 4. - sindex
               rhodusti(i) = 1./power*(grid(i+1)**power-grid(i)**power)
            ENDIF
         ENDDO

         print *,'Hutchison min, max size ',smin,smax
         print *,'Hutchison grid bounds: ',grid(1:ndusttypes+1)

         !--Sum the contributions from each cell for total dust content
         rhodtot = SUM(rhodusti)

         !--Calculate the total dust fraction from the dust-to-gas ratio
         dustfrac_tot = dust_to_gas_tot/(1. + dust_to_gas_tot)

         !--Calculate the normalisation factor (\propto 1/rhotot) & scale
         !  the dust fractions. Note: dust density and dust fraction
         !  have the same power-law dependence on s.
         norm         = dustfrac_tot/rhodtot
         dustfrac(:)  = norm*rhodusti(:)

         !--Check to make sure the integral determining the
         !  contributions is correct
         IF (sindex == 4.) THEN
            exact = LOG(grid(ndusttypes+1)/grid(1))
         ELSE
            exact = 1./power*(grid(ndusttypes+1)**power-grid(1)**power)
         ENDIF
         IF (ABS(rhodtot-exact)/exact.GT.tol) THEN
            PRINT*,'Piecewise integration of MRN distribution not',
     &             ' matching the exact solution!'
            CALL quit(1)
         ENDIF
      ENDIF

      END SUBROUTINE set_dustfrac_power_law

c***************************************************************
c
c  Subroutine to set the grain size if spread of sizes, 
c     optionally returns an array of size bins in log space
c
c***************************************************************
      SUBROUTINE set_grainsize(smin,smax,grid)
      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/dustfluidgrains'

      REAL, INTENT(IN)  :: smin,smax
      REAL, OPTIONAL, INTENT(OUT) :: grid(:)
      INTEGER :: i
      REAL :: log_ds
      REAL :: log_grid(ndusttypes+1)
      REAL :: grainsizecgs(ndusttypes) = 0.

      smincgs = smin
      smaxcgs = smax

      ! check whether grid is passed in, and if so, that is large enough
      IF (PRESENT(grid)) THEN
         IF (SIZE(grid) < ndusttypes+1) THEN
            PRINT*, 'error trying to pass grid of insufficient size',
     &              ' to set_grainsize()'
         ENDIF
      ENDIF

      IF (ndusttypes==1) THEN
         !--This case should never occur, but for completeness...
         DO i = 1, ndusttypes
            WRITE(*,'(A)') ' Enter grain size for dust (cm)'
            READ (*,*) grainsizecgs(i)
         END DO
      ELSEIF (smax==smin .and. ndusttypes>1) THEN
         !--If all the same grain size, then just scale dustfrac
         grainsizecgs(:) = smax
      ELSE
         !--Create a uniform grid with N+1 points between smax and
         !  smin (inclusive)
         log_ds = LOG10(smax/smin)/real(ndusttypes)
         DO i = 1,ndusttypes+1
            log_grid(i) = LOG10(smin) + (i-1)*log_ds
         END DO

         !--Convert grid coordinates back to real space
         log_grid = 10.**log_grid

         !--Find representative s for each cell
         !  (skewed towards small grains because there are more small
         !  grains than large grains)
         DO i = 1,ndusttypes
            grainsizecgs(i) = SQRT(log_grid(i)*log_grid(i+1))
         END DO

         ! if we supplied grid, then return it
         IF (PRESENT(grid)) THEN
            !--this is no longer log of the grid at this point
            grid = log_grid 
         ENDIF
      ENDIF

      !--Set the grain properties relating to grain size
      sgrain(:)    = grainsizecgs(:)/udist
      !densgrain(:) = graindenscgs(:)/unit_density
      !massgrain(:) = 4./3.*pi*densgrain(:)*sgrain(:)**3

      END SUBROUTINE set_grainsize

c******************************************************************
c
c  Subroutine to number strings for N dust species
c
c******************************************************************
      SUBROUTINE nduststrings(pre_string,post_string,complete_string)
      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/dustfluidgrains'

      CHARACTER(LEN=*),   INTENT(IN)  :: pre_string,post_string
      CHARACTER(LEN=120), INTENT(OUT) :: complete_string(ndusttypes)
      
      INTEGER :: i,total_len,int_len
      CHARACTER(LEN=20) :: num_string

      int_len = 0
      IF (ndusttypes.GT.1) int_len = FLOOR(LOG10(REAL(ndusttypes)
     &                               + tiny)) + 1
      total_len = len(pre_string) + int_len  + len(post_string)
      IF (len(complete_string) < total_len) THEN
         PRINT*,'N dust string is not long enough!'
         CALL quit(0)
      ENDIF

      IF (ndusttypes.GT.1) THEN
         DO i = 1,ndusttypes
            WRITE(num_string,'(I0)') i
            WRITE(complete_string(i),'(A)') pre_string//
     &                       TRIM(ADJUSTL(num_string))//post_string
         ENDDO
      ELSE
         WRITE(complete_string,'(A)') pre_string//post_string
      ENDIF

      RETURN
      END SUBROUTINE nduststrings

