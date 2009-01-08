      SUBROUTINE montecarloRT(npart, xyzmh,vxyzu,trho,time)
      
      USE torus_mod, only : torus

      IMPLICIT NONE

      INCLUDE 'idim'
#ifdef MPI
      INCLUDE 'mpif.h'
#endif      

      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/files'

! Arguments
      INTEGER, INTENT(IN) :: npart
      REAL, INTENT(IN)    :: xyzmh(5,idim)
      REAL, INTENT(INOUT) :: vxyzu(4,idim)
      REAL*4, INTENT(IN)  :: trho(idim)
      REAL, INTENT(IN)    :: time

! Local variables
      INTEGER         :: i, ierr, num_gas
      REAL            :: temperature(idim)
      REAL            :: totalgasmass
      REAL, PARAMETER :: temp_min=10.0 ! Minimum allowed temperature
      LOGICAL, save   :: ifirst = .true.

! For conversion between internal energy and temperature
      REAL :: TempToU, uToTemp

! Time interval between torus calls
      REAL, parameter :: radTimeStep = 0.01
! Time of previous torus call
      REAL, save      :: lastRadTime      
! No. of pt mass on prev torus call
      INTEGER, save   :: lastNptmass       
      INTEGER :: nptmass_all
     
! fixed source parameters
      logical, parameter :: use_fixed_source=.true. 
      real, parameter    :: fix_source_R = 2.0    ! solar radii
      real, parameter    :: fix_source_L = 1.78   ! l_sun
      real, parameter    :: fix_source_T = 4000   ! kelvin

      TempToU = 3.0/2.0 * Rg/gmw/uergg
      uToTemp =  1.0 / TempToU

c--Only perform radiative transfer if ptmass exist, otherwise don't modify u()
     
#ifdef MPI
      CALL MPI_ALLREDUCE(nptmass, nptmass_all, 1, MPI_INTEGER, 
     &                   MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      nptmass_all = nptmass
#endif
! If this is the first step then remember the time
      IF (ifirst) THEN
         lastRadTime = time
         lastNptmass = nptmass_all
         ifirst = .false.
      END IF

! Call torus if time since last call exceeds limit or new pt mass created
      IF ( (nptmass_all > lastNptmass .or. 
     &     time-lastRadTime >= radTimeStep) .and. nptmass_all > 0 ) THEN

        WRITE (*,*) 'nptmass, lastNptmass, time, lastRadTime=', 
     &               nptmass_all, lastNptmass, time, lastRadTime
        lastRadTime = time
        lastNptmass = nptmass_all

c Find total gas mass and number of gas particles. 
         num_gas = 0
         totalgasmass = 0.0
         DO i=1,npart
            IF ( iphase(i) == 0 ) THEN 
               num_gas      = num_gas + 1
               totalgasmass = totalgasmass + xyzmh(4,i)
            END IF
         END DO
         WRITE(*,*) 'Total mass of gas particles=   ', totalgasmass

c Pass array of sph temperatures to torus to use as initial conditions
        DO i=1,npart
           temperature(i) = vxyzu(4,i) * uToTemp
        END DO

c--Call torus for Monte-Carlo radiative transfer     

        if ( use_fixed_source ) then 

           CALL torus(idim,  npart, nptmass, num_gas, 
     &          xyzmh, trho,  iphase, 
     &          udist, umass, utime, 
     &          time,  temperature, temp_min,
     &          totalgasmass, file1,              
     &          fix_source_R=fix_source_R, fix_source_L=fix_source_L,
     &          fix_source_T=fix_source_T, remove_radius=haccall*udist )

        else

           CALL torus(idim,  npart, nptmass, num_gas, 
     &          xyzmh, trho,  iphase, 
     &          udist, umass, utime, 
     &          time,  temperature, temp_min,
     &          totalgasmass, file1 )

        end if

c     
c--Set thermal energy of gas particles using torus's temperatures     
c
         DO i=1,npart
            vxyzu(4,i) =  max(temperature(i), temp_min) * TempToU 
         END DO

      ENDIF 

      END SUBROUTINE montecarloRT

