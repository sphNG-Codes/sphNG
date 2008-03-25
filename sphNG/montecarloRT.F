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

! Arguments
      INTEGER, INTENT(IN) :: npart
      REAL, INTENT(IN)    :: xyzmh(5,idim)
      REAL, INTENT(INOUT) :: vxyzu(4,idim)
      REAL*4, INTENT(IN)  :: trho(idim)
      REAL, INTENT(IN)    :: time

! Local variables
      INTEGER       :: i, igas, ierr
      REAL          :: temperature(npart), gaspartmass(npart)
      REAL          :: totalgasmass
      LOGICAL, save :: ifirst = .true.


      REAL :: TempToU, uToTemp

! Time interval between torus calls
      REAL, parameter :: radTimeStep = 0.05
! Time of previous torus call
      REAL, save      :: lastRadTime      
! No. of pt mass on prev torus call
      INTEGER, save   :: lastNptmass       
      INTEGER :: nptmass_all
      INTEGER :: num_gas
     
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

c Pass array of sph temperatures to torus to use as initial conditions

         IF ( encal /= 'v' ) THEN 
            print *, "sphtorus needs encal=v at present"
            STOP
         END IF

         igas = 0
         totalgasmass = 0.0
         DO i=1,npart
            IF ( iphase(i) == 0 ) THEN 
               igas = igas + 1
               temperature(igas) = vxyzu(4,i) * uToTemp
               gaspartmass(igas) = xyzmh(4,i)
               totalgasmass = totalgasmass + gaspartmass(igas)
            END IF
         END DO
         num_gas = igas

         WRITE(*,*) 'Total number of gas particles= ', num_gas
         WRITE(*,*) 'Total mass of gas particles=   ', totalgasmass

c     
c--Call torus for Monte-Carlo radiative transfer     
c

! Don't pass the whole array, just pass the bit that we're using
! What does trho do?
! Why not use mass from xyzmh variable?
! listpm is unused

         CALL torus(idim,npart,nactive,xyzmh,trho,iphase,nptmass,
     &       listpm,udist,umass,utime,time,gaspartmass(1),
     &       num_gas, temperature(1:num_gas))
c     
c--Set thermal energy of gas particles using torus's temperatures     
c     
         DO i=1, npart
            IF ( iphase(i) == 0 ) THEN 
               vxyzu(4,i) =  MAX(temperature(igas),10.0) * TempToU
            END IF
         END DO
         
      ENDIF 

      END SUBROUTINE montecarloRT

