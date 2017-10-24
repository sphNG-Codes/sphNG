      SUBROUTINE montecarloRT(npart, xyzmh,vxyzu,trho,time)
      
      USE torus_mod, only : torus

#ifdef MPIALL
#include "mpi_sup.h"
#endif

!      IMPLICIT NONE

      INCLUDE 'idim'
#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
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
      REAL, INTENT(IN)    :: xyzmh(5,mmax2)
      REAL, INTENT(INOUT) :: vxyzu(4,idim2)
      REAL*4, INTENT(IN)  :: trho(idim)
      REAL, INTENT(IN)    :: time

! Local variables
      INTEGER         :: i, num_gas, num_gas_all !, ierr
      REAL            :: temperature(idim)
      REAL            :: sum_tem, sum_tem_all 
      REAL            :: totalgasmass
      REAL, PARAMETER :: temp_min=10.0 ! Minimum allowed temperature
      LOGICAL, save   :: ifirst = .true.

! For conversion between internal energy and temperature
      REAL :: TempToU, uToTemp

! Time interval between torus calls
      REAL, parameter :: radTimeStep = 1.0
! Time of previous torus call
      REAL, save      :: lastRadTime      
! No. of pt mass on prev torus call
      INTEGER, save   :: lastNptmass       
      INTEGER :: nptmass_all
     
! fixed source parameters
      logical, parameter :: use_fixed_source=.false. 
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

        lastRadTime = time
        lastNptmass = nptmass_all

! Find total gas mass and number of gas particles. 
        num_gas = 0
        totalgasmass = 0.0
        DO i=1,npart
           IF ( iphase(i) == 0 ) THEN 
              num_gas      = num_gas + 1
              totalgasmass = totalgasmass + xyzmh(4,i)
           END IF
        END DO

! Pass array of sph temperatures to torus to use as initial conditions
        sum_tem = 0.0
        DO i=1,npart
           IF ( iphase(i) == 0 ) THEN
              temperature(i) = vxyzu(4,i) * uToTemp
              sum_tem = sum_tem + temperature(i)
           ELSE
              temperature(i) = temp_min
           END IF
        END DO

! Calculate mean temperature prior to Torus call
#ifdef MPI
        CALL MPI_ALLREDUCE(sum_tem, sum_tem_all, 1, MPI_REAL8, 
     &       MPI_SUM, MPI_COMM_WORLD, ierr) 

        CALL MPI_ALLREDUCE(num_gas, num_gas_all, 1, MPI_INTEGER, 
     &       MPI_SUM, MPI_COMM_WORLD, ierr) 

#else
        sum_tem_all = sum_tem
        num_gas_all = num_gas
#endif

#ifdef MPI
        if (iproc == 0) then
#endif
           write(*,*) "Mean gas temperature before torus=", 
     &          sum_tem_all / real(num_gas_all) 
#ifdef MPI
        end if
#endif

c--Call torus for Monte-Carlo radiative transfer     

        if ( use_fixed_source ) then 

c$$$           CALL torus(idim,  npart, nptmass, num_gas, 
c$$$     &          xyzmh, trho,  iphase, 
c$$$     &          udist, umass, utime, 
c$$$     &          time,  temperature,
c$$$     &          totalgasmass, file1,              
c$$$     &          fix_source_R=fix_source_R, fix_source_L=fix_source_L,
c$$$     &          fix_source_T=fix_source_T, remove_radius=haccall*udist )

           write(*,*) 
     &      "ERROR: Fixed source not available. See montecarloRT.F"

        else

           CALL torus(idim,  npart, nptmass, num_gas, 
     &          xyzmh, trho,  iphase, 
     &          udist, umass, utime, 
     &          time,  temperature,
     &          totalgasmass, file1 )

        end if

c$$$#ifdef MPI
c$$$        if (iproc == 0) then
c$$$#endif
c$$$        write(*,*) "Fixing T=20, u=", 20.0*TempToU
c$$$#ifdef MPI
c$$$        end if
c$$$#endif
c$$$        temperature(:) = 20.0


c     
c--Set thermal energy of gas particles using torus's temperatures     
c
        sum_tem = 0.0
         DO i=1,npart
            IF ( iphase(i) == 0 ) THEN 
               vxyzu(4,i) =  max(temperature(i), temp_min) * TempToU 
               sum_tem = sum_tem + temperature(i)
            END IF
         END DO

! Calculate mean temperature after Torus call. 
#ifdef MPI
        CALL MPI_ALLREDUCE(sum_tem, sum_tem_all, 1, MPI_REAL8, 
     &       MPI_SUM, MPI_COMM_WORLD, ierr) 

#else
        sum_tem_all = sum_tem
#endif

#ifdef MPI
        if (iproc == 0) then
#endif
           write(*,*) "Mean gas temperature after torus=", 
     &          sum_tem_all / real(num_gas_all) 
#ifdef MPI
        end if
#endif


      ENDIF 

      END SUBROUTINE montecarloRT

