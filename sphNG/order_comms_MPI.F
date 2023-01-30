c--initialise send and recieve lists for MPI
c     when not using domain_checking these lists allow densityiterate
c     to loop over all other MPI ranks.
c     if using domain_checking these lists are built each timestep and
c     do not need to be initialised.
      SUBROUTINE init_mpi_loops

#include "mpi_sup.h"

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'igrape'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidomains'
      INCLUDE 'COMMONS/logun'

      INTEGER ii,jj

      IF (domain_checking) THEN
         CONTINUE
      ELSE
         IF (ideniteratebal) THEN
            nloopmpi = numproc - 1
            DO ii = 1, numproc - 1
               ipsend(ii) = MOD(iproc+ii,numproc)
               iprec(ii) = MOD(numproc+iproc-ii,numproc)
            ENDDO
         ELSE
            nloopmpi = numproc
            nprocsend = numproc - 1
            isendproc = iproc
            jj = 0
            DO ii = 0, numproc - 1
               IF (iproc .NE. ii) THEN
                  jj = jj + 1
                  ipsend(jj) = ii
               ENDIF
               iprec(ii+1) = ii
            ENDDO
         ENDIF
      ENDIF
      
      END SUBROUTINE init_mpi_loops
