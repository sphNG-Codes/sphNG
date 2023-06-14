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

      INTEGER ii,jj,i,ipsendsetup(numproc)
      INTEGER ipsendall(numproc,numproc)

      IF (domain_checking) THEN
c--build pairwise list so that all procs can pair communicate in loop
         i = 0
         DO ii=1,numproc
            IF (ii .NE. iproc+1) THEN
               i = i + 1
               ipsendsetup(i) = ii
            ENDIF
         ENDDO
         ipsendsetup(numproc) = 0
         CALL MPI_ALLGATHER(ipsendsetup,numproc,MPI_INTEGER,
     &        ipsendall,numproc,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         nprocsend = numproc-1
         nprocsendall(1:numproc) = numproc-1
         CALL order_mpi_comms_bal(ipsendall)
         ipairwiseall(1:numproc) = ipsend(1:numproc)
c         WRITE(iprint,*) ipairwiseall
c--domain checking must find comms order in first call to densityiterate
         nprocsendlast = 0
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

      
c--communication ordering for SENDRECV implementation
c     pairs all ranks that need to communicate while ensuring that
c     no iteration in nloopmpi contains a rank more than once
      SUBROUTINE order_mpi_comms_bal(ipsendall)

#include "mpi_sup.h"
      
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidomains'
#ifdef MPIDEBUGDC
      INCLUDE 'COMMONS/logun'
#endif
      
      INTEGER ipsendall(numproc,numproc),ipsendinit(numproc,numproc)
      INTEGER indnrs(numproc),iplist(numproc),iflist(numproc)
      INTEGER iplocstart,iplocend
      INTEGER i,ii,j,jj,k
      INTEGER icompcount,npresent,imark,nfilled,iswap,ipst
      INTEGER iholder,ipsij
      LOGICAL rswb

      ipsendinit = ipsendall
      CALL insert_sort_index(numproc,nprocsendall,indnrs)

#ifdef MPIDEBUGDC
      WRITE (iprint,*) "START order_mpi_comms_bal"
      CALL debug_order_comms(ipsendall,indnrs)
#endif

c--This loop sorts the 2D array of rank send lists such that all      
c     ranks are paired with no duplicates for any given cycle
      icompcount = 0
      DO i=1,numproc
         IF (icompcount .EQ. numproc) EXIT         
         icompcount = 0
         npresent = 0
         iplist = 0
         iflist = 0
         nfilled = 0
         DO jj=1,numproc
            j = indnrs(jj)
c--if ij zero move to next row
            IF (ipsendall(i,j) .EQ. 0) THEN
               icompcount = icompcount + 1
               CYCLE
            ENDIF
            imark = 0
c--check that column has not already been filled by its pair            
            DO k=1,nfilled
               IF (iflist(k) .EQ. j) THEN
                  imark = 1
                  EXIT
               ENDIF
            ENDDO
            IF (imark .EQ. 1) THEN
               CYCLE
            ENDIF
c--check if send rank is not already present in this column
c     if it is re-order row to find rank that is not already present
c     if this is not possible set ij to 0
            DO ii=i,nprocsendall(j)
               iswap = 0
               DO k=1,npresent
                  ipst = ipsendall(ii,j)
                  IF (iplist(k) .EQ. ipst) THEN
                     iswap = 1
                  ELSEIF (iplist(k) .EQ. j) THEN
                     PRINT*,iproc,"current ij in iplist"
                     CALL quit(0)
                  ELSEIF (ipst .EQ. 0) THEN
                     PRINT*,iproc,"ipst 0",ii,i,j
                     CALL quit(0)
                  ENDIF
               ENDDO
               
               IF (iswap .EQ. 0) THEN
                  IF (i .EQ. ii) THEN
                     EXIT
                  ELSE
                     iholder = ipsendall(i,j)
                     ipsendall(i,j) = ipsendall(ii,j)
                     ipsendall(ii,j) = iholder
                     EXIT
                  ENDIF
               ELSEIF (ii .EQ. nprocsendall(j)) THEN
                  nprocsendall(j) = nprocsendall(j) + 1
                  DO k=nprocsendall(j),i+1,-1
                     ipsendall(k,j) = ipsendall(k-1,j)
                  ENDDO
                  ipsendall(i,j) = 0
                  EXIT
               ENDIF
            ENDDO

c--if ij was set to zero move to next column
            ipsij = ipsendall(i,j)
            IF (ipsij .EQ. 0)THEN 
               CYCLE
            ELSEIF (ipsendall(i,ipsij) .EQ. 0) THEN
               if (iproc.eq.0) PRINT*,"RANK SORT Warning A1",
     &              ipsij,ipsendall(i,ipsij)
               CALL debug_order_comms(ipsendinit,indnrs)
               CALL debug_order_comms(ipsendall,indnrs)
               CALL quit(0)
            ENDIF

c--add to list of ranks present in column
            IF (npresent .GT. numproc) PRINT*,"filled list overflow",
     &           iproc,npresent
            npresent = npresent + 1
            iplist(npresent) = ipsij
            IF (ipsij .LE. 0) THEN
               PRINT*,iproc,"RANK SORT Warning A2",ipsij,npresent,
     &              i,j,nprocsendall(j)
               PRINT *,"rowjb",i,j,nfilled,npresent
c               PRINT 99002,ipsendall(1:numproc,j)
               CALL debug_order_comms(ipsendinit,indnrs)
               CALL debug_order_comms(ipsendall,indnrs)
               CALL quit(0)
            ENDIF
c--set pair of current row directly
            rswb = .TRUE.
            IF (ipsendall(i,ipsij) .NE. j) THEN
               DO k=i+1,nprocsendall(ipsij)
                  IF (ipsendall(k,ipsij) .EQ. j) THEN
                     ipsendall(k,ipsij) = ipsendall(i,ipsij)
                     ipsendall(i,ipsij) = j
                     rswb = .FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF (rswb) THEN
                  PRINT*,"RANK SORT Warning B",
     &                 k,j,ipsij,nprocsendall(ipsij)
                  CALL debug_order_comms(ipsendinit,indnrs)
                  CALL debug_order_comms(ipsendall,indnrs)
                  CALL quit(0)
               ENDIF
            ENDIF
            
            npresent = npresent + 1
            iplist(npresent) = j
            DO k=1,nfilled
               IF (iflist(k) .EQ. ipsij) THEN
                  PRINT*,iproc,"RANK SORT Warning C",
     &                 k,ipsij,iflist(k)
                  CALL debug_order_comms(ipsendinit,indnrs)
                  CALL debug_order_comms(ipsendall,indnrs)
                  CALL quit(0)
               ENDIF
            ENDDO
            nfilled = nfilled + 1
            iflist(nfilled) = ipsij
            if (npresent .eq. numproc+1) THEN
               PRINT*,iproc,"RANK SORT Warning D",
     &              npresent,nfilled
               CALL debug_order_comms(ipsendinit,indnrs)
               CALL debug_order_comms(ipsendall,indnrs)
               CALL quit(0)
            ENDIF
         ENDDO
      ENDDO
      CONTINUE

c--build send list and check it's identical to recieve list
c     in general there is no need for them to be identical but this
c     implementation currently requires it
      iprec = 0
      iplocstart = 0
      iplocend = numproc + 1
      DO i=1,numproc
         ipsij = ipsendall(i,iproc+1)
         IF (ipsij .NE. 0) THEN
            iplocstart = iplocstart + 1
            ipsend(iplocstart) = ipsij-1
            IF (ipsendall(i,ipsij) .EQ. iproc+1) THEN
               iprec(iplocstart) = ipsij-1
            ELSE
               PRINT*, "ERROR: sendrecv mismatch A",iproc,
     &              iprec(i), ipsend(i)
               CALL quit(0)
            ENDIF
         ELSE
            iplocend = iplocend - 1
            ipsend(iplocend) = ipsij-1
            iprec(iplocend) = ipsij-1
         ENDIF
      ENDDO

#ifdef MPIDEBUGDC
      WRITE(iprint,*) " "
      WRITE (iprint,*) "END order_mpi_comms_bal"
      CALL debug_order_comms(ipsendall,indnrs)
      WRITE(iprint,*) " "
99002 FORMAT (I4)
      DO j=1,nprocsend
         WRITE(iprint,99002,advance="no") ipsend(j)
      ENDDO
      WRITE(iprint,*) " "
      DO j=1,nprocsend
         WRITE(iprint,99002,advance="no") iprec(j)
      ENDDO
      WRITE(iprint,*) " "
#endif
      DO i=1,numproc
         IF (iprec(i) .NE. ipsend(i)) THEN
            PRINT*, "ERROR: sendrecv mismatch B",iproc,
     &           iprec(i), ipsend(i)
            CALL quit(0)
         ENDIF
      ENDDO
      nloopmpi = nprocsend
      
      END SUBROUTINE order_mpi_comms_bal

      
c--communication ordering for SEND implementation
c     make groups of ranks that have no overlap in their communication. 
c     these groups can send during the same iteration over nloopmpi.
      SUBROUTINE order_mpi_comms_else(ipsendall,ipsendindx)

#include "mpi_sup.h"
      
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidomains'
#ifdef MPIDEBUGDC
      INCLUDE 'COMMONS/logun'
#endif
      
      INTEGER ipsendall(numproc,numproc)
      INTEGER ipsendindx(numproc),indnrs(numproc),iplist(numproc)
      INTEGER i,i1,i2,j,jj,k,kk
      INTEGER isendorder
      INTEGER iplisti1,niplist,niplistnew,isendproctemp
      LOGICAL complete(numproc),break
      
c--sort ranks high to low over number of ranks to communicate with         
      CALL insert_sort_index(numproc,nprocsendall,indnrs)
#ifdef MPIDEBUGDC
      WRITE(iprint,*) " "
      WRITE (iprint,*) "ipsend",nprocsend
      WRITE(iprint,99002) ipsend(1:nprocsend)
      WRITE (iprint,*) "iprec"
      WRITE(iprint,99002) iprec(1:nprocsend)
      WRITE (iprint,*) "indnrs"
      WRITE (iprint,99002) indnrs
      WRITE (iprint,*) "nprocsendall"
      WRITE (iprint,99002) nprocsendall(1:numproc)
      WRITE(iprint,*) " "
      DO i=1,numproc
         WRITE(iprint,99004) i,ipsendall(1:numproc,i)
      ENDDO
99004 FORMAT (1X,1(I3),"  ",1X,32(I3))
99002 FORMAT (1X, 32(I3))
#endif      
      complete = .FALSE.
      nloopmpi = 0
      isendorder = 0
      DO WHILE (isendorder .LT. numproc)
c--find an unassigned rank         
         DO jj=1,numproc
            j = indnrs(jj)
            IF (.NOT. complete(j)) EXIT
         ENDDO
c--mark it as complete add send ranks to group list         
         isendorder = isendorder + 1
         complete(j) = .TRUE.
         niplist = nprocsendall(j)
         iplist(1:niplist) = ipsendall(1:niplist,j)
#ifdef MPIDEBUGDC      
         WRITE(iprint,*) ""
         WRITE(iprint,99002,advance="no") iplist(1:niplist)
#endif
c--is this rank sending or receiving during this cycle over nloopmpi
c--if sending
         IF (ANY(ipsendall(1:nprocsendall(j),j) .EQ.
     &        iproc+1)) THEN
            nloopmpi = nloopmpi + 1
            iprec(nloopmpi) = j-1
            IF (j .EQ. iproc+1) THEN
               PRINT*,iproc,"ERROR: sending and receiving comms"
               CALL quit(0)
            ENDIF
c--elseif receiveing            
         ELSEIF (j .EQ. iproc+1) THEN
            nloopmpi = nloopmpi + 1
            isendproc = nloopmpi - 1
            iprec(nloopmpi) = -1
         ENDIF
c--loop over rows that have not yet been assigned (complete)
c     or test previously for this cycle
         DO kk=1,numproc
            k = indnrs(kk)
            IF (complete(k)) THEN
               CYCLE
            ENDIF
            break = .FALSE.
c--can rank k be added to the current group?
            DO i1=1,niplist
               iplisti1 = iplist(i1) 
               DO i2=1,nprocsendall(k)
                  IF (iplisti1 .EQ. ipsendall(i2,k)) THEN
                     break = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (break) EXIT
            ENDDO
            IF (break) THEN
c--rank k cannot be added to this group so iterate 
               CYCLE
            ELSE
c--add rank k to current group and its send ranks to this group's list
               niplistnew = niplist + nprocsendall(k)
               IF (niplistnew .GT. numproc) THEN
                  WRITE(*,*) iproc,"ERROR: niplist too big",niplistnew
                  CALL quit(0)
               ENDIF
               iplist(niplist+1:niplistnew) =
     &              ipsendall(1:nprocsendall(k),k)
#ifdef MPIDEBUGDC               
               WRITE(iprint,99002,advance="no")
     &              ipsendall(1:nprocsendall(k),k)
#endif
               niplist = niplistnew
c--mark k as an assigned rank          
               isendorder = isendorder + 1
               complete(k) = .TRUE.
c--if k communicates with this rank add to receive list
               IF (k .EQ. iproc+1) THEN
                  nloopmpi = nloopmpi + 1
                  isendproc = nloopmpi - 1
                  iprec(nloopmpi) = -1
               ENDIF
               IF (ANY(ipsendall(1:nprocsendall(k),k) .EQ.
     &              iproc+1)) THEN
                  nloopmpi = nloopmpi + 1
                  iprec(nloopmpi) = k-1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
      IF (isendorder .GT. numproc) THEN
         PRINT*,iproc,"ERROR: isendorder",isendorder,numproc
         CALL quit(0)
c      ELSEIF (nloopmpi .NE. nprocsend + 1) THEN
c         PRINT*,iproc,"ERROR: nloopmpi",nloopmpi,nprocsend
c         CALL quit(0)
      ENDIF

c--convert from fortran indexes to MPI ranks      
      DO j=1,nprocsend
         ipsend(j) = ipsendindx(j) - 1
      ENDDO
      
#ifdef MPIDEBUGDC
      WRITE(iprint,*) "order "
      WRITE(iprint,*) "isp,nloo,nps",isendproc,nloopmpi,nprocsend
      WRITE(iprint,99002) ipsend(1:nprocsend)
      WRITE(iprint,99002) iprec(1:nloopmpi)
#endif

      END SUBROUTINE order_mpi_comms_else



      SUBROUTINE debug_order_comms(ipsendall,indnrs)

#include "mpi_sup.h"
      
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidomains'
      INCLUDE 'COMMONS/logun'
      
      INTEGER ipsendall(numproc,numproc)
      INTEGER indnrs(numproc)
      INTEGER i,j
      
99004 FORMAT (1X,1(I4),"  ")
99002 FORMAT (I4)
      WRITE (iprint,*)"nprocsendall - number of comms for each MPI rank"
      DO j=1,numproc
         WRITE (iprint,99002,advance="no") nprocsendall(j)
      ENDDO
      WRITE(iprint,*) " "
      WRITE (iprint,*)"indnrs - ranks sorted by nprocsendall,",
     &     "most comms to least"
      DO j=1,numproc
         WRITE (iprint,99002,advance="no") indnrs(j)
      ENDDO
      WRITE(iprint,*) " "
      WRITE(iprint,*) " "
      DO i=1,numproc
         WRITE(iprint,99004,advance="no") i
         DO j=1,numproc
            WRITE(iprint,99002,advance="no") ipsendall(j,i)
         ENDDO
         WRITE(iprint,*) " "
      ENDDO

      END SUBROUTINE debug_order_comms
