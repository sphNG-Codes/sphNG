      SUBROUTINE domain_selection(xyzmh,nlstdo,ldolist,
     &     nldistantgn,listdistantgn,nsinknode,ntot,i5REAL8,i7REAL4)

#include "mpi_sup.h"
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidomains'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/treecom_P'

      INTEGER nlstdo,ldolist(idim),ntot
      INTEGER nldistantgn,listdistantgn(ngrav_nodes)
      INTEGER i5REAL8,i7REAL4
      REAL xyzmh(5,mmax2)

      INTEGER n,i,j,k,l,ipart,nld,ipair,ii
      INTEGER inow,iindex,ipsendi
      REAL xyzki,rcs,qrr,rr,acc

      INTEGER indexMPI_qrad,indexMPI5gn

      INTEGER nbothsend,nbothrec,nlstdopairrec,nlstdopair,npartsend
      INTEGER nprocrec,ipscount,isinkrec,nreccheck,isendsink,nsinknode
      INTEGER instart,isinknend,isendparts,ipstemp
      REAL point_cuboid_d2,cuboid_dist2,cdist2

      INTEGER nlstdoall(numproc),listsendgn(ngrav_nodes)
      INTEGER idispranks(numproc),ipsendbuf(numproc*numproc)
      INTEGER ipsendall(numproc,numproc),ipsendindx(numproc)
      INTEGER ipcount(numproc),ipcountorder(numproc)
      REAL de_all(6,numproc),csgas_all(6,numproc)
      REAL de_act(6),csgas_act(6),deact_all(6,numproc)
      LOGICAL order_comms,igravreverse,igravrevrec

      REAL t1,t2

      de_act = 1E30
      csgas_act = 1E30

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstdo,ldolist,radkernel)
C$OMP& shared(xyzmh,iphase)
C$OMP& private(n,k,l,ipart,rcs,xyzki)
C$OMP& reduction(MIN:de_act,csgas_act)
      DO n=1,nlstdo
         ipart = ldolist(n)
         IF (domain_checking) THEN
            rcs = csgrowthlim * radkernel * xyzmh(5,ipart)
            DO k = 1, 3
               l = k + 3
               xyzki = xyzmh(k,ipart)
               de_act(k) = MIN(de_act(k),xyzki)
               de_act(l) = MIN(de_act(l),-xyzki)
               IF (iphase(ipart) .EQ. 0 .OR.
     &              iphase(ipart) .GT. 10) THEN
                  csgas_act(k) = MIN(csgas_act(k),xyzki-rcs)
                  csgas_act(l) = MIN(csgas_act(l),rcs-xyzki)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C$OMP END PARALLEL DO

c
c--domain checking optimisations
c     before iteration begins need to establish which
c     MPI ranks this rank will send particles to
c     others will receive gravity nodes instead (if igphi.EQ.1)
c
      CALL MPI_ALLGATHER(domain_extent,6,MPI_REAL8,
     &     de_all,6,MPI_REAL8,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLGATHER(de_cs_gas,6,MPI_REAL8,
     &     csgas_all,6,MPI_REAL8,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLGATHER(de_act,6,MPI_REAL8,
     &     deact_all,6,MPI_REAL8,MPI_COMM_WORLD,ierr)

      CALL MPI_ALLGATHER(nlstdo,1,MPI_INTEGER,
     &     nlstdoall,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      nprocsend = 0
      ipsendindx = 0
      nldistantgn = 0
      instart = 2*ntot+3
      nsinknode = 0
      isinknend = mmax + ngrav_nodes
c
c--pairwise sendrecv
c
      DO ii=1,numproc-1
         ipair = ipairwiseall(ii)
         i = ipair + 1
         nlstdopair = nlstdoall(i)
c
c--if > 0 active particles, does target rank have potential neighbours?
c
         IF (nlstdo .GT. 0) THEN
            cdist2 = MIN(cuboid_dist2(de_act,csgas_all(1:6,i)),
     &           cuboid_dist2(csgas_act,de_all(1:6,i)))
c
c--will need to send gravity nodes if not recieving particles
c
         ELSEIF (nlstdopair .GT. 0) THEN
            cdist2 = 1E30
         ELSE
c
c--no particles active locally or on target rank so nothing to do
c
            CYCLE
         ENDIF
         IF (cdist2 .LT. tiny) THEN
c
c--potential neighbours so add rank to list to send particles to
c
            nprocsend = nprocsend + 1
            ipsendindx(nprocsend) = i
            igravreverse = .FALSE.
         ELSE
c
c--no potential neighbours so flag for possible reverse gravity
c
            igravreverse = .TRUE.
         ENDIF

         IF (igphi .EQ. 1) THEN
            isendparts = 0
            CALL MPI_SENDRECV(igravreverse,1,MPI_LOGICAL,ipair,29,
     &           igravrevrec,1,MPI_LOGICAL,ipair,29,
     &           MPI_COMM_WORLD,istatus,ierr)
c
c--if > 0 active particles on target rank and not receiving particles,
c     find gravity nodes that need to be sent
            IF (nlstdopair .GT. 0 .AND. igravrevrec) THEN
               CALL grav_only_treef(deact_all(1:6,i),xyzmh,acc,
     &              npartsend,nbothsend,listsendgn(:),0.0)
c
c--if gravity nodes to be sent include a tree leaf (particle)
c
               IF(npartsend.GT.0) THEN
c
c--set flag to receive all particles from rank i instead
c
                  nlstdopair = -1
                  nbothsend = 0
c elseif # of gravity nodes is large compared to # of particles
c               ELSEIF(nlstdopair .LT. 10*nbothsend) THEN
c                  nlstdopair = -1
c                  nbothsend = 0
               ENDIF
            ELSE
c
c--if no foreign active particles, no gravity nodes need be sent
c
               nbothsend = 0
            ENDIF

            CALL MPI_SENDRECV(nlstdopair,1,MPI_INTEGER,ipair,30,
     &           nlstdopairrec,1,MPI_INTEGER,ipair,30,
     &           MPI_COMM_WORLD,istatus,ierr)
c
c--not receiving gravity nodes therefore send all particles
c
            IF (nlstdopairrec .EQ. -1) THEN
               nprocsend = nprocsend + 1
               ipsendindx(nprocsend) = i
#ifdef MPIDEBUGDC
               PRINT*,iproc,": WARNING - sending extra rank",i,
     &              npartsend,nbothsend
#endif
            ENDIF

            IF (nbothsend .LT. 0) THEN
               PRINT*,iproc,": ERROR - must send at least one node",
     &              npartsend,nbothsend
               CALL quit(1)
            ELSEIF (instart+nbothsend+nsinknode.GE.isinknend-1) THEN
               PRINT*,': ERROR - ngrav_nodes too small',
     &              ntot,nbothrec
               CALL quit(1)
c            ELSEIF (nbothsend.GE.30) THEN
c               WRITE(iprint,*)"sending 30+gnodes",nbothsend
            ENDIF
c
c--don't send gravity if no active paticles at target or
c     particles being sent
c
            IF (nlstdopair .LE. 0 .AND. nlstdopairrec .LE. 0) THEN
               CYCLE
            ELSE
c
c--Check if standalone tree node would ideally have been sub-divided
c
               IF (nptmass .GT. 0) THEN
                  rr = point_cuboid_d2(xyzmh(1:3,nroot+1),
     &                 deact_all(1:6,i))
                  qrr = qrad(1,nroot+1)**2
                  IF (0.7**2*rr.LT.qrr) THEN
                     WRITE(iprint,*) "Normally would open node here",
     &                    0.7**2*rr,qrr,rr
                  ENDIF
                  isendsink = 1
               ELSE
                  isendsink = 0
               ENDIF
               CALL insert_sort(nbothsend,listsendgn(1:nbothsend))

               CALL MPI_TYPE_INDEXED(nbothsend,
     &              lblocklengths,listsendgn,
     &              i5REAL8,indexMPI5gn,ierr)
               CALL MPI_TYPE_COMMIT(indexMPI5gn,ierr)

               CALL MPI_TYPE_INDEXED(nbothsend,
     &              lblocklengths,listsendgn,
     &              i7REAL4,indexMPI_qrad,ierr)
               CALL MPI_TYPE_COMMIT(indexMPI_qrad,ierr)

               CALL MPI_SENDRECV(xyzmh,1,indexMPI5gn,ipair,28,
     &              xyzmh(1,instart), ngrav_nodes, i5REAL8,
     &              ipair, 28, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,i5REAL8,nbothrec,ierr)

               CALL MPI_SENDRECV(qrad,1,indexMPI_qrad,ipair,27,
     &              qrad(1,instart), ngrav_nodes, i7REAL4,
     &              ipair, 27, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,i7REAL4,nreccheck,ierr)
               IF (nbothrec.NE.nreccheck) THEN
                  PRINT*,iproc,': ERROR - nbothrec.NE.nreccheck',
     &                 nbothrec,nreccheck
                  CALL quit(1)
               ENDIF
c
c--receive sink node at end of qrad array.
c
               CALL MPI_SENDRECV(xyzmh(1,nroot+1),isendsink,
     &              i5REAL8,ipair,29,
     &              xyzmh(1,isinknend-nsinknode), 1, i5REAL8,
     &              ipair, 29, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,i5REAL8,isinkrec,ierr)
               IF (isinkrec .NE. 0 .AND. isinkrec .NE. 1) THEN
                  PRINT*,iproc,': ERROR - isinkrec NE 1or0',isinkrec
                  CALL quit(1)
               ENDIF

               CALL MPI_SENDRECV(qrad(1,nroot+1),isendsink,
     &              i7REAL4,ipair,30,
     &              qrad(1,isinknend-nsinknode), 1, i7real4,
     &              ipair, 30, MPI_COMM_WORLD, istatus, ierr)

               CALL MPI_TYPE_FREE(indexMPI5gn,ierr)
               CALL MPI_TYPE_FREE(indexMPI_qrad,ierr)

               DO nld=nldistantgn+1,nldistantgn+nbothrec
                  listdistantgn(nld) = instart
                  instart = instart + 1
               ENDDO
               nldistantgn = nldistantgn + nbothrec
               nsinknode = nsinknode + isinkrec
            ENDIF
         ENDIF
      ENDDO
      instart = nldistantgn
      DO k=isinknend-nsinknode,isinknend
         instart = instart + 1
         listdistantgn(instart) = k
      ENDDO
      CALL insert_sort(nprocsend,ipsendindx(1:nprocsend))
c
c--sort send order such that common ranks are prioritised
c
      IF (ideniteratebal) THEN
         ipcount = 0
         DO i=1,nprocsend
            ipcount(ipsendindx(i)) = 1
         ENDDO
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,ipcount,numproc,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
         CALL insert_sort_index(numproc,ipcount,ipcountorder)
         nloopmpi = ipcount(ipcountorder(1))
         inow = 1
         DO j=1,numproc
            iindex = ipcountorder(j)
            DO i=1,nprocsend
               ipsendi = ipsendindx(i)
               IF (iindex .EQ. ipsendi) THEN
                  ipsendindx(i) = ipsendindx(inow)
                  ipsendindx(inow) = ipsendi
                  inow = inow + 1
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      nprocsendmax = MAX(nprocsendmax,nprocsend)
      nprocsendmin = MIN(nprocsendmin,nprocsend)
c
c--have to re-order MPI comms if any rank has changed its send list
c
      order_comms = .TRUE.
      IF (nprocsend .EQ. nprocsendlast) THEN
         order_comms = .FALSE.
         DO i=1,nprocsend
            IF (ipsendindx(i) .NE. ipsendlast(i)) THEN
               order_comms = .TRUE.
               EXIT
            ENDIF
         ENDDO
      ENDIF
      DO i=1,nprocsend
         ipsendlast(i) = ipsendindx(i)
      ENDDO
      nprocsendlast = nprocsend

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,order_comms,1,
     &     MPI_INTEGER,MPI_LOR,MPI_COMM_WORLD,ierr)

      IF (order_comms) THEN
         CALL MPI_ALLGATHER(nprocsend,1,MPI_INTEGER,
     &        nprocsendall,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         idispranks(1) = 0
         DO i=1,numproc-1
            idispranks(i+1) = idispranks(i) + nprocsendall(i)
         ENDDO
         CALL MPI_ALLGATHERV(ipsendindx,nprocsend,MPI_INTEGER,
     &        ipsendbuf,nprocsendall,idispranks,
     &        MPI_INTEGER,MPI_COMM_WORLD,ierr)
c
c--build array of which ranks are sending to which
c
         ipscount = 0
         nprocrec = 0
         iprec = 0
         DO j=1,numproc
            DO i=1,numproc
               IF (i .LE. nprocsendall(j)) THEN
                  ipscount = ipscount + 1
                  ipstemp = ipsendbuf(ipscount)
                  ipsendall(i,j) = ipstemp
                  IF (ipstemp .EQ. iproc+1) THEN
                     nprocrec = nprocrec + 1
                     iprec(nprocrec) = j
                  ENDIF
               ELSE
                  ipsendall(i,j) = 0
               ENDIF
            ENDDO
         ENDDO
         IF (ideniteratebal) THEN
            CALL order_sendrecv(ipsendall)
         ELSE
            CALL order_mpi_comms_else(ipsendall,ipsendindx)
         ENDIF
      ENDIF

      END SUBROUTINE domain_selection
