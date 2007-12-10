      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &     dvxyzu,dha,npart,ntot,ireal,alphaMM,ekcle,Bxyz,dBxyz)
c************************************************************
c                                                           *
c  This subroutine drives the computation of the forces on  *
c     every particle on the list.                           *
c                                                           *
c************************************************************
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION xyzmh(5,mmax),vxyzu(4,idim),dvxyzu(4,idim)
      REAL*4 dha(1+isizealphaMM,idim),alphaMM(isizealphaMM,idim)
      DIMENSION ireal(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bxyz(3,imhd),dBxyz(3,imhd)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/hagain'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/avail'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/active'

c      INCLUDE 'COMMONS/timeextra'

#ifdef MPI
      PARAMETER (ineighproc = idim/2)
      DIMENSION llistsend(idim), llistrec(idim),
     &     nneighsentback(nummaxproc),
     &     nneightogetback(nummaxproc),
     &     llistsentback(ineighproc,nummaxproc),
     &     inumberindiv(nummaxproc),
     &     inumbercumm(nummaxproc),
     &     inumberproc(nummaxproc)
      DIMENSION nneighrec(ineighproc)
      DIMENSION lblocklengths(idim)
      REAL*4 rhorec(ineighproc), divvrec(ineighproc)
#endif

      CHARACTER*7 where
      DIMENSION dedxyz(3,iradtrans)

      DATA where/'derivi'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine derivi')
c
c--Set constants first time around
c
      uradconst = radconst/uergcc
      nlst_in = 1
      nlst_end = nlst
      IF (itrace.EQ.'all') WRITE (iprint, 99002) nlst_in, nlst_end
99002 FORMAT(' derivi ',I8, I8)
c
c--Compute the neighbor indexes & gravitational forces of the distant 
c     particles for all the particles in the list
c
      IF (igrape.EQ.0) THEN
         IF (nlst_end.GT.nptmass .OR. itreeupdate .OR. 
     &        iptintree.EQ.2) THEN
            CALL insulate(3, ntot, npart, xyzmh, dvxyzu)
c
c--Keep using tree until all sinks have been done after accretion event
c     (important if sinks have individual timesteps)
c
            IF (nlst_end.GE.nptmass) itreeupdate = .FALSE.
         ELSE
c
c--Don't bother to update gravity from gas particles acting on sinks
c    Also means potential energy and neighbours of sinks are not updated
            DO i = nlst_in, nlst_end
               ipart = llist(i)
               dvxyzu(1,ipart) = gravxyzstore(1,ipart)
               dvxyzu(2,ipart) = gravxyzstore(2,ipart)
               dvxyzu(3,ipart) = gravxyzstore(3,ipart)
               poten(ipart) = potenstore(ipart)
            END DO
         ENDIF
      ELSEIF (igrape.EQ.1) THEN
         CALL insulate(4, ntot, npart, xyzmh, dvxyzu)
      ELSE
         CALL error(where,1)
      ENDIF
c
c--Find neighbours and calculate gravity on and from point masses
c     The point masses are no longer in the TREE/GRAPE - done separately for
c     higher accuracy
c
      IF (nptmass.GT.0) THEN
         IF (iptintree.EQ.0) THEN
            CALL gptall(xyzmh,npart,dvxyzu)
         ELSEIF (iptintree.EQ.1) THEN
            CALL gforspt(nlst_end,llist,npart,ntot,xyzmh,dvxyzu)
         ENDIF
      ENDIF
c
c--Only want to do density and updates for non-sinks
c
c      IF (nlst_end.GT.nptmass) THEN
c
c--Compute the pressure, divv, etc. on list particles but
c     do not search twice the list particles
c
         IF (itiming) CALL getused(tdens1)

         CALL densityi(npart,ntot,xyzmh,vxyzu,ekcle,
     &                 nlst_in,nlst_end,llist,itime)

         IF (itiming) THEN
            CALL getused(tdens2)
            tdens = tdens + (tdens2 - tdens1)
         ENDIF
c
c--Predict the pressure, divv, etc. on list-particle neighbors
c
         IF (itiming) CALL getused(td11)
         nlistavail = 0
         DO i = 1, nlst
            ipart = llist(i)
            DO j = 1, nneigh(ipart)
               IF (j.GE.nlmax) THEN
                  jpart = neighover(j-nlmax+1,ABS(neighb(nlmax,ipart))) 
               ELSE
                  jpart = neighb(j,ipart)
               ENDIF
               IF (jpart.GT.idim) THEN
                  WRITE (*,*) 'jpart.GT.idim'
                  WRITE (iprint,*) 'jpart.GT.idim'
                  CALL quit
               ENDIF
               IF (jpart.LT.1) THEN
                  WRITE (*,*) 'jpart.LT.1'
                  WRITE (iprint,*) 'jpart.LT.1'
                  CALL quit
               ENDIF
               IF (iavail(jpart).EQ.0 .AND.(.NOT.iscurrent(jpart))) THEN
                  iavail(jpart) = 1
                  nlistavail = nlistavail + 1
                  IF (nlistavail+nlst.GT.idim) THEN
                     WRITE (*,*) 'nlistavail+nlst.GT.idim'
                     WRITE (iprint,*) 'nlistavail+nlst.GT.idim'
                     CALL quit
                  ENDIF
                  llist(nlst+nlistavail) = jpart
               ENDIF
            END DO
         END DO
         IF (itiming) THEN
            CALL getused(td12)
            td1 = td1 + (td12 - td11)
         ENDIF

#ifdef MPI
#ifdef MPIDEBUG
         print *,iproc,': Start MPI'
#endif
         CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
         CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
         DO i = 1, idim
            lblocklengths(i) = 1
         END DO
         DO i = 1, nlst
            ipos = llist(i)
            llistsend(i) = ipos - 1
            nneighfull(ipos) = nneigh(ipos)
         END DO
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlst, 1, llistsend, i5REAL8,
c     &        indexMPI5, ierr)
         CALL MPI_TYPE_INDEXED(nlst, lblocklengths, llistsend, i5REAL8,
     &        indexMPI5, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlst, 1, llistsend, i4REAL8,
c     &        indexMPI4, ierr)
         CALL MPI_TYPE_INDEXED(nlst, lblocklengths, llistsend, i4REAL8,
     &        indexMPI4, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

         inumofreturns = 0
         inumofsends = 0
         maxnneighsentback = 0
c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
         DO i = 1, numproc - 1
            iahead = MOD(iproc+i,numproc)
            ibehind = MOD(numproc+iproc-i,numproc)
c
c--Send to node ahead, receive from node behind
c
#ifdef MPIDEBUG
            print *,iproc,': Sending xyzmh to ',iahead,' rec from ',
     &           ibehind
#endif

            CALL MPI_SENDRECV(xyzmh,1,indexMPI5,iahead,0,
     &           xyzmh(1,2*ntot+3), idim, i5REAL8, ibehind,
     &           0, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)

            CALL MPI_SENDRECV(vxyzu,1,indexMPI4,iahead,1,
     &           vxyzu(1,ntot+1), idim, i4REAL8, ibehind,
     &           1, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i4REAL8, icheck, ierr)
            IF (icheck.NE.inumber) THEN
               WRITE (*,*) 'ERROR - icheck.NE.inumber 1 ',iproc,
     &              ibehind,icheck,inumber
               CALL quit
            ENDIF

#ifdef MPIDEBUG
            print *,iproc,': Received ',inumber
#endif
c
c--This should return neighbours lists for the particles as well as 
c     gravitational forces from nodes and leaves
c
            inumreturn = 0
            DO j = 1, inumber
               jpart = j + ntot
               IF (jpart.GT.idim) THEN
                  WRITE (*,*) 'jpart.GT.idim'
                  CALL quit
               ENDIF
               iphase(jpart) = 0
               CALL treef(jpart,npart,ntot,xyzmh,acc,igphi,
     &              fsx,fsy,fsz,epot)
               dvxyzu(1,jpart) = fsx
               dvxyzu(2,jpart) = fsy
               dvxyzu(3,jpart) = fsz
               dvxyzu(4,jpart) = epot
c
c--Make a list of indices to return densities
c
               IF (nneigh(jpart).GE.1) THEN
                  inumreturn = inumreturn + 1
                  llistrec(inumreturn) = jpart
               ENDIF
            END DO
            ntotplusinumber = ntot + inumber
#ifdef MPIDEBUG
            print *,iproc,': number found neigh = ',inumreturn
#endif
            nneightogetback(ibehind+1) = inumreturn
            IF (inumreturn.GT.0) THEN
               inumofreturns = inumofreturns + 1
c
c--The neighbours lists should then be used to calculate contributions to
c     rho and divv (and curlv if required) and add them onto (rho,divv)passed
c     -1 is flag for using passed quantities
c     Only need to do this if there are neighbours!
c
               CALL densityi(npart,ntot,xyzmh,vxyzu,ekcle,1,
     &              inumreturn,llistrec,itime)
#ifdef MPIDEBUG
               print *,iproc,' done density'
#endif
c
c--Need to mark neighbours of remote particles to have their rho etc updated
c
               IF (itiming) CALL getused(td11)
               DO l = 1, inumreturn
                  ipart = llistrec(l)
                  DO j = 1, nneigh(ipart)
                     IF (j.GE.nlmax) THEN
                        jpart = neighover(j-nlmax+1,
     &                       ABS(neighb(nlmax,ipart))) 
                     ELSE
                        jpart = neighb(j,ipart)
                     ENDIF
                     IF (iavail(jpart).EQ.0 .AND.
     &                    (.NOT.iscurrent(jpart))) THEN
                        iavail(jpart) = 1
                        nlistavail = nlistavail + 1
                        IF (nlistavail+nlst.GT.idim) THEN
                           WRITE (*,*) 'nlistavail+nlst.GT.idim'
                           WRITE (iprint,*) 'nlistavail+nlst.GT.idim'
                           CALL quit
                        ENDIF
                        llist(nlst+nlistavail) = jpart
                     ENDIF
                  END DO
               END DO
               IF (itiming) THEN
                  CALL getused(td12)
                  td1 = td1 + (td12 - td11)
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': set iavail'
#endif
c
c--Return indices
c
#ifdef MPIDEBUG
               print *,iproc,': modifying llistrec'
#endif
               DO j = 1, inumreturn
                  llistrec(j) = llistrec(j) - ntot
               END DO
c
c--End of work to do on neighbours
c
            ENDIF
#ifdef MPIDEBUG
            print *,iproc,': sending indices to ',ibehind,
     &           inumreturn
#endif

            CALL MPI_SENDRECV(llistrec,inumreturn, MPI_INTEGER, ibehind,
     &           10, llistrec(inumreturn+1), nlst, MPI_INTEGER, iahead,
     &           10, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_INTEGER,inumberneigh,ierr)
#ifdef MPIDEBUG
            print *,iproc,': returned indices to ',ibehind,' and got ',
     &           inumberneigh,' from ',iahead
#endif
c
c--Return densities of particles with density contributions (if any)
c
c--Note - need to change values of llistrec because in C array indices
c     start from 0 not 1 and MPI assumes C-type indexing
c
            DO j = 1, inumreturn
               llistrec(j) = llistrec(j) + ntot - 1
            END DO

c            CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1,
c     &           llistrec, MPI_REAL4, indexMPI1return, ierr)
            CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths,
     &           llistrec, MPI_REAL4, indexMPI1return, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI1return,ierr)

c            CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1, 
c     &           llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
            CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths, 
     &           llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
            CALL MPI_TYPE_COMMIT(indexMPI_I1return,ierr)

            CALL MPI_SENDRECV(rho, 1, indexMPI1return, ibehind, 11, 
     &           rhorec, inumberneigh, MPI_REAL4, iahead, 11, 
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
            IF (icheck.NE.inumberneigh) THEN
               WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &              iproc
               CALL quit
            ENDIF
            
#ifdef MPIDEBUG
            print *,iproc,': sent rho'
#endif

            CALL MPI_SENDRECV(divv, 1, indexMPI1return, ibehind, 12, 
     &           divvrec, inumberneigh, MPI_REAL4, iahead, 12, 
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
            IF (icheck.NE.inumberneigh) THEN
               WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &              iproc
               CALL quit
            ENDIF

#ifdef MPIDEBUG
            print *,iproc,': sent divv'
#endif

            CALL MPI_SENDRECV(nneigh, 1 ,indexMPI_I1return, ibehind, 13,
     &           nneighrec, inumberneigh, MPI_INTEGER, iahead, 13, 
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
            IF (icheck.NE.inumberneigh) THEN
               WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &              iproc
               CALL quit
            ENDIF
#ifdef MPIDEBUG
            print *,iproc,': sent nneigh'
#endif

            CALL MPI_TYPE_FREE(indexMPI_I1return,ierr)
            CALL MPI_TYPE_FREE(indexMPI1return,ierr)
            IF (inumberneigh.GT.0) inumofsends = inumofsends + 1
#ifdef MPIDEBUG
            print *,iproc,': all sent'
#endif
            DO k = 1, inumberneigh
               ipos = llist(llistrec(k + inumreturn))
               rho(ipos) = rho(ipos) + rhorec(k)
               divv(ipos) = divv(ipos) + divvrec(k)
c--Note that cannot change nneigh() because this still needs to have the number
c     of particles in the *local* neighbours list only.  nneighfull needed for
c     hdot.
               nneighfull(ipos) = nneighfull(ipos) + nneighrec(k)
               llistsentback(k,iahead+1) = ipos - 1
#ifdef MPIDEBUG
#endif
            END DO

#ifdef MPIDEBUG
            print *,iproc,': received neighbour data from ',iprocrec
#endif
            nneighsentback(iahead+1) = inumberneigh
            maxnneighsentback = MAX(maxnneighsentback,
     &           inumberneigh)
c
c--Send gravity forces back to sending process
c
#ifdef MPIDEBUG
            print *,iproc,' sending forces back to ',ibehind
#endif
            CALL MPI_SENDRECV(dvxyzu(1,ntot+1),4*inumber,
     &           MPI_REAL8,ibehind,14,
     &           dvxyzu(1,ntotplusinumber+1),4*nlst,MPI_REAL8,
     &                 iahead, 14, 
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
            IF (ireturned.NE.4*nlst) THEN
               WRITE (*,*) 'ERROR - ireturned.NE.4*nlst dvxyz'
               CALL quit
            ENDIF

            DO l = 1, nlst
               ipos = llist(l)
               jpos = ntotplusinumber + l
               DO k = 1, 3
                  dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
               END DO
               poten(ipos) = poten(ipos) + dvxyzu(4,jpos)
            END DO

#ifdef MPIDEBUG
            print *,iproc,' sent forces'
#endif

         END DO



c
c--IMPLEMENTATION USING SEPARATE SEND/RECV's FARMED OUT BY 1 PROCESS AT A TIME
c
         GOTO 4000

         DO i = 0, numproc - 1
c            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
            IF (iproc.EQ.i) THEN

c
c--Send active node data to be processed by other processes.  Unfortunately
c     this means that the local process doesn't do much, but as N->infinity
c     this will not matter.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     CALL MPI_SEND(xyzmh,1,indexMPI5,j,0,
     &                    MPI_COMM_WORLD, ierr)
                     CALL MPI_SEND(vxyzu,1,indexMPI4,j,1,
     &                    MPI_COMM_WORLD, ierr)
                  ENDIF
               END DO
c
c--Receive back density, divv, nneigh contributions to be added on to 
c     local value
c
               inumofsends = 0
               DO j = 1, numproc - 1
#ifdef MPIDEBUG
                  print *,iproc,': receiving neighbour data'
#endif
                  CALL MPI_RECV(llistrec,nlst,MPI_INTEGER,
     &                 MPI_ANY_SOURCE,10,MPI_COMM_WORLD,istatus,ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_INTEGER,inumberneigh,
     &                 ierr)
                  iprocrec = istatus(MPI_SOURCE)

                  IF (inumberneigh.GT.ineighproc) THEN
                     WRITE (*,*) 'ERROR - inumberneigh.GT.ineighproc',
     &                    inumberneigh,ineighproc,' from ',iprocrec
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': received list of neighbours from ',
     &                 iprocrec, inumberneigh
#endif

                  IF (inumberneigh.GT.0) THEN
                     inumofsends = inumofsends + 1
                     CALL MPI_RECV(rhorec,inumberneigh,MPI_REAL4,
     &                    iprocrec, 11, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,
     &                 ierr)
                     IF (icheck.NE.inumberneigh) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': received rho ',icheck
#endif

                     CALL MPI_RECV(divvrec,inumberneigh,MPI_REAL4,
     &                    iprocrec, 12, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,
     &                 ierr)
                     IF (icheck.NE.inumberneigh) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': received divv ',icheck
#endif

                     CALL MPI_RECV(nneighrec,inumberneigh,MPI_INTEGER,
     &                    iprocrec, 13, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_INTEGER,icheck,
     &                 ierr)
                     IF (icheck.NE.inumberneigh) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': received nneigh ',icheck,iprocrec
#endif

c                  goto 300

                     DO k = 1, inumberneigh
                        ipos = llist(llistrec(k))
                        rho(ipos) = rho(ipos) + rhorec(k)
                        divv(ipos) = divv(ipos) + divvrec(k)
c--Note that cannot change nneigh() because this still needs to have the number
c     of particles in the *local* neighbours list only.  nneighfull needed for
c     hdot.
                     nneighfull(ipos) = nneighfull(ipos) + nneighrec(k)
                        llistsentback(k,iprocrec+1) = ipos - 1
#ifdef MPIDEBUG
                        IF (iproc.EQ.0 .AND. ipos.EQ.1300) THEN
                           print *,iproc,' CEN nneighrec ',iprocrec,
     &                          nneighrec(k),k
                        ENDIF
#endif
                     END DO
 300                 continue
#ifdef MPIDEBUG
               print *,iproc,': received neighbour data from ',iprocrec
#endif
                  ENDIF
                  nneighsentback(iprocrec+1) = inumberneigh
                  maxnneighsentback = MAX(maxnneighsentback,
     &                 inumberneigh)
c
c--Receive back gravity-only forces from remote process
c
#ifdef MPIDEBUG
                  print *,iproc,': receiving forces into ',ntot+1,nlst
#endif
                  CALL MPI_RECV(dvxyzu(1,ntot+1),4*nlst,MPI_REAL8,
     &                 MPI_ANY_SOURCE, 14, MPI_COMM_WORLD, istatus,
     &                 ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,
     &                 ierr)
                  IF (ireturned.NE.4*nlst) THEN
                     WRITE (*,*) 'ERROR - ireturned.NE.4*nlst dvxyzprec'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': received forces ',istatus(MPI_SOURCE)
#endif

                  DO l = 1, nlst
                     ipos = llist(l)
                     jpos = ntot + l
                     DO k = 1, 3
                        dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
                     END DO
                     poten(ipos) = poten(ipos) + dvxyzu(4,jpos)
                  END DO
#ifdef MPIDEBUG
                  print *,iproc,': stored forces '
#endif
               END DO
c               icount = 0
c               DO j = 1, nlst
c                  IF (nneigh(llist(j))*4+3.NE.nneighfull(llist(j))) THEN
c                     icount = icount + 1
c
c                     IF (icount.LT.100) THEN
c                        l = llist(j)
c                     print *,j,l,nneigh(l),nneighfull(l),
c     &                sqrt(xyzmh(1,l)**2+xyzmh(2,l)**2+xyzmh(3,l)**2),
c     &                       xyzmh(5,l)
c                     ENDIF
c                  ENDIF
c               END DO
#ifdef MPIDEBUG
               print *,iproc,': Numbers of neighbours ',nlst,icount
#endif
c
c--Density, divv, nneigh are now complete for this process, do for next one
c     Also know how many and which particles have neighbours on the remote
c     process and need to be sent back to get additional forces
c
c--Other processes do work on the particles being sent
c
            ELSE
c
c--Note that because xyzmh is also used to store the tree structure (which 
c     extends up to ntot*2 + 2) we need to read these in differently than
c     the vxyzu() etc values.
c
               CALL MPI_RECV(xyzmh(1,2*ntot+3),idim,i5REAL8,
     &              MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

               CALL MPI_RECV(vxyzu(1,ntot+1), idim, i4REAL8, iprocrec,
     &              1, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i4REAL8, icheck, ierr)
               IF (icheck.NE.inumber) THEN
                  WRITE (*,*) 'ERROR - icheck.NE.inumber 1 ',iproc,
     &                 iprocrec,icheck,inumber
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got ',inumber,' particles'
#endif

c
c--This should return neighbours lists for the particles as well as 
c     gravitational forces from nodes and leaves
c
               inumreturn = 0
               DO j = 1, inumber
c                  if (itime.NE.0 .AND. MOD(j,1000).EQ.0) 
c     &                 print *,iproc,': treef ',j
                  jpart = j + ntot
                  iphase(jpart) = 0
                  IF (jpart.GT.idim) THEN
                     WRITE (*,*) 'jpart.GT.idim'
                     CALL quit
                  ENDIF
                  CALL treef(jpart,npart,ntot,xyzmh,acc,igphi,
     &                 fsx,fsy,fsz,epot)
                  dvxyzu(1,jpart) = fsx
                  dvxyzu(2,jpart) = fsy
                  dvxyzu(3,jpart) = fsz
                  dvxyzu(4,jpart) = epot
c                   if (iproc.EQ.1 .AND. MOD(j,1000).EQ.0 .AND. 
c     &                 itime.NE.0)
c     &                 print *,iproc,': F ',j,fsx,epot,nneigh(jpart)
c
c--Make a list of indices to return densities
c
                  IF (nneigh(jpart).GE.1) THEN
                     inumreturn = inumreturn + 1
                     llistrec(inumreturn) = jpart
                  ENDIF
               END DO
#ifdef MPIDEBUG
               print *,iproc,': number returned neigh = ',inumreturn
#endif
               nneightogetback(iprocrec+1) = inumreturn
               IF (inumreturn.GT.0) THEN
c
c--The neighbours lists should then be used to calculate contributions to
c     rho and divv (and curlv if required) and add them onto (rho,divv)passed
c     -1 is flag for using passed quantities
c     Only need to do this if there are neighbours!
c
                  CALL densityi(npart,ntot,xyzmh,vxyzu,ekcle,1,
     &                 inumreturn,llistrec,itime)
#ifdef MPIDEBUG
                  print *,iproc,' done density'
#endif
c
c--Need to mark neighbours of remote particles to have their rho etc updated
c
                  IF (itiming) CALL getused(td11)
                  DO l = 1, inumreturn
                     ipart = llistrec(l)
                     DO j = 1, nneigh(ipart)
                        IF (j.GE.nlmax) THEN
                           jpart = neighover(j-nlmax+1,
     &                          ABS(neighb(nlmax,ipart))) 
                        ELSE
                           jpart = neighb(j,ipart)
                        ENDIF
                        IF (iavail(jpart).EQ.0 .AND.
     &                       (.NOT.iscurrent(jpart))) THEN
                           iavail(jpart) = 1
                           nlistavail = nlistavail + 1
                           IF (nlistavail+nlst.GT.idim) THEN
                              WRITE (*,*) 'nlistavail+nlst.GT.idim'
                              WRITE (iprint,*) 'nlistavail+nlst.GT.idim'
                              CALL quit
                           ENDIF
                           llist(nlst+nlistavail) = jpart
                        ENDIF
                     END DO
                  END DO
                  IF (itiming) THEN
                     CALL getused(td12)
                     td1 = td1 + (td12 - td11)
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': set iavail'
#endif
c
c--Return indices
c
#ifdef MPIDEBUG
                  print *,iproc,': modifying llistrec'
#endif
                  DO j = 1, inumreturn
                     llistrec(j) = llistrec(j) - ntot
                  END DO
#ifdef MPIDEBUG
                  print *,iproc,': sending indices to ',iprocrec,
     &                 inumreturn
#endif

                  CALL MPI_SEND(llistrec, inumreturn, MPI_INTEGER,
     &                 iprocrec, 10, MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': returned indices to ',iprocrec
#endif
c
c--Return densities of particles with density contributions (if any)
c
c--Note - need to change values of llistrec because in C array indices
c     start from 0 not 1 and MPI assumes C-type indexing
c
                  DO j = 1, inumreturn
                     llistrec(j) = llistrec(j) + ntot - 1
                  END DO

c                  CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1,
c     &                 llistrec, MPI_REAL4, indexMPI1return, ierr)
                  CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths,
     &                 llistrec, MPI_REAL4, indexMPI1return, ierr)
                  CALL MPI_TYPE_COMMIT(indexMPI1return,ierr)

c                  CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1, 
c     &                 llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
                  CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths, 
     &                 llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
                  CALL MPI_TYPE_COMMIT(indexMPI_I1return,ierr)

                  CALL MPI_SEND(rho, 1, indexMPI1return,
     &                 iprocrec, 11, MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent rho'
#endif

                  CALL MPI_SEND(divv, 1, indexMPI1return,
     &                 iprocrec, 12, MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent divv'
#endif

                  CALL MPI_SEND(nneigh, 1 ,indexMPI_I1return,
     &                 iprocrec, 13, MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent nneigh'
#endif

                  CALL MPI_TYPE_FREE(indexMPI_I1return,ierr)
                  CALL MPI_TYPE_FREE(indexMPI1return,ierr)
                  inumofreturns = inumofreturns + 1
#ifdef MPIDEBUG
                  print *,iproc,': all sent'
#endif
               ELSE
c
c--Return fact that there are no nneighbours
c
                  CALL MPI_SEND(llistrec, 0, MPI_INTEGER, iprocrec,
     &                 10, MPI_COMM_WORLD, ierr)
               ENDIF
c
c--Send gravity forces back to sending process
c
#ifdef MPIDEBUG
               print *,iproc,' sending forces'
#endif
               CALL MPI_SEND(dvxyzu(1,ntot+1),4*inumber,MPI_REAL8,
     &              iprocrec,14,MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
               print *,iproc,' sent forces'
#endif
            ENDIF
         END DO

 4000    CONTINUE

         CALL MPI_TYPE_FREE(indexMPI4,ierr)
         CALL MPI_TYPE_FREE(indexMPI5,ierr)

#ifdef MPIDEBUG
         print *,iproc,': FINISHED PART 1 ',itime
         IF (iproc.EQ.0) print *,iproc,':A',dvxyzu(1,625),
     &        dvxyzu(2,625),dvxyzu(3,625),xyzmh(1,625),
     &        xyzmh(2,625),xyzmh(3,625),
     &        rho(625),divv(625),nneighfull(625)
         IF (iproc.EQ.0) print *,iproc,':AA',dvxyzu(1,1300),
     &        dvxyzu(2,1300),dvxyzu(3,1300),xyzmh(1,1300),
     &        xyzmh(2,1300),xyzmh(3,1300),
     &        rho(1300),divv(1300),nneighfull(1300)
         IF (iproc.EQ.1) print *,iproc,':A',dvxyzu(1,25786),
     &        dvxyzu(2,25786),dvxyzu(3,25786),xyzmh(1,25786),
     &        xyzmh(2,25786),xyzmh(3,25786),
     &        rho(25786),divv(25786),nneighfull(25786)
         IF (iproc.EQ.2) print *,iproc,':A',dvxyzu(1,460),
     &        dvxyzu(2,460),dvxyzu(3,460),xyzmh(1,460),
     &        xyzmh(2,460),xyzmh(3,460),
     &        rho(460),divv(460),nneighfull(460)
         IF (iproc.EQ.3) print *,iproc,':A',dvxyzu(1,25646),
     &        dvxyzu(2,25646),dvxyzu(3,25646),xyzmh(1,25646),
     &        xyzmh(2,25646),xyzmh(3,25646),
     &        rho(25646),divv(25646),nneighfull(25646)
         print *,' '
#endif
#endif

         nlstall = nlst + nlistavail

         IF (itiming) CALL getused(td21)
C$OMP PARALLEL default(none)
C$OMP& shared(npart,it1,imax,itime,it0,isteps,dt,imaxstep,ekcle)
C$OMP& shared(divv,rho,dumrho,vxyzu,pr,vsound,ntot,ireal,ibound,icall)
C$OMP& shared(iphase,nlst,llist,nlistavail,iavail,encal,uradconst)
C$OMP& private(i,j,deltat,deltarho,ipart)
C$OMP DO SCHEDULE (runtime)
         DO i = 1, nlst+nlistavail
            ipart = llist(i)
            IF (i.LE.nlst) THEN
               dumrho(ipart) = rho(ipart)
c
c--Set e(i) for the first time around because density only set here
c     This sets radiation and matter to have same initial temperature
c
               IF (icall.EQ.1 .AND. encal.EQ.'r') THEN
                  IF (ekcle(1,ipart).EQ.0.0) THEN
                     ekcle(3,ipart) = getcv(rho(ipart),vxyzu(4,ipart))
                     ekcle(1,ipart) = uradconst*(vxyzu(4,ipart)/
     &                    ekcle(3,ipart))**4/rho(ipart)
                  ENDIF
               ENDIF
#ifdef MPI
c
c--Need to calculate full pressure, sound speed from FULL density
c
               CALL eospg(ipart,vxyzu,dumrho,pr,vsound,ekcle)
#endif
            ELSE
               iavail(ipart) = 0
            
               IF (ipart.LE.npart) THEN
                  IF (it1(ipart).EQ.imax) THEN
            deltat = dt*(itime - it0(ipart) - isteps(ipart)/2)/imaxstep
                  ELSE
                     deltat = dt*(itime - it0(ipart))/imaxstep
                  ENDIF
c     
c--Update the density value at neighbor's locations
c--Avoid, though, abrupt changes in density
c
                  deltarho = -deltat*divv(ipart)
                  IF (ABS(deltarho).GT.rho(ipart)/2.) THEN
                     deltarho = SIGN(1.0,deltarho)*rho(ipart)/2.0
                  ENDIF
                  dumrho(ipart) = rho(ipart) + deltarho
               
                  CALL eospg(ipart,vxyzu,dumrho,pr,vsound,ekcle)
               ENDIF
            ENDIF
         END DO
C$OMP END DO

C$OMP DO SCHEDULE (runtime)
         DO i = npart + 1, ntot
            j = ireal(i)
            dumrho(i) = dumrho(j)
            pr(i) = pr(j)
            vsound(i) = vsound(j)
            divv(i) = 0.
            vxyzu(4,i) = vxyzu(4,j)
            IF (encal.EQ.'r' .AND. ibound.EQ.100) 
     &           ekcle(1,i) = ekcle(1,j)
         END DO
C$OMP END DO
C$OMP END PARALLEL

         IF (itiming) THEN
            CALL getused(td22)
            td2 = td2 + (td22 - td21)
         ENDIF
c
c--End if for nlst>nptmass
c
c      ENDIF
c
c--Compute implicit radiative transfer
c
      IF (itiming) CALL getused(tass1)

      IF(encal.EQ.'r' .OR. encal.EQ.'m') THEN
         WRITE (*,*) 'Calling RT at realtime ',dt*itime/imaxstep+gt
         IF (encal.EQ.'r') THEN
c           CALL ASS(nlst_in,nlst_end,nlstall,llist,dt,itime,npart,
c     &        xyzmh,vxyzu,ekcle,dumrho,dedxyz)
         ELSE
c            CALL montecarloRT(nlst_in,nlst_end,nlstall,llist,npart,
c     &         xyzmh,vxyzu,dumrho)
         ENDIF

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstall,vxyzu,dumrho,pr,vsound,llist,ekcle)
C$OMP& private(i,ipart)
         DO i = 1, nlstall
            ipart = llist(i)
            CALL eospg(ipart,vxyzu,dumrho,pr,vsound,ekcle)
         END DO
C$OMP END PARALLEL DO
      END IF

      IF (itiming) THEN
        CALL getused(tass2)
        tass = tass + (tass2 - tass1)
      ENDIF

#ifdef MPI
c
c--Now need to do forcei contributions from neighbours on remote nodes.
c     ONLY sends particles known to have neighbours on the remote nodes
c     from the previous send for all nlst particles which calculated
c     gravity forces (except those from neighbours), rho, divv and neighbours.
c
c     NOTE: The TOTAL number of particles being sent to a process from ALL
c     other processes for neighbours calculations MUST BE LESS THAN idim.
c     In addition, the number being sent to any individual process must be
c     less than ineighproc = idim/10.
c     These should be okay since only boarder particles need to be sent.
c
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
      print *,iproc,': START2 ',maxnneighsentback,nneighsentback(1),
     &     nneighsentback(2),nneighsentback(3),nneighsentback(4),
     &     inumofreturns
#endif
      inumbertotal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (maxnneighsentback.GT.0) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done (unlike for
c     the above MPI calls where data is processed as it is received).  This
c     is because here it is assumed that the total number of neighbours on
c     which forces need to be calculated is less than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (nneighsentback(j+1).GT.0) THEN
#ifdef MPIDEBUG
                        print *,iproc,': sending neighbour data to ',j,
     &                       ' starting at ',llistsentback(1,j+1),'+1'
#endif
c
c--Else does not need to send any particle back to this particular process.
c
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c                       1,llistsentback(1,j+1),i5REAL8,indexMPI5,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i5REAL8,indexMPI5,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

                        CALL MPI_SEND(xyzmh,1,indexMPI5,j,20,
     &                       MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
                        print *,iproc,' sent xyzmh to ',j
#endif

                        CALL MPI_SEND(vxyzu,1,indexMPI4,j,21,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent vxyzu'
#endif

                        CALL MPI_SEND(rho,1,indexMPI1,j,22,
     &                       MPI_COMM_WORLD, ierr)
c                  do ll = 1, 10
c                     print *,'S   ',rho(llistsentback(ll,j+1)+1)
c                  end do
c                  do ll = 10, 1, -1
c                     print *,'S   ',rho(llistsentback(
c     &                    nneighsentback(j+1)-ll+1,j+1)+1)
c                  end do
#ifdef MPIDEBUG
                        print *,iproc,' sent rho'
#endif

                        CALL MPI_SEND(divv,1,indexMPI1,j,23,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent divv'
#endif

                        CALL MPI_SEND(pr,1,indexMPI1,j,24,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent pr to ',j
#endif

                        CALL MPI_SEND(vsound,1,indexMPI1,j,25,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent vsound to ',j
#endif

                        CALL MPI_TYPE_FREE(indexMPI1,ierr)
                        CALL MPI_TYPE_FREE(indexMPI4,ierr)
                        CALL MPI_TYPE_FREE(indexMPI5,ierr)
#ifdef MPIDEBUG
                        print *,iproc,': sent neighbour data to ',j
#endif
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned.LT.inumofreturns .AND. 
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
               print *,iproc,': receiving neigh data ',inumberreturned,
     &              inumofreturns
#endif
               istart = ntot + inumbertotal + 1
               CALL MPI_RECV(xyzmh(1,istart+ntot+2), idim,i5REAL8,
     &              MPI_ANY_SOURCE, 20, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
      IF (iproc.EQ.1 .AND. iprocrec.EQ.0) 
     &     print *,iproc,':E',xyzmh(1,istart+ntot+2+625-1),
     &     xyzmh(2,istart+ntot+2+625-1),xyzmh(3,istart+ntot+2+625-1),
     &     istart+ntot+2+625-1

      IF (iproc.EQ.2 .AND. iprocrec.EQ.3) 
     &     print *,iproc,':E',xyzmh(1,istart+ntot+2+2054-1),
     &     xyzmh(2,istart+ntot+2+2054-1),xyzmh(3,istart+ntot+2+2054-1),
     &     istart+ntot+2+2054-1

               print *,iproc,': got xyzmh from ',iprocrec,inumber,
     &     ' to be put into ',istart
#endif

               IF (inumber.GT.0) THEN
                  CALL MPI_RECV(vxyzu(1,istart),inumber,i4REAL8,
     &                 iprocrec, 21, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i4REAL8, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 1 ',
     &                    iproc
                     CALL quit
                  ENDIF

                  CALL MPI_RECV(dumrho(istart),inumber,MPI_REAL4,
     &                 iprocrec, 22, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 2 ',
     &                    iproc
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  IF (iproc.EQ.1 .and. iprocrec.EQ.0) 
     &                 print *,iproc,': ??? ',iprocrec,
     &                 dumrho(istart+1300-1),xyzmh(1,istart+ntot+2+
     &                 1300-1),xyzmh(2,istart+ntot+2+
     &                 1300-1),xyzmh(3,istart+ntot+2+
     &                 1300-1),istart+1300-1
c                  do ll = 1, 10
c                     print *,'R   ',dumrho(istart+ll-1)
c                  end do
c                  do ll = 10, 1, -1
c                     print *,'R   ',dumrho(istart+inumber-ll)
c                  end do
                  print *,iproc,': got rho from ',iprocrec,icheck
#endif

                  CALL MPI_RECV(divv(istart),inumber,MPI_REAL4,iprocrec,
     &                 23, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 3 ',
     &                    iproc
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': got divv from ',iprocrec,icheck
#endif

                  CALL MPI_RECV(pr(istart), inumber, MPI_REAL4,iprocrec,
     &                 24, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 4 ',
     &                    iproc
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': got pr from ',iprocrec,icheck
#endif

                CALL MPI_RECV(vsound(istart),inumber,MPI_REAL4,iprocrec,
     &                 25, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 5 ',
     &                    iproc
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': got vsound from ',iprocrec,icheck
#endif
               ELSE
                  WRITE (*,*) 'ERROR - number returned <=0: ',inumber
                  CALL quit
               ENDIF
               inumbertotal = inumbertotal + inumber
               inumberindiv(inumberreturned) = inumber
               inumbercumm(inumberreturned) = istart
               inumberproc(inumberreturned) = iprocrec
#ifdef MPIDEBUG
               print *,iproc,' set numbers '
#endif
            ENDIF
         ENDIF
      END DO
      IF (inumbertotal.GT.idim) THEN
         WRITE (*,*) 'ERROR - inumbertotal.GT.idim ', iproc
         CALL quit
      ENDIF
      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR - inumberreturned.NE.inumofreturns ', iproc
         CALL quit
      ENDIF
c     
c--Calculate ONLY neighbours lists for remote particles.  Also need to ensure
c     that forces and poten are zeroed before call of forcei.
c
      DO j = 1, inumbertotal
         jpart = ntot + j
         llist(nlst_end + j) = jpart
         iphase(jpart) = 0
         CALL treef(jpart,npart,ntot,xyzmh,acc,0,fsx,fsy,fsz,epot)
#ifdef MPIDEBUG
c         IF (iproc.EQ.1 .AND. jpart.EQ.33629) print *,iproc,': NN ',
c     &        dumrho(jpart),nneigh(jpart),j,jpart,xyzmh(1,jpart+
c     &        ntot+2),xyzmh(2,jpart+ntot+2),xyzmh(3,jpart+ntot+2)
         IF (iproc.EQ.1 .AND. jpart.EQ.33630) print *,iproc,': NN ',
     &        dumrho(jpart),nneigh(jpart),j,jpart,xyzmh(1,jpart+
     &        ntot+2),xyzmh(2,jpart+ntot+2),xyzmh(3,jpart+ntot+2)
c         IF (iproc.EQ.1 .AND. jpart.EQ.33631) print *,iproc,': NN ',
c     &        dumrho(jpart),nneigh(jpart),j,jpart,xyzmh(1,jpart+
c     &        ntot+2),xyzmh(2,jpart+ntot+2),xyzmh(3,jpart+ntot+2)
#endif
         DO k = 1, 3
            dvxyzu(k,jpart) = 0.0
         END DO
         poten(jpart) = 0.0
      END DO
      nlst_tot = nlst_end + inumbertotal

#ifdef MPIDEBUG
      print *,iproc,': Received all particles with remote neighbours'
      print *,' '
      print *,' '
#endif
      istart = ntot + 1 + ntot + 2
#else
      nlst_tot = nlst_end
#endif
c
c--Compute forces on EACH particle
c
      IF (itiming) CALL getused(tforce1)

#ifdef MPIDEBUG
      IF (iproc.EQ.0) print *,iproc,': NLOCAL ',nneigh(1300)
#endif

      CALL forcei(nlst_in,nlst_end,nlst_tot,llist,dt,itime,npart,ntot,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,dedxyz)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

#ifdef MPI




c            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)





#ifdef MPIDEBUG
         IF (iproc.EQ.2) print *,iproc,':G',dvxyzu(1,ntot + 1+2054-1),
     &              dvxyzu(2,ntot+1+2054-1),dvxyzu(3,ntot+1+2054-1),
     &     ntot+1+2054-1,ntot
#endif
c
c--Now need to return forcei contributions from neighbours on remote nodes.
c
      istart = ntot + 1
      istartrec = istart + inumbertotal
      IF (istartrec + maxnneighsentback.GT.2*idim) THEN
         WRITE (*,*) 'ERROR - istartrec + maxnneighsentback.GT.2*idim'
         CALL quit
      ENDIF
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back forces, du, potential energy to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' forces ',j
               IF (iproc.EQ.14) THEN
                  DO iii = 0, numproc -1
                     print *,iproc,': LIST ',iii,nneighsentback(iii+1)
                  END DO
               ENDIF
#endif
               CALL MPI_RECV(dvxyzu(1,istartrec),4*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,30,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.4*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.4*nnsentback dvxyzu'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
c               print *,iproc,': got forces from ',iprocrec,' put into ',
c     &              llistsentback(1,iprocrec+1)+1
#endif

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1                  
                  DO k = 1, 4
                     dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
                  END DO
               END DO

               CALL MPI_RECV(poten(istartrec),maxnneighsentback,
     &              MPI_REAL4,iprocrec,31,MPI_COMM_WORLD,istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsendback poten'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got poten from ',iprocrec
#endif

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  poten(ipos) = poten(ipos) + poten(jpos)
               END DO
            END DO
c
c--Other processes send the forces back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                  print *,iproc,': sending forces to ',iprocsend,
     &                 ' istart ',istart,' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(dvxyzu(1,istart),4*inumber,MPI_REAL8,
     &                 iprocsend,30,MPI_COMM_WORLD, ierr)
                  CALL MPI_SEND(poten(istart),inumber,MPI_REAL4,
     &                 iprocsend,31,MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent forces and poten to ',iprocsend
#endif
                  GOTO 7765
               ENDIF
            END DO
 7765       CONTINUE

            GOTO 8888

            IF (inumberreturned.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned = inumberreturned + 1
               istart = inumbercumm(inumberreturned)
               inumber = inumberindiv(inumberreturned)
               iprocsend = inumberproc(inumberreturned)
#ifdef MPIDEBUG
               print *,iproc,': sending forces to ',iprocsend,
     &              ' istart ',istart,' returned ',inumberreturned,
     &              ' of ',inumofreturns
#endif
               CALL MPI_SEND(dvxyzu(1,istart),4*inumber,MPI_REAL8,
     &              iprocsend,30,MPI_COMM_WORLD, ierr)
               CALL MPI_SEND(poten(istart),inumber,MPI_REAL4,
     &              iprocsend,31,MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
               print *,iproc,': sent forces and poten to ',iprocsend
c         IF (iproc.EQ.2) print *,iproc,':D',dvxyzu(1,istart+2054-1),
c     &              dvxyzu(2,istart+2054-1),dvxyzu(3,istart+2054-1)
#endif
               istart = istart + inumber
            ENDIF

 8888       CONTINUE

         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND 
c     because none will be expected)
c
      END DO

      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR - inumberreturned.NE.inumofreturns ', iproc
         CALL quit
      ENDIF

      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)

#ifdef MPIDEBUG
      print *,iproc,': FINISHED ',itime
         IF (iproc.EQ.0) print *,iproc,':B',dvxyzu(1,625),
     &        dvxyzu(2,625),dvxyzu(3,625),xyzmh(1,625),
     &        xyzmh(2,625),xyzmh(3,625)
         IF (iproc.EQ.1) print *,iproc,':B',dvxyzu(1,25786),
     &        dvxyzu(2,25786),dvxyzu(3,25786),xyzmh(1,25786),
     &        xyzmh(2,25786),xyzmh(3,25786)
         IF (iproc.EQ.2) print *,iproc,':B',dvxyzu(1,460),
     &        dvxyzu(2,460),dvxyzu(3,460),xyzmh(1,460),
     &        xyzmh(2,460),xyzmh(3,460)
         IF (iproc.EQ.3) print *,iproc,':B',dvxyzu(1,25646),
     &        dvxyzu(2,25646),dvxyzu(3,25646),xyzmh(1,25646),
     &        xyzmh(2,25646),xyzmh(3,25646)
      print *,' '
#endif
#endif

      RETURN
      END
