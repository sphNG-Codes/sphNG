      SUBROUTINE dust_implicit (dt,npart,ntot,xyzmh,vxyzu,
     &            nlst_in,nlst_end,itime)
c************************************************************
c                                                           *
c  Subroutine to compute the density and smoothing lengths  *
c     self-consistently using iteration if necessary.       *
c     It also calculates the velocity divergence, the curl, *
c     and the gravity softening term with variable h, and   *
c     interpolates the density for particles which are      *
c     neighbours of particles that are currently in the     *
c     list.  This subroutine uses the binary tree algorithm *
c     to locate neighbours.  Neighbours are not stored in   *
c     a list.                                               *
c                                                           *
c     Code written by MRB and DJP (14/12/2005).             *
c     MPI version written by MRB (26/06/2007).
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2)

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/dustimplicit'

#ifdef MPICOPY
      INCLUDE 'COMMONS/mpicopy'
#endif

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine dust_implicit')
c
c--Initialise
c
#ifdef MPIALL
#ifdef MPIDEBUG
      print *,iproc,': Entered dust_implicit ',itime,nlst_end
#endif
#else
c      print *,': Entered dust_implicit ',itime
#endif

#ifdef MPI
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)

      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
#endif
c
c--SEMI-IMPLICIT DUST DRAG needs to be done before the velocities
c     and internal energies are needed to calculate other quantities.
c     The velocities and internal energies are modified.
c     Need to:
c        1) Calculate norm() and rhoG() for dust particles.  The norm()
c           must be calculated after densities and smoothing lengths
c           are set because it needs to know the gas particle densities.
c        2) They will need to be circulated to other MPI processes
c        3) Then need to calculate new vxyzu() looping over dust and
c           changing its velocity and those gas particles that it
c           interacts with.  Need to calculate both local and foreign
c           contributions to changes.
c
c--Calculate local contributions to dustnorm and rhoG
c
      CALL dust_norm_rhoG(nlst_in,nlst_end,npart,ntot,
     &     xyzmh,SIZE(xyzmh,2))
c
#ifdef MPI
c
c--For MPI job, need to add contributions to dustnorm and rhoG
c     from other MPI processes
c
c      print *,'MPI: dust_norm_rhoG ',ncompactlocal,inumbertotal,
c     &     ncompact,nlst_end

      CALL dust_norm_rhoG(ncompactlocal+1,ncompactlocal+inumbertotal,
     &     npart,ntot,xyzmh,SIZE(xyzmh,2))
c
c--Now need to return contributions from neighbours on remote nodes.
c
      istart = ntot + 1
      istartrec = istart + inumbertotal
      IF (istartrec + maxnneighsentback.GT.2*idim) THEN
         WRITE (*,*) 'ERROR - istartrec + maxnneighsentback.GT.2*idim'
         CALL quit(1)
      ENDIF
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back dustnorm & rhoG to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' dustnorm contributions ',j,
     &              ' with maxnneighsentback ',maxnneighsentback
#endif
               CALL MPI_RECV(dustnorm(istartrec),maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,70,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               IF (istartrec+ireturned.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istartrec+ireturned.GT.idim2 ',
     &                 istartrec,ireturned
                  CALL quit(1)
               ENDIF
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsentback dust ',
     &                 ireturned,nneighsentback(iprocrec+1)
                  CALL quit(1)
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got dustnorm from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
               print *,iproc,': received values ',
     &              (dustnorm(j3),j3=istartrec,
     &              istartrec+1)
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,dustnorm)
C$OMP& shared(istartrec,iphase)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1

                  IF (iphase(ipos).LT.11 .AND. 
     &                 dustnorm(ipos) + dustnorm(jpos).GT.0.) THEN
                     WRITE (*,*) 'ERROR - non-dust has non-zero',
     &                    ' dustnorm ',ipos,jpos,dustnorm(ipos),
     &                    dustnorm(jpos)
                     CALL quit(0)
                  ENDIF

                  dustnorm(ipos) = dustnorm(ipos) + dustnorm(jpos)
               END DO
C$OMP END PARALLEL DO

               CALL MPI_RECV(rhoG(istartrec),maxnneighsentback,
     &              MPI_REAL8,iprocrec,71,MPI_COMM_WORLD,istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsendback rhoG'
                  CALL quit(1)
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got rhoG from ',iprocrec
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,rhoG,istartrec)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  rhoG(ipos) = rhoG(ipos) + rhoG(jpos)
               END DO
C$OMP END PARALLEL DO
            END DO
c
c--Other processes send the dustnorm & rhoG values back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                   print *,iproc,': sending dust to ',iprocsend,
     &                 ' istart ',istart,' number ',inumber,
     &                 ' returned ',inumberreturned,
     &                 ' of ',inumofreturns
                   print *,iproc,': sending values ',
     &                  (dustnorm(j2),j2=istart,istart+1)
#endif
                  CALL MPI_SEND(dustnorm(istart),inumber,MPI_REAL8,
     &                 iprocsend,70,MPI_COMM_WORLD, ierr)
                  CALL MPI_SEND(rhoG(istart),inumber,MPI_REAL8,
     &                 iprocsend,71,MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
                       print *,iproc,': sent dustnorm,rhoG to ',
     &                 iprocsend,' and has sent ',inumber,' values'
#endif
                  GOTO 8765
               ENDIF
            END DO
 8765       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO

      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR2 - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit(1)
      ENDIF
c
c--Now need to transfer resulting dustnorm() & isteps() values for dust
c     particles that are interacting with gas particles on other MPI
c     processes!
c
      inumbertotal_2 = 0
      inumberreturned_2 = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done (unlike for
c     the above MPI calls where data is processed as it is received).  This
c     is because here it is assumed that the total number of neighbours for
c     which other quantities need to be calculated is less than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (nneighsentany(j+1)) THEN
#ifdef MPIDEBUG
                       print *,iproc,': sending neighbour data to ',j,
     &                       ' starting at ',llistsentback(1,j+1),'+1',
     &                       ' and sending ',nneighsentback(j+1),
     &                       ' elements'
#endif
c
c--Else does not need to send any particle back to this particular process.
c
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL8,indexMPI1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL8,indexMPI1,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

                           CALL MPI_SEND(dustnorm,1,indexMPI1,j,73,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent dustnorm'
#endif
                           CALL MPI_TYPE_FREE(indexMPI1,ierr)

                        IF (nneighsentback(j+1).GT.0) THEN

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_INTEGER4,
c     &                  indexMPI_INT1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                          llistsentback(1,j+1),MPI_INTEGER4,
     &                          indexMPI_INT1,ierr)
                           CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)

                           CALL MPI_SEND(isteps,1,indexMPI_INT1,j,74,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent isteps'
#endif
                           CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)

#ifdef MPIDEBUG
                       print *,iproc,': sent neighbour data to ',j,
     &                       radkernel,xyzmh(1,1)
#endif
                        ENDIF
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned_2.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned_2 = inumberreturned_2 + 1
#ifdef MPIDEBUG
            print *,iproc,': receiving neigh data ',inumberreturned_2,
     &              inumofreturns
#endif
               istart = ntot + inumbertotal_2 + 1
               CALL MPI_RECV(dustnorm(istart),idim,MPI_REAL8,
     &              i, 73, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               IF (istart+inumber.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istart+inumber.GT.idim2'
                  CALL quit(1)
               ENDIF
               iprocrec = istatus(MPI_SOURCE)

               DO iii2 = istart, istart+inumber-1
                  IF (dustnorm(iii2).GT.0. .AND. iphase(iii2).LT.11)
     &                 THEN
                     print *,'ERROR - xxx ',iii2,dustnorm(iii2),
     &                    iphase(iii2)
                     CALL quit(0)
                  ENDIF
               END DO

#ifdef MPIDEBUG
               print *,iproc,': got dustnorm from ',iprocrec,inumber,
     &     ' to be put into ',istart
#endif
               IF (inumber.GT.0) THEN
                  CALL MPI_RECV(isteps(istart),inumber,MPI_INTEGER4,
     &                 iprocrec, 74, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER4,icheck,ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber s1 ',
     &                    iproc
                     CALL quit(1)
                  ENDIF

               ENDIF
               inumbertotal_2 = inumbertotal_2 + inumber
#ifdef MPIDEBUG
               print *,iproc,' set numbers ',radkernel,xyzmh(1,1)
#endif
            ENDIF
         ENDIF
      END DO
#ifdef MPIDEBUG
      print *,iproc,' inumbertotal_2 is ',inumbertotal_2
#endif
      IF (inumbertotal_2.NE.inumbertotal) THEN
         WRITE (*,*) 'ERROR - inumbertotal_2.NE.inumbertotal '
         CALL quit(1)
      ENDIF
      IF (inumberreturned_2.NE.inumberreturned) THEN
         WRITE (*,*) 'ERROR1 - inumberreturned_2.NE.inumberreturned ',
     &        iproc,inumberreturned_2,inumberreturned,inumofreturns
         CALL quit(1)
      ENDIF

#ifdef MPIDEBUG
      print *,iproc,': Received all particles with remote neighbours ',
     &     inumbertotal_2
      print *,' '
      print *,' '
#endif
c
c--END -- Have transferred quantities to other MPI processes
c
#endif
c
c--Calculate new velocities due to drag interaction between dust & gas
c
      DO n = 1, ncompact
         ipart = ivar(3,n)
         DO k = 1, 4
            vxyzu_change(k,ipart) = 0.
         END DO
      END DO

      CALL dust_drag(nlst_in,nlst_end,npart,ntot,
     &     xyzmh,SIZE(xyzmh,2),vxyzu,SIZE(vxyzu,2),itime,dt)
#ifdef MPI
c
c--For MPI job, need to add contributions to dustnorm and rhoG
c     from other MPI processes
c
c      print *,'MPI: dust_drag ',ncompactlocal,inumbertotal,
c     &     ncompact

      DO n = ncompactlocal+1,ncompactlocal+inumbertotal
         ipart = ivar(3,n)
         DO k = 1, 4
            vxyzu_change(k,ipart) = 0.
         END DO
      END DO

      CALL dust_drag(ncompactlocal+1,
     &     ncompactlocal+inumbertotal,
     &     npart,ntot,xyzmh,SIZE(xyzmh,2),vxyzu,SIZE(vxyzu,2),
     &     itime,dt)
c
c--Now need to return contributions from neighbours on remote nodes.
c
      istart = ntot + 1
      istartrec = istart + inumbertotal
      IF (istartrec + maxnneighsentback.GT.2*idim) THEN
         WRITE (*,*) 'ERROR - istartrec + maxnneighsentback.GT.2*idim'
         CALL quit(1)
      ENDIF
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back vxyzu_change() to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' vxyzu_change contributions ',j,
     &              ' with maxnneighsentback ',maxnneighsentback
#endif
               CALL MPI_RECV(vxyzu_change(1,istartrec),
     &              maxnneighsentback,
     &              i4REAL8,MPI_ANY_SOURCE,72,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, i4REAL8, ireturned, ierr)
               IF (istartrec+ireturned.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istartrec+ireturned.GT.idim2 ',
     &                 istartrec,ireturned
                  CALL quit(1)
               ENDIF
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsentback dust ',
     &                 ireturned,nneighsentback(iprocrec+1)
                  CALL quit(1)
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got vxyzu_change from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,vxyzu_change)
C$OMP& shared(istartrec)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  DO k = 1, 4
                     vxyzu_change(k,ipos) = vxyzu_change(k,ipos) +
     &                    vxyzu_change(k,jpos)
                  END DO
               END DO
C$OMP END PARALLEL DO

            END DO
c
c--Other processes send the dustnorm & rhoG values back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                   print *,iproc,': sending vxyzu_change to ',iprocsend,
     &                 ' istart ',istart,' number ',inumber,
     &                 ' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(vxyzu_change(1,istart),inumber,i4REAL8,
     &                 iprocsend,72,MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
                       print *,iproc,': sent vxyzu_change to ',
     &                 iprocsend,' and has sent ',inumber,' values'
#endif
                  GOTO 9765
               ENDIF
            END DO
 9765       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO

      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR2 - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit(1)
      ENDIF
#endif
c
c--Now update vxyzu for full changes from dust-gas drag
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,ivar,vxyzu,vxyzu_change)
C$OMP& private(n,ipart,j)
      DO n = nlst_in, nlst_end
         ipart = ivar(3,n)
         DO j = 1, 4
            vxyzu(j,ipart) = vxyzu(j,ipart) + vxyzu_change(j,ipart)
         END DO
      END DO
C$OMP END PARALLEL DO

#ifdef MPI
c
c--Now need to send updated vxyzu to other MPI processes for 
c     particles that interact with particles on other nodes.
c
      inumbertotal_2 = 0
      inumberreturned_2 = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done (unlike for
c     the above MPI calls where data is processed as it is received).  This
c     is because here it is assumed that the total number of neighbours for
c     which other quantities need to be calculated is less than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (nneighsentany(j+1)) THEN
#ifdef MPIDEBUG
                       print *,iproc,': sending neighbour data to ',j,
     &                       ' starting at ',llistsentback(1,j+1),'+1',
     &                       ' and sending ',nneighsentback(j+1),
     &                       ' elements'
#endif
c
c--Else does not need to send any particle back to this particular process.
c
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

                           CALL MPI_SEND(vxyzu,1,indexMPI4,j,75,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent vxyzu'
#endif
                           CALL MPI_TYPE_FREE(indexMPI4,ierr)
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned_2.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned_2 = inumberreturned_2 + 1
#ifdef MPIDEBUG
            print *,iproc,': receiving neigh data ',inumberreturned_2,
     &              inumofreturns
#endif
               istart = ntot + inumbertotal_2 + 1
               CALL MPI_RECV(vxyzu(1,istart),idim,i4REAL8,
     &              i, 75, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
               IF (istart+inumber.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istart+inumber.GT.idim2 ',
     &                 istart,inumber
                  CALL quit(1)
               ENDIF
               iprocrec = istatus(MPI_SOURCE)

               inumbertotal_2 = inumbertotal_2 + inumber
#ifdef MPIDEBUG
               print *,iproc,' set numbers ',radkernel,xyzmh(1,1)
#endif
            ENDIF
         ENDIF
      END DO
#ifdef MPIDEBUG
      print *,iproc,' inumbertotal_2 is ',inumbertotal_2
#endif
      IF (inumbertotal_2.NE.inumbertotal) THEN
         WRITE (*,*) 'ERROR - inumbertotal_2.NE.inumbertotal ', iproc,
     &        inumbertotal_2,inumbertotal,idim,(inumberindiv(ix),ix=1,8)
         CALL quit(1)
      ENDIF
      IF (inumberreturned_2.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR1 - inumberreturned_2.NE.inumofreturns ',
     &        iproc,inumberreturned_2,inumofreturns
         CALL quit(1)
      ENDIF

#ifdef MPIDEBUG
      print *,iproc,': Received all particles with remote neighbours ',
     &     inumbertotal_2
      print *,' '
      print *,' '
#endif
#endif
c
c--END MPI SECTION that transfers new values of vxyzu()
c

#ifdef MPI
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
#endif

      IF (itrace.EQ.'all') WRITE (iprint,300)
 300  FORMAT ('exit subroutine dust implicit')

c      print *,'Exited dust implicit ',itime

      RETURN
      END

c=============================================================================
