      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &     dvxyzu,dha,npart,ntot,ireal,alphaMM,ekcle,Bevolxyz,dBevolxyz)
c************************************************************
c                                                           *
c  This subroutine drives the computation of the forces on  *
c     every particle on the list.                           *
c                                                           *
c  MRB 14/12/2005:                                          *
c                                                           *
c  The non-gradh version calculates neighbours and gravity  *
c     with neighbours stored in a list, then computes a     *
c     list of particles whose values are required to be     *
c     interpolated (density, pressure, sound speed) and     *
c     does the interpolation, then updates ghosts, the      *
c     does implicit radiative transfer, and finally calls   *
c     forcei to calculate the rest of the forces.           *
c                                                           *
c  In this grad-h version of the code, the structure is     *
c     very different because of the interations required to *
c     set density and the fact that enormous numbers of     *
c     neighbours can be obtained to make density and h      *
c     consistent (overflowing any global neighbour store).  *
c     Neighbours are no longer stored. Rather, they are got *
c     each time for density and forcei. Note that this will *
c     not work for radiative transfer because of the        *
c     iterative solve - would need to calculate neighbours  *
c     each iteration !!!!  In the grad-h version of the     *
c     code, density and neighbours (h's) are calculated     *
c     simulataneously in densityiterate_gradh.f.  This code *
c     also calculates divv() and other quantities as in the *
c     normal density.f code.  It also interpolates the      *
c     density, pressure and sound speed of neighbouring     *
c     particles that are required for forces at the same    *
c     time, rather than constructing a list and doing it    *
c     later.  Thus, it combines density.f and the middle    *
c     part of derivi.  Forcei_gradh then calculates all the *
c     forces, including all gravity forces when it re-finds *
c     the neighbours of each particle.  Thus, it combines   *
c     the first part of derivi.f and forcei.f .             *
c                                                           *
c  MHD note: (DJP 6.1.06)                                   *
c   The MHD quantities passed to this routine               *
c     are the *evolved* MHD variables. These could be       *
c     B, B/rho (usual option) or the Euler potentials.      *
c     However, we send just the magnetic field into the     *
c     force routines. The dBevolxyz returned by derivi      *
c     is the derivative required for evolving the magnetic  *
c     field.                                                *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif

      DIMENSION xyzmh(5,mmax),vxyzu(4,idim),dvxyzu(4,idim)
      REAL*4 dha(1+isizealphaMM,idim),alphaMM(isizealphaMM,idim)
      DIMENSION ireal(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevolxyz(3,imhd),dBevolxyz(3,imhd)

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
c     Bxyz is stored here for calculation of energy & writing to dump file
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/vsmooth'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/compact'

      CHARACTER*7 where
      DIMENSION dedxyz(3,iradtrans),listex(idim)

      DATA where/'derivi'/
c
c--Allow for tracing flow
c
c      IF (gt.GT.50.) print *,iproc,': ENTERED derivi'

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine derivi_gradh')
      IF (igrape.NE.0) THEN
         WRITE (iprint,*) 'ERROR: derivi_P_gradh must have igrape.EQ.0'
         CALL quit
      ENDIF
      IF (nlmax.NE.1) THEN
         WRITE (iprint,*) 'ERROR: derivi_P_gradh must have nlmax.EQ.1'
         CALL quit
      ENDIF
c     
c--Set constants first time around
c
      uradconst = radconst/uergcc
      nlst_in = 1
      nlst_end = nlst

      IF (itrace.EQ.'all') WRITE (iprint, 99002) nlst_in, nlst_end
99002 FORMAT(' derivi_gradh ',I8, I8)
c
c--Find self-consistent density and smoothing length for all the particles
c     in the list.  Only calculate density etc for non-sinks.
c     for Euler potentials, this also means get B from grad alpha x grad beta
c

c
c--For MPI may need to call even if nothing but sinks being done because
c     other processes may be doing gas particles
c
#ifdef MPI
#else
c      IF (nlst_end.GT.nptmass) THEN
#endif
         IF (itiming) CALL getused(tdens1)

#ifdef MPIDEBUG
         print *,'DERIVI: calling densityiterate_gradh'
#endif

         CALL densityiterate_gradh(dt,npart,ntot,xyzmh,vxyzu,dvxyzu,
     &        nlst_in,nlst_end,llist,itime,ekcle,Bevolxyz,Bxyz,
     &        dBevolxyz)

#ifdef MPIDEBUG
         print *,'DERIVI: called densityiterate_gradh'
#endif

         IF (itiming) THEN
            CALL getused(tdens2)
            tdens = tdens + (tdens2 - tdens1)
         ENDIF

c         IF (encal.EQ.'v' .OR. encal.EQ.'i') THEN
c            CALL montecarloRT(npart, xyzmh, vxyzu, dumrho, 
c     &                        dt*itime/imaxstep+gt)
c         END IF

#ifdef MPI
#else
c      ENDIF
#endif
c
c--Implicit hyperdiffusion of div B
c
      IF (itiming) CALL getused(tass1)

c      IF(imhd.EQ.idim) THEN
      IF(.FALSE.) THEN
         WRITE (*,*) 'Calling Hyper at realtime ',dt*itime/imaxstep+gt
c         CALL divBdiffuse(dt,nlst_in,nlst_end,npart,llist,
c     &        xyzmh,dumrho,Bxyz)

         IF (varmhd.EQ.'Bvol') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,Bevolxyz,Bxyz,llist)
C$OMP& private(i,ipart)
            DO i = nlst_in,nlst_end
               ipart = llist(i)
               Bevolxyz(1,ipart) = Bxyz(1,ipart)
               Bevolxyz(2,ipart) = Bxyz(2,ipart)
               Bevolxyz(3,ipart) = Bxyz(3,ipart)
            END DO
C$OMP END PARALLEL DO
         ELSEIF (varmhd.EQ.'Brho') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,Bevolxyz,Bxyz,llist,dumrho)
C$OMP& private(i,ipart)
            DO i = nlst_in,nlst_end
               ipart = llist(i)
               Bevolxyz(1,ipart) = Bxyz(1,ipart)/dumrho(ipart)
               Bevolxyz(2,ipart) = Bxyz(2,ipart)/dumrho(ipart)
               Bevolxyz(3,ipart) = Bxyz(3,ipart)/dumrho(ipart)
            END DO
C$OMP END PARALLEL DO
         ENDIF
      END IF
c
c--div B projection
c
c      IF (imhd.EQ.idim) THEN
c         IF (varmhd.EQ.'Brho' .OR. varmhd.EQ.'Bvol') THEN
c            divBmax = 0.
c            DO i=1,npart
c               divBmax = max(divcurlB(1,i),divBmax)
c            ENDDO
c            WRITE(iprint,*) 'div B max = ',divBmax
c            CALL divBclean(nlst_in,nlst_end,npart,ntot,llist,
c     &                     xyzmh,rho,Bevolxyz)
c         ENDIF
c      ENDIF
      
      IF (itiming) THEN
        CALL getused(tass2)
        tass = tass + (tass2 - tass1)
      ENDIF

#ifdef MPI
c
c--Now need to do forcei contributions from neighbours on remote nodes.
c     ONLY sends particles known to have neighbours on the remote nodes
c     from the previous send (in densityiterate_gradh) for all nlst 
c     particles (excluding sinks) which calculated gravity forces (except 
c     those from neighbours), rho, divv and number of neighbours (list of
c     neighbours is not stored in gradh code).
c
c     NOTE: The TOTAL number of particles being sent to a process from ALL
c     other processes for neighbours calculations MUST BE LESS THAN idim.
c     In addition, the number being sent to any individual process must be
c     less than ineighproc = idim/10.
c     These limits should be okay since only boarder particles need to be sent.
c
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(isizealphaMM,MPI_REAL4,ialphaREAL4,ierr)

      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)
      CALL MPI_TYPE_COMMIT(ialphaREAL4,ierr)

#ifdef MPIDEBUG
      print *,iproc,': START2 ',maxnneighsentback,nneighsentback(1),
     &     nneighsentback(2),nneighsentback(3),nneighsentback(4),
     &     inumofreturns
#endif
      inumbertotal = 0
      inumberreturned = 0
      inumofsends = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done (unlike for
c     the above MPI calls where data is processed as it is received).  This
c     is because here it is assumed that the total number of remote particles
c     on which forces need to be calculated by this MPI process is less 
c     than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
c
c--Only sends anything at all to another MPI process if it needs to (the
c     other MPI process knows when it doesn't need to expect anything)
c
                     IF (nneighsentany(j+1)) THEN
                        inumofsends = inumofsends + 1
#ifdef MPIDEBUG
                        print *,iproc,': sending neighbour data to ',j,
     &                       ' starting at ',llistsentback(1,j+1),'+1'
#endif
c
c--Else does not need to send any particle back to this particular process.
c
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

c                        IF (iproc.EQ.0) THEN
c                           print *,' SentbackDER ',llistsentback(1,j+1),
c     &                       iunique(iorig(llistsentback(1,j+1)+1)),
c     &                       i,j,iunique(iorig(llistsentback(2,j+1)+1))
c
c                           DO ijk = 1, nneighsentback(j+1)
c                              print *,iorig(llistsentback(ijk,j+1)+1)-1
c                           END DO
c                        ENDIF

                        DO ijk = 1, nneighsentback(j+1)
                        listex(ijk) = iorig(llistsentback(ijk,j+1)+1)-1
                        END DO

                        CALL MPI_SEND(divv,1,indexMPI1,j,23,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent divv'
#endif
                        IF (nneighsentback(j+1).GT.0) THEN

                        IF (imhd.EQ.idim) THEN
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)
                        ENDIF

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i2REAL4,indexMPI2,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i2REAL4,indexMPI2,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI2,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),ialphaREAL4,
c     &                  indexMPIalpha,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),ialphaREAL4,
     &                       indexMPIalpha,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIalpha,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_INTEGER8,indexMPI_I8,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &    listex,MPI_INTEGER8,indexMPI_I8,ierr)
c     &    iorig(llistsentback(1,j+1)+1),MPI_INTEGER8,indexMPI_I8,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI_I8,ierr)

                        CALL MPI_SEND(iunique,1,indexMPI_I8,j,200,
     &                       MPI_COMM_WORLD, ierr)

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

                        CALL MPI_SEND(gradhs,1,indexMPI2,j,26,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent gradhs to ',j
#endif

                        IF (imhd.EQ.idim) THEN
                           CALL MPI_SEND(Bxyz,1,indexMPI3,j,27,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent Bxyz to ',j
#endif

                           CALL MPI_SEND(Bevolxyz,1,indexMPI3,j,28,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent Bevolxyz to ',j
#endif

                           CALL MPI_TYPE_FREE(indexMPI3,ierr)
                        ENDIF

                        IF (ifsvi.EQ.6 .OR. imhd.EQ.idim) THEN
                           CALL MPI_SEND(alphaMM,1,indexMPIalpha,j,29,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent alphaMM to ',j
#endif                           
                        ENDIF

                        CALL MPI_TYPE_FREE(indexMPI_I8,ierr)
                        CALL MPI_TYPE_FREE(indexMPIalpha,ierr)
                        CALL MPI_TYPE_FREE(indexMPI2,ierr)
                        CALL MPI_TYPE_FREE(indexMPI1,ierr)
                        ENDIF

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
c
c--Need to store this particles beyond the end of the tree nodes, hence the
c     +ntot+2 in the line below (this only applies to xyzmh, not the other
c     arrays such as vxyzu).
c
               CALL MPI_RECV(divv(istart), idim, MPI_REAL4,
     &              i, 23, MPI_COMM_WORLD, istatus, ierr)
c     &              MPI_ANY_SOURCE, 23, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
                  print *,iproc,': got divv from ',iprocrec,inumber
#endif

               IF (inumber.GT.0) THEN
                  CALL MPI_RECV(iunique(istart), inumber, MPI_INTEGER8,
     &                 iprocrec,200, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER8,icheck,ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 200 ',
     &                    iproc
                     CALL quit
                  ENDIF

c                  print *,iproc,': Recved ',(iunique(ijk),ijk=istart,
c     &                 istart+inumber)

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

                CALL MPI_RECV(gradhs(1,istart),inumber,i2REAL4,
     &                 iprocrec, 26, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i2REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 8 ',
     &                    iproc
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': got gradhs from ',iprocrec,icheck
#endif

                  IF (imhd.EQ.idim) THEN
                     CALL MPI_RECV(Bxyz(1,istart),inumber,i3REAL8,
     &                  iprocrec, 27, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,i3REAL8,icheck,ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 6 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got Bxyz from ',iprocrec,icheck
#endif

                     CALL MPI_RECV(Bevolxyz(1,istart),inumber,i3REAL8,
     &                  iprocrec, 28, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,i3REAL8,icheck,ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 7 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got Bxyz from ',iprocrec,icheck
#endif
                  ENDIF

                  IF (ifsvi.EQ.6 .OR. imhd.EQ.idim) THEN
                     CALL MPI_RECV(alphaMM(1,istart), inumber,
     &                    ialphaREAL4, iprocrec, 29, MPI_COMM_WORLD, 
     &                    istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, ialphaREAL4, icheck, 
     &                    ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 8 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got alphaMM from ',iprocrec,icheck
#endif
                  ENDIF
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
         iorig(jpart) = jpart

         llist(nlst_end + j) = jpart
         iphase(jpart) = 0

ccc         CALL treef(jpart,npart,ntot,xyzmh,acc,0,fsx,fsy,fsz,epot)

         DO k = 1, 3
            dvxyzu(k,jpart) = 0.0
         END DO
         poten(jpart) = 0.0
         dq(jpart) = 0.0
         IF (imhd.EQ.idim) THEN
            DO k = 1, 3
               dBevolxyz(k,jpart) = 0.
            END DO
         ENDIF
      END DO
      nlst_tot = nlst_end + inumbertotal

      print *,nlst_tot,nlst_end,inumbertotal

#ifdef MPIDEBUG
      print *,iproc,': Received all particles with remote neighbours'
      print *,' '
      print *,' '
#endif

#else
      nlst_tot = nlst_end
#endif
      IF (nlst_tot.NE.ncompact) THEN
         WRITE (*,*) 'ERROR - ncompact.NE.nlst_tot ',ncompact,nlst_tot,
     &        nlst_end
         CALL quit
      ENDIF
c     
c--Compute implicit radiative transfer
c     
      IF (itiming) CALL getused(tass1)

      IF(encal.EQ.'r') THEN
         WRITE (*,*) 'Calling ass at realtime ',dt*itime/imaxstep+gt
         WRITE(*,*) 'radiation/gradh not yet implemented'
         STOP

c        WRITE (*,*) 'nlstall is ',nlstall

c        CALL ASS(nlst_in,nlst_end,llist,dt,itime,npart,ntot,
c     &           xyzmh,vxyzu,ekcle,dumrho,dedxyz,alphaMM)

cC$OMP PARALLEL DO SCHEDULE(runtime) default(none)
cC$OMP& shared(nlstall,vxyzu,dumrho,pr,vsound,llist,iphase,ekcle)
cC$OMP& private(i,ipart)
c            DO i = 1, nlstall
c               ipart = llist(i)
c               IF (iphase(ipart).EQ.0) 
c     &              CALL eospg(ipart,vxyzu,dumrho,pr,vsound,ekcle)
c            END DO
cC$OMP END PARALLEL DO
      END IF

      IF (itiming) THEN
         CALL getused(tass2)
         tass = tass + (tass2 - tass1)
      ENDIF
c
c--Compute forces on EACH particle
c
      IF (itiming) CALL getused(tforce1)

#ifdef MPIDEBUG
      print *,iproc,': calling forcei'
#endif


c      DO i = 1, nlst_tot
c         ipart = llist(i)
c         DO k = 1, 3
c            dvxyzu(k,ipart) = 0.0
c         END DO
c      END DO



c      CALL forcei(nlst_in,nlst_tot,llist,dt,itime,npart,ntot,
c     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,
c     &     dedxyz,Bxyz,dBevolxyz,Bevolxyz)
      CALL forcei(nlst_in,nlst_end,nlst_tot,llist,dt,itime,npart,ntot,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,
     &     dedxyz,Bxyz,dBevolxyz,Bevolxyz)

#ifdef MPIDEBUG
      print *,iproc,': called forcei '
#endif

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

#ifdef MPI
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
c--Receive back forces, du, potential energy, dBevolxyz to be added on 
c     to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' forces ',j
#endif
               CALL MPI_RECV(iready,1,
     &              MPI_INTEGER,MPI_ANY_SOURCE,98,MPI_COMM_WORLD,
     &              istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)

               CALL MPI_RECV(dvxyzu(1,istartrec),4*maxnneighsentback,
     &              MPI_REAL8,iprocrec,30,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
c               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.4*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.4*nnsentback dvxyzu'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got forces from ',iprocrec,' put into ',
     &              llistsentback(1,iprocrec+1)+1
#endif
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
c                  IF (ipos.EQ.1164) THEN
c                     print *,iproc,ipos,jpos,' extra ',dvxyzu(1,ipos),
c     &                    dvxyzu(1,jpos)
c                  ENDIF
                  DO k = 1, 4
                     dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
c                     IF (iunique(iorig(ipos)).EQ.2757) THEN
c                        print *,iproc,': FOR ',k,dvxyzu(k,ipos),
c     &                       dvxyzu(k,jpos),iprocrec,
c     &                       ipos,jpos
c                     ENDIF
                  END DO
               END DO

#ifdef MPIDEBUG
c      IF (iproc.EQ.0) 
c     &     print *,'called forcei PRE3 ',dvxyzu(1,68),
c     &              dvxyzu(2,68),dvxyzu(3,68),dvxyzu(1,2325),
c     &              dvxyzu(2,2325),dvxyzu(3,2325)
#endif

               CALL MPI_RECV(poten(istartrec),maxnneighsentback,
     &              MPI_REAL4,iprocrec,31,MPI_COMM_WORLD,istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsendback poten'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
c               print *,iproc,': got poten from ',iprocrec
#endif

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  poten(ipos) = poten(ipos) + poten(jpos)
               END DO

               IF (imhd.EQ.idim) THEN
               CALL MPI_RECV(dBevolxyz(1,istartrec),3*maxnneighsentback,
     &              MPI_REAL8,iprocrec,32,MPI_COMM_WORLD,istatus,
     &              ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
                  IF (ireturned.NE.3*nneighsentback(iprocrec+1)) THEN
               WRITE (*,*) 'ERROR - ireturned.NE.3*nnsentback dBevolxyz'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got dBevolxyz from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif
                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1
                     DO k = 1, 3
                        dBevolxyz(k,ipos) = dBevolxyz(k,ipos) + 
     &                       dBevolxyz(k,jpos)
c                        IF (iunique(iorig(ipos)).EQ.2757) THEN
c                           print *,iproc,': dB ',k,dBevolxyz(k,ipos),
c     &                          dBevolxyz(k,jpos),iprocrec,
c     &                          ipos,jpos
c                        ENDIF
                     END DO
                  END DO
               ENDIF

#ifdef MPIDEBUG
               print *,iproc,': received forces, poten ',j,' of ',
     &              inumofsends
#endif

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
                  iready = 1
                  CALL MPI_SEND(iready,1,MPI_INTEGER,
     &                 iprocsend,98,MPI_COMM_WORLD, ierr)
                  CALL MPI_SEND(dvxyzu(1,istart),4*inumber,MPI_REAL8,
     &                 iprocsend,30,MPI_COMM_WORLD, ierr)
                  CALL MPI_SEND(poten(istart),inumber,MPI_REAL4,
     &                 iprocsend,31,MPI_COMM_WORLD, ierr)

                  IF (imhd.EQ.idim) THEN
                     CALL MPI_SEND(dBevolxyz(1,istart),3*inumber,
     &                 MPI_REAL8,iprocsend,32,MPI_COMM_WORLD, ierr)
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': sent forces and poten to ',iprocsend
#endif
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                  GOTO 7765
               ENDIF
            END DO
 7765       CONTINUE

         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO

      GOTO 779

      IF (iproc.EQ.0) THEN
         IF (gt.EQ.0.) THEN
            DO n = 1, nlst
               ipart = llist(n)
               IF (xyzmh(1,ipart).LT.0.2 .AND. xyzmh(1,ipart).GT.0.0
     &              .AND. xyzmh(2,ipart).GT.1.0
     &              .AND. xyzmh(2,ipart).LT.1.2 
     &              .AND. xyzmh(3,ipart).LT.0.2
     &              .AND. xyzmh(3,ipart).GT.0.) THEN                  
                  print *,'Particle near x-bound ',ipart
c                  write (88,99004) ipart
99004             FORMAT('# ',I8)
                  iparticlex = ipart
               ENDIF
               IF (xyzmh(1,ipart).LT.1.2 .AND. xyzmh(1,ipart).GT.1.0
     &              .AND. xyzmh(2,ipart).GT.0.0
     &              .AND. xyzmh(2,ipart).LT.0.2 
     &              .AND. xyzmh(3,ipart).LT.0.2
     &              .AND. xyzmh(3,ipart).GT.0.) THEN
                  print *,'Particle near y-bound ',ipart
c                  write (89,99004) ipart
                  iparticley = ipart
               ENDIF
            END DO
         ENDIF
         write (88,99903) xyzmh(1,iparticlex),xyzmh(2,iparticlex),
     &        xyzmh(3,iparticlex),xyzmh(5,iparticlex),
     &        rho(iparticlex),divv(iparticlex),
     &        gradhs(1,iparticlex),gradhs(2,iparticlex),
     &        dvxyzu(1,iparticlex),dvxyzu(2,iparticlex),
     &        dvxyzu(3,iparticlex),
     &        dvxyzu(4,iparticlex),gt,
     &        nneigh(iparticlex)
         write (89,99903) xyzmh(1,iparticley),xyzmh(2,iparticley),
     &        xyzmh(3,iparticley),xyzmh(5,iparticley),
     &        rho(iparticley),divv(iparticley),
     &        gradhs(1,iparticley),gradhs(2,iparticley),
     &        dvxyzu(1,iparticley),dvxyzu(2,iparticley),
     &        dvxyzu(3,iparticley),
     &        dvxyzu(4,iparticley),gt,
     &        nneigh(iparticley)
         CALL FLUSH(88)
         CALL FLUSH(89)
      ENDIF


      GOTO 779

      IF (iproc.EQ.0) THEN
         DO n = 1, nlst
            ipart = llist(n)
            write (88,99902) xyzmh(1,ipart),xyzmh(2,ipart),
     &           xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
         END DO
      ELSE
         DO n = 1, nlst
            ipart = llist(n)
            write (89,99902) xyzmh(1,ipart),xyzmh(2,ipart),
     &           xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
         END DO
      ENDIF

 779  CONTINUE


      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR - inumberreturned.NE.inumofreturns ', iproc
         CALL quit
      ENDIF

      CALL MPI_TYPE_FREE(ialphaREAL4,ierr)
      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
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

c         DO n = 1, nlst
c            ipart = llist(n)
c            write (88,99901) xyzmh(1,ipart),xyzmh(2,ipart),
c     &           xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
cc     &           gradhs(1,ipart),gradhs(2,ipart),
c     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
c     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
c         END DO
99901 FORMAT (11(1PE12.5,1X),I8)
99902 FORMAT (13(1PE12.5,1X),I8)
99903 FORMAT (13(1PE17.10,1X),I8)

c      STOP

c      DO ipos = 1, npart
c         IF (iunique(iorig(ipos)).EQ.2757) THEN
c            DO k = 1,4
c               print *,iproc,': FOR-FINAL ',k,dvxyzu(k,ipos)
c            END DO
c            DO k = 1,3
c               print *,iproc,': dB-FINAL ',k,Bxyz(k,ipos),
c    &              Bevolxyz(k,ipos),dBevolxyz(k,ipos)
c            END DO
c         ENDIF             
c      END DO
                   
c      CALL quit
                   
      RETURN
      END
