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

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION xyzmh(5,mmax2),vxyzu(4,idim2),dvxyzu(4,idim3)
      REAL*4 dha(1+isizealphaMM,idim2),alphaMM(isizealphaMM,idim2)
      DIMENSION ireal(idim)
      DIMENSION ekcle(5,iradtrans2)
      DIMENSION Bevolxyz(imhdevol,imhd2),dBevolxyz(imhdevol,imhd3)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
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
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/planetesimal'

      DIMENSION dedxyz(3,iradtrans2)
#ifdef MPI
      DIMENSION listex(idim)
#endif

c
c--Allow for tracing flow
c
c      IF (gt.GT.50.) 
#ifdef MPIDEBUG
      print *,iproc,': ENTERED derivi ',nlst
#endif

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

      IF (iexf.EQ.7 .OR. ibound.EQ.102) THEN
         realtime = dt*itime/imaxstep+gt
         IF (rplanet.GT.1.0E-10) THEN
            pradfac(1) = (rplanet + (0.01*exp(-4.*realtime/pi)))/rplanet
c--Double shrink time.
c            pradfac(1) = (rplanet + (0.01*exp(-2.*realtime/pi)))/rplanet
            IF (pradfac(1).LT.1.001) pradfac(1) = 1.0
         ELSEIF (initialptm.EQ.5) THEN
            DO i = 1, nptmass
               pradfac(i) = (xyzmh(5,listpm(i)) + (0.01*exp(-4.*
     &              realtime/pi)))/xyzmh(5,listpm(i))
            IF (pradfac(i).LT.1.001) pradfac(i) = 1.0
            ENDDO
         ELSE
            pradfac(:) = 0.0
         ENDIF
      ENDIF


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
         print *,iproc,' DERIVI: calling densityiterate_gradh'
#endif

         CALL densityiterate_gradh(dt,npart,ntot,xyzmh,vxyzu,dvxyzu,
     &        nlst_in,nlst_end,llist,itime,ekcle,Bevolxyz,Bxyz,
     &        dBevolxyz)

#ifdef MPIDEBUG
         print *,iproc,' DERIVI: called densityiterate_gradh'
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
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, imhdevolREAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(isizealphaMM,MPI_REAL4,ialphaREAL4,ierr)

      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(imhdevolREAL8,ierr)
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
     &                       ' starting at ',llistsentback(1,j+1),'+1 ',
     &                       nneighsentback(j+1),' list '
c     &               ,(llistsentback(kkk,j+1),kkk=1,nneighsentback(j+1))
#endif
c
c--Else does not need to send any particle back to this particular process.
c
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

                        DO ijk = 1, nneighsentback(j+1)
                        listex(ijk) = iorig(llistsentback(ijk,j+1)+1)-1
                        END DO

                        CALL MPI_SEND(divv,1,indexMPI1,j,23,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent divv'
#endif
                        IF (nneighsentback(j+1).GT.0) THEN

                        IF (encal.EQ.'r') THEN
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i5REAL8,indexMPI5,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i5REAL8,indexMPI5,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI5,ierr)
                        ENDIF

                        IF (imhd.EQ.idim) THEN

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)

                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),imhdevolREAL8,
     &                       indexMPIimhdevol,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIimhdevol,ierr)

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

                           CALL MPI_SEND(Bevolxyz,1,indexMPIimhdevol,
     &                          j,28,MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent Bevolxyz to ',j
#endif

                           CALL MPI_TYPE_FREE(indexMPI3,ierr)
                           CALL MPI_TYPE_FREE(indexMPIimhdevol,ierr)
                        ENDIF

                        IF (ifsvi.GE.6 .OR. imhd.EQ.idim) THEN
                           CALL MPI_SEND(alphaMM,1,indexMPIalpha,j,29,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent alphaMM to ',j
#endif                           
                        ENDIF

                        IF (encal.EQ.'r') THEN
                           CALL MPI_SEND(ekcle,1,indexMPI5,j,230,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent ekcle to ',j
#endif
                           CALL MPI_TYPE_FREE(indexMPI5,ierr)
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
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (istart+inumber.GE.idim2) THEN
                  WRITE (*,*) iproc,' ERROR - istart+inumber.GE.idim ',
     &                 istart,inumber,idim,ntot,inumbertotal
                  CALL quit
               ENDIF
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

#ifdef MPIDEBUG
c                  print *,iproc,': Recved ',(iunique(ijk),ijk=istart,
c     &                 istart+inumber)
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

                     CALL MPI_RECV(Bevolxyz(1,istart),inumber,
     &                    imhdevolREAL8, iprocrec, 28, 
     &                    MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,imhdevolREAL8,
     &                    icheck,ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 7 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got Bxyz from ',iprocrec,icheck
#endif
                  ENDIF

                  IF (ifsvi.GE.6 .OR. imhd.EQ.idim) THEN
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

                  IF (encal.EQ.'r') THEN
                     CALL MPI_RECV(ekcle(1,istart),inumber,i5REAL8,
     &                    iprocrec, 230, MPI_COMM_WORLD,istatus,ierr)
                     CALL MPI_GET_COUNT(istatus,i5REAL8,icheck,ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 9 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got ekcle from ',iprocrec,
     &                    icheck,istart
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
c--Need to set up lists for remote particles and ensure
c     that forces and poten are zeroed before call of forcei.
c
      DO j = 1, inumbertotal
         jpart = ntot + j
         iorig(jpart) = jpart

         llist(nlst_end + j) = jpart

         DO k = 1, 3
            dvxyzu(k,jpart) = 0.0
         END DO
         poten(jpart) = 0.0
         dq(jpart) = 0.0
         IF (imhd.EQ.idim) THEN
            DO k = 1, imhdevol
               dBevolxyz(k,jpart) = 0.
            END DO
            DO k = 1, 5
               divcurlB(k,jpart) = 0.
            END DO
            DO k = 1, 9
               gradB(k,jpart) = 0.
            END DO
         ENDIF
      END DO
      nlst_tot = nlst_end + inumbertotal

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
      IF (nlst_end.NE.ncompactlocal) THEN
         WRITE (*,*) 'ERROR - ncompactlocal.NE.nlst_end ',ncompactlocal,
     &        nlst_end,nlst_tot,ncompact
         CALL quit
      ENDIF
c
c--Compute implicit magnetic resistivity
c
      IF (imhd.EQ.idim .AND. (iresist.LT.0 .OR. 
     &     isubcycle_divB.GT.0 .AND. ndivBsubcycles.GT.0)) THEN
       IF (iresist.LT.0) THEN
         CALL resistivity_implicit(nlst_in,nlst_end,llist,dt,itime,
     &        npart,ntot,Bxyz,xyzmh,vxyzu,dumrho,vsound,alphaMM,nit,
     &        xerror)

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
       ENDIF
c
c--Compute divBcleaning via subcycling (NOTE: If used with implicit resistivity
c     there MAY NEED to be an MPI transfer between here and the above implicit
c     resistivity call)
c
         IF (imhdevol.EQ.4 .AND. 
     &      isubcycle_divB.GT.0 .AND. ndivBsubcycles.GT.0) THEN
            IF (itiming) CALL getused(tdivBsub1)

            CALL divBclean_subcycle(dt,npart,ntot,xyzmh,dumrho,Bevolxyz)

            IF (varmhd.EQ.'Bvol') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,Bevolxyz,Bxyz,llist)
C$OMP& private(i,ipart)
               DO i = nlst_in,nlst_end
                  ipart = llist(i)
                  Bxyz(1,ipart) = Bevolxyz(1,ipart)
                  Bxyz(2,ipart) = Bevolxyz(2,ipart)
                  Bxyz(3,ipart) = Bevolxyz(3,ipart)
               END DO
C$OMP END PARALLEL DO
            ELSEIF (varmhd.EQ.'Brho') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,Bevolxyz,Bxyz,llist,dumrho)
C$OMP& private(i,ipart)
               DO i = nlst_in,nlst_end
                  ipart = llist(i)
                  Bxyz(1,ipart) = Bevolxyz(1,ipart)*dumrho(ipart)
                  Bxyz(2,ipart) = Bevolxyz(2,ipart)*dumrho(ipart)
                  Bxyz(3,ipart) = Bevolxyz(3,ipart)*dumrho(ipart)
               END DO
C$OMP END PARALLEL DO
            ENDIF

            IF (itiming) THEN
               CALL getused(tdivBsub2)
               tdivBsub = tdivBsub + (tdivBsub2 - tdivBsub1)
            ENDIF
         ENDIF

#ifdef MPI
c
c--Transfer new values of Bxyz, Bevolxyz to other MPI processes
c
      inumbertotallocal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
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

                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)

                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),imhdevolREAL8,
     &                  indexMPIimhdevol,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIimhdevol,ierr)

                        CALL MPI_SEND(Bxyz,1,indexMPI3,j,27,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent Bxyz to ',j
#endif

                        CALL MPI_SEND(Bevolxyz,1,indexMPIimhdevol,j,28,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent Bevolxyz to ',j
#endif

                        CALL MPI_TYPE_FREE(indexMPIimhdevol,ierr)
                        CALL MPI_TYPE_FREE(indexMPI3,ierr)
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
               print *,iproc,': receiving Bxyz, Bevolxyz data ',
     &              inumberreturned,inumofreturns
#endif
               istart = ntot + inumbertotallocal + 1

               CALL MPI_RECV(Bxyz(1,istart),idim,i3REAL8,
     &              i, 27, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,i3REAL8,inumber,ierr)
               icheck = inumber
               IF (istart + icheck.GE.idim2) THEN
                  WRITE (*,*) iproc,': ERROR - istart+icheck.GE.idim2 ',
     &                 istart,icheck,idim2,ntot,inumbertotallocal,
     &                 inumberreturned,inumofreturns,
     &                 nneightogetback(i+1)
                  CALL quit
               ENDIF
               iprocrec = istatus(MPI_SOURCE)
#ifdef MPIDEBUG
               print *,iproc,': got Bxyz from ',iprocrec,inumber
#endif

               CALL MPI_RECV(Bevolxyz(1,istart),idim,imhdevolREAL8,
     &              iprocrec, 28, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,imhdevolREAL8,inumber,ierr)
               IF (icheck.NE.inumber) THEN
                  WRITE (*,*) 'ERROR - icheck.NE.inumber 17 ',
     &                 iproc
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got Bxyz from ',iprocrec,inumber
#endif

               inumbertotallocal = inumbertotallocal + inumber
               IF (inumberindiv(inumberreturned).NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv B ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit
               ENDIF
               IF (inumbercumm(inumberreturned).NE.istart) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm B ',istart,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberprocB ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit
      ENDIF

#endif

      ENDIF   ! Endif for imhd.EQ.idim and iresist<0 or subcycling
c
c--Compute implicit radiative transfer
c     
      IF (itiming) CALL getused(tass1)

      IF(encal.EQ.'r' .OR. (encal.EQ.'i' .AND. ibound.EQ.102) 
     &     .OR. (encal.EQ.'i' .AND. ibound.EQ.100)) THEN

ccc         thermalenergy = 0.
ccc         radenergy = 0.

c         DO i = 1, npart
c            IF (iphase(i).EQ.0) THEN
c               thermalenergy = thermalenergy + xyzmh(4,i)*vxyzu(4,i)
c               radenergy = radenergy + xyzmh(4,i)*ekcle(1,i)
c            ENDIF
c         END DO

ccc         DO n = 1, ncompactlocal
ccc            i = ivar(3,n)

c            IF (iunique(iorig(i)).EQ.79225) THEN
c               print *,'COM ',i,iunique(iorig(i)),iproc,ivar(1,n)

c            print *,' TRY ',iproc,icall,i,ekcle(1,i)

ccc               thermalenergy = thermalenergy + xyzmh(4,i)*vxyzu(4,i)
ccc               radenergy = radenergy + xyzmh(4,i)*ekcle(1,i)

c               DO k = 1, ivar(1,n)
c                  icompact = ivar(2,n) + k
c                  j = ijvar(icompact)
c
c                  print *,'COM-N ',j
c
c                  thermalenergy = thermalenergy + xyzmh(4,j)*vxyzu(4,j)
c                  radenergy = radenergy + xyzmh(4,j)*ekcle(1,j)
c               END DO
c            ELSEIF (iproc.EQ.0) THEN
c               DO k = 1, ivar(1,n)
c                  icompact = ivar(2,n) + k
c                  j = ijvar(icompact)

c                  IF (iunique(iorig(j)).EQ.79225) THEN

c                     DO nnn = 1, ncompactlocal
c                        IF (iunique(iorig(ivar(3,nnn))).EQ.79225) THEN
c                           DO kkk = 1, ivar(1,nnn)
c                              icompact2 = ivar(2,nnn) + kkk
c                              IF (i.EQ.ijvar(icompact2)) GOTO 8911
c                           END DO
c                        ENDIF
c                     END DO

c                     print *,'COM-J ',j,i

c                  thermalenergy = thermalenergy + xyzmh(4,i)*vxyzu(4,i)
c                     radenergy = radenergy + xyzmh(4,i)*ekcle(1,i)

c 8911                CONTINUE
c                  ENDIF
c               END DO

c            ENDIF
ccc         END DO
            
ccc         print 99555,iproc,icall,thermalenergy,radenergy,
ccc     &        thermalenergy+radenergy
ccc99555    FORMAT(I2,': ENER-B ',I1,1PE14.7,1PE14.7,1PE14.7)

c         CALL quit

         IF (encal.EQ.'r') THEN
c            WRITE (*,*) 'Calling ass at realtime ',dt*itime/imaxstep+gt
            CALL ASS(nlst_in,nlst_end,llist,dt,itime,npart,ntot,
     &           xyzmh,vxyzu,ekcle,dumrho,dedxyz,alphaMM)
         ELSEIF (ibound.EQ.102 .AND. use_tprof) THEN
c-- Added for Isothermal whole disc case with t-profile
            rho0 = 0.2
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,llist,xyzmh,gamma,gmw,iphase)
C$OMP& shared(hoverr,tprof,Rg,uergg,vxyzu,rho0,rho)
C$OMP& private(i,j,radius,boundtempl)
            DO j = nlst_in, nlst_end
               i = llist(j)
               IF (iphase(i).EQ.0) THEN
                  radius = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
                  boundtempl =  gmw*hoverr**2*radius**
     &                 (tprof+1)/((Rg/uergg)*gamma*radius)
                  vxyzu(4,i) = getu(rho(i), boundtempl)
                  IF (rho(i).GT.rho0) THEN
                     vxyzu(4,i) = vxyzu(4,i)*(rho(i)/rho0)**0.8
                  ENDIF
               ENDIF
            END DO
C$OMP END PARALLEL DO
         ELSEIF (ibound.EQ.100 .OR. ibound.EQ.102) THEN
c--Added for Isothermal disc section case and whole disc with u profile
            rho0 = 0.2
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,llist,xyzmh,gamma,iphase)
C$OMP& shared(hoverr,vxyzu,rho0,rho)
C$OMP& private(i,j)
            DO j = nlst_in, nlst_end
               i = llist(j)
               IF (iphase(i).EQ.0) THEN
                  vxyzu(4,i)=hoverr**2/(SQRT(xyzmh(1,i)**2+
     &                 xyzmh(2,i)**2)*gamma*(gamma-1.0))
                  IF (rho(i).GT.rho0) THEN
                     vxyzu(4,i) = vxyzu(4,i)*(rho(i)/rho0)**0.8
                  ENDIF
               ENDIF
            END DO
C$OMP END PARALLEL DO
         ENDIF
         
ccc         thermalenergy = 0.
ccc         radenergy = 0.
c         DO i = 1, npart
c            IF (iphase(i).EQ.0) THEN
c               thermalenergy = thermalenergy + xyzmh(4,i)*vxyzu(4,i)
c               radenergy = radenergy + xyzmh(4,i)*ekcle(1,i)
c            ENDIF
c         END DO

ccc         DO n = 1, ncompactlocal
ccc            i = ivar(3,n)

c            IF (iunique(iorig(i)).EQ.79225) THEN
c               print *,'COM ',i,iunique(iorig(i)),iproc,ivar(1,n)

c            print *,' TRY ',iproc,icall,i,ekcle(1,i)

ccc               thermalenergy = thermalenergy + xyzmh(4,i)*vxyzu(4,i)
ccc               radenergy = radenergy + xyzmh(4,i)*ekcle(1,i)

c               DO k = 1, ivar(1,n)
c                  icompact = ivar(2,n) + k
c                  j = ijvar(icompact)
c
c                  print *,'COM-N ',j
c
c                  thermalenergy = thermalenergy + xyzmh(4,j)*vxyzu(4,j)
c                  radenergy = radenergy + xyzmh(4,j)*ekcle(1,j)
c               END DO
c            ELSEIF (iproc.EQ.0) THEN
c               DO k = 1, ivar(1,n)
c                  icompact = ivar(2,n) + k
c                  j = ijvar(icompact)

c                  IF (iunique(iorig(j)).EQ.79225) THEN

c                     DO nnn = 1, ncompactlocal
c                        IF (iunique(iorig(ivar(3,nnn))).EQ.79225) THEN
c                           DO kkk = 1, ivar(1,nnn)
c                              icompact2 = ivar(2,nnn) + kkk
c                              IF (i.EQ.ijvar(icompact2)) GOTO 8910
c                           END DO
c                        ENDIF
c                     END DO

c                     print *,'COM-J ',j,i

c                  thermalenergy = thermalenergy + xyzmh(4,i)*vxyzu(4,i)
c                     radenergy = radenergy + xyzmh(4,i)*ekcle(1,i)

c 8910                CONTINUE
c                  ENDIF
c               END DO

c            ENDIF
ccc         END DO


ccc         print 99556,iproc,icall,thermalenergy,radenergy,thermalenergy+
ccc     &        radenergy
ccc99556    FORMAT(I2,': ENER-A ',I1,1PE14.7,1PE14.7,1PE14.7)

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,vxyzu,dumrho,pr,vsound,llist,iphase)
C$OMP& shared(ekcle)
C$OMP& shared(gasdist,gasdrag)
C$OMP& private(i,ipart)

            DO i = nlst_in, nlst_end
               ipart = llist(i)
               IF (iphase(ipart).EQ.0) 
     &              CALL eospg(ipart,vxyzu,dumrho,pr,vsound,ekcle)
               IF (iphase(ipart).EQ.11 .AND. gasdrag)
     &              gasdist(ipart) = 1.0E6
            END DO

C$OMP END PARALLEL DO





#ifdef MPI
c
c--Transfer new values of vxyzu, ekcle, pr, vsound
c     to other MPI processes
c
      inumbertotallocal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
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

                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

                        CALL MPI_SEND(vxyzu,1,indexMPI4,j,210,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexMPI4,ierr)


                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i5REAL8,indexekcle,ierr)
                        CALL MPI_TYPE_COMMIT(indexekcle,ierr)

                        CALL MPI_SEND(ekcle,1,indexekcle,j,211,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexekcle,ierr)
#ifdef MPIDEBUGRT
                        print *,iproc,' sent ekcle'
#endif

                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),MPI_REAL4,indexMPI_R4,ierr)
                CALL MPI_TYPE_COMMIT(indexMPI_R4,ierr)

                        CALL MPI_SEND(pr,1,indexMPI_R4,j,212,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(vsound,1,indexMPI_R4,j,213,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexMPI_R4,ierr)

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
#ifdef MPIDEBUGRT
               print *,iproc,': receiving vxyzu,ekcle,pr,vsound data ',
     &              inumberreturned,inumofreturns
#endif
               istart = ntot + inumbertotallocal + 1
               CALL MPI_RECV(vxyzu(1,istart), idim, MPI_REAL8,
     &              i, 210, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               icheck = inumber/4
               IF (istart + icheck.GE.idim2) THEN
                  WRITE (*,*) iproc,': ERROR - istart+icheck.GE.idim2 ',
     &                 istart,icheck,idim2,ntot,inumbertotallocal,
     &                 inumberreturned,inumofreturns,
     &                 nneightogetback(i+1)
                  CALL quit
               ENDIF
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUGRT
               print *,iproc,': got vxyzu from ',iprocrec,inumber
#endif
               CALL MPI_RECV(ekcle(1,istart), idim, MPI_REAL8,
     &              i, 211, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (icheck.NE.inumber/5) THEN
                  WRITE (*,*) iproc,': ERROR - icheck.NE.inumber/5'
                  CALL quit
               ENDIF

#ifdef MPIDEBUGRT
               print *,iproc,': got ekcle from ',iprocrec,inumber
#endif
               CALL MPI_RECV(pr(istart), idim, MPI_REAL4,
     &              i, 212, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (icheck.NE.inumber) THEN
                  WRITE (*,*) iproc,': ERROR - icheck.NE.inumber A'
                  CALL quit
               ENDIF

               CALL MPI_RECV(vsound(istart), idim, MPI_REAL4,
     &              i, 213, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (icheck.NE.inumber) THEN
                  WRITE (*,*) iproc,': ERROR - icheck.NE.inumber B'
                  CALL quit
               ENDIF

               inumbertotallocal = inumbertotallocal + inumber
               IF (inumberindiv(inumberreturned).NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv RT ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit
               ENDIF
               IF (inumbercumm(inumberreturned).NE.istart) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm RT ',istart,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberproc RT ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit
      ENDIF

#endif


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
      print *,iproc,': calling forcei ',itime
#endif


      CALL forcei(nlst_in,nlst_end,nlst_tot,llist,dt,itime,npart,ntot,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,
     &     dedxyz,Bxyz,dBevolxyz,Bevolxyz)


#ifdef MPIDEBUG
      print *,iproc,': called forcei ',itime
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
      IF (istartrec + maxnneighsentback.GT.idim3) THEN
         WRITE (*,*) 'ERROR - istartrec + maxnneighsentback.GT.idim3'
         CALL quit
      ENDIF
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back forces, du, potential energy, dBevolxyz, divcurlB 
c     to be added on to local values
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
               IF (istartrec+ireturned.GT.idim3) THEN
                  WRITE (*,*) 'ERROR - istartrec+ireturned.GT.idim3 ',
     &                 istartrec,ireturned,idim3
                  CALL quit
               ENDIF
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
                  DO k = 1, 4
                     dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
                  END DO
               END DO

               IF (istartrec+nneighsentback(iprocrec+1).GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istartrec+nneighsentback',
     &                 '(iprocrec+1).GT.idim2 ',istartrec,
     &                 nneighsentback(iprocrec+1),idim2
                  CALL quit
               ENDIF

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
               CALL MPI_RECV(dBevolxyz(1,istartrec),
     &                 imhdevol*maxnneighsentback,
     &                 MPI_REAL8,iprocrec,32,MPI_COMM_WORLD,istatus,
     &                 ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
                  IF (ireturned.NE.imhdevol*nneighsentback(iprocrec+1))
     &                 THEN
               WRITE (*,*) 'ERROR - ireturned.NE.imhdevol',
     &                    '*nnsentback dBevolxyz ',imhdevol,
     &                    nneighsentback(iprocrec+1),ireturned
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got dBevolxyz from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif
                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1
                     DO k = 1, imhdevol
                        dBevolxyz(k,ipos) = dBevolxyz(k,ipos) + 
     &                       dBevolxyz(k,jpos)
                     END DO
                  END DO

               CALL MPI_RECV(divcurlB(1,istartrec),5*maxnneighsentback,
     &              MPI_REAL4,iprocrec,33,MPI_COMM_WORLD,istatus,
     &              ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL4,ireturned,ierr)
                  IF (ireturned.NE.5*nneighsentback(iprocrec+1)) THEN
               WRITE (*,*) 'ERROR - ireturned.NE.5*nnsentback divcurlB'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got divcurlB from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif
                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1
                     DO k = 1, 5
                        divcurlB(k,ipos) = divcurlB(k,ipos) + 
     &                       divcurlB(k,jpos)
                     END DO
                  END DO

               CALL MPI_RECV(gradB(1,istartrec),9*maxnneighsentback,
     &              MPI_REAL4,iprocrec,34,MPI_COMM_WORLD,istatus,
     &              ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL4,ireturned,ierr)
                  IF (ireturned.NE.9*nneighsentback(iprocrec+1)) THEN
               WRITE (*,*) 'ERROR - ireturned.NE.9*nnsentback gradB'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got gradB from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif
                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1
                     DO k = 1, 9
                        gradB(k,ipos) = gradB(k,ipos) + 
     &                       gradB(k,jpos)
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
                     CALL MPI_SEND(dBevolxyz(1,istart),imhdevol*inumber,
     &                 MPI_REAL8,iprocsend,32,MPI_COMM_WORLD, ierr)
                     CALL MPI_SEND(divcurlB(1,istart),5*inumber,
     &                 MPI_REAL4,iprocsend,33,MPI_COMM_WORLD, ierr)
                     CALL MPI_SEND(gradB(1,istart),9*inumber,
     &                 MPI_REAL4,iprocsend,34,MPI_COMM_WORLD, ierr)
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

      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR - inumberreturned.NE.inumofreturns ', iproc
         CALL quit
      ENDIF

      CALL MPI_TYPE_FREE(ialphaREAL4,ierr)
      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(imhdevolREAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)

#ifdef MPIDEBUG
      print *,iproc,': FINISHED MPI force transfer ',itime
      print *,' '
#endif

#endif

c
c--Calculate source terms for variable viscosity and resistivity
c     This used to be done in forcei, but with MPI could not be because
c     divcurlB (e.g. with B/rho) cannot be calculated until the end of
c     forcei.  Therefore, the source term calculation for both the 
c     viscosity and the resistivity has been moved to this new routine
c     (which is essentially a very small fraction of the original forcei
c     routine).
c
      CALL dissipative_source(nlst_in,nlst_end,nlst_tot,llist,itime,
     &     npart,ntot,xyzmh,dha,dumrho,vsound,alphaMM,Bxyz,Bevolxyz,
     &     dBevolxyz)

#ifdef MPIDEBUG
      print *,iproc,': Exited derivi ',itime
#endif
                   
      RETURN
      END
