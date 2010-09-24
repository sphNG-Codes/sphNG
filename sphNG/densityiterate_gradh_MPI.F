      SUBROUTINE densityiterate_gradh (dt,npart,ntot,xyzmh,vxyzu,dvxyzu,
     &            nlst_in,nlst_end,listp,itime,ekcle,Bevol,Bxyz,dBxyz)
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

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2), dvxyzu(4,idim3)
      DIMENSION listp(idim2)
      DIMENSION ekcle(5,iradtrans2)
      DIMENSION Bevol(3,imhd3),Bxyz(3,imhd2),dBxyz(3,imhd2)

c--Neides=INT(4./3.*pi*8.*hfact**3))
      PARAMETER (hfact = 1.2)

      PARAMETER (htol = 1.e-3)
      PARAMETER (hstretch = 1.01)
      PARAMETER (maxiterations = 500)
      PARAMETER (maxitsnr = 30)
c--this weight is equivalent to m/(rho*h^3) in the grad h version
      PARAMETER (weight = 1./hfact**3)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/ener1'
c--treecom_P is included to use listparents to make list of inactive particles
c     that are neighbours of active particles for updating their densities etc
      INCLUDE 'COMMONS/treecom_P'
c      INCLUDE 'COMMONS/vsmooth'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/gravi'

      INTEGER*1 isucceed(idim)
      REAL*4 hminbisec(idim),hmaxbisec(idim)

      DIMENSION ldolist(idim),hi_old(idim)
      DIMENSION gradmhd(15,imhd2)
#ifdef MPI
      DIMENSION nneighrec(ineighproc),nneighsentbackold(nummaxproc)
      REAL*4 rhorec(ineighproc), dumrhorec(ineighproc), rhonextmax
      LOGICAL*1 isend(idim)
      INTEGER*8 iuniquesend(iptdim)
#endif

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine densityiterate')
c
c--Initialise
c
      iproblem = itime
#ifdef MPI
c      print *,iproc,': Entered densityiterate ',itime
#else
c      print *,': Entered densityiterate ',itime
#endif

      uradconst = radconst/uergcc
      third = 1./3.
      rhonext = 0.
      icreate = 0
      icreatetot = 0
      radcrit2 = radcrit*radcrit
      numparticlesdone = numparticlesdone + nlst_end
      nwarnup = 0
      nwarndown = 0
      nwarnroundoff = 0
      stressmax = 0.
      nbisection = 0

      ncompact = 0
      icompact = 0
      nlistupdated = 0
c
c--For constant pressure boundaries, use a minimum density
c  equal to the external density
c
c      IF (ibound.EQ.7) THEN
c         rhomin = 0.25*rhozero
c         print*,'rhomin = ',rhomin
c      ELSE
      rhomin = 0.      
c      ENDIF

c
c--Need to leave in sink particles.  Even though they do not need to be 
c     iterated, they need to be included for the gravitational forces from
c     gas particles on other MPI processes to be included (for iptintree=1,2).
c     This list is reused for the subsequent iterations when it is expected 
c     that it will contain fewer and fewer particles as gas particles have
c     their correct smoothing length and density set.
c
      nlstdo = 0
      nptmasssend = 0
      DO n = nlst_in, nlst_end
         ipart = listp(n)
         IF (iphase(ipart).GE.0) THEN
            nlstdo = nlstdo + 1
            ldolist(nlstdo) = listp(n)
#ifdef MPI
            IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
               nptmasssend = nptmasssend + 1
               iuniquesend(nptmasssend) = iunique(iorig(ipart))
            ENDIF
#endif
         ENDIF
      END DO
#ifdef MPI
      IF (individualtimesteps.NE.2) THEN
         IF (nptmasssend.NE.nptmass) THEN
            WRITE (*,*) 'ERROR - nptmasssend.NE.nptmass ',nptmasssend,
     &           nptmass
            CALL quit
         ENDIF
      ENDIF
#endif

      nlstsend = nlstdo

C$OMP PARALLEL DO SCHEDULE(runtime) default(none) 
C$OMP& shared(nlstdo,ldolist,isucceed,xyzmh,hi_old,it1,isteps,imaxstep)
C$OMP& shared(divv,rho,dumrho,nneigh,dvxyzu,rhomin,dt,imax,poten,dq)
C$OMP& shared(iphase,dBxyz,it0,itime,curlv,gradhs)
#ifdef MPI
C$OMP& shared(lsendlist,isend)
#endif
C$OMP& private(n,ipart,pmassi,hi,dhdrhoi,deltat,k,hchange)
      DO n = 1, nlstdo
         ipart = ldolist(n)
         isucceed(ipart) = 0

#ifdef MPI
         lsendlist(n) = ipart - 1
         isend(n) = .TRUE.
#endif
         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN
            pmassi = xyzmh(4,ipart)
            hi = xyzmh(5,ipart)
            hi_old(ipart) = hi
c
c--Predict h
c
            dhdrhoi = - hi/(3.*(pmassi*(hfact/hi)**3 + rhomin))
            IF (it1(ipart).EQ.itime) THEN
               deltat = (dt*isteps(ipart)/2)/imaxstep
            ELSEIF (it0(ipart).EQ.itime) THEN
               deltat = (dt*isteps(ipart))/imaxstep
            ELSE
               WRITE (*,*) 'ERROR - it1 AND it2 .NE. itime'
               CALL quit
            ENDIF

            hchange = - dhdrhoi*divv(ipart)*deltat

c            IF (ABS(hchange/xyzmh(5,ipart)).GT.0.5) THEN
c               print *,' LARGE ',ipart,xyzmh(5,ipart),dhdrhoi,
c     &              divv(ipart),deltat,it1(ipart),isteps(ipart),
c     &              imaxstep,dt,itime,xyzmh(1,ipart),xyzmh(2,ipart),
c     &              xyzmh(3,ipart)
c            ENDIF

            IF (ABS(hchange/xyzmh(5,ipart)).LT.0.5) THEN
               xyzmh(5,ipart) = xyzmh(5,ipart) + hchange
            ENDIF

            rho(ipart) = 0.
            dumrho(ipart) = 0.
            divv(ipart) = 0.
            curlv(ipart) = 0.
            gradhs(2,ipart) = 0.
         ENDIF
         nneigh(ipart) = 0
         DO k = 1, 4
            dvxyzu(k,ipart) = 0.
         END DO
         poten(ipart) = 0.
         dq(ipart) = 0.
         IF (imhd.EQ.idim) THEN
            DO k = 1, 3
               dBxyz(k,ipart) = 0.
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
c
c--Iterate density calculation for particle ipart
c
#ifdef MPI

      IF (itiming) CALL getused(td11)

      inumofreturns = 0
      inumofsends = 0
      maxnneighsentback = 0
      nneighsentanyatall = .FALSE.
C$OMP PARALLEL DO SCHEDULE(static) default(none)
C$OMP& shared(numproc,nneighsentany)
C$OMP& private(i)
      DO i = 1, numproc
         nneighsentany(i) = .FALSE.
      END DO
C$OMP END PARALLEL DO

      CALL MPI_TYPE_CONTIGUOUS(15, MPI_REAL8, i15REAL8, ierr)
c      CALL MPI_TYPE_CONTIGUOUS(12, MPI_REAL8, i12REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)

      CALL MPI_TYPE_COMMIT(i15REAL8,ierr)
c      CALL MPI_TYPE_COMMIT(i12REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)

      CALL MPI_ALLREDUCE(nlstdo,nlstdo_tot,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
#else
      nlstdo_tot = nlstdo
#endif

      DO iteration = 1, maxiterations

#ifdef MPI
#ifdef MPIDEBUG
c      IF (itime.EQ.2097152)
              print *,iproc,'DEN_ITERATE: iteration ',iteration
#endif
c
c--MPI: Need to get rho, gradh, numneigh contributions from other processes
c
#ifdef MPIDEBUG
         print *,iproc,': Start MPI densityiterate ',numproc



c         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


#endif
         iprocdebug = iproc

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlstsend, 1, lsendlist, 
c     &       i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_INDEXED(nlstsend, lblocklengths, lsendlist, 
     &        i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlstsend, 1, lsendlist, 
c     &       MPI_INTEGER1, indexMPI_INT1, ierr)
         CALL MPI_TYPE_INDEXED(nlstsend, lblocklengths, lsendlist, 
     &        MPI_INTEGER1, indexMPI_INT1, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--There are two implementations, one for balanced jobs, one for unbalanced
c
         nptmasslocal = nptmass
         nptmasslocal2 = nptmass
         DO ii = 1, numproc
            IF (ideniteratebal) THEN
               IF (ii.EQ.numproc) GOTO 444
c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
c              IF (nptmasstot.GE.1) print *,iproc,' CIRCULAR ',ii
               iahead = MOD(iproc+ii,numproc)
               ibehind = MOD(numproc+iproc-ii,numproc)
c
c--Send to node ahead, receive from node behind
c
#ifdef MPIDEBUG
               print *,iproc,': Sending xyzmh to ',iahead,' rec from ',
     &              ibehind,xyzmh(1,lsendlist(1)+1),numproc
#endif

               CALL MPI_SENDRECV(xyzmh,1,indexMPI5,iahead,0,
     &              xyzmh(1,2*ntot+3), idim, i5REAL8, ibehind,
     &              0, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
               IF (2*ntot+3+inumber.GE.mmax2) THEN
                  WRITE (*,*) iproc,': ERROR - 2*ntot+3+inumber',
     &                 '.GE.mmax2'
                  CALL quit
               ENDIF

#ifdef MPIDEBUG
               print *,iproc,': Received ',inumber,radkernel,
     &              xyzmh(1,2*ntot+3),numproc
#endif

               CALL MPI_SENDRECV(iphase,1,indexMPI_INT1,iahead,1,
     &              iphase(ntot+1), idim, MPI_INTEGER1, ibehind,
     &              1, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_INTEGER1,igotnumber,ierr)
               IF (ntot+1+inumber.GE.idim2) THEN
                  WRITE (*,*) iproc,': ERROR - ntot+1+inumber.GE.idim'
                  CALL quit
               ENDIF
               IF (igotnumber.NE.inumber) THEN
                  WRITE (iprint, *) 'ERROR - igotnumber.NE.inumber:',
     &                 'iphase'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               IF (itime.EQ.2097152)
     &              print *,iproc,': Received iphase ',inumber
#endif

               IF (iteration.EQ.1) THEN
                  CALL MPI_SENDRECV(iuniquesend, nptmasssend, 
     &                 MPI_INTEGER8, iahead, 9,
     &                 iuniquestore(nptmasslocal2+1), iptdim, 
     &                 MPI_INTEGER8, ibehind, 9, 
     &                 MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_INTEGER8,igotnumber,
     &                 ierr)
                  nptmasslocal2 = nptmasslocal2 + igotnumber
                  IF (nptmasslocal2.GT.iptdim) THEN
                     WRITE (*,*) 'ERROR - nptmasslocal2.GT.iptdim ',
     &                    nptmasslocal2, igotnumber
                     CALL quit
                  ENDIF
               ENDIF

            ELSE
c
c--IMPLEMENTATION USING SEND
c
               IF (ii.EQ.iproc+1) THEN
                  DO jj = 0, numproc - 1
                     IF (iproc.NE.jj) THEN
#ifdef MPIDEBUG
                        print *,iproc,': Sending ',nlstsend,
     &                       ' xyzmh to ',jj,xyzmh(1,lsendlist(1)+1),
     &                       iphase(lsendlist(1)+1)
#endif

                        CALL MPI_SEND(xyzmh,1,indexMPI5,jj,0,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(iphase,1,indexMPI_INT1,jj,1,
     &                       MPI_COMM_WORLD, ierr)

                        IF (iteration.EQ.1) THEN
                           CALL MPI_SEND(iuniquesend,nptmasssend,
     &                          MPI_INTEGER8,jj,9,MPI_COMM_WORLD, ierr)
                        ENDIF
                     ENDIF
                  END DO
                  inumber = 0
                  ibehind = ii-1
               ELSE
                  ibehind = ii-1
c
c--Receive data
c
                  CALL MPI_RECV(xyzmh(1,2*ntot+3),idim,i5REAL8,ibehind,
     &                 0,MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
                  IF (2*ntot+3+inumber.GE.mmax2) THEN
                 WRITE (*,*) iproc,': ERROR - 2*ntot+3+inumber.GE.mmax2'
                     CALL quit
                  ENDIF

#ifdef MPIDEBUG
                  print *,iproc,': Received ',inumber,' from ',ibehind,
     &              xyzmh(1,2*ntot+3),numproc
#endif

                  CALL MPI_RECV(iphase(ntot+1),idim,MPI_INTEGER1,
     &                 ibehind,1,MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_INTEGER1,igotnumber,
     &                 ierr)
                  IF (ntot+1+inumber.GE.idim2) THEN
                  WRITE (*,*) iproc,': ERROR - ntot+1+inumber.GE.idim2'
                     CALL quit
                  ENDIF
                  IF (igotnumber.NE.inumber) THEN
                     WRITE (iprint, *) 'ERROR - igotnumber.NE.inumber:',
     &                    ' iphase'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': Received iphase ',inumber,
     &                 ' from ',ibehind,iphase(ntot+1)
#endif

                  IF (iteration.EQ.1) THEN
                     CALL MPI_RECV(iuniquestore(nptmasslocal2+1), 
     &                 iptdim, MPI_INTEGER8,
     &                 ibehind, 9, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_INTEGER8,igotnumber,
     &                    ierr)
                     nptmasslocal2 = nptmasslocal2 + igotnumber
                     IF (nptmasslocal2.GT.iptdim) THEN
                        WRITE (*,*) 'ERROR - nptmasslocal2.GT.iptdim ',
     &                       nptmasslocal2, igotnumber
                        CALL quit
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
c
c--This calculates the contributions to rhoi, gradhi, numneighreal, and 
c     gravity forces from non-neighbours
c
            DO i = 1, inumber
               ipart = i + ntot
               IF (ipart.GT.idim2) THEN
                  WRITE (*,*) 'ipart.GT.idim2'
                  CALL quit
               ENDIF
               IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
                  IF (iteration.NE.1) THEN
                     WRITE (*,*) 'ERROR - iteration.NE.1'
                     CALL quit
                  ENDIF
                  nptmasslocal = nptmasslocal + 1
                  listpm(nptmasslocal) = ipart
                  listrealpm(ipart) = nptmasslocal
                  IF (iuniquestore(nptmasslocal).NE.
     &                 iuniquestoreold(nptmasslocal)) nptmasstotlast=-1

c                  WRITE (*,*) iproc,' Trans ',ipart,
c     &                 iuniquestore(nptmasslocal),
c     &                 iuniquestoreold(nptmasslocal),nptmasstotlast

                  iuniquestoreold(nptmasslocal) = 
     &                 iuniquestore(nptmasslocal)
               ENDIF
            END DO
            IF (iteration.EQ.1) THEN
               IF (nptmasslocal2.NE.nptmasslocal) THEN
                  WRITE (*,*) 'ERROR - nptmasslocal2.NE.nptmasslocal ',
     &                 nptmasslocal2,nptmasslocal
                  CALL quit
               ENDIF
            ENDIF

            inumreturn = 0
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumber,ntot,npart,xyzmh,dvxyzu,inumreturn,llistrec)
C$OMP& shared(nlstdo_tot)
C$OMP& private(i,ipart,nneighlocal)
            DO i = 1, inumber
               ipart = i + ntot
               IF (ipart.GT.idim2) THEN
                  WRITE (*,*) 'ipart.GT.idim2'
                  CALL quit
               ENDIF
               CALL densitygradh(ipart,npart,ntot,nlstdo_tot,xyzmh,
     &              nneighlocal,dvxyzu)
c
c--Make a list of indices to return densities
c
               IF (nneighlocal.GT.0) THEN
C$OMP CRITICAL (addtollistrec)
                  inumreturn = inumreturn + 1
                  llistrec(inumreturn) = i
C$OMP END CRITICAL (addtollistrec)
               ENDIF
            END DO
C$OMP END PARALLEL DO

            ntotplusinumber = ntot + inumber

#ifdef MPIDEBUG
c            IF (nptmasstot.GE.1) 
      IF (itime.EQ.2097152)
     &    print *,iproc,': number found neigh = ',inumreturn,numproc,
     &           inumber
#endif

            IF (iteration.EQ.1) THEN
c
c--nneightogetback does not record actual number, just 0 or 1.
c     Does not need the actual number, and getting the actual number would be
c     difficult because of the iterations to set h's and, thus, neighbours.
c
c--Also records whether any particles had neighbours during *ANY* of the
c     iterations, regardless of whether or not the particles had neighbours
c     on this process when the particle's h had eventually converged.
c
               IF (inumreturn.GT.0) THEN
                  inumofreturns = inumofreturns + 1
                  nneightogetback(ibehind+1) = 1
               ELSE
                  nneightogetback(ibehind+1) = 0
               ENDIF
            ELSE
               IF (inumreturn.GT.0 .AND. 
     &              nneightogetback(ibehind+1).EQ.0) THEN
                  inumofreturns = inumofreturns + 1
                  nneightogetback(ibehind+1) = 1
               ENDIF
            ENDIF
            IF (ideniteratebal) THEN
c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
#ifdef MPIDEBUG
               print *,iproc,': sending indices to ',ibehind,
     &              inumreturn,numproc
#endif

               CALL MPI_SENDRECV(llistrec, inumreturn, MPI_INTEGER,
     &           ibehind, 10, llistrec(inumreturn+1), nlstsend, 
     &           MPI_INTEGER, iahead, 10, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_INTEGER,inumberneigh,ierr)
               IF (inumreturn+inumberneigh.GT.ineighproc) THEN
               WRITE (*,*) 'ERROR - inumreturn+inumberneigh.GT.ineighp'
                  CALL quit
               ENDIF

#ifdef MPIDEBUG
               print *,iproc,': returned indices to ',ibehind,
     &              ' and got ',inumberneigh,' from ',iahead,numproc
#endif
c
c--Return densities of particles with density contributions (if any)
c
c--Note - need to change values of llistrec because in C array indices
c     start from 0 not 1 and MPI assumes C-type indexing
c
C$OMP PARALLEL DO SCHEDULE(static) default(none)
C$OMP& shared(inumreturn,llistrec,ntot)
C$OMP& private(j)
               DO j = 1, inumreturn
                  llistrec(j) = llistrec(j) + ntot - 1
               END DO
C$OMP END PARALLEL DO

c            CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1,
c     &           llistrec, MPI_REAL4, indexMPI1return, ierr)
               CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths,
     &              llistrec, MPI_REAL4, indexMPI1return, ierr)
               CALL MPI_TYPE_COMMIT(indexMPI1return,ierr)

c            CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1,
c     &           llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
               CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths,
     &              llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
               CALL MPI_TYPE_COMMIT(indexMPI_I1return,ierr)

               IF (inumberneigh.GT.ineighproc) THEN
                  WRITE (*,*) 'ERROR - inumberneigh.GT.ineighproc'
                  CALL quit
               ENDIF

               CALL MPI_SENDRECV(rho, 1, indexMPI1return, ibehind, 11,
     &              rhorec, inumberneigh, MPI_REAL4, iahead, 11,
     &              MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
               IF (icheck.GT.ineighproc) THEN
                  WRITE (*,*) 'ERROR - icheck.GT.ineighproc 1 ',icheck
                  CALL quit
               ENDIF
               IF (icheck.NE.inumberneigh) THEN
                  WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                 iproc
                  CALL quit
               ENDIF

#ifdef MPIDEBUG
               print *,iproc,': sent rho ',numproc
#endif

               CALL MPI_SENDRECV(dumrho, 1, indexMPI1return,ibehind,12,
     &           dumrhorec, inumberneigh, MPI_REAL4, iahead, 12,
     &           MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
               IF (icheck.NE.inumberneigh) THEN
                  WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                 iproc
                  CALL quit
               ENDIF

#ifdef MPIDEBUG
               print *,iproc,': sent gradh'
#endif

               CALL MPI_SENDRECV(nneigh,1,indexMPI_I1return,ibehind,13,
     &              nneighrec, inumberneigh, MPI_INTEGER, iahead, 13,
     &              MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_INTEGER,icheck,ierr)
               IF (icheck.NE.inumberneigh) THEN
                  WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                 iproc
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': sent nneigh'
#endif

               CALL MPI_TYPE_FREE(indexMPI_I1return,ierr)
               CALL MPI_TYPE_FREE(indexMPI1return,ierr)
               IF (iteration.EQ.1) THEN
                  IF (inumberneigh.GT.0) THEN
                     inumofsends = inumofsends + 1
                     nneighsentany(iahead+1) = .TRUE.
                  ENDIF
               ELSE
                  IF (inumberneigh.GT.0 .AND. 
     &                 .NOT.nneighsentany(iahead+1)) THEN
c     &              nneighsentback(iahead+1).EQ.0) THEN
                     inumofsends = inumofsends + 1
                     nneighsentany(iahead+1) = .TRUE.
                  ENDIF
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': all sent'
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumberneigh,lsendlist,llistrec,inumreturn)
C$OMP& shared(rho,rhorec,dumrho,dumrhorec,nneigh,nneighrec)
C$OMP& shared(iteration,llistsentback,iahead,nneighsentback)
C$OMP& private(k,ipos)
               DO k = 1, inumberneigh
                  ipos = lsendlist(llistrec(k + inumreturn)) + 1

                  rho(ipos) = rho(ipos) + rhorec(k)
                  dumrho(ipos) = dumrho(ipos) + dumrhorec(k)
                  nneigh(ipos) = nneigh(ipos) + nneighrec(k)
                  IF (iteration.EQ.1) THEN
                     llistsentback(k,iahead+1) = ipos - 1
                  ELSE
                     llistsentback(k+nneighsentback(iahead+1),
     &                    iahead+1) = ipos - 1
                  ENDIF
               END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
               print *,iproc,': received neighbour data ',numproc
#endif
               IF (iteration.EQ.1) THEN
#ifdef MPIDEBUG
                  print *,iproc,': reports ',inumberneigh,
     &                 ' foreign neighbours on process ',iahead+1
#endif
                  nneighsentbackold(iahead+1) = 0
                  nneighsentback(iahead+1) = inumberneigh
               ELSE
                  nneighsentbackold(iahead+1) = nneighsentback(iahead+1)
                  nneighsentback(iahead+1) = nneighsentback(iahead+1) +
     &                 inumberneigh
               ENDIF
               maxnneighsentback = MAX(maxnneighsentback,
     &              nneighsentback(iahead+1))
c
c--Send gravity forces back to sending process (includes forces from gas
c	particles on sinks if iptintree=1 and also sink particles if
c       iptintree=2).
c
#ifdef MPIDEBUG
               print *,iproc,' sending forces back to ',ibehind,numproc,
     &              ntot,inumber,ntotplusinumber,nlstsend,
     &              ntot+nlstsend,idim3
#endif
               IF (ntotplusinumber+nlstsend.GE.idim3) THEN
                  WRITE (*,*) iproc,
     &                 ': ERROR - ntotplusinumber+nlstsend.GE.idim3 ',
     &                 ntotplusinumber,nlstsend,idim3,npart,ntot
                  CALL quit
               ENDIF
               CALL MPI_SENDRECV(dvxyzu(1,ntot+1),4*inumber,
     &              MPI_REAL8,ibehind,15,dvxyzu(1,ntotplusinumber+1),
     &              4*nlstsend,MPI_REAL8, iahead, 15, 
     &              MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
               IF (ireturned.NE.4*nlstsend) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.4*nlstsend dvxyz'
                  CALL quit
               ENDIF

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstsend,lsendlist,ntotplusinumber,dvxyzu)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nlstsend
                  ipos = lsendlist(l) + 1
                  jpos = ntotplusinumber + l

                  DO k = 1, 4
                     dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
                  END DO
               END DO
C$OMP END PARALLEL DO
c
c--End circular send/receive
c
            ELSE
c
c--IMPLEMENTATION USING SEND
c
c
c--Receive data
c
               IF (ii.EQ.iproc+1) THEN
                  DO jj = 0 , numproc - 1
                     IF (iproc.NE.jj) THEN
                        CALL MPI_RECV(llistrec, nlstsend, 
     &                       MPI_INTEGER, MPI_ANY_SOURCE, 10, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                        CALL MPI_GET_COUNT(istatus,MPI_INTEGER,
     &                       inumberneigh,ierr)
                        iahead = istatus(MPI_SOURCE)
                        IF (inumberneigh.GT.ineighproc) THEN
                           WRITE (*,*) 'ERROR - ',
     &                          'inumberneigh.GT.ineighp'
                           CALL quit
                        ENDIF
#ifdef MPIDEBUG
                        print *,iproc,': got indices ',inumberneigh,
     &                       ' from ',iahead,numproc
#endif
                        CALL MPI_RECV(rhorec,inumberneigh,MPI_REAL4,
     &                       iahead,11,MPI_COMM_WORLD, istatus, ierr)
                        CALL MPI_GET_COUNT(istatus,MPI_REAL4,
     &                       icheck,ierr)
                        IF (icheck.GT.ineighproc) THEN
                           WRITE (*,*) 'ERROR - icheck.GT.',
     &                          'ineighproc 1 ',icheck
                           CALL quit
                        ENDIF
                        IF (icheck.NE.inumberneigh) THEN
                           WRITE (*,*) 'ERROR - icheck.NE.',
     &                          'inumberneigh 1 ',iproc
                           CALL quit
                        ENDIF
#ifdef MPIDEBUG
                        print *,iproc,': got rho ',icheck
#endif

                        CALL MPI_RECV(dumrhorec,inumberneigh,MPI_REAL4,
     &                       iahead,12,MPI_COMM_WORLD,istatus,ierr)
                        CALL MPI_GET_COUNT(istatus,MPI_REAL4,
     &                       icheck,ierr)
                        IF (icheck.NE.inumberneigh) THEN
                           WRITE (*,*) 'ERROR - icheck.NE.',
     &                          'inumberneigh 1 ',iproc
                           CALL quit
                        ENDIF
#ifdef MPIDEBUG
                        print *,iproc,': got gradh ',icheck
#endif

                        CALL MPI_RECV(nneighrec,inumberneigh, 
     &                       MPI_INTEGER,iahead,13,
     &                       MPI_COMM_WORLD,istatus,ierr)
#ifdef MPIDEBUG
                        print *,iproc,': rec '
#endif
                        CALL MPI_GET_COUNT(istatus,MPI_INTEGER,
     &                       icheck,ierr)
#ifdef MPIDEBUG
                        print *,iproc,': count '
#endif
                        IF (icheck.NE.inumberneigh) THEN
                           WRITE (*,*) 'ERROR - icheck.NE.',
     &                          'inumberneigh 1 ',iproc
                           CALL quit
                        ENDIF
#ifdef MPIDEBUG
                        print *,iproc,': sent nneigh'
#endif
                        IF (iteration.EQ.1) THEN
                           IF (inumberneigh.GT.0) THEN
                              inumofsends = inumofsends + 1
                              nneighsentany(iahead+1) = .TRUE.
                           ENDIF
                        ELSE
                           IF (inumberneigh.GT.0 .AND. 
     &                          .NOT.nneighsentany(iahead+1)) THEN
                              inumofsends = inumofsends + 1
                              nneighsentany(iahead+1) = .TRUE.
                           ENDIF
                        ENDIF
#ifdef MPIDEBUG
                        print *,iproc,': all received'
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumberneigh,lsendlist,llistrec)
C$OMP& shared(rho,rhorec,dumrho,dumrhorec,nneigh,nneighrec)
C$OMP& shared(iteration,llistsentback,iahead,nneighsentback)
C$OMP& private(k,ipos)
                        DO k = 1, inumberneigh
                           ipos = lsendlist(llistrec(k))+1

                           rho(ipos) = rho(ipos) + rhorec(k)
                           dumrho(ipos) = dumrho(ipos) + dumrhorec(k)
                           nneigh(ipos) = nneigh(ipos) + nneighrec(k)
                           IF (iteration.EQ.1) THEN
                              llistsentback(k,iahead+1) = ipos - 1
                           ELSE
                              llistsentback(k+nneighsentback(iahead+1),
     &                             iahead+1) = ipos - 1
                           ENDIF
                        END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
                        print *,iproc,': received neighbour data'
#endif
                        IF (iteration.EQ.1) THEN
#ifdef MPIDEBUG
                           print *,iproc,': reports ',inumberneigh,
     &                          ' foreign neighbours on process ',
     &                          iahead
#endif
                           nneighsentbackold(iahead+1) = 0
                           nneighsentback(iahead+1) = inumberneigh
                        ELSE
                           nneighsentbackold(iahead+1) = 
     &                          nneighsentback(iahead+1)
                           nneighsentback(iahead+1) = 
     &                          nneighsentback(iahead+1) +
     &                          inumberneigh
                        ENDIF
                        maxnneighsentback = MAX(maxnneighsentback,
     &                       nneighsentback(iahead+1))
c
c--Receive gravity forces back (includes forces from gas
c	particles on sinks if iptintree=1 and also sink particles if
c       iptintree=2).
c
                        IF (ntotplusinumber+nlstsend.GE.idim3) THEN
                           WRITE (*,*) iproc,
     &                  ': ERROR - ntotplusinumber+nlstsend.GE.idim ',
     &                  ntotplusinumber,nlstsend,idim,npart,ntot
                           CALL quit
                        ENDIF
                        CALL MPI_RECV(dvxyzu(1,ntotplusinumber+1),
     &                       4*nlstsend, MPI_REAL8, iahead, 15, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                        CALL MPI_GET_COUNT(istatus, MPI_REAL8,
     &                       ireturned, ierr)
                        IF (ireturned.NE.4*nlstsend) THEN
                           WRITE (*,*) 'ERROR - ireturned.NE.',
     &                          '4*nlstsend dvxyz'
                           CALL quit
                        ENDIF

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstsend,lsendlist,ntotplusinumber,dvxyzu)
C$OMP& private(l,ipos,jpos)
                        DO l = 1, nlstsend
                           ipos = lsendlist(l) + 1
                           jpos = ntotplusinumber + l

                           DO k = 1, 4
                              dvxyzu(k,ipos) = dvxyzu(k,ipos) + 
     &                             dvxyzu(k,jpos)
                           END DO
                        END DO
C$OMP END PARALLEL DO
                     ENDIF
                  END DO
               ELSE
c
c--Send data
c
                  ibehind = ii - 1

                  CALL MPI_SEND(llistrec,inumreturn,MPI_INTEGER,ibehind,
     &                 10, MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent indices to ',ibehind
#endif
c
c--Return densities of particles with density contributions (if any)
c
c--Note - need to change values of llistrec because in C array indices
c     start from 0 not 1 and MPI assumes C-type indexing
c
C$OMP PARALLEL DO SCHEDULE(static) default(none)
C$OMP& shared(inumreturn,llistrec,ntot)
C$OMP& private(j)
                  DO j = 1, inumreturn
                     llistrec(j) = llistrec(j) + ntot - 1
                  END DO
C$OMP END PARALLEL DO

c            CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1,
c     &           llistrec, MPI_REAL4, indexMPI1return, ierr)
                  CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths,
     &                 llistrec, MPI_REAL4, indexMPI1return, ierr)
                  CALL MPI_TYPE_COMMIT(indexMPI1return,ierr)

c            CALL MPI_TYPE_CREATE_INDEXED_BLOCK(inumreturn, 1,
c     &           llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
                  CALL MPI_TYPE_INDEXED(inumreturn, lblocklengths,
     &                 llistrec,MPI_INTEGER,indexMPI_I1return,ierr)
                  CALL MPI_TYPE_COMMIT(indexMPI_I1return,ierr)

                  CALL MPI_SEND(rho, 1, indexMPI1return, ibehind, 11,
     &                 MPI_COMM_WORLD, istatus, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent rho to ',ibehind
#endif
                  CALL MPI_SEND(dumrho, 1, indexMPI1return, ibehind, 12,
     &                 MPI_COMM_WORLD, istatus, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent gradh to ',ibehind,inumreturn,
     &                 llistrec(1)
#endif
                  CALL MPI_SEND(nneigh, 1 ,indexMPI_I1return, 
     &                 ibehind, 13, MPI_COMM_WORLD, istatus, ierr)
#ifdef MPIDEBUG
                  print *,iproc,': sent nneigh to ',ibehind
#endif
                  CALL MPI_TYPE_FREE(indexMPI_I1return,ierr)
                  CALL MPI_TYPE_FREE(indexMPI1return,ierr)
#ifdef MPIDEBUG
                  print *,iproc,': all sent to ',ibehind
#endif
c
c--Send gravity forces back to sending process (includes forces from gas
c       particles on sinks if iptintree=1 and also sink particles if
c       iptintree=2).
c
#ifdef MPIDEBUG
                  print *,iproc,' sending forces back to ',ibehind,
     &                 numproc,ntot,nlstsend,ntot+nlstsend,idim
#endif
                  CALL MPI_SEND(dvxyzu(1,ntot+1),4*inumber,
     &                 MPI_REAL8,ibehind,15,
     &                 MPI_COMM_WORLD, istatus, ierr)
               ENDIF
c
c--End send/receive of data
c
            ENDIF
c
c--END OF FOREIGN CONTRIBUTIONS
c
         END DO
 444     CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
         CALL MPI_TYPE_FREE(indexMPI5,ierr)
         IF (individualtimesteps.NE.2) THEN
            IF (iteration.EQ.1) THEN
               IF (nptmasslocal2.NE.nptmasstot) THEN
                  WRITE (*,*) 'ERROR - nptmasslocal2.NE.nptmasstot ',
     &                 nptmasslocal2, nptmasstot
                  CALL quit
               ENDIF
            ENDIF
         ENDIF
c
c--End of MPI-only section
c
#endif

c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      IF (itime.EQ.000)
c         print *,iproc,': DOING LOCAL'
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

         numsucceeded = 0
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstdo_tot,npart,ntot,listp,divv,curlv,gradhs)
C$OMP& shared(xyzmh,vxyzu,pr,vsound,rho,ekcle,rhomin)
C$OMP& shared(nneigh,selfnormkernel)
C$OMP& shared(cnormk,radkernel,dvtable,ddvtable,wij,grwij)
C$OMP& shared(listpm,iphase,dphidh,uradconst,icall,encal)
C$OMP& shared(iprint,nptmass,iptmass,iorig,third)
C$OMP& shared(dumrho,iscurrent,ibound,ireal)
C$OMP& shared(isteps,it0,it1,imax,imaxstep,dt,itime)
C$OMP& shared(varmhd,Bevol,Bxyz,vsmooth,Bsmooth,nlstdo,ldolist)
C$OMP& shared(dvxyzu,hi_old,isucceed,poten,isend,numproc,iproc,gradmhd)
C$OMP& shared(ivar,ijvar,ncompact,icompact,iteration)
C$OMP& shared(hminbisec,hmaxbisec,iupdated,nlistupdated,listparents)
C$OMP& private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi)
C$OMP& private(divvi,curlvxi,curlvyi,curlvzi,gradhi,gradsofti)
C$OMP& private(pmassj,hi_oldi,hi1,hi21,hi31,hi41,hneigh)
C$OMP& private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2,rcut)
C$OMP& private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij)
C$OMP& private(projv,procurlvx,procurlvy,procurlvz)
C$OMP& private(l,iptcur,dphi,dwdhi,dpotdh,numneighi)
C$OMP& private(numneighreal,rhohi,dhdrhoi,omegai,func,dfdh1)
C$OMP& private(hnew,deltarho)
C$OMP& private(nneighlocal,nneighi,ncompacthere,icompacthere)
C$OMP& private(numberupdated)
C$OMP& reduction(MAX:rhonext,imaxit)
C$OMP& reduction(+:inumit,inumfixed,inumrecalc,nwarnup,nwarndown)
C$OMP& reduction(+:nwarnroundoff,nbisection,numsucceeded)
      DO n = 1, nlstdo
         ipart = ldolist(n)
c
c--Adds on local contributions to rho, gradh, and nneigh and long-range
c	 gravity
c
         CALL densitygradh(ipart,npart,ntot,nlstdo_tot,xyzmh,
     &        nneighlocal,dvxyzu)

         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN

         IF (nneighlocal.GT.nneighmax) THEN
            WRITE (iprint,*) 'ERROR: nneighlocal exceeds nneighmax'
            WRITE (iprint,*) iorig(ipart),nneighlocal,nneighmax
            CALL quit
         ENDIF
c
c--gradhi value is stored in dumrho temporarily
c
         gradhi = dumrho(ipart)
         rhoi = rho(ipart)

         pmassi = xyzmh(4,ipart)
         hi = xyzmh(5,ipart)
         hi1 = 1.0/hi
         hi_oldi = hi_old(ipart)
         nneighi = nneigh(ipart)

c
c--Iteration business
c  These lines define the relationship between h and rho
c     omega is the term in the denominator as in Monaghan(2001)
c
         rhohi = pmassi*(hfact*hi1)**3 - rhomin
         dhdrhoi = -hi/(3.*(rhoi + rhomin))
         omegai = 1. - dhdrhoi*gradhi
         func = rhohi - rhoi

         IF (isucceed(ipart).NE.-1) THEN
c
c--Newton-Raphson iteration
c
            dfdh1 = dhdrhoi/omegai
            hnew = hi - func*dfdh1
c               write(iprint,*)
c     &          'newton raphson (',iteration,'): hnew = ',hnew
         ELSE
c
c--Bisection iteration
c
            IF (func.lt. 0.) THEN
               hmaxbisec(ipart) = hi
            ELSE
               hminbisec(ipart) = hi
            ENDIF
            hnew = 0.5*(hminbisec(ipart) + hmaxbisec(ipart))
            write(iprint,*) 'bisection (',iteration,'): hnew =',
     &           hnew,ipart,hi,func,hmaxbisec(ipart),hminbisec(ipart)
c            write(*,*) 'bisection (',iteration,'): hnew =',hnew,ipart,
c     &           hi,func,hmaxbisec(ipart),hminbisec(ipart),nneighi,
c     &           gradhi,rhoi,rhohi,dhdrhoi,omegai
         ENDIF
c
c--Don't allow sudden jumps to huge numbers of neighbours
c  (Newton-Raphson only)
c
         IF (isucceed(ipart).NE.-1) THEN
            IF (hnew.GT.1.2*hi) THEN
               nwarnup = nwarnup + 1
               hnew = 1.2*hi
            ELSEIF (hnew.LT.0.8*hi) THEN
               nwarndown = nwarndown + 1
               hnew = 0.8*hi
            ENDIF
         ENDIF

c         IF (itime.GE.142868480)
c     &  write(*,*)'H-stuff ',hnew,hi,hi_oldi,omegai,nneighi,
c     &        isucceed(ipart),iphase(ipart),ipart

         IF ((hnew.LE.0. .OR. omegai.LE.tiny .OR. 
     &        iteration.EQ.maxitsnr .OR. nneighi.LE.0) 
     &        .AND. isucceed(ipart).NE.-1) THEN
c
c--Switch to Bisection if not converging or running into trouble
c
            WRITE(iprint,*) 'WARNING: switching to bisection on'//
     &           ' particle ',iorig(ipart),ipart,'(',ipart,') hi = ',hi,
     &           ' hnew = ',hnew
c            WRITE(*,*) 'WARNING: switching to bisection on'//
c     &           ' particle ',iorig(ipart),'(',ipart,') hi = ',hi,
c     &           ' hnew = ',hnew,nneighi,gradhi,rhoi,rhohi,dhdrhoi,
c     &           omegai,iteration
            IF (nneighi.LE.0) THEN
               WRITE(iprint,*) '(particle has no neighbours) ',nneighi
            ENDIF
            IF (omegai.LE.tiny) THEN
               WRITE(iprint,*) '(omega < tiny) ',omegai
            ENDIF
            IF (iteration.EQ.maxitsnr) THEN
               WRITE(iprint,*) '(more than ',maxitsnr,' iterations)'
            ENDIF

            nbisection = nbisection + 1
            isucceed(ipart) = -1
            hminbisec(ipart) = 0.
            hmaxbisec(ipart) = 1.e6
!--Don't have to start with ridiculous h, 
!  just something reasonably big between above limits
            hnew = 2.*hi_oldi
!               hnew = 0.5*(hminbisec(ipart) + hmaxbisec(ipart))
c
c--Otherwise check for convergence
c
         ELSEIF (ABS((hnew-hi)/hi_oldi).LT.htol.AND.isucceed(ipart).EQ.0
     &           .AND. ABS((hnew-hi)/hi).LT.htol
     &           .OR. ABS((hnew-hi)/hi_oldi).LT.1.0E-6
     &           .AND. ABS((hnew-hi)/hi).LT.htol
     &           ) THEN

            imaxit = MAX(imaxit, iteration)
            inumit = inumit + iteration
            IF (isucceed(ipart).EQ.-1) THEN
               write(iprint,*) 'SWITCH BACK TO N-R ',ipart,hnew,hi
c               write(*,*) 'SWITCH BACK TO N-R ',ipart,hnew,hi
               isucceed(ipart) = 0  ! switch back to N-R
            ELSE
c
c--Store quantities if converged
c
               isucceed(ipart) = 1
               numsucceeded = numsucceeded + 1

               rho(ipart) = rhoi
               dumrho(ipart) = rhoi

               xyzmh(5,ipart) = hi
               gradhs(1,ipart) = 1./omegai

               poten(ipart) = dvxyzu(4,ipart)

               IF (iphase(ipart).EQ.0) THEN
c
c--Set e(i) for the first time around because density only set here
c     This sets radiation and matter to have same initial temperature
c
                  IF (icall.EQ.1 .AND. encal.EQ.'r') THEN
                     IF (ekcle(1,ipart).EQ.0.0) THEN
                       ekcle(3,ipart) = getcv(rho(ipart),vxyzu(4,ipart))
                        ekcle(1,ipart) = uradconst*(vxyzu(4,ipart)/
     &                      ekcle(3,ipart))**4/rho(ipart)
                        ekcle(2,ipart) = getkappa(vxyzu(4,ipart),
     &                       ekcle(3,ipart),rho(ipart))
                     ENDIF
                  ENDIF
c
c--Pressure and sound velocity from ideal gas law...
c
                  CALL eospg(ipart, vxyzu, rho, pr, vsound, ekcle)
               ENDIF
c
c--Calculate other quantities using the interaction list from local process
c
c               CALL divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
c     &              nneighlocal,neighlist,gradmhd)
c
c--Set up compact list of neighbours within h_i and non-active neighbours
c     within h_j
c
C$OMP CRITICAL (listcompact)
               ncompact = ncompact + 1
               ncompacthere = ncompact
               icompacthere = icompact
               icompact = icompact + nneighlocal
C$OMP END CRITICAL (listcompact)
               IF (icompacthere+nneighlocal.GT.icompactmax) THEN
                  WRITE (*,*) 'ERROR - compact not large enough DI ',
     &                 icompacthere, nneighlocal, icompactmax
                  CALL quit
               ENDIF
               ivar(1,ncompacthere) = nneighlocal
               ivar(2,ncompacthere) = icompacthere
               ivar(3,ncompacthere) = ipart
               DO k = 1, nneighlocal
                  j = neighlist(k)
                  ijvar(icompacthere + k) = j
                  IF (.NOT.iscurrent(j) .AND. .NOT.iupdated(j) 
     &                 .AND. j.LE.npart) THEN
                     iupdated(j) = .TRUE.
C$OMP FLUSH(iupdated)
C$OMP CRITICAL (listupdated)
                     nlistupdated = nlistupdated + 1
                     numberupdated = nlistupdated
C$OMP END CRITICAL (listupdated)
                     IF (numberupdated.GT.idim) THEN
                        WRITE (iprint,*) 'ERROR - numberupdated.GT.idim'
                        CALL quit
                     ENDIF
                     listparents(numberupdated) = j
                  ENDIF
               END DO

            ENDIF
         ENDIF

         ELSE
c
c--Else for sink particles
c
         isucceed(ipart) = 1
         numsucceeded = numsucceeded + 1
         poten(ipart) = dvxyzu(4,ipart)
c
c--Set up compact list of neighbours within h_i and non-active neighbours
c     within h_j
c
C$OMP CRITICAL (listcompact)
         ncompact = ncompact + 1
         ncompacthere = ncompact
C$OMP END CRITICAL (listcompact)
         ivar(1,ncompacthere) = 0
         ivar(2,ncompacthere) = 1
         ivar(3,ncompacthere) = ipart
c
c--End of sink particles
c
         ENDIF

#ifdef MPI
         IF (isucceed(ipart).EQ.1 .OR. 
     &        (hnew.LT.hi .AND. nneigh(ipart).EQ.0)) THEN
            isend(ipart) = .FALSE.
         ELSE
            isend(ipart) = .TRUE.
         ENDIF
#endif

c         hi_old(ipart) = hi
         IF (isucceed(ipart).NE.1) xyzmh(5,ipart) = hnew

      END DO
C$OMP END PARALLEL DO
c
c--Check to see whether all particles have converged
c
#ifdef MPIDEBUG
c      IF (itime.EQ.000)
           print *,'numsucceeded ',numsucceeded,nlstdo,htol,numproc
#endif

      IF (numsucceeded.EQ.nlstdo) THEN
#ifdef MPI
         iscore = 1
#else
         GOTO 30
#endif
#ifdef MPI
      ELSE
         iscore = 0
#endif
      ENDIF
#ifdef MPI
c
c--Test other processes to see if they have finished also
c
      CALL MPI_ALLREDUCE(iscore,ifinished,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)

#ifdef MPIDEBUG
c      IF (itime.EQ.000) THEN
      print *,iproc,': ifinished is ',ifinished,iscore,numproc
      print *,iproc,': maxnneighsentback is ',maxnneighsentback

      print *,iproc,': INFORMATION ',iteration
      DO i = 1, numproc
         print *,iproc,': ',i,nneighsentback(i),nneighsentbackold(i),
     &        nneighsentany(i)
      END DO
      print *,iproc,': inumofsends ',inumofsends
c      ENDIF
#endif

      IF (ifinished.EQ.numproc) GOTO 30
#endif
c
c--Otherwise, make new list of particles that are still to converge
c
      nlstsend = 0
      nlstdo_old = nlstdo
      nlstdo = 0
      DO n = 1, nlstdo_old
         ipart = ldolist(n)
         IF (isucceed(ipart).NE.1) THEN
            nlstdo = nlstdo + 1
            ldolist(nlstdo) = ipart
#ifdef MPI
            IF (isend(ipart)) THEN
               nlstsend = nlstsend + 1
               lsendlist(nlstsend) = ipart - 1
            ENDIF
#endif
         ENDIF
      END DO
c
c--Re-zero quantities
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstdo,ldolist,rho,dumrho,divv,curlv,gradhs,nneigh)
C$OMP& shared(dvxyzu,poten)
C$OMP& private(n,ipart,k)
      DO n = 1, nlstdo
         ipart = ldolist(n)

         rho(ipart) = 0.
         dumrho(ipart) = 0.
         divv(ipart) = 0.
         curlv(ipart) = 0.
         gradhs(2,ipart) = 0.
         nneigh(ipart) = 0
         DO k = 1, 4
            dvxyzu(k,ipart) = 0.
         END DO
         poten(ipart) = 0.
      END DO
C$OMP END PARALLEL DO

#ifdef MPI
c
c--Need to remove those particles that are being sent back now from the lists
c     of particles that will be sent back to each MPI process to calculate
c     forces since they will be added back onto those lists if they 
c     have neighbours anyway (and we don't want them done more than once!)
c
      DO i = 1, numproc
         IF (iproc.NE.i-1) THEN
#ifdef MPIDEBUG
            print *,iproc,': Trying remove ',nneighsentback(i),
     &           nneighsentbackold(i),nlstsend,nlstdo
#endif
            IF (nneighsentback(i).GE.nneighsentbackold(i)+1) THEN
               newvalue = nneighsentbackold(i)
               DO k = nneighsentbackold(i)+1, nneighsentback(i)
                  ifound = 0
                  DO n = 1, nlstsend
                     IF (llistsentback(k,i).EQ.lsendlist(n)) THEN
                        ifound = 1
                        GOTO 29
                     ENDIF
                  END DO
 29               IF (ifound.EQ.0) THEN
                     newvalue = newvalue + 1
                     llistsentback(newvalue,i) = llistsentback(k,i)
                  ENDIF
               END DO
               nneighsentback(i) = newvalue
c               IF (nneighsentback(i).EQ.0) inumofsends = inumofsends-1
            ENDIF
         ENDIF
      END DO

#ifdef MPIDEBUG
      print *,iproc,': REMOVED neighbours being sent back ',iteration
      DO i = 1, numproc
         print *,iproc,': ',i,nneighsentback(i),nneighsentbackold(i)
      END DO
#endif

#endif
c
c--End of iteration loop
c
      END DO
      WRITE (iprint,*) 'ERROR: density iteration failed for ',nlstdo,
     &     ' particles ',ifinished,iproc
      WRITE (*,*) 'ERROR: density iteration failed for ',nlstdo,
     &     ' particles ',ifinished,iproc
#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      CALL quit

 30   CONTINUE
      ncompactlocal = ncompact
      nptmasstotlast = nptmasstot

#ifdef MPI
      IF (itiming) THEN
         CALL getused(td12)
         td1 = td1 + (td12 - td11)
      ENDIF
#endif
c
c--FINISHED CALCULATION OF NEIGHBOURS, DENSITY, AND LONG-RANGE GRAVITY
c
c
c--Now calculated local contributions to divv, curlv, etc
c
      IF (itiming) CALL getused(td31)

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end)
C$OMP& shared(xyzmh,vxyzu,pr,vsound,rho,ekcle,Bevol)
C$OMP& shared(iphase,uradconst,icall,encal)
C$OMP& shared(numproc,iproc,gradmhd,npart,ntot)
C$OMP& shared(ivar,ijvar)
C$OMP& private(n,ipart,k)
C$OMP& private(nneighlocal)
      DO n = nlst_in, nlst_end
         ipart = ivar(3,n)
         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN
c
c--Calculate other quantities using the interaction list from local process
c
            nneighlocal = ivar(1,n)

            CALL divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
     &           nneighlocal,ijvar(ivar(2,n)+1),gradmhd)
         ENDIF
      END DO
C$OMP END PARALLEL DO

      IF (itiming) THEN
         CALL getused(td32)
         td3 = td3 + (td32 - td31)
      ENDIF
c
c--For MPI job, need to add contributions to other quantities (such as divv,
c     curlv, etc) from other MPI processes
c
#ifdef MPI

      IF (itiming) CALL getused(td21)

#ifdef MPIDEBUG
c      IF (itime.EQ.000) THEN
         print *,iproc,': HERE '
         DO i = 1, numproc
         print *,iproc,': VALUES ',i-1,nneightogetback(i),inumofreturns
         END DO
         DO i = 1, numproc
         print *,iproc,': SEND ',i-1,nneighsentback(i),maxnneighsentback
         END DO
c      ENDIF
#endif

      nneighsentanyatall = .FALSE.
      DO i = 1, numproc
         nneighsentanyatall = nneighsentanyatall .OR. nneighsentany(i)
      ENDDO
      inumbertotal = 0
      inumberreturned = 0
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
c      IF (itime.EQ.000)
                       print *,iproc,': sending neighbour data to ',j,
     &                       ' starting at ',llistsentback(1,j+1),'+1',
     &                       ' and sending ',nneighsentback(j+1),
     &                       ' elements'
#endif
c
c--Else does not need to send any particle back to this particular process.
c
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c                       1,llistsentback(1,j+1),i5REAL8,indexMPI5,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i5REAL8,indexMPI5,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

c                        IF (iproc.EQ.0) THEN
c                           print *,' Sent back ',llistsentback(1,j+1),
c     &                       iunique(iorig(llistsentback(1,j+1)+1)),
c     &                          i,j
c                        ENDIF

                        CALL MPI_SEND(xyzmh,1,indexMPI5,j,40,
     &                       MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
                        print *,iproc,' sent xyzmh to ',j
#endif
                        CALL MPI_TYPE_FREE(indexMPI5,ierr)

                        IF (nneighsentback(j+1).GT.0) THEN

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i4REAL8,indexMPI4,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

                           CALL MPI_SEND(vxyzu,1,indexMPI4,j,41,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent vxyzu'
#endif
                           CALL MPI_TYPE_FREE(indexMPI4,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

                           CALL MPI_SEND(dumrho,1,indexMPI1,j,42,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent dumrho'
#endif
                           CALL MPI_TYPE_FREE(indexMPI1,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_INTEGER1,
c     &                  indexMPI_INT1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                          llistsentback(1,j+1),MPI_INTEGER1,
     &                          indexMPI_INT1,ierr)
                           CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)

                           CALL MPI_SEND(iphase,1,indexMPI_INT1,j,43,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent iphase'
#endif
                           CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)

                           IF (imhd.EQ.idim) THEN
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)

                              CALL MPI_SEND(Bevol,1,indexMPI3,j,44,
     &                             MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                              print *,iproc,' sent Bevol'
#endif
                              CALL MPI_TYPE_FREE(indexMPI3,ierr)
                           ENDIF
                        ENDIF

#ifdef MPIDEBUG
c      IF (itime.EQ.000)
                       print *,iproc,': sent neighbour data to ',j,
     &                       radkernel,xyzmh(1,1)
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
c      IF (itime.EQ.000)
            print *,iproc,': receiving neigh data ',inumberreturned,
     &              inumofreturns,radkernel,xyzmh(1,1)
#endif
               istart = ntot + inumbertotal + 1
               CALL MPI_RECV(xyzmh(1,istart+ntot+2), idim,i5REAL8,
     &              i, 40, MPI_COMM_WORLD, istatus, ierr)
c     &              MPI_ANY_SOURCE, 40, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
               IF (istart+ntot+2+inumber.GT.mmax2) THEN
                  WRITE (*,*) 'ERROR - istart+ntot+2+inumber.GT.mmax2'
                  CALL quit
               ENDIF
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
c      IF (itime.EQ.000) THEN
               print *,iproc,': got xyzmh from ',iprocrec,inumber,
     &     ' to be put into ',istart,radkernel,xyzmh(1,1)

               print *,iproc,': POSITION BACK ',xyzmh(1,istart+ntot+2),
     &              xyzmh(2,istart+ntot+2),xyzmh(3,istart+ntot+2)
c      ENDIF
#endif
               IF (istart+inumber.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istart+inumber.GT.idim2 ',
     &                 istart,inumber
                  CALL quit
               ENDIF

               IF (inumber.GT.0) THEN
                  CALL MPI_RECV(vxyzu(1,istart),inumber,i4REAL8,
     &                 iprocrec, 41, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i4REAL8, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 1 ',
     &                    iproc
                     CALL quit
                  ENDIF

                  CALL MPI_RECV(dumrho(istart),inumber,MPI_REAL4,
     &                 iprocrec, 42, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 2 ',
     &                    iproc
                     CALL quit
                  ENDIF

                  CALL MPI_RECV(iphase(istart),inumber,MPI_INTEGER1,
     &                 iprocrec, 43, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER1,icheck,ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 3 ',
     &                    iproc
                     CALL quit
                  ENDIF

                  IF (imhd.EQ.idim) THEN
                     CALL MPI_RECV(Bevol(1,istart),inumber,i3REAL8,
     &                    iprocrec, 44, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, i3REAL8, icheck, ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 4 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got Bevol from ',iprocrec,icheck,
     &                    radkernel,xyzmh(1,1)
#endif
                  ENDIF
               ENDIF
               inumbertotal = inumbertotal + inumber
               inumberindiv(inumberreturned) = inumber
               inumbercumm(inumberreturned) = istart
               inumberproc(inumberreturned) = iprocrec
#ifdef MPIDEBUG
               print *,iproc,' set numbers ',radkernel,xyzmh(1,1)
#endif
            ENDIF
         ENDIF
      END DO
#ifdef MPIDEBUG
      print *,iproc,' inumbertotal is ',inumbertotal
#endif
      IF (inumbertotal.GT.idim) THEN
         WRITE (*,*) 'ERROR - inumbertotal.GT.idim ', iproc,
     &        inumbertotal,idim,(inumberindiv(ix),ix=1,8)
         CALL quit
      ENDIF
      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR1 - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit
      ENDIF

#ifdef MPIDEBUG
c      IF (itime.EQ.000) THEN
      print *,iproc,': Received all particles with remote neighbours ',
     &     inumbertotal,radkernel,xyzmh(1,1)
      print *,' '
      print *,' '
c      ENDIF
#endif
c
c--Calculate ONLY neighbours lists for remote particles, then contributions
c     to divv, curlv, and Euler potential alpha and beta, etc
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumbertotal,ntot,listp,nlst_end,iphase,hasghost)
C$OMP& shared(npart,xyzmh,acc,nneigh,vxyzu,Bevol,gradmhd)
C$OMP& shared(ncompact,icompact,ivar,ijvar,iscurrent)
C$OMP& shared(iupdated,nlistupdated,listparents,iprint)
C$OMP& private(fsx,fsy,fsz,epot,nneighlocal,ncompacthere,icompacthere)
C$OMP& private(i,k,ipart,jpart,numberupdated)
      DO i = 1, inumbertotal

         ipart = ntot + i
         listp(nlst_end + i) = ipart
         hasghost(ipart) = .FALSE.

         CALL treef(ipart,npart,ntot,xyzmh,acc,0,fsx,fsy,fsz,epot)

         nneighlocal = nneigh(ipart)

         CALL divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
     &        nneighlocal,neighlist,gradmhd)
c
c--Set up compact list of neighbours within h_i and non-active neighbours
c     within h_j
c
C$OMP CRITICAL (listcompact2)
         ncompact = ncompact + 1
         ncompacthere = ncompact
         icompacthere = icompact
         icompact = icompact + nneighlocal
C$OMP END CRITICAL (listcompact2)
         IF (icompacthere+nneighlocal.GT.icompactmax) THEN
            WRITE (*,*) 'ERROR - compact not large enough DI ',
     &           icompacthere, nneighlocal, icompactmax
            CALL quit
         ENDIF
         ivar(1,ncompacthere) = nneighlocal
         ivar(2,ncompacthere) = icompacthere
         ivar(3,ncompacthere) = ipart
         DO k = 1, nneighlocal
            jpart = neighlist(k)
            ijvar(icompacthere + k) = jpart
            IF (.NOT.iscurrent(jpart) .AND. .NOT.iupdated(jpart) 
     &           .AND. jpart.LE.npart) THEN
               iupdated(jpart) = .TRUE.
C$OMP FLUSH(iupdated)
C$OMP CRITICAL (listupdated)
               nlistupdated = nlistupdated + 1
               numberupdated = nlistupdated
C$OMP END CRITICAL (listupdated)
               IF (numberupdated.GT.idim) THEN
                  WRITE (iprint,*) 'ERROR - numberupdated.GT.idim'
                  CALL quit
               ENDIF
               listparents(numberupdated) = jpart
            ENDIF
         END DO
      END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
      print *,iproc,': has calculated all divvcurlgrad ',inumbertotal
c         IF (iproc.EQ.2) print *,iproc,':G',dvxyzu(1,ntot + 1+2054-1),
c     &              dvxyzu(2,ntot+1+2054-1),dvxyzu(3,ntot+1+2054-1),
c     &     ntot+1+2054-1,ntot
#endif

c      STOP

c      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
c--Now need to return contributions from neighbours on remote nodes.
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
c      IF (itime.EQ.000)
                    print *,iproc,': expecting to get ',inumofsends,
     &              ' divv contributions ',j,
     &              ' with maxnneighsentback ',maxnneighsentback
c               IF (iproc.EQ.14) THEN
c                  DO iii = 0, numproc -1
c                     print *,iproc,': LIST ',iii,nneighsentback(iii+1)
c                  END DO
c               ENDIF
#endif
               CALL MPI_RECV(divv(istartrec),maxnneighsentback,
     &              MPI_REAL4,MPI_ANY_SOURCE,50,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (istartrec+ireturned.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istartrec+ireturned.GT.idim2 ',
     &                 istartrec,ireturned
                  CALL quit
               ENDIF
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsentback divv ',
     &                 ireturned,nneighsentback(iprocrec+1)
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got divv from ',iprocrec,' put into ',
     &              llistsentback(1,iprocrec+1)+1
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,divv,istartrec)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  divv(ipos) = divv(ipos) + divv(jpos)
               END DO
C$OMP END PARALLEL DO

               CALL MPI_RECV(curlv(istartrec),maxnneighsentback,
     &              MPI_REAL4,iprocrec,51,MPI_COMM_WORLD,istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsendback curlv'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got curlv from ',iprocrec,radkernel,
     &              xyzmh(1,1)
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,curlv,istartrec)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  curlv(ipos) = curlv(ipos) + curlv(jpos)
               END DO
C$OMP END PARALLEL DO

               CALL MPI_RECV(gradhs(1,istartrec),2*maxnneighsentback,
     &              MPI_REAL4,iprocrec,52,MPI_COMM_WORLD,istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (ireturned.NE.2*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsendback gradhs'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
c      IF (itime.EQ.000)
                    print *,iproc,': got gradhs2 from ',iprocrec
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,gradhs,istartrec)
C$OMP& private(l,ipos,jpos)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  gradhs(2,ipos) = gradhs(2,ipos) + gradhs(2,jpos)
               END DO
C$OMP END PARALLEL DO

c
c--Receive back quantities for MHD gradients (if required)
c
               IF (imhd.EQ.idim) THEN

                  CALL MPI_RECV(gradmhd(1,istartrec),
     &                 15*maxnneighsentback,MPI_REAL8,
     &                 iprocrec,53,MPI_COMM_WORLD,istatus,ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
                  IF (ireturned.NE.15*nneighsentback(iprocrec+1)) THEN
                     WRITE (*,*) 'ERROR - ireturned.NE.nnsendbk gradmhd'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': got gradMHD from ',iprocrec
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,gradmhd,istartrec)
C$OMP& private(l,ipos,jpos)
                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1
                     DO k = 1, 15
                        gradmhd(k,ipos)=gradmhd(k,ipos)+gradmhd(k,jpos)
                     END DO
                  END DO
C$OMP END PARALLEL DO
               ENDIF

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
c      IF (itime.EQ.000)
                   print *,iproc,': sending forces to ',iprocsend,
     &                 ' istart ',istart,' number ',inumber,
     &                 ' returned ',inumberreturned,
     &                 ' of ',inumofreturns,xyzmh(1,1)
#endif
c                  iloop = iloop + 1
c                  print *,'SETTING iloop ',iloop
c                  IF (iloop.EQ.2 .and. inumber.EQ.5360) THEN
c                     inumber = 1
c                  ENDIF

                  CALL MPI_SEND(divv(istart),inumber,MPI_REAL4,
     &                 iprocsend,50,MPI_COMM_WORLD, ierr)
                  CALL MPI_SEND(curlv(istart),inumber,MPI_REAL4,
     &                 iprocsend,51,MPI_COMM_WORLD, ierr)
                  CALL MPI_SEND(gradhs(1,istart),inumber,i2REAL4,
     &                 iprocsend,52,MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
c      IF (itime.EQ.000)
                       print *,iproc,': sent divv,curlv,gradhs to ',
     &                 iprocsend,' and has sent ',inumber,' values'
#endif
c
c--Send back quantities for MHD gradients (if required)
c
                  IF (imhd.EQ.idim) THEN
c                     CALL MPI_SEND(gradmhd(1,istart),inumber,i12REAL8,
c     &                    iprocsend,53,MPI_COMM_WORLD, ierr)
                     CALL MPI_SEND(gradmhd(1,istart),inumber,i15REAL8,
     &                    iprocsend,53,MPI_COMM_WORLD, ierr)

                  ENDIF

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
         WRITE (*,*) 'ERROR2 - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit
      ENDIF

      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)
c      CALL MPI_TYPE_FREE(i12REAL8,ierr)
      CALL MPI_TYPE_FREE(i15REAL8,ierr)

#ifdef MPIDEBUG
c      IF (itime.EQ.000)
           print *,iproc,': FINISHED ',itime
c         IF (iproc.EQ.0) print *,iproc,':B',dvxyzu(1,625),
c     &        dvxyzu(2,625),dvxyzu(3,625),xyzmh(1,625),
c     &        xyzmh(2,625),xyzmh(3,625)
c         IF (iproc.EQ.1) print *,iproc,':B',dvxyzu(1,25786),
c     &        dvxyzu(2,25786),dvxyzu(3,25786),xyzmh(1,25786),
c     &        xyzmh(2,25786),xyzmh(3,25786)
c         IF (iproc.EQ.2) print *,iproc,':B',dvxyzu(1,460),
c     &        dvxyzu(2,460),dvxyzu(3,460),xyzmh(1,460),
c     &        xyzmh(2,460),xyzmh(3,460)
c         IF (iproc.EQ.3) print *,iproc,':B',dvxyzu(1,25646),
c     &        dvxyzu(2,25646),dvxyzu(3,25646),xyzmh(1,25646),
c     &        xyzmh(2,25646),xyzmh(3,25646)
c      print *,' '
#endif

           IF (itiming) THEN
              CALL getused(td22)
              td2 = td2 + (td22 - td21)
           ENDIF
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_end,listp,iphase,rho,divv,curlv,gradhs,xyzmh,cnormk)
C$OMP& shared(nptmass,iptmass,listpm,gradmhd,varmhd,Bxyz,rhomin)
C$OMP& shared(Bextx,Bexty,Bextz,radcrit2)
C$OMP& private(j,l,ipart,rhoi,dhdrhoi,iptcur)
C$OMP& private(rxxi,rxyi,rxzi,ryyi,ryzi,rzzi)
C$OMP& private(dalphaxi,dalphayi,dalphazi,dbetaxi,dbetayi,dbetazi)
C$OMP& private(denom,ddenom,gradalphaxi,gradalphayi,gradalphazi)
C$OMP& private(gradbetaxi,gradbetayi,gradbetazi,term)
C$OMP& private(dgammaxi,dgammayi,dgammazi)
C$OMP& private(gradgammaxi,gradgammayi,gradgammazi)
C$OMP& reduction(MAX:rhonext)
C$OMP& reduction(+:nwarnroundoff)
      DO j = 1, nlst_end
         ipart = listp(j)
         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN
c
c--Need to complete calculations of divv, curlv, and gradhs(2,*)
c
            rhoi = rho(ipart)
            divv(ipart) = divv(ipart)*gradhs(1,ipart)
            curlv(ipart) = curlv(ipart)*gradhs(1,ipart)

            dhdrhoi = - xyzmh(5,ipart)/(3.*(rhoi + rhomin))

            gradhs(2,ipart) = gradhs(2,ipart)*dhdrhoi
c
c--Find particle with highest density outside radcrit of point mass
c
            IF (iptmass.NE.0) THEN
               DO l = 1, nptmass
                  iptcur = listpm(l)
                  IF ( (xyzmh(1,ipart) - xyzmh(1,iptcur))**2 +
     &                 (xyzmh(2,ipart) - xyzmh(2,iptcur))**2 +
     &                 (xyzmh(3,ipart) - xyzmh(3,iptcur))**2 
     &                 .LT.radcrit2) GOTO 50
               END DO
               rhonext = MAX(rhonext, rho(ipart))
            ENDIF
 50         CONTINUE

c
c--calculate B from the evolved Euler potentials
c
            IF (imhd.EQ.idim .AND. iphase(ipart).EQ.0) THEN
               IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
                  rxxi = gradmhd(1,ipart)
                  rxyi = gradmhd(2,ipart)
                  rxzi = gradmhd(3,ipart)
                  ryyi = gradmhd(4,ipart)
                  ryzi = gradmhd(5,ipart)
                  rzzi = gradmhd(6,ipart)
                  dalphaxi = gradmhd(7,ipart)
                  dalphayi = gradmhd(8,ipart)
                  dalphazi = gradmhd(9,ipart)
                  dbetaxi = gradmhd(10,ipart)
                  dbetayi = gradmhd(11,ipart)
                  dbetazi = gradmhd(12,ipart)
c
c--compute grad alpha and grad beta using exact linear interpolation
c  (see Price 2004)
c
                  denom = rxxi*ryyi*rzzi + 2.*rxyi*rxzi*ryzi
     &               - rxxi*ryzi*ryzi - ryyi*rxzi*rxzi - rzzi*rxyi*rxyi

                  IF (abs(denom).GT.tiny) THEN
                     ddenom = 1./denom

                     gradalphaxi =(dalphaxi*(ryyi*rzzi - ryzi*ryzi)
     &                        + dalphayi*(rxzi*ryzi - rzzi*rxyi)
     &                        + dalphazi*(rxyi*ryzi - rxzi*ryyi))*ddenom
                     gradalphayi =(dalphaxi*(ryzi*rxzi - rxyi*rzzi)
     &                        + dalphayi*(rzzi*rxxi - rxzi*rxzi)
     &                        + dalphazi*(rxyi*rxzi - rxxi*ryzi))*ddenom
                     gradalphazi =(dalphaxi*(rxyi*ryzi - rxzi*ryyi)
     &                        + dalphayi*(rxyi*rxzi - rxxi*ryzi)
     &                        + dalphazi*(rxxi*ryyi - rxyi*rxyi))*ddenom
                     gradbetaxi = (dbetaxi*(ryyi*rzzi - ryzi*ryzi)
     &                        + dbetayi*(rxzi*ryzi - rzzi*rxyi)
     &                        + dbetazi*(rxyi*ryzi - rxzi*ryyi))*ddenom
                     gradbetayi = (dbetaxi*(ryzi*rxzi - rxyi*rzzi)
     &                        + dbetayi*(rzzi*rxxi - rxzi*rxzi)
     &                        + dbetazi*(rxyi*rxzi - rxxi*ryzi))*ddenom
                     gradbetazi = (dbetaxi*(rxyi*ryzi - rxzi*ryyi)
     &                        + dbetayi*(rxyi*rxzi - rxxi*ryzi)
     &                        + dbetazi*(rxxi*ryyi - rxyi*rxyi))*ddenom
                  ELSE
c
c--standard first derivative (use in case of round-off error problems)
c
                     nwarnroundoff = nwarnroundoff + 1
                     term = cnormk*gradhs(1,ipart)/rhoi
                     gradalphaxi = dalphaxi*term
                     gradalphayi = dalphayi*term
                     gradalphazi = dalphazi*term
                     gradbetaxi = dbetaxi*term
                     gradbetayi = dbetayi*term
                     gradbetazi = dbetazi*term
                  ENDIF
                  IF (varmhd.EQ.'vecp') THEN
c
c--vector potential: B = curl A
c
                     dgammaxi = gradmhd(13,ipart)
                     dgammayi = gradmhd(14,ipart)
                     dgammazi = gradmhd(15,ipart)
                     IF (abs(denom).GT.tiny) THEN
                        gradgammaxi = (dgammaxi*(ryyi*rzzi - ryzi*ryzi)
     &                        + dgammayi*(rxzi*ryzi - rzzi*rxyi)
     &                        + dgammazi*(rxyi*ryzi - rxzi*ryyi))*ddenom
                        gradgammayi = (dgammaxi*(ryzi*rxzi - rxyi*rzzi)
     &                        + dgammayi*(rzzi*rxxi - rxzi*rxzi)
     &                        + dgammazi*(rxyi*rxzi - rxxi*ryzi))*ddenom
                        gradgammazi = (dgammaxi*(rxyi*ryzi - rxzi*ryyi)
     &                        + dgammayi*(rxyi*rxzi - rxxi*ryzi)
     &                        + dgammazi*(rxxi*ryyi - rxyi*rxyi))*ddenom
                     ELSE
                        gradgammaxi = dgammaxi*term
                        gradgammayi = dgammayi*term
                        gradgammazi = dgammazi*term
                     ENDIF
                     term = gradgammayi - gradbetazi
                     Bxyz(1,ipart) = term + Bextx
                     term = gradalphazi - gradgammaxi
                     Bxyz(2,ipart) = term + Bexty
                     term = gradbetaxi - gradalphayi
                     Bxyz(3,ipart) = term + Bextz
                  ELSE
c
c--Euler potentials: B = grad alpha cross grad beta
c
                     term= gradalphayi*gradbetazi-gradalphazi*gradbetayi
                     Bxyz(1,ipart)= term
                     term= gradalphazi*gradbetaxi-gradalphaxi*gradbetazi
                     Bxyz(2,ipart)= term
                     term= gradalphaxi*gradbetayi-gradalphayi*gradbetaxi
                     Bxyz(3,ipart)= term
                  ENDIF
c             ELSE
c
c--add self contribution and store smoothed velocity
c
c               vsmooth(1,ipart) = cnormk*(vbarxi
c     &                                  + weight*vxyzu(1,ipart)*wij(0))
c               vsmooth(2,ipart) = cnormk*(vbaryi
c     &                                  + weight*vxyzu(2,ipart)*wij(0))
c               vsmooth(3,ipart) = cnormk*(vbarzi
c     &                                  + weight*vxyzu(3,ipart)*wij(0))
c               Bsmooth(1,ipart) = cnormk*(Bbarxi
c     &                                  + weight*Bevol(1,ipart)*wij(0))
c               Bsmooth(2,ipart) = cnormk*(Bbaryi
c     &                                  + weight*Bevol(2,ipart)*wij(0))
c               Bsmooth(3,ipart) = cnormk*(Bbarzi
c     &                                  + weight*Bevol(3,ipart)*wij(0))
c               print*,ipart,'v       = ',vxyzu(1:3,ipart)
c               print*,ipart,'vsmooth = ',vsmooth(:,ipart)
               ENDIF
            ENDIF
         ENDIF
      END DO

C$OMP PARALLEL default(none)
C$OMP& shared(npart,ntot,ireal,rho,dumrho,xyzmh,pr,vsound,divv,curlv)
C$OMP& shared(vxyzu,gradhs,Bxyz,encal,ibound,ekcle,varmhd,listp,Bevol)
C$OMP& shared(nlst_in,nlst_end,listparents,nlistupdated,iupdated,icall)
C$OMP& shared(dt,itime,iphase)
C$OMP& private(i,j,k,ipart)
c
c--Update density and quantities depending on density for inactive particles
c     that are neighbours of active particles
c
C$OMP DO SCHEDULE (runtime)
      DO i = 1, nlistupdated
         j = listparents(i)
         iupdated(j) = .FALSE.
         CALL extrapolate(j,dt,itime,dumrho,vxyzu,pr,vsound,
     &        ekcle,Bxyz,Bevol)
      END DO
C$OMP END DO
c
c--Copy changed values onto ghost particles
c
C$OMP DO SCHEDULE (runtime)
      DO i = npart + 1, ntot
         j = ireal(i)
         rho(i) = rho(j)
         dumrho(i) = dumrho(j)
         xyzmh(5,i) = xyzmh(5,j)
         pr(i) = pr(j)
         vsound(i) = vsound(j)
         divv(i) = divv(j)
         curlv(i) = curlv(j)
         vxyzu(4,i) = vxyzu(4,j)
         gradhs(1,i) = gradhs(1,j)
         gradhs(2,i) = gradhs(2,j)
         IF (imhd.EQ.idim) THEN
            IF (varmhd.EQ.'eulr' .or. varmhd.EQ.'vecp') THEN
               Bxyz(1,i) = Bxyz(1,j)
               Bxyz(2,i) = Bxyz(2,j)
               Bxyz(3,i) = Bxyz(3,j)            
            ENDIF
         ENDIF
         IF (encal.EQ.'r' .AND. ibound/10.EQ.10) THEN
            IF(icall.EQ.1) THEN
               DO k=1,5
                  ekcle(k,i) = ekcle(k,j)
               END DO
            ELSE
               ekcle(1,i) = ekcle(1,j)
            ENDIF
         ENDIF
      END DO
C$OMP END DO
c
c--Calculate B from the evolved magnetic field variable
c
      IF (imhd.EQ.idim) THEN
         IF (varmhd.eq.'Bvol') THEN
C$OMP DO SCHEDULE (runtime)
            DO i=nlst_in,nlst_end
               ipart = listp(i)
               Bxyz(1,ipart) = Bevol(1,ipart)
               Bxyz(2,ipart) = Bevol(2,ipart)
               Bxyz(3,ipart) = Bevol(3,ipart)
            ENDDO
C$OMP END DO
C$OMP DO SCHEDULE (runtime)
            DO ipart=npart+1,ntot
               Bxyz(1,ipart) = Bevol(1,ipart)
               Bxyz(2,ipart) = Bevol(2,ipart)
               Bxyz(3,ipart) = Bevol(3,ipart)
            ENDDO
C$OMP END DO
         ELSEIF (varmhd.EQ.'Brho') THEN
C$OMP DO SCHEDULE (runtime)
            DO i=nlst_in,nlst_end
               ipart = listp(i)
               Bxyz(1,ipart) = Bevol(1,ipart)*dumrho(ipart)
               Bxyz(2,ipart) = Bevol(2,ipart)*dumrho(ipart)
               Bxyz(3,ipart) = Bevol(3,ipart)*dumrho(ipart)
            ENDDO
C$OMP END DO
C$OMP DO SCHEDULE (runtime)
            DO ipart=npart+1,ntot
               Bxyz(1,ipart) = Bevol(1,ipart)*dumrho(ipart)
               Bxyz(2,ipart) = Bevol(2,ipart)*dumrho(ipart)
               Bxyz(3,ipart) = Bevol(3,ipart)*dumrho(ipart)
            ENDDO      
C$OMP END DO
         ELSEIF (varmhd.NE.'eulr' .AND. varmhd.NE.'vecp') THEN
            STOP 'unknown MHD variable in Bevol->Bxyz conversion'
         ENDIF
      ENDIF

C$OMP END PARALLEL
ccC$OMP END PARALLEL DO

      IF (nwarnup.GT.0) THEN
         WRITE (iprint,*) 'WARNING: restricted h jump (up) ',
     &        nwarnup,' times'
      ENDIF
      IF (nwarndown.GT.0) THEN
         WRITE (iprint,*) 'WARNING: restricted h jump (down) ',
     &        nwarndown,' times'
      ENDIF
      IF (nwarnroundoff.GT.0) THEN
         IF (varmhd.EQ.'vecp') THEN
           WRITE (iprint,*) 'WARNING: denom in vecp gradients zero on ',
     &        nwarnroundoff,' particles'
         ELSEIF (varmhd.EQ.'eulr') THEN
          WRITE (iprint,*) 'WARNING: denom in euler gradients zero on ',
     &        nwarnroundoff,' particles'
         ENDIF
      ENDIF
      IF (nbisection.GT.0) THEN
         WRITE(iprint,*) 'WARNING: used bisection on ',nbisection,
     &       ' particles'
      ENDIF
c
c--Possible to create a point mass
c
      IF (rhonext.GT.rhocrea .AND. iptmass.NE.0 .AND. icall.EQ.3) THEN
c
c--Find particle with highest density outside radcrit of point mass
c
         irhonex = 0
         DO n = nlst_in, nlst_end
            ipart = listp(n)
            IF (rho(ipart).EQ.rhonext) THEN
               irhonex = ipart
            ENDIF
         END DO
         IF (irhonex.EQ.0) THEN
            WRITE(iprint,*)'Failed to find densest particle ',rhonext
            CALL quit
         ENDIF
c
c--Make sure that all neighbours of point mass candidate are being
c     done on this time step. Otherwise, not possible to accrete
c     them to form a point mass and it may create a point mass without
c     accreting many particles!
c
c--NOTE: For MPI job, still need to check for neighbours on another process
c	For this, need to check within forcei when particles are sent back
c	Only need to check particle that has passed all local tests and
c	also has neighbours on other process.
c	Need to check that all its neighbours are currently being evolved.
c	Also need to check that it is the candidate sink particle with the
c	greatest density on all MPI processes (only form 1 sink at a time).
c	These things only need to be checked if there is a particle for
c	which icreate=1 on any MPI process.
c
         IF ((2.0*xyzmh(5,irhonex)).LT.hacc) THEN
            WRITE(iprint,*)'Ptmass creation passed h ',xyzmh(5,irhonex)

	    nlist = 0
            CALL getneigh(irhonex,npart,xyzmh(5,irhonex),xyzmh,idim,
     &           nlist,iptneigh,nearl)

            iokay = 1
            DO n = 1, nlist
               j = nearl(n)
               IF (it0(j).NE.itime) iokay = 0
            END DO
c
c--Set creation flag to true. Other tests done in accrete.f
c
            IF (iokay.EQ.1) THEN
               icreate = 1 
               WRITE(iprint,*) ' and all local particles on step'
            ELSE
               WRITE(iprint,*) ' but not all local particles on step'
            ENDIF
         ELSE
            WRITE(iprint,*) 'Ptmass creation failed on h ',
     &           xyzmh(5,irhonex)
         ENDIF
      ENDIF

#ifdef MPI
  250 CALL MPI_ALLREDUCE(icreate,icreatetot,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
#ifdef MPIDEBUGS
      print *,iproc,': icreatetot ',icreatetot,icreate
#endif
      IF (icreatetot.GE.1) THEN
c
c--Check most dense particle first
c
	 IF (icreatetot.GT.1) THEN
   	    IF (icreate.EQ.0) rhonext = 0.
	    CALL MPI_ALLREDUCE(rhonext,rhonextmax,1,MPI_REAL4,MPI_MAX,
     &           MPI_COMM_WORLD,ierr)
#ifdef MPIDEBUGS
      print *,iproc,': rhonextmax ',rhonextmax,rhonext
#endif
            IF (rhonext.EQ.rhonextmax) THEN
               itest = iproc
            ELSE
               itest = 0
            ENDIF
c
c--If more than one process has a particle with density=rhonextmax, then
c	take the process with the largest rank first
c
            CALL MPI_ALLREDUCE(itest,itestmax,1,MPI_INTEGER,MPI_MAX,
     &           MPI_COMM_WORLD,ierr)
            iproccreate = itestmax
#ifdef MPIDEBUGS
      print *,iproc,': iproccreate ',iproccreate
#endif
	 ELSE
            itest = icreate*iproc
            CALL MPI_ALLREDUCE(itest,itestmax,1,MPI_INTEGER,MPI_MAX,
     &           MPI_COMM_WORLD,ierr)
            iproccreate = itestmax
#ifdef MPIDEBUGS
      print *,iproc,': iproccreate ',iproccreate
#endif
         ENDIF
c
c--iproccreate is the process whose particle is being tested.  All MPI
c	processes know this.
c
	 IF (iproccreate.EQ.iproc) THEN
c
c--Check for neighbours on other processes - only need to loop over "sendback"
c	lists to see which processes need to check
c
	    IF (nneighsentanyatall) THEN
	       DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (nneighsentany(j+1)) THEN
                        CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8,
     &	                     ierr)
                        CALL MPI_TYPE_COMMIT(i5REAL8,ierr)

                        CALL MPI_SEND(xyzmh(1,irhonex),1,i5REAL8,j,60,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUGS
      print *,iproc,': sent irhonex data to ',j
#endif

	                CALL MPI_TYPE_FREE(i5REAL8,ierr)
c
c--Receive back information on whether neighbouring particles are all on 
c	current timestep or not
c
                        CALL MPI_RECV(iokay,1,MPI_INTEGER,j,61,
     &                       MPI_COMM_WORLD, istatus, ierr)
#ifdef MPIDEBUGS
      print *,iproc,': was told that iokay=',iokay,' by ',j
#endif
                        IF (iokay.EQ.0) icreate = 0
                     ENDIF
                  ENDIF
               END DO
            ENDIF
	 ELSE
            IF (nneightogetback(iproccreate+1).GT.0) THEN
c
c--Expect to get particle which has neighbours on this process
c
               itest = 2*ntot+inumbertotal+3
               IF (itest.GT.mmax2) THEN
                  WRITE (*,*) 'ERROR - itest.GT.mmax2 ',itest,mmax2
                  CALL quit
               ENDIF
               CALL MPI_RECV(xyzmh(1,itest), 5, MPI_REAL8,
     &              MPI_ANY_SOURCE, 60, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (inumber.NE.5) THEN
                  WRITE (iprint,*) 'ERROR - Particle test ',inumber
                  CALL quit
               ENDIF
c
c--Find neighbours of this particle and test to see whether they are all
c	being evolved on the current step
c
               nlist = 0
               CALL getneigh(itest,npart,xyzmh(5,itest),xyzmh,
     &              mmax2,nlist,iptneigh,nearl)

               iokay = 1
               DO n = 1, nlist
                  j = nearl(n)
                  IF (it0(j).NE.itime) iokay = 0
               END DO
c
c--Send back the result
c
	       CALL MPI_SEND(iokay,1,MPI_INTEGER,iprocrec,61,
     &              MPI_COMM_WORLD, ierr)
	    ENDIF
	 ENDIF
c
c--Everyone needs to find out the overall result
c
         iokay = icreate
         CALL MPI_BCAST(iokay,1,MPI_INTEGER,iproccreate,
     &        MPI_COMM_WORLD, ierr)         
#ifdef MPIDEBUGS
      print *,iproc,': has been told that creation is ',iokay,
     &        iproccreate
#endif
         IF (iokay.EQ.0) GOTO 250
         icreatetot = 1
         IF (iproc.NE.iproccreate) icreate = 0
#ifdef MPIDEBUGS
      print *,iproc,': icreate, icreatetot ',icreate,
     &        icreatetot
#endif
      ENDIF
#else
      icreatetot = icreate
#endif
      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine densityiterate')

#ifdef MPI
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      IF (itime.EQ.000)
c           print *,iproc,': Exited densityiterate'
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      CALL quit
#endif

c      print *,'Exited densityiterate'

      RETURN
      END



c=============================================================================
c
c--Subroutine for calculating rho, gradh and number of neighbours
c
      SUBROUTINE densitygradh(ipart,npart,ntot,nlst_doing,
     &     xyzmh,numneighi,dvxyzu)

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax2), dvxyzu(4,idim2)

      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/initpt'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
c
c--Needed for MPI code
c
#ifdef MPI
#ifdef MPIDEBUG
      IF (iproblem.EQ.2097152) print *,iproc,' Entered ',
     &     ' densitygrad ',ipart,npart,ntot,iphase(ipart),xyzmh(1,ipart)
#endif
#endif

      IF (ipart.GT.npart) THEN
         iparttree = ipart + ntot + 2
      ELSE
         iparttree = ipart
         nneighprior = nneigh(ipart)
      ENDIF

      IF (nlst_doing.GT.nptmasstot .OR. nptmasstot.NE.nptmasstotlast
     &     .OR. iptintree.EQ.2 .OR. initialptm.EQ.5) THEN
         CALL treef(ipart,npart,ntot,xyzmh,acc,igphi,fsx,fsy,fsz,epot)

         IF (iptintree.NE.2) THEN
            IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
               iptcur = listrealpm(ipart)
               IF (listpm(iptcur).NE.ipart) THEN
                  WRITE (*,*) 'ERROR - listpm(iptcur).NE.ipart'
                  CALL quit
               ENDIF
               gravxyzpstore(1,iptcur) = fsx
               gravxyzpstore(2,iptcur) = fsy
               gravxyzpstore(3,iptcur) = fsz
               gravxyzpstore(4,iptcur) = epot
            ENDIF
         ENDIF
      ELSE
         iptcur = listrealpm(ipart)
         IF (listpm(iptcur).NE.ipart) THEN
            WRITE (*,*) 'ERROR - listpm(iptcur).NE.ipart'
            CALL quit
         ENDIF
         fsx = gravxyzpstore(1,iptcur)
         fsy = gravxyzpstore(2,iptcur)
         fsz = gravxyzpstore(3,iptcur)
         epot = gravxyzpstore(4,iptcur)
      ENDIF

      rhoi = 0.
      gradhi = 0.
      numneighi = 0
      IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN
         numneighi = nneigh(ipart)
c
c--Calculate density contribution by looping over interacting neighbors
c
         xi = xyzmh(1,iparttree)
         yi = xyzmh(2,iparttree)
         zi = xyzmh(3,iparttree)
         pmassi = xyzmh(4,iparttree)
         hi = xyzmh(5,iparttree)
         hi1 = 1./hi
         hi21 = hi1*hi1
         hi31 = hi21*hi1
         hi41 = hi21*hi21

         IF (numneighi.GT.0) THEN
            DO k = 1, numneighi
               j = neighlist(k)

               IF (iphase(j).EQ.iphase(ipart)) THEN

                  dx = xi - xyzmh(1,j)
                  dy = yi - xyzmh(2,j)
                  dz = zi - xyzmh(3,j)
#ifdef PERIODIC_NO_GHOSTS
                  CALL modbound(dx,dy,dz)
#endif
                  pmassj = xyzmh(4,j)

                  rij2 = dx*dx + dy*dy + dz*dz + tiny
                  v2 = rij2*hi21

                  IF (v2.LT.radkernel**2) THEN
                     rij1 = SQRT(rij2)
c
c--Get kernel quantities from interpolation in table
c
                     index = v2*ddvtable
                     dxx = v2 - index*dvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dwdx = (wij(index1) - wij(index))*ddvtable
                     wtij = (wij(index) + dwdx*dxx)*hi31
                     dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                     grwtij = (grwij(index) + dgrwdx*dxx)*hi41/rij1
c
c--Derivative w.r.t. h for grad h correction terms (and dhdrho)
c
                     dwdhi = (-rij2*grwtij - 3.*wtij)*hi1
                     gradhi = gradhi + pmassj*dwdhi
c
c--Compute density
c
                     rhoi = rhoi + pmassj*wtij
                  ENDIF
               ENDIF
            END DO
         ENDIF
      ENDIF
c
c--Needed for MPI code
c
      IF (ipart.GT.npart) THEN
         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN
            rho(ipart) = cnormk*rhoi
            dumrho(ipart) = cnormk*gradhi
            nneigh(ipart) = numneighi
         ENDIF
         dvxyzu(1,ipart) = fsx
         dvxyzu(2,ipart) = fsy
         dvxyzu(3,ipart) = fsz
         dvxyzu(4,ipart) = epot
      ELSE
         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN
            rho(ipart) = rho(ipart) +
     &           cnormk*(rhoi + selfnormkernel*pmassi*hi31)
            dumrho(ipart) = dumrho(ipart) + cnormk*(gradhi +
     &           selfnormkernel*pmassi*(-3.*hi41))
            nneigh(ipart) = nneighprior + numneighi
         ENDIF
         dvxyzu(1,ipart) = dvxyzu(1,ipart) + fsx
         dvxyzu(2,ipart) = dvxyzu(2,ipart) + fsy
         dvxyzu(3,ipart) = dvxyzu(3,ipart) + fsz
         dvxyzu(4,ipart) = dvxyzu(4,ipart) + epot
      ENDIF
c
c--Add gravity from sink particles
c
c      IF (nptmass.GE.1) print *,iproc,': ',ipart,dvxyzu(1,ipart)

      IF (nptmass.GT.0) THEN
c
c--Add gravity on gas and sinks from sink particles
c
	 IF (iptintree.EQ.0) THEN
            CALL gpti(ipart, npart, ntot, xyzmh, dvxyzu)
c
c--On sinks and gas but without gravity from gas on sinks
c
	 ELSEIF (iptintree.EQ.1) THEN
            CALL gforspt(ipart,npart,ntot,xyzmh,dvxyzu)
         ENDIF
      ENDIF
c
c--Add gravity on sink particles from gas
c
      IF (iptintree.EQ.0 .AND. 
     &     iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
         CALL gptfromgas(ipart, npart, ntot, xyzmh, dvxyzu)
      ENDIF

      RETURN
      END

c=============================================================================
c
c--Subroutine for calculating other quantities such as divv, curlv, alpha,
c     and beta for Euler potentials, etc
c
c--Also now updates the densities
c
c
      SUBROUTINE divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
     &     numneighi,neighbourslist,gradmhd)

      INCLUDE 'idim'

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2), Bevol(3,imhd3),
     &     gradmhd(15,imhd2)
      DIMENSION neighbourslist(numneighi)

      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/phase'
c
c--Needed for MPI code
c
      IF (ipart.GT.npart) THEN
         iparttree = ipart + ntot + 2
      ELSE
         iparttree = ipart
      ENDIF

      xi = xyzmh(1,iparttree)
      yi = xyzmh(2,iparttree)
      zi = xyzmh(3,iparttree)
      pmassi = xyzmh(4,iparttree)
      hi = xyzmh(5,iparttree)

      hi1 = 1./hi
      hi21 = hi1*hi1
      hi31 = hi21*hi1
      hi41 = hi21*hi21

      divvi = 0.
      curlvxi = 0.
      curlvyi = 0.
      curlvzi = 0.
      gradsofti = 0.

      IF (imhd.EQ.idim) THEN
         DO k = 1, 15
            gradmhd(k,ipart) = 0.
         END DO
      ENDIF

c      vbarxi = 0.
c      vbaryi = 0.
c      vbarzi = 0.
c      Bbarxi = 0.
c      Bbaryi = 0.
c      Bbarzi = 0.

      vxi = vxyzu(1,ipart)
      vyi = vxyzu(2,ipart)
      vzi = vxyzu(3,ipart)
c
c--Calculate B from the evolved magnetic field variable
c
      IF (imhd.EQ.idim .AND. iphase(ipart).EQ.0) THEN
         IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
            alphai= Bevol(1,ipart)
            betai= Bevol(2,ipart)
            gammai = Bevol(3,ipart)
         ENDIF
      ENDIF

      DO k = 1, numneighi
         j = neighbourslist(k)

         IF (iphase(j).EQ.iphase(ipart)) THEN

         dx = xi - xyzmh(1,j)
         dy = yi - xyzmh(2,j)
         dz = zi - xyzmh(3,j)
#ifdef PERIODIC_NO_GHOSTS
         IF (varmhd.EQ.'eulr') THEN
            dalpha= alphai - Bevol(1,j)
            dbeta= betai - Bevol(2,j)
            CALL modboundeulr(dx,dy,dz,dalpha,dbeta)
         ELSEIF (varmhd.EQ.'vecp') THEN
            write(iprint,*) 
     &           'densityiterate: modbound not implemented for vecp'
            CALL quit
         ELSE
            CALL modbound(dx,dy,dz)
         ENDIF
#endif
         pmassj = xyzmh(4,j)
         rij2 = dx*dx + dy*dy + dz*dz + tiny
         v2 = rij2*hi21

         IF (v2.LT.radkernel**2) THEN
            rij1 = SQRT(rij2)

            dvx = vxi - vxyzu(1,j)
            dvy = vyi - vxyzu(2,j)
            dvz = vzi - vxyzu(3,j)
c
c--Get kernel quantities from interpolation in table
c
            index = v2*ddvtable
            dxx = v2 - index*dvtable
            index1 = index + 1
            IF (index.GE.itable) THEN
               index = itable
               index1 = itable
            ENDIF
            dwdx = (wij(index1) - wij(index))*ddvtable
            wkern = (wij(index) + dwdx*dxx)
            wtij = wkern*hi31
            dgrwdx = (grwij(index1) - grwij(index))*ddvtable
            grwtij = (grwij(index) + dgrwdx*dxx)*hi41/rij1
            dpotdh = (dphidh(index1) - dphidh(index))*ddvtable
            dphi = (dphidh(index) + dpotdh*dxx)*hi21
c
c--Derivative of gravitational potential w.r.t. h
c
            gradsofti = gradsofti - pmassj*dphi
c
c--Velocity divergence times density
c
            projv = grwtij*(dvx*dx + dvy*dy + dvz*dz)
            divvi = divvi - pmassj*projv

            IF (iphase(ipart).EQ.0) THEN
c
c--Velocity curl in 3D times density
c
               procurlvz = grwtij*(dvy*dx - dvx*dy)
               procurlvy = grwtij*(dvx*dz - dvz*dx)
               procurlvx = grwtij*(dvz*dy - dvy*dz)

               curlvxi = curlvxi - pmassj*procurlvx
               curlvyi = curlvyi - pmassj*procurlvy
               curlvzi = curlvzi - pmassj*procurlvz
c
c--get B from the evolved Euler potentials
c
c
c--grad alpha and grad beta (Euler potentials)
c
               IF (imhd.EQ.idim) THEN
                  IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
                     grpmi= pmassj*grwtij
c
c--rxx, rxy, rxz, ryy, ryz, rzz
c
                     gradmhd(1,ipart) = gradmhd(1,ipart) - grpmi*dx*dx
                     gradmhd(2,ipart) = gradmhd(2,ipart) - grpmi*dx*dy
                     gradmhd(3,ipart) = gradmhd(3,ipart) - grpmi*dx*dz
                     gradmhd(4,ipart) = gradmhd(4,ipart) - grpmi*dy*dy
                     gradmhd(5,ipart) = gradmhd(5,ipart) - grpmi*dy*dz
                     gradmhd(6,ipart) = gradmhd(6,ipart) - grpmi*dz*dz
#ifndef PERIODIC_NO_GHOSTS
                     dalpha= alphai - Bevol(1,j)
                     dbeta= betai - Bevol(2,j)
                     dgamma = gammai - Bevol(3,j)
#endif
c
c--dalpha(x,y,z) and dbeta(x,y,z)
c
                  gradmhd(7,ipart) = gradmhd(7,ipart) -grpmi*dalpha*dx
                  gradmhd(8,ipart) = gradmhd(8,ipart) -grpmi*dalpha*dy
                  gradmhd(9,ipart) = gradmhd(9,ipart) -grpmi*dalpha*dz

                  gradmhd(10,ipart) = gradmhd(10,ipart)-grpmi*dbeta*dx
                  gradmhd(11,ipart) = gradmhd(11,ipart)-grpmi*dbeta*dy
                  gradmhd(12,ipart) = gradmhd(12,ipart)-grpmi*dbeta*dz

                  !--strictly only for vecp but should not matter
                  gradmhd(13,ipart) = gradmhd(13,ipart)-grpmi*dgamma*dx
                  gradmhd(14,ipart) = gradmhd(14,ipart)-grpmi*dgamma*dy
                  gradmhd(15,ipart) = gradmhd(15,ipart)-grpmi*dgamma*dz

c               ELSE
c
c--smoothed velocity for use in the B or B/rho evolution
c
c                  vbarxi = vbarxi + weight*vxyzu(1,j)*wkern
c                  vbaryi = vbaryi + weight*vxyzu(2,j)*wkern
c                  vbarzi = vbarzi + weight*vxyzu(3,j)*wkern
c
c--smoothed Bevol for div B reduction
c
c                  Bbarxi = Bbarxi + weight*Bevol(1,j)*wkern
c                  Bbaryi = Bbaryi + weight*Bevol(2,j)*wkern
c                  Bbarzi = Bbarzi + weight*Bevol(3,j)*wkern
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         ENDIF

      END DO
c
c--Need for MPI
c
      IF (ipart.GT.npart) THEN
         divv(ipart) = cnormk*divvi
         curlv(ipart) = cnormk*SQRT(curlvxi**2+curlvyi**2+curlvzi**2)
         gradhs(2,ipart) = gradsofti
      ELSE
         divv(ipart) = cnormk*divvi
         curlv(ipart) = cnormk*SQRT(curlvxi**2+curlvyi**2+curlvzi**2)
c
c--Add self contribution
c
         gradhs(2,ipart) = gradsofti - pmassi*dphidh(0)*hi21
      ENDIF

      RETURN
      END
