      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &     dvxyzu,dha,npart,ntot,ireal,alphaMM,ekcle,Bxyz,dBxyz)
c************************************************************
c                                                           *
c  This subroutine drives the computation of the forces on  *
c     every particle on the list.                           *
c                                                           *
c************************************************************
#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'idim'
      INCLUDE 'igrape'

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
      DIMENSION xyzmhpassed(5,idim), fxyzpext(4,idim),
     &     fxyzprec(4,idim), llistsend(idim)
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
            CALL gforspt(xyzmh,dvxyzu)
         ENDIF
      ENDIF
c
c--Only want to do density and updates for non-sinks
c
      IF (nlst_end.GT.nptmass) THEN
c
c--Compute the pressure, divv, etc. on list particles but
c     do not search twice the list particles
c
         IF (itiming) CALL getused(tdens1)

         CALL densityi(npart,xyzmh,vxyzu,ekcle,
     &        nlst_in,nlst_end,llist,itime)

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
         nlstall = nlst + nlistavail
         IF (itiming) THEN
            CALL getused(td12)
            td1 = td1 + (td12 - td11)
         ENDIF

#ifdef MPI
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      DO i = 1, nlst
         llistsend(i) = llist(i) - 1
      END DO
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlst, 1, llistsend, i5REAL8,
     &     indexMPI5, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlst, 1, llistsend, i4REAL8,
     &     indexMPI4, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlst, 1, llistsend, MPI_REAL4,
     &     indexMPI1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlst, 1, llistsend,MPI_INTEGER,
     &     indexMPI_I1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_I1,ierr)

      DO i = 1, numproc - 1
         iahead = MOD(iproc+i,numproc)
         ibehind = MOD(numproc+iproc-i,numproc)
c
c--Send to node ahead, receive from node behind
c
c         print *,iproc,': Sending xyzmh to ',iahead,' rec from ',ibehind
c         DO j = 1, 10
c            print *,iproc,': S   ',j,xyzmh(1,llist(j)),xyzmh(3,llist(j))
c         END DO

         CALL MPI_SENDRECV(xyzmh,1,indexMPI5,iahead,0,
     &        xyzmhpassed, idim, i5REAL8, ibehind,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)

         CALL MPI_SENDRECV(vxyzu,1,indexMPI4,iahead,0,
     &        vxyzupassed, idim, i4REAL8, ibehind,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)

         CALL MPI_SENDRECV(rho,1,indexMPI1,iahead,0,
     &        rhopassed, idim, MPI_REAL4, ibehind,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)

         CALL MPI_SENDRECV(divv,1,indexMPI1,iahead,0,
     &        divvpassed, idim, MPI_REAL4, ibehind,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)

         CALL MPI_SENDRECV(nneigh,1,indexMPI_I1,iahead,0,
     &        nneighpassed, idim, MPI_INTEGER, ibehind,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)

c         print *,iproc,': Done Sendrecv'
c         DO j = 1, 10
c            print *,iproc,': R   ',j,xyzmhpassed(1,j),xyzmhpassed(3,j)
c         END DO
c         print *,iproc,': Received ',inumber
c
c--Loop over particles from other node to calculate forces on them
c     and to get neighbours of particles.
c     Then needs to calculate contributions to rho, divv, curlv,
c     nneigh.  From this and partial rho passed from other node, can
c     determine full rho, pressure, divv, curlv and then determines
c     contributions to gas forces and gravity from neighbours.
c     The returned dataset then consists of contributions to
c     fxyzp, rho, divv, curlv, nneigh.
c     NOTE: fxyzp also needs to include forces from sink particles.
c     Question: What happens with accretion!!
c
         DO j = 1, inumber
            DO k = 1, 4
               fxyzpext(k,j) = 0.0
            END DO
c
c--This should return neighbours lists for the particles as well as node
c     and leaf forces in fxyzpext, and the contribution to nneigh(ext)
c
c            CALL treefext(j,xyzmhpassed,xyzmh,acc,igphi,fxyzpext)
c
c--The neighbours lists should then be used to calculate contributions to
c     rho and divv (and curlv if required) and add them onto (rho,divv)passed
c
c            CALL densityext()
c
c--Then calculate additional forces
c
c            CALL forceiext(fxyzpext)
         END DO
c
c--Send quantities back to node behind, receive forces from node ahead
c
c         print *,'Sending forces to ',ibehind,' rec from ',iahead
c         DO j = 1, 10
c            print *,iproc,': S   ',j,fxyzpext(1,j)
c         END DO

         CALL MPI_SENDRECV(fxyzpext,4*inumber,MPI_REAL8,ibehind,1,
     &        fxyzprec,4*nlst,MPI_REAL8,iahead,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
         IF (ireturned.NE.4*nlst) THEN
            WRITE (*,*) 'ERROR - ireturned.NE.4*nlst'
            CALL quit
         ENDIF
c
c--CANNOT read totals directly back into rho() because multiple nodes may be
c     contributing to the same particle's values (need to add contributions
c     together).
c
         CALL MPI_SENDRECV(rhoext,inumber,MPI_REAL4,ibehind,1,
     &        rhorec,nlst,MPI_REAL4,iahead,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
         IF (ireturned.NE.4*nlst) THEN
            WRITE (*,*) 'ERROR - ireturned.NE.4*nlst'
            CALL quit
         ENDIF

         CALL MPI_SENDRECV(divvext,inumber,MPI_REAL4,ibehind,1,
     &        divvrec,nlst,MPI_REAL4,iahead,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
         IF (ireturned.NE.4*nlst) THEN
            WRITE (*,*) 'ERROR - ireturned.NE.4*nlst'
            CALL quit
         ENDIF

         CALL MPI_SENDRECV(rhoext,inumber,MPI_REAL4,ibehind,1,
     &        rhorec,nlst,MPI_REAL4,iahead,
     &        MPI_ANY_TAG, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
         IF (ireturned.NE.4*nlst) THEN
            WRITE (*,*) 'ERROR - ireturned.NE.4*nlst'
            CALL quit
         ENDIF

c         DO j = 1, 10
c            print *,iproc,': R   ',j,fxyzprec(1,j)
c         END DO


         DO j = 1, nlst
            DO k = 1, 3
               fxyzu(k,llist(j)) = fxyzu(k,llist(j)) + fxyzprec(k,j)
               
            END DO
            poten(llist(j)) = poten(llist(j)) + fxyzprec(4,j)
         END DO
      END DO

c      print *,iproc,': Received all forces'
c
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      CALL MPI_TYPE_FREE(indexMPI1,ierr)
      CALL MPI_TYPE_FREE(indexMPI4,ierr)
      CALL MPI_TYPE_FREE(indexMPI5,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)

      print *,iproc,': FINISHED'
      print *,' '
#endif


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
      ENDIF
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
c
c--Compute forces on EACH particle
c
      IF (itiming) CALL getused(tforce1)

      CALL forcei(nlst_in,nlst_end,llist,dt,itime,npart,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,dedxyz)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

      RETURN
      END
