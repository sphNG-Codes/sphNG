      SUBROUTINE densityiterate_gradh (dt,npart,ntot,xyzmh,vxyzu,dvxyzu,
     &            nlst_in,nlst_end,list,itime,ekcle,Bevol,Bxyz)
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

      DIMENSION xyzmh(5,mmax), vxyzu(4,idim), dvxyzu(4,idim)
      DIMENSION list(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevol(3,imhd),Bxyz(3,imhd)

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
c      INCLUDE 'COMMONS/vsmooth'

      LOGICAL*1 Bisection
      INTEGER isucceed(idim),hminbisec(idim),hmaxbisec(idim)

      DIMENSION ldolist(idim),hi_old(idim)
      DIMENSION gradmhd(12,imhd)
#ifdef MPI
      DIMENSION nneighrec(ineighproc),nneighsentbackold(nummaxproc)
      REAL*4 rhorec(ineighproc), dumrhorec(ineighproc)
      LOGICAL*1 isend(idim)
#endif

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine densityiterate')
c
c--Initialise
c
      uradconst = radconst/uergcc
      third = 1./3.
      rhonext = 0.
      icreate = 0
      radcrit2 = radcrit*radcrit
      numparticlesdone = numparticlesdone + nlst_end
      nwarnup = 0
      nwarndown = 0
      nwarnroundoff = 0
      stressmax = 0.
      nbisection = 0
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
c--Get rid of sink particles (they do not need to be iterated).  This list
c     is reused for the subsequent iterations when it is expected that it will 
c     contain fewer and fewer particles.
c
      nlstdo = 0
      DO n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).EQ.0) THEN
            nlstdo = nlstdo + 1
            ldolist(nlstdo) = list(n)
         ENDIF
      END DO

      nlstsend = nlstdo

C$OMP PARALLEL DO SCHEDULE(runtime) default(none) 
C$OMP& shared(nlstdo,ldolist,isucceed,xyzmh,hi_old,it1,isteps,imaxstep)
C$OMP& shared(divv,rho,dumrho,nneigh,dvxyzu,rhomin,dt,imax,poten)
#ifdef MPI
C$OMP& shared(lsendlist,isend)
#endif
C$OMP& private(n,ipart,pmassi,hi,dhdrhoi,deltat,k)
      DO n = 1, nlstdo
         ipart = ldolist(n)
         isucceed(ipart) = 0

#ifdef MPI
         lsendlist(n) = ipart - 1
         isend(n) = .TRUE.
#endif
         pmassi = xyzmh(4,ipart)
         hi = xyzmh(5,ipart)
         hi_old(ipart) = hi
c
c--Predict h
c
         dhdrhoi = - hi/(3.*(pmassi*(hfact/hi)**3 + rhomin))
         IF (it1(ipart).EQ.imax) THEN
            deltat = (dt*isteps(ipart)/2)/imaxstep
         ELSE
            deltat = (dt*isteps(ipart))/imaxstep
         ENDIF
         xyzmh(5,ipart) = xyzmh(5,ipart) - dhdrhoi*divv(ipart)*deltat

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
c
c--Iterate density calculation for particle ipart
c
#ifdef MPI
      inumofreturns = 0
      inumofsends = 0
      maxnneighsentback = 0
#endif

      DO iteration = 1, maxiterations

#ifdef MPIDEBUG
         print *,'DEN_ITERATE: iteration ',iteration
#endif

#ifdef MPI
c
c--MPI: Need to get rho, gradh, numneigh contributions from other processes
c
#ifdef MPIDEBUG
         print *,iproc,': Start MPI densityiterate ',numproc
#endif
         iprocdebug = iproc

         CALL MPI_TYPE_CONTIGUOUS(12, MPI_REAL8, i12REAL8, ierr)
         CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
         CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
         CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
         CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)

         CALL MPI_TYPE_COMMIT(i12REAL8,ierr)
         CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
         CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
         CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
         CALL MPI_TYPE_COMMIT(i2REAL4,ierr)

         DO i = 1, idim
            lblocklengths(i) = 1
         END DO
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlstsend, 1, lsendlist, 
c     &       i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_INDEXED(nlstsend, lblocklengths, lsendlist, 
     &        i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI5,ierr)
c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
         DO ii = 1, numproc - 1
            iahead = MOD(iproc+ii,numproc)
            ibehind = MOD(numproc+iproc-ii,numproc)
c
c--Send to node ahead, receive from node behind
c
#ifdef MPIDEBUG
            print *,iproc,': Sending xyzmh to ',iahead,' rec from ',
     &           ibehind,radkernel,xyzmh(1,1),numproc
#endif

            CALL MPI_SENDRECV(xyzmh,1,indexMPI5,iahead,0,
     &           xyzmh(1,2*ntot+3), idim, i5REAL8, ibehind,
     &           0, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)

#ifdef MPIDEBUG
            print *,iproc,': Received ',inumber,radkernel,xyzmh(1,1),
     &           numproc
#endif
c
c--This calculates the contributions to rhoi, gradhi, numneighreal, and 
c     gravity forces from non-neighbours
c
            inumreturn = 0
            DO i = 1, inumber

               ipart = i + ntot
               IF (ipart.GT.idim) THEN
                  WRITE (*,*) 'ipart.GT.idim'
                  CALL quit
               ENDIF
               iphase(ipart) = 0

               CALL densitygradh(ipart,npart,ntot,xyzmh,nneighlocal,
     &              dvxyzu)
c
c--Make a list of indices to return densities
c
               IF (nneighlocal.GT.0) THEN
                  inumreturn = inumreturn + 1
                  llistrec(inumreturn) = i
               ENDIF
            END DO
            ntotplusinumber = ntot + inumber

#ifdef MPIDEBUG
            print *,iproc,': number found neigh = ',inumreturn,numproc
#endif

            IF (iteration.EQ.1) THEN
c
c--nneightogetback does not record actual number, just 0 or 1.
c     Does not need the actual number, and getting the actual number would be
c     difficult because of the iterations to set h's and, thus, neighbours.
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

#ifdef MPIDEBUG
            print *,iproc,': sending indices to ',ibehind,inumreturn,
     &           numproc
#endif

            CALL MPI_SENDRECV(llistrec,inumreturn, MPI_INTEGER, ibehind,
     &           10, llistrec(inumreturn+1), nlstsend, MPI_INTEGER, 
     &           iahead, 10, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_INTEGER,inumberneigh,ierr)

#ifdef MPIDEBUG
            print *,iproc,': returned indices to ',ibehind,' and got ',
     &           inumberneigh,' from ',iahead,numproc
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
            print *,iproc,': sent rho ',numproc
#endif

            CALL MPI_SENDRECV(dumrho, 1, indexMPI1return, ibehind, 12,
     &           dumrhorec, inumberneigh, MPI_REAL4, iahead, 12,
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
            IF (icheck.NE.inumberneigh) THEN
               WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &              iproc
               CALL quit
            ENDIF

#ifdef MPIDEBUG
            print *,iproc,': sent gradh'
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
            IF (iteration.EQ.1) THEN
               IF (inumberneigh.GT.0) inumofsends = inumofsends + 1
            ELSE
               IF (inumberneigh.GT.0 .AND. 
     &              nneighsentback(iahead+1).EQ.0) THEN
                  inumofsends = inumofsends + 1
               ENDIF
            ENDIF
#ifdef MPIDEBUG
            print *,iproc,': all sent'
#endif
            DO k = 1, inumberneigh
               ipos = lsendlist(llistrec(k + inumreturn)) + 1
               rho(ipos) = rho(ipos) + rhorec(k)
               dumrho(ipos) = dumrho(ipos) + dumrhorec(k)
               nneigh(ipos) = nneigh(ipos) + nneighrec(k)
               IF (iteration.EQ.1) THEN
                  llistsentback(k,iahead+1) = ipos - 1
               ELSE
                  llistsentback(k+nneighsentback(iahead+1),
     &                 iahead+1) = ipos - 1
               ENDIF
#ifdef MPIDEBUG
               IF (k.EQ.1) print *,iproc,': POSITION ',xyzmh(1,ipos),
     &              xyzmh(2,ipos),xyzmh(3,ipos)
#endif
            END DO

#ifdef MPIDEBUG
            print *,iproc,': received neighbour data ',numproc
#endif
            IF (iteration.EQ.1) THEN
#ifdef MPIDEBUG
               print *,iproc,': reports ',inumberneigh,
     &              ' foreign neighbours on process ',iahead+1
#endif
               nneighsentbackold(iahead+1) = 0
               nneighsentback(iahead+1) = inumberneigh
            ELSE
               nneighsentbackold(iahead+1) = nneighsentback(iahead+1)
               nneighsentback(iahead+1) = nneighsentback(iahead+1) +
     &              inumberneigh
            ENDIF
            maxnneighsentback = MAX(maxnneighsentback,
     &           nneighsentback(iahead+1))
c
c--Send gravity forces back to sending process
c
#ifdef MPIDEBUG
            print *,iproc,' sending forces back to ',ibehind,numproc
#endif
            CALL MPI_SENDRECV(dvxyzu(1,ntot+1),4*inumber,
     &           MPI_REAL8,ibehind,14,
     &           dvxyzu(1,ntotplusinumber+1),4*nlstsend,MPI_REAL8,
     &                 iahead, 14, 
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
            IF (ireturned.NE.4*nlstsend) THEN
               WRITE (*,*) 'ERROR - ireturned.NE.4*nlstsend dvxyz'
               CALL quit
            ENDIF

            DO l = 1, nlstsend
               ipos = lsendlist(l) + 1
               jpos = ntotplusinumber + l

               DO k = 1, 4
                  dvxyzu(k,ipos) = dvxyzu(k,ipos) + dvxyzu(k,jpos)
               END DO

c               poten(ipos) = poten(ipos) + dvxyzu(4,jpos)
            END DO

#ifdef MPIDEBUG
            print *,iproc,' sent forces ',radkernel,xyzmh(1,1),numproc
#endif

         END DO
c
c--End circular send/receive
c
c
c--End of MPI-only section
c
#endif


         GOTO 777

      IF (iproc.EQ.0) THEN
         DO n = ntot+1, ntot+inumber
            ipart2 = n + ntot + 2
            ipart = n
            write (90,99901) xyzmh(1,ipart2),xyzmh(2,ipart2),
     &          xyzmh(3,ipart2),xyzmh(5,ipart2),rho(ipart),divv(ipart),
c     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
c99901                     FORMAT (11(1PE12.5,1X),I8)
         END DO
      ELSE
         DO n = ntot+1, ntot+inumber
            ipart2 = n + ntot + 2
            ipart = n
            write (91,99901) xyzmh(1,ipart2),xyzmh(2,ipart2),
     &          xyzmh(3,ipart2),xyzmh(5,ipart2),rho(ipart),divv(ipart),
c     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
         END DO
      ENDIF


      IF (iproc.EQ.0) THEN
         DO n = 1, nlstsend
            ipart = lsendlist(n) + 1
            write (92,99901) xyzmh(1,ipart),xyzmh(2,ipart),
     &          xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
c     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
c99901                     FORMAT (11(1PE12.5,1X),I8)
         END DO
      ELSE
         DO n = 1, nlstsend
            ipart = lsendlist(n) + 1
            write (93,99901) xyzmh(1,ipart),xyzmh(2,ipart),
     &          xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
c     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
         END DO
      ENDIF


 777  CONTINUE



         numsucceeded = 0
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,npart,ntot,list,divv,curlv,gradhs)
C$OMP& shared(hmax,xyzmh,vxyzu,pr,vsound,rho,ekcle,rhomin)
C$OMP& shared(nneigh,neighb,neighover,selfnormkernel)
C$OMP& shared(cnormk,radkernel,dvtable,ddvtable,wij,grwij)
C$OMP& shared(listpm,iphase,dphidh,uradconst,icall,encal)
C$OMP& shared(iprint,nptmass,iptmass,radcrit2,iorig,third)
C$OMP& shared(dumrho,iscurrent,ibound,ireal)
C$OMP& shared(isteps,it0,it1,imax,imaxstep,dt,itime)
C$OMP& shared(varmhd,Bevol,Bxyz,vsmooth,Bsmooth,nlstdo,ldolist)
C$OMP& shared(dvxyzu,hi_old,isucceed,poten,isend,numproc,iproc,gradmhd)
C$OMP& private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi)
C$OMP& private(divvi,curlvxi,curlvyi,curlvzi,gradhi,gradsofti)
C$OMP& private(pmassj,hi_oldi,hi1,hi21,hi31,hi41,hneigh)
C$OMP& private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2,rcut)
C$OMP& private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij)
C$OMP& private(projv,procurlvx,procurlvy,procurlvz)
C$OMP& private(l,iptcur,dphi,dwdhi,dpotdh,iteration,numneighi)
C$OMP& private(numneighreal,rhohi,dhdrhoi,omegai,func,dfdh1)
C$OMP& private(hnew,deltat,deltarho)
C$OMP& private(dalphaxi,dalphayi,dalphazi,dbetaxi,dbetayi,dbetazi)
C$OMP& private(denom,ddenom,bisection)
C$OMP& private(rxxi,rxyi,rxzi,ryyi,ryzi,rzzi)
C$OMP& private(gradalphaxi,gradalphayi,gradalphazi)
C$OMP& private(gradbetaxi,gradbetayi,gradbetazi,term)
C$OMP& private(vbarxi,vbaryi,vbarzi,hminbisec,hmaxbisec)
C$OMP& private(alphai,betai,wkern,dalpha,dbeta,grpmi)
C$OMP& private(nneighlocal,nneighi)
C$OMP& reduction(MAX:rhonext,imaxit)
C$OMP& reduction(+:inumit,inumfixed,inumrecalc,nwarnup,nwarndown)
C$OMP& reduction(+:nwarnroundoff,nbisection,numsucceeded)
      DO n = 1, nlstdo
         ipart = ldolist(n)
c
c--Don't need to exclude sink particles because already excluded
c
ccc         IF (iphase(ipart).GE.1) GOTO 50
c
c--Adds on local contributions to rho, gradh, and nneigh
c
         CALL densitygradh(ipart,npart,ntot,xyzmh,nneighlocal,
     &        dvxyzu)

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
c               write(iprint,*) 'bisection (',iteration,'): hnew =',hnew
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

         IF ((hnew.LE.0. .OR. omegai.LE.tiny .OR. 
     &        iteration.EQ.maxitsnr .OR. nneighi.LE.0) 
     &        .AND. isucceed(ipart).NE.-1) THEN
c
c--Switch to Bisection if not converging or running into trouble
c
            WRITE(iprint,*) 'WARNING: switching to bisection on'//
     &           ' particle ',iorig(ipart),'(',ipart,') hi = ',hi,
     &           ' hnew = ',hnew
            IF (nneighi.LE.0) THEN
               WRITE(iprint,*) '(particle has no neighbours)'
            ENDIF
            IF (omegai.LE.tiny) THEN
               WRITE(iprint,*) '(omega < tiny)'
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
         ELSEIF (ABS(hnew-hi)/hi_oldi.LT.htol) THEN
            imaxit = MAX(imaxit, iteration)
            inumit = inumit + iteration
            IF (isucceed(ipart).EQ.-1) THEN
               write(iprint,*) 'SWITCH BACK TO N-R ',ipart,hnew,hi
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
c
c--Pressure and sound velocity from ideal gas law...
c
               CALL eospg(ipart, vxyzu, rho, pr, vsound, ekcle)
c
c--Calculate other quantities using the interaction list from local process
c
               CALL divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
     &              nneighlocal,gradmhd)

c
c--Find particle with highest density outside radcrit of point mass
c
c               IF (iptmass.NE.0) THEN
c                  DO l = 1, nptmass
c                     iptcur = listpm(l)
c                     IF ( (xi - xyzmh(1,iptcur))**2 + 
c     &                    (yi - xyzmh(2,iptcur))**2 + 
c     &               (zi - xyzmh(3,iptcur))**2 .LT.radcrit2) GOTO 50
c                  END DO
c                  rhonext = MAX(rhonext, rho(ipart))
c               ENDIF
            ENDIF
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
      print *,iproc,': ifinished is ',ifinished,iscore,numproc
      print *,iproc,': maxnneighsentback is ',maxnneighsentback

      print *,iproc,': INFORMATION ',iteration
      DO i = 1, numproc
         print *,iproc,': ',i,nneighsentback(i),nneighsentbackold(i)
      END DO
      print *,iproc,': inumofsends ',inumofsends
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

#ifdef MPI
            IF (isend(ipart)) THEN
               nlstsend = nlstsend + 1
               lsendlist(nlstsend) = ipart - 1
            ENDIF
#endif
         ENDIF
      END DO
#ifdef MPI
c
c--Need to remove those particles that are being sent back from the lists of
c     particles being sent back to each MPI process since they will be added
c     back onto those lists if they have neighbours anyway (and we don't want
c     them done more than once!)
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
               IF (nneighsentback(i).EQ.0) inumofsends = inumofsends-1
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
     &     ' particles'
      CALL quit

 30   CONTINUE


      GOTO 778

      IF (iproc.EQ.0) THEN
         DO n = 1, nlstdo
            ipart = ldolist(n)
            write (88,99901) xyzmh(1,ipart),xyzmh(2,ipart),
     &           xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
c     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
99901              FORMAT (11(1PE12.5,1X),I8)
         END DO
      ELSE
         DO n = 1, nlstdo
            ipart = ldolist(n)
            write (89,99901) xyzmh(1,ipart),xyzmh(2,ipart),
     &           xyzmh(3,ipart),xyzmh(5,ipart),rho(ipart),divv(ipart),
c     &           gradhs(1,ipart),gradhs(2,ipart),
     &           dvxyzu(1,ipart),dvxyzu(2,ipart),dvxyzu(3,ipart),
     &           dvxyzu(4,ipart),poten(ipart),nneigh(ipart)
         END DO
      ENDIF



      STOP

 778  CONTINUE




c
c--For MPI job, need to add contributions to other quantities (such as divv,
c     curlv, etc) from other MPI processes
c
#ifdef MPI
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
c     is because here it is assumed that the total number of neighbours for
c     which other quantities need to be calculated is less than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (nneighsentback(j+1).GT.0) THEN
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
c     &                  1,llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),i3REAL8,indexMPI3,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)

                        CALL MPI_SEND(xyzmh,1,indexMPI5,j,40,
     &                       MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
                        print *,iproc,' sent xyzmh to ',j
#endif

                        CALL MPI_SEND(vxyzu,1,indexMPI4,j,41,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent vxyzu'
#endif
                        IF (imhd.EQ.idim) THEN
                           CALL MPI_SEND(Bevol,1,indexMPI3,j,43,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent Bevol'
#endif
                        ENDIF

                        CALL MPI_TYPE_FREE(indexMPI3,ierr)
                        CALL MPI_TYPE_FREE(indexMPI4,ierr)
                        CALL MPI_TYPE_FREE(indexMPI5,ierr)
#ifdef MPIDEBUG
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
               print *,iproc,': receiving neigh data ',inumberreturned,
     &              inumofreturns,radkernel,xyzmh(1,1)
#endif
               istart = ntot + inumbertotal + 1
               CALL MPI_RECV(xyzmh(1,istart+ntot+2), idim,i5REAL8,
     &              MPI_ANY_SOURCE, 40, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
               print *,iproc,': got xyzmh from ',iprocrec,inumber,
     &     ' to be put into ',istart,radkernel,xyzmh(1,1)

               print *,iproc,': POSITION BACK ',xyzmh(1,istart+ntot+2),
     &              xyzmh(2,istart+ntot+2),xyzmh(3,istart+ntot+2)
#endif

               IF (inumber.GT.0) THEN
                  CALL MPI_RECV(vxyzu(1,istart),inumber,i4REAL8,
     &                 iprocrec, 41, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i4REAL8, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 1 ',
     &                    iproc
                     CALL quit
                  ENDIF

                  IF (imhd.EQ.idim) THEN
                     CALL MPI_RECV(Bevol(1,istart),inumber,i3REAL8,
     &                    iprocrec, 43, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, i3REAL8, icheck, ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 3 ',
     &                       iproc
                        CALL quit
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got Bevol from ',iprocrec,icheck,
     &                    radkernel,xyzmh(1,1)
#endif
                  ENDIF
               ELSE
                  WRITE (*,*) 'ERROR - number returned <=0: ',inumber
                  CALL quit
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
      IF (inumbertotal.GT.idim) THEN
         WRITE (*,*) 'ERROR - inumbertotal.GT.idim ', iproc
         CALL quit
      ENDIF
      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR1 - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit
      ENDIF

#ifdef MPIDEBUG
      print *,iproc,': Received all particles with remote neighbours ',
     &     inumbertotal,radkernel,xyzmh(1,1)
      print *,' '
      print *,' '
#endif
c
c--Calculate ONLY neighbours lists for remote particles, then contributions
c     to divv, curlv, and Euler potential alpha and beta, etc
c
      DO j = 1, inumbertotal

         jpart = ntot + j
         list(nlst_end + j) = jpart
         iphase(jpart) = 0

         CALL treef(jpart,npart,ntot,xyzmh,acc,0,fsx,fsy,fsz,epot)

         numneighlocal = nneigh(jpart)

         CALL divvcurlvgrad(jpart,npart,ntot,xyzmh,vxyzu,Bevol,
     &        numneighlocal,gradmhd)
      END DO

#ifdef MPIDEBUG
      print *,iproc,': has calculated all divvcurlgrad ',inumbertotal
c         IF (iproc.EQ.2) print *,iproc,':G',dvxyzu(1,ntot + 1+2054-1),
c     &              dvxyzu(2,ntot+1+2054-1),dvxyzu(3,ntot+1+2054-1),
c     &     ntot+1+2054-1,ntot
#endif

c      STOP

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
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

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  divv(ipos) = divv(ipos) + divv(jpos)
               END DO

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

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  curlv(ipos) = curlv(ipos) + curlv(jpos)
               END DO

               CALL MPI_RECV(gradhs(1,istartrec),2*maxnneighsentback,
     &              MPI_REAL4,iprocrec,52,MPI_COMM_WORLD,istatus,ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (ireturned.NE.2*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsendback gradhs'
                  CALL quit
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got gradhs2 from ',iprocrec
#endif

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  gradhs(2,ipos) = gradhs(2,ipos) + gradhs(2,jpos)
               END DO
c
c--Receive back quantities for MHD gradients (if required)
c
               IF (imhd.EQ.idim) THEN

                  CALL MPI_RECV(gradmhd(1,istartrec),
     &                 12*maxnneighsentback,MPI_REAL8,
     &                 iprocrec,53,MPI_COMM_WORLD,istatus,ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
                  IF (ireturned.NE.12*nneighsentback(iprocrec+1)) THEN
                     WRITE (*,*) 'ERROR - ireturned.NE.nnsendbk gradmhd'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': got gradMHD from ',iprocrec
#endif

                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1
                     DO k = 1, 12
                        gradmhd(k,ipos)=gradmhd(k,ipos)+gradmhd(k,jpos)
                     END DO
                  END DO
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
                  print *,iproc,': sent divv,curlv,gradhs to ',
     &                 iprocsend,' and has sent ',inumber,' values'
#endif
c
c--Send back quantities for MHD gradients (if required)
c
                  IF (imhd.EQ.idim) THEN
                     CALL MPI_SEND(gradmhd(1,istart),inumber,i12REAL8,
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
      CALL MPI_TYPE_FREE(i12REAL8,ierr)

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

      DO j = 1, nlst_end
         ipart = list(j)
         IF (iphase(ipart).EQ.0) THEN
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
                  IF ( (xi - xyzmh(1,iptcur))**2 +
     &                 (yi - xyzmh(2,iptcur))**2 +
     &                 (zi - xyzmh(3,iptcur))**2 .LT.radcrit2) GOTO 50
               END DO
               rhonext = MAX(rhonext, rho(ipart))
            ENDIF
 50         CONTINUE

c
c--calculate B from the evolved Euler potentials
c
            IF (imhd.EQ.idim) THEN
               IF (varmhd.EQ.'eulr') THEN
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
c
c--grad alpha cross grad beta
c
                  term= gradalphayi*gradbetazi - gradalphazi*gradbetayi
                  Bxyz(1,ipart)= term
                  term= gradalphazi*gradbetaxi - gradalphaxi*gradbetazi
                  Bxyz(2,ipart)= term
                  term= gradalphaxi*gradbetayi - gradalphayi*gradbetaxi
                  Bxyz(3,ipart)= term
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
C$OMP& shared(vxyzu,gradhs,Bxyz,encal,ibound,ekcle,varmhd,list,Bevol)
C$OMP& shared(nlst_in,nlst_end)
C$OMP& private(i,j,ipart)
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
            IF (imhd.EQ.idim .AND. varmhd.EQ.'eulr') THEN
               Bxyz(1,i) = Bxyz(1,j)
               Bxyz(2,i) = Bxyz(2,j)
               Bxyz(3,i) = Bxyz(3,j)
            ENDIF
            IF (encal.EQ.'r' .AND. ibound.EQ.100) 
     &           ekcle(1,i) = ekcle(1,j)
         END DO
C$OMP END DO
c
c--Calculate B from the evolved magnetic field variable
c
      IF (imhd.EQ.idim) THEN
         IF (varmhd.eq.'Bvol') THEN
C$OMP DO SCHEDULE (runtime)
            DO i=nlst_in,nlst_end
               ipart = list(i)
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
               ipart = list(i)
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
         ELSEIF (varmhd.NE.'eulr') THEN
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
         WRITE (iprint,*) 'WARNING: denom in euler gradients zero on ',
     &        nwarnroundoff,' particles'
      ENDIF
      IF (nbisection.GT.0) THEN
         WRITE(iprint,*) 'WARNING: used bisection on ',nbisection,
     &       ' particles'
      ENDIF

c
c--Possible to create a point mass
c
      IF (rhonext.GT.rhocrea .AND. nptmass.LT.iptdim .AND.
     &                     iptmass.NE.0 .AND. icall.EQ.3) THEN
c
c--Find particle with highest density outside radcrit of point mass
c
         irhonex = 0
         DO n = nlst_in, nlst_end
            ipart = list(n)
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
         IF ((2.0*xyzmh(5,irhonex)).LT.hacc) THEN
            WRITE(iprint,*)'Ptmass creation passed h ',xyzmh(5,irhonex)

            CALL getneigh(irhonex,npart,xyzmh(5,irhonex),xyzmh,nlist,
     &           iptneigh,nearl)

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
               WRITE(iprint,*) ' and all particles on step'
            ELSE
               WRITE(iprint,*) ' but not all particles on step'
            ENDIF
         ELSE
            WRITE(iprint,*) 'Ptmass creation failed on h ',
     &           xyzmh(5,irhonex)
         ENDIF
      ENDIF

      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine densityiterate')

      RETURN
      END



c=============================================================================
c
c--Subroutine for calculating rho, gradh and number of neighbours
c
      SUBROUTINE densitygradh(ipart,npart,ntot,xyzmh,numneighi,dvxyzu)

      INCLUDE 'idim'

      DIMENSION xyzmh(5,mmax), dvxyzu(4,idim)

      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/dumderivi'

      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
c
c--Needed for MPI code
c
      IF (ipart.GT.npart) THEN
         iparttree = ipart + ntot + 2
      ELSE
         iparttree = ipart
         nneighprior = nneigh(ipart)
      ENDIF

      xi = xyzmh(1,iparttree)
      yi = xyzmh(2,iparttree)
      zi = xyzmh(3,iparttree)
      pmassi = xyzmh(4,iparttree)
      hi = xyzmh(5,iparttree)

      CALL treef(ipart,npart,ntot,xyzmh,acc,igphi,fsx,fsy,fsz,epot)
c
c--Calculate density contribution by looping over interacting neighbors
c
      numneighi = nneigh(ipart)
      rhoi = 0.
      gradhi = 0.
      nneighi = 0
      hi1 = 1./hi
      hi21 = hi1*hi1
      hi31 = hi21*hi1
      hi41 = hi21*hi21

      IF (numneighi.GT.0) THEN

         DO k = 1, numneighi
            j = neighlist(k)

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
         END DO
      ENDIF
c
c--Needed for MPI code
c
      IF (ipart.GT.npart) THEN
         rho(ipart) = cnormk*rhoi
         dumrho(ipart) = cnormk*gradhi
         nneigh(ipart) = numneighi
         dvxyzu(1,ipart) = fsx
         dvxyzu(2,ipart) = fsy
         dvxyzu(3,ipart) = fsz
         dvxyzu(4,ipart) = epot
      ELSE
         rho(ipart) = rho(ipart) +
     &        cnormk*(rhoi + selfnormkernel*pmassi*hi31)
         dumrho(ipart) = dumrho(ipart) + cnormk*(gradhi +
     &        selfnormkernel*pmassi*(-3.*hi41))
         nneigh(ipart) = nneighprior + numneighi
         dvxyzu(1,ipart) = dvxyzu(1,ipart) + fsx
         dvxyzu(2,ipart) = dvxyzu(2,ipart) + fsy
         dvxyzu(3,ipart) = dvxyzu(3,ipart) + fsz
         dvxyzu(4,ipart) = dvxyzu(4,ipart) + epot
      ENDIF

      RETURN
      END

c=============================================================================
c
c--Subroutine for calculating other quantities such as divv, curlv, alpha,
c     and beta for Euler potentials, etc
c
      SUBROUTINE divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
     &     numneighi,gradmhd)

      INCLUDE 'idim'

      DIMENSION xyzmh(5,mmax), vxyzu(4,idim), Bevol(3,imhd),
     &     gradmhd(12,imhd)

      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
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
         DO k = 1, 12
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
      IF (imhd.EQ.idim) THEN
         IF (varmhd.EQ.'eulr') THEN
            alphai= Bevol(1,ipart)
            betai= Bevol(2,ipart)
         ENDIF
      ENDIF
                  
      DO k = 1, numneighi
         j = neighlist(k)

         dx = xi - xyzmh(1,j)
         dy = yi - xyzmh(2,j)
         dz = zi - xyzmh(3,j)
#ifdef PERIODIC_NO_GHOSTS
         IF (varmhd.EQ.'eulr') THEN
            dalpha= alphai - Bevol(1,j)
            dbeta= betai - Bevol(2,j)
            CALL modboundeulr(dx,dy,dz,dalpha,dbeta)
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
c--Velocity divergence times density
c
            projv = grwtij*(dvx*dx + dvy*dy + dvz*dz)
            divvi = divvi - pmassj*projv
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
c--Derivative of gravitational potential w.r.t. h
c
            gradsofti = gradsofti - pmassj*dphi
c
c--get B from the evolved Euler potentials
c
c
c--grad alpha and grad beta (Euler potentials)
c
            IF (imhd.EQ.idim) THEN
               IF (varmhd.EQ.'eulr') THEN
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
      END DO
c
c--Need for MPI
c
      IF (ipart.GT.npart) THEN
         divv(ipart) = cnormk*divvi
         curlv(ipart) = cnormk*SQRT(curlvxi**2+curlvyi**2+curlvzi**2)
         gradhs(2,ipart) = gradsofti
      ELSE
         divv(ipart) = divv(ipart) + cnormk*divvi
         curlv(ipart) = curlv(ipart) + 
     &        cnormk*SQRT(curlvxi**2+curlvyi**2+curlvzi**2)
c
c--Add self contribution
c
         gradhs(2,ipart) = gradhs(2,ipart) +  gradsofti 
     &        - pmassi*dphidh(0)*hi21
      ENDIF

      RETURN
      END
