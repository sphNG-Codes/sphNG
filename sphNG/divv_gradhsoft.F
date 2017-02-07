      SUBROUTINE divv_gradhsoft (dt,npart,ntot,xyzmh,vxyzu,
     &            nlst_in,nlst_end,listp,itime,ekcle,Bevol,Bxyz)
c************************************************************
c                                                           *
c  Subroutine to compute the velocity divergence, the curl, *
c     and the gravity softening term with variable h, and   *
c     interpolates the density for particles which are      *
c     neighbours of particles that are currently in the     *
c     list.                                                 *
c                                                           *
c     NOTE: Only the magnitude of curlv is calculated, and  *
c     for MPI jobs it is NOT calculated correctly !!!       *
c                                                           *
c     Essentially, it calculates any quantity except        *
c     density and smoothing length and dh/dt = gradhs(1,*)  *
c     that is required BEFORE the forces are calculated     *
c                                                           *
c     This used to all be done in densityiterate, but was   *
c     split off into a separate subroutine by MRB on        *
c     12/05/2015                                            *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2)
      DIMENSION listp(idim2)
      DIMENSION ekcle(5,iradtrans2)
      DIMENSION Bevol(imhdevol,imhd3),Bxyz(3,imhd2)

      INCLUDE 'COMMONS/astrcon'
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
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/sightlines'

#ifdef MPICOPY
      INCLUDE 'COMMONS/mpicopy'
#endif

      COMMON /rhominval/ rhomin

      DIMENSION gradmhd(15,imhd2)

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine divv_gradhsoft')
c
c--Initialise
c
#ifdef MPIALL
#ifdef MPIDEBUG
      print *,iproc,': Entered divv_gradhsoft ',itime,nlst_end
#endif
#else
c      print *,': Entered divv_gradhsoft ',itime
#endif

      IF (itiming) CALL getused(td31)

#ifdef MPI
      CALL MPI_TYPE_CONTIGUOUS(15, MPI_REAL8, i15REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)

      CALL MPI_TYPE_COMMIT(i15REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)
#endif
c
c--Now calculate local contributions to divv, curlv, etc
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end)
C$OMP& shared(xyzmh,vxyzu,pr,vsound,rho,ekcle,Bevol)
C$OMP& shared(iphase,icall,encal)
C$OMP& shared(numproc,iproc,gradmhd,npart,ntot)
C$OMP& shared(ivar,ijvar)
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend)
#endif
C$OMP& private(n,ipart,k)
C$OMP& private(nneighlocal)
#ifdef MPICOPY
      DO n = numberstart, numberend
#else
      DO n = nlst_in, nlst_end
#endif
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

#ifdef MPICOPY
c
c--For MPICOPY job, need to transfer gradhs(1-2) between processes.  This 
c     needs to be done here because gradhs(1) is required for calculating
c     the final values of divv and curlv
c
c     Also need to transfer the other products of divvcurlvgrad which
c     are: divv(), curlv(), and gradmhd(1-15,i)
c
c--Transfer gradhs(1-2)
c
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)

      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
     &     lblocklengths, llisttrans(numberstart),
     &     i2REAL4, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their gradhs(1-2) data
c
      CALL MPI_ALLGATHERV(gradhs,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,i2REAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put gradhs(1-2) data into correct places in gradhs(1-2)
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,real4transfer2)
C$OMP&shared(gradhs,numproc,iproc)
C$OMP&private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 2
                  gradhs(k,llisttrans(i) + 1) = 
     &                 real4transfer2((i-1)*2 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': gradhs(1-2) move '
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
      CALL MPI_TYPE_FREE(i2REAL4,ierr)
c
c--Transfer divv() and curlv()
c
      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
     &     lblocklengths, llisttrans(numberstart),
     &     MPI_REAL4, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their divv() data
c
      CALL MPI_ALLGATHERV(divv,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,MPI_REAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put divv() data into correct places in divv()
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,real4transfer2)
C$OMP&shared(divv,numproc,iproc)
C$OMP&private(i,j)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               divv(llisttrans(i) + 1) = 
     &              real4transfer2(i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': divv() move '
#endif
c
c--All processes transfer their curlv() data
c
      CALL MPI_ALLGATHERV(curlv,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,MPI_REAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put curlv() data into correct places in curlv()
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,real4transfer2)
C$OMP&shared(curlv,numproc,iproc)
C$OMP&private(i,j)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               curlv(llisttrans(i) + 1) = 
     &              real4transfer2(i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': curlv() move '
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)

      IF (imhd.EQ.idim) THEN
         IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
c
c--Transfer gradmhd(1-15)
c
      CALL MPI_TYPE_CONTIGUOUS(15, MPI_REAL8, i15REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i15REAL8,ierr)

      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
     &     lblocklengths, llisttrans(numberstart),
     &     i15REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their gradmhd(1-15) data
c
      CALL MPI_ALLGATHERV(gradmhd,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,i15REAL8,
     &     MPI_COMM_WORLD,ierr)
c
c--Put gradmhd(1-15) data into correct places in gradmhd(1-15)
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,realtransfer5to15)
C$OMP&shared(gradmhd,numproc,iproc)
C$OMP&private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 15
                  gradmhd(k,llisttrans(i) + 1) = 
     &                 realtransfer5to15((i-1)*15 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': gradmhd(1-15) move '
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
      CALL MPI_TYPE_FREE(i15REAL8,ierr)

         ENDIF
      ENDIF
#endif
c
c--END OF MPICOPY SECTION
c

#ifdef MPI
c
c--For MPI job, need to add contributions to other quantities (such as divv,
c     curlv, etc) from other MPI processes
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumbertotal,ntot)
C$OMP& shared(npart,xyzmh,vxyzu,Bevol,gradmhd)
C$OMP& shared(ncompactlocal,ivar,ijvar)
C$OMP& private(nneighlocal,ncompacthere)
C$OMP& private(i,ipart)
      DO i = 1, inumbertotal
         ncompacthere = ncompactlocal + i
         ipart = ivar(3,ncompacthere)
         nneighlocal = ivar(1,ncompacthere)

         CALL divvcurlvgrad(ipart,npart,ntot,xyzmh,vxyzu,Bevol,
     &        nneighlocal,ijvar(ivar(2,ncompacthere)+1),gradmhd)
      END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
      print *,iproc,': has calculated all divvcurlgrad ',inumbertotal
#endif
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
c--Receive back forces, du, potential energy to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' divv contributions ',j,
     &              ' with maxnneighsentback ',maxnneighsentback
#endif
               CALL MPI_RECV(divv(istartrec),maxnneighsentback,
     &              MPI_REAL4,MPI_ANY_SOURCE,50,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, ireturned, ierr)
               IF (istartrec+ireturned.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istartrec+ireturned.GT.idim2 ',
     &                 istartrec,ireturned
                  CALL quit(1)
               ENDIF
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsentback divv ',
     &                 ireturned,nneighsentback(iprocrec+1)
                  CALL quit(1)
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
                  CALL quit(1)
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
                  CALL quit(1)
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
                     CALL quit(1)
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
                   print *,iproc,': sending forces to ',iprocsend,
     &                 ' istart ',istart,' number ',inumber,
     &                 ' returned ',inumberreturned,
     &                 ' of ',inumofreturns,xyzmh(1,1)
#endif
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
         CALL quit(1)
      ENDIF

      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i15REAL8,ierr)

#ifdef MPIDEBUG
           print *,iproc,': FINISHED ',itime
#endif
#endif

c Initialise nwarnroundoff counter to zero
      nwarnroundoff = 0


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
               ENDIF
            ENDIF
         ENDIF
      END DO
C$OMP END PARALLEL DO

c Emit warning if nwarnroundoff > 0
c Used to be done in densityiterate
      IF (nwarnroundoff.GT.0) THEN
         IF (varmhd.EQ.'vecp') THEN
           WRITE (iprint,*) 'WARNING: denom in vecp gradients zero on ',
     &        nwarnroundoff,' particles'
         ELSEIF (varmhd.EQ.'eulr') THEN
          WRITE (iprint,*) 'WARNING: denom in euler gradients zero on ',
     &        nwarnroundoff,' particles'
         ENDIF
      ENDIF


C$OMP PARALLEL default(none)
C$OMP& shared(npart,ntot,ireal,rho,dumrho,xyzmh,pr,vsound,divv,curlv)
C$OMP& shared(vxyzu,gradhs,Bxyz,encal,ibound,ekcle,varmhd,listp,Bevol)
C$OMP& shared(nlst_in,nlst_end,listparents,nlistupdated,iupdated,icall)
C$OMP& shared(dt,itime,iphase,dust_tk)
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
         dumrho(i) = dumrho(j)
         pr(i) = pr(j)
         vsound(i) = vsound(j)
         divv(i) = divv(j)
         curlv(i) = curlv(j)
         vxyzu(4,i) = vxyzu(4,j)
         gradhs(2,i) = gradhs(2,j)
         IF (imhd.EQ.idim) THEN
            IF (varmhd.EQ.'eulr' .or. varmhd.EQ.'vecp') THEN
               Bxyz(1,i) = Bxyz(1,j)
               Bxyz(2,i) = Bxyz(2,j)
               Bxyz(3,i) = Bxyz(3,j)            
            ENDIF
         ENDIF
         IF (encal.EQ.'r') THEN
            IF (ibound/10.EQ.10) THEN
               IF(icall.EQ.1) THEN
                  DO k=1,5
                     ekcle(k,i) = ekcle(k,j)
                  END DO
               ELSE
                  ekcle(1,i) = ekcle(1,j)
               ENDIF
            ENDIF
            IF (idustRT.GT.0) THEN
               dust_tk(1,i) = dust_tk(1,j)
               dust_tk(2,i) = dust_tk(2,j)
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

      IF (itiming) THEN
         CALL getused(td32)
         td3 = td3 + (td32 - td31)
      ENDIF

      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine divv_gradhsoft')

c      print *,'Exited divv_gradhsoft ',itime

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

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2), Bevol(imhdevol,imhd3),
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
            CALL quit(1)
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
