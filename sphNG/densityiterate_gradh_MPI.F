      SUBROUTINE densityiterate_gradh (dt,npart,ntot,xyzmh,vxyzu,dvxyzu,
     &            nlst_in,nlst_end,listp,itime,ekcle,Bevol,Bxyz,dBxyz,
     &            dustvar,ddust,rgrain)
c************************************************************
c                                                           *
c  Subroutine to compute the density and smoothing lengths  *
c     self-consistently using iteration if necessary.       *
c     This subroutine uses the binary tree algorithm        *
c     to locate neighbours.                                 *
c                                                           *
c     Code written by MRB and DJP (14/12/2005).             *
c     MPI version written by MRB (26/06/2007).              *
c                                                           *
c     On 12/05/2015, original version was split into two    *
c     subroutines, with the second half going in to         *
c     divv_gradhsoft.F                                      *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif

#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidomains'
#endif

      INCLUDE 'COMMONS/accstat'

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2), dvxyzu(4,idim3)
      DIMENSION listp(idim2)
      DIMENSION ekcle(5,iradtrans2)
      DIMENSION Bevol(imhdevol,imhd3),Bxyz(3,imhd2),
     &     dBxyz(imhdevol,imhd2),dustvar(ndusttypes,idim_dustFluid2),
     &     ddust(ndusttypes,idim_dustFluid2)
      DIMENSION rgrain(idim_grow2)

      PARAMETER (htol = 1.e-3)
      PARAMETER (hstretch = 1.01)
      PARAMETER (maxiterations = 500)
      PARAMETER (maxitsnr = 30)

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
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/sightlines'
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/dustfluidvelu'
      INCLUDE 'COMMONS/tstopvar'
      INCLUDE 'COMMONS/HY09accel'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif

      COMMON /rhominval/ rhomin

#ifdef MPICOPY
      INCLUDE 'COMMONS/mpicopy'
#endif

      INTEGER*1 isucceed(idim)
      REAL*4 hminbisec(idim),hmaxbisec(idim)

      DIMENSION ldolist(idim),hi_old(idim),hi_2back(idim)
      LOGICAL*1 iconvergence1,iconvergence2,iconvergence3
      LOGICAL ifound
#ifdef MPI
      DIMENSION nneighrec(ineighproc),nneighsentbackold(nummaxproc)
      REAL*4 rhorec(ineighproc), dumrhorec(ineighproc), rhonextmax
      LOGICAL*1 isend(idim)
      INTEGER*8 iuniquesend(iptdim)

      REAL de_all(6,numproc),csgas_all(6,numproc)
      REAL csbuff(numproc),dvxyzuinit(4,idim)
      INTEGER idispranks(numproc),ipsendbuf(numproc*numproc)
      INTEGER ipsendall(numproc,numproc),ipsendindx(numproc)
      INTEGER listdistantgn(ngrav_nodes),listsendgn(ngrav_nodes)
      LOGICAL order_comms
#endif

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine densityiterate')
c
c--Initialise
c
      iproblem = itime
#ifdef MPIALL
#ifdef MPIDEBUG
      print *,iproc,': Entered densityiterate ',itime,nlst_end
#endif
#else
c      print *,': Entered densityiterate ',itime
#endif

      IF (itiming) CALL getused(td11)

      distancemax = MAX(rmax, xmax, ymax, zmax, rcyl)

      uradconst = radconst/uergcc
      third = 1./3.
      rhonext = 0.
      icreate = 0
      icreatetot = 0
      radcrit2 = radcrit*radcrit
      numparticlesdone = numparticlesdone + nlst_end
      nwarnup = 0
      nwarndown = 0
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
            CALL quit(1)
         ENDIF
      ENDIF
#endif
#ifdef MPICOPY
#ifdef MPIDEBUG
      print *,iproc,' nlstdo ',nlstdo
#endif

      numbertodo = nlstdo/numproc
      numbertodohere = numbertodo
      IF (MOD(nlstdo,numproc).GT.iproc) 
     &     numbertodohere = numbertodohere + 1
      idisplacementsi(1) = 0
      DO i = 1, numproc
         irecvcounti(i) = numbertodo
         IF (MOD(nlstdo,numproc).GT.i-1)
     &        irecvcounti(i) = irecvcounti(i) + 1
      END DO
      idisplacementsi(1) = 0
      DO i = 1, numproc - 1
         idisplacementsi(i+1) = idisplacementsi(i) + irecvcounti(i)
      END DO

#ifdef MPIDEBUG
      print *,iproc,': recv ',(irecvcounti(i),i=1,numproc)
      print *,iproc,': disp ',(idisplacementsi(i),i=1,numproc)
#endif

      numberstart = idisplacementsi(iproc+1) + 1
      numberend = numberstart + numbertodohere - 1
      numberendkeep = numberend

      IF (MOD(nlstdo,numproc).NE.0) THEN
         numbertodomax = numbertodo + 1
      ELSE
         numbertodomax = numbertodo
      ENDIF

#ifdef MPIDEBUG
      print *,iproc,': numbertodomax ',numbertodo,numbertodomax,
     &     numberstart,numberend,numbertodohere
#endif

      numbertodoherekeep = numbertodohere
C$OMP PARALLEL DO SCHEDULE(static) default(none)
C$OMP& shared(nlstdo,ldolist,llisttrans)
C$OMP& private(i)
      DO i = 1, nlstdo
         llisttrans(i) = ldolist(i) - 1
      END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
      print *,iproc,': nlstdo,todo,start,end ',
     &     nlstdo,numbertodo,numberstart,numberend
#endif
      ncompact = ncompact + numberstart - 1
#endif

      nlstsend = nlstdo

#ifdef MPI
      CALL MPI_TYPE_CONTIGUOUS(15, MPI_REAL8, i15REAL8, ierr)
c      CALL MPI_TYPE_CONTIGUOUS(12, MPI_REAL8, i12REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(8, MPI_REAL8, i8REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(imhdevol, MPI_REAL8, imhdevolREAL8,ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(7, MPI_REAL4, i7REAL4, ierr)

      CALL MPI_TYPE_COMMIT(i15REAL8,ierr)
c      CALL MPI_TYPE_COMMIT(i12REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i8REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(imhdevolREAL8,ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)
      CALL MPI_TYPE_COMMIT(i7REAL4,ierr)

c--domain checking optimisations
c     before iteration begins need to establish which
c     MPI ranks this rank will send particles to
c     others will receive gravity nodes instead (if igphi.EQ.1)
      IF (domain_checking) THEN
c--need to reverse the sign of the second half of domain extent arrays
c     they are stored this way so that C$OMP& reduction(MIN:
c     applies to the entirety of each array in step.
         DO k=4,6
            domain_extent(k) = -domain_extent(k)
            de_cs_gas(k) = -de_cs_gas(k)
         ENDDO
         
         CALL MPI_ALLGATHER(domain_extent,6,MPI_REAL8,
     &        de_all,6,MPI_REAL8,MPI_COMM_WORLD,ierr)
         CALL MPI_ALLGATHER(de_cs_gas,6,MPI_REAL8,
     &        csgas_all,6,MPI_REAL8,MPI_COMM_WORLD,ierr)

         nprocsend = 0
         ipsendindx = 0
         nldistantgn = 0
         instart = 2*ntot+3
c--pairwise sendrecv
         DO ii=1,numproc-1
            ipair = ipairwiseall(ii)
            i = ipair + 1
            cdist2 = MIN(cuboid_dist2(de_cs_gas,de_all(1:6,i)),
     &           cuboid_dist2(domain_extent,csgas_all(1:6,i)))
            IF (cdist2 .LT. tiny) THEN
c--neighbouring rank so add it to list to send particles to
               nprocsend = nprocsend + 1
               ipsendindx(nprocsend) = i
               CYCLE
            ELSEIF (igphi .EQ. 1) THEN
               idatsend = 0
c--otherwise send gravity nodes in other direction
               CALL grav_only_treef(de_all(1:6,i),xyzmh,acc,
     &              npartsend,nbothsend,listsendgn(:),0.0)
               IF (npartsend .GT. 0) THEN
c--cannot send particle data so send all particless on rank
c     this case probably arise with diagonally adjacent domains
                  idatsend = 1
#ifdef MPIDEBUGDC
                  PRINT*,iproc,": WARNING - sending extra rank",i
     &                 npartsend,nbothsend
#endif
               ELSEIF (nbothsend .LE. 0) THEN
                  PRINT*,iproc,": ERROR - must send at least one node",
     &                 npartsend,nbothsend
                  CALL quit(1)
               ELSEIF (instart+nbothsend.GE.mmax+ngrav_nodes) THEN
                  PRINT*,': ERROR - ngrav_nodes too small',
     &                 ntot,nbothrec
                  CALL quit(1)
               ELSEIF (nbothsend.GE.30) THEN
                  WRITE(iprint,*)"sending 30+gnodes",nbothsend
               ENDIF

               CALL MPI_SENDRECV(idatsend,1,MPI_INTEGER,ipair,29,
     &              idatrec, MPI_INTEGER, i5REAL8,
     &              ipair, 29, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,i5REAL8,nbothrec,ierr)

               IF (idatsend .EQ. 1 .OR. idatrec .EQ. 1) THEN
                  nprocsend = nprocsend + 1
                  ipsendindx(nprocsend) = i
                  CYCLE
               ELSE
                  CALL insert_sort(nbothsend,listsendgn(1:nbothsend))

                  CALL MPI_TYPE_INDEXED(nbothsend,
     &                 lblocklengths,listsendgn,
     &                 i5REAL8,indexMPI5gn,ierr)
                  CALL MPI_TYPE_COMMIT(indexMPI5gn,ierr)

                  CALL MPI_TYPE_INDEXED(nbothsend,
     &                 lblocklengths,listsendgn,
     &                 i7REAL4,indexMPI_qrad,ierr)
                  CALL MPI_TYPE_COMMIT(indexMPI_qrad,ierr)

                  CALL MPI_SENDRECV(xyzmh,1,indexMPI5gn,ipair,28,
     &                 xyzmh(1,instart), ngrav_nodes, i5REAL8,
     &                 ipair, 28, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,i5REAL8,nbothrec,ierr)

                  CALL MPI_SENDRECV(qrad,1,indexMPI_qrad,ipair,27,
     &                 qrad(1,instart), ngrav_nodes, i7REAL4,
     &                 ipair, 27, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,i7REAL4,nreccheck,ierr)
                  IF (nbothrec.NE.nreccheck) THEN
                     PRINT*,iproc,': ERROR - nbothrec.NE.nreccheck',
     &                    nbothrec,nreccheck
                     CALL quit(1)
                  ENDIF

                  CALL MPI_TYPE_FREE(indexMPI5gn,ierr)
                  CALL MPI_TYPE_FREE(indexMPI_qrad,ierr)

                  DO nld=nldistantgn+1,nldistantgn+nbothrec
                     listdistantgn(nld) = instart
                     instart = instart + 1
                  ENDDO
                  nldistantgn = nldistantgn + nbothrec
               ENDIF
            ENDIF
         ENDDO
         CALL insert_sort(nprocsend,ipsendindx(1:nprocsend))

c--reverse second half of domain extent arrays again         
         DO k=4,6
            domain_extent(k) = -domain_extent(k)
            de_cs_gas(k) = -de_cs_gas(k)
         ENDDO

c--have to re-order MPI comms if any rank has changed its send list         
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
     &        MPI_INTEGER,MPI_LOR,MPI_COMM_WORLD,ierr)

         IF (order_comms) THEN
            CALL MPI_ALLGATHER(nprocsend,1,MPI_INTEGER,
     &        nprocsendall,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            idispranks(1) = 0
            DO i=1,numproc-1
               idispranks(i+1) = idispranks(i) + nprocsendall(i)
            ENDDO
            CALL MPI_ALLGATHERV(ipsendindx,nprocsend,MPI_INTEGER,
     &           ipsendbuf,nprocsendall,idispranks,
     &           MPI_INTEGER,MPI_COMM_WORLD,ierr)

c--build array of which ranks are sending to which
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
               CALL order_mpi_comms_bal(ipsendall)
            ELSE
               CALL order_mpi_comms_else(ipsendall,ipsendindx)
            ENDIF
         ENDIF
      ENDIF
#endif

>>>>>>> da3df6b... Removed infrastructure for treating domain gas and star/dust domain extents separately.
C$OMP PARALLEL DO SCHEDULE(runtime) default(none) 
C$OMP& shared(nlstdo,ldolist,isucceed,xyzmh,hi_old,it1,isteps,imaxstep)
C$OMP& shared(divv,rho,dumrho,nneigh,dvxyzu,rhomin,dt,imax,poten,dq)
C$OMP& shared(iphase,dBxyz,divcurlB,gradB,it0,itime,hfact)
C$OMP& shared(curlv,gradhs,sightcolumns,hi_2back)
C$OMP& shared(planetesimaltimestep,planetesimalnorm,iplanetesimals)
C$OMP& shared(stellarrad,nptmasstot,ddust,dustfluiddu,dustfluidvxyz)
C$OMP& shared(stoppingtime,gas_accel)
#ifdef MPI
C$OMP& shared(lsendlist,isend)
#endif
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend)
#endif
C$OMP& private(n,ipart,pmassi,hi,dhdrhoi,deltat,k,hchange)
c
c--Loop is the same for MPICOPY and normal code because some variables
c     need to be zeroed both those particles done remotely and by this
c     MPI process.
c
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
            hi_2back(ipart) = hi
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
               CALL quit(1)
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
            DO k = 1, imhdevol
               dBxyz(k,ipart) = 0.
            END DO
            DO k = 1, 5
               divcurlB(k,ipart) = 0.
            END DO
            DO k = 1, 9
               gradB(k,ipart) = 0.
            END DO
         ENDIF
         IF (idustFluid.NE.0) THEN
            ddust(:,ipart) = 0.
            dustfluiddu(ipart) = 0.
            stoppingtime(:,ipart) = 0.
         ENDIF
         IF (idustFluid.NE.0 .OR. imakedust) THEN
            dustfluidvxyz(1:3,:,ipart) = 0.
         ENDIF
         IF (idimHY09.EQ.idim) THEN
            gas_accel(:,ipart) = 0.
         ENDIF
         IF (idustRT.GT.0) THEN
#ifdef MPI
            sightcolumns(1:nsightlinesmax,1:2,ipart) = 0.
#endif
         ENDIF
         IF (istellar_radiation) THEN
            stellarrad(1:4,1:nptmasstot,ipart) = 0.
         ENDIF
         IF (iplanetesimals.GT.0) THEN
            planetesimaltimestep(ipart) = 0.0
            planetesimalnorm(ipart) = 0.0
         ENDIF
      END DO
C$OMP END PARALLEL DO

#ifdef MPI
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
      CALL MPI_TYPE_CONTIGUOUS(8, MPI_REAL8, i8REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(imhdevol, MPI_REAL8, imhdevolREAL8,ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)

      CALL MPI_TYPE_COMMIT(i15REAL8,ierr)
c      CALL MPI_TYPE_COMMIT(i12REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i8REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(imhdevolREAL8,ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)

      CALL MPI_ALLREDUCE(nlstdo,nlstdo_tot,1,MPI_INTEGER,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
#else
      nlstdo_tot = nlstdo
#endif

      IF (istellarfeedback.GT.0 .OR. istellar_radiation) THEN
c
c--For stellar radiation from sink particles, need to know about all
c     sink particles on all MPI processes.  For non-MPI, still keep
c     this in the same array for consistency.
c      
         DO ipt = 1, nptmass
            stellar_xyzmhrti(1:5,ipt) = xyzmh(1:5,listpm(ipt))
            IF (istellar_radiation) THEN
               stellar_xyzmhrti(6:8,ipt) = stellar_radtempion(1:3,ipt)
            ELSE
               stellar_xyzmhrti(6:8,ipt) = 0.
            ENDIF
         END DO
#ifdef MPI
c
c--Get iproc=0 to gather all sink information, then broadcast.  This has
c     the advantage that all MPI processes have the same order
c
         IF (iproc.EQ.0) THEN
            ipos = nptmass + 1
            IF (ipos.GT.iptdim) THEN
               WRITE (*,*) 'ERROR - ipos.GT.iptdim'
               CALL quit(1)
            ENDIF
            DO i = 1, numproc-1
               CALL MPI_RECV(stellar_xyzmhrti(1,ipos), 8*iptdim,
     &              MPI_REAL8, MPI_ANY_SOURCE, 218, MPI_COMM_WORLD,
     &              istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               IF (MOD(inumber,8).NE.0) THEN
                  WRITE (*,*) iproc,': ERROR - 218 MOD(inumber,8).NE.0',
     &              i,inumber,iptdim,ipos
                  CALL quit(1)
               ENDIF
               ipos = ipos + inumber/8
            END DO
            IF (ipos-1.NE.nptmasstot) THEN
               WRITE (*,*) iproc,': ERROR - ipos-1.NE.nptmasstot',
     &              ipos-1,nptmasstot
               CALL quit(1)
            ENDIF
         ELSE
            CALL MPI_SEND(stellar_xyzmhrti(1,1),nptmass,i8REAL8,0,218,
     &                       MPI_COMM_WORLD, ierr)
         ENDIF
         CALL MPI_BCAST(stellar_xyzmhrti(1,1),nptmasstot,i8REAL8,0,
     &        MPI_COMM_WORLD, ierr)
c
c--Find local sink particle in the global list
c
         DO ipt = 1, nptmass
            DO jpt = 1, nptmasstot
               ifound = .TRUE.
               DO k = 1, 3
                  ifound = (ifound .AND. 
     &      (ABS(xyzmh(k,listpm(ipt))-stellar_xyzmhrti(k,jpt)).LT.tiny))
               END DO
               IF (ifound) THEN
                  liststellarID(ipt) = jpt
                  GOTO 400
               ENDIF
            END DO
 400        CONTINUE
         END DO
#endif
      ENDIF
c
c--Iterate density calculation for particle ipart
c
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

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nlstsend, 1, lsendlist, 
c     &       MPI_INTEGER, indexMPI_INT, ierr)
         CALL MPI_TYPE_INDEXED(nlstsend, lblocklengths, lsendlist, 
     &        MPI_INTEGER, indexMPI_INT, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI_INT,ierr)
c
c--There are two implementations, one for balanced jobs, one for unbalanced
c
         nptmasslocal = nptmass
         nptmasslocal2 = nptmass
         DO ii = 1, nloopmpi
            IF (ideniteratebal) THEN
c               IF (ii.EQ.numproc) GOTO 444
c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
c              IF (nptmasstot.GE.1) print *,iproc,' CIRCULAR ',ii
c               iahead = MOD(iproc+ii,numproc)
c               ibehind = MOD(numproc+iproc-ii,numproc)
               iahead = ipsend(ii)
               ibehind = iprec(ii)
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
     &                 '.GE.mmax2',ntot,inumber
                  CALL quit(1)
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
                  CALL quit(1)
               ENDIF
               IF (igotnumber.NE.inumber) THEN
                  WRITE (iprint, *) 'ERROR - igotnumber.NE.inumber:',
     &                 'iphase'
                  CALL quit(1)
               ENDIF
#ifdef MPIDEBUG
               IF (itime.EQ.2097152)
     &              print *,iproc,': Received iphase ',inumber
#endif

               IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) THEN
                  CALL MPI_SENDRECV(icolumnnext,1,indexMPI_INT,
     &                 iahead, 103,
     &                 icolumnnext(ntot+1), idim, MPI_INTEGER, 
     &                 ibehind, 103, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_INTEGER,
     &                 igotnumber,ierr)
                  IF (ntot+1+inumber.GE.idim2) THEN
                     WRITE (*,*) iproc,
     &                    ': ERROR - ntot+1+inumber.GE.idim (2)'
                     CALL quit(1)
                  ENDIF
                  IF (igotnumber.NE.inumber) THEN
                     WRITE (iprint, *) 
     &                    'ERROR - igotnumber.NE.inumber: colnext'
                     CALL quit(1)
                  ENDIF
               ENDIF

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
                     CALL quit(1)
                  ENDIF
               ENDIF

            ELSE
c
c--IMPLEMENTATION USING SEND
c
               IF (ii.EQ.isendproc+1) THEN
                  DO jjj = 1, nprocsend
                     jj = ipsend(jjj)
                     IF (iproc.EQ.jj) THEN
                        PRINT*,iproc,": ERROR - rank sending to itself"
                        CALL quit(0)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': Sending ',nlstsend,
     &                    ' xyzmh to ',jj,xyzmh(1,lsendlist(1)+1),
     &                    iphase(lsendlist(1)+1)
#endif

                     CALL MPI_SEND(xyzmh,1,indexMPI5,jj,0,
     &                    MPI_COMM_WORLD, ierr)

                     CALL MPI_SEND(iphase,1,indexMPI_INT1,jj,1,
     &                    MPI_COMM_WORLD, ierr)

                     IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) 
     &                    THEN
                        CALL MPI_SEND(icolumnnext,1,indexMPI_INT,
     &                       jj,103,MPI_COMM_WORLD, ierr)
                     ENDIF

                     IF (iteration.EQ.1) THEN
                        CALL MPI_SEND(iuniquesend,nptmasssend,
     &                       MPI_INTEGER8,jj,9,MPI_COMM_WORLD, ierr)
                     ENDIF
                  END DO
                  inumber = 0
                  ibehind = isendproc
               ELSE
                  ibehind = iprec(ii)
c
c--Receive data
c
                  CALL MPI_RECV(xyzmh(1,2*ntot+3),idim,i5REAL8,ibehind,
     &                 0,MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
                  IF (2*ntot+3+inumber.GE.mmax2) THEN
                WRITE (*,*) iproc,': ERROR - 2*ntot+3+inumber.GE.mmax2',
     &                    ntot,inumber,mmax2
                     CALL quit(1)
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
                     CALL quit(1)
                  ENDIF
                  IF (igotnumber.NE.inumber) THEN
                     WRITE (iprint, *) 'ERROR - igotnumber.NE.inumber:',
     &                    ' iphase'
                     CALL quit(1)
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': Received iphase ',inumber,
     &                 ' from ',ibehind,iphase(ntot+1)
#endif

                  IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) THEN
                     CALL MPI_RECV(icolumnnext(ntot+1),idim,
     &                    MPI_INTEGER, ibehind, 103, 
     &                    MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_INTEGER,igotnumber,
     &                    ierr)
                     IF (ntot+1+inumber.GE.idim2) THEN
                        WRITE (*,*) iproc,
     &                       ': ERROR - ntot+1+inumber.GE.idim2'
                        CALL quit(1)
                     ENDIF
                     IF (igotnumber.NE.inumber) THEN
                        WRITE (iprint, *) 
     &                       'ERROR - igotnumber.NE.inumber: colnext'
                        CALL quit(1)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': Received icolumnnext ',inumber,
     &                    ' from ',ibehind,iphase(ntot+1)
#endif
                  ENDIF

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
                        CALL quit(1)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF   ! ideniteratebal
c
c--This calculates the contributions to rhoi, gradhi, numneighreal, and 
c     gravity forces from non-neighbours
c
            DO i = 1, inumber
               ipart = i + ntot
               IF (ipart.GT.idim2) THEN
                  WRITE (*,*) 'ipart.GT.idim2'
                  CALL quit(1)
               ENDIF
               IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
                  IF (iteration.NE.1) THEN
                     WRITE (*,*) 'ERROR - iteration.NE.1'
                     CALL quit(1)
                  ENDIF
                  nptmasslocal = nptmasslocal + 1
                  listpm(nptmasslocal) = ipart
                  listrealpm(ipart) = nptmasslocal
                  IF (iuniquestore(nptmasslocal).NE.
     &                 iuniquestoreold(nptmasslocal)) nptmasstotlast=-1

c                  WRITE (*,*) iproc,': Trans ',ipart,
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
                  CALL quit(1)
               ENDIF
            ENDIF

            inumreturn = 0
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumber,ntot,npart,xyzmh,dvxyzu,inumreturn,llistrec)
C$OMP& shared(nlstdo_tot,itime,iproc,iphase)
C$OMP& private(i,ipart,nneighlocal)
            DO i = 1, inumber
               ipart = i + ntot
               IF (ipart.GT.idim2) THEN
                  WRITE (*,*) 'ipart.GT.idim2'
                  CALL quit(1)
               ENDIF

               CALL densitygradh(ipart,itime,npart,ntot,nlstdo_tot,
     &              xyzmh,nneighlocal,dvxyzu)

               IF (istellar_radiation) THEN
                  IF (iphase(ipart).EQ.0)
     &                 CALL column_los(ipart,itime,ntot,xyzmh)
               ENDIF
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
                  CALL quit(1)
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
                  CALL quit(1)
               ENDIF

               CALL MPI_SENDRECV(rho, 1, indexMPI1return, ibehind, 11,
     &              rhorec, inumberneigh, MPI_REAL4, iahead, 11,
     &              MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_REAL4,icheck,ierr)
               IF (icheck.GT.ineighproc) THEN
                  WRITE (*,*) 'ERROR - icheck.GT.ineighproc 1 ',icheck
                  CALL quit(1)
               ENDIF
               IF (icheck.NE.inumberneigh) THEN
                  WRITE (*,*) 'ERROR - icheck.NE.inumberneigh 1 ',
     &                 iproc
                  CALL quit(1)
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
                  CALL quit(1)
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
                  CALL quit(1)
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
c      particles on sinks if iptintree=1 and also sink particles if
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
                  CALL quit(1)
               ENDIF
               CALL MPI_SENDRECV(dvxyzu(1,ntot+1),4*inumber,
     &              MPI_REAL8,ibehind,15,dvxyzu(1,ntotplusinumber+1),
     &              4*nlstsend,MPI_REAL8, iahead, 15, 
     &              MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
               IF (ireturned.NE.4*nlstsend) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.4*nlstsend dvxyz'
                  CALL quit(1)
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
c--Receive back quantities for column densities on sightlines (if req)
c
               IF (idustRT.GT.0) THEN

                  CALL MPI_SENDRECV(sightcolumns(1,1,ntot+1),
     &                 nsightlinesmax*2*inumber,
     &                 MPI_REAL4,ibehind,16,
     &                 sightcolumns(1,1,ntotplusinumber+1),
     &                 nsightlinesmax*2*nlstsend,MPI_REAL4, iahead, 16,
     &                 MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL4,ireturned,ierr)
                  IF (ireturned.NE.nsightlinesmax*2*nlstsend) THEN
                     WRITE (*,*) 'ERROR - ireturned.NE.',
     &                    'nsightlinesmax*2*nlstsend'
                     CALL quit
                  ENDIF

#ifdef MPIDEBUG
                  print *,iproc,': transferred sightcolumns ',ibehind,
     &                 iahead
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstsend,lsendlist,ntotplusinumber,sightcolumns)
C$OMP& private(l,k,ipos,jpos)
                  DO l = 1, nlstsend
                     ipos = lsendlist(l) + 1
                     jpos = ntotplusinumber + l

                     DO k = 1, nsightlinesmax
                        sightcolumns(k,1,ipos) =sightcolumns(k,1,ipos) +
     &                       sightcolumns(k,1,jpos)
                        sightcolumns(k,2,ipos) =sightcolumns(k,2,ipos) +
     &                       sightcolumns(k,2,jpos)
                     END DO
                  END DO
C$OMP END PARALLEL DO
               ENDIF   ! idustRT>0
               IF (istellar_radiation) THEN
c
c--Transfer stellarrad contributions from other MPI processes
c     NOTE:  This transfers the data for iptdim sink particles all at
c     once (so iptdim should not be too much larger than the actual
c     number of sinks or the total amount of data transferred will be
c     much larger than necessary).
c
                  CALL MPI_SENDRECV(stellarrad(1,1,ntot+1),
     &                 4*iptdim*inumber,
     &                 MPI_REAL8,ibehind,201,
     &                 stellarrad(1,1,ntotplusinumber+1),
     &                 4*iptdim*nlstsend,MPI_REAL8, iahead, 201,
     &                 MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,ireturned,ierr)
                  IF (ireturned.NE.4*iptdim*nlstsend) THEN
                     WRITE (*,*) 'ERROR - ireturned.NE.',
     &                    '4*iptdim*nlstsend'
                     CALL quit
                  ENDIF
#ifdef MPIDEBUG
                  print *,iproc,': transferred stellarrad ',ibehind,
     &                 iahead
#endif
C$OMP PARALLEL default(none)
C$OMP& shared(nlstsend,lsendlist,ntotplusinumber,nptmasstot,stellarrad)
C$OMP& private(l,k,ipos,jpos)
C$OMP DO SCHEDULE(runtime)
                  DO l = 1, nlstsend
                     ipos = lsendlist(l) + 1
                     jpos = ntotplusinumber + l

                     DO k = 1, nptmasstot
                        stellarrad(1:4,k,ipos) =
     &                       stellarrad(1:4,k,ipos) +
     &                       stellarrad(1:4,k,jpos)
                     END DO
                  END DO
C$OMP END DO
C$OMP END PARALLEL
               ENDIF   ! istellar_radiation
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
               IF (ii.EQ.isendproc+1) THEN
                  DO jjj = 1, nprocsend
                     jj = ipsend(jjj)
                     CALL MPI_RECV(llistrec, nlstsend, 
     &                    MPI_INTEGER, MPI_ANY_SOURCE, 10, 
     &                    MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_INTEGER,
     &                    inumberneigh,ierr)
                     iahead = istatus(MPI_SOURCE)
                     IF (inumberneigh.GT.ineighproc) THEN
                        WRITE (*,*) 'ERROR - ',
     &                       'inumberneigh.GT.ineighp'
                        CALL quit(1)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got indices ',inumberneigh,
     &                    ' from ',iahead,numproc
#endif
                     CALL MPI_RECV(rhorec,inumberneigh,MPI_REAL4,
     &                    iahead,11,MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL4,
     &                    icheck,ierr)
                     IF (icheck.GT.ineighproc) THEN
                        WRITE (*,*) 'ERROR - icheck.GT.',
     &                       'ineighproc 1 ',icheck
                        CALL quit(1)
                     ENDIF
                     IF (icheck.NE.inumberneigh) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.',
     &                       'inumberneigh 1 ',iproc
                        CALL quit(1)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got rho ',icheck
#endif

                     CALL MPI_RECV(dumrhorec,inumberneigh,MPI_REAL4,
     &                    iahead,12,MPI_COMM_WORLD,istatus,ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL4,
     &                    icheck,ierr)
                     IF (icheck.NE.inumberneigh) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.',
     &                       'inumberneigh 1 ',iproc
                        CALL quit(1)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got gradh ',icheck
#endif

                     CALL MPI_RECV(nneighrec,inumberneigh, 
     &                    MPI_INTEGER,iahead,13,
     &                    MPI_COMM_WORLD,istatus,ierr)
#ifdef MPIDEBUG
                     print *,iproc,': rec '
#endif
                     CALL MPI_GET_COUNT(istatus,MPI_INTEGER,
     &                    icheck,ierr)
#ifdef MPIDEBUG
                     print *,iproc,': count '
#endif
                     IF (icheck.NE.inumberneigh) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.',
     &                       'inumberneigh 1 ',iproc
                        CALL quit(1)
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
     &                       .NOT.nneighsentany(iahead+1)) THEN
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
     &                          iahead+1) = ipos - 1
                        ENDIF
                     END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
                     print *,iproc,': received neighbour data'
#endif
                     IF (iteration.EQ.1) THEN
#ifdef MPIDEBUG
                        print *,iproc,': reports ',inumberneigh,
     &                       ' foreign neighbours on process ',
     &                       iahead
#endif
                        nneighsentbackold(iahead+1) = 0
                        nneighsentback(iahead+1) = inumberneigh
                     ELSE
                        nneighsentbackold(iahead+1) = 
     &                       nneighsentback(iahead+1)
                        nneighsentback(iahead+1) = 
     &                       nneighsentback(iahead+1) +
     &                       inumberneigh
                     ENDIF
                     maxnneighsentback = MAX(maxnneighsentback,
     &                    nneighsentback(iahead+1))
c
c--Receive gravity forces back (includes forces from gas
c      particles on sinks if iptintree=1 and also sink particles if
c       iptintree=2).
c
                     IF (ntotplusinumber+nlstsend.GE.idim3) THEN
                        WRITE (*,*) iproc,
     &                  ': ERROR - ntotplusinumber+nlstsend.GE.idim ',
     &                  ntotplusinumber,nlstsend,idim,npart,ntot
                        CALL quit(1)
                     ENDIF
                     CALL MPI_RECV(dvxyzu(1,ntotplusinumber+1),
     &                    4*nlstsend, MPI_REAL8, iahead, 15, 
     &                    MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, MPI_REAL8,
     &                    ireturned, ierr)
                     IF (ireturned.NE.4*nlstsend) THEN
                        WRITE (*,*) 'ERROR - ireturned.NE.',
     &                       '4*nlstsend dvxyz'
                        CALL quit(1)
                     ENDIF

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstsend,lsendlist,ntotplusinumber,dvxyzu)
C$OMP& private(l,ipos,jpos)
                     DO l = 1, nlstsend
                        ipos = lsendlist(l) + 1
                        jpos = ntotplusinumber + l
                        
                        DO k = 1, 4
                           dvxyzu(k,ipos) = dvxyzu(k,ipos) + 
     &                          dvxyzu(k,jpos)
                        END DO
                     END DO
C$OMP END PARALLEL DO

c
c--Receive back quantities for column densities on sightlines (if req)
c
                     IF (idustRT.GT.0) THEN

                        CALL MPI_RECV(sightcolumns(1,1,
     &                     ntotplusinumber+1),
     &                     nsightlinesmax*2*nlstsend, MPI_REAL4, iahead, 
     &                     16, MPI_COMM_WORLD,istatus,ierr)
                        CALL MPI_GET_COUNT(istatus,MPI_REAL4,
     &                     ireturned,ierr)
                        IF (ireturned.NE.nsightlinesmax*2*nlstsend) THEN
                           WRITE (*,*) 'ERROR - ireturned.NE.',
     &                          'nsightlinesmax*2*nlstsend'
                           CALL quit
                        ENDIF
#ifdef MPIDEBUG
                        print *,iproc,': got sightcolumns from ',iahead
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstsend,lsendlist,ntotplusinumber,sightcolumns)
C$OMP& private(l,ipos,jpos)
                        DO l = 1, nlstsend
                           ipos = lsendlist(l) + 1
                           jpos = ntotplusinumber + l

                           DO k = 1, nsightlinesmax
                              sightcolumns(k,1,ipos) = 
     &                             sightcolumns(k,1,ipos) +
     &                             sightcolumns(k,1,jpos)
                              sightcolumns(k,2,ipos) = 
     &                             sightcolumns(k,2,ipos) +
     &                             sightcolumns(k,2,jpos)
                           END DO
                        END DO
C$OMP END PARALLEL DO
                        ENDIF

                  END DO
               ELSE
c
c--Send data
c
                  ibehind = iprec(ii)

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

                  IF (idustRT.GT.0) THEN
                     CALL MPI_SEND(sightcolumns(1,1,ntot+1),
     &                    nsightlinesmax*2*inumber,
     &                    MPI_REAL4,ibehind,16,
     &                    MPI_COMM_WORLD, istatus, ierr)
                  ENDIF
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
         CALL MPI_TYPE_FREE(indexMPI_INT,ierr)
         CALL MPI_TYPE_FREE(indexMPI5,ierr)
c
c--If individualtimesteps=2, then not all sink particles may have been
c     transferred because some might not be being evolved this timestep.
c
         IF (individualtimesteps.NE.2) THEN 
            IF (iteration.EQ.1) THEN
               IF (nptmasslocal2.NE.nptmasstot) THEN
                  WRITE (*,*) 'ERROR - nptmasslocal2.NE.nptmasstot ',
     &                 nptmasslocal2, nptmasstot
                  CALL quit(1)
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
C$OMP& shared(nneigh,selfnormkernel,hfact,distancemax)
C$OMP& shared(cnormk,radkernel,dvtable,ddvtable,wij,grwij)
C$OMP& shared(listpm,iphase,dphidh,uradconst,icall,encal)
C$OMP& shared(iprint,nptmass,iptmass,iorig,third)
C$OMP& shared(dumrho,iscurrent,ibound,ireal,hi_2back)
C$OMP& shared(isteps,it0,it1,imax,imaxstep,dt,itime)
C$OMP& shared(varmhd,Bevol,Bxyz,nlstdo,ldolist)
C$OMP& shared(dvxyzu,hi_old,isucceed,poten,isend,numproc,iproc)
C$OMP& shared(ivar,ijvar,ncompact,icompact,iteration,dust_tk)
C$OMP& shared(hminbisec,hmaxbisec,iupdated,nlistupdated,listparents)
C$OMP& shared(heatingISR,sightcolumns,nsightlines,umass,udist)
C$OMP& shared(chemistry,icolumnnext,icolumnsteps,bounddens)
C$OMP& shared(dustvar)
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend,ijvartemp,listparentstemp)
#endif
C$OMP& private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi)
C$OMP& private(divvi,curlvxi,curlvyi,curlvzi,gradhi,gradsofti)
C$OMP& private(pmassj,hi_oldi,hi1,hi21,hi31,hi41,hneigh)
C$OMP& private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2,rcut)
C$OMP& private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij)
C$OMP& private(projv,procurlvx,procurlvy,procurlvz)
C$OMP& private(l,iptcur,dphi,dwdhi,dpotdh,numneighi)
C$OMP& private(numneighreal,rhohi,dhdrhoi,omegai,func,dfdh1)
C$OMP& private(hnew,deltarho,heatingISRival,dust_gas)
C$OMP& private(nneighlocal,nneighi,ncompacthere,icompacthere)
C$OMP& private(numberupdated,dust_kappai,photo,dust_cooling)
C$OMP& private(cooling,gas_dust,cosmic_ray,func_gas,derivative)
C$OMP& private(heatingISRi,Gphotoelectrici,ext,exp_Av,itemp)
C$OMP& private(u_found,cv_val,gas_temp,ijk,frac_change,h2column)
C$OMP& private(columnnumberdensity,expext,selfval,selffac,selfshield)
C$OMP& private(exp_selfval,Qabs550)
C$OMP& private(iconvergence1,iconvergence2,iconvergence3)
C$OMP& reduction(MAX:rhonext,imaxit)
C$OMP& reduction(+:inumit,inumfixed,inumrecalc,nwarnup,nwarndown)
C$OMP& reduction(+:nbisection,numsucceeded)
#ifdef MPICOPY
      DO n = numberstart, numberend
#else
      DO n = 1, nlstdo
#endif
         ipart = ldolist(n)
c
c--Adds on local contributions to rho, gradh, and nneigh and long-range
c       gravity
c
         CALL densitygradh(ipart,itime,npart,ntot,nlstdo_tot,xyzmh,
     &        nneighlocal,dvxyzu)
c
c--Non-sink particles
c
         IF (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10) THEN

            IF (nneighlocal.GT.nneighmax) THEN
               WRITE (iprint,*) 'ERROR: nneighlocal exceeds nneighmax'
               WRITE (iprint,*) iorig(ipart),nneighlocal,nneighmax
               CALL quit(1)
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
               write(iprint,99100) iteration,hnew,ipart,nneighi,
     &              hi,func,hmaxbisec(ipart),hminbisec(ipart)
99100          FORMAT ('bisection (',I4,'): hnew =',1PE12.5,1X,I9,1X,I5,
     &              4(1X,1PE12.5))
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
c
c--Output data if a single particle is failing to converge
c
c            IF (nlstdo.EQ.1) THEN
c               itemp = ldolist(nlstdo)
c               WRITE (66,77666) itemp,xyzmh(1,itemp),xyzmh(2,itemp),
c     &              hnew,omegai,isucceed(itemp),iteration,
c     &              maxitsnr,nneighi,hi,hi_oldi,hi_2back(itemp)
c77666          FORMAT(I8,1X,4(1PE12.5,1X),4(I5,1X),3(1PE12.5,1X))
c            ENDIF

            IF ((hnew.LE.0. .OR. omegai.LE.tiny .OR. 
     &           iteration.EQ.maxitsnr .OR. nneighi.LE.0) 
     &           .AND. isucceed(ipart).NE.-1) THEN
c
c--Switch to Bisection if not converging or running into trouble
c
               WRITE(iprint,*) 'WARNING: switching to bisection on'//
     &              ' particle ',iorig(ipart),ipart,'(',ipart,') hi = ',
     &              hi,' hnew = ',hnew
c               WRITE(*,*) 'WARNING: switching to bisection on'//
c     &              ' particle ',iorig(ipart),'(',ipart,') hi = ',hi,
c     &              ' hnew = ',hnew,nneighi,gradhi,rhoi,rhohi,dhdrhoi,
c     &              omegai,iteration
               IF (nneighi.LE.0) THEN
                  WRITE(iprint,*) '(particle has no neighbours) ',
     &                 nneighi
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
               hmaxbisec(ipart) = 2.0*distancemax
!--Don't have to start with ridiculous h, 
!  just something reasonably big between above limits
               hnew = 2.*hi_oldi
!               hnew = 0.5*(hminbisec(ipart) + hmaxbisec(ipart))
c
c--Otherwise check for convergence
c
            ELSE
               iconvergence1 = (ABS((hnew-hi)/hi_oldi).LT.htol
     &              .AND.isucceed(ipart).EQ.0
     &              .AND. ABS((hnew-hi)/hi).LT.htol)
               iconvergence2 = ABS((hnew-hi)/hi_oldi).LT.1.0E-6
     &              .AND. ABS((hnew-hi)/hi).LT.htol
c
c--Type 3 convergence was added because very ocassionally, a single 
c     particle would get stuck oscillating between two solutions,
c     for example, between solutions with 63 or 64 neighbours where
c     h differed by slightly more than the required tolerance.
c
               iconvergence3 = (iteration.GT.maxitsnr .AND.
     &              isucceed(ipart).EQ.0 .AND.
     &              ABS((hi_2back(ipart)-hnew)/hi_oldi).LT.1.0E-8)

               IF (iconvergence1 .OR. iconvergence2 .OR.
     &              iconvergence3
c     &              .OR. (iteration.GE.maxiterations-1)
     &              ) THEN

                  imaxit = MAX(imaxit, iteration)
                  inumit = inumit + iteration
                  IF (isucceed(ipart).EQ.-1) THEN
                     write(iprint,*) 'SWITCH BACK TO N-R ',ipart,hnew,hi
c                     write(*,*) 'SWITCH BACK TO N-R ',ipart,hnew,hi
                     isucceed(ipart) = 0 ! switch back to N-R
                  ELSE
                     IF (iconvergence3 .AND. .NOT. iconvergence1) THEN
                        WRITE (iprint,11003) ipart,hnew,hi,
     &                       hi_2back(ipart),hi_oldi,xyzmh(ipart,1),
     &                       xyzmh(ipart,2),xyzmh(ipart,3)
                        WRITE (*,11003) ipart,hnew,hi,
     &                       hi_2back(ipart),hi_oldi,xyzmh(ipart,1),
     &                       xyzmh(ipart,2),xyzmh(ipart,3)
11003                   FORMAT ('N-R got stuck oscillating ',I9,
     &                       7(1X,1PE12.5))
                     ENDIF
c
c--Store quantities if converged
c
                     isucceed(ipart) = 1
                     numsucceeded = numsucceeded + 1

                     rho(ipart) = rhoi
                     dumrho(ipart) = rhoi

                     xyzmh(5,ipart) = hi
                     gradhs(1,ipart) = 1./omegai
c
c--The 4th element is used to compute the potential energy while in
c     densityiterate
c
                     poten(ipart) = dvxyzu(4,ipart)

                     IF (iphase(ipart).EQ.0) THEN
                        IF (iradtrans.EQ.idim) THEN
c
c--Calculating the ISR attenuation factors for a non-MPI calculation can
c     be done in column_density.F.  However, it needs to be done here for
c     the MPI code (see comments in column_density.F).
c
#ifdef MPI
#ifdef RT
                           IF (idustRT.GT.0) THEN
                              IF (ioptimise_column.EQ.1) THEN

                              IF (itime.LT.icolumnnext(ipart)) GOTO 500
                              IF (itime.GT.icolumnnext(ipart)) THEN
                                 WRITE (*,*) 
     &                           'ERROR - itime.GT.icolumnnext(ipart) ',
     &                                itime,icolumnnext(ipart),ipart
                              ENDIF
                              ENDIF

                              heatingISRi = 0.
                              Gphotoelectrici = 0.
                              exp_Av = 0.
                              selfshield = 0.
                              DO k = 1, nsightlines
c
c--Based on Zucconi et al. (2001)
c
                                 columnnumberdensity = 
     &                                sightcolumns(k,1,ipart)*
     &                                umass/udist**2 /(gmw*1.67E-24)
     &                                + bounddens
                                 h2column = sightcolumns(k,2,ipart)*
     &                                umass/udist**2 /(gmw*1.67E-24)
     &                                + bounddens

                                 ext = columnnumberdensity * 
     &                                Qv(c/0.0000550,Qabs550)

                                 heatingISRi = heatingISRi + 
     &                                get_heatingISR(ext)
                                 Gphotoelectrici = Gphotoelectrici +
     &                                get_Gphotoelectric(ext)
                                 expext = EXP(-ext)
                                 exp_Av = exp_Av + expext

                                 selfval = SQRT(1.0+h2column/(5.0E+14))
                                 exp_selfval = EXP(MAX(-8.5E-4*selfval,
     &                                -60.0))
                                 selffac = 
     &                                0.965/(1.0+h2column/(5.0E+14))**2+
     &                                0.035*exp_selfval/selfval
                                 IF (expext.LT.1.0E-8) expext = 1.0E-8
                                 selfshield = selfshield + 
     &                                expext**3.74 * selffac
                              END DO
                              heatingISR(1,ipart) = 
     &                             heatingISRi/nsightlines
                              heatingISR(2,ipart) = 
     &                             Gphotoelectrici/nsightlines
                              heatingISR(3,ipart) = exp_Av/nsightlines
                              heatingISR(4,ipart) = 
     &                             selfshield/nsightlines
                           ENDIF
#endif
#endif
c
c--Set e(i) for the first time around because density only set here
c     This sets radiation and matter to have same initial temperature
c
#ifdef MPI
#ifdef RT
 500                       CONTINUE
#endif
#endif
                           IF (icall.EQ.1 .AND. encal.EQ.'r') THEN
                              IF (ekcle(1,ipart).EQ.0.0) THEN
#ifdef RT
                                 IF (idustRT.GT.0) THEN
c
c--Set temperature of the radiation to be 1 K initially (i.e. negligible)
c
                                    ekcle(1,ipart) =uradconst/rho(ipart)
c
c--If starting from a dump file without radiative transfer information,
c     need to define a starting value of U/T, otherwise the call to 
c     dust_temperature.f below will fail.
c
                                    IF (ekcle(3,ipart).LE.0.0) THEN
                                       ekcle(3,ipart) = 
     &                                 getcv(rho(ipart),vxyzu(4,ipart))
                                    ENDIF

                                    gas_temp = 
     &                                   gas_temperature_equilibrium(
     &                          ipart,dust_temperature(ipart,ntot,
     &                          uradconst,0.0,vxyzu(4,ipart),ekcle,
     &                          rhoi,dust_kappai,dust_cooling,
     &                          heatingISRival,dust_gas),rhoi,photo,
     &                          cooling,gas_dust,cosmic_ray,func_gas,
     &                          derivative)

                                    vxyzu(4,ipart) = 
     &                                   getu(rho(ipart),gas_temp)
                                    ekcle(3,ipart) = 
     &                                getcv(rho(ipart),vxyzu(4,ipart))

                                    dust_tk(1,ipart) = 
     &                                   dust_temperature(ipart,ntot,
     &                          uradconst,0.0,vxyzu(4,ipart),ekcle,
     &                          rhoi,dust_kappai,dust_cooling,
     &                          heatingISRival,dust_gas)
                                    dust_tk(2,ipart) = dust_kappai

c                        IF (iproc.EQ.0) THEN
cC$OMP CRITICAL (write19)
ccc                           IF (rhoi.GT.0.05 .AND. rhoi.LT.0.08)
c                       WRITE (19,445) SQRT(xyzmh(1,ipart)**2+
c     &                       xyzmh(2,ipart)**2+xyzmh(3,ipart)**2),
c     &                       vxyzu(4,ipart)/ekcle(3,ipart),
c     &                       dust_temperature(ipart,ntot,uradconst,
c     &                       0.0,vxyzu(4,ipart),ekcle,
c     &                       rhoi,dust_kappai,dust_cooling,
c     &                       heatingISRival,dust_gas),
c     &                       dust_kappai,rhoi,photo,
c     &                       dust_cooling,heatingISRival,dust_gas,
c     &                       cooling,gas_dust,cosmic_ray,func_gas,
c     &                       derivative,gas_temp,vxyzu(4,ipart),
c     &                       ekcle(3,ipart),chemistry(1,ipart),
c     &                       chemistry(2,ipart),chemistry(3,ipart),
c     &                       heatingISR(3,ipart),ipart
cC$OMP END CRITICAL (write19)
c 445                   FORMAT(21(1PE12.5,1X),I8)
c                        ENDIF
                                 ELSE
#endif
c
c--Set specific radiation energy to be in thermal equilibrium with the
c     matter initially
c
                                    ekcle(3,ipart) =
     &                                getcv(rho(ipart),vxyzu(4,ipart))
                                    ekcle(1,ipart) = 
     &                                   uradconst*(vxyzu(4,ipart)/
     &                                   ekcle(3,ipart))**4/rho(ipart)
#ifdef RT
                                 ENDIF
#endif
c
c--Set opacities
c
                                 ekcle(2,ipart) = 
     &                                getkappa(vxyzu(4,ipart),
     &                                ekcle(3,ipart),rho(ipart))
                              ENDIF
                           ENDIF
                        ENDIF
c
c--Pressure and sound velocity from ideal gas law...
c
                        CALL eospg(ipart, vxyzu, rho, pr, vsound, ekcle,
     &                       dustvar)
                     ENDIF
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
                        WRITE (*,*) 
     &                       'ERROR - compact not large enough DI ',
     &                       icompacthere, nneighlocal, icompactmax
                        CALL quit(1)
                     ENDIF
                     ivar(1,ncompacthere) = nneighlocal
                     ivar(2,ncompacthere) = icompacthere
                     ivar(3,ncompacthere) = ipart

                     DO k = 1, nneighlocal
                        j = neighlist(k)
#ifdef MPICOPY
c
c--For MPICOPY need to store neighbours in temp array that is combined
c     into ijvar() across all processes later
c
                        ijvartemp(icompacthere + k) = j
#else
                        ijvar(icompacthere + k) = j
#endif
                        IF (.NOT.iscurrent(j) .AND. .NOT.iupdated(j) 
     &                       .AND. j.LE.npart) THEN
                           iupdated(j) = .TRUE.
C$OMP FLUSH(iupdated)
C$OMP CRITICAL (listupdated)
                           nlistupdated = nlistupdated + 1
                           numberupdated = nlistupdated
C$OMP END CRITICAL (listupdated)
                           IF (numberupdated.GT.idim) THEN
                              WRITE (iprint,*) 
     &                             'ERROR - numberupdated.GT.idim'
                              CALL quit(1)
                           ENDIF
#ifdef MPICOPY
c
c--For MPICOPY need to store parents list in temp array that is combined
c     across all processes later
c
                           listparentstemp(numberupdated) = j
#else
                           listparents(numberupdated) = j
#endif
                        ENDIF
                     END DO

                  ENDIF
               ENDIF
            ENDIF

         ELSE
c
c--Else for sink particles
c
         isucceed(ipart) = 1
         numsucceeded = numsucceeded + 1
c
c--The 4th element is used to compute the potential energy while in
c     densityiterate
c
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
         IF (isucceed(ipart).NE.1) THEN
            hi_2back(ipart) = hi
            xyzmh(5,ipart) = hnew
         ENDIF

      END DO
C$OMP END PARALLEL DO
c
c--Check to see whether all particles have converged
c
#ifdef MPIDEBUG
#ifdef MPI
      print *,'numsucceeded ',numsucceeded,nlstdo,htol,numproc
#else
      print *,'numsucc ',numsucceeded,htol,iproc,numbertodohere
#endif
#endif

#ifdef MPICOPY
      IF (numsucceeded.EQ.numbertodohere) THEN
#else
      IF (numsucceeded.EQ.nlstdo) THEN
#endif
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
#ifdef MPICOPY
      numberend_old = numberend
      numberend = numberstart - 1
      DO n = numberstart, numberend_old
         ipart = ldolist(n)
         IF (isucceed(ipart).NE.1) THEN
            numberend = numberend + 1
            ldolist(numberend) = ipart
         ENDIF
      END DO
      numbertodohere = numberend - numberstart + 1
#else
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
#endif
c
c--Re-zero quantities
c
c--Doesn't zero quantities if dens_iterate has failed, allowing printing.
c
      IF (iteration.NE.maxiterations) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstdo,ldolist,rho,dumrho,divv,curlv,gradhs,nneigh)
C$OMP& shared(dvxyzu,poten,dust_tk,stellarrad,nptmasstot)
#ifdef MPI
C$OMP& shared(sightcolumns)
#endif
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend)
#endif
C$OMP& private(n,ipart,k)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, nlstdo
#endif
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
            IF (idustRT.GT.0) THEN
               DO k = 1, 2
                  dust_tk(k,ipart) = 0.
               END DO
#ifdef MPI
               DO k = 1, nsightlinesmax
                  sightcolumns(k,1,ipart) = 0.
                  sightcolumns(k,2,ipart) = 0.
               END DO
#endif
            ENDIF
            IF (istellar_radiation) THEN
               stellarrad(1:4,1:nptmasstot,ipart) = 0.
            ENDIF
            poten(ipart) = 0.
         END DO
C$OMP END PARALLEL DO
      ENDIF
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
#ifdef MPICOPY
      WRITE (iprint,*) 'ERROR: density iteration failed for ',
     &     numberend-numberstart,' particles ',ifinished,iproc
      WRITE (*,*) 'ERROR: density iteration failed for ',
     &     numberend-numberstart,' particles ',ifinished,iproc
#else
#ifdef MPI
      IF (nlstdo.GT.0) THEN
#endif
      WRITE (iprint,*) 'ERROR: density iteration failed for ',nlstdo,
     &     ' particles ',ifinished,iproc
c
c--Added to give more information about failed particle. (fort.65)
c
      itemp = ldolist(nlstdo)

11001 FORMAT (I8,X,I4,X,I8,6(X,1PE15.8),X,L1)
      WRITE (65,11001) itemp, iphase(itemp), nneigh(itemp), rho(itemp),
     &     xyzmh(1,itemp), xyzmh(2,itemp), xyzmh(3,itemp),
     &     xyzmh(5,itemp), 0, iscurrent(itemp)
      WRITE (65,*)
      DO k = 1, nneigh(itemp)
         iitemp = neighlist(k)
         rtemp = SQRT((xyzmh(1,itemp)-xyzmh(1,iitemp))**2 + 
     &        (xyzmh(2,itemp)-xyzmh(2,iitemp))**2 + 
     &        (xyzmh(3,itemp)-xyzmh(3,iitemp))**2)
         WRITE (65,11001) iitemp, iphase(iitemp), nneigh(iitemp),
     &        rho(iitemp), xyzmh(1,iitemp), xyzmh(2,iitemp),
     &        xyzmh(3,iitemp), xyzmh(5,iitemp), rtemp,
     &        iscurrent(iitemp)
      ENDDO

      WRITE (*,*) 'ERROR: density iteration failed for ',nlstdo,
     &     ' particles ',ifinished,iproc
#ifdef MPI
      ENDIF
#endif
#endif
#ifdef MPIALL
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      CALL quit(1)

 30   CONTINUE

c
c--Now for MPICOPY parallelisation, need to share all the neighbours lists
c     and gravitational forces, gradhs, etc between MPI processes.
c     Need to share: rho(),dumrho(),nneigh(),
c     dumvxyzu(1-3),xyzmh(5),gradhs(1),poten()
c
c     Need to call eospg for shared particles (not done locally)
c
c     Need to set ivar(1-3), store neighbours in ijvar() and need
c     to get listparents()
c
c     For RT:  Need to share start values of ekcle(1-3) if undefined
c
#ifdef MPICOPY
      numberend = numberendkeep
c
c--Send around rho, dumrho, and poten
c
      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
     &     lblocklengths, llisttrans(numberstart),
     &     MPI_REAL4, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their rho data
c
      CALL MPI_ALLGATHERV(rho,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,MPI_REAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put rho data into correct places in rho
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,real4transfer2)
C$OMP&shared(rho,numproc,iproc)
C$OMP&private(i,j)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1) 
               rho(llisttrans(i) + 1) = real4transfer2(i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
      print *,iproc,': rho move ',1,numbertodo*iproc,
     &     numbertodo*(iproc+1)+1, nlstdo_tot
#endif
c
c--All processes transfer their dumrho data
c
      CALL MPI_ALLGATHERV(dumrho,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,MPI_REAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put rho data into correct places in rho
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,real4transfer2)
C$OMP&shared(dumrho,numproc,iproc)
C$OMP&private(i,j)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               dumrho(llisttrans(i) + 1) = real4transfer2(i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
      print *,iproc,': dumrho move ',1,numbertodo*iproc,
     &     numbertodo*(iproc+1)+1, nlstdo_tot
#endif
c
c--All processes transfer their poten data
c
c     NOTE: ONLY NEEDED FOR TESTING densityiterate
c
cTEST      CALL MPI_ALLGATHERV(poten,1,indexMPI_INT1,
cTEST     &     real4transfer2,irecvcounti,idisplacementsi,MPI_REAL4,
cTEST     &     MPI_COMM_WORLD,ierr)
c
c--Put rho data into correct places in rho
c
cTEST      DO j = 0, numproc - 1
cTEST         IF (j.NE.iproc) THEN
cTEST            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
cTEST     &           irecvcounti(j+1)
cTEST               poten(llisttrans(i) + 1) = real4transfer2(i)
cTEST            END DO
cTEST         ENDIF
cTEST      END DO
#ifdef MPIDEBUG
cTEST      print *,iproc,': poten move ',1,numbertodo*iproc,
cTEST     &     numbertodo*(iproc+1)+1, nlstdo_tot
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
c
c--Send around nneigh
c
      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
     &     lblocklengths, llisttrans(numberstart),
     &     MPI_INTEGER, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their nneigh data
c
      CALL MPI_ALLGATHERV(nneigh,1,indexMPI_INT1,
     &     itransfer3,irecvcounti,idisplacementsi,MPI_INTEGER,
     &     MPI_COMM_WORLD,ierr)
c
c--Put nneigh data into correct places in nneigh
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,itransfer3)
C$OMP&shared(nneigh,numproc,iproc)
C$OMP&private(i,j)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               nneigh(llisttrans(i) + 1) = itransfer3(i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': nneigh move ',1,numbertodo*iproc,
     &     numbertodo*(iproc+1)+1, nlstdo_tot
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
c
c--Send around xyzmh(5)
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(numberstart,numbertodoherekeep,llisttrans)
C$OMP&private(i)
      DO i = numberstart, numberstart + numbertodoherekeep - 1
         llisttrans(i) = llisttrans(i)*5 + 4
      END DO
C$OMP END PARALLEL DO
      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
     &     lblocklengths, llisttrans(numberstart),
     &     MPI_REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(numberstart,numbertodoherekeep,llisttrans)
C$OMP&private(i)
      DO i = numberstart, numberstart + numbertodoherekeep - 1
         llisttrans(i) = llisttrans(i)/5
      END DO
C$OMP END PARALLEL DO
c
c--All processes transfer their xyzmh(5) data
c
      CALL MPI_ALLGATHERV(xyzmh,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,MPI_REAL8,
     &     MPI_COMM_WORLD,ierr)
c
c--Put xyzmh(5) data into correct places in xyzmh(5)
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans,realtransfer5to15)
C$OMP&shared(xyzmh,numproc,iproc)
C$OMP&private(i,j)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               xyzmh(5,llisttrans(i) + 1) = realtransfer5to15(i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO

#ifdef MPIDEBUG
      print *,iproc,': xyzmh(5) move ',1,numbertodo*iproc,
     &     numbertodo*(iproc+1)+1, nlstdo_tot
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
c
c--Send around dvxyzu(1-4)
c
c     NOTE: ONLY NEEDED FOR TESTING densityiterate
c
cTEST      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
cTEST      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)

cTEST      CALL MPI_TYPE_INDEXED(numbertodoherekeep, 
cTEST     &     lblocklengths, llisttrans(numberstart),
cTEST     &     i4REAL8, indexMPI_INT1, ierr)
cTEST      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their dvxyzu(1-4) data
c
cTEST      CALL MPI_ALLGATHERV(dvxyzu,1,indexMPI_INT1,
cTEST     &     realtransfer5to15,irecvcounti,idisplacementsi,i4REAL8,
cTEST     &     MPI_COMM_WORLD,ierr)
c
c--Put dvxyzu(1-4) data into correct places in dvxyzu(1-4)
c
cTEST      DO j = 0, numproc - 1
cTEST         IF (j.NE.iproc) THEN
cTEST            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
cTEST     &           irecvcounti(j+1)
cTEST               DO k = 1, 4
cTEST                  dvxyzu(k,llisttrans(i) + 1) = 
cTEST     &                 realtransfer5to15((i-1)*4 + k)
cTEST               END DO
cTEST            END DO
cTEST         ENDIF
cTEST      END DO

#ifdef MPIDEBUG
cTEST      print *,iproc,': dvxyzu(1-4) move ',1,numbertodo*iproc,
cTEST     &     numbertodo*(iproc+1)+1, nlstdo_tot
#endif

cTEST      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
cTEST      CALL MPI_TYPE_FREE(i4REAL8,ierr)
c
c--Send around ivar(1-3) values
c
      DO i = 1, numproc
         irecvcounti(i) = irecvcounti(i)*3
         idisplacementsi(i) = idisplacementsi(i)*3

#ifdef MPIDEBUG
         print *,iproc,': mod irecv ',i,irecvcounti(i),
     &        idisplacementsi(i)
#endif
      END DO

#ifdef MPIDEBUG
      print *,iproc,': check ncompact ',ncompact,numbertodo*numproc
      print *,iproc,': numbertodokeephere ',numbertodoherekeep
#endif

      CALL MPI_ALLGATHERV(ivar(1,numberstart),numbertodoherekeep*3,
     &     MPI_INTEGER,itransfer3,irecvcounti,idisplacementsi,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
c
c--Put ivar(1-3) data into correct places in ivar(1-3)
c
      DO i = 1, numproc
         irecvcounti(i) = irecvcounti(i)/3
         idisplacementsi(i) = idisplacementsi(i)/3
      END DO

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,itransfer3)
C$OMP&shared(ivar,numproc,iproc)
C$OMP&private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 3
                  ivar(k,i) = itransfer3((i-1)*3 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO

      ncompact = nlstdo_tot

#ifdef MPIDEBUG
      print *,iproc,': ivar(1-3) move ',1,numbertodo*iproc,
     &     numbertodo*(iproc+1)+1, nlstdo_tot, ncompact
      print *,iproc,': icompact ',icompact
#endif
c
c--Transfer neighbours arrary ijvar.  Need to use ALLGATHERV since
c     each group of particles will have a different number of neighbours.
c
c--First need each process to determine how many neighbours it is expecting
c     to receive from other nodes
c
c
      irecvcountn(iproc) = icompact
      CALL MPI_ALLGATHER(irecvcountn(iproc),1,
     &     MPI_INTEGER,irecvcountn,1,MPI_INTEGER,
     &     MPI_COMM_WORLD,ierr)

#ifdef MPIDEBUG
      print *,iproc,': irecvcountn ',(irecvcountn(i),i=1,numproc)
#endif

      idisplacementsn(1) = 0
      DO i = 1, numproc-1
         idisplacementsn(i+1) = idisplacementsn(i) + irecvcountn(i)
      END DO

      CALL MPI_ALLGATHERV(ijvartemp,icompact,MPI_INTEGER,
     &        ijvar,irecvcountn,idisplacementsn,MPI_INTEGER,
     &        MPI_COMM_WORLD,ierr)

#ifdef MPIDEBUG
      print *,' recvcount ',(irecvcountn(i),i=1,numproc)
      print *,' idisplacements ',(idisplacementsn(i),i=1,numproc)
#endif
c
c--Now need to update ivar(2) to point to the correct neighbours by adding
c     displacements
c
      DO i = 1, numproc
         DO j = 1, irecvcounti(i)
            jlocal = idisplacementsi(i) + j
            ivar(2,jlocal) = ivar(2,jlocal) + idisplacementsn(i)
         END DO
      END DO
c
c--Now need to transfer the list of parents of updated particles
c     Need to set iupdated=.TRUE. for all in list as well.
c
c
c--First need each process to determine how many parents it is expecting
c     to receive from other nodes
c
c
      irecvcountn(iproc) = nlistupdated
      CALL MPI_ALLGATHER(irecvcountn(iproc),1,
     &     MPI_INTEGER,irecvcountn,1,MPI_INTEGER,
     &     MPI_COMM_WORLD,ierr)

#ifdef MPIDEBUG
      print *,iproc,': irecvcountn updated',
     &     (irecvcountn(i),i=1,numproc)
#endif

      idisplacementsn(1) = 0
      DO i = 1, numproc-1
         idisplacementsn(i+1) = idisplacementsn(i) + irecvcountn(i)
      END DO

      CALL MPI_ALLGATHERV(listparentstemp(1),nlistupdated,MPI_INTEGER,
     &        listparents,irecvcountn,idisplacementsn,MPI_INTEGER,
     &        MPI_COMM_WORLD,ierr)
c
c--Unset local
c
#ifdef MPIDEBUG
      print *,iproc,': nlistupdated ',nlistupdated
#endif

      DO i = 1, nlistupdated
         iupdated(listparentstemp(i)) = .FALSE.
      END DO
      nlistupdated = idisplacementsn(numproc) + irecvcountn(numproc)
c
c--Need to make sure no duplicates and all are set
c
      itemp = 0
      DO i = 1, nlistupdated
         IF (.NOT.iupdated(listparents(i))) THEN
            iupdated(listparents(i))=.TRUE.
            itemp = itemp + 1
            listparents(itemp) = listparents(i)
         ENDIF
      END DO
      nlistupdated = itemp
c
c--Calculate eospg values for particles done on other processes
c
c--Pressure and sound velocity from ideal gas law...
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP&shared(idisplacementsi,irecvcounti,llisttrans)
C$OMP&shared(vxyzu,rho,pr,vsound,ekcle,numproc,iproc,iphase,dustvar)
C$OMP&private(i,j,ipart)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               ipart = llisttrans(i) + 1

               IF (iphase(ipart).EQ.0)
     &              CALL eospg(ipart, vxyzu, rho, pr, vsound, ekcle,
     &              dustvar)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#endif

      IF (itiming) THEN
         CALL getused(td12)
         td1 = td1 + (td12 - td11)
      ENDIF

      ncompactlocal = ncompact
      nptmasstotlast = nptmasstot
c
c--FINISHED CALCULATION OF NEIGHBOURS, DENSITY, AND LONG-RANGE GRAVITY
c
c
c--BEGIN MPI SECTION
c
c--For MPI job, need to transfer quantities of particles that interact
c     with other MPI processes to those processes.
c
      IF (itiming) CALL getused(td21)
#ifdef MPI
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
c     &                  1,llistsentback(1,j+1),imhdevolREAL8,indexMPI3,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &             llistsentback(1,j+1),imhdevolREAL8,indexMPI3,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)

                              CALL MPI_SEND(Bevol,1,indexMPI3,j,44,
     &                             MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                              print *,iproc,' sent Bevol'
#endif
                              CALL MPI_TYPE_FREE(indexMPI3,ierr)
                           ENDIF

                           IF (idim_grow.EQ.idim) THEN
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL8,indexMPI,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL8,indexMPI,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI,ierr)

                              CALL MPI_SEND(rgrain,1,indexMPI,j,45,
     &                             MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                              print *,iproc,' sent rgrain'
#endif
                              CALL MPI_TYPE_FREE(indexMPI,ierr)

c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &                  1,llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                  llistsentback(1,j+1),MPI_REAL4,indexMPI1,ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

                           CALL MPI_SEND(vsound,1,indexMPI1,j,46,
     &                          MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent vsound'
#endif
                           CALL MPI_TYPE_FREE(indexMPI1,ierr)

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
                  CALL quit(1)
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
                  CALL quit(1)
               ENDIF

               IF (inumber.GT.0) THEN
                  CALL MPI_RECV(vxyzu(1,istart),inumber,i4REAL8,
     &                 iprocrec, 41, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i4REAL8, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 1 ',
     &                    iproc
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(dumrho(istart),inumber,MPI_REAL4,
     &                 iprocrec, 42, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck, ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 2 ',
     &                    iproc
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(iphase(istart),inumber,MPI_INTEGER1,
     &                 iprocrec, 43, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER1,icheck,ierr)
                  IF (icheck.NE.inumber) THEN
                     WRITE (*,*) 'ERROR - icheck.NE.inumber 3 ',
     &                    iproc
                     CALL quit(1)
                  ENDIF

                  IF (imhd.EQ.idim) THEN
                     CALL MPI_RECV(Bevol(1,istart),inumber,
     &                    imhdevolREAL8, iprocrec, 44, MPI_COMM_WORLD, 
     &                    istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, imhdevolREAL8, 
     &                    icheck, ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 4 ',
     &                       iproc
                        CALL quit(1)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got Bevol from ',iprocrec,icheck,
     &                    radkernel,xyzmh(1,1)
#endif
                  ENDIF

                  IF (idim_grow.EQ.idim) THEN
                     CALL MPI_RECV(rgrain(istart),inumber,
     &                    MPI_REAL8, iprocrec, 45, MPI_COMM_WORLD,
     &                    istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, MPI_REAL8,
     &                    icheck, ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 5 ',
     &                       iproc
                        CALL quit(1)
                     ENDIF
#ifdef MPIDEBUG
                     print *,iproc,': got rgrain from ',iprocrec,
     &                    icheck,istart
#endif                     

                     CALL MPI_RECV(vsound(istart),inumber,MPI_REAL4,
     &                 iprocrec, 46, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, MPI_REAL4, icheck,ierr)
                     IF (icheck.NE.inumber) THEN
                        WRITE (*,*) 'ERROR - icheck.NE.inumber 6 ',
     &                       iproc
                        CALL quit(1)
                     ENDIF

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
         CALL quit(1)
      ENDIF
      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERROR1 - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit(1)
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
c--END -- Have transferred quantities to other MPI processes
c
c
c--Calculate ONLY neighbours lists for remote particles
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(inumbertotal,ntot,listp,nlst_end,iphase,hasghost)
C$OMP& shared(npart,xyzmh,acc,nneigh)
C$OMP& shared(ncompact,icompact,ivar,ijvar,iscurrent)
C$OMP& shared(iupdated,nlistupdated,listparents,iprint,itime)
C$OMP& private(fsx,fsy,fsz,epot,nneighlocal,ncompacthere,icompacthere)
C$OMP& private(i,k,ipart,jpart,numberupdated)
      DO i = 1, inumbertotal

         ipart = ntot + i

         listp(nlst_end + i) = ipart
         hasghost(ipart) = .FALSE.

         CALL treef(ipart,itime,npart,ntot,xyzmh,acc,0,fsx,fsy,fsz,epot)
         nneighlocal = nneigh(ipart)
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
            CALL quit(1)
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
                  CALL quit(1)
               ENDIF
               listparents(numberupdated) = jpart
            ENDIF
         END DO
      END DO
C$OMP END PARALLEL DO

      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(imhdevolREAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)
      CALL MPI_TYPE_FREE(i8REAL8,ierr)
c      CALL MPI_TYPE_FREE(i12REAL8,ierr)
      CALL MPI_TYPE_FREE(i15REAL8,ierr)

#endif
c
c--END MPI SECTION
c
c
c--Copy changed values onto ghost particles
c
C$OMP PARALLEL DO SCHEDULE (runtime) default(none)
C$OMP& shared(npart,ntot,ireal,rho,dumrho,xyzmh,gradhs)
#ifdef NONIDEAL
C$OMP& shared(eta_nimhd,nden_nimhd)
#endif
C$OMP& private(i,j)
      DO i = npart + 1, ntot
         j = ireal(i)
         rho(i) = rho(j)
         dumrho(i) = dumrho(j)
         xyzmh(5,i) = xyzmh(5,j)
         gradhs(1,i) = gradhs(1,j)
#ifdef NONIDEAL
         eta_nimhd(:,i) = eta_nimhd(:,j)
         nden_nimhd(:,i) = nden_nimhd(:,j)
#endif
      END DO
C$OMP END PARALLEL DO
c
c--BEGIN FIND CANDIDATE SINK PARTICLE
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_end,listp,iphase,rho,divv,curlv,gradhs,xyzmh,cnormk)
C$OMP& shared(nptmass,iptmass,listpm,varmhd,Bxyz,rhomin)
C$OMP& shared(Bextx,Bexty,Bextz,radcrit2)
C$OMP& shared(grav_accel,dvxyzu)
C$OMP& private(j,l,ipart,rhoi,dhdrhoi,iptcur)
C$OMP& private(rxxi,rxyi,rxzi,ryyi,ryzi,rzzi)
C$OMP& private(dalphaxi,dalphayi,dalphazi,dbetaxi,dbetayi,dbetazi)
C$OMP& private(denom,ddenom,gradalphaxi,gradalphayi,gradalphazi)
C$OMP& private(gradbetaxi,gradbetayi,gradbetazi,term)
C$OMP& private(dgammaxi,dgammayi,dgammazi)
C$OMP& private(gradgammaxi,gradgammayi,gradgammazi)
C$OMP& reduction(MAX:rhonext)
      DO j = 1, nlst_end
         ipart = listp(j)
         IF (iphase(ipart).EQ.0) THEN
            rhoi = rho(ipart)
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
         ENDIF
c
c--Also store gravitational acceleration for dust growth
c
         IF (idimHY09.EQ.idim) grav_accel(1:3,ipart) = dvxyzu(1:3,ipart)
      END DO
C$OMP END PARALLEL DO
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
            CALL quit(1)
         ENDIF
c
c--Make sure that all neighbours of point mass candidate are being
c     done on this time step. Otherwise, not possible to accrete
c     them to form a point mass and it may create a point mass without
c     accreting many particles!
c
c--NOTE: For MPI job, still need to check for neighbours on another process
c      For this, need to check within forcei when particles are sent back
c      Only need to check particle that has passed all local tests and
c      also has neighbours on other process.
c      Need to check that all its neighbours are currently being evolved.
c      Also need to check that it is the candidate sink particle with the
c      greatest density on all MPI processes (only form 1 sink at a time).
c      These things only need to be checked if there is a particle for
c      which icreate=1 on any MPI process.
c
         IF ((2.0*xyzmh(5,irhonex)).LT.hacc) THEN
            IF (iaccnoisy) THEN
             WRITE(iprint,*)'Ptmass creation passed h ',xyzmh(5,irhonex)
            ENDIF
            iaccph1 = iaccph1 + 1 
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
               IF (iaccnoisy) THEN 
                 WRITE(iprint,*) ' and all local particles on step'
               ENDIF
               iaccph2 = iaccph2 + 1
            ELSE
               IF (iaccnoisy) THEN
                 WRITE(iprint,*) ' but not all local particles on step'
               ENDIF
               iaccf5 = iaccf5 + 1
            ENDIF
         ELSE
            IF (iaccnoisy) THEN
              WRITE(iprint,*) 'Ptmass creation failed on h ',
     &             xyzmh(5,irhonex)
            ENDIF
            iaccf6 = iaccf6 + 1
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
c      take the process with the largest rank first
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
c      processes know this.
c
       IF (iproccreate.EQ.iproc) THEN
c
c--Check for neighbours on other processes - only need to loop over "sendback"
c      lists to see which processes need to check
c
          IF (nneighsentanyatall) THEN
             DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (nneighsentany(j+1)) THEN
                        CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8,
     &                           ierr)
                        CALL MPI_TYPE_COMMIT(i5REAL8,ierr)

                        CALL MPI_SEND(xyzmh(1,irhonex),1,i5REAL8,j,60,
     &                       MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUGS
      print *,iproc,': sent irhonex data to ',j
#endif

                      CALL MPI_TYPE_FREE(i5REAL8,ierr)
c
c--Receive back information on whether neighbouring particles are all on 
c      current timestep or not
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
                  CALL quit(1)
               ENDIF
               CALL MPI_RECV(xyzmh(1,itest), 5, MPI_REAL8,
     &              MPI_ANY_SOURCE, 60, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (inumber.NE.5) THEN
                  WRITE (iprint,*) 'ERROR - Particle test ',inumber
                  CALL quit(1)
               ENDIF
c
c--Find neighbours of this particle and test to see whether they are all
c      being evolved on the current step
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
c
c--END FIND CANDIDATE SINK PARTICLE
c
      IF (nwarnup.GT.0) THEN
         WRITE (iprint,*) 'WARNING: restricted h jump (up) ',
     &        nwarnup,' times'
      ENDIF
      IF (nwarndown.GT.0) THEN
         WRITE (iprint,*) 'WARNING: restricted h jump (down) ',
     &        nwarndown,' times'
      ENDIF
      IF (nbisection.GT.0) THEN
         WRITE(iprint,*) 'WARNING: used bisection on ',nbisection,
     &       ' particles'
      ENDIF

      IF (itiming) THEN
         CALL getused(td22)
         td2 = td2 + (td22 - td21)
      ENDIF

      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine densityiterate')

#ifdef MPIDEBUG
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      IF (itime.EQ.000)
c           print *,iproc,': Exited densityiterate'
c      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      CALL quit(1)
#endif

c      print *,'Exited densityiterate ',itime

      RETURN
      END



c=============================================================================
c
c--Subroutine for calculating rho, gradh and number of neighbours
c
      SUBROUTINE densitygradh(ipart,itime,npart,ntot,nlst_doing,
     &     xyzmh,numneighi,dvxyzu)

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax2), dvxyzu(4,idim3)

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
      INCLUDE 'COMMONS/call'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
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

      IF (icall.NE.4) THEN
         IF (nlst_doing.GT.nptmasstot .OR. nptmasstot.NE.nptmasstotlast
     &        .OR. iptintree.EQ.2 .OR. initialptm.EQ.5) THEN
            CALL treef(ipart,itime,npart,ntot,xyzmh,acc,igphi,
     &           fsx,fsy,fsz,epot)

            IF (iptintree.NE.2) THEN
               IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
                  iptcur = listrealpm(ipart)
                  IF (listpm(iptcur).NE.ipart) THEN
                     WRITE (*,*) 'ERROR - listpm(iptcur).NE.ipart'
                     CALL quit(1)
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
               CALL quit(1)
            ENDIF
            fsx = gravxyzpstore(1,iptcur)
            fsy = gravxyzpstore(2,iptcur)
            fsz = gravxyzpstore(3,iptcur)
            epot = gravxyzpstore(4,iptcur)
         ENDIF
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
      IF (iptintree.EQ.0 .AND. igphi.NE.-1 .AND.
     &     iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
         CALL gptfromgas(ipart, npart, ntot, xyzmh, dvxyzu)
      ENDIF

      RETURN
      END
