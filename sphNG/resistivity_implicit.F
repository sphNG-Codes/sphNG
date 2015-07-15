      SUBROUTINE resistivity_implicit(nlst_in,nlst_end,list,dtdo,
     &     itime,npart,ntot,Bxyz,xyzmh,vxyzu,rho,vsound,alphaMMpass,
     &     nit,error)
  
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION list(idim2)
      DIMENSION Bxyz(3,imhd2)
      DIMENSION xyzmh(5,mmax2)
      DIMENSION vxyzu(4,idim2)
      REAL*4 rho(idim2),vsound(idim2),alphaMMpass(isizealphaMM,idim2)

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/fullstep'

#ifdef MPICOPY
      INCLUDE 'COMMONS/mpicopy'
#endif

      PARAMETER (nswmax = 1000)
      PARAMETER (tolmhd = 1.0e-6)
c
c--This factor is 1.0 for Backwards Euler implicit method and 0.5 for 
c     Crank-Nicolson implicit method
c
      PARAMETER (fac = 1.0)

      INTEGER nosweep
      LOGICAL moresweep
      REAL rhoi, rhoj
      LOGICAL*1 ido,ido_j

      DIMENSION BxyzNew(3,imhd2)
      DIMENSION ido(idim2),ido_j(idim2)

      REAL maxerrB,maxerrB2,maxerrBlast,maxerrB2last
#ifdef MPIALL
      REAL maxerrBtot, maxerrB2tot
      LOGICAL moresweeptot
#endif

#ifdef MPIALL
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8, ierr)
#endif

      WRITE (*,*) 'Entered Resistivity_Implicit ',nlst_end,ncompact
c
c--Set errors to zero for iteration start
c
      dtimax = dtdo/imaxstep
      moresweep = .FALSE.

      numcomp = 0
      ipos = 1
      print *,' ENTERED resistivity_implicit ',itime

C$OMP PARALLEL default(none)
C$OMP& shared(xmaxerr1old,xmaxerr2old,ncompact,ivar)
C$OMP& shared(varinew,ido,ido_j)
C$OMP& private(i,n,k)
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, ntests
         xmaxerr1old(i) = 0.0
         xmaxerr2old(i) = 0.0
      END DO
C$OMP END DO
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)

         ido(i) = .FALSE.
         ido_j(i) = .FALSE.

         DO k = 1, 7
            varinew(k,i) = 0.
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL

      ihasghostcount = 0
      ihasghost = 0

C$OMP PARALLEL default(none)
C$OMP& shared(BxyzNew,Bxyz)
C$OMP& shared(ddvtable,vxyzu)
C$OMP& shared(icall,dtimax,dtdo,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij)
C$OMP& shared(xyzmh,dvtable,grwij,cnormk)
C$OMP& shared(ihasghost,iphase,ifullstep,iresist,etamhd)
C$OMP& shared(ntot,ncompactlocal,radkernel,gradhs,ncompact)
C$OMP& shared(ido,ido_j,vsound,alphaMMpass,alphamin,ifsvi)
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend)
#endif
C$OMP& private(n,i,j,k,l,rhoi,icompact,pmi,iparttree)
C$OMP& private(dti,dx,dy,dz)
C$OMP& private(rij2,rij,rij1,pmj,rhoj,hi,hj,hi21,hj21,hi41,hj41)
C$OMP& private(v2i,vi,v2j,vj,index,dxx,index1,dgrwdx,grwtij)
C$OMP& private(dWi,dWj,dWidrmj,dWidrmi,dWjdrmj)
C$OMP& private(etai,etaj,etaij,dWidrmjrhojrhoieta,dWidrmirhojrhoieta)
C$OMP& private(dWjdrmjrhojrhoieta)
C$OMP& private(B2i,B2j,vsigi,vsigj,etaart)
C$OMP& reduction(+:ihasghostcount)
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
      DO n = numberstart, numberend
#else
      DO n = 1, ncompact
#endif
         i = ivar(3,n)
c
c--Only do for gas
c
         IF (iphase(i).EQ.0) THEN
c
c--Needed for MPI code
c
            IF (i.GT.ntot) THEN
               iparttree = i + ntot + 2
            ELSE
               iparttree = i
            ENDIF
c
c--Set original and new values of B
c
            DO j = 1, 3
               BxyzNew(j,i) = Bxyz(j,i)
            END DO
c
c--Particle is on this MPI process and needs to be iterated over
c
            IF (n.LE.ncompactlocal) THEN
               IF(hasghost(i)) THEN
                  ihasghostcount = ihasghostcount + 1
               ENDIF
c
c--Set the time interval over which the particle needs to be evolved
c
               IF (icall.EQ.1) THEN
                  IF (isteps(i).EQ.0) THEN
                     dti = dtdo*1.0d-12
                  ELSE
                     dti = dtimax*isteps(i)*1.0d-12
                  ENDIF
               ELSEIF (ifullstep(i)) THEN
                  dti = dtimax*isteps(i)
               ELSE
                  dti = dtimax*isteps(i)/2.0
               ENDIF
            ENDIF
c
c--Calculate quantities which do not change during iterations
c
            pmi = xyzmh(4,iparttree)
            hi = xyzmh(5,iparttree)
            hi21 = 1./(hi*hi)
            hi41 = hi21*hi21
            rhoi = rho(i)

            IF (ABS(iresist).GE.2) THEN
               etai = etafunc(rhoi,vxyzu(4,i),etamhd,iresist)
c               print *,i,rhoi,vxyzu(4,i),etai
            ELSEIF (ABS(iresist).EQ.1) THEN
               etai = etamhd
            ELSE
               etai = 0.
            ENDIF

            B2i = Bxyz(1,i)**2 + Bxyz(2,i)**2 + Bxyz(3,i)**2
            vsigi = SQRT(vsound(i)**2 + B2i/rhoi)
            IF (ifsvi.GE.6) THEN
               etaart = alphaMMpass(2,i)*vsigi*hi
            ELSE
               etaart = alphamin(2)*vsigi*hi
            ENDIF
            IF (etai/etaart.GT.0.01) ido(i) = .TRUE.

            DO k = 1, ivar(1,n)
               icompact = ivar(2,n) + k
               j = ijvar(icompact)
c
c--Need to make sure that E and U values are loaded for non-active neighbours
c
               DO l = 1, 3
                  BxyzNew(l,j) = Bxyz(l,j)
               END DO
c
c--Calculate other quantities
c
               dx = xyzmh(1,iparttree) - xyzmh(1,j)
               dy = xyzmh(2,iparttree) - xyzmh(2,j)
               dz = xyzmh(3,iparttree) - xyzmh(3,j)
               rij2 = dx*dx + dy*dy + dz*dz + tiny
               rij = SQRT(rij2)
               rij1 = 1./rij

               pmj = xyzmh(4,j)
               rhoj = rho(j)

               hj = xyzmh(5,j)
               hj21 = 1./(hj*hj)
               hj41 = hj21*hj21
c
c--Set kernel quantities
c
               v2i = rij2*hi21
               vi = rij/hi

               v2j = rij2*hj21
               vj = rij/hj

               IF (vi.LT.radkernel) THEN
                  index = v2i*ddvtable
                  dxx = v2i - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                  grwtij = (grwij(index) + dgrwdx*dxx)*hi41
                  dWi = grwtij * cnormk
               ELSE
                  dWi = 0.
               ENDIF

               IF (vj.LT.radkernel) THEN
                  index = v2j*ddvtable
                  dxx = v2j - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                  grwtij = (grwij(index) + dgrwdx*dxx)*hj41
                  dWj = grwtij * cnormk
               ELSE
                  dWj = 0.
               ENDIF

               dWidrmj = 0.5*pmj*dWi*rij1*gradhs(1,i)
               dWidrmi = 0.5*pmi*dWi*rij1*gradhs(1,i)
               dWjdrmj = 0.5*pmj*dWj*rij1*gradhs(1,j)
c
c--Set resistivity terms
c
               IF (ABS(iresist).GE.2) THEN
                  etaj = etafunc(rhoj,vxyzu(4,j),etamhd,iresist)
               ELSEIF (ABS(iresist).EQ.1) THEN
                  etaj = etamhd
               ELSE
                  etaj = 0.
               ENDIF

               B2j = Bxyz(1,j)**2 + Bxyz(2,j)**2 + Bxyz(3,j)**2
               vsigj = SQRT(vsound(j)**2 + B2j/rhoj)
               IF (ifsvi.GE.6) THEN
                  etaart = alphaMMpass(2,j)*vsigj*hj
               ELSE
                  etaart = alphamin(2)*vsigj*hj
               ENDIF
               IF (etaj/etaart.GT.0.01) THEN
                  ido_j(i) = .TRUE.
               ENDIF

               etaij = etai + etaj

               dWidrmjrhojrhoieta = dWidrmj/(rhoj*rhoi)*etaij
               dWidrmirhojrhoieta = dWidrmi/(rhoj*rhoi)*etaij
               dWjdrmjrhojrhoieta = dWjdrmj/(rhoj*rhoi)*etaij

               varij(1,icompact) = rhoj
               varij(2,icompact) = dWidrmjrhojrhoieta
               varij(3,icompact) = dWidrmirhojrhoieta
               varij(4,icompact) = dWjdrmjrhojrhoieta
            END DO
            vari(1,n) = dti
            vari(2,n) = rhoi
         ENDIF
      END DO
C$OMP END DO
C$OMP END PARALLEL

      IF (ihasghostcount.GE.1) ihasghost = 1
c
c--Begin iterating
c
c      print *,iproc, 'begin iterations '

      maxerrBlast = 1.0E+30
      maxerrB2last = 1.0E+30

      DO nosweep = 1, nswmax
c         print *,iproc, 'it ',nosweep
c
c--Set error to zero for this iteration    
c
         maxerrB = 0.0
         maxerrB2 = 0.0
c
c--Calculate fluxlimiter values without using separate subroutine
c
C$OMP PARALLEL default(none)
C$OMP& shared(vari,ivar,ijvar,varij,varinew)
C$OMP& shared(BxyzNew,Bxyz,rho,itime)
C$OMP& shared(nosweep)
C$OMP& shared(npart,nghost,ihasghost,ireal,ncompact)
C$OMP& shared(xyzmh)
C$OMP& shared(icall)
C$OMP& shared(iphase,maxerrBlast,maxerrB2last)
C$OMP& shared(ncompactlocal,iscurrent,ido,ido_j)
#ifdef MPI
C$OMP& shared(istart,ntot,istartrec,inumbertotal,inumberreturned)
C$OMP& shared(numproc,iproc,inumofsends,maxnneighsentback,istatus)
C$OMP& shared(ierr,ireturned,iprocrec,nneighsentback,ipos,jpos,inumber)
C$OMP& shared(llistsentback,inumofreturns,jjj,inumberproc,iprocsend)
C$OMP& shared(inumbertotallocal,nneighsentanyatall,nneighsentany)
C$OMP& shared(transfer,nneightogetback,inumberindiv,inumbercumm)
#endif
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend,numbertodoherekeep,llisttrans)
C$OMP& shared(numproc,iproc,lblocklengths,indexMPI_INT1)
C$OMP& shared(irecvcounti,idisplacementsi,realtransfer5to15,i3REAL8)
#endif
C$OMP& private(i,j,k,n,rhoi,rhoj)
C$OMP& private(icompact,dti,Bnumeratorx,Bnumeratory,Bnumeratorz)
C$OMP& private(BnumeratorOldx,BnumeratorOldy,BnumeratorOldz)
C$OMP& private(Bdenominator,tempval)
C$OMP& private(dWidrmjrhojrhoieta,dWidrmirhojrhoieta,dWjdrmjrhojrhoieta)
C$OMP& private(Bnewxi,Bnewyi,Bnewzi,dBcombx,dBcomby,dBcombz)
#ifdef MPI
C$OMP& private(l)
#endif
#ifdef MPICOPY
C$OMP& private(ierr)
#endif
C$OMP& reduction(MAX:maxerrB)
C$OMP& reduction(MAX:maxerrB2)

C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO j = 1, 7
               varinew(j,i) = 0.
            END DO
         END DO
C$OMP END DO

C$OMP DO SCHEDULE(runtime)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            DO k = 1, 3
               BxyzNew(k,i) = BxyzNew(k,j)
            END DO
         END DO
C$OMP END DO

c
c--Particle I loop for calculating I-J quantities
c

C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, ncompact
#endif
            i = ivar(3,n)

            IF (iphase(i).EQ.0 .AND. (ido(i) .OR. ido_j(i))) THEN

               rhoi = vari(2,n)
c
c--NOTE: Needs to do this loop even for boundaryparticles because active
c     boundary particles will need to contribute to the varinew() 
c     quantities (i.e. diffusion terms) of particle j due to the way that
c     particle j only finds neighbours inside h_j or non-active particles
c     inside h_i.  The varinew() quantities of a boundaryparticle are
c     not used, but its contributions to j are.
c
c--Initialising counters to zero for this particle
c
               Bnumeratorx = 0.0
               Bnumeratory = 0.0
               Bnumeratorz = 0.0
               BnumeratorOldx = 0.0
               BnumeratorOldy = 0.0
               BnumeratorOldz = 0.0
               Bdenominator = 0.0
c
c--All the neighbours loop
c
               DO k = 1, ivar(1,n)
                  icompact = ivar(2,n) + k
                  j = ijvar(icompact)

                  rhoj = varij(1,icompact)
                  dWidrmjrhojrhoieta = varij(2,icompact)
                  dWidrmirhojrhoieta = varij(3,icompact)
                  dWjdrmjrhojrhoieta = varij(4,icompact)
c	    
c--Calculate numerator and denominator
c
                  Bnumeratorx = Bnumeratorx +
     &                 dWidrmjrhojrhoieta*BxyzNew(1,j)

                  Bnumeratory = Bnumeratory +
     &                 dWidrmjrhojrhoieta*BxyzNew(2,j)

                  Bnumeratorz = Bnumeratorz +
     &                 dWidrmjrhojrhoieta*BxyzNew(3,j)

                  BnumeratorOldx = BnumeratorOldx +
     &                 dWidrmjrhojrhoieta*Bxyz(1,j)

                  BnumeratorOldy = BnumeratorOldy +
     &                 dWidrmjrhojrhoieta*Bxyz(2,j)

                  BnumeratorOldz = BnumeratorOldz +
     &                 dWidrmjrhojrhoieta*Bxyz(3,j)

                  Bdenominator = Bdenominator +
     &                 dWidrmjrhojrhoieta
c
c--For current(j) need to add contribution due to i for h_j
c
                  IF (iscurrent(j) .AND. ido(j)) THEN
                     tempval = dWidrmirhojrhoieta*BxyzNew(1,i)
C$OMP ATOMIC
                     varinew(1,j) = varinew(1,j) + tempval

                     tempval = dWidrmirhojrhoieta*BxyzNew(2,i)
C$OMP ATOMIC
                     varinew(2,j) = varinew(2,j) + tempval

                     tempval = dWidrmirhojrhoieta*BxyzNew(3,i)
C$OMP ATOMIC
                     varinew(3,j) = varinew(3,j) + tempval

                     tempval = dWidrmirhojrhoieta*Bxyz(1,i)
C$OMP ATOMIC
                     varinew(4,j) = varinew(4,j) + tempval

                     tempval = dWidrmirhojrhoieta*Bxyz(2,i)
C$OMP ATOMIC
                     varinew(5,j) = varinew(5,j) + tempval

                     tempval = dWidrmirhojrhoieta*Bxyz(3,i)
C$OMP ATOMIC
                     varinew(6,j) = varinew(6,j) + tempval

C$OMP ATOMIC
                     varinew(7,j) = varinew(7,j) + dWidrmirhojrhoieta

                  ELSE
                     Bnumeratorx = Bnumeratorx +
     &                    dWjdrmjrhojrhoieta*BxyzNew(1,j)

                     Bnumeratory = Bnumeratory +
     &                    dWjdrmjrhojrhoieta*BxyzNew(2,j)

                     Bnumeratorz = Bnumeratorz +
     &                    dWjdrmjrhojrhoieta*BxyzNew(3,j)

                     BnumeratorOldx = BnumeratorOldx +
     &                    dWjdrmjrhojrhoieta*Bxyz(1,j)

                     BnumeratorOldy = BnumeratorOldy +
     &                    dWjdrmjrhojrhoieta*Bxyz(2,j)

                     BnumeratorOldz = BnumeratorOldz +
     &                    dWjdrmjrhojrhoieta*Bxyz(3,j)

                     Bdenominator = Bdenominator +
     &                    dWjdrmjrhojrhoieta

                  ENDIF
               END DO           !J-loop

C$OMP ATOMIC
               varinew(1,i) = varinew(1,i) + Bnumeratorx
C$OMP ATOMIC
               varinew(2,i) = varinew(2,i) + Bnumeratory
C$OMP ATOMIC
               varinew(3,i) = varinew(3,i) + Bnumeratorz
C$OMP ATOMIC
               varinew(4,i) = varinew(4,i) + BnumeratorOldx
C$OMP ATOMIC
               varinew(5,i) = varinew(5,i) + BnumeratorOldy
C$OMP ATOMIC
               varinew(6,i) = varinew(6,i) + BnumeratorOldz
C$OMP ATOMIC
               varinew(7,i) = varinew(7,i) + Bdenominator
            ENDIF
         END DO
C$OMP END DO

#ifdef MPICOPY
c
c--For MPICOPY job, unfortunately need to do MPI reduce because in the above
c     loop particle 'i' gives varinew contributions to particle 'j'.
c
c--Transfer varinew(1-7,i)
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO k = 1, 7
               realtransfer5to15((n-1)*7+k) = varinew(k,i)
               realtransfer5to15(7*idim+(n-1)*7+k) = 0.
            END DO
         END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
         print *,iproc,': about to do varinew(1-7) reduce ',ncompact
#endif
c
c--All processes reduce the varinew(1-7) data
c
         CALL MPI_ALLREDUCE(realtransfer5to15,
     &        realtransfer5to15(7*idim+1),7*ncompact,
     &        MPI_REAL8,MPI_SUM,
     &        MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Now do reduction of varinew(1-7)
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO k = 1, 7
               varinew(k,i) = realtransfer5to15(7*idim+(n-1)*7 + k)
            END DO
         END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
         print *,iproc,': varinew(1-7) reduce ',ncompact
#endif

C$OMP END SINGLE
#endif

#ifdef MPI
c
c--Need to add contributions to varinew(1-7) from other MPI processes
c
C$OMP SINGLE
      istart = ntot + 1
      istartrec = istart + inumbertotal
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back varinew(1-7) to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' varinew ',j
#endif
               CALL MPI_RECV(varinew(1,istartrec),7*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,153,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.7*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.7*nnsentback varinw'
                  CALL quit(1)
               ENDIF

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1

                  DO k = 1, 7
                     varinew(k,ipos) = varinew(k,ipos) + varinew(k,jpos)
                  END DO
               END DO

#ifdef MPIDEBUG
               print *,iproc,': got varinew from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif

            END DO
c
c--Other processes send the values back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                  print *,iproc,': sending varinew to ',iprocsend,
     &                 ' istart ',istart,' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(varinew(1,istart),7*inumber,MPI_REAL8,
     &                 iprocsend,153,MPI_COMM_WORLD,ierr)
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                  GOTO 7703
               ENDIF
            END DO
 7703       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO
C$OMP END SINGLE
#endif

c
c--Particle I loop for calculating I only quantities
c
C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, ncompactlocal
#endif
            i = ivar(3,n)
            rhoi = vari(2,n)

            IF (iphase(i).EQ.0 .AND. ido(i)) THEN

               dti = vari(1,n)

               Bnumeratorx = dti*varinew(1,i)*rhoi
               Bnumeratory = dti*varinew(2,i)*rhoi
               Bnumeratorz = dti*varinew(3,i)*rhoi
               BnumeratorOldx = dti*varinew(4,i)*rhoi
               BnumeratorOldy = dti*varinew(5,i)*rhoi
               BnumeratorOldz = dti*varinew(6,i)*rhoi
               Bdenominator = dti*varinew(7,i)*rhoi

               Bnewxi = (Bxyz(1,i)*(1.0 + (1.0-fac)*Bdenominator) 
     &              - fac*Bnumeratorx - (1.0-fac)*BnumeratorOldx)/
     &              (1.0 - fac*Bdenominator)
               Bnewyi = (Bxyz(2,i)*(1.0 + (1.0-fac)*Bdenominator) 
     &              - fac*Bnumeratory - (1.0-fac)*BnumeratorOldy)/
     &              (1.0 - fac*Bdenominator)
               Bnewzi = (Bxyz(3,i)*(1.0 + (1.0-fac)*Bdenominator) 
     &              - fac*Bnumeratorz - (1.0-fac)*BnumeratorOldz)/
     &              (1.0 - fac*Bdenominator)

               dBcombx = fac*(Bdenominator*Bnewxi - Bnumeratorx)
     &              +(1.0-fac)*(Bdenominator*Bxyz(1,i) - BnumeratorOldx)
               dBcomby = fac*(Bdenominator*Bnewyi - Bnumeratory)
     &              +(1.0-fac)*(Bdenominator*Bxyz(2,i) - BnumeratorOldy)
               dBcombz = fac*(Bdenominator*Bnewzi - Bnumeratorz)
     &              +(1.0-fac)*(Bdenominator*Bxyz(3,i) - BnumeratorOldz)
 667           CONTINUE

c               Bnewxi = (Bxyz(1,i) - Bnumeratorx)/
c     &              (1.0 - Bdenominator)
c               Bnewyi = (Bxyz(2,i) - Bnumeratory)/
c     &              (1.0 - Bdenominator)
c               Bnewzi = (Bxyz(3,i) - Bnumeratorz)/
c     &              (1.0 - Bdenominator)

c               dBcombx = Bdenominator*Bnewxi - Bnumeratorx
c               dBcomby = Bdenominator*Bnewyi - Bnumeratory
c               dBcombz = Bdenominator*Bnewzi - Bnumeratorz
c
c--And the error is...
c
c--Relative erro
c
               maxerrB2 = MAX(maxerrB2,ABS(BxyzNew(1,i) - Bnewxi))
               maxerrB2 = MAX(maxerrB2,ABS(BxyzNew(2,i) - Bnewyi))
               maxerrB2 = MAX(maxerrB2,ABS(BxyzNew(3,i) - Bnewzi))
c
c--Function error
c
               maxerrB = MAX(maxerrB,Bxyz(1,i) + dBcombx - BxyzNew(1,i))
               maxerrB = MAX(maxerrB,Bxyz(2,i) + dBcomby - BxyzNew(2,i))
               maxerrB = MAX(maxerrB,Bxyz(3,i) + dBcombz - BxyzNew(3,i))
c
c--Copy values
c
               BxyzNew(1,i) = Bnewxi
               BxyzNew(2,i) = Bnewyi
               BxyzNew(3,i) = Bnewzi

            ENDIF
         END DO ! I-loop
C$OMP END DO

c         CALL quit(1)

#ifdef MPICOPY
C$OMP SINGLE
c
c--Need to transfer BxyzNew(1-3,i) so that when exits, can load
c         into Bxyz, or, if another iteration is do, so that
c         the latest values for all particles are used.
c
c--Transfer BxyzNew(1-3,i)
c
      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     i3REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their BxyzNew(1-3) data
c
      CALL MPI_ALLGATHERV(BxyzNew,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,i3REAL8,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Put BxyzNew(1-3) data into correct places in BxyzNew
c
C$OMP DO SCHEDULE(runtime)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 3
                 BxyzNew(k,llisttrans(i) + 1) =
     &                 realtransfer5to15((i-1)*3 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': BxyzNew(1-3) move ',numbertodoherekeep
#endif
      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
C$OMP END SINGLE
#endif

C$OMP END PARALLEL 

#ifdef MPIALL
c
c--Need to do global check for maxerrB, maxerrB2 to check
c     for errors and convergence across all MPI processes
c
c--Applies to BOTH MPI and MPICOPY
c
         CALL MPI_ALLREDUCE(maxerrB,maxerrBtot,1,MPI_REAL8,MPI_MAX,
     &        MPI_COMM_WORLD,ierr)
         maxerrB = maxerrBtot
         CALL MPI_ALLREDUCE(maxerrB2,maxerrB2tot,1,MPI_REAL8,MPI_MAX,
     &        MPI_COMM_WORLD,ierr)
         maxerrB2 = maxerrB2tot
#endif

#ifdef MPIALL
         IF (iproc.EQ.0) THEN
#endif
         IF (MOD(nosweep,1).EQ.0) THEN
           PRINT 97001,nosweep,nswmax,maxerrB,maxerrB2
97001      FORMAT('GSIMPL: Finished iteration ',I4,' of ',I4,' (',
     &        1PE12.5,1PE12.5,')')
         ENDIF
#ifdef MPIALL
         ENDIF
#endif
c
c--The actual test
c
         IF(maxerrB.LE.tolmhd.AND.maxerrB2.LE.tolmhd) THEN
            PRINT *,"Complete with ",nosweep," iterations"
            GOTO 150
         ENDIF

         maxerrBlast = maxerrB
         maxerrB2last = maxerrB2
c
c--Test for convergence to non-zero value (incorrect minimum)
c     Must have equal value at least twice to stop detecting up and down
c     as non-convergence
c
c         GOTO 333

         DO itest = ntests,1,-1

c            GOTO 332

            IF (nosweep.GT.10) THEN

c               GOTO 330

               IF (maxerrB.GT.0.99999*xmaxerr1old(itest).AND.
     &              maxerrB.LT.1.00001*xmaxerr1old(itest)) THEN
                  DO iii = 1, numcomp
                     IF (maxerrB.GT.0.99999*xmaxerrcomp1(iii).AND.
     &                    maxerrB.LT.1.00001*xmaxerrcomp1(iii)) THEN
                        PRINT *,"GSIMPL: Non-convergence1 ",maxerrB,
     &                  xmaxerrcomp1(iii),(xmaxerr1old(ii),ii=1,ntests)

                        moresweep = .TRUE.
                        GOTO 333

                     ENDIF
                  END DO
                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
                  xmaxerrcomp1(numcomp) = maxerrB
                  PRINT *,"GSIMPL: Almost Non-convergence ",
     &                 maxerrB,(xmaxerr1old(ii),ii=1,ntests)
               ENDIF
 330           CONTINUE
               IF (maxerrB2.GT.0.99999*xmaxerr2old(itest).AND.
     &              maxerrB2.LT.1.00001*xmaxerr2old(itest)) THEN
                  DO iii = 1, numcomp
                     IF (maxerrB2.GT.0.99999*xmaxerrcomp2(iii).AND.
     &                    maxerrB2.LT.1.00001*xmaxerrcomp2(iii)) THEN
                        PRINT *,"GSIMPL: Non-convergence1-2 ",maxerrB2,
     &                  xmaxerrcomp2(iii),(xmaxerr2old(ii),ii=1,ntests)


                        moresweep = .TRUE.
                        GOTO 333

                     ENDIF
                  END DO
                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
                  xmaxerrcomp2(numcomp) = maxerrB2
                  PRINT *,"GSIMPL: Almost Non-convergence ",
     &                 maxerrB2,(xmaxerr2old(ii),ii=1,ntests)
               ENDIF
            ENDIF
 332        CONTINUE
            IF (itest.NE.1) THEN
               xmaxerr1old(itest) = xmaxerr1old(itest-1)
               xmaxerr2old(itest) = xmaxerr2old(itest-1)
            ELSE
               xmaxerr1old(itest) = maxerrB
               xmaxerr2old(itest) = maxerrB2
            ENDIF
         END DO
 333     CONTINUE

#ifdef MPIALL
c
c--Need to do global check for moresweep across all MPI processes
c
c--Applies to BOTH MPI and MPICOPY
c
         CALL MPI_ALLREDUCE(moresweep,moresweeptot,1,MPI_LOGICAL,
     &        MPI_LAND,MPI_COMM_WORLD,ierr)
         moresweep = moresweeptot
#endif

         IF (moresweep) GOTO 998

#ifdef MPI
c
c--If another iteration is required, transfer new values of EU0, ekcle(2-3)
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
#ifdef MPIDEBUG
                        print *,iproc,': sending EU0 to ',j,
     &                     ' starting at ',llistsentback(1,j+1),'+1 ',
     &                       nneighsentback(j+1),' list '
c     &             ,(llistsentback(kkk,j+1),kkk=1,nneighsentback(j+1))
#endif
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i3REAL8,indexEU,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i3REAL8,indexEU,ierr)
                CALL MPI_TYPE_COMMIT(indexEU,ierr)

                        CALL MPI_SEND(BxyzNew,1,indexEU,j,154,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexEU,ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent BxyzNew'
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
               print *,iproc,': receiving BxyzNew data ',
     &              inumberreturned,inumofreturns
#endif
               istart = ntot + inumbertotallocal + 1
               CALL MPI_RECV(BxyzNew(1,istart), idim, MPI_REAL8,
     &              i, 154, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
               print *,iproc,': got BxyzNew from ',iprocrec,inumber
#endif
               inumbertotallocal = inumbertotallocal + inumber/3
               IF (inumberindiv(inumberreturned)*3.NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumbercumm(inumberreturned).NE.istart) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm ',istart,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberproc ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit(1)
      ENDIF

#endif

      END DO ! Iterations loop
c
c--Maximum number of iterations reached
c
      WRITE (iprint,99001) nswmax, maxerrB, maxerrB2
99001 FORMAT('RES-IMPL: Warning. Maximum iterations reached ',
     &     I5,2(1X,1PE12.5))
      PRINT *,"GSIMPL: Warning. Maximum iterations reached"
      moresweep = .TRUE.
      GOTO 998
c
c--Output success
c
 150  CONTINUE
#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
      PRINT *,"Succeeded with ",nosweep," iterations ",maxerrB,maxerrB2
#ifdef MPI
      ENDIF
#endif

 998  CONTINUE

      nit = nosweep
      error = MAX(maxerrB,maxerrB2)
c
c--And that done, return everything to ASS
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompactlocal,ivar)
C$OMP& shared(iphase,Bxyz,BxyzNew,ido)
C$OMP& private(n,i,j)
      DO n = 1, ncompactlocal
         i = ivar(3,n)
         IF (iphase(i).EQ.0 .AND. ido(i)) THEN
            DO j = 1, 3
               Bxyz(j,i) = BxyzNew(j,i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO

 999  CONTINUE

#ifdef MPIALL
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
#endif

      WRITE (*,*) 'Exiting Resistivity_Implicit ',nlst_end,ncompact

      RETURN

      END !SUBROUTINE GSIMPL
