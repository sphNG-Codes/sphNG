      SUBROUTINE dust1_implicit(dtdo,itime,npart,ntot,xyzmh,rho,
     &      dustvar,nit,errorS,moresweep)
c***********************************************************
c                                                          *
c   Subroutine to implicitly solve the dust variable       *
c   evolution equation (Ballabio et al. 2018) using        *
c   backward Euler implicit scheme and Gauss Seidel        *
c   iterative method to solve for interactions between     *
c   particles. The implicit scheme forms a quartic which   *
c   is solved analytically.                                *
c                                                          *
c   DE 24/02/22                                            *
c***********************************************************

#ifdef MPIALL
#define IMPLICITNONE
#include "mpi_sup.h"
#else
      IMPLICIT NONE
#endif

      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif
  
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/fullstep'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/dustfluidgrains'
      INCLUDE 'COMMONS/dustimplstat'

      INTEGER nswmax
      PARAMETER (nswmax = 100)
      REAL tol
      PARAMETER (tol = 0.001)
      REAL lincrit
      PARAMETER (lincrit = 0.0005)
 
      LOGICAL moresweep, moresweep2
      REAL, DIMENSION(ndusttypes,idim_dustFluid2) :: dustvar
      REAL, DIMENSION(ndusttypes,idim_dustFluid2) :: dustvarNew
      REAL, DIMENSION(5,mmax2) :: xyzmh
      REAL, DIMENSION(4,idim2) :: vxyzu
      REAL*4, DIMENSION(idim2) :: rho

      REAL, DIMENSION(ndusttypes,idim_dustFluid2) :: ts_store

      REAL, DIMENSION(ndusttypes,2,idim_dustFluid2) :: varidust

      REAL, DIMENSION(ndusttypes) :: dustfraci, dustfracj
      SAVE dustfraci, dustfracj
C$OMP THREADPRIVATE(dustfraci, dustfracj)
      REAL, DIMENSION(ndusttypes) :: term1, term2
      SAVE term1, term2
C$OMP THREADPRIVATE(term1, term2)

      REAL ts_temp
      REAL rhoi, rhoj
      REAL pri, prj
      REAL dtimax, dtdo, dti
      REAL dtmax
      REAL rij, rij1, rij2
      REAL hi, hi21, hi41
      REAL hj, hj21, hj41
      REAL hij
      REAL vi, v2i
      REAL vj, v2j
      REAL grwtij
      REAL index, index1
      REAL dxx, dx, dy, dz
      REAL dgrwdx
      REAL dWi, dWj
      REAL dWidrmjrhoj, dWjdrmjrhoj, dWidrmirhoj
      REAL dWidrmjrhojpr, dWidrmirhojpr, dWjdrmjrhojpr
      REAL pmi, pmj
      REAL tempval1, tempval2
      REAL s4term, s3term, s2term, s1term, s0term
      REAL S1i
      REAL maxerrS, maxerrSlast
      REAL errorS
      REAL gradhi, gradhj
      REAL vsoundi, vsoundj
      REAL dustfracisum, dustfracjsum

      INTEGER itime
      INTEGER nosweep
      INTEGER numcomp
      INTEGER i, j, k, kk, n, l
      INTEGER ihasghost, ihasghostcount
      INTEGER ntot, npart
      INTEGER iparttree
      INTEGER icompact
      INTEGER itest
      INTEGER ii, iii
      INTEGER iprint
      INTEGER nit
      INTEGER iregime

#ifdef MPIALL
      REAL maxerrStot
      LOGICAL moresweeptot
#endif

#ifdef MPI
      INTEGER istart, istartrec, inumberreturned
      INTEGER iprocrec, ireturned
      INTEGER ipos, jpos, inumber, iprocsend, jjj
      INTEGER inumbertotallocal
      INTEGER idustREAL8, idust2REAL8,indexMPIDustVar
#endif

c-------------------------------------------------------
c                        VARIABLES                     |
c------------------------------------------------------|
c                                                      |
c dustvar = Dust variable                              |
c dustvarNew = Dust variable after implicit scheme     |
c rhoi, rhoj = Density for i and j                     |
c pri, prj = Pressure for i and j                      |
c                                                      |
c-------------------------------------------------------

#ifdef MPIALL
      CALL MPI_TYPE_CONTIGUOUS(ndusttypes, MPI_REAL8, idustREAL8, ierr)
      CALL MPI_TYPE_COMMIT(idustREAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(ndusttypes*2, MPI_REAL8,idust2REAL8,ierr)
      CALL MPI_TYPE_COMMIT(idust2REAL8, ierr)
#endif

c
c--This routine only needs to be called during evolution
c     (i.e. icall=2 or icall=3).
c
      IF (icall.NE.2 .AND. icall.NE.3) THEN
         WRITE (*,*) 'ERROR - dust1_implicit called ',icall
         CALL quit(0)
      ENDIF
c
c--Set max time step
c
      dtimax = dtdo/imaxstep
c
c--Set more sweep logical
c
      moresweep = .FALSE.
c
c--Things for convergence test
c
      numcomp = 0
      xmaxerr1old(:) = 0.0
c
c--Initialise ghost counters
c
      ihasghostcount = 0
      ihasghost = 0

C$OMP PARALLEL default(none)
C$OMP& shared(dustvar,dustvarNew)
C$OMP& shared(xyzmh,rho,pr)
C$OMP& shared(dtimax,dtdo,isteps,icall,ihasghost,hasghost)
C$OMP& shared(ntot,npart,nghost,ncompactlocal,ncompact)
C$OMP& shared(vari,ivar,ijvar,varij,varij2)
C$OMP& shared(grwij,dvtable,cnormk,radkernel)
C$OMP& shared(ddvtable,ifullstep,iphase)
C$OMP& shared(gradhs)
C$OMP& shared(gamma,vsound,ts_store,sgrain,densgrain)
C$OMP& private(rij,rij2,rij1,hi,hj,hij,hi21,hj21,hi41,hj41)
C$OMP& private(vi,v2i,vj,v2j,index,dxx,index1,dgrwdx)
C$OMP& private(grwtij)
C$OMP& private(dWi,dWj,dWidrmjrhoj,dWjdrmjrhoj,dWidrmirhoj)
C$OMP& private(dWidrmjrhojpr,dWidrmirhojpr,dWjdrmjrhojpr)
C$OMP& private(iparttree,dx,dy,dz)
C$OMP& private(pmi,pmj,rhoi,rhoj,pri,prj)
C$OMP& private(dti)
C$OMP$ private(n,i,j,k,l,icompact)
C$OMP& private(gradhi,gradhj)
C$OMP& private(dustfracisum,dustfracjsum,ts_temp)
C$OMP& private(iregime,vsoundi,vsoundj)

C$OMP& reduction(+:ihasghostcount)

C$OMP DO SCHEDULE(runtime)
c
c--Set up values that don't change during sweeps
c
      DO n = 1, ncompact
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
            ENDIF    ! i.GT.ntot
c
c--Set original and new values of dustvar
c
            DO l = 1, ndusttypes
               dustvarNew(l,i) = dustvar(l,i)
            END DO
c
c--Particle on this MPI process needs to be iterated over
c
            IF (n.LE.ncompactlocal) THEN
               IF (hasghost(i)) THEN
                  ihasghostcount = ihasghostcount + 1
               ENDIF ! hasghost(i)
c
c--Set time interval over which particle needs to be evolved
c     This routine only needs to be called during evolution
c     (i.e. icall=2 or icall=3).
c
               IF (ifullstep(i)) THEN
                  dti = dtimax*isteps(i)
               ELSE
                  dti = dtimax*isteps(i)/2.0
               ENDIF
            ENDIF    ! n.LE.ncompactlocal
c
c--Calculate quantites which do not change during iterations
c
            pmi = xyzmh(4,iparttree)
            hi = xyzmh(5,iparttree)
            hi21 = 1./(hi*hi)
            hi41 = hi21*hi21
            rhoi = rho(i)
            pri = pr(i)
            vsoundi = vsound(i)
c            gradhi = gradhs(1,i)

            !--calculate the stopping time for the dust
c            tsi(:) = stoppingtime_old(i,rhoi,dustvarNew(:,i))
            CALL dustvar2dustfrac(1.0/rhoi,dustvar(1,i),dustfraci,
     &           dustfracisum)
            DO l = 1, ndusttypes
c            DO l = ndusttypes, ndusttypes
               CALL get_ts(ts_temp,sgrain(l),densgrain(l),
     &              rhoi*(1.-dustfracisum),rhoi*dustfracisum,
     &              vsoundi,0.,iregime)            
C$OMP ATOMIC WRITE
               ts_store(l,i) = ts_temp
            END DO

c            print *,'sgrain ',i,sgrain(:)
c            print *,'densgrain ',densgrain(:)
c            print *,'rhoi,dfracisum,vsoundi ',rhoi,dustfracisum,vsoundi

            DO k = 1, ivar(1,n)
               icompact = ivar(2,n) + k
               j = ijvar(icompact)
c
c--For gass only
c        
               IF (iphase(j).EQ.0) THEN
c            
c--Need to make sure s values are loaded for non-active neighbours
c
                  DO l = 1, ndusttypes
                     dustvarNew(l,j) = dustvar(l,j)
                  END DO
c
c--Calculate other quantites
c
                  dx = xyzmh(1,iparttree) - xyzmh(1,j)
                  dy = xyzmh(2,iparttree) - xyzmh(2,j)
                  dz = xyzmh(3,iparttree) - xyzmh(3,j)
                  rij2 = dx*dx + dy*dy + dz*dz + 1.0d-30
                  rij = SQRT(rij2)
                  rij1 = 1./rij

                  pmj = xyzmh(4,j)
                  rhoj = rho(j)
                  prj = pr(j)
                  vsoundj = vsound(j)
                  hj = xyzmh(5,j)
                  hj21 = 1./(hj*hj)
                  hj41 = hj21*hj21
c                  gradhj = gradhs(1,j)

c
c--Save the stopping times of the dust grains.  It is assumed they 
c   won't change between iterations.  This is not quite true because
c   the stopping time does depend on the dust density which can
c   slowly change due to dust diffusion, but we assume this is small.
c
c                  tsj(:) = stoppingtime_old(j,rhoj,dustvarNew(:,j))
                  CALL dustvar2dustfrac(1.0/rhoj,dustvarNew(1,j),
     &                 dustfracj,dustfracjsum)
                  DO l = 1, ndusttypes
c                  DO l = ndusttypes, ndusttypes
                     CALL get_ts(ts_temp,sgrain(l),densgrain(l),
     &                    rhoj*(1.-dustfracjsum),rhoj*dustfracjsum,
     &                    vsoundj,0.,iregime)
C$OMP ATOMIC WRITE
                     ts_store(l,j) = ts_temp
                  END DO
c
c--Set kernel quantities
c
                  v2i = rij2*hi21
                  vi = rij/hi

                  v2j = rij2*hj21
                  vj = rij/hj

                  IF (vi.LT.radkernel) THEN
                     index = v2i*ddvtable
                     index1 = index + 1
                     dxx = v2i - index*dvtable
                     IF (index1.GT.itable) index1 = itable
                     dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                     grwtij = (grwij(index) + dgrwdx*dxx)*cnormk
                     dWi = grwtij*hi41!*gradhi
                  ELSE
                     dWi = 0.
                  ENDIF         ! vi.LT.radkernel
c
c--Do similar for the j particles... 
c
                  IF (vj.LT.radkernel) THEN
                     index = v2j*ddvtable
                     index1 = index + 1
                     dxx = v2j - index*dvtable 
                     IF (index1.GT.itable) index1 = itable 
                     dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                     grwtij = (grwij(index) + dgrwdx*dxx)*cnormk
                     dWj = grwtij*hj41!*gradhj
                     ELSE
                     dWj = 0. 
                  ENDIF         ! vj.LT.radkernel
c
c--Set useful quantities, i.e. the Lij quants
c
                  dWidrmjrhoj = 0.5*pmj*dWi*rij1/(rhoj*rhoi)
                  dWidrmirhoj = 0.5*pmi*dWi*rij1/(rhoj*rhoi)
                  dWjdrmjrhoj = 0.5*pmj*dWj*rij1/(rhoj*rhoi)
c
c--Include pressure
c        
                  dWidrmjrhojpr = 0.5*dWidrmjrhoj*(pri-prj)
                  dWidrmirhojpr = 0.5*dWidrmirhoj*(pri-prj)
                  dWjdrmjrhojpr = 0.5*dWjdrmjrhoj*(pri-prj)
c
c--Store these quantities in an array so they are more readily 
c   available in the cache (this will speed things up)
c
                  varij(1,icompact) = rhoj
                  varij(2,icompact) = dWidrmjrhojpr
                  varij(3,icompact) = dWidrmirhojpr
                  varij(4,icompact) = dWjdrmjrhojpr
               ENDIF            ! iphase(j).EQ.0
            END DO              ! k = 1, ivar(1,n)
c
c--Store these values in an array so they are more readily available
c   in the cache (this will speed things up)
c
            vari(1,n) = dti
            vari(2,n) = rhoi
         ENDIF          ! iphase(i).EQ.0
      END DO            ! n = 1, ncompact
C$OMP END DO
C$OMP END PARALLEL
c
c--Reset ghost counter
c
      IF (ihasghostcount.GE.1) ihasghost = 1
c
c--Begin iterating
c
      maxerrSlast = 1.0E+30

      DO nosweep = 1, nswmax
c
c--Set error to zero for this iteration
c



c         varidust(:,:,:) = 0.




         maxerrS = 0.0
C$OMP PARALLEL default(none)
C$OMP& shared(vari,ivar,ijvar,varij,varij2,varidust)
C$OMP& shared(dustvarNew,dustvar)
C$OMP& shared(npart,ntot,nghost,ncompact,ncompactlocal)
C$OMP& shared(xyzmh,rho)
C$OMP& shared(itime,ireal,icall,iphase,ihasghost)
C$OMP& shared(iscurrent)
C$OMP& shared(maxerrSlast)
C$OMP& shared(nosweep)
C$OMP& shared(moresweep)
C$OMP& shared(ts_store)
#ifdef MPI
C$OMP& shared(inumbertotal,numproc,iproc,inumofsends,maxnneighsentback)
C$OMP& shared(istatus,ierr,nneighsentback,ipos,jpos,inumber)
C$OMP& shared(istart,istartrec,inumberreturned,ireturned,iprocrec)
C$OMP& shared(llistsentback,inumofreturns,inumberproc,inumbercumm)
C$OMP& shared(inumberindiv,iprocsend)
C$OMP& shared(idustREAL8,idust2REAL8,indexMPIDustVar)
#endif
C$OMP& private(i,j,k,kk,l,n,icompact)
C$OMP& private(rhoi,rhoj,dti)
C$OMP& private(dWidrmjrhojpr,dWidrmirhojpr,dWjdrmjrhojpr)
C$OMP& private(tempval1,tempval2)
C$OMP& private(moresweep2)
C$OMP& private(s4term,s3term,s2term,s1term,s0term)
C$OMP& private(S1i)
C$OMP& reduction(MAX:maxerrS)
c
c--Set varidust to zero
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO j = 1,2
               DO l = 1, ndusttypes
                  varidust(l,j,i) = 0.0
               END DO
            END DO
         END DO ! n = 1, ncompact
C$OMP END DO
c
c--Set ghost values
c
C$OMP DO SCHEDULE(runtime)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            DO l = 1, ndusttypes
               dustvarNew(l,i) = dustvarNew(l,j)
            END DO
         END DO ! i = npart + 1, npart + nghost*ihasghost            
C$OMP END DO
c
c--Loop for i, calculating ij quantites
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)

            IF (iphase(i).EQ.0) THEN
               rhoi = vari(2,n)
c               tsi(:) = ts_store(:,i)

c--term1 is from the tsj sums: term1 = tsj*Lij*sj/(1+sj**2)
c   term2 is from the tsi sums: term2 = Lij*sj
c
c   where Lij = mj*(pri-prj)*Fij/rij/rhoj
c  
c   from the quartic: a4 = dt*term1/2/rhoi
c                     a3=0
c                     a2 = dt*term1/rhoi + dt*tsi*term2/2/rhoi
c                     a1 = 1
c                     a0 = dt*term1/2/rhoi + dt*tsi*term2/2/rhoi + si(n)
c
               DO l = 1, ndusttypes
                  term1(l) = 0.0
                  term2(l) = 0.0
               END DO
               DO k = 1, ivar(1,n)
                  icompact = ivar(2,n) + k
                  j = ijvar(icompact)
                  IF (iphase(j).EQ.0) THEN
                     rhoj = varij(1,icompact)
                     dWidrmjrhojpr = varij(2,icompact)
                     dWidrmirhojpr = varij(3,icompact)
                     dWjdrmjrhojpr = varij(4,icompact)
c                     tsj(:) = ts_store(:,j)

                     DO l = 1, ndusttypes

                      term1(l) = term1(l) + dWidrmjrhojpr*ts_store(l,j)*
     &                       dustvarNew(l,j) /(1.0 + dustvarNew(l,j)**2)
                        term2(l) = term2(l) + dWidrmjrhojpr*
     &                       dustvarNew(l,j)*ts_store(l,i)
                        IF (iscurrent(j)) THEN
                           tempval1 = dWidrmirhojpr*ts_store(l,i)*
     &                          dustvarNew(l,i)/(1.0+dustvarNew(l,i)**2)
                           tempval2 = dWidrmirhojpr*
     &                          dustvarNew(l,i)*ts_store(l,j)
C$OMP ATOMIC
                           varidust(l,1,j) = varidust(l,1,j)-tempval1
C$OMP ATOMIC
                           varidust(l,2,j) = varidust(l,2,j)-tempval2
                        ELSE
                        term1(l) = term1(l)+dWjdrmjrhojpr*ts_store(l,j)*
     &                       dustvarNew(l,j) /(1.0 + dustvarNew(l,j)**2)
                           term2(l) = term2(l) + dWjdrmjrhojpr*
     &                       dustvarNew(l,j)*ts_store(l,i)
                        ENDIF   ! iscurrent(j)

                     END DO

                  ENDIF         ! iphase(j).EQ.0
               END DO           ! k = 1, ivar(1,n) (j loop)
               DO l = 1, ndusttypes
c               DO l = ndusttypes, ndusttypes
C$OMP ATOMIC
                  varidust(l,1,i) = varidust(l,1,i) + term1(l)
C$OMP ATOMIC
                  varidust(l,2,i) = varidust(l,2,i) + term2(l)
               END DO
c--At this point still need the factor of dt(/2)/rhoi
            ENDIF               ! iphase(i).EQ.0
         END DO                 ! n = 1, ncompact (i loop)
C$OMP END DO

#ifdef MPI
c
c--Need to add contributions to varidust() from other MPI processes
c
C$OMP SINGLE
         istart = ntot + 1
         istartrec = istart + inumbertotal
         inumberreturned = 0
         DO i = 0, numproc - 1
            IF (iproc.EQ.i) THEN
c
c--Receive back varidust to be added on to local values
c
               DO j = 1, inumofsends
#ifdef MPIDEBUG
                  print *,iproc,': expecting to get ',inumofsends,
     &                 ' varidust ',j
#endif
                  CALL MPI_RECV(varidust(1,1,istartrec),
     &                 maxnneighsentback,idust2REAL8,
     &                 MPI_ANY_SOURCE,153,MPI_COMM_WORLD,istatus,ierr)
                  CALL MPI_GET_COUNT(istatus,idust2REAL8,ireturned,ierr)
                  iprocrec = istatus(MPI_SOURCE)
                  IF (ireturned.NE.nneighsentback(iprocrec+1)) THEN
                     WRITE (*,*) 'ERROR - ireturned.NE.nelements',
     &                    '*nnsentback varidust'
                     CALL quit(1)
                  ENDIF

                  DO l = 1, nneighsentback(iprocrec+1)
                     ipos = llistsentback(l,iprocrec+1)+1
                     jpos = istartrec + l - 1

                     DO k = 1, 2
                        DO kk = 1, ndusttypes
                           varidust(kk,k,ipos) = varidust(kk,k,ipos) + 
     &                          varidust(kk,k,jpos)
                        END DO
                     END DO
                  END DO

#ifdef MPIDEBUG
                  print *,iproc,': got varidust from ',iprocrec,
     &                 ' put into ',llistsentback(1,iprocrec+1)+1
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
                     print *,iproc,': sending varidust to ',iprocsend,
     &                   ' istart ',istart,' returned ',inumberreturned,
     &                   ' of ',inumofreturns
#endif
                     CALL MPI_SEND(varidust(1,1,istart),
     &                    inumber,idust2REAL8,
     &                    iprocsend,153,MPI_COMM_WORLD,ierr)
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                     GOTO 7703
                  ENDIF
               END DO
 7703          CONTINUE
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
         DO n = 1, ncompactlocal
            IF (moresweep) CYCLE
            i = ivar(3,n)
            IF (iphase(i).EQ.0) THEN
               rhoi = vari(2,n)
               dti = vari(1,n)

               DO l = 1, ndusttypes
                  term1(l) = varidust(l,1,i)*dti
                  term2(l) = varidust(l,2,i)*dti
               END DO
c               PRINT*, 'term1 = ',i,dti,term1(:)
c               PRINT*, 'term2 = ',i,dti,term2(:)
c     
c--Set quartic coefficients
c     
c               IF (term1.EQ.0) THEN
c                  PRINT*, 'ERROR - leading term in dust1_implicit is
c     &                     zero'
c                  PRINT*, 'dti = ',dti,i
c                  PRINT*, 'varidust(1,1,i) = ',varidust(1,1,i), i
c                  PRINT*, 'rhoi = ',rhoi,i
c                  CALL quit(1)
c               ENDIF
               DO l = 1, ndusttypes
                  IF (moresweep) CYCLE
c               DO l = ndusttypes, ndusttypes
c
c--Set the moresweep input
c
                  moresweep2 = .FALSE.

                  s4term = term1(l)
                  s3term = 0.0
                  s2term = 2*term1(l) + term2(l)
                  s1term = 1.0
                  s0term = term1(l) + term2(l) - dustvar(l,i)
c               PRINT*, 's4term = ',i,s4term
c               PRINT*, 's3term = ',i,s3term
c               PRINT*, 's2term = ',i,s2term
c               PRINT*, 's1term = ',i,s1term
c               PRINT*, 's0term = ',i,s0term
c
c--Rescale so leading coefficient is 1
c
                  IF (dustvarNew(l,i).GT.lincrit) THEN
                     IF (s4term.NE.0.) THEN   
                        s2term = s2term/s4term
                        s1term = s1term/s4term
                        s0term = s0term/s4term
                        CALL quartic_dust1(s2term,s1term,s0term,
     &                       dustvarNew(l,i),S1i,moresweep2,i,1)
c
c--Solve the quartic, S1i is the solution
c
                        IF (moresweep2) THEN                    
cC$OMP CRITICAL (isucceedprint)
c                        PRINT *,"QUARTIC FAILED, RETURNED TO IMPLICIT"
c                        PRINT *,"PRINTING PARTICLE QUANTITIES"
c                        PRINT *,"Particle index = ",i
c                         
c                        PRINT *,"Particle positions:"
c                        PRINT *,"x = ",xyzmh(1,i)
c                        PRINT* ,"y = ",xyzmh(2,i)
c                        PRINT *,"z = ",xyzmh(3,i)
c                        PRINT *,"dti = ",dti
c                        PRINT *,"tsi = ",tsi
c                        PRINT *,"rhoi = ",rhoi
c                        PRINT *,"Dust variable = ",dustvar(:,i)
c                        PRINT *,"Dust variable new = ",dustvarNew(i)
c                        PRINT *,"S1i = ",S1i 
cC$OMP END CRITICAL (isucceedprint)
c                        CALL quit(1)
cC$OMP CRITICAL (moresweepset)
C$OMP ATOMIC WRITE                           
                           moresweep =.TRUE.
cC$OMP END CRITICAL (moresweepset)
                           CYCLE

c                        dustvarNew(:,i) = dustvar(:,i)         
c                           RETURN
c                           EXIT
                        ENDIF
                     ELSE
                        CALL quadratic(s2term,s1term,s0term,
     &                    dustvarNew(l,i),S1i,moresweep2,i)
                        IF (moresweep2) THEN
c                        PRINT *,"QUADRATIC FAILED, RETURN TO IMPLICIT"
c                        PRINT *,"PRINTING PARTICLE QUANTITIES"
c                        PRINT *,"Particle index = ",i
c                         
c                        PRINT *,"Particle positions:"
c                        PRINT *,"x = ",xyzmh(1,i)
c                        PRINT* ,"y = ",xyzmh(2,i)
c                        PRINT *,"z = ",xyzmh(3,i)
c                        PRINT *,"dti = ",dti
c                        PRINT *,"tsi = ",tsi
c                        PRINT *,"rhoi = ",rhoi
c                        PRINT *,"Dust variable = ",dustvar(:,i)
c                        PRINT *,"Dust variable new = ",dustvarNew(i)
c                        PRINT *,"S1i = ",S1i
cC$OMP CRITICAL (moresweepset)
C$OMP ATOMIC WRITE                           
                           moresweep=.TRUE.
cC$OMP END CRITICAL (moresweepset)
                           CYCLE
                        ENDIF
                     ENDIF   
                  ELSEIF (dustvarNew(l,i).LE.lincrit) THEN
                     S1i = MAX(-s0term,0.) 
                  ENDIF
c
c--Tests for negativity
c
                  IF (S1i.LT.0.0) THEN
C$OMP CRITICAL (moresweepset)
                     PRINT*,"dust1_implicit: dustvar has gone negative",
     &                    l,i,xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)
                     PRINT*,"dustvarNew(l,i) = ",dustvarNew(l,i)
                     PRINT*,"S1i = ",S1i
C$OMP END CRITICAL (moresweepset)
C$OMP ATOMIC WRITE                           
                     moresweep=.TRUE.
                     CYCLE
                  ENDIF         ! S1i.LE.0.0
c
c--Error
c
                  IF (S1i.GE.0.) THEN
                     IF (dustvarNew(l,i).GT.lincrit) THEN
                        maxerrS = MAX(maxerrS,
     &                       ABS((dustvarNew(l,i)-S1i)/S1i))
                        dustvarNew(l,i) = S1i
                     ELSE
                        maxerrS = MAX(maxerrS,
     &                       ABS((dustvarNew(l,i)-S1i)/S1i))
                        dustvarNew(l,i) = S1i
                     ENDIF
                  ENDIF
               END DO ! l = 1, ndusttypes

c               PRINT*, 'maxerrS = ',maxerrS   
            ENDIF ! iphase(i).EQ.0
 900        CONTINUE
          END DO ! n = 1, compactlocal
C$OMP END DO

          IF (.NOT.moresweep) THEN

C$OMP DO SCHEDULE(runtime)
             DO i = npart + 1, npart + nghost*ihasghost
                j = ireal(i)
                DO l = 1, ndusttypes
                   dustvarNew(l,i) = dustvarNew(l,j)
                END DO
             END DO
C$OMP END DO

          ENDIF
C$OMP END PARALLEL

#ifdef MPIALL
c
c--Need to do global check for moresweep, and for maxerrE2, maxerrU2 to check
c     for errors and convergence across all MPI processes
c
         CALL MPI_ALLREDUCE(moresweep,moresweeptot,1,MPI_LOGICAL,
     &        MPI_LAND,MPI_COMM_WORLD,ierr)
         moresweep = moresweeptot
         CALL MPI_ALLREDUCE(maxerrS,maxerrStot,1,MPI_REAL8,MPI_MAX,
     &        MPI_COMM_WORLD,ierr)
         maxerrS = maxerrStot
#endif

         IF (moresweep) GOTO 999
c
c--Convergence test
c
         IF (maxerrS.LE.tol) THEN
            GOTO 150
         ENDIF
         maxerrSlast = maxerrS
c
c--Test for convergence to non-zero value (incorrect minimum)
c     Must have equal value at least twice to stop detecting up and down
c     as non-convergence
c
         DO itest = ntests,1,-1
            IF (nosweep.GT.10) THEN
               IF (maxerrS.GT.0.99999*xmaxerr1old(itest).AND.
     &              maxerrS.LT.1.00001*xmaxerr1old(itest)) THEN
                  DO iii = 1, numcomp
                     IF (maxerrS.GT.0.99999*xmaxerrcomp1(iii).AND.
     &                    maxerrS.LT.1.00001*xmaxerrcomp1(iii)) THEN
                         PRINT *,"GSIMPL: Non-convergence1 ",maxerrS,
     &                   xmaxerrcomp1(iii),(xmaxerr1old(ii),ii=1,ntests)
                         moresweep = .TRUE.
                        GOTO 998
                     ENDIF 
                  END DO ! iii = 1, numcomp
                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
                  xmaxerrcomp1(numcomp) = maxerrS
                  PRINT *,"GSIMPL: Almost Non-convergence ",
     &                maxerrS,(xmaxerr1old(ii),ii=1,ntests)
               ENDIF
            ENDIF ! nosweep.GT.10
            IF (itest.NE.1) THEN
               xmaxerr1old(itest) = xmaxerr1old(itest-1)
            ELSE
               xmaxerr1old(itest) = maxerrS
            ENDIF ! itest.NE.1
         END DO ! itest = ntests, 1, -1

#ifdef MPI
c
c--If another iteration is required, transfer new values of dustvarNew
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
c--Send active node data to be processed by other processes.  Data from al
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
                           print *,iproc,': sending dustvarNew to ',j,
     &                        ' starting at ',llistsentback(1,j+1),'+1 ',
     &                        nneighsentback(j+1),' list '
c     &             ,(llistsentback(kkk,j+1),kkk=1,nneighsentback(j+1))
#endif
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),idustREAL8,indexEU,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &             llistsentback(1,j+1),idustREAL8,indexMPIDustVar,ierr)
                CALL MPI_TYPE_COMMIT(indexMPIDustVar,ierr)

                           CALL MPI_SEND(dustvarNew,1,indexMPIDustVar,j,
     &                        154,MPI_COMM_WORLD, ierr)

                           CALL MPI_TYPE_FREE(indexMPIDustVar,ierr)
#ifdef MPIDEBUG
                           print *,iproc,' sent dustvarNew'
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
     &              nneightogetback(i+1).GT.0) THEN
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                  print *,iproc,': receiving dustvarNew data ',
     &                 inumberreturned,inumofreturns
#endif
                  istart = ntot + inumbertotallocal + 1
                  CALL MPI_RECV(dustvarNew(1,istart), idim, idustREAL8,
     &                 i, 154, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, idustREAL8, inumber, ierr)
                  iprocrec = istatus(MPI_SOURCE)
#ifdef MPIDEBUG
                  print *,iproc,': got dustvarN from ',iprocrec,inumber
#endif
                  inumbertotallocal = inumbertotallocal + inumber
                  IF (inumberindiv(inumberreturned).NE.inumber) THEN
                     WRITE (*,*) iproc,': ERR - inumberindiv dust ',
     &                    inumber,inumberindiv(inumberreturned),
     &                    inumberreturned
                     CALL quit(1)
                  ENDIF

                  IF (inumbercumm(inumberreturned).NE.istart) THEN
                     WRITE (*,*) iproc,': ERR - inumbercum dust',istart,
     &                    inumbercumm(inumberreturned),inumberreturned,
     &                    ntot,inumbertotallocal
                     CALL quit(1)
                  ENDIF
                  IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                     WRITE (*,*) iproc,': ERR - inumberproc dust',
     &                    iprocrec,
     &                    inumberproc(inumberreturned),inumberreturned
                     CALL quit(1)
                  ENDIF
               ENDIF
            ENDIF
         END DO
         IF (inumbertotallocal.NE.inumbertotal) THEN
            WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &           inumbertotallocal,inumbertotal
            CALL quit(1)
         ENDIF

#endif

      END DO ! nosweep = 1, nswmax  (i.e. the iterations loop)
c
c--Maximum number of iterations reached
c
      WRITE (iprint,99001) nswmax, maxerrS
99001 FORMAT('DUST1-IMPL: Warning. Maximum iterations reached ',
     &     I5,(1X,1PE12.5))
      PRINT *,"GSIMPL: Warning. Maximum iterations reached"
      moresweep = .TRUE.
      GOTO 998
c
c--Output success
c
150   CONTINUE
c
c--If convergence reached before max iterations
c
      ndustimpl_suc = ndustimpl_suc + 1
      ndustimpl_tot = ndustimpl_tot + nosweep

#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
      PRINT 99100,itime,nosweep,maxerrS,moresweep
99100 FORMAT('Succeeded at ',I10,' with ',I3,
     &     ' iterations, error is ',1PE12.5,' ',L)
#ifdef MPI
      ENDIF
#endif

c
c--If convergence max iterations reached
c
998   CONTINUE

      nit = nosweep
      errorS = maxerrS
c
c--Return everything
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompactlocal,ivar)
C$OMP& shared(iphase,dustvar,dustvarNew)
C$OMP& private(n,i,l)
      DO n = 1, ncompactlocal
         i = ivar(3,n)
         IF (iphase(i).EQ.0) THEN
            DO l = 1, ndusttypes
               dustvar(l,i) = dustvarNew(l,i)
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO

999   CONTINUE ! Return old values

#ifdef MPIALL
      CALL MPI_TYPE_FREE(idust2REAL8,ierr)
      CALL MPI_TYPE_FREE(idustREAL8,ierr)
#endif

      RETURN

      END ! SUBROUTINE dust1_implicit
