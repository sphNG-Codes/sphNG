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

      IMPLICIT NONE

      INCLUDE 'idim'
  
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

      INTEGER nswmax
      PARAMETER (nswmax = 100)
      REAL tol
      PARAMETER (tol = 0.001)
      REAL lincrit
      PARAMETER (lincrit = 0.0005)
 
      LOGICAL moresweep, moresweep2
      REAL, DIMENSION(idim_dustFluid2) :: dustvar
      REAL, DIMENSION(idim_dustFluid2) :: dustvarNew
      REAL, DIMENSION(5,mmax2) :: xyzmh
      REAL, DIMENSION(4,idim2) :: vxyzu
      REAL*4, DIMENSION(idim2) :: rho

      REAL rhoi, rhoj
      REAL pri, prj
      REAL tsi, tsj
      REAL maxerrS, maxerrSlast
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
      REAL itime
      REAL stoppingtime
      REAL term1, term2
      REAL tempval1, tempval2
      REAL s4term, s3term, s2term, s1term, s0term
      REAL S1i
      REAL maxerrE2
      REAL errorS
      REAL gradhi, gradhj
      REAL s

      INTEGER nosweep
      INTEGER numcomp
      INTEGER ipos
      INTEGER i, j, k, n, l
      INTEGER ihasghost, ihasghostcount
      INTEGER ntot, npart
      INTEGER iparttree
      INTEGER icompact
      INTEGER itest
      INTEGER ii, iii
      INTEGER iprint
      INTEGER nit
c-------------------------------------------------------
c                        VARIABLES                     |
c------------------------------------------------------|
c                                                      |
c dustvar = Dust variable                              |
c dustvarNew = Dust variable after implicit scheme     |
c rhoi, rhoj = Density for i and j                     |
c pri, prj = Pressure for i and j                      |
c tsi, tsj = Stopping time for i and j                 |
c                                                      |
c-------------------------------------------------------

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
      ipos = 1
C$OMP PARALLEL default(none)
C$OMP& shared(xmaxerr1old,ncompact)
C$OMP& shared(varinew,ivar)
C$OMP& shared(rho)
C$OMP& private(i,n,k)
C$OMP& private(rhoi,rhoj,tsi,tsj)
c
c--Set up values that won't change during sweeps
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, ntests
         xmaxerr1old(i) = 0.0
      END DO ! i = 1, ntests
C$OMP END DO
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)
c         ido(i) = .FALSE.
c         ido_j(i) = .FALSE.
         DO k = 1, 2
            varinew(k,i) = 0.0
         END DO ! k = 1, 2
      END DO ! n = 1, ncompact
C$OMP END DO
c C$OMP DO SCHEDULE(runtime)
c
c--Doing a loop to get the stopping time here (with rho) might be 
c   faster?
c
c   Need to set up an array for the stopping time
c
c      DO n = 1, ncompact
c         i = ivar(3,n)
c         IF (iphase(i).EQ.0) THEN
c            rhoi = rho(i)
c            tsi(i) = stoppingtime(i,rhoi) ! need to get rhoi in here?  
c            DO k = 1, ivar(1,n)
c               icompact = ivar(2,n) + k
c               j = ijvar(icompact)
c                  IF (iphase(j).EQ.0) THEN
c                     rhoj = rho(j)
c                     tsj(j) = stoppingtime(j,rhoj)
c                  ENDIF ! iphase(j).EQ.0
c            END DO ! k = 1, ivar(1,n)
c         ENDIF ! iphase(i).EQ.0
c      END DO ! n = 1, ncompact
c
c C$OMP END DO
C$OMP END PARALLEL
c
c - Initialise ghost counters
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
C$OMP& shared(gamma)
C$OMP& private(rij,rij2,rij1,hi,hj,hij,hi21,hj21,hi41,hj41)
C$OMP& private(vi,v2i,vj,v2j,index,dxx,index1,dgrwdx)
C$OMP& private(grwtij)
C$OMP& private(dWi,dWj,dWidrmjrhoj,dWjdrmjrhoj,dWidrmirhoj)
C$OMP& private(dWidrmjrhojpr,dWidrmirhojpr,dWjdrmjrhojpr)
C$OMP& private(iparttree,dx,dy,dz)
C$OMP& private(pmi,pmj,rhoi,rhoj,pri,prj)
C$OMP& private(dti,tsi,tsj)
C$OMP$ private(n,i,j,k,l,icompact)
C$OMP& private(gradhi,gradhj)

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
            ENDIF ! i.GT.ntot
c
c--Set original and new values of dustvar
c
            dustvarNew(i) = dustvar(i)
c
c--Particle on this MPI process needs to be iterated over
c
            IF (n.LE.ncompactlocal) THEN
               IF (hasghost(i)) THEN
                  ihasghostcount = ihasghostcount + 1
               ENDIF ! hasghost(i)
c
c--Set time interval over which particle needs to be evolved
c
               IF (icall.EQ.1) THEN
                  IF (isteps(i).EQ.0) THEN
                     dti = dtdo*1.0d-12
                  ELSE
                     dti = dtimax*isteps(i)*1.0d-12
                  ENDIF ! isteps(i).EQ.0
               ELSEIF (ifullstep(i)) THEN
                  dti = dtimax*isteps(i)
               ELSE
                  dti = dtimax*isteps(i)/2.0
               ENDIF    ! icall.EQ.1
            ENDIF       ! n.LE.ncompactlocal
c
c--Calculate quantites which do not change during iterations
c
            pmi = xyzmh(4,iparttree)
            hi = xyzmh(5,iparttree)
            hi21 = 1./(hi*hi)
            hi41 = hi21*hi21
            rhoi = rho(i)
            pri = pr(i)
            tsi = stoppingtime(i,rhoi,dustvarNew(i))
c            gradhi = gradhs(1,i)

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
                  dustvarNew(j) = dustvar(j)
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
                  tsj = stoppingtime(j,rhoj,dustvarNew(j))
                  hj = xyzmh(5,j)
                  hj21 = 1./(hj*hj)
                  hj41 = hj21*hj21
c                  gradhj = gradhs(1,j)
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
c
c--Save the stopping times of the dust grains as they won't change
c   between iterations as dustvar is not involved in their calculation
c                  
                  varij2(1,icompact) = tsi
                  varij2(2,icompact) = tsj
               ENDIF            ! iphase(j).EQ.0
            END DO              ! k = 1, ivar(1,n)
c
c--Store these values in an array so they are more readily available
c   in the cache (this will speed things up)
c
            vari(1,n) = dti
c            PRINT*,'dti = ',vari(1,n),n
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
         maxerrS = 0.0
C$OMP PARALLEL default(none)
C$OMP& shared(vari,ivar,ijvar,varij,varij2,varinew)
C$OMP& shared(dustvarNew,dustvar)
C$OMP& shared(npart,nghost,ncompact,ncompactlocal)
C$OMP& shared(xyzmh,rho)
C$OMP& shared(itime,ireal,icall,iphase,ihasghost)
C$OMP& shared(iscurrent)
C$OMP& shared(maxerrSlast)
C$OMP& shared(nosweep)
C$OMP& shared(moresweep)
c--ifdef MPI
c--ifdef MPICOPY
C$OMP& private(i,j,k,n,icompact)
C$OMP& private(rhoi,rhoj,tsi,tsj,dti)
C$OMP& private(term1,term2)
C$OMP& private(dWidrmjrhojpr,dWidrmirhojpr,dWjdrmjrhojpr)
C$OMP& private(tempval1,tempval2)
C$OMP& private(moresweep2)
C$OMP& private(s4term,s3term,s2term,s1term,s0term)
C$OMP& private(S1i)

C$OMP& reduction(MAX:maxerrS)
c
c--Set varinew to zero
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO j = 1, 2
               varinew(j,i) = 0.0
            END DO
         END DO ! n = 1, ncompact
C$OMP END DO
c
c--Set ghost values
c
C$OMP DO SCHEDULE(runtime)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            dustvarNew(i) = dustvarNew(j)
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
               term1 = 0.0
               term2 = 0.0
               DO k = 1, ivar(1,n)
                  icompact = ivar(2,n) + k
                  j = ijvar(icompact)
                  IF (iphase(j).EQ.0) THEN
                     rhoj = varij(1,icompact)
                     dWidrmjrhojpr = varij(2,icompact)
                     dWidrmirhojpr = varij(3,icompact)
                     dWjdrmjrhojpr = varij(4,icompact)
                     tsi = varij2(1,icompact)
                     tsj = varij2(2,icompact)

                     term1 = term1 + dWidrmjrhojpr*tsj*dustvarNew(j) /
     &                          (1 + dustvarNew(j)**2)
                     term2 = term2 + dWidrmjrhojpr*dustvarNew(j)*tsi
                     IF (iscurrent(j)) THEN
                        tempval1 = dWidrmirhojpr*tsi*dustvarNew(i)
     &                          /(1 + dustvarNew(i)**2)
                        tempval2 = dWidrmirhojpr*dustvarNew(i)*tsj
C$OMP ATOMIC
                        varinew(1,j) = varinew(1,j) - tempval1
C$OMP ATOMIC
                        varinew(2,j) = varinew(2,j) - tempval2
                     ELSE
                        term1 = term1 + dWjdrmjrhojpr*tsj*dustvarNew(j)
     &                          /(1 + dustvarNew(j)**2)
                        term2 = term2 + dWjdrmjrhojpr*dustvarNew(j)*tsi
                     ENDIF      ! iscurrent(j)
                  ENDIF         ! iphase(j).EQ.0
               END DO           ! k = 1, ivar(1,n) (j loop)
C$OMP ATOMIC
            varinew(1,i) = varinew(1,i) + term1
C$OMP ATOMIC
            varinew(2,i) = varinew(2,i) + term2
c               PRINT*, 'term1,i = ',term1,i
c--At this point still need the factor of dt(/2)/rhoi
            ENDIF               ! iphase(i).EQ.0
         END DO                 ! n = 1, ncompact (i loop)
C$OMP END DO

c
c--Particle I loop for calculating I only quantities
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompactlocal
            i = ivar(3,n)
            IF (iphase(i).EQ.0) THEN
               rhoi = vari(2,n)
               dti = vari(1,n)
               tsi = varij2(1,n)
               tsj = varij2(2,n)

               term1 = varinew(1,i)*dti
               term2 = varinew(2,i)*dti
c               PRINT*, 'term1 = ',term1,i
c               PRINT*, 'term2 = ',term2,i
c     
c--Set quartic coefficients
c     
c               IF (term1.EQ.0) THEN
c                  PRINT*, 'ERROR - leading term in dust1_implicit is
c     &                     zero'
c                  PRINT*, 'dti = ',dti,i
c                  PRINT*, 'varinew(1,i) = ',varinew(1,i), i
c                  PRINT*, 'rhoi = ',rhoi,i
c                  CALL quit(1)
c               ENDIF    
               s4term = term1
               s3term = 0.0
               s2term = 2*term1 + term2
               s1term = 1.0
               s0term = term1 + term2 - dustvar(i)
c               PRINT*, 's4term = ',s4term
c               PRINT*, 's3term = ',s3term
c               PRINT*, 's2term = ',s2term
c               PRINT*, 's1term = ',s1term
c               PRINT*, 's0term = ',s0term
c
c--Set the moresweep input
c
               moresweep2 = .FALSE.
c
c--Rescale so leading coefficient is 1
c
               IF (dustvarNew(i).GT.lincrit) THEN
                  IF (s4term.NE.0.) THEN   
                     s2term = s2term/s4term
                     s1term = s1term/s4term
                     s0term = s0term/s4term
                     CALL quartic_dust1(s2term,s1term,s0term,
     &                    dustvarNew(i),S1i,moresweep2,i,1)
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
c                        PRINT *,"Dust variable = ",dustvar(i)
c                        PRINT *,"Dust variable new = ",dustvarNew(i)
c                        PRINT *,"S1i = ",S1i 
cC$OMP END CRITICAL (isucceedprint)
c                        CALL quit(1)
C$OMP CRITICAL (moresweepset)
                        moresweep =.TRUE.
C$OMP END CRITICAL (moresweepset)

c                        dustvarNew(i) = dustvar(i)         
c                           RETURN
c                           EXIT
                        CYCLE
                     ENDIF
                  ELSE
                     CALL quadratic(s2term,s1term,s0term,dustvarNew(i),
     &                    S1i,moresweep2,i)
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
c                        PRINT *,"Dust variable = ",dustvar(i)
c                        PRINT *,"Dust variable new = ",dustvarNew(i)
c                        PRINT *,"S1i = ",S1i
C$OMP CRITICAL (moresweepset)
                        moresweep=.TRUE.
C$OMP END CRITICAL (moresweepset)
                        CYCLE
                     ENDIF
                  ENDIF   
               ELSEIF (dustvarNew(i).LE.lincrit) THEN
                  S1i = MAX(-s0term,0.) 
               ENDIF
c
c--Tests for negativity
c
               IF (S1i.LT.0.0) THEN
C$OMP CRITICAL (moresweepset)
                  PRINT*,"dust1_implicit: dustvar has gone negative",
     &                 i,xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)
                  PRINT*,"dustvarNew(i) = ",dustvarNew(i)
                  PRINT*,"S1i = ",S1i
C$OMP END CRITICAL (moresweepset)
                  moresweep=.TRUE.
                  CYCLE
               ENDIF ! S1i.LE.0.0
c
c--Error
c
               IF (S1i.GE.0.) THEN
                  IF (dustvarNew(i).GT.lincrit) THEN
                     maxerrS = MAX(maxerrS,ABS((dustvarNew(i)-S1i)/S1i))
                     dustvarNew(i) = S1i
                  ELSE
                     maxerrS = MAX(maxerrS,ABS((dustvarNew(i)-S1i)/S1i))
                     dustvarNew(i) = S1i
                  ENDIF
               ENDIF
c               PRINT*, 'maxerrS = ',maxerrS   
            ENDIF ! iphase(i).EQ.0
          END DO ! n = 1, compactlocal
C$OMP END DO

C$OMP DO SCHEDULE(runtime)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            dustvarNew(i) = dustvarNew(j)
         END DO
C$OMP END DO

C$OMP END PARALLEL

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
     &                maxerrE2,(xmaxerr1old(ii),ii=1,ntests)
               ENDIF
            ENDIF ! nosweep.GT.10
            IF (itest.NE.1) THEN
               xmaxerr1old(itest) = xmaxerr1old(itest-1)
            ELSE
               xmaxerr1old(itest) = maxerrS
            ENDIF ! itest.NE.1
         END DO ! itest = ntests, 1, -1
c
c--After this are MPI processes, then the end of
c   iterations
c
      END DO ! nosweep = 1, nswmax
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
      PRINT*,'Suceeded with ',nosweep,' iterations error is ',maxerrS
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
C$OMP& shared(xyzmh,vari)
C$OMP& private(n,i,j)
      DO n = 1, ncompactlocal
         i = ivar(3,n)
         IF (iphase(i).EQ.0) THEN
            dustvar(i) = dustvarNew(i)
         ENDIF
      END DO
C$OMP END PARALLEL DO

999   CONTINUE ! Return old values
      RETURN
      
      END ! SUBROUTINE dust1_implicit
