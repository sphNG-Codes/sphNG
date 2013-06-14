      SUBROUTINE divBsubcycle(nsubsteps,dtmax,npart,xyzmh,rho,
     &           Bxyzpnew,divBmax)
  
      INCLUDE 'idim'

      REAL dtmax,divBmax
      INTEGER npart
      DIMENSION xyzmh(5,idim2)
      REAL*4 rho(idim2)
      DIMENSION Bxyzpnew(imhdevol,imhd2)
c      DIMENSION dBevolxyz(imhdevol,imhd2)

      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/eosq'

      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/gradhterms'

      REAL Bxyzphalf(imhdevol,imhd2)
      REAL dBxyzp(imhdevol,imhd2)

      INTEGER nosweep
      REAL dW,dr
      LOGICAL*1 moresweep
      REAL dx,dy,dz

      COMMON /implicit3/ vari2(3,idim2)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine divBiterate')
c
c--Set up constants in Code units
c
      tolerance = 1.0e-3
      dtimax = dtmax/imaxstep
c
c--Set errors to zero for iteration start
c
      numoscillations = 0
      numequal = 0
      numcomp = 0
      ipos = 1

      ihasghostcount = 0
      ihasghost = 0
      icompact = 0
      Bxyzmax = 0.0

C$OMP PARALLEL default(none)
C$OMP& shared(icall,dtimax,dtmax,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij,varij2,vari2)
C$OMP& shared(xyzmh,dvtable,ddvtable,grwij,cnormk)
C$OMP& shared(nlst0,ihasghost,radkernel,gradhs)
C$OMP& shared(Bxyzpnew,Bxyzp,Bxyzphalf,dBxyzp,vsound,ncompact)
C$OMP& private(n,i,j,k,icompact,pmi,hj,hj21,hj41,grkerni,grkernj)
C$OMP& private(rxyi,ryzi,rxzi,rxxi,ryyi,rzzi)
C$OMP& private(dti,dx,dy,dz)
C$OMP& private(rij2,rij,rij1,dr,pmj,rhoi,rhoj,hi,hi21,hi41)
C$OMP& private(v2,v,index,dxx,index1,dgrwdx,grwtij,dW)
C$OMP& private(pmjdWrij1rhoj,runix,runiy,runiz,denom1)
C$OMP& private(dBxi,dByi,dBzi,projdB,B2i,vsigi,vs2i,rho1i)
C$OMP& reduction(+:ihasghostcount)
C$OMP& reduction(+:totalmagenergy)
C$OMP& reduction(MAX:Bxyzmax)
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(static)
      DO n = 1, ncompact
         i = ivar(3,n)

         IF(hasghost(i)) THEN
            ihasghostcount = ihasghostcount + 1
         ENDIF

         DO j = 1, imhdevol
            dBxyzp(j,i) = 0.
         END DO

         dBxi = 0.
         dByi = 0.
         dBzi = 0.

         pmi = xyzmh(4,i)
         hi = xyzmh(5,i)
         hi21 = 1./(hi*hi)
         hi41 = hi21*hi21
         rhoi = rho(i)
         rho1i = 1.0/rhoi

         B2i = Bxyzpnew(1,i)**2 + Bxyzpnew(2,i)**2 + Bxyzpnew(3,i)**2
         vs2i = vsound(i)**2 + B2i*rhoi
         ! just use regular fast MHD wave speed
c         vs2i = 10.* vs2i * nsubsteps !**2

         vsigi = sqrt(vs2i)

         IF (icall.EQ.1) THEN
            IF (isteps(i).EQ.0) THEN
               dti = dtmax*1.0d-12
            ELSE
               dti = dtimax*isteps(i)*1.0d-12
               dti = dtimax*isteps(i)
            ENDIF
         ELSEIF (n.LE.nlst0) THEN
            dti = dtimax*isteps(i)
         ELSE
            dti = dtimax/2.0*isteps(i)
         ENDIF

c         dti = dti/nsubsteps

         DO k = 1, ivar(1,n)
            icompact = ivar(2,n) + k
            j = ijvar(icompact)

            dx = xyzmh(1,i) - xyzmh(1,j)
            dy = xyzmh(2,i) - xyzmh(2,j)
            dz = xyzmh(3,i) - xyzmh(3,j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            rij = SQRT(rij2)
            rij1 = 1./rij
            dr = rij

            pmj = xyzmh(4,j)
            hj = xyzmh(5,j)
            hj21 = 1.0/(hj*hj)
            hj41 = hj21*hj21
            rhoj = rho(j)

            v2 = rij2*hi21
            v = rij/hi
            IF (v.LT.radkernel) THEN
               index = v2*ddvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dgrwdx = (grwij(index1) - grwij(index))*ddvtable
               grkerni = (grwij(index) + dgrwdx*dxx)*cnormk*hi41*
     &              gradhs(1,i)
            ELSE
               grkerni = 0.
            ENDIF

            v2 = rij2*hj21
            v = rij/hj
            IF (v.LT.radkernel) THEN
               index = v2*ddvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dgrwdx = (grwij(index1) - grwij(index))*ddvtable
               grkernj = (grwij(index) + dgrwdx*dxx)*cnormk*hj41*
     &              gradhs(1,j)
            ELSE
               grkernj = 0.
            ENDIF

            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1

c            B2j = Bxyz(1,j)**2 + Bxyz(2,j)**2 + Bxyz(3,j)**2
c--rough signal velocity only (no div v term)
c            vsigij = 0.5*(vsigi + sqrt(vsound(j)**2 + B2j/rhoj))
c            pmjdWrij1rhoj = pmj*dW*vsigij/rhoj

            varij(1,icompact) = grkerni * pmj * rho1i * rho1i
            varij(2,icompact) = grkerni * pmi * rho1i * rho1i
            varij(3,icompact) = grkernj * pmj / (rhoj * rhoj)
            varij(4,icompact) = rhoj

            varij2(1,icompact) = runix
            varij2(2,icompact) = runiy
            varij2(3,icompact) = runiz
         END DO

         vari(1,n) = rhoi

         ! don't need 0.5 divv term during substepping, only when part of full MHD equations
         vari2(1,n) = vs2i
         vari2(2,n) = 0.8*vsigi/hi !+ 0.5*divv(i)*rho1i
         vari2(3,n) = dti
      END DO
C$OMP END DO
C$OMP SINGLE

      IF (ihasghostcount.GE.1) ihasghost = 1
C$OMP END SINGLE
C$OMP DO SCHEDULE(static)
      DO i = npart + 1, npart + nghost*ihasghost
         j = ireal(i)
         DO k = 1, 7
c            dvdx(k,i) = dvdx(k,j)
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL
c
c--Begin iterating
c
 50   divBinitial = 0.

      DO nosweep = 1, nsubsteps
c      print *, 'it ',nosweep
c
c--Calculate div B and dBpsi
c
         CALL divB_dBpsi(Bxyzpnew)

c         DO n = 1, ncompact
c            i = ivar(3,n)           
c            DO j = 1, 3
c               dBevolxyz(j,i) = dBevolxyz(j,i) + varinew(j,i)
c            END DO
c            dpsidt = -vari2(1,n)*varinew(4,i) - vari2(2,n)*Bxyzpnew(4,i)
c            dBevolxyz(4,i) = dpsidt
c         END DO
c
c         GOTO 333



c         STOP
c
c--Can now calculate half-timestep values of B and psi
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompact,ivar,vari2,Bxyzphalf,Bxyzpnew,varinew,nosweep)
C$OMP& shared(dBxyzp)
C$OMP& private(n,i,dti,dpsidt)
C$OMP& reduction(MAX:divBinitial)
         DO n = 1, ncompact
            i = ivar(3,n)
c
c--Do 1st bit of full timestep integration
c
            dti = vari2(3,n) /256.0

            Bxyzphalf(1,i) = Bxyzpnew(1,i) + varinew(1,i)*dti

c            IF (i.EQ.1) THEN
c            print *,'Bx/rho ',Bxyzpnew(1,i),varinew(1,i),dti
c            ENDIF

            Bxyzphalf(2,i) = Bxyzpnew(2,i) + varinew(2,i)*dti
            Bxyzphalf(3,i) = Bxyzpnew(3,i) + varinew(3,i)*dti
            dpsidt = -vari2(1,n)*varinew(4,i) - vari2(2,n)*Bxyzpnew(4,i)
            Bxyzphalf(4,i) = Bxyzpnew(4,i) + dpsidt*dti

            DO j = 1,3
               dBxyzp(j,i) = dBxyzp(j,i) + varinew(j,i)*dti
            END DO
            dBxyzp(4,i) = dBxyzp(4,i) + dpsidt*dti


            IF (nosweep.EQ.1) THEN
               divBinitial = MAX(divBinitial,ABS(varinew(4,i)))
            ENDIF
c
c--Integrate to 1/2 timestep
c
            dti = 0.5 * vari2(3,n)

            Bxyzpnew(1,i) = Bxyzpnew(1,i) + varinew(1,i)*dti
            Bxyzpnew(2,i) = Bxyzpnew(2,i) + varinew(2,i)*dti
            Bxyzpnew(3,i) = Bxyzpnew(3,i) + varinew(3,i)*dti
            Bxyzpnew(4,i) = Bxyzpnew(4,i) + dpsidt*dti
         END DO
C$OMP END PARALLEL DO
c
c--Modify values on ghosts for half-timestep
c
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            DO k = 1, 4
               Bxyzpnew(k,i) = Bxyzpnew(k,j)
            END DO
         END DO
c
c--Then calculate derivatives of B and psi at half-timestep
c
         CALL divB_dBpsi(Bxyzpnew)
c
c--Them calculate full-timestep values of B and psi
c
         f22 = 255./256.

         divBmax = 0.
         divBmean = 0.
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompact,ivar,vari2,Bxyzphalf,Bxyzpnew,dBxyzp,varinew,f22)
C$OMP& private(n,i,dti,dpsidt)
C$OMP& reduction(MAX:divBmax)
C$OMP& reduction(+:divBmean)
         DO n = 1, ncompact
            i = ivar(3,n)

            divBmax = MAX(divBmax,ABS(varinew(4,i)))
            divBmean = divBmean + ABS(varinew(4,i))
c
c--Do the 2nd bit of full timestep integration
c
            dti = f22 * vari2(3,n)

            Bxyzpnew(1,i) = Bxyzphalf(1,i) + varinew(1,i)*dti
            Bxyzpnew(2,i) = Bxyzphalf(2,i) + varinew(2,i)*dti
            Bxyzpnew(3,i) = Bxyzphalf(3,i) + varinew(3,i)*dti
            dpsidt = -vari2(1,n)*varinew(4,i) -vari2(2,n)*Bxyzpnew(4,i)
            Bxyzpnew(4,i) = Bxyzphalf(4,i) + dpsidt*dti

            DO j = 1,3
               dBxyzp(j,i) = dBxyzp(j,i) + varinew(j,i)*dti
            END DO
            dBxyzp(4,i) = dBxyzp(4,i) + dpsidt*dti
         END DO
C$OMP END PARALLEL DO
c
c--Modify values on ghosts for next timestep
c
         IF (nosweep.LT.nsubsteps) THEN
            DO i = npart + 1, npart + nghost*ihasghost
               j = ireal(i)
               DO k = 1, 4
                  Bxyzpnew(k,i) = Bxyzpnew(k,j)
               END DO
            END DO
         ENDIF

c         print *,'it ',nosweep,divBmax,divBmean/ncompact,ncompact

 333     CONTINUE

      END DO ! End loop over nsubsteps

c
c--Then calculate derivatives of B and psi at the end
c
      CALL divB_dBpsi(Bxyzpnew)

c      divBmax = 0.
c      DO n = 1, ncompact
c         i=ivar(3,n)

c         DO j = 1, 3
c            dBevolxyz(j,i) = dBevolxyz(j,i) + 
c     &           dBxyzp(j,i)/(vari2(3,n)*nsubsteps)
c         END DO
c         dBevolxyz(4,i) = dBxyzp(4,i)/(vari2(3,n)*nsubsteps)
c      END DO

      print *,'Maximum divergence is ',divBmax,nsubsteps,divBinitial

c      GOTO 50

cC$OMP DO SCHEDULE(static)
c         DO i = npart + 1, npart + nghost*ihasghost
c            j = ireal(i)
c         END DO
cC$OMP END DO
cC$OMP END PARALLEL
c
c--And that done, return everything to ASS
c

c      print *,'Initial magnetic energy = ',totalmagenergy
      totalmagenergy = 0.

cC$OMP PARALLEL DO SCHEDULE(static) default(none)
cC$OMP& shared(Bxyznew,Bxyz,xyzmh,rho,dBxyz)
cC$OMP& private(n,i,j,k,icompact)
cC$OMP& reduction(+:totalmagenergy)
c
c--Copy arrays for all particles
c
c      DO n = 1, ncompact
c         i = ivar(3,n)
         
c         totalmagenergy = totalmagenergy + xyzmh(4,i)*(Bxyz(1,i)**2 +
c     &        Bxyz(2,i)**2 + Bxyz(3,i)**2)/rho(i)
c      END DO
cC$OMP END PARALLEL DO

c      print *,'Final magnetic energy   = ',totalmagenergy

c      stop

      RETURN

      END





      SUBROUTINE divB_dBpsi(Bxyzpnew)

      INCLUDE 'idim'

      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/current'

      REAL Bxyzpnew(imhdevol,imhd2)
c
c--Zero arrays for storing div B and dBpsi
c
C$OMP PARALLEL default(none)
C$OMP& shared(ivar,ijvar,vari,varij,varij2)
C$OMP& shared(Bxyzpnew,iscurrent,varinew,ncompact)
C$OMP& private(i,j,k,n,divBi,dBpsix,dBpsiy,dBpsiz,rhoi,rhoj)
C$OMP& private(icompact,pmjdWirho2i,pmidWirho2i,pmjdWjrho2j)
C$OMP& private(runix,runiy,runiz,dBx,dBy,dBz,projdB,tempval)
C$OMP& private(pmjdWirho1i)
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)
         DO k = 1, 4
            varinew(k,i) = 0.
         END DO
      END DO
C$OMP END DO
c
c--Calculate divB values without using separate subroutine
c
c--Particle I loop
c
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)

         divBi = 0.
         dBpsix = 0.
         dBpsiy = 0.
         dBpsiz = 0.

         rhoi = vari(1,n)
c
c--All the neighbours loop
c
         DO k = 1, ivar(1,n)
            icompact = ivar(2,n) + k
            j = ijvar(icompact)

            pmjdWirho2i = varij(1,icompact)
            pmjdWirho1i = pmjdWirho2i * rhoi
            pmidWirho2i = varij(2,icompact)
            pmjdWjrho2j = varij(3,icompact)
            rhoj = varij(4,icompact)

            runix = varij2(1,icompact)
            runiy = varij2(2,icompact)
            runiz = varij2(3,icompact)
c
c--Calculate div B (only contributions from within hi)
c
            dBx = Bxyzpnew(1,i)*rhoi - Bxyzpnew(1,j)*rhoj
            dBy = Bxyzpnew(2,i)*rhoi - Bxyzpnew(2,j)*rhoj
            dBz = Bxyzpnew(3,i)*rhoi - Bxyzpnew(3,j)*rhoj
            projdB = dBx*runix + dBy*runiy + dBz*runiz

c            print *,'projdB-n ',i,j,projdB,dBx,runix,
c     &           pmjdWirho1i*rhoi

            divBi = divBi - pmjdWirho1i*projdB
c
c--Calculate dBpsi
c
c--Add terms from particle j within hi
c
            tempval = pmjdWirho2i*Bxyzpnew(4,i)
            dBpsix = dBpsix - tempval*runix
            dBpsiy = dBpsiy - tempval*runiy
            dBpsiz = dBpsiz - tempval*runiz
c
c--Add terms to particle j (if j is active)
c
            IF (iscurrent(j)) THEN
               tempval = pmidWirho2i*Bxyzpnew(4,i)
C$OMP ATOMIC
               varinew(1,j) = varinew(1,j) + tempval*runix
C$OMP ATOMIC
               varinew(2,j) = varinew(2,j) + tempval*runiy
C$OMP ATOMIC
               varinew(3,j) = varinew(3,j) + tempval*runiz
c
c--Add terms to particle i from inactive j particle
c
            ELSE
               tempval = pmjdWjrho2j*Bxyzpnew(4,j)
               dBpsix = dBpsix - tempval*runix
               dBpsiy = dBpsiy - tempval*runiy
               dBpsiz = dBpsiz - tempval*runiz
            ENDIF
         END DO                 !J-loop
c
c--Store div B and dBpsi
c
C$OMP ATOMIC
         varinew(1,i) = varinew(1,i) + dBpsix
C$OMP ATOMIC
         varinew(2,i) = varinew(2,i) + dBpsiy
C$OMP ATOMIC
         varinew(3,i) = varinew(3,i) + dBpsiz

         varinew(4,i) = divBi

c      print *,'divB total ',divBi

      END DO
C$OMP END DO
C$OMP END PARALLEL

      RETURN

      END
