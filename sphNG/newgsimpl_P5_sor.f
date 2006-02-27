      SUBROUTINE GSIMPL(dtmax,nlst_in,nlst_end,nlstall,npart,list,
     &     ekcle,xyzmh,vxyzu,dedxyz,rho,moresweep,nit,error)
  
      INCLUDE 'idim'

      DIMENSION list(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION dedxyz(3,iradtrans)
      REAL*4 rho(idim)

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/neighbor_P'
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

      PARAMETER (nswmax = 200000)

      COMMON /getkap/ iflag
      
      REAL gradE(idim),flux(idim)
      REAL EU0(2,idim)
      REAL U1i,E1i
      REAL maxerrE,maxerrU
      REAL origEU(2,idim)

      DIMENSION dvdx(9,idim)
      
      INTEGER nosweep,lwst,hgst,stepa,sw
      REAL dW,vpi,dv,b2,b3,b4,b1,dr,hmean,E0ij,E1ij,Rh,Q
      REAL rhomean,cs,U0i,U0j,a,eta,vmu,dtsweeps
      LOGICAL*1 moresweep
      LOGICAL moresweep2,thirdsweepfail
      REAL dx,dy,dz,dvz,dvx,dvy,lightspeed
      LOGICAL qfailmoresw,trapez

      PARAMETER (icompactmax=60*idim)
      DIMENSION vari(3,idim),varij(7,icompactmax),varij2(3,icompactmax)
      DIMENSION ivar(2,idim),ijvar(icompactmax)      
      DIMENSION oneovermu(idim)

      PARAMETER (ntests=10)
      REAL xmaxerrold(ntests),xmaxerrcomp(ntests)

      REAL maxerrE2,maxerrU2
c
c--Set up constants in Code units
c
      lightspeed = c / udist * utime
      uradconst = radconst / uergcc
      dtimax = dtmax/imaxstep

c      print *,'entry gsimpl',dtmax,nlst_in,nlst_end,nlstall
c     &  ,npart,list(1),ekcle(1,1),ekcle(2,1),xyzmh(5,1),vxyzu(1,1),
c     &     isteps(1)


      boundtemp=12.037
c
c--Set errors to zero for iteration start
c
      numoscillations = 0
      numequal = 0
      numcomp = 0
      ipos = 1
      DO i = 1, ntests
         xmaxerrold(i) = 0.0
      END DO
c
c--Make compact list of neighbours
c
      ihasghostcount = 0
      ihasghost = 0
      icompact = 0
      DO n = nlst_in, nlst_end
        i = list(n)

        ivar(1,n) = nneigh(i)
        ivar(2,n) = icompact
c
c--NOTE: n loop cannot be parallised because of this loop:
c
        DO k = 1, nneigh(i)
           icompact = icompact + 1
           IF (icompact.GT.icompactmax) THEN
              WRITE (*,*) 'ERROR - compact not large enough'
              STOP
           ENDIF
           IF (k.GE.nlmax) THEN
              j = neighover(k-nlmax+1,ABS(neighb(nlmax,i)))
           ELSE
              j = neighb(k,i)
           ENDIF
           IF(i.EQ.j) THEN
              PRINT *,"TRAPIMPL: Particle interacting with itself"
              PRINT *,"TRAPIMPL: Error. Bye."
              STOP
           END IF
           ijvar(icompact) = j
        END DO
      END DO
c      print *, 'done here'
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,nlstall,list,EU0,uergg,ekcle)
C$OMP& shared(vxyzu,lightspeed,uradconst,oneovermu,iprint)
C$OMP& shared(icall,dtimax,dtmax,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij,ijvar2,varij2)
C$OMP& shared(xyzmh,dvtable,grwij,cnormk,alpha,beta)
C$OMP& shared(dvdx,nlst0,ihasghost,iflag,origEU)
C$OMP& private(n,i,j,k,pmi,rhoi,icompact)
C$OMP& private(dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi)
C$OMP& private(dvzdxi,dvzdyi,dvzdzi,dti,dx,dy,dz)
C$OMP& private(rij2,rij,rij1,dr,pmj,rhoj,hmean,hmean21,hmean41)
C$OMP& private(v2,v,index,dxx,index1,dgrwdx,grwtij,dW,dvx,dvy,dvz)
C$OMP& private(rhomean,dvdotdr,dv,vmu,dvdWmj05)
C$OMP& private(dWdrlightrhorhom,pmjdWrijrhoi)
C$OMP& private(pmjdWrunix,pmjdWruniy,pmjdWruniz)
C$OMP& reduction(+:ihasghostcount)
c
c--Copy arrays for all particles
c
C$OMP DO SCHEDULE(static)
      DO n = nlst_in, nlstall
      	 i = list(n)

         EU0(1,i) = ekcle(1,i)
         EU0(2,i) = vxyzu(4,i)
c
c--Note that CV and Kappa have already been done in ASS
c
         oneovermu(i) = GET1OVERMU(rho(i),vxyzu(4,i))
      END DO
C$OMP END DO
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(static)
      DO n = nlst_in, nlst_end
         i = list(n)

         origEU(1,i) = ekcle(1,i)
         origEU(2,i) = vxyzu(4,i)

         IF(hasghost(i)) THEN
            ihasghostcount = ihasghostcount + 1
         ENDIF

         dvxdxi = 0.
         dvxdyi = 0.
         dvxdzi = 0.
         dvydxi = 0.
         dvydyi = 0.
         dvydzi = 0.
         dvzdxi = 0.
         dvzdyi = 0.
         dvzdzi = 0.

         IF (icall.EQ.1) THEN
            IF (isteps(i).EQ.0) THEN
               dti = dtmax*1.0d-12
            ELSE
               dti = dtimax*isteps(i)*1.0d-12
            ENDIF
         ELSEIF (n.LE.nlst0) THEN
            dti = dtimax*isteps(i)
         ELSE
            dti = dtimax/2.0*isteps(i)
         ENDIF

         pmi = xyzmh(4,i)
         rhoi = rho(i)

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
            rhoj = rho(j)

            hmean = 0.5*(xyzmh(5,i) + xyzmh(5,j))
            hmean21 = 1./(hmean*hmean)
            hmean41 = hmean21*hmean21

c         print *,'loop1b ',k,j

            v2 = rij2*hmean21
            v = rij/hmean

c         print *,'loop1c ',k,j,v2

            index = v2/dvtable
            dxx = v2 - index*dvtable
            index1 = index + 1
            IF (index1.GT.itable) index1 = itable
            dgrwdx = (grwij(index1) - grwij(index))/dvtable
            grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
            dW = grwtij * cnormk

c         print *,'loop2 ',k,j

            dvx = vxyzu(1,i) - vxyzu(1,j)
            dvy = vxyzu(2,i) - vxyzu(2,j)
            dvz = vxyzu(3,i) - vxyzu(3,j)

            rhomean = 0.5*(rhoi+rhoj)

            dvdotdr = dvx*dx + dvy*dy + dvz*dz
            dv = dvdotdr/dr

            IF(dvdotdr.GT.0.0) THEN
               vmu = 0.0
            ELSE
ccc               vmu = hmean*dvdotdr/(rij**2+0.01*hmean**2)
               vmu = dv
            END IF

            dvdWmj05 = 0.5*pmj*dv*dW

            dWdrlightrhorhom = lightspeed*dW/dr*pmj/(rhoi*rhoj)

            pmjdWrijrhoi = pmj*dW*rij1/rhoi
            pmjdWrunix = pmjdWrijrhoi*dx
            pmjdWruniy = pmjdWrijrhoi*dy
            pmjdWruniz = pmjdWrijrhoi*dz
c
c--Calculates density(i) times the gradient of velocity
c
            dvxdxi = dvxdxi - dvx*pmjdWrunix
            dvxdyi = dvxdyi - dvx*pmjdWruniy
            dvxdzi = dvxdzi - dvx*pmjdWruniz
            dvydxi = dvydxi - dvy*pmjdWrunix
            dvydyi = dvydyi - dvy*pmjdWruniy
            dvydzi = dvydzi - dvy*pmjdWruniz
            dvzdxi = dvzdxi - dvz*pmjdWrunix
            dvzdyi = dvzdyi - dvz*pmjdWruniy
            dvzdzi = dvzdzi - dvz*pmjdWruniz

c         print *,'loop3 ',k,j

            varij(1,icompact) = pmj
            varij(2,icompact) = rhoj
            varij(3,icompact) = rhomean
            varij(4,icompact) = dvdWmj05
            varij(5,icompact) = vmu
            varij(6,icompact) = dWdrlightrhorhom
            varij(7,icompact) = pmj*dW

            varij2(1,icompact) = pmjdWrunix
            varij2(2,icompact) = pmjdWruniy
            varij2(3,icompact) = pmjdWruniz
         END DO

         dvdx(1,i) = dvxdxi
         dvdx(2,i) = dvxdyi
         dvdx(3,i) = dvxdzi
         dvdx(4,i) = dvydxi
         dvdx(5,i) = dvydyi
         dvdx(6,i) = dvydzi
         dvdx(7,i) = dvzdxi
         dvdx(8,i) = dvzdyi
         dvdx(9,i) = dvzdzi

         vari(1,n) = dti
         vari(2,n) = pmi
         vari(3,n) = rhoi
      END DO
C$OMP END DO
C$OMP SINGLE
c      print *, 'entry single',ekcle(2,22),npart,nghost,ihasghostcount

      IF (ihasghostcount.GE.1) ihasghost = 1
C$OMP END SINGLE
C$OMP DO SCHEDULE(static)
      DO i = npart + 1, npart + nghost*ihasghost
         j = ireal(i)
         DO k = 1, 9
            dvdx(k,i) = dvdx(k,j)
         END DO
      END DO
C$OMP END PARALLEL
c
c--Begin iterating
c
c      print *, 'begin iterations',ekcle(2,22)
      DO nosweep = 1, nswmax
c      print *, 'it ',nosweep
c
c--Set error to zero for this iteration    
c
         maxerrE2 = 0.0
         maxerrU2 = 0.0
c
c--Calculate fluxlimiter values without using separate subroutine
c
C$OMP PARALLEL default(none)
C$OMP& shared(lwst,hgst,stepa,list,vari,ivar,varij2,ijvar,varij,uergg)
C$OMP& shared(dedxyz,ekcle,EU0,uradconst,lightspeed,rho)
C$OMP& shared(oneovermu,nosweep,ipos,boundtemp)
C$OMP& shared(npart,nghost,ihasghost,ireal,dvdx,nlst_in,nlst_end)
C$OMP& shared(origEU,moresweep,xyzmh,gamma,alpha,beta,Rg,gmw)
C$OMP& shared(udens,radconst,nlstall,iflag,icall,ifsvi)
C$OMP& private(i,j,k,n,dedxi,dedyi,dedzi,pmi,rhoi,pmj,rhoj)
C$OMP& private(icompact,dti,pres_denominator,diffusion_numerator)
C$OMP& private(diffusion_denominator,pres_numerator,tfour)
C$OMP& private(radpresdenom,rhomean,dvdWmj05,vmu,dWdrlightrhorhom)
C$OMP& private(cs,vpi,bi,bj,b1,gradEi2,rpdiag,rpall,gradvPi,vsig)
C$OMP& private(pmjdWrunix,pmjdWruniy,pmjdWruniz,Eij1)
C$OMP& private(betaval,chival,gammaval,u4term,u1term,u0term)
C$OMP& private(gradE1i,tsr1i,moresweep2,U1i,E1i,dUcomb,dEcomb)
C$OMP& private(presioverrhoi2,presjoverrhoj2,pmjdW)
C$OMP& private(cvcold,ucold,ecold,residualE,residualU)
C$OMP& reduction(MAX:maxerrU2)
C$OMP& reduction(MAX:maxerrE2)
C$OMP DO SCHEDULE(static)
         DO n = nlst_in, nlst_end
            i = list(n)

c            if (icall.EQ.1 .and. n.EQ.1000) print *,i,EU0(1,i),EU0(2,i),
c     &           ekcle(2,i),ekcle(3,i)
            
            dedxi = 0.
            dedyi = 0.
            dedzi = 0.

            pmi = vari(2,n)
            rhoi = vari(3,n)

            DO k = 1, ivar(1,n)
               icompact = ivar(2,n) + k
               j = ijvar(icompact)

               pmj = varij(1,icompact)
               rhoj = varij(2,icompact)
               pmjdWrunix = varij2(1,icompact)
               pmjdWruniy = varij2(2,icompact)
               pmjdWruniz = varij2(3,icompact)
c
c--Calculates the gradient of E (where E=rho*e)
c
               Eij1 = rhoi*EU0(1,i) - rhoj*EU0(1,j)

               dedxi = dedxi - Eij1*pmjdWrunix
               dedyi = dedyi - Eij1*pmjdWruniy
               dedzi = dedzi - Eij1*pmjdWruniz
            END DO
            dedxyz(1,i) = dedxi
            dedxyz(2,i) = dedyi
            dedxyz(3,i) = dedzi

            gradE1i=SQRT(dedxi**2+dedyi**2+dedzi**2)

            tsr1i = ABS(gradE1i)/(EU0(1,i)*(rhoi**2)*ekcle(2,i))

            ekcle(4,i) = (2. + tsr1i ) / (6. + 3.0*tsr1i + tsr1i**2)
!          ekcle(4,i) = (120. + tsr1i ) / (360. + 3.0*tsr1i + tsr1i**2)
!          ekcle(4,i) = 1.0/3.0

            ekcle(5,i) = ekcle(4,i) + ekcle(4,i)**2 * tsr1i**2
!          eddington(i) = 1.0/3.0
c
c--NOTE: ***** Forcing lambda and eddington to be 1/3 *****
c
c            IF (rhoi.GT.10000.0) THEN
c               ekcle(4,i) = 1.0/3.0
c               ekcle(5,i) = 1.0/3.0
c            ENDIF

         END DO
C$OMP END DO
C$OMP DO SCHEDULE(static)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            ekcle(4,i) = ekcle(4,j)
            ekcle(5,i) = ekcle(5,j)
            dedxyz(1,i) = dedxyz(1,j)
            dedxyz(2,i) = dedxyz(2,j)
            dedxyz(3,i) = dedxyz(3,j)
         END DO
C$OMP END DO
c
c--Particle I loop
c
C$OMP DO SCHEDULE(static)
         DO n = nlst_in, nlst_end
            i = list(n)

c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) i,n
            
            dti = vari(1,n)
            pmi = vari(2,n)
            rhoi = vari(3,n)
c
c--Initialising counters to zero for this particle
c
            diffusion_numerator = 0.0
            diffusion_denominator = 0.0
            pres_numerator = 0.0
            pres_denominator = 0.0
c
c--All the neighbours loop
c
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl 1',
c     &           ivar(1,n)
            
            DO k = 1, ivar(1,n)
               icompact = ivar(2,n) + k
               j = ijvar(icompact)

               pmj = varij(1,icompact)
               rhoj = varij(2,icompact)
               rhomean = varij(3,icompact)
               dvdWmj05 = varij(4,icompact)
               vmu = varij(5,icompact)
               pmjdW = varij(7,icompact)
               dWdrlightrhorhom = varij(6,icompact)
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'loop',
c     &           k,j,dWdrlightrhorhom
c
c--Current mean sound speed
c
               presioverrhoi2 = Rg*EU0(2,i)/ekcle(3,i)*
     &              oneovermu(i)/rhoi/uergg
               presjoverrhoj2 = Rg*EU0(2,j)/ekcle(3,j)*
     &              oneovermu(j)/rhoj/uergg
               cs = (SQRT(gamma*presioverrhoi2*rhoi) + 
     &              SQRT(gamma*presjoverrhoj2*rhoj))/2.0
               vpi = ((-alpha*cs*vmu)+(beta*vmu**2))/rhomean
c
c--Work out numerator for pressure
c
c               pres_numerator = pres_numerator+
c     &              dvdWmj05*(vpi+presjoverrhoj2)
ccc     &              + pmjdW*cs/rhomean*1.0*(EU0(2,i)-EU0(2,j))
c
c               pres_denominator = pres_denominator+
c     &              dvdWmj05*presioverrhoi2/EU0(2,i)

               pres_numerator = pres_numerator+
     &              dvdWmj05*(vpi)

               pres_denominator = pres_denominator+
     &              2.0 * dvdWmj05*presioverrhoi2/EU0(2,i)
c
c--Add thermal conductivity
c
               IF (ifsvi.EQ.7) THEN
                  vsig = cs - 2.0*dvdWmj05/pmjdW
                  IF (vsig.GT.0.0) THEN
                     pres_numerator = pres_numerator -
     &                    0.05*pmjdW*vsig/rhomean*EU0(2,j)
                     pres_denominator = pres_denominator +
     &                    0.05*pmjdW*vsig/rhomean
                  ENDIF
               ENDIF
c
c--Set c*lambda/kappa*rho term for current quantities
c
               bi = ekcle(4,i)/(ekcle(2,i)*rhoi)
               bj = ekcle(4,j)/(ekcle(2,j)*rhoj)
               b1 = (4.0*bi*bj)/(bi+bj)
c	    
c--Diffusion numerator and denominator
c
               diffusion_numerator = diffusion_numerator -
     &              dWdrlightrhorhom*b1*EU0(1,J)*rhoj
               diffusion_denominator = diffusion_denominator +
     &              dWdrlightrhorhom*b1*rhoi
c               if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'la',
c     &              k,j,diffusion_denominator,dWdrlightrhorhom,b1,rhoi,
c     &              ekcle(4,i),ekcle(2,i),ekcle(4,j),ekcle(2,j),rhoj
            END DO              !J-loop
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl 2'
c
c--Radiation pressure...
c
            gradEi2 = (dedxyz(1,i)**2+dedxyz(2,i)**2+dedxyz(3,i)**2)
          
            IF (gradEi2.EQ.0.0) THEN
               gradvPi = 0.0
            ELSE
               rpdiag=0.5*(1.0-ekcle(5,i))
               rpall=0.5*(3.0-ekcle(5,i))/gradEi2
               gradvPi=(((rpdiag+rpall*dedxyz(1,i)**2)*dvdx(1,i))+
     $              ((rpall*dedxyz(1,i)*dedxyz(2,i))*dvdx(2,i))+
     $              ((rpall*dedxyz(1,i)*dedxyz(3,i))*dvdx(3,i))+
     $              ((rpall*dedxyz(2,i)*dedxyz(1,i))*dvdx(4,i))+
     $              ((rpdiag+rpall*dedxyz(2,i)**2)*dvdx(5,i))+
     $              ((rpall*dedxyz(2,i)*dedxyz(3,i))*dvdx(6,i))+
     $              ((rpall*dedxyz(3,i)*dedxyz(1,i))*dvdx(7,i))+
     $              ((rpall*dedxyz(3,i)*dedxyz(2,i))*dvdx(8,i))+
     $              ((rpdiag+rpall*dedxyz(3,i)**2)*dvdx(9,i)))
            ENDIF

            radpresdenom = gradvPi/rhoi * EU0(1,i)

c         radpresdenom=0.0

            tfour=uradconst*lightspeed*ekcle(2,i)* 
     &         ((rhoi*EU0(1,i)/uradconst) - ((EU0(2,i)/ekcle(3,i))**4))
c
c--Now solve those equations...
c
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl 3'
            betaval = lightspeed*ekcle(2,i)*rhoi*dti
            chival = dti*(diffusion_denominator-radpresdenom/EU0(1,I))-
     &         betaval
            gammaval = uradconst*lightspeed*ekcle(2,i)/(ekcle(3,i))**4
            u4term = gammaval*dti*(chival + betaval - 1.0)
            u1term = (chival-1.0)*(1.0-dti*pres_denominator)
            u0term = betaval*origEU(1,i) + (chival-1.0)*(-origEU(2,i)-
     &        dti*pres_numerator) + dti*diffusion_numerator*betaval

c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl3a ',
c     &           u0term,u1term,u4term,gammaval,chival,betaval
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl3b ',
c     &           diffusion_denominator,radpresdenom,EU0(1,I),betaval

            IF (u1term.GT.0.0 .AND. u0term.GT.0.0 .OR. u1term.LT.0.0 
     &           .AND. u0term.LT.0.0) THEN
C$OMP CRITICAL(quart)
            print *,"ngs ",u4term,u1term,u0term,betaval,chival,gammaval
            print *,"    ",ekcle(2,i),rhoi,dti
            print *,"    ",diffusion_denominator,diffusion_numerator
            print *,"    ",pres_denominator,pres_numerator,uradconst
            print *,"    ",radpresdenom,EU0(1,I),EU0(2,I),ekcle(3,i)
            print *,"    ",lightspeed,origEU(1,i),origEU(2,i)
C$OMP END CRITICAL(quart)
C$OMP CRITICAL (moresweepset)
               moresweep=.TRUE.
C$OMP END CRITICAL (moresweepset)
               GOTO 200
            ENDIF

            u1term = u1term/u4term
            u0term = u0term/u4term

            moresweep2 = .FALSE.

c            u1term = u1term/1.0E+12
c            u0term = u0term/1.0E+16

c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl 4',
c     &           u1term,u0term,EU0(2,I),U1i,u4term
c            print *,i,u1term,u0term,EU0(2,I),U1i,chival,pres_denominator
c     &           ,dti,diffusion_denominator,radpresdenom,EU0(1,I),
c     &           betaval
            CALL quartic_gs1t(u1term,u0term,EU0(2,I),U1i,moresweep2,
     &           nlst_end,nlstall)
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl 5'

c            U1i = U1i*1.0E+4

            IF (moresweep2) THEN
C$OMP CRITICAL (moresweepset)
               moresweep=.TRUE.
               PRINT *,"Info: ",EU0(2,I)/ekcle(3,i),xyzmh(1,i),
     &              xyzmh(2,i),xyzmh(3,i)
C$OMP END CRITICAL (moresweepset)               
            ENDIF

            E1i = (origEU(1,i)+dti*diffusion_numerator+gammaval*dti*
     &           U1i**4)/(1.0-chival)

c            IF (EU0(2,I)/ekcle(3,i).GT.5000.0) THEN
c            IF (ABS(((rhoi*EU0(1,i)/uradconst) - ((EU0(2,i)/ekcle(3,i))**4))/
c     &              ((EU0(2,i)/ekcle(3,i))**4)) THEN
c               tfour = 0.0
c               tfour = (origEU(1,i) - EU0(1,I))/dti + (diffusion_numerator + 
c     &              diffusion_denominator*EU0(1,I) - radpresdenom)
c            ENDIF
            dUcomb = pres_numerator + pres_denominator*EU0(2,I) + tfour
            dEcomb = diffusion_numerator + diffusion_denominator*
     &           EU0(1,I) - tfour - radpresdenom
c
c--Tests for negativity
c
            IF(U1i.LE.0.0) THEN
C$OMP CRITICAL (moresweepset)
               moresweep=.TRUE.
               PRINT *,"GSIMPL: U has gone negative ",i,xyzmh(1,i),
     &              xyzmh(2,i),xyzmh(3,i)
C$OMP END CRITICAL (moresweepset)
            ENDIF

            IF(E1i.LE.0.0) THEN
C$OMP CRITICAL (moresweepset)
               moresweep=.TRUE.
               PRINT *,"GSIMPL: E has gone negative ",i,xyzmh(1,i),
     &              xyzmh(2,i),xyzmh(3,i)
C$OMP END CRITICAL (moresweepset)
            END IF
c
c--And the error is...
c
            IF (EU0(2,I)/ekcle(3,i).GT.0.0) THEN
               maxerrE2 = MAX(maxerrE2,1.0*ABS((EU0(1,I) - E1i) /E1i) 
c     &              *10.0
     &              )



               residualE = 0.0
            ELSE
               maxerrE2 = MAX(maxerrE2,ABS((origEU(1,i) + (dEcomb)*
     &              dti - E1i) /E1i))
               residualE = origEU(1,i) + (dEcomb)*dti - E1i
            ENDIF
c            IF (maxerrE2.EQ.ABS((origEU(1,i) + (dEcomb)*
c     &           dti - E1i) /E1i)) THEN
cC$OMP CRITICAL (ipostest)
c               ipos = i
c               PRINT *,"        ",i,sqrt(xyzmh(1,i)**2+xyzmh(2,i)**2+
c     &              xyzmh(3,i)**2),rhoi,U1i/ekcle(3,i),
c     &              (U1i/ekcle(3,i))**4,rhoi*E1i/uradconst,
c     &              (EU0(2,I)/ekcle(3,i))**4,rhoi*EU0(1,I)/uradconst,
c     &              origEU(1,i),E1i,(dEcomb),dti,
c     &              diffusion_numerator,diffusion_denominator,
c     &              tfour,radpresdenom,uradconst,lightspeed,ekcle(2,i)
cC$OMP END CRITICAL (ipostest)
c            ENDIF


c            maxerrE2 = MAX(maxerrE2,ABS((origEU(1,i) + (dEcomb)*
c     &           dti - EU0(1,I)) /EU0(1,I)))

c            maxerrU2 = MAX(maxerrU2,ABS((origEU(2,i)+(dUcomb)* 
c     &           dti - EU0(2,I))/EU0(2,I)))
            IF (EU0(2,I)/ekcle(3,i).GT.0.0) THEN
               maxerrU2 = MAX(maxerrU2,1.0*ABS((EU0(2,I) - U1i) /U1i)
c     &              *10.0
     &              )



               residualU = 0.0
            ELSE
               maxerrU2 = MAX(maxerrU2,ABS((origEU(2,i)+(dUcomb)* 
     &              dti - U1i)/U1i))
               residualU = origEU(2,i)+(dUcomb)*dti - U1i
            ENDIF
c         
c--Copy values
c
            EU0(1,I) = E1i
            EU0(2,I) = U1i
c            EU0(1,I) = E1i - residualE/200000000.0*(1.0+nosweep/100.0)**1
c            EU0(2,I) = U1i - residualU/200000000.0*(1.0+nosweep/100.0)**1

c            cvcold = 1.5*Rg/(gmw*uergg)
c            ucold = boundtemp*cvcold
c            ecold = uradconst*(boundtemp)**4/rhoi
c            IF (EU0(1,I).LT.ecold .OR. U0(I).LT.ucold) THEN
c               EU0(1,I) = ecold
c               EU0(2,I) = ucold
c            ENDIF
c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71) write (*,*) 'pl 6'

            ekcle(3,i) = GETCV(rho(i),EU0(2,i))
            oneovermu(i) = GET1OVERMU(rho(i),EU0(2,i))

c            if (nlst_end.EQ.34 .AND. nlstall.EQ.71 .AND. i.EQ.22) THEN
c            if (i.EQ.22) THEN
c               iflag = 1
c            else
c               iflag = 0
c            endif

            ekcle(2,i) = GETKAPPA(EU0(2,i),ekcle(3,i),rho(i))

 200        CONTINUE
         END DO ! I-loop
C$OMP END DO
C$OMP DO SCHEDULE(static)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            ekcle(2,i) = ekcle(2,j)
         END DO
C$OMP END DO
C$OMP END PARALLEL
c         write (*,*) 'place 4'
         IF (MOD(nosweep,1).EQ.0) THEN
c            PRINT *,"GSIMPL: Finished iteration ",
c     $        nosweep," of ",nswmax," (",maxerrE2,maxerrU2,")"
         ENDIF

         IF (moresweep) RETURN
c
c--The actual test
c
         IF(maxerrE2.LE.tolerance.AND.maxerrU2.LE.tolerance) THEN
c            PRINT *,"Complete with ",nosweep," iterations"
            GOTO 150
         ENDIF
c
c--Test for oscillations
c
         xmaxerrtot = MAX(maxerrE2,maxerrU2)
         IF (xmaxerrtot.GT.xmaxerrold(1) .AND. nosweep.GT.10) THEN
            numoscillations = numoscillations + 1
            IF (numoscillations.GT.100) THEN
               PRINT *,"GSIMPL: Oscillating ",xmaxerrtot,
     &              (xmaxerrold(ii),ii=1,ntests)
               moresweep = .TRUE.
               RETURN
            ENDIF
         ENDIF
c
c--Test for convergence to non-zero value (incorrect minimum)
c     Must have equal value at least twice to stop detecting up and down
c     as non-convergence
c
c         GOTO 333

         DO itest = ntests,1,-1

c            GOTO 332

            IF (nosweep.GT.10) THEN
               IF (xmaxerrtot.GT.0.99999*xmaxerrold(itest).AND.
     &              xmaxerrtot.LT.1.00001*xmaxerrold(itest)) THEN
                  DO iii = 1, numcomp
                     IF (xmaxerrtot.GT.0.99999*xmaxerrcomp(iii).AND.
     &                    xmaxerrtot.LT.1.00001*xmaxerrcomp(iii)) THEN
                        PRINT *,"GSIMPL: Non-convergence ",xmaxerrtot,
     &                    xmaxerrcomp(iii),(xmaxerrold(ii),ii=1,ntests)
                        moresweep = .TRUE.
                        RETURN
                     ENDIF
                  END DO
                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
                  xmaxerrcomp(numcomp) = xmaxerrtot
                  PRINT *,"GSIMPL: Almost Non-convergence ",
     &                 xmaxerrtot,(xmaxerrold(ii),ii=1,ntests)
               ENDIF
            ENDIF
 332        CONTINUE
            IF (itest.NE.1) THEN
               xmaxerrold(itest) = xmaxerrold(itest-1)
            ELSE
               xmaxerrold(itest) = xmaxerrtot
            ENDIF
         END DO
 333     CONTINUE
      END DO ! Iterations loop
c
c--Maximum number of iterations reached
c
      PRINT *,"GSIMPL: Warning. Maximum iterations reached"
      moresweep = .TRUE.
      RETURN
      STOP
c
c--Output success
c
 150  PRINT *,"Succeeded with ",nosweep," iterations ",maxerrE2,maxerrU2
      nit = nosweep
      error = MAX(maxerrE2,maxerrU2)
c
c--And that done, return everything to ASS
c
C$OMP PARALLEL DO SCHEDULE(static) default(none)
C$OMP& shared(nlst_in,nlst_end,list,vxyzu,EU0,iflag)
C$OMP& shared(Rg,gmw,uergg,boundtemp,uradconst,rho,ekcle,oneovermu)
C$OMP& private(n,i)
C$OMP& private(cvcold,ucold,ecold)
      DO n = nlst_in,nlst_end
         i = list(n)
         ekcle(1,i) = EU0(1,i)
         vxyzu(4,i) = EU0(2,i)

         cvcold = 1.5*Rg/(gmw*uergg)
         ucold = 0.5*boundtemp*cvcold
         ecold = uradconst*(0.5*boundtemp)**4/rho(i)
         IF (ekcle(1,i).LT.ecold .OR. vxyzu(4,i).LT.ucold) THEN
c            ekcle(1,i) = ecold
c            vxyzu(4,i) = ucold
c            ekcle(3,i) = GETCV(rho(i),vxyzu(4,i))
c            oneovermu(i) = GET1OVERMU(rho(i),vxyzu(4,i))

c            if (i.EQ.22) THEN
c              iflag = 1
c            else
c               iflag = 0
c            endif

c            ekcle(2,i) = GETKAPPA(vxyzu(4,i),ekcle(3,i),rho(i))
         ENDIF
      END DO
C$OMP END PARALLEL DO

c      write (*,*) 'exit gsimpl',ekcle(2,22)

      RETURN
      END !SUBROUTINE GSIMPL
