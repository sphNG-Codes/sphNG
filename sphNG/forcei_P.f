      SUBROUTINE forcei(nlst_in,nlst_end,list,dt,itime,
     &      npart,xyzmh,vxyzu,fxyzu,dha,
     &      trho,pr,vsound,alphaMMpass,ekcle,dedxyz)
c************************************************************
c                                                           *
c  This subroutine computes the forces on particle ipart    *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,idim), vxyzu(4,idim)
      REAL*4 trho(idim),pr(idim),vsound(idim),dha(1+isizealphaMM,idim)
      REAL*4 alphaMMpass(isizealphaMM,idim)
      DIMENSION list(idim)
      DIMENSION fxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION dedxyz(3,iradtrans)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/vsmooth'

      REAL*4 ddvxyz(3,idim)

      CHARACTER*7 where

      DATA where/'forcei'/
      DATA epsil/1.E-2/
      DATA epsil2/1.E-4/
      DATA epsilvbar/0.5/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine forcei')

      realtime = dt*itime/imaxstep + gt

      IF (nlst_end.GT.nptmass) THEN

      IF (ifsvi.EQ.6 .OR. XSPH) THEN
c
c--Calculates ddvx, ddvy, ddvz only if viscosity=6
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,list,nneigh,iphase,neighover,neighb)
C$OMP& shared(xyzmh,radkernel,trho,grwij,vxyzu,dvtable,ddvtable)
C$OMP& shared(ddvxyz,cnormk)
C$OMP& shared(epsilvbar,vsmooth,wij)
C$OMP& private(n,ipart,ddvxi,ddvyi,ddvzi,k,j,dx,dy,dz,rij2,rij,rij1)
C$OMP& private(hmean,hmean21,hmean41,v,v2,index,dxx,index1,rhoj,grwtij)
C$OMP& private(dgrwdx,grpm,dvx,dvy,dvz,edotv,termx,termy,termz)
C$OMP& private(ddvscalar,vsmoothxi,vsmoothyi,vsmoothzi,hmean1,hmean31)
C$OMP& private(rhoi,dwdx,wtij,wtijmrho,temp)
      DO 180 n = nlst_in, nlst_end
         ipart = list(n)

         vsmooth(1,ipart) = vxyzu(1,ipart)
         vsmooth(2,ipart) = vxyzu(2,ipart)
         vsmooth(3,ipart) = vxyzu(3,ipart)

         IF (iphase(ipart).GE.1) GOTO 180

         rhoi = trho(ipart)

         ddvxi = 0.
         ddvyi = 0.
         ddvzi = 0.

         vsmoothxi = 0.
         vsmoothyi = 0.
         vsmoothzi = 0.
c
c--Loop over neighbors
c
         DO 170 k = 1, nneigh(ipart)
            IF (k.GE.nlmax) THEN
               j = neighover(k-nlmax+1,ABS(neighb(nlmax,ipart)))
            ELSE
               j = neighb(k,ipart)
            ENDIF

            IF (iphase(j).GE.1) GOTO 170
c
c--Gravity and potential energy
c
            dx = xyzmh(1,ipart) - xyzmh(1,j)
            dy = xyzmh(2,ipart) - xyzmh(2,j)
            dz = xyzmh(3,ipart) - xyzmh(3,j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            rij = SQRT(rij2)
            rij1 = 1./rij
c
c--Define mean h
c
            hmean = 0.5*(xyzmh(5,ipart) + xyzmh(5,j))
            hmean1 = 1.0/hmean
            hmean21 = hmean1*hmean1
            hmean31 = hmean21*hmean1
            hmean41 = hmean21*hmean21

            v2 = rij2*hmean21
            v = rij/hmean

            index = v2*ddvtable
            dxx = v2 - index*dvtable
            index1 = index + 1
            IF (index1.GT.itable) index1 = itable

            IF (v.LT.radkernel) THEN
               rhoj = trho(j)
c 
c--Get kernel quantities from interpolation in table
c
               dwdx = (wij(index1) - wij(index))*ddvtable
               wtij = (wij(index) + dwdx*dxx)*hmean31
               dgrwdx = (grwij(index1) - grwij(index))*ddvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
               wtijmrho = xyzmh(4,j)*wtij/(rhoi+rhoj)
               grpm = xyzmh(4,j)*grwtij

               dvx = vxyzu(1,ipart) - vxyzu(1,j)
               dvy = vxyzu(2,ipart) - vxyzu(2,j)
               dvz = vxyzu(3,ipart) - vxyzu(3,j)
c
c--Smoothed velocity
c
               vsmoothxi = vsmoothxi - dvx*wtijmrho
               vsmoothyi = vsmoothyi - dvy*wtijmrho
               vsmoothzi = vsmoothzi - dvz*wtijmrho
c
c--Calculates grad of div velocity (straigh, not times density)
c
               edotv = (dx*dvx + dy*dvy + dz*dvz)/rij
               termx = - grpm/rhoj*(5.0*dx/rij*edotv - dvx)/rij
               termy = - grpm/rhoj*(5.0*dy/rij*edotv - dvy)/rij
               termz = - grpm/rhoj*(5.0*dz/rij*edotv - dvz)/rij
               ddvxi = ddvxi + termx
               ddvyi = ddvyi + termy
               ddvzi = ddvzi + termz
               ddvscalar = (termx*dvx + termy*dvy + termz*dvz)/
     &              SQRT(dvx**2 + dvy**2 + dvz**2)
c               print *,"edotv",edotv,ddvxi,ddvyi,ddvzi
            ENDIF
 170     CONTINUE
c
c--Store quantities
c
         ddvxyz(1,ipart)=cnormk*ddvxi
         ddvxyz(2,ipart)=cnormk*ddvyi
         ddvxyz(3,ipart)=cnormk*ddvzi

         temp = 2.0*cnormk*epsilvbar

         vsmooth(1,ipart) = vsmooth(1,ipart) + temp*vsmoothxi
         vsmooth(2,ipart) = vsmooth(2,ipart) + temp*vsmoothyi
         vsmooth(3,ipart) = vsmooth(3,ipart) + temp*vsmoothzi

 180   CONTINUE
C$OMP END PARALLEL DO
       ENDIF
c
c--Initialize
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,npart,hmin,list,nneigh,neimin,neimax)
C$OMP& shared(icall,dt,imaxstep,isteps)
C$OMP& shared(xyzmh,vxyzu,dha,fxyzu,trho,pr,vsound,dq)
C$OMP& shared(neighb,neighover,dvtable,ddvtable,psoft)
C$OMP& shared(fmass,fpoten,dphidh,part2kernel,part1kernel,radkernel)
C$OMP& shared(part2potenkernel,part1potenkernel,grwij)
C$OMP& shared(divv,curlv,beta,alpha,poten,dgrav,nlst0)
C$OMP& shared(cnormk,epsil,epsil2,where,pext)
C$OMP& shared(iphase,listpm,iprint,nptmass,iorig)
C$OMP& shared(alphaMMpass,alphamin,ddv,ddvxyz)
C$OMP& shared(igrp,igphi,ifsvi,iexf)
C$OMP& shared(ifcor,iexpan,iener,damp)
C$OMP& shared(realtime,ekcle,encal,dedxyz,vsmooth)
C$OMP& private(n,ipart,stepsi,numneigh,vsig)
C$OMP& private(xi,yi,zi,vxi,vyi,vzi,pmassi,dhi,hi,gravxi,gravyi,gravzi)
C$OMP& private(poteni,dphiti,gradxi,gradyi,gradzi,artxi,artyi,artzi)
C$OMP& private(pdvi,dqi,rhoi,pro2i,vsoundi,k,j,hj,dx,dy,dz)
C$OMP& private(rij2,rij,rij1,pmassj,runix,runiy,runiz,hmean,hmean21)
C$OMP& private(hmean41,dhmean,v2,v,index,dxx,index1,rij2grav,rijgrav,fm)
C$OMP& private(phi,dphi,dfmassdx,dfptdx,dpotdh,xmasj,rhoj,robar)
C$OMP& private(dgrwdx,grwtij,grpm,poro2j,dvx,dvy,dvz,projv,vsbar)
C$OMP& private(f,adivi,acurlvi,fi,adivj,acurlvj,fj,t12j)
C$OMP& private(ddvxi,ddvyi,ddvzi,edotv)
C$OMP& private(vlowcorrection,qi,qj)
C$OMP& private(ii,iptcurv,xii,yii,zii,vpos)
C$OMP& private(alphamean,projddv,termx,termy,termz,ddvscalar)
C$OMP& reduction(+:ioutmin,ioutsup,ioutinf)
C$OMP& reduction(MIN:inmin,inminsy)
C$OMP& reduction(MAX:inmax,inmaxsy)

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit
         ELSEIF (iphase(ipart).EQ.0) THEN
c
c--Derivative of smoothing length
c       
            IF (icall.EQ.3) THEN
               numneigh = nneigh(ipart)
               inmin = MIN(inmin,numneigh)
               inmax = MAX(inmax,numneigh)
               inminsy = MIN(inminsy,numneigh)
               inmaxsy = MAX(inmaxsy,numneigh)
               IF (xyzmh(5,ipart).LT.hmin .AND. numneigh.GT.neimin)
     &              ioutmin = ioutmin + 1
               IF (numneigh.GT.neimax) ioutsup = ioutsup + 1
               IF (numneigh.LT.neimin) ioutinf = ioutinf + 1
            ENDIF
            stepsi = dt*isteps(ipart)/imaxstep
            CALL hdot(npart, ipart, stepsi, xyzmh, dha)
         ENDIF
      END DO
C$OMP END DO

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).GE.1) GOTO 80
c
c--Compute forces on particle ipart
c
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         pmassi = xyzmh(4,ipart)
         hi = xyzmh(5,ipart)
         IF (XSPH) THEN
            vxi = vsmooth(1,ipart)
            vyi = vsmooth(2,ipart)
            vzi = vsmooth(3,ipart)
         ELSE
            vxi = vxyzu(1,ipart)
            vyi = vxyzu(2,ipart)
            vzi = vxyzu(3,ipart)
         ENDIF
         dhi = dha(1,ipart)

         gravxi = 0.
         gravyi = 0.
         gravzi = 0.
         poteni = 0.
         dphiti = 0.

         gradxi = 0.
         gradyi = 0.
         gradzi = 0.

         artxi = 0.
         artyi = 0.
         artzi = 0.

         pdvi = 0.
         dqi = 0.

         ddvxi = 0.
         ddvyi = 0.
         ddvzi = 0.
         ddvscalar = 0.

         rhoi = trho(ipart)
         IF (iphase(ipart).GE.1) THEN
            pro2i = 0.0
         ELSE
            pro2i = (pr(ipart) - pext)/(rhoi*rhoi)
         ENDIF
         vsoundi = vsound(ipart)

         stepsi = dt*isteps(ipart)/imaxstep
c
c--Loop over neighbors
c
         DO 70 k = 1, nneigh(ipart)
            IF (k.GE.nlmax) THEN
               j = neighover(k-nlmax+1,ABS(neighb(nlmax,ipart)))
            ELSE
               j = neighb(k,ipart)
            ENDIF

            IF (iphase(j).GE.1) GOTO 70

            IF (iphase(j).EQ.-1) THEN
               WRITE(iprint,*)'ERROR - Accreted particle as neighbour!'
               WRITE(iprint,*) j,iorig(j),xyzmh(1,j),xyzmh(2,j),
     &              xyzmh(3,j),vxyzu(1,j),vxyzu(2,j),vxyzu(3,j)
               WRITE(iprint,*) ipart,iorig(ipart),icall,xi,yi,zi
               CALL quit
            ENDIF
c
c--Gravity and potential energy
c
            dx = xi - xyzmh(1,j)
            dy = yi - xyzmh(2,j)
            dz = zi - xyzmh(3,j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            rij = SQRT(rij2)
            rij1 = 1./rij
            pmassj = xyzmh(4,j)
            hj = xyzmh(5,j)
c
c--Unit vectors
c
            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1
c
c--Define mean h
c
            IF (iphase(ipart).GE.1) THEN
               hmean = hj/2.0
            ELSEIF (iphase(j).GE.1) THEN
               hmean = hi/2.0
            ELSE
               hmean = 0.5*(hi + hj)
            ENDIF
            hmean21 = 1./(hmean*hmean)
            hmean41 = hmean21*hmean21
c
c--dhmean uses old dh(j) if particle j is not currently being evaluated
c
            IF (iphase(ipart).GE.1) THEN
               dhmean = dha(1,j)*stepsi
            ELSE IF (iphase(j).GE.1) THEN
               dhmean = dhi*stepsi
            ELSE
               dhmean = 0.5*(dhi + dha(1,j))*stepsi
            ENDIF
            v2 = rij2*hmean21
            v = rij/hmean

            index = v2*ddvtable
            dxx = v2 - index*dvtable
            index1 = index + 1
            IF (index1.GT.itable) index1 = itable

            IF (igrape.EQ.0 .AND. igphi.NE.0) THEN
               IF (isoft.EQ.1) THEN
                  rij2grav = dx*dx + dy*dy + dz*dz + psoft**2
                  rijgrav = SQRT(rij2grav)
                  fm = 1.0
                  phi = - 1./rijgrav
                  dphi = 0.0
               ELSEIF (isoft.EQ.0) THEN
                  rij2grav = rij2
                  rijgrav = rij
                  IF (v.GE.radkernel) THEN
                     fm = 1.0
                     phi = -rij1
                     dphi = 0.0
                  ELSE
                     dfmassdx = (fmass(index1) - fmass(index))*ddvtable
                     fm = (fmass(index) + dfmassdx*dxx)
                     dfptdx = (fpoten(index1) - fpoten(index))*ddvtable
                     phi = (fpoten(index) + dfptdx*dxx)/hmean
                     dpotdh = (dphidh(index1) - dphidh(index))*ddvtable
                     dphi = (dphidh(index) + dpotdh*dxx)*hmean21*dhmean
                     IF (v.GT.part2kernel) THEN
                        phi = phi + rij1*part2potenkernel
                     ELSEIF (v.GT.part1kernel) THEN
                        phi = phi + rij1*part1potenkernel
                     ENDIF
                  ENDIF
               ELSE
                  CALL error(where,1)
               ENDIF
c
c--Gravitational force calculation
c
               IF (j.LE.npart) THEN
                  xmasj = fm*pmassj/(rij2grav*rijgrav)
                  gravxi = gravxi - xmasj*dx
                  gravyi = gravyi - xmasj*dy
                  gravzi = gravzi - xmasj*dz
                  poteni = poteni + phi*pmassj
                  dphiti = dphiti + pmassj*dphi
               ENDIF
            ENDIF
c
c--Pressure and artificial viscosity
c     There is no pressure between point masses and particles
c     There is no viscosity between point masses and particles
c     There is no viscosity or pressure between two point masses
c
c--No artificial viscosity or pressure between particles across a point mass
c
            DO ii = 1, nptmass
               iptcurv = listpm(ii)
               xii = xyzmh(1,iptcurv)
               yii = xyzmh(2,iptcurv)
               zii = xyzmh(3,iptcurv)
               vpos = (xii-xi)*(xii-xyzmh(1,j)) + 
     &              (yii-yi)*(yii-xyzmh(2,j)) +
     &              (zii-zi)*(zii-xyzmh(3,j))
               IF (vpos.LT.0.0) GOTO 70
            END DO 

            IF(iphase(ipart).GE.1 .OR. iphase(j).GE.1) GOTO 70

            IF (v.LT.radkernel) THEN
               rhoj = trho(j)
               robar = 0.5*(rhoi + rhoj)
c 
c--Get kernel quantities from interpolation in table
c
               dgrwdx = (grwij(index1) - grwij(index))*ddvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
               grpm = pmassj*grwtij
c
c--Pressure gradient and pdv
c
               poro2j = grpm*(pro2i + (pr(j) - pext)/(rhoj**2))
               gradxi = gradxi + poro2j*runix
               gradyi = gradyi + poro2j*runiy
               gradzi = gradzi + poro2j*runiz

               IF (XSPH) THEN
                  dvx = vxi - vsmooth(1,j)
                  dvy = vyi - vsmooth(2,j)
                  dvz = vzi - vsmooth(3,j)
               ELSE
                  dvx = vxi - vxyzu(1,j)
                  dvy = vyi - vxyzu(2,j)
                  dvz = vzi - vxyzu(3,j)
               ENDIF
               projv = dvx*runix + dvy*runiy + dvz*runiz

               pdvi = pdvi + grpm*projv
c
c--Calculates grad of div velocity (straight, not times density)
c
               IF (ifsvi.EQ.6) THEN
                  edotv = (dx*dvx + dy*dvy + dz*dvz)/rij
                  termx = - grpm/rhoj*(5.0*dx/rij*edotv - dvx)/rij
                  termy = - grpm/rhoj*(5.0*dy/rij*edotv - dvy)/rij
                  termz = - grpm/rhoj*(5.0*dz/rij*edotv - dvz)/rij
                  ddvxi = ddvxi + termx
                  ddvyi = ddvyi + termy
                  ddvzi = ddvzi + termz
c                    ddvscalar = (termx*dvx + termy*dvy + termz*dvz)/
cc     &              SQRT(dvx**2 + dvy**2 + dvz**2)
ccc               ddvscalar = ddvscalar+hi*projv*rij/(rij**2+epsil2*hi*hi)
                  projddv=ddvxyz(1,ipart)*runix+ddvxyz(2,ipart)*runiy+
     &                 ddvxyz(3,ipart)*runiz 
             ddvscalar=ddvscalar+hi*projddv*rij/(rij**2+epsil2*hi*hi)
               ENDIF
c
c--Artificial viscosity and energy dissipation
c
               vsbar = 0.5*(vsoundi + vsound(j))

c             IF (ifsvi.NE.0 .AND. projv.LT.0.) THEN
               IF (ifsvi.NE.0 .AND. projv.LT.0. .AND. j.LE.npart) THEN
c
c--Calculate artificial viscosity:
c     If ifsvi=1 then normal viscosity
c     If ifsvi=2 then divv/curl weighted viscosity
c     If ifsvi=3 then viscosity reduced linearly to zero below vsound/2
c     If ifsvi=4 then Balsara viscosity (divv/curl weighted,but times pressure)
c     If ifsvi=5 then viscosity in Hernquist and Katz, ApJS 70, 424
c     If ifsvi=6 then viscosity switch in Morris and Monaghan, 1997, 
c                  J. Comp. Phys. 136, 41-50.  This formulation does not use
c                  beta - it sets beta to be 2*alpha automatically.
c     If ifsvi=7 then use standard viscosity (type 1) but with thermal 
c                 conductivity term as well
c
                  IF (ifsvi.NE.5) THEN
ccc                     f = projv*v/(v2 + epsil)
                     f = projv
                     IF (ifsvi.EQ.2 .OR. ifsvi.EQ.4) THEN
                        adivi = ABS(divv(ipart)/rhoi)
                        acurlvi = ABS(curlv(ipart)/rhoi)
                        fi = adivi/(adivi+acurlvi+epsil2*vsoundi/hi)
                        adivj = ABS(divv(j)/rhoj)
                        acurlvj = ABS(curlv(j)/rhoj)
                        fj = adivj/(adivj+acurlvj+epsil2*vsound(j)/hj)

                        IF (ifsvi.EQ.2) THEN
                           f = f*(fi+fj)/2.0

                           t12j = grpm*f*(beta*f - alpha*vsbar)/robar
                        ELSEIF (ifsvi.EQ.4) THEN
                           f = f*(fi+fj)/(vsoundi+vsound(j))

                           t12j = poro2j*f*(beta*f - alpha)
                        ENDIF
                     ELSEIF (ifsvi.EQ.6) THEN
                        alphamean = (alphaMMpass(1,ipart) + 
     &                       alphaMMpass(1,j))/2.0
                        t12j = alphamean*grpm*f*(2.0*f - vsbar)/robar
                     ELSE
                        t12j = grpm*f*(beta*f - alpha*vsbar)/robar
                     ENDIF
                     IF (ifsvi.EQ.3) THEN
                        vlowcorrection = ABS(projv/vsbar)
                        IF (vlowcorrection.LT.0.5) 
     &                       t12j = 2.0*vlowcorrection*t12j
                     ENDIF
                  ELSE
c
c--Hernquist and Katz
c
                     IF (divv(ipart).LT.0) THEN
                        adivi = ABS(divv(ipart)/rhoi)
                        qi = hi*rhoi*adivi*(alpha*vsoundi + 
     &                       beta*hi*adivi)
                     ELSE
                        qi = 0.0
                     ENDIF
                     IF (divv(j).LT.0) THEN
                        adivj = ABS(divv(j)/rhoj)
                        hj = xyzmh(5,j)
                        qj = hj*rhoj*adivj*(alpha*vsound(j) + 
     &                       beta*hj*adivj)
                     ELSE
                        qj = 0.0
                     ENDIF
                     t12j = grpm*(qi/(rhoi**2) + qj/(rhoj**2))
                  ENDIF

                  artxi = artxi + t12j*runix
                  artyi = artyi + t12j*runiy
                  artzi = artzi + t12j*runiz
                  dqi = dqi + t12j*projv
               ENDIF
c
c--Add thermal conductivity
c
               IF (ifsvi.EQ.7) THEN
                  vsig = vsbar - projv
                  IF (vsig.GT.0.0) THEN
                     dqi = dqi + 0.05*
     &                    grpm*vsig/robar*(vxyzu(4,ipart)-vxyzu(4,j))
                  ENDIF
               ENDIF
            ENDIF
 70      CONTINUE
c
c--Store quantities
c 
        IF (igrape.EQ.0 .AND. igphi.NE.0) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) + gravxi
            fxyzu(2,ipart) = fxyzu(2,ipart) + gravyi
            fxyzu(3,ipart) = fxyzu(3,ipart) + gravzi
            poten(ipart) = poten(ipart) + poteni
c
c--Save correction on gravitational potential
c
            IF (n.LE.nlst0) THEN
               IF (icall.EQ.3 .AND. igphi.NE.0) THEN
                  IF (isoft.EQ.0 .AND. iphase(ipart).EQ.0) 
     &                 dgrav(ipart) = dgrav(ipart) + dphiti
               ENDIF
            ENDIF
         ENDIF
c
c--Pressure gradients
c
         IF (igrp.NE.0) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) - gradxi*cnormk
            fxyzu(2,ipart) = fxyzu(2,ipart) - gradyi*cnormk
            fxyzu(3,ipart) = fxyzu(3,ipart) - gradzi*cnormk
         ENDIF
c
c--Artificial viscosity
c
         IF (ifsvi.NE.0) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) - artxi*cnormk
            fxyzu(2,ipart) = fxyzu(2,ipart) - artyi*cnormk
            fxyzu(3,ipart) = fxyzu(3,ipart) - artzi*cnormk
         ENDIF
         fxyzu(4,ipart) = pdvi
         dq(ipart) = dqi

         IF (ifsvi.EQ.6) THEN
cc         dha(2,ipart) = 0.2*vsoundi*(alphamin-alphaMMpass(1,ipart))/hi-
cc     &           MIN(divv(ipart)/rhoi+0.5*vsoundi/hi,0.0)
c          dha(2,ipart) = 0.2*vsoundi*(alphamin-alphaMMpass(1,ipart))/hi-
c     &           MIN(divv(ipart)/rhoi,0.0)
           dha(2,ipart)= 0.05*vsoundi*(alphamin-alphaMMpass(1,ipart))/hi
            ddv(ipart) = cnormk*SQRT(ddvxi**2 + ddvyi**2 + ddvzi**2)
            ddv(ipart) = ddvscalar
        IF (divv(ipart).LT.0.0.AND.hi*ddv(ipart).LT.
     &           -2.0*vsoundi/hi) THEN
               dha(2,ipart) = dha(2,ipart) - 
     &                  2.0*(hi*ddv(ipart)+2.0*vsoundi/hi)
c     &              1.0*SQRT(ABS(hi*ddv(ipart)*divv(ipart))
            ENDIF
         ENDIF
c
c--Radiation pressure force
c
         IF (encal.EQ.'r') THEN
            DO j = 1, 3               
               fxyzu(j,ipart) = fxyzu(j,ipart) - 
     &              ekcle(4,ipart)/trho(ipart)*dedxyz(j,ipart)
            END DO
         ENDIF
c
c--Damp velocities if appropiate
c
         IF (damp.NE.0.) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) - damp*vxyzu(1,ipart)
            fxyzu(2,ipart) = fxyzu(2,ipart) - damp*vxyzu(2,ipart)
            fxyzu(3,ipart) = fxyzu(3,ipart) - damp*vxyzu(3,ipart)
         ENDIF
c
c--Sink particles jump in here
c
 80      CONTINUE
c
c--Energy conservation
c  
         IF (iphase(ipart).EQ.0) THEN
            CALL energ(ipart,realtime,vxyzu, fxyzu, xyzmh)  
         ELSE
            dha(1,ipart) = 0.0
            fxyzu(4,ipart) = 0.0
         ENDIF
c
c--External forces
c
         IF (iexf.GE.1) CALL externf(ipart,realtime,xyzmh,fxyzu,iexf)
c
c--Coriolis and centrifugal forces
c
         IF (ifcor.NE.0) CALL coriol(ipart,realtime,xyzmh,vxyzu,fxyzu) 
c
c--Homologous expansion or contraction
c
         IF (iexpan.GT.0) CALL homexp(ipart,realtime,vxyzu,fxyzu) 

      END DO
C$OMP END DO
C$OMP END PARALLEL

      ELSE
c
c--Only sink particles being moved
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,list,dha,fxyzu,xyzmh,iexf,ifcor)
C$OMP& shared(iexpan,realtime,vxyzu,vsmooth)
C$OMP& private(n,ipart)
      DO n = nlst_in, nlst_end
         ipart = list(n)
c
c--Set changes in h and energy to zero
c
         dha(1,ipart) = 0.0
         fxyzu(4,ipart) = 0.0

         IF (XSPH) THEN
            vsmooth(1,ipart) = vxyzu(1,ipart)
            vsmooth(2,ipart) = vxyzu(2,ipart)
            vsmooth(3,ipart) = vxyzu(3,ipart)
         ENDIF
c
c--External forces
c
         IF (iexf.GE.1) CALL externf(ipart,realtime,xyzmh,fxyzu,iexf)
c
c--Coriolis and centrifugal forces
c
         IF (ifcor.NE.0) CALL coriol(ipart,realtime,xyzmh,vxyzu,fxyzu) 
c
c--Homologous expansion or contraction
c
         IF (iexpan.GT.0) CALL homexp(ipart,realtime,vxyzu,fxyzu) 
      END DO
C$OMP END PARALLEL DO
      ENDIF
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine forcei')

      RETURN
      END











