      SUBROUTINE forcei(nlst_in,nlst_end,list,dt,itime,
     &      npart,xyzmh,vxyzu,fxyzu,dha,
     &      trho,pr,vsound,alphaMMpass,ekcle,dedxyz,Bxyz,dBxyz)
c************************************************************
c                                                           *
c  This subroutine computes the forces on particle ipart    *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,idim), vxyzu(4,idim)
      REAL*4 trho(idim), pr(idim), vsound(idim), dha(2,idim)
      REAL*4 alphaMMpass(idim)
      DIMENSION list(idim)
      DIMENSION fxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION dedxyz(3,iradtrans)
      DIMENSION Bxyz(3,imhd),dBxyz(3,imhd)

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
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/gradhterms'

      REAL*4 ddvxyz(3,idim)

      CHARACTER*7 where

      DATA where/'forcei'/
      DATA epsil/1.E-2/
      DATA epsil2/1.E-4/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine forcei')

      realtime = dt*itime/imaxstep + gt

      IF (ifsvi.EQ.6) THEN
c
c--Calculates ddvx, ddvy, ddvz only if viscosity=6
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,list,nneigh,iphase,neighover,neighb)
C$OMP& shared(xyzmh,radkernel,trho,grwij,vxyzu,dvtable)
C$OMP& shared(ddvxyz,cnormk)
C$OMP& private(n,ipart,ddvxi,ddvyi,ddvzi,k,j,dx,dy,dz,rij2,rij,rij1)
C$OMP& private(hi,hi21,hi41,vi,v2i,gradhi,index,dxx,index1,rhoj,grwtij)
C$OMP& private(dgrwdx,grpm,dvx,dvy,dvz,edotv,termx,termy,termz)
C$OMP& private(ddvscalar)
      DO 180 n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).GE.1) GOTO 180

         ddvxi = 0.
         ddvyi = 0.
         ddvzi = 0.
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
c--Use hi
c
            hi = xyzmh(5,ipart)
            hi21 = 1./(hi*hi)
            hi41 = hi21*hi21
            gradhi = gradhs(1,ipart)

            v2i = rij2*hi21
            vi = rij/hi

            index = v2i/dvtable
            dxx = v2i - index*dvtable
            index1 = index + 1
            IF (index1.GT.itable) index1 = itable

            IF (vi.LT.radkernel) THEN
               rhoj = trho(j)
c 
c--Get kernel quantities from interpolation in table
c
               dgrwdx = (grwij(index1) - grwij(index))/dvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hi41
               grpm = xyzmh(4,j)*grwtij*gradhi

               dvx = vxyzu(1,ipart) - vxyzu(1,j)
               dvy = vxyzu(2,ipart) - vxyzu(2,j)
               dvz = vxyzu(3,ipart) - vxyzu(3,j)
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

 180   CONTINUE
C$OMP END PARALLEL DO
       ENDIF
c
c--Initialize
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,npart,hmin,list,nneigh,neimin,neimax)
C$OMP& shared(icall,dt,imaxstep,isteps)
C$OMP& shared(xyzmh,vxyzu,dha,fxyzu,trho,pr,vsound,dq,gradhs)
C$OMP& shared(neighb,neighover,dvtable,psoft)
C$OMP& shared(fmass,fpoten,dphidh,part2kernel,part1kernel,radkernel)
C$OMP& shared(part2potenkernel,part1potenkernel,grwij)
C$OMP& shared(divv,curlv,beta,alpha,poten,dgrav,nlst0)
C$OMP& shared(cnormk,epsil,epsil2,where,pext)
C$OMP& shared(iphase,listpm,iprint,nptmass,iorig)
C$OMP& shared(alphaMMpass,alphamin,ddv,ddvxyz)
C$OMP& shared(igrp,igphi,ifsvi,iexf)
C$OMP& shared(ifcor,iexpan,iener,damp)
C$OMP& shared(realtime,ekcle,encal,dedxyz)
C$OMP& shared(Bxyz,dBxyz,Bextx,Bexty,Bextz)
C$OMP& private(n,ipart,stepsi,numneigh)
C$OMP& private(xi,yi,zi,vxi,vyi,vzi,pmassi,dhi,hi,gravxi,gravyi,gravzi)
C$OMP& private(poteni,dphiti,gradxi,gradyi,gradzi,artxi,artyi,artzi)
C$OMP& private(pdvi,dqi,rhoi,pro2i,vsoundi,k,j,hj,dx,dy,dz)
C$OMP& private(rho1i,rho21i,sqrtrho1i,rho1j)
C$OMP& private(rij2,rij,rij1,pmassj,runix,runiy,runiz,hmean,hmean21)
C$OMP& private(hi1,hi21,hi41,v2i,vi)
C$OMP& private(hj1,hj21,hj41,v2j,vj)
C$OMP& private(index,dxx,index1,rij2grav,rijgrav,fm)
C$OMP& private(phi,dphi,dfmassdx,dfptdx,dpotdh,xmasj,rhoj,robar)
C$OMP& private(dgrwdx,grwtij,grpm,poro2j,dvx,dvy,dvz,projv,vsbar)
C$OMP& private(f,adivi,acurlvi,fi,adivj,acurlvj,fj,t12j)
C$OMP& private(ddvxi,ddvyi,ddvzi,edotv)
C$OMP& private(vlowcorrection,qi,qj)
C$OMP& private(ii,iptcurv,xii,yii,zii,vpos)
C$OMP& private(alphamean,projddv,termx,termy,termz,ddvscalar)
C$OMP& private(gradhi,gradsofti,gradpi,gradpj)
C$OMP& private(Bxi,Byi,Bzi,B2i,B2ext)
C$OMP& private(Bxj,Byj,Bzj,B2j,projBi,projBj,dBx,dBy,dBz)
C$OMP& private(dBxideali,dByideali,dBzideali,dBxdissi,dBydissi,dBzdissi)
C$OMP& private(divBi,curlBxi,curlByi,curlBzi)
C$OMP& reduction(+:ioutmin,ioutsup,ioutinf)
C$OMP& reduction(MIN:inmin,inminsy)
C$OMP& reduction(MAX:inmax,inmaxsy)

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)

      END DO
C$OMP END DO

      IF (imhd.EQ.idim) THEN
         B2ext = Bextx**2 + Bexty**2 + Bextz**2
      ELSE
         B2ext = 0.
         B2i = 0.
         B2j = 0.
      ENDIF

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit
         ELSEIF (iphase(ipart).GE.1) THEN
            GOTO 80
         ENDIF  
c
c--Compute forces on particle ipart
c
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         hi = xyzmh(5,ipart)
         hi1 = 1./hi
         hi21 = hi1*hi1
         hi41 = hi21*hi21
         
         vxi = vxyzu(1,ipart)
         vyi = vxyzu(2,ipart)
         vzi = vxyzu(3,ipart)
         gradhi = gradhs(1,ipart)
         gradsofti = gradhs(2,ipart)
         IF (imhd.EQ.idim) THEN
            Bxi = Bxyz(1,ipart)
            Byi = Bxyz(2,ipart)
            Bzi = Bxyz(3,ipart)
            B2i = Bxi**2 + Byi**2 + Bzi**2
	    dBxideali = 0.
	    dByideali = 0.
	    dBzideali = 0.
	    dBxdissi = 0.
	    dBydissi = 0.
	    dBzdissi = 0.
            divBi = 0.
            curlBxi = 0.
            curlByi = 0.
            curlBzi = 0.
         ENDIF

         gravxi = 0.
         gravyi = 0.
         gravzi = 0.
         poteni = 0.
         dsoftxi = 0.
         dsoftyi = 0.
         dsoftzi = 0.

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
         rho1i = 1./rhoi
         sqrtrho1i = SQRT(rho1i)
         rho21i = rho1i*rho1i
         IF (iphase(ipart).GE.1) THEN
            pro2i = 0.0
         ELSE
            pro2i = (pr(ipart) - pext + 0.5*(B2i-B2ext))*rho21i
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
            hj1 = 1./hj
            hj21 = hj1*hj1
            hj41 = hj21*hj21
c
c--Unit vectors
c
            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1
            dvx = vxi - vxyzu(1,j)
            dvy = vyi - vxyzu(2,j)
            dvz = vzi - vxyzu(3,j)
            projv = dvx*runix + dvy*runiy + dvz*runiz

            v2i = rij2*hi21
            vi = rij*hi1
            v2j = rij2*hj21
            vj = rij*hj1
            rhoj = trho(j)
            rho1j = 1./rhoj
            robar1 = 2./(rhoi + rhoj)
            IF (imhd.EQ.idim) THEN
               Bxj = Bxyz(1,j)
               Byj = Bxyz(2,j)
               Bzj = Bxyz(3,j)
	       dBx = Bxi - Bxj
	       dBy = Byi - Byj
	       dBz = Bzi - Bzj
               B2j = Bxj**2 + Byj**2 + Bzj**2
               projBi = Bxi*runix + Byi*runiy + Bzi*runiz
               projBj = Bxj*runix + Byj*runiy + Bzj*runiz
            ENDIF
c
c--Using hi
c
            IF (vi.LT.radkernel) THEN
               index = v2i/dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dxx = v2i - index*dvtable	       
               dgrwdx = (grwij(index1)-grwij(index))/dvtable ! slope
c              (note that kernel gradient is multiplied by gradhi)
	       grkerni = (grwij(index)+ dgrwdx*dxx)*hi41*gradhi
               grpmi = grkerni*pmassj
c
c--i contribution to pressure gradient and pdv
c
               gradpi = grpmi*pro2i
               pdvi = pdvi + grpmi*projv
c
c--Calculates grad of div velocity (straight, not times density)
c
               IF (ifsvi.EQ.6) THEN
                  edotv = (dx*dvx + dy*dvy + dz*dvz)/rij
                  termx = - grpmi/rhoj*(5.0*dx/rij*edotv - dvx)/rij
                  termy = - grpmi/rhoj*(5.0*dy/rij*edotv - dvy)/rij
                  termz = - grpmi/rhoj*(5.0*dz/rij*edotv - dvz)/rij
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
c--i contribution to force softening (including pseudo-pressure term)
c
               IF (isoft.EQ.0 .OR. isoft.EQ.2) THEN
                  dfmassdx = (fmass(index1) - fmass(index))/dvtable
                  fmi = (fmass(index) + dfmassdx*dxx)
                  dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                  phii = (fpoten(index) + dfptdx*dxx)*hi1
                  IF (vi.GT.part2kernel) THEN
                     phii = phii + rij1*part2potenkernel
                  ELSEIF (vi.GT.part1kernel) THEN
                     phii = phii + rij1*part1potenkernel
                  ENDIF
                  dsofttermi = 0.5*grpmi*gradsofti
                  dsoftxi = dsoftxi - dsofttermi*runix
                  dsoftyi = dsoftyi - dsofttermi*runiy
                  dsoftzi = dsoftzi - dsofttermi*runiz
               ENDIF
               
               IF (imhd.EQ.idim) THEN
c		  
c--time derivative terms
c
	          dBxideali = dBxideali - grpmi*dvx*projBi
	          dByideali = dByideali - grpmi*dvy*projBi
	          dBzideali = dBzideali - grpmi*dvz*projBi
c		  
c--compute divB
c
		  projdB = dBx*runix + dBy*runiy + dBz*runiz
	          divBi = divBi - grpmi*projdB
c
c--compute current
c	          
                  curlBxi = curlBxi + grpmi*(dBy*runiz - dBz*runiy)
		  curlByi = curlByi + grpmi*(dBz*runix - dBx*runiz)
		  curlBzi = curlBzi + grpmi*(dBx*runiy - dBy*runix)
               ENDIF
            ELSE
               grkerni = 0.
               gradpi = 0.
               fmi = 1.
               phii = -rij1
            ENDIF
c
c--Using hj
c
            IF (vj.LT.radkernel) THEN
               index = v2j/dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dxx = v2j - index*dvtable
	       dgrwdx = (grwij(index1)-grwij(index))/dvtable ! slope
c              (note that kernel gradient is multiplied by gradhj)
	       grkernj = (grwij(index)+ dgrwdx*dxx)*hj41*gradhs(1,j)
               grpmj = grkernj*pmassj
c
c--j contribution to pressure gradient
c
               poro2j = (pr(j) - pext + 0.5*(B2j-B2ext))*rho1j*rho1j
               gradpj = grpmj*poro2j
c
c--j contribution to force softening (including pseudo-pressure term)
c
               IF (isoft.EQ.0 .OR. isoft.EQ.2) THEN
                  dfmassdx = (fmass(index1) - fmass(index))/dvtable
                  fmj = (fmass(index) + dfmassdx*dxx)
                  dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                  phij = (fpoten(index) + dfptdx*dxx)*hj1
                  IF (vj.GT.part2kernel) THEN
                     phij = phij + rij1*part2potenkernel
                  ELSEIF (vj.GT.part1kernel) THEN
                     phij = phij + rij1*part1potenkernel
                  ENDIF
                  dsofttermj = 0.5*grpmj*gradhs(2,j)
                  dsoftxi = dsoftxi - dsofttermj*runix
                  dsoftyi = dsoftyi - dsofttermj*runiy
                  dsoftzi = dsoftzi - dsofttermj*runiz
               ENDIF
            ELSE
               grkernj = 0.
               gradpj = 0.
               fmj = 1.
               phij = -rij1
            ENDIF

c
c--Artificial viscosity and pressure forces
c
c--No artificial viscosity or pressure between particles across a point mass
c
c            DO ii = 1, nptmass
c               iptcurv = listpm(ii)
c               xii = xyzmh(1,iptcurv)
c               yii = xyzmh(2,iptcurv)
c               zii = xyzmh(3,iptcurv)
c               vpos = (xii-xi)*(xii-xyzmh(1,j)) + 
c     &              (yii-yi)*(yii-xyzmh(2,j)) +
c     &              (zii-zi)*(zii-xyzmh(3,j))
c               IF (vpos.LT.0.0) GOTO 60
c            END DO 

            IF (vi.LT.radkernel .OR. vj.LT.radkernel) THEN
c
c--add pressure term
c
               gradxi = gradxi + (gradpi + gradpj)*runix
               gradyi = gradyi + (gradpi + gradpj)*runiy
               gradzi = gradzi + (gradpi + gradpj)*runiz
c
c--calculate average of kernel gradients
c
	       grwtij = 0.5*(grkerni + grkernj)
               grpm = pmassj*grwtij

               IF (imhd.EQ.idim) THEN
c
c--anisotropic magnetic force (Morris formalism)
c
                  rhoij1 = rho1i*rho1j
                  projBext = Bextx*runix + Bexty*runiy + Bextz*runiz
		  fanisoxi = fanisoxi
     &                     + grpm*((Bxj-Bextx)*(projBj-projBext)
     &                           - (Bxi-Bextx)*(projBi-projBext))*rhoij1
		  fanisoyi = fanisoyi 
     &                     + grpm*((Byj-Bexty)*(projBj-projBext)
     &                           - (Byi-Bexty)*(projBi-projBext))*rhoij1
		  fanisozi = fanisozi
     &                     + grpm*((Bzj-Bextz)*(projBj-projBext)
     &                           - (Bzi-Bextz)*(projBi-projBext))*rhoij1
c
c--signal velocity (MHD)
c
                  vsoundj = vsound(j)
		  vs2i = vsoundi**2 + B2i*rho1i
		  vs2j = vsoundj**2 + B2j*rho1j
		  vsproji = 2.*vsoundi*projBi*sqrtrho1i
		  vsprojj = 2.*vsoundj*projBj*SQRT(rho1j)		     
		  vsigi = 0.5*(SQRT(vs2i - vsproji)
     &                        +SQRT(vs2i + vsproji))
 		  vsigj = 0.5*(SQRT(vs2j - vsprojj)
     &                        +SQRT(vs2j + vsprojj))
                  vsbar = 0.5*(vsigi + vsigj)
c
c--artificial resistivity
c
		  IF (j.LE.npart) THEN
                     alphaB = 1.0
                     termB = alphaB*grpm*(vsbar + abs(projv))*robar1
c                    dBxdissi = dBxdissi + termB*(dBx - runix*projdB)*robar1
c                    dBydissi = dBydissi + termB*(dBy - runiy*projdB)*robar1
c                    dBzdissi = dBzdissi + termB*(dBz - runiz*projdB)*robar1
		     dBxdissi = dBxdissi + termB*dBx*robar1
		     dBydissi = dBydissi + termB*dBy*robar1
		     dBzdissi = dBzdissi + termB*dBz*robar1

		     !!--add contribution to thermal energy
		     dB2 = dBx*dBx + dBy*dBy + dBz*dBz
c                    dqi = dqi - termB*(dB2 - projdB**2)*robar1
                     dqi = dqi - termB*dB2*robar1
                  ENDIF
               ELSE
c
c--signal velocity (hydro)
c
                  vsbar = 0.5*(vsoundi + vsound(j))
               ENDIF
c
c--Artificial viscosity and energy dissipation
c
               IF (ifsvi.NE.0 .AND. projv.LT.0. .AND. j.LE.npart) THEN
c
c--Calculate artificial viscosity:
c     If ifsvi=1 then normal viscosity
c     If ifsvi=2 then divv/curl weighted viscosity
c     If ifsvi=6 then viscosity switch in Morris and Monaghan, 1997, 
c                  J. Comp. Phys. 136, 41-50.  This formulation does not use
c                  beta - it sets beta to be 2*alpha automatically.
c
                  f = projv
                  IF (ifsvi.EQ.2) THEN
                     adivi = ABS(divv(ipart)*rho1i)
                     acurlvi = ABS(curlv(ipart)*rho1i)
                     fi = adivi/(adivi+acurlvi+epsil2*vsoundi/hi)
                     adivj = ABS(divv(j)*rho1j)
                     acurlvj = ABS(curlv(j)*rho1j)
                     fj = adivj/(adivj+acurlvj+epsil2*vsound(j)/hj)
                     f = f*(fi+fj)/2.0
                     t12j = grpm*f*(beta*f - alpha*vsbar)*robar1
                  ELSEIF (ifsvi.EQ.6) THEN
                     alphamean = (alphaMMpass(ipart) + 
     &                       alphaMMpass(j))/2.0
                     t12j = alphamean*grpm*f*(2.0*f - vsbar)*robar1
                  ELSE
                     t12j = grpm*f*(beta*f - alpha*vsbar)*robar1
                  ENDIF

                  artxi = artxi + t12j*runix
                  artyi = artyi + t12j*runiy
                  artzi = artzi + t12j*runiz
                  dqi = dqi + t12j*projv
               ENDIF
            ENDIF
            
60          CONTINUE
c
c--gravitational force softening
c            
            IF (igrape.EQ.0 .AND. igphi.NE.0) THEN
c--Plummer
               IF (isoft.EQ.1) THEN
                  rij2grav = dx*dx + dy*dy + dz*dz + psoft**2
                  rijgrav = SQRT(rij2grav)
                  fm = 1.0
                  phi = - 1./rijgrav
c--Average softening kernel
               ELSEIF (isoft.EQ.0) THEN
                  rij2grav = rij2
                  rijgrav = rij
                  fm = 0.5*(fmi + fmj)
                  phi = 0.5*(phii + phij)
c--Average softening lengths
               ELSEIF (isoft.EQ.2) THEN
                  rij2grav = rij2
                  rijgrav = rij
                  IF (iphase(ipart).GE.1) THEN
                     hmean = hj/2.0
                  ELSEIF (iphase(j).GE.1) THEN
                     hmean = hi/2.0
                  ELSE
                     hmean = 0.5*(hi + hj)
                  ENDIF
                  hmean21 = 1./(hmean*hmean)
                  v2 = rij2*hmean21
                  v = SQRT(v2)
                  IF (v2.GE.radkernel*radkernel) THEN
                     fm = 1.0
                     phi = -rij1
                     dphi = 0.0
                  ELSE
                     index = v2/dvtable
                     index1 = index + 1
                     IF (index.GT.itable) index = itable
                     IF (index1.GT.itable) index1 = itable
                     dxx = v2 - index*dvtable
                     dfmassdx = (fmass(index1) - fmass(index))/dvtable
                     fm = (fmass(index) + dfmassdx*dxx)
                     dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                     phi = (fpoten(index) + dfptdx*dxx)/hmean
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
               ENDIF
            ENDIF

 70      CONTINUE
c
c--Store quantities
c 
         IF (igrape.EQ.0 .AND. igphi.NE.0) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) + gravxi + dsoftxi*cnormk
            fxyzu(2,ipart) = fxyzu(2,ipart) + gravyi + dsoftyi*cnormk
            fxyzu(3,ipart) = fxyzu(3,ipart) + gravzi + dsoftzi*cnormk
            poten(ipart) = poten(ipart) + poteni
            dgrav(ipart) = 0.
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
c--Anisotropic magnetic force, time derivative of B/rho, div/curl B
c
         IF (imhd.EQ.idim) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) + fanisoxi*cnormk
            fxyzu(2,ipart) = fxyzu(2,ipart) + fanisoyi*cnormk
            fxyzu(3,ipart) = fxyzu(3,ipart) + fanisozi*cnormk
            dBxyz(1,ipart) = cnormk*(dBxideali*rho21i + dBxdissi)
            dBxyz(2,ipart) = cnormk*(dByideali*rho21i + dBydissi)
            dBxyz(3,ipart) = cnormk*(dBzideali*rho21i + dBzdissi)
            divcurlB(1,ipart) = cnormk*divBi*rho1i
            divcurlB(2,ipart) = cnormk*curlBxi*rho1i
            divcurlB(3,ipart) = cnormk*curlByi*rho1i
            divcurlB(4,ipart) = cnormk*curlBzi*rho1i
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
cc         dha(2,ipart) = 0.2*vsoundi*(alphamin-alphaMMpass(ipart))/hi-
cc     &           MIN(divv(ipart)/rhoi+0.5*vsoundi/hi,0.0)
c          dha(2,ipart) = 0.2*vsoundi*(alphamin-alphaMMpass(ipart))/hi-
c     &           MIN(divv(ipart)/rhoi,0.0)
           dha(2,ipart) = 0.05*vsoundi*(alphamin-alphaMMpass(ipart))/hi
            ddv(ipart) = cnormk*SQRT(ddvxi**2 + ddvyi**2 + ddvzi**2)
            ddv(ipart) = ddvscalar
        IF (divv(ipart).LT.0.0.AND.hi*ddv(ipart).LT.
     &           -2.0*vsoundi/hi) THEN
               dha(2,ipart) = dha(2,ipart) - 
     &                  2.0*(hi*ddv(ipart)+2.0*vsoundi/hi)
c     &              1.0*SQRT(ABS(hi*ddv(ipart)*divv(ipart))
            ENDIF
         ENDIF

 80      CONTINUE
c
c--Energy conservation
c  
         IF (iphase(ipart).EQ.0) THEN
            CALL energ(ipart,realtime,vxyzu, fxyzu)  
         ELSE
            dha(1,ipart) = 0.0
            fxyzu(4,ipart) = 0.0
         ENDIF
c
c--External forces
c
         IF (iexf.GE.1) CALL externf(ipart,xyzmh,fxyzu,iexf)
c
c--Coriolis and centrifugal forces
c
         IF (ifcor.NE.0) CALL coriol(ipart,realtime,xyzmh,vxyzu,fxyzu) 
c
c--Homologous expansion or contraction
c
         IF (iexpan.GT.0) CALL homexp(ipart,realtime,vxyzu,fxyzu) 
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

c      IF (ipart.EQ.1) THEN
c         write(iprint,99203) vxyzu(4,ipart)
c99203    FORMAT('pdv ',1F15.10)
c         write(iprint,99204) dq(ipart)
c99204    FORMAT('dq ',1F15.10)
c         write(iprint,99207) fxyzu(1,ipart), fxyzu(2,ipart), 
c     &           fxyzu(3,ipart)
c99207    FORMAT('fx,fy,fz ',1F15.10,1F15.10,1F15.10)
c      END IF

      END DO
C$OMP END DO
C$OMP END PARALLEL
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine forcei')

      RETURN
      END











