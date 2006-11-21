      SUBROUTINE forcei(nlst_in,nlst_end,listp,dt,itime,
     &      npart,xyzmh,vxyzu,fxyzu,dha,trho,pr,vsound,alphaMMpass,
     &      ekcle,dedxyz,Bxyz,dBxyz,Bevolxyz)
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
      DIMENSION listp(idim)
      DIMENSION fxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION dedxyz(3,iradtrans)
      DIMENSION Bxyz(3,imhd),dBxyz(3,imhd)
      DIMENSION Bevolxyz(3,imhd) ! needed for prediction only
c--this weight is equivalent to m/(rho*h^3) in the grad h version
      PARAMETER (hfact = 1.2)
      PARAMETER (weight = 1./hfact**3)

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
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/vsmooth'
c
c--Used for listparents list to keep a list of iupdated
c
      INCLUDE 'COMMONS/treecom_P'

      CHARACTER*7 where

      DATA where/'forcei'/
      DATA epsil/1.E-2/
      DATA epsil2/1.E-4/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine forcei')

      nlistupdated = 0
      realtime = dt*itime/imaxstep + gt

      IF (imhd.EQ.idim .AND. ibound.EQ.7) THEN
c--external magnetic pressure boundaries -- use dummy variables 
c  as can still have external fields even if ibound.ne.7
         Bextxi = 0 !Bextx
         Bextyi = 0 !Bexty
         Bextzi = 0 !Bextz
         B2ext = Bextx**2 + Bexty**2 + Bextz**2
      ELSE
         Bextxi = 0.
         Bextyi = 0.
         Bextzi = 0.
         B2ext = 0.
      ENDIF
c
c--also subtract maximum stress
c      
c      stressmax = 0.
c      IF (imhd.EQ.idim) THEN
c         DO i=1,npart
c            IF (iphase(i).EQ.0) THEN
c               B2i = Bxyz(1,i)**2 + Bxyz(2,i)**2 + Bxyz(3,i)**2
c               stressterm = max(0.5*B2i - pr(i),0.)
c               stressmax = max(stressterm,stressmax)
cc               print*,0.5*B2i-pr(i),stressterm,stressmax
c            ENDIF
c         ENDDO
c         print*,'stressmax = ',stressmax
c      ENDIF
c
c--Initialize
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,npart,listp,nneigh,neimin,neimax)
C$OMP& shared(icall,dt,itime,imaxstep,isteps)
C$OMP& shared(xyzmh,vxyzu,dha,fxyzu,trho,pr,vsound,dq,gradhs)
C$OMP& shared(neighb,neighover,dvtable,ddvtable,psoft)
C$OMP& shared(fmass,fpoten,part2kernel,part1kernel,radkernel)
C$OMP& shared(part2potenkernel,part1potenkernel,grwij)
C$OMP& shared(divv,curlv,beta,alpha,poten,dgrav,nlst0)
C$OMP& shared(cnormk,epsil,epsil2,where,pext)
C$OMP& shared(iphase,listpm,iprint,nptmass,iorig)
C$OMP& shared(alphaMMpass,alphamin)
C$OMP& shared(igrp,igphi,ifsvi,iexf)
C$OMP& shared(ifcor,iexpan,iener,damp)
C$OMP& shared(realtime,ekcle,encal,dedxyz,acc)
C$OMP& shared(Bxyz,dBxyz,Bextxi,Bextyi,Bextzi,B2ext,divcurlB,Bevolxyz)
C$OMP& shared(gravxyzstore,potenstore,varmhd,iupdated,listparents)
C$OMP& shared(iscurrent,ireal,ibound,nlistupdated)
C$OMP& private(n,ipart,stepsi)
C$OMP& private(xi,yi,zi,vxi,vyi,vzi,pmassi,dhi,hi,gravxi,gravyi,gravzi)
C$OMP& private(fxi,fyi,fzi,numneigh,hmin)
C$OMP& private(poteni,dphiti,gradxi,gradyi,gradzi,artxi,artyi,artzi)
C$OMP& private(pdvi,dqi,rhoi,pro2i,pro2j,vsoundi,k,j,hj,dx,dy,dz)
C$OMP& private(rho1i,rho21i,sqrtrho1i,rho1j)
C$OMP& private(rij2,rij,rij1,pmassj,runix,runiy,runiz,hmean,hmean21)
C$OMP& private(hi1,hi21,hi41,v2i,vi)
C$OMP& private(hj1,hj21,hj41,v2j,vj)
C$OMP& private(index,dxx,index1,rij2grav,rijgrav,fm)
C$OMP& private(phi,dphi,dfmassdx,dfptdx,dpotdh,xmasj,rhoj,robar)
C$OMP& private(dgrwdx,grwtij,grpm,dvx,dvy,dvz,projv,vsbar)
C$OMP& private(f,adivi,acurlvi,fi,adivj,acurlvj,fj,t12j)
C$OMP& private(ddvxi,ddvyi,ddvzi,edotv)
C$OMP& private(vlowcorrection,qi,qj)
C$OMP& private(tdecay1,source)
C$OMP& private(ii,iptcurv,xii,yii,zii,vpos)
C$OMP& private(alphamean,projddv,termx,termy,termz,ddvscalar)
C$OMP& private(gradhi,gradsofti,gradpi,gradpj)
C$OMP& private(Bxi,Byi,Bzi,B2i,Bevolxi,Bevolyi,Bevolzi)
C$OMP& private(Bxj,Byj,Bzj,B2j,projBi,projBj,dBx,dBy,dBz)
C$OMP& private(dBxideali,dByideali,dBzideali,dBxdissi,dBydissi,dBzdissi)
C$OMP& private(divBi,curlBxi,curlByi,curlBzi,fanisoxi,fanisoyi,fanisozi)
C$OMP& private(dsoftxi,dsoftyi,dsoftzi,vsigi,vsigj,alphaB)
C$OMP& private(vsproji,vsprojj,termb,vs2i,vs2j,dB2,robar1,fmi)
C$OMP& private(grkerni,grpmi,phii,dsofttermi,projdB,projbext)
C$OMP& private(grkernj,grpmj,phij,dsofttermj,rhoij1,fmj,vsoundj)
C$OMP& private(grkerntable)
C$OMP& private(irealj)
C$OMP& reduction(+:ioutmin,ioutsup,ioutinf)
C$OMP& reduction(MIN:inmin,inminsy)
C$OMP& reduction(MAX:inmax,inmaxsy)

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = listp(n)
c
c--Zero forces and change in thermal energy, and 
c     rate of change of h (the latter so that integrating h in step does
c     not matter).
c
         IF (nlst_end.GT.nptmass) THEN
c
c--Get neighbour lists and calculate gravity on the tree
c     number of neighbours is returned in nneigh() as usual
c     list of neighbours is returned in neighlist(nneigh())
c     which is a threadprivate variable (one copy for each thread)
c     NOTE: Assumes that nlmax=1 is set to save memory for grad-h code.
c     If not, will still work but cache re-use will be VERY BAD !!!!
c     
            IF (iphase(ipart).EQ.0 .OR. iptintree.GE.1) THEN
              CALL treef(ipart,npart,xyzmh,acc,igphi,fxi,fyi,fzi,poteni)

               fxyzu(1,ipart) = fxi
               fxyzu(2,ipart) = fyi
               fxyzu(3,ipart) = fzi
               poten(ipart) = poteni

               gravxyzstore(1,ipart) = fxi
               gravxyzstore(2,ipart) = fyi
               gravxyzstore(3,ipart) = fzi
               potenstore(ipart) = poteni
            ENDIF
         ELSE
c
c--Don't bother to update gravity from gas particles acting on sinks
c    Also means potential energy and neighbours of sinks are not updated
c
            fxyzu(1,ipart) = gravxyzstore(1,ipart)
            fxyzu(2,ipart) = gravxyzstore(2,ipart)
            fxyzu(3,ipart) = gravxyzstore(3,ipart)
            poten(ipart) = potenstore(ipart)
         ENDIF
            
         fxyzu(4,ipart) = 0.
         dha(1,ipart) = 0.0

         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit
         ELSEIF (iphase(ipart).GE.1) THEN
            GOTO 80
         ENDIF

         IF (icall.EQ.3) THEN
            numneigh = nneigh(ipart)
            inmin = MIN(inmin,numneigh)
            inmax = MAX(inmax,numneigh)
            inminsy = MIN(inminsy,numneigh)
            inmaxsy = MAX(inmaxsy,numneigh)
            IF (xyzmh(5,ipart).LT.hmin .AND. numneigh.GT.neimin)
     &           ioutmin = ioutmin + 1
            IF (numneigh.GT.neimax) ioutsup = ioutsup + 1
            IF (numneigh.LT.neimin) ioutinf = ioutinf + 1
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
            Bevolxi = Bevolxyz(1,ipart)
            Bevolyi = Bevolxyz(2,ipart)
            Bevolzi = Bevolxyz(3,ipart)
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
            fanisoxi = 0.
            fanisoyi = 0.
            fanisozi = 0.
         ELSE
            B2i = 0.
            B2j = 0.
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
c--note that pressure term includes isotropic magnetic pressure         
c         pro2i = (pr(ipart) - pext)*rho21i
         pro2i = max((pr(ipart) - pext + 0.5*(B2i-B2ext)),0.0)*rho21i
         vsoundi = vsound(ipart)
	 vs2i = vsoundi**2 + B2i*rho1i

         stepsi = dt*isteps(ipart)/imaxstep
c
c--Loop over neighbors
c
         DO 70 k = 1, nneigh(ipart)
            j = neighlist(k)

            IF (iphase(j).GE.1) GOTO 70

            IF (iphase(j).EQ.-1) THEN
               WRITE(iprint,*)'ERROR - Accreted particle as neighbour!'
               WRITE(iprint,*) j,iorig(j),xyzmh(1,j),xyzmh(2,j),
     &              xyzmh(3,j),vxyzu(1,j),vxyzu(2,j),vxyzu(3,j)
               WRITE(iprint,*) ipart,iorig(ipart),icall,xi,yi,zi
               CALL quit
            ENDIF
c
c--Interpolate density, pressure, vsound and B of neighbour that is not 
c    being evolved. This may have already been done in radiative 
c    transfer routines (flag iupdated set if so).
c
            IF (.NOT.iscurrent(j) .AND. .NOT.iupdated(j)) THEN
               iupdated(j) = .TRUE.
C$OMP FLUSH(iupdated)
C$OMP CRITICAL (listupdated)
               nlistupdated = nlistupdated + 1
               IF (nlistupdated.GT.idim) THEN
                  WRITE (iprint,*) 'ERROR - nlistupdate.GT.idim'
                  CALL quit
               ENDIF
               listparents(nlistupdated) = j
C$OMP END CRITICAL (listupdated)
               IF (j.LE.npart) THEN
                  CALL extrapolate(j,dt,itime,trho,vxyzu,pr,vsound,
     &                             ekcle,Bxyz,Bevolxyz)
               ELSEIF (iupdated(ireal(j))) THEN
                  iupdated(j) = .TRUE.
C$OMP FLUSH(iupdated)
C$OMP CRITICAL (listupdated)
                  nlistupdated = nlistupdated + 1
                  IF (nlistupdated.GT.idim) THEN
                     WRITE (iprint,*) 'ERROR - nlistupdate.GT.idim'
                     CALL quit
                  ENDIF
                  listparents(nlistupdated) = j
C$OMP END CRITICAL (listupdated)
                  irealj = ireal(j)
                  trho(j) = trho(irealj)
                  pr(j) = pr(irealj)
                  vsound(j) = vsound(irealj)
                  divv(j) = 0.
                  IF (encal.EQ.'r' .AND. ibound.EQ.100)
     &                ekcle(1,j) = ekcle(1,irealj)
                  IF (imhd.EQ.idim) THEN
                     IF (varmhd.EQ.'Brho') THEN
                        Bxyz(1,j) = Bxyz(1,irealj)
                        Bxyz(2,j) = Bxyz(2,irealj)
                        Bxyz(3,j) = Bxyz(3,irealj)
                     ENDIF
                  ENDIF
               ELSE
                  irealj = ireal(j)
                  iupdated(irealj) = .TRUE.
                  iupdated(j) = .TRUE.
C$OMP FLUSH(iupdated)
C$OMP CRITICAL (listupdated)
                  nlistupdated = nlistupdated + 2
                  IF (nlistupdated.GT.idim) THEN
                     WRITE (iprint,*) 'ERROR - nlistupdate.GT.idim'
                     CALL quit
                  ENDIF
                  listparents(nlistupdated-1) = j
                  listparents(nlistupdated) = ireal(j)
C$OMP END CRITICAL (listupdated)
                  CALL extrapolate(irealj,dt,itime,trho,vxyzu,pr,vsound,
     &                             ekcle,Bxyz,Bevolxyz)
                  
                  trho(j) = trho(irealj)
                  pr(j) = pr(irealj)
                  vsound(j) = vsound(irealj)
                  divv(j) = 0.
                  IF (encal.EQ.'r' .AND. ibound.EQ.100)
     &                ekcle(1,j) = ekcle(1,irealj)
                  IF (imhd.EQ.idim) THEN
                     IF (varmhd.EQ.'Brho') THEN
                        Bxyz(1,j) = Bxyz(1,irealj)
                        Bxyz(2,j) = Bxyz(2,irealj)
                        Bxyz(3,j) = Bxyz(3,irealj)
                     ENDIF
                  ENDIF

               ENDIF
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
c               projBsi = Bsmooth(1,ipart)*runix + Bsmooth(2,ipart)*runiy
c     &                + Bsmooth(3,ipart)*runiz
c               projBsj = Bsmooth(1,j)*runix + Bsmooth(2,j)*runiy
c     &                + Bsmooth(3,j)*runiz
               projBi = Bxi*runix + Byi*runiy + Bzi*runiz
               projBj = Bxj*runix + Byj*runiy + Bzj*runiz
	       projdB = dBx*runix + dBy*runiy + dBz*runiz
               
               projBext = Bextxi*runix + Bextyi*runiy + Bextzi*runiz
c	       projBrho2i = (projBi - projBext)*rho21i
c	       projBrho2j = (projBj - projBext)*rho1j*rho1j

            ENDIF
c
c--Using hi
c
            IF (vi.LT.radkernel) THEN
               index = v2i*ddvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dxx = v2i - index*dvtable	       
               dgrwdx = (grwij(index1)-grwij(index))*ddvtable ! slope
c              (note that kernel gradient is multiplied by gradhi)
	       grkerntable = (grwij(index)+ dgrwdx*dxx)
               grkerni = grkerntable*hi41*gradhi
               grpmi = grkerni*pmassj
c
c--i contribution to pressure gradient and pdv
c
               gradpi = grpmi*pro2i
               pdvi = pdvi + grpmi*projv
c
c--i contribution to force softening (including pseudo-pressure term)
c
               IF (isoft.EQ.0) THEN
                  dfmassdx = (fmass(index1) - fmass(index))*ddvtable
                  fmi = (fmass(index) + dfmassdx*dxx)
                  dfptdx = (fpoten(index1) - fpoten(index))*ddvtable
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
c	          dBxideali = dBxideali 
c     &                - grpmi*(vsmooth(1,ipart)-vsmooth(1,j))*projBsi
c	          dByideali = dByideali
c     &                - grpmi*(vsmooth(2,ipart)-vsmooth(2,j))*projBsi
c	          dBzideali = dBzideali
c     &                - grpmi*(vsmooth(3,ipart)-vsmooth(3,j))*projBsi
c		  
c--compute divB
c
c	          divBi = divBi - grpmi*projdB
	          divBi = divBi - weight*grkerntable*hi1*projdB
c
c--compute current
c	          
                  curlBxi = curlBxi + grpmi*(dBy*runiz - dBz*runiy)
		  curlByi = curlByi + grpmi*(dBz*runix - dBx*runiz)
		  curlBzi = curlBzi + grpmi*(dBx*runiy - dBy*runix)
c                  curlBxi = curlBxi + weight*grkerntable*hi1*
c     &                               (dBy*runiz - dBz*runiy)
c		  curlByi = curlByi + weight*grkerntable*hi1*
c     &                               (dBz*runix - dBx*runiz)
c		  curlBzi = curlBzi + weight*grkerntable*hi1*
c     &                               (dBx*runiy - dBy*runix)
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
               index = v2j*ddvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dxx = v2j - index*dvtable
	       dgrwdx = (grwij(index1)-grwij(index))*ddvtable ! slope
c              (note that kernel gradient is multiplied by gradhj)
	       grkernj = (grwij(index)+ dgrwdx*dxx)*hj41*gradhs(1,j)
               grpmj = grkernj*pmassj
c
c--j contribution to pressure gradient and isotropic mag force
c
c               pro2j = (pr(j) - pext)*rho1j*rho1j
               pro2j = max((pr(j) - pext + 0.5*(B2j-B2ext)),0.0)
     &                     *rho1j*rho1j
               gradpj = grpmj*pro2j
c
c--j contribution to force softening (including pseudo-pressure term)
c
               IF (isoft.EQ.0 .OR. isoft.EQ.2) THEN
                  dfmassdx = (fmass(index1) - fmass(index))*ddvtable
                  fmj = (fmass(index) + dfmassdx*dxx)
                  dfptdx = (fpoten(index1) - fpoten(index))*ddvtable
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
c--No artificial viscosity, resistivity, pressure, or MHD between particles 
c     across a point mass
c
            DO ii = 1, nptmass
               iptcurv = listpm(ii)
               xii = xyzmh(1,iptcurv)
               yii = xyzmh(2,iptcurv)
               zii = xyzmh(3,iptcurv)
               vpos = (xii-xi)*(xii-xyzmh(1,j)) + 
     &              (yii-yi)*(yii-xyzmh(2,j)) +
     &              (zii-zi)*(zii-xyzmh(3,j))
               IF (vpos.LT.0.0) GOTO 60
            END DO

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
                  projBext = Bextxi*runix + Bextyi*runiy + Bextzi*runiz
		  fanisoxi = fanisoxi
     &                   + grpm*((Bxj-Bextxi)*(projBj-projBext)
     &                         - (Bxi-Bextxi)*(projBi-projBext))*rhoij1
		  fanisoyi = fanisoyi 
     &                   + grpm*((Byj-Bextyi)*(projBj-projBext)
     &                         - (Byi-Bextyi)*(projBi-projBext))*rhoij1
		  fanisozi = fanisozi
     &                   + grpm*((Bzj-Bextzi)*(projBj-projBext)
     &                         - (Bzi-Bextzi)*(projBi-projBext))*rhoij1
c
c--exactly momentum-conserving form (no external boundaries yet)
c
c                  sxxi = -Bxi*Bxi + stressmax
c                  sxyi = -Bxi*Byi + stressmax
c                  sxzi = -Bxi*Bzi + stressmax
c                  syyi = -Byi*Byi + stressmax
c                  syzi = -Byi*Bzi + stressmax
c                  szzi = -Bzi*Bzi + stressmax
c                  
c                  sxxj = -Bxj*Bxj + stressmax
c                  sxyj = -Bxj*Byj + stressmax
c                  sxzj = -Bxj*Bzj + stressmax
c                  syyj = -Byj*Byj + stressmax
c                  syzj = -Byj*Bzj + stressmax
c                  szzj = -Bzj*Bzj + stressmax
c
c                  fanisoxi = pmassj*(
c                  (sxxi*runix + sxyi*runiy + sxzi*runiz)*rho21i*grkerni
c                  (sxxj*runix + sxyj*runiy + sxzj*runiz)*rho21j*grkernj)
c                 
c                  fanisoyi = pmassj*(
c                  (sxyi*runix + syyi*runiy + syzi*runiz)*rho21i*grkerni
c                  (sxyj*runix + syyj*runiy + syzj*runiz)*rho21j*grkernj)
c
c                  fanisozi = pmassj*(
c                  (sxzi*runix + syzi*runiy + szzi*runiz)*rho21i*grkerni
c                  (sxzj*runix + syzj*runiy + szzj*runiz)*rho21j*grkernj)
c
c--signal velocity (MHD)
c
                  vsoundj = vsound(j)
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
		  IF (j.LE.npart .OR. ibound.EQ.11) THEN
                     alphaB = 0.5*(alphaMMpass(2,ipart)
     &                           + alphaMMpass(2,j))
c                     alphaB = 1.0
c                     termB = alphaB*grpm*MAX(vsbar - projv,0.0)*robar1
                     termB = alphaB*grpm*(vsbar + 2.0*abs(projv))*robar1
c                    dBxdissi = dBxdissi + termB*(dBx - runix*projdB)*robar1
c                    dBydissi = dBydissi + termB*(dBy - runiy*projdB)*robar1
c                    dBzdissi = dBzdissi + termB*(dBz - runiz*projdB)*robar1

                     IF (varmhd.EQ.'eulr') THEN
                     dBxdissi = dBxdissi + termB*(Bevolxi-Bevolxyz(1,j))
                     dBydissi = dBydissi + termB*(Bevolyi-Bevolxyz(2,j))
                     dBzdissi = dBzdissi + termB*(Bevolzi-Bevolxyz(3,j))
		     ELSE
                     dBxdissi = dBxdissi + termB*dBx*robar1
		     dBydissi = dBydissi + termB*dBy*robar1
		     dBzdissi = dBzdissi + termB*dBz*robar1
                     ENDIF
                     
c
c--this is grad(div B) for diffusion of divergence errors
c
c                     dBxdissi = dBxdissi 
c     &                  - 10.*termB*(0.2*dBx - runix*projdB)*robar1
c                     dBydissi = dBydissi
c     &                  - 10.*termB*(0.2*dBy - runiy*projdB)*robar1
c                     dBzdissi = dBzdissi
c     &                  - 10.*termB*(0.2*dBz - runiz*projdB)*robar1

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
               IF (ifsvi.NE.0 .AND. projv.LT.0. .AND. 
     &                        (j.LE.npart.OR.ibound.EQ.11)) THEN
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
                     alphamean = (alphaMMpass(1,ipart) + 
     &                       alphaMMpass(1,j))/2.0
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
c           add self contribution to potential
            poten(ipart) = poten(ipart) + poteni 
     &                                  + xyzmh(4,ipart)*fpoten(0)*hi1
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
c--Anisotropic magnetic force, div/curl B
c
         IF (imhd.EQ.idim) THEN
            fxyzu(1,ipart) = fxyzu(1,ipart) + fanisoxi*cnormk
            fxyzu(2,ipart) = fxyzu(2,ipart) + fanisoyi*cnormk
            fxyzu(3,ipart) = fxyzu(3,ipart) + fanisozi*cnormk
            divcurlB(1,ipart) = cnormk*divBi !!*rho1i
            divcurlB(2,ipart) = cnormk*curlBxi*rho1i
            divcurlB(3,ipart) = cnormk*curlByi*rho1i
            divcurlB(4,ipart) = cnormk*curlBzi*rho1i
c
c--time derivative of B/rho, B (constructed from B/rho and rho derivatives) 
c  or the Euler potentials (dissipation only)
c
            IF (varmhd(1:1).EQ.'B') THEN
               dBxyz(1,ipart) = cnormk*(dBxideali*rho21i + dBxdissi)
               dBxyz(2,ipart) = cnormk*(dByideali*rho21i + dBydissi)
               dBxyz(3,ipart) = cnormk*(dBzideali*rho21i + dBzdissi)
               IF (varmhd.EQ.'Bvol') THEN
                  dBxyz(1,ipart) = rhoi*dBxyz(1,ipart) 
     &                            + Bxyz(1,ipart)*rho1i*divv(ipart)
                  dBxyz(2,ipart) = rhoi*dBxyz(2,ipart) 
     &                            + Bxyz(2,ipart)*rho1i*divv(ipart)
                  dBxyz(3,ipart) = rhoi*dBxyz(3,ipart) 
     &                            + Bxyz(3,ipart)*rho1i*divv(ipart)
               ENDIF
            ELSE
               dBxyz(1,ipart) = cnormk*dBxdissi
               dBxyz(2,ipart) = cnormk*dBydissi
               dBxyz(3,ipart) = cnormk*dBzdissi
            ENDIF
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
c
c--Morris & Monaghan switch source and decay terms for both artificial viscosity
c  and artificial resistivity (see Price & Monaghan 2005)
c
         IF (ifsvi.EQ.6) THEN
            vsigi = SQRT(vs2i)
            tdecay1 = 0.2*vsigi/hi
c--balsara factor
c            adivi = ABS(divv(ipart)*rho1i)
c            acurlvi = ABS(curlv(ipart)*rho1i)
c            fi = adivi/(adivi+acurlvi+epsil2*vsigi/hi)

cc         dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*tdecay1 -
cc     &           MIN(divv(ipart)/rhoi+0.5*vsoundi/hi,0.0)
            dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*tdecay1 -
     &             MIN(divv(ipart)/rhoi,0.0) !!*fi
            IF (imhd.EQ.idim) THEN
               dB2 = divcurlB(1,ipart)**2 + divcurlB(2,ipart)**2
     &             + divcurlB(3,ipart)**2 + divcurlB(4,ipart)**2
               source = SQRT(MAX(dB2*rho1i,
     &                           vs2i*divcurlB(1,ipart)**2/B2i))
               dha(3,ipart) = (alphamin(2)-alphaMMpass(2,ipart))*tdecay1
     &                        + source
            ENDIF
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
c--Energy conservation
c

         CALL energ(ipart,realtime,vxyzu, fxyzu, xyzmh)  
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
c--Forces that operate on both sink particles and gas from here on:
c
 80      CONTINUE
c
c--External forces
c
         IF (iexf.GE.1) CALL externf(ipart,realtime,xyzmh,fxyzu,
     &      trho,iexf)
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
c
c--Set flag for whether a particle has been updated or not back to false
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlistupdated
         iupdated(listparents(i)) = .FALSE.
      END DO
C$OMP END DO 
C$OMP END PARALLEL
c
c--Calculate gravity on and from point masses (when sink particles done outside
c     or partially outside the tree) - done separately for higher accuracy
c
      IF (nptmass.GT.0) THEN
         IF (iptintree.EQ.0) THEN
            CALL gptall(xyzmh,npart,fxyzu)
         ELSEIF (iptintree.EQ.1) THEN
            CALL gforspt(xyzmh,fxyzu)
         ENDIF
      ENDIF
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine forcei')

      RETURN
      END


      SUBROUTINE extrapolate(j,dt,itime,trho,vxyzu,pr,vsound,ekcle,
     &                       Bxyz,Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine extrapolates the density and quantities  *
c  which depend on the density for non-active particles     *
c                                                           *
c************************************************************
      INCLUDE 'idim'
      
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/timei'
      
      DIMENSION vxyzu(4,idim)
      REAL*4 trho(idim), pr(idim), vsound(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bxyz(3,imhd)
      DIMENSION Bevolxyz(3,imhd)
      
      IF (it1(j).EQ.imax) THEN
         deltat = dt*(itime - it0(j) - isteps(j)/2)/imaxstep
      ELSE
         deltat = dt*(itime - it0(j))/imaxstep
      ENDIF
c
c--Update the density value at neighbor's locations
c--Avoid, though, abrupt changes in density
c
      deltarho = -deltat*divv(j)
      IF (ABS(deltarho).GT.rho(j)/2.) THEN
         deltarho = SIGN(1.0,deltarho)*rho(j)/2.0
      ENDIF
C$OMP CRITICAL (dumrhoj)
      trhoj = rho(j) + deltarho
      CALL eospg(j, vxyzu, trho, pr, vsound, ekcle)
      
      IF (imhd.EQ.idim) THEN
c
c--use interpolated density to update B from B/rho if necessary
c  (NB the equivalent is not done for the Euler potentials as this
c   would involve too much work - so they are slightly wrong but
c   hopefully not much)
c
         IF (varmhd.EQ.'Brho') THEN
            Bxyz(1,j) = Bevolxyz(1,j)*trho(j)
            Bxyz(2,j) = Bevolxyz(2,j)*trho(j)
            Bxyz(3,j) = Bevolxyz(3,j)*trho(j)
         ENDIF
      ENDIF
C$OMP END CRITICAL (dumrhoj)
      RETURN
      
      END
