      SUBROUTINE densityiterate_gradh (dt,npart,ntot,xyzmh,vxyzu,
     &            nlst_in,nlst_end,list,itime,ekcle,Bevol,Bxyz)
c************************************************************
c                                                           *
c  Subroutine to compute the density and smoothing lengths  *
c     self-consistently using iteration if necessary.       *
c     It also calculates the velocity divergence, the curl, *
c     and the gravity softening term with variable h, and   *
c     interpolates the density for particles which are      *
c     neighbours of particles that are currently in the     *
c     list.  This subroutine uses the binary tree algorithm *
c     to locate neighbours.  Neighbours are not stored in   *
c     a list.                                               *
c                                                           *
c     Code written by MRB and DJP (14/12/2005).             *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,idim), vxyzu(4,idim), list(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevol(3,imhd),Bxyz(3,imhd)

c--Neides=INT(4./3.*pi*8.*hfact**3))
      PARAMETER (hfact = 1.2)

      PARAMETER (htol = 1.e-3)
      PARAMETER (hstretch = 1.01)
      PARAMETER (maxiterations = 500)
c--this weight is equivalent to m/(rho*h^3) in the grad h version
      PARAMETER (weight = 1./hfact**3)

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
c      INCLUDE 'COMMONS/vsmooth'

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine densityiterate')
c
c--Initialise
c
      uradconst = radconst/uergcc
      third = 1./3.
      rhonext = 0.
      icreate = 0
      radcrit2 = radcrit*radcrit
      numparticlesdone = numparticlesdone + nlst_end
      nwarnup = 0
      nwarndown = 0
      stressmax = 0.
c
c--for constant pressure boundaries, use a minimum density
c  equal to the external density
c
c      IF (ibound.EQ.7) THEN
c         rhomin = 0.25*rhozero
c         print*,'rhomin = ',rhomin
c      ELSE
         rhomin = 0.      
c      ENDIF

C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,list,divv,curlv,gradhs)
C$OMP& shared(hmax,xyzmh,vxyzu,pr,vsound,rho,ekcle,rhomin)
C$OMP& shared(nneigh,neighb,neighover,selfnormkernel)
C$OMP& shared(cnormk,radkernel,dvtable,ddvtable,wij,grwij)
C$OMP& shared(listpm,iphase,dphidh,uradconst,icall,encal)
C$OMP& shared(iprint,nptmass,iptmass,radcrit2,iorig,third)
C$OMP& shared(dumrho,iscurrent,npart,ibound,ntot,ireal)
C$OMP& shared(isteps,it0,it1,imax,imaxstep,dt,itime)
C$OMP& shared(varmhd,Bevol,Bxyz,vsmooth,Bsmooth)
C$OMP& private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi)
C$OMP& private(divvi,curlvxi,curlvyi,curlvzi,gradhi,gradsofti)
C$OMP& private(pmassj,hi_old,hi1,hi21,hi31,hi41,hneigh)
C$OMP& private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2,rcut)
C$OMP& private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij)
C$OMP& private(projv,procurlvx,procurlvy,procurlvz)
C$OMP& private(l,iptcur,dphi,dwdhi,dpotdh,iteration,numneighi)
C$OMP& private(numneighreal,rhohi,dhdrhoi,omegai,func,dfdh1)
C$OMP& private(hnew,deltat,deltarho)
C$OMP& private(dalphaxi,dalphayi,dalphazi,dbetaxi,dbetayi,dbetazi)
C$OMP& private(rxxi,rxyi,rxzi,ryyi,ryzi,rzzi)
C$OMP& private(gradalphaxi,gradalphayi,gradalphazi)
C$OMP& private(gradbetaxi,gradbetayi,gradbetazi,term)
C$OMP& private(vbarxi,vbaryi,vbarzi)
C$OMP& private(alphai,betai,wkern,dalpha,dbeta,grpmi)
C$OMP& reduction(MAX:rhonext,imaxit)
C$OMP& reduction(+:inumit,inumfixed,inumrecalc,nwarnup,nwarndown)

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)

         IF (iphase(ipart).GE.1) GOTO 50

         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         pmassi = xyzmh(4,ipart)
         hi = xyzmh(5,ipart)
         hi_old = hi
         hneigh = 0.
c
c--Predict h
c
         dhdrhoi = - hi/(3.*(pmassi*(hfact/hi)**3 + rhomin))
         IF (it1(ipart).EQ.imax) THEN
            deltat = (dt*isteps(ipart)/2)/imaxstep
         ELSE
            deltat = (dt*isteps(ipart))/imaxstep
         ENDIF
         hi = hi - dhdrhoi*divv(ipart)*deltat
c
c--Iterate density calculation for particle ipart
c
         DO iteration = 1, maxiterations
            hi1 = 1./hi
            hi21 = hi1*hi1
            hi31 = hi21*hi1
            hi41 = hi21*hi21

            IF (hi.GT.hneigh) THEN
               hneigh = hstretch*hi
               rcut = hneigh*radkernel
               CALL getneighi(ipart,xi,yi,zi,rcut,
     &              numneighi,neighlist,xyzmh)
               inumrecalc = inumrecalc + 1
               IF (numneighi.EQ.0) WRITE(*,*) 'ZERO NEIGHBOURS!'
            ENDIF
c
c--Calculate density by looping over interacting neighbors
c
            rhoi = 0.
            gradhi = 0.
            numneighreal = 0
            DO k = 1, numneighi
               j = neighlist(k)

               dx = xi - xyzmh(1,j)
               dy = yi - xyzmh(2,j)
               dz = zi - xyzmh(3,j)
               pmassj = xyzmh(4,j)

               rij2 = dx*dx + dy*dy + dz*dz + tiny
               v2 = rij2*hi21

               IF (v2.LT.radkernel**2) THEN
                  numneighreal = numneighreal + 1
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
            END DO
c
c--Add self contribution
c
            rhoi = cnormk*(rhoi + selfnormkernel*pmassi*hi31)
            gradhi = cnormk*(gradhi +
     &           selfnormkernel*pmassi*(-3.*hi41))  
c
c--Iteration business
c  These lines define the relationship between h and rho
c     omega is the term in the denominator as in Monaghan(2001)
c
            rhohi = pmassi*(hfact*hi1)**3 - rhomin
            dhdrhoi = -hi/(3.*(rhoi + rhomin))
            omegai = 1. - dhdrhoi*gradhi
c
c--Newton-Raphson iteration
c
            func = rhohi - rhoi
            dfdh1 = dhdrhoi/omegai
            hnew = hi - func*dfdh1
c
c--Don't allow sudden jumps to huge numbers of neighbours
c
            IF (hnew.GT.1.2*hi) THEN
               nwarnup = nwarnup + 1
c               WRITE (*,*) 'restricting h jump (up) on particle ',
c     &               iorig(ipart),hi,hnew
               hnew = 1.2*hi !!hi - 0.5*func*dfdh1
            ELSEIF (hnew.LT.0.8*hi) THEN
               nwarndown = nwarndown + 1
c               WRITE (*,*) 'restricting h jump (down) on particle ',
c     &               iorig(ipart),hi,hnew
               hnew = 0.8*hi
            ELSEIF (hnew.LE.0. .OR. (omegai.LE.tiny)) THEN
c
c--Take fixed point if Newton-Raphson running into trouble
c  (i.e. if gradients are very wrong)
c
               WRITE (*,*) 'doing fixed point',hnew,gradhi,omegai,hfact,
     &              pmassi,rhoi,hi,hi31
               hnew = hfact*(pmassi/rho(ipart))**third
               inumfixed = inumfixed + 1           
            ENDIF

c            IF (numneighreal.GT.500) THEN
c               WRITE(iprint,*) 'part: ',iorig(ipart),' has ',numneighi,
c     &              numneighreal,' neighbours '
c            ENDIF
            IF (numneighreal.LE.0) THEN
               WRITE (iprint,*) 'WARNING: particle ',ipart,
     &            ' has no neighbours h=',hi,hi_old,'setting h=',
     &              hfact*(pmassi/rho(ipart))**third
               hnew = max(hfact*(pmassi/rho(ipart))**third,1.1*hnew)
            ENDIF
            
            IF (ABS(hnew-hi)/hi_old.LT.htol .AND. omegai.GT.0.) THEN
               imaxit = MAX(imaxit, iteration)
               inumit = inumit + iteration
               GOTO 30 
            ENDIF

            hi = hnew
         END DO
         WRITE (iprint,*) 'ERROR: density iteration failed'
         WRITE (iprint,*) iorig(ipart),numneighi,numneighreal,hnew,
     &        hi,hi_old,omegai,dhdrhoi,rhoi
         CALL quit
c
c--Store quantities if converged
c
 30      rho(ipart) = rhoi
         dumrho(ipart) = rhoi
         IF (numneighreal.GT.nneighmax) THEN
            WRITE (iprint,*) 'ERROR: numneighreal exceeds nneighmax'
            WRITE (iprint,*) iorig(ipart),numneighreal,nneighmax
            CALL quit
         ENDIF

         nneigh(ipart) = numneighreal
         xyzmh(5,ipart) = hi
         gradhs(1,ipart) = 1./omegai
c
c--Set e(i) for the first time around because density only set here
c     This sets radiation and matter to have same initial temperature
c
         IF (icall.EQ.1 .AND. encal.EQ.'r') THEN
            IF (ekcle(1,ipart).EQ.0.0) THEN
               ekcle(3,ipart) = getcv(rho(ipart),vxyzu(4,ipart))
               ekcle(1,ipart) = uradconst*(vxyzu(4,ipart)/
     &              ekcle(3,ipart))**4/rho(ipart)
            ENDIF
         ENDIF
c
c--Pressure and sound velocity from ideal gas law...
c
         CALL eospg(ipart, vxyzu, rho, pr, vsound, ekcle)
c
c--Calculate other quantities by looping over interacting neighbours.
c     Simultaneously interpolate neighbours to current time for those
c     particles that are not active for this call to derivi.
c
         divvi = 0.
         curlvxi = 0.
         curlvyi = 0.
         curlvzi = 0.
         gradsofti = 0.
         dalphaxi= 0.
         dalphayi= 0.
         dalphazi= 0.
         dbetaxi= 0.
         dbetayi= 0.
         dbetazi= 0.
         rxxi = 0.
         rxyi = 0.
         rxzi = 0.
         ryyi = 0.
         ryzi = 0.
         rzzi = 0.
         
c         vbarxi = 0.
c         vbaryi = 0.
c         vbarzi = 0.
c         Bbarxi = 0.
c         Bbaryi = 0.
c         Bbarzi = 0.

         vxi = vxyzu(1,ipart)
         vyi = vxyzu(2,ipart)
         vzi = vxyzu(3,ipart)
c
c--Calculate B from the evolved magnetic field variable
c
         IF (imhd.EQ.idim) THEN
            IF (varmhd.EQ.'eulr') THEN
               alphai= Bevol(1,ipart)
               betai= Bevol(2,ipart)
            ENDIF
         ENDIF
                  
         DO k = 1, numneighi
            j = neighlist(k)

            dx = xi - xyzmh(1,j)
            dy = yi - xyzmh(2,j)
            dz = zi - xyzmh(3,j)
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
c--Velocity divergence times density
c
               projv = grwtij*(dvx*dx + dvy*dy + dvz*dz)
               divvi = divvi - pmassj*projv
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
c--Derivative of gravitational potential w.r.t. h
c
               gradsofti = gradsofti - pmassj*dphi
c
c--get B from the evolved Euler potentials
c
c
c--grad alpha and grad beta (Euler potentials)
c
               IF (imhd.EQ.idim) THEN
                  IF (varmhd.EQ.'eulr') THEN
                     grpmi= pmassj*grwtij

                     rxxi = rxxi - grpmi*dx*dx
                     rxyi = rxyi - grpmi*dx*dy
                     rxzi = rxzi - grpmi*dx*dz
                     ryyi = ryyi - grpmi*dy*dy
                     ryzi = ryzi - grpmi*dy*dz
                     rzzi = rzzi - grpmi*dz*dz
                     dalpha= alphai - Bevol(1,j)
                     dbeta= betai - Bevol(2,j)

                     dalphaxi= dalphaxi - grpmi*dalpha*dx
                     dalphayi= dalphayi - grpmi*dalpha*dy
                     dalphazi= dalphazi - grpmi*dalpha*dz

                     dbetaxi= dbetaxi - grpmi*dbeta*dx
                     dbetayi= dbetayi - grpmi*dbeta*dy
                     dbetazi= dbetazi - grpmi*dbeta*dz
c                  ELSE
c
c--smoothed velocity for use in the B or B/rho evolution
c
c                     vbarxi = vbarxi + weight*vxyzu(1,j)*wkern
c                     vbaryi = vbaryi + weight*vxyzu(2,j)*wkern
c                     vbarzi = vbarzi + weight*vxyzu(3,j)*wkern
c
c--smoothed Bevol for div B reduction
c
c                     Bbarxi = Bbarxi + weight*Bevol(1,j)*wkern
c                     Bbaryi = Bbaryi + weight*Bevol(2,j)*wkern
c                     Bbarzi = Bbarzi + weight*Bevol(3,j)*wkern
                  ENDIF
               ENDIF
            ENDIF
         END DO
c
c--Add self contribution
c
         divv(ipart) = cnormk*divvi*gradhs(1,ipart)
         curlv(ipart) = cnormk*SQRT(curlvxi**2+curlvyi**2+curlvzi**2)*
     &        gradhs(1,ipart)
         gradhs(2,ipart) = dhdrhoi*(gradsofti - pmassi*dphidh(0)*hi21) 
c
c--Find particle with highest density outside radcrit of point mass
c
         IF (iptmass.NE.0) THEN
            DO l = 1, nptmass
               iptcur = listpm(l)
               IF ( (xi - xyzmh(1,iptcur))**2 + 
     &              (yi - xyzmh(2,iptcur))**2 + 
     &              (zi - xyzmh(3,iptcur))**2 .LT.radcrit2) GOTO 50
            END DO
            rhonext = MAX(rhonext, rho(ipart))
         ENDIF
c
c--calculate B from the evolved Euler potentials
c
         IF (imhd.EQ.idim) THEN
            IF (varmhd.EQ.'eulr') THEN
c
c--compute grad alpha and grad beta using exact linear interpolation
c  (see Price 2004)
c
               ddenom = 1./(rxxi*ryyi*rzzi + 2.*rxyi*rxzi*ryzi
     &                    - rxxi*ryzi**2 - ryyi*rxzi**2 - rzzi*rxyi**2)
     
               gradalphaxi = (dalphaxi*(ryyi*rzzi - ryzi**2)
     &                     +  dalphayi*(rxzi*ryzi - rzzi*rxyi)
     &                     +  dalphazi*(rxyi*ryzi - rxzi*ryyi))*ddenom
               gradalphayi = (dalphaxi*(ryzi*rxzi - rxyi*rzzi)
     &                     +  dalphayi*(rzzi*rxxi - rxzi**2)
     &                     +  dalphazi*(rxyi*rxzi - rxxi*ryzi))*ddenom
               gradalphazi = (dalphaxi*(rxyi*ryzi - rxzi*ryyi)
     &                     +  dalphayi*(rxyi*rxzi - rxxi*ryzi)
     &                     +  dalphazi*(rxxi*ryyi - rxyi**2))*ddenom
               gradbetaxi = (dbetaxi*(ryyi*rzzi - ryzi**2)
     &                    +  dbetayi*(rxzi*ryzi - rzzi*rxyi)
     &                    +  dbetazi*(rxyi*ryzi - rxzi*ryyi))*ddenom
               gradbetayi = (dbetaxi*(ryzi*rxzi - rxyi*rzzi)
     &                    +  dbetayi*(rzzi*rxxi - rxzi**2)
     &                    +  dbetazi*(rxyi*rxzi - rxxi*ryzi))*ddenom
               gradbetazi = (dbetaxi*(rxyi*ryzi - rxzi*ryyi)
     &                    +  dbetayi*(rxyi*rxzi - rxxi*ryzi)
     &                    +  dbetazi*(rxxi*ryyi - rxyi**2))*ddenom
c
c--uncomment the following lines for the standard first derivative
c
c               term = cnormk*gradhs(1,ipart)/rhoi
c               gradalphaxi = dalphaxi*term
c               gradalphayi = dalphayi*term
c               gradalphazi = dalphazi*term
c               gradbetaxi = dbetaxi*term
c               gradbetayi = dbetayi*term
c               gradbetazi = dbetazi*term
c
c--grad alpha cross grad beta
c
               term= gradalphayi*gradbetazi - gradalphazi*gradbetayi
               Bxyz(1,ipart)= term
               term= gradalphazi*gradbetaxi - gradalphaxi*gradbetazi
               Bxyz(2,ipart)= term
               term= gradalphaxi*gradbetayi - gradalphayi*gradbetaxi
               Bxyz(3,ipart)= term
c            ELSE
c
c--add self contribution and store smoothed velocity
c
c               vsmooth(1,ipart) = cnormk*(vbarxi 
c     &                                  + weight*vxyzu(1,ipart)*wij(0))
c               vsmooth(2,ipart) = cnormk*(vbaryi
c     &                                  + weight*vxyzu(2,ipart)*wij(0))
c               vsmooth(3,ipart) = cnormk*(vbarzi
c     &                                  + weight*vxyzu(3,ipart)*wij(0))
c               Bsmooth(1,ipart) = cnormk*(Bbarxi 
c     &                                  + weight*Bevol(1,ipart)*wij(0))
c               Bsmooth(2,ipart) = cnormk*(Bbaryi
c     &                                  + weight*Bevol(2,ipart)*wij(0))
c               Bsmooth(3,ipart) = cnormk*(Bbarzi
c     &                                  + weight*Bevol(3,ipart)*wij(0))
c               print*,ipart,'v       = ',vxyzu(1:3,ipart)
c               print*,ipart,'vsmooth = ',vsmooth(:,ipart)
            ENDIF
         ENDIF
 50   CONTINUE
      END DO
C$OMP END DO
      IF (nwarnup.GT.0) THEN
         WRITE (iprint,*) 'WARNING: restricted h jump (up) ',
     &        nwarnup,' times'
      ENDIF
      IF (nwarndown.GT.0) THEN
         WRITE (iprint,*) 'WARNING: restricted h jump (down) ',
     &        nwarndown,' times'
      ENDIF

c
c--copy changed values onto ghost particles
c
C$OMP DO SCHEDULE (runtime)
         DO i = npart + 1, ntot
            j = ireal(i)
            rho(i) = rho(j)
            dumrho(i) = dumrho(j)
            xyzmh(5,i) = xyzmh(5,j)
            pr(i) = pr(j)
            vsound(i) = vsound(j)
            divv(i) = divv(j)
            curlv(i) = curlv(j)
            vxyzu(4,i) = vxyzu(4,j)
            gradhs(1,i) = gradhs(1,j)
            gradhs(2,i) = gradhs(2,j)
            IF (imhd.EQ.idim) THEN
               Bxyz(1,i) = Bxyz(1,j)
               Bxyz(2,i) = Bxyz(2,j)
               Bxyz(3,i) = Bxyz(3,j)
            ENDIF
            IF (encal.EQ.'r' .AND. ibound.EQ.100) 
     &           ekcle(1,i) = ekcle(1,j)
         END DO
C$OMP END DO
C$OMP END PARALLEL
ccC$OMP END PARALLEL DO
c
c--Possible to create a point mass
c
      IF (rhonext.GT.rhocrea .AND. nptmass.LT.iptdim .AND.
     &                     iptmass.NE.0 .AND. icall.EQ.3) THEN
c
c--Find particle with highest density outside radcrit of point mass
c
         DO n = nlst_in, nlst_end
            ipart = list(n)
            IF (rho(ipart).EQ.rhonext) THEN
               irhonex = ipart
            ENDIF
         END DO
c
c--Make sure that all neighbours of point mass candidate are being
c     done on this time step. Otherwise, not possible to accrete
c     them to form a point mass and it may create a point mass without
c     accreting many particles!
c
         IF ((2.0*xyzmh(5,irhonex)).LT.hacc) THEN
            WRITE(iprint,*)'Ptmass creation passed h ',xyzmh(5,irhonex)

            CALL getneigh(irhonex,npart,xyzmh(5,irhonex),xyzmh,nlist,
     &           iptneigh,nearl)

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
               WRITE(iprint,*) ' and all particles on step'
            ELSE
               WRITE(iprint,*) ' but not all particles on step'
            ENDIF
         ELSE
            WRITE(iprint,*) 'Ptmass creation failed on h ',
     &           xyzmh(5,irhonex)
         ENDIF
      ENDIF

      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine densityiterate')

      RETURN
      END
