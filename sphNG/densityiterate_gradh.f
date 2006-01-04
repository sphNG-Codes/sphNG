      SUBROUTINE densityiterate_gradh (dt,npart,ntot,xyzmh,vxyzu,
     &            nlst_in,nlst_end,list,itime,ekcle)
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

c--Neides=INT(4./3.*pi*8.*hfact**3))
      PARAMETER (hfact = 1.2)

      PARAMETER (rhomin = 0.0)
      PARAMETER (htol = 1.e-3)
      PARAMETER (hstretch = 1.01)
      PARAMETER (maxiterations = 500)

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

      LOGICAL*1 iupdated(idim)

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
      ihasghostcount = 0
      numparticlesdone = numparticlesdone + nlst_end

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,list,divv,curlv,gradhs)
C$OMP& shared(hmax,xyzmh,vxyzu,pr,vsound,rho,ekcle)
C$OMP& shared(nneigh,neighb,neighover,selfnormkernel)
C$OMP& shared(cnormk,radkernel,dvtable,wij,grwij)
C$OMP& shared(listpm,iphase,dphidh,uradconst,icall,encal)
C$OMP& shared(iprint,nptmass,iptmass,radcrit2,iorig,third)
C$OMP& shared(dumrho,iscurrent,npart,hasghost)
C$OMP& shared(isteps,it0,it1,imax,imaxstep,dt,itime,iupdated)
C$OMP& private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi)
C$OMP& private(divvi,curlvxi,curlvyi,curlvzi,gradhi,gradsofti)
C$OMP& private(pmassj,hi_old,hi1,hi21,hi31,hi41,hneigh)
C$OMP& private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2)
C$OMP& private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij)
C$OMP& private(projv,procurlvx,procurlvy,procurlvz)
C$OMP& private(l,iptcur,dphi,dwdhi,dpotdh,iteration,numneighi)
C$OMP& private(numneighreal,rhohi,dhdrhoi,omegai,func,dfdh1)
C$OMP& private(hnew,deltat,deltarho)
C$OMP& reduction(MAX:rhonext,imaxit)
C$OMP& reduction(+:ihasghostcount,inumit,inumfixed,inumrecalc)
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
               CALL getneighi(ipart,xi,yi,zi,hneigh,numneighi,
     &              neighlist,xyzmh)
               inumrecalc = inumrecalc + 1
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
                  iupdated(j) = .FALSE.
                  IF (hasghost(j)) ihasghostcount = ihasghostcount + 1
                  rij1 = SQRT(rij2)
c
c--Get kernel quantities from interpolation in table
c
                  index = v2/dvtable
                  dxx = v2 - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dwdx = (wij(index1) - wij(index))/dvtable
                  wtij = (wij(index) + dwdx*dxx)*hi31
                  dgrwdx = (grwij(index1) - grwij(index))/dvtable
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
            rhohi = pmassi*(hfact*hi1)**3
            dhdrhoi = -hi/(3.*(rhohi + rhomin))
            omegai = 1. - dhdrhoi*gradhi
c
c--Newton-Raphson iteration
c
            func = rhohi - rhoi
            dfdh1 = dhdrhoi/omegai
            hnew = hi - func*dfdh1
c
c--Take fixed point if Newton-Raphson running into trouble
c  (i.e. if gradients are very wrong)
c
            IF (hnew.LE.0. .OR. (omegai.LE.tiny)) THEN
c               WRITE (*,*) 'doing fixed point',gradhi,omegai,hfact,
c     &              pmassi,rhoi,third
               hnew = hfact*(pmassi/rhoi)**third
               inumfixed = inumfixed + 1
c
c--Don't allow sudden jumps to huge numbers of neighbours
c
            ELSEIF (hnew.GT.1.2*hi .OR. hnew.LT.0.8*hi) THEN
c               WRITE (*,*) 'large h jump on particle ',iorig(ipart)
               hnew = hi - 0.5*func*dfdh1
            ENDIF

c            IF (numneighreal.GT.500) THEN
c               WRITE(iprint,*) 'part: ',iorig(ipart),' has ',numneighi,
c     &              numneighreal,' neighbours '
c            ENDIF

            IF (ABS(hnew-hi)/hi_old.LT.htol .AND. omegai.GT.0) THEN
               imaxit = MAX(imaxit, iteration)
               inumit = inumit + iteration
               GOTO 30 
            ENDIF

            hi = hnew
         END DO
         WRITE (iprint,*) 'ERROR: density iteration failed'
         WRITE (iprint,*) iorig(ipart),numneighi,numneighreal,hnew,
     &        hi,hi_old,omegai
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
c--Calculate other quantities by looping over interacting neighbours.
c     Simultaneously interpolate neighbours to current time for those
c     particles that are not active for this call to derivi.
c
         divvi = 0.
         curlvxi = 0.
         curlvyi = 0.
         curlvzi = 0.
         gradsofti = 0.

         vxi = vxyzu(1,ipart)
         vyi = vxyzu(2,ipart)
         vzi = vxyzu(3,ipart)

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
               index = v2/dvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index.GE.itable) THEN
                  index = itable
                  index1 = itable
               ENDIF
               dwdx = (wij(index1) - wij(index))/dvtable
               wtij = (wij(index) + dwdx*dxx)*hi31
               dgrwdx = (grwij(index1) - grwij(index))/dvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hi41/rij1
               dpotdh = (dphidh(index1) - dphidh(index))/dvtable
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
c--Interpolate density, pressure and vsound of neighbour that is not 
c     being evolved
c
               IF (.NOT.iscurrent(j) .AND. .NOT.iupdated(j)) THEN
                  iupdated(j) = .TRUE.
C$OMP FLUSH(iupdated)
                  IF (j.LE.npart) THEN
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
                     dumrho(j) = rho(j) + deltarho
                     CALL eospg(j, vxyzu, dumrho, pr, vsound, ekcle)
C$OMP END CRITICAL (dumrhoj)
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
c--Pressure and sound velocity from ideal gas law...
c
         CALL eospg(ipart, vxyzu, rho, pr, vsound, ekcle)
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

 50   CONTINUE
      END DO
C$OMP END PARALLEL DO
      IF (ihasghostcount.GT.0) THEN
C$OMP PARALLEL DO SCHEDULE (static) default(none)
C$OMP&shared(npart,ntot,ireal,dumrho,pr,vsound,divv,vxyzu,encal,ibound)
C$OMP&shared(ekcle)
C$OMP&private(i,j) 
         DO i = npart + 1, ntot
            j = ireal(i)
            dumrho(i) = dumrho(j)
            pr(i) = pr(j)
            vsound(i) = vsound(j)
            divv(i) = 0.
            vxyzu(4,i) = vxyzu(4,j)
            IF (encal.EQ.'r' .AND. ibound.EQ.100)
     &           ekcle(1,i) = ekcle(1,j)
         END DO
C$OMP END PARALLEL DO
      ENDIF
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
