      SUBROUTINE accrete(dt, realtime, isave)
c************************************************************
c                                                           *
c  This routine does accretion for pt masses                *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/accnum'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/accurpt'
      INCLUDE 'COMMONS/binfile'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/accrem'
      INCLUDE 'COMMONS/timeextra'

      DIMENSION numberacc(iptdim)
      REAL*4 ptminner(iptdim)
      CHARACTER*7 where
      
      DATA where/'accrete'/
c
c--Only allow accretion of GAS particles evaluated at the CURRENT timestep
c     iremove is initialised in evol.f to -1
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst0,llist,iphase,iremove,nptmass,numberacc,ptminner)
C$OMP& private(i,j,iii)

C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst0
         j = llist(i)
         IF (iphase(j).EQ.0) iremove(j) = 0
      END DO
C$OMP END DO

C$OMP DO SCHEDULE(runtime)
      DO iii = 1, nptmass
         numberacc(iii) = 0
         ptminner(iii) = 0.
      END DO
C$OMP END DO
C$OMP END PARALLEL

c
c--CREATION OF NEW POINT MASS, by accretion of particles to form it
c
      IF (icreate.EQ.1 .AND. iremove(irhonex).EQ.0) THEN
c
c--Test before creating a new point mass:
c
c     (a) whether the particle(irhonex) contains a jeans mass within 2h
c           =>calculate alpha
c
c     (b) calculate beta ratio of rotational to gravitational pot. energy
c
c     (c) calculate the divergence of the acceleration
c           -ve => self gravitating/collapsing
c           +ve => in process of tidal/disruption or core bounce
c 
c
c--Fraction of hacc for accreting a particle regardless of tests
c
         hacc2  = hacc*hacc
         IF(iptmass.EQ.1) THEN
            haccmin2 = hacc2
         ELSE
            haccmin2 = haccall*haccall
         ENDIF

         xi = xyzmh(1,irhonex)
         yi = xyzmh(2,irhonex)
         zi = xyzmh(3,irhonex)
         pmassi = xyzmh(4,irhonex)
         hi = xyzmh(5,irhonex)
         vxi = vxyzu(1,irhonex)
         vyi = vxyzu(2,irhonex)
         vzi = vxyzu(3,irhonex)
         ui  =  vxyzu(4,irhonex)
         f1vxi = f1vxyzu(1,irhonex)
         f1vyi = f1vxyzu(2,irhonex)
         f1vzi = f1vxyzu(3,irhonex)
         f1ui = f1vxyzu(4,irhonex)
         f1hi = f1ha(1,irhonex)
         rhoi = rho(irhonex)
         dgravi = dgrav(irhonex)
         poteni = poten(irhonex)
         pri = pr(irhonex)
         divvi = divv(irhonex)
         it0i = it0(irhonex)
         it1i = it1(irhonex)
         it2i = it2(irhonex)
         istepsi = isteps(irhonex)
         spinxi = 0.
         spinyi = 0.
         spinzi = 0.

         tkin = 0.
         trotx = 0.
         troty = 0.
         trotz = 0.
         tgrav = 0.
         divai = 0.
         gama1 = gamma - 1.0
         IF ( varsta.NE.'entropy' ) THEN
            tterm = pmassi*ui
         ELSEIF (gama1.EQ.0.) THEN
            tterm = 1.5*ui
         ELSE 
            tterm = pmassi*ui*rhoi**gama1/gama1
         ENDIF

         CALL getneigh(irhonex,npart,hi,xyzmh,nlist,iptneigh,nearl)
         write (iprint,*) 'getneigh ',iorig(irhonex),npart,hi,nlist

         nptlist(nptmass+1) = nlist

         nvalid = 0
         DO k = 1, nlist
            j = nearl(k)
            nearpt(k,nptmass+1) = j
            dx = xyzmh(1,j) - xi
            dy = xyzmh(2,j) - yi
            dz = xyzmh(3,j) - zi
            pmassj = xyzmh(4,j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            IF (iremove(j).EQ.0 .AND. rij2.LT.hacc2) THEN
               nvalid = nvalid + 1
               hmean = 0.5*(hi + xyzmh(5,j))
               hmean21 = 1./(hmean*hmean)
               hmean31 = hmean21/hmean
               hmean41 = hmean21*hmean21
               dvx = vxyzu(1,j) - vxi
               dvy = vxyzu(2,j) - vyi
               dvz = vxyzu(3,j) - vzi
c
c--Relative kinetic energy, tkin
c
               vtot2 = dvx*dvx + dvy*dvy + dvz*dvz
               tkin = tkin + pmassj*vtot2
c
c--Relative rotational energy around x, trotx
c
               r2yz = dz*dz + dy*dy
               rvx = dy*dvz - dz*dvy
               IF(r2yz.NE.0.) trotx = trotx + pmassj*rvx*rvx/r2yz
c
c--Relative rotational energy around y, troty
c            
               r2xz = dx*dx + dz*dz
               rvy = dz*dvx - dx*dvz
               IF(r2xz.NE.0.) troty = troty + pmassj*rvy*rvy/r2xz
c
c--Relative rotational energy around z, trotz
c            
               r2xy = dx*dx + dy*dy
               rvz = dx*dvy - dy*dvx
               IF(r2xy.NE.0.) trotz = trotz + pmassj*rvz*rvz/r2xy

               v2 = rij2*hmean21
               rij = SQRT(rij2)
               rij1 = 1.0/rij 
               v = rij/hmean
c
c--Get kernel quantities from interpolation in table
c
               IF (v.LT.radkernel) THEN
                  index = v2/dvtable
                  dxx = v2 - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dgrwdx = (grwij(index1) - grwij(index))/dvtable
                  grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
                  dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                  phi = (fpoten(index) + dfptdx*dxx)/hmean
                  IF (v2.GT.1.) phi = phi + rij1/15.0
               ELSE
                  grwtij = 0.
                  phi = -rij1
               ENDIF
c
c--Acceleration divergence times density...
c
               dax = f1vxyzu(1,j) - f1vxi
               day = f1vxyzu(2,j) - f1vyi
               daz = f1vxyzu(3,j) - f1vzi
               proja = grwtij*(dax*dx + day*dy + daz*dz)/rij
               divai = divai - pmassj*proja
c
c--Gravitational energy of particles
c
               tgrav = tgrav + phi*pmassj*pmassi
c
c--Thermal energy, tterm
c     
               IF ( varsta.NE.'entropy' ) THEN
                  tterm = tterm + pmassj*vxyzu(4,j)
               ELSEIF (gama1.EQ.0.) THEN
                  tterm = tterm + 1.5*vxyzu(4,j)
               ELSE 
                  tterm = tterm + pmassj*vxyzu(4,j)*rho(j)**gama1/gama1
               ENDIF
            ENDIF 
         END DO
c
c--Normalise acceleration divergence, divai
c
         WRITE(iprint,*) 'Ptmass nvalid, nlist ',nvalid, nlist,
     &        iorig(irhonex),hi,rho(irhonex),iphase(irhonex),
     &        nneigh(irhonex)
         divai = cnormk*divai
         IF (divai.GE.0) THEN
            WRITE(iprint,*)'Divai +ve => no pt mass creation yet ',
     &           divai
            icreate = 0
            GOTO 100
         ENDIF
c
c--Other normalisations
c
         tkin = 0.5*tkin
         trotx = 0.5*trotx
         troty = 0.5*troty
         trotz = 0.5*trotz
         trot = SQRT(trotx*trotx + troty*troty + trotz*trotz)
c
c--Now calculate tgrav, the potential energy => alpha
c
         DO k1 = 1, nlist
            j1 = nearl(k1)
            xj1 = xyzmh(1,j1)
            yj1 = xyzmh(2,j1)
            zj1 = xyzmh(3,j1)
            pmassj1 = xyzmh(4,j1)

            rx = xj1 - xi
            ry = yj1 - yi
            rz = zj1 - zi
            r2 = rx*rx + ry*ry + rz*rz
            IF (iremove(j1).NE.0 .OR. r2.GE.hacc2) GOTO 20

            DO k2 = k1+1, nlist
               j2 = nearl(k2)

               rx = xyzmh(1,j2) - xi
               ry = xyzmh(2,j2) - yi
               rz = xyzmh(3,j2) - zi
               r2 = rx*rx + ry*ry + rz*rz
               IF (iremove(j2).NE.0 .OR. r2.GE.hacc2) GOTO 10

               dx = xyzmh(1,j2) - xj1
               dy = xyzmh(2,j2) - yj1
               dz = xyzmh(3,j2) - zj1
               pmassj2 = xyzmh(4,j2)
               rij2 = dx*dx + dy*dy + dz*dz + tiny
               rij = SQRT(rij2)
               rij1 = 1./rij
c
c--Define mean h ...
c                        
               hmean = 0.5*(xyzmh(5,j1) + xyzmh(5,j2))
               hmean21 = 1./(hmean*hmean)
               v2 = rij2*hmean21
               IF (v2.LT.radkernel*radkernel) THEN
                  index = v2/dvtable
                  dxx = v2 - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                  phi = (fpoten(index) + dfptdx*dxx)/hmean
                  IF (v2.GT.1.) phi = phi + rij1/15.0
               ELSE
                  phi = -rij1
               ENDIF
               tgrav = tgrav + phi*pmassj1*pmassj2
 10         END DO
 20      END DO
c
c--Now test to see if the particles to be turned into a point mass 
c     satisfy all of the criteria
c
         xmjeans = ABS(tgrav/tterm)
         alphapt = 1./xmjeans
         betatotal = ABS(trot/tgrav)
         alphatot = alphapt + betatotal
         total = tterm + tgrav + tkin
         IF(alphapt.GT.0.5) THEN
            WRITE(iprint,99001)
            WRITE(iprint,*)'Ptmass failed on alpha = ',alphapt
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans,
     &           alphapt,betatotal,rho(irhonex),hi,
     &           iorig(irhonex)
            icreate = 0
            CALL FLUSH (iprint)
            GOTO 100
         ELSE IF (alphatot.GT.1.0) THEN
            WRITE(iprint,99001)
            WRITE(iprint,*)'Ptmass failed on alpha + beta = ',
     &                     alphapt,betatotal
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans,
     &           alphapt,betatotal,rho(irhonex),hi,
     &           iorig(irhonex)
            icreate = 0
            CALL FLUSH (iprint)
            GOTO 100
         ELSE IF(total.GE.0) THEN
            WRITE(iprint,99001)
            WRITE(iprint,*)'Ptmass failed on total energy (pos) = ', 
     &           total
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans,
     &           alphapt,betatotal,rho(irhonex),hi,
     &           iorig(irhonex)
            icreate = 0
            CALL FLUSH (iprint)
            GOTO 100
         ELSE
            WRITE(iprint,99002)
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans,
     &           alphapt,betatotal,rho(irhonex),hi
     &           ,iorig(irhonex)
            CALL FLUSH (iprint)
99001       FORMAT(' PROTOSTAR FORMATION UNSUCCESSFUL !!!')          
99002       FORMAT(' PROTOSTAR FORMATION SUCCESSFUL !!!')          
99003       FORMAT(' total energy                   :',1PE14.5,/,
     &     ' thermal energy                 :',1PE14.5,/,
     &     ' gravitational potential energy :',1PE14.5,/,
     &     ' kinetic energy                 :',1PE14.5,/,
     &     ' Jeans no.                      :',1PE14.5,/,
     &     ' alpha                          :',1PE14.5,/,
     &     ' beta total                     :',1PE14.5,/,
     &     ' rho(irhonex)                   :',1PE14.5,/,
     &     ' h(irhonex)                     :',1PE14.5,/,
     &     ' irhonex                        :',I8)
         ENDIF
c
c--Create point mass from particles
c
c--Make a dump next time save is called
c
         iptcreat = 1
         nptmass = nptmass + 1
         IF (nptmass.GE.iptdim) CALL error(where,3)

         numberacc(nptmass) = 0
         ptminner(nptmass) = 0.0
         DO jj = 1,nptlist(nptmass)
            j = nearpt(jj,nptmass)
            rx = xyzmh(1,j)-xi
            ry = xyzmh(2,j)-yi
            rz = xyzmh(3,j)-zi
            pmassj = xyzmh(4,j)
            dvx = vxyzu(1,j)-vxi
            dvy = vxyzu(2,j)-vyi
            dvz = vxyzu(3,j)-vzi
              
            r2 = rx*rx+ry*ry+rz*rz
            IF (r2.LT.hacc2.AND.iremove(j).EQ.0) THEN
               iremove(j) = 1
               iphase(j) = -1
               iaccr = 1
               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR nlstacc'
                  CALL quit
               ENDIF
               listacc(nlstacc) = j
               numberacc(nptmass) = numberacc(nptmass) + 1
               ptminner(nptmass) = ptminner(nptmass) + pmassj
c
c--Accrete particle if it lies within hacc
c      accrete mass
c      angular momentum
c      linear momentum
c
               totalmass = pmassi + pmassj
               spinm = pmassj*pmassi/totalmass
               spinxi = spinxi + spinm*(ry*dvz - dvy*rz)
               spinyi = spinyi + spinm*(dvx*rz - rx*dvz)
               spinzi = spinzi + spinm*(rx*dvy - dvx*ry)
               vxi = (pmassi*vxi + pmassj*vxyzu(1,j))/totalmass
               vyi = (pmassi*vyi + pmassj*vxyzu(2,j))/totalmass
               vzi = (pmassi*vzi + pmassj*vxyzu(3,j))/totalmass
               ui  = (pmassi*ui  + pmassj*vxyzu(4,j))/totalmass
               xi = (pmassi*xi + pmassj*xyzmh(1,j))/totalmass
               yi = (pmassi*yi + pmassj*xyzmh(2,j))/totalmass
               zi = (pmassi*zi + pmassj*xyzmh(3,j))/totalmass
               f1vxi = (pmassi*f1vxi + pmassj*f1vxyzu(1,j))/totalmass
               f1vyi = (pmassi*f1vyi + pmassj*f1vxyzu(2,j))/totalmass
               f1vzi = (pmassi*f1vzi + pmassj*f1vxyzu(3,j))/totalmass
               WRITE(iprint,*)'add = ', jj, r2, pmassi
               pmassi = totalmass
            ENDIF
         END DO
         WRITE(iprint,77001) pmassi, realtime
77001        FORMAT('PROTOSTAR CREATION, mass = ',1PE12.5,
     &        ' time = ',1PE14.7)
         WRITE(iprint,77002) hacc
77002        FORMAT('   accretion radius = ',1PE12.5)
c
c--Set spin arrays for ptmass and change the number of point masses
c
         spinx(nptmass) = spinxi
         spiny(nptmass) = spinyi
         spinz(nptmass) = spinzi
         spinadx(nptmass) = spinxi
         spinady(nptmass) = spinyi
         spinadz(nptmass) = spinzi
c
c--Finally set new point mass's properties
c
         xyzmh(1,irhonex) = xi
         xyzmh(2,irhonex) = yi
         xyzmh(3,irhonex) = zi
         xyzmh(4,irhonex) = pmassi

         xyzmh(5,irhonex) = hacc

         vxyzu(1,irhonex) = vxi
         vxyzu(2,irhonex) = vyi
         vxyzu(3,irhonex) = vzi

         vxyzu(4,irhonex) = ui         

         xmomsyn(nptmass) = pmassi*vxi
         ymomsyn(nptmass) = pmassi*vyi
         zmomsyn(nptmass) = pmassi*vzi
         xmomadd(nptmass) = 0.0
         ymomadd(nptmass) = 0.0
         zmomadd(nptmass) = 0.0

         f1vxyzu(1,irhonex) = f1vxi
         f1vxyzu(2,irhonex) = f1vyi
         f1vxyzu(3,irhonex) = f1vzi
         f1vxyzu(4,irhonex) = f1ui
         f1ha(1,irhonex) = f1hi
         ptmsyn(nptmass) = pmassi
         ptmadd(nptmass) = 0.0
         angaddx(nptmass) = 0.0
         angaddy(nptmass) = 0.0
         angaddz(nptmass) = 0.0

         rho(irhonex) = rhoi
         dgrav(irhonex) = 0.
         poten(irhonex) = poteni
         pr(irhonex) = pri
         divv(irhonex) = divvi
         IF (istepmin.LT.istepsi) THEN
            it0(irhonex) = it0i
            it1(irhonex) = it0(irhonex) + istepmin/2
            it2(irhonex) = it0(irhonex) + istepmin
            isteps(irhonex) = istepmin
         ELSE
            it0(irhonex) = it0i
            it1(irhonex) = it1i
            it2(irhonex) = it2i
            isteps(irhonex) = istepsi
         ENDIF

         CALL getneigh(irhonex,npart,xyzmh(5,irhonex),xyzmh,nlist,
     &        iptneigh,nearl)

         nptlist(nptmass) = nlist

         DO k = 1, nlist
            nearpt(k,nptmass) = nearl(k)
         END DO

         iphase(irhonex) = iptmass
         IF (initialptm.EQ.0) initialptm = iptmass
         listpm(nptmass) = irhonex
         hasghost(irhonex) = .FALSE.
      ENDIF
c
c--ACCRETION OF PARTICLES NEAR AN EXISTING POINT MASS
c
c--Method for accreting a particle
c     iphase = 1  Point mass that accretes everything regardless of tests.
c              2  Point mass that accretes particles if they pass tests.
c              3  Point mass that accretes part of a particle until all gone,
c                    they must also pass tests.
c              4  Point mass with accretion radius boundary corrections
c                   (a) smoothing length corrections
c                   (b) local density gradient correction of density
c                   (c) local pressure gradient correction to pressure force
c                   (d) local shear viscosity correction
c
c     hacc    = outer radius at which particle's accretion begins
c     haccall = radius at which all particles are accreted without test
c
c
 100  DO iii = 1, nptmass
         IF (nptlist(iii).GT.0) THEN

         i = listpm(iii)
         
         hacccur = xyzmh(5,i)
         IF (iphase(i).EQ.3) THEN
            haccmin = haccall
         ELSE
            haccmin = hacccur
         ENDIF
         hacccur2 = hacccur*hacccur
         haccmin2 = haccmin*haccmin
         haccall2 = haccall*haccall

         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         pmassi = xyzmh(4,i)
         hi = xyzmh(5,i)
         xiold = xi
         yiold = yi
         ziold = zi
         vxi = vxyzu(1,i)
         vyi = vxyzu(2,i)
         vzi = vxyzu(3,i)
         ui = vxyzu(4,i)
         f1vxi = f1vxyzu(1,i)
         f1vyi = f1vxyzu(2,i)
         f1vzi = f1vxyzu(3,i)

         xtemp = 0.0
         ytemp = 0.0
         ztemp = 0.0
         f1vxtemp = 0.0
         f1vytemp = 0.0
         f1vztemp = 0.0
         d2vxtemp = 0.0
         d2vytemp = 0.0
         d2vztemp = 0.0
         utemp = 0.0
         spinxtemp = 0.0
         spinytemp = 0.0
         spinztemp = 0.0

         DO jj = 1, nptlist(iii)
            j = nearpt(jj,iii)
            rx = xyzmh(1,j) - xi
            ry = xyzmh(2,j) - yi
            rz = xyzmh(3,j) - zi
            r2 = rx*rx + ry*ry + rz*rz
c
c--Is point mass neighbour inside accretion radius and accretable?
c
            IF (r2.LT.hacccur2 .AND. iremove(j).EQ.0) THEN
c
c--Check to ensure that the particle is actually bound to the 
c     point mass and not just passing through its neighbourhood
c      (a) particle must be bound
c      (b) particle must be more bound to current point mass than any other
c      (c) specific angular momentum of particle must be less than that
c             required for it to form a circular orbit at hacc (8/9/94)
c
               dvx = vxyzu(1,j) - vxi
               dvy = vxyzu(2,j) - vyi
               dvz = vxyzu(3,j) - vzi

               r1 = SQRT(r2)
               divvr2 = dvx*dvx + dvy*dvy + dvz*dvz
               vrad = (dvx*rx + dvy*ry + dvz*rz)/r1
               vrad2 = vrad*vrad
               vkep2 = pmassi/r1
               vpotdif = -vkep2 + divvr2/2.
               vtan2 = divvr2 - vrad2
               specangmom2 = vtan2*r2
               specangmomhacc2 = pmassi*hacccur
c
c--Check which point mass the particle is MOST bound to (6/8/94)
c
               nboundto = iii
               DO kk = 1, nptmass
                  IF (kk.NE.iii) THEN
                     ikk = listpm(kk)
                     rxtemp = xyzmh(1,j) - xyzmh(1,ikk)
                     rytemp = xyzmh(2,j) - xyzmh(2,ikk)
                     rztemp = xyzmh(3,j) - xyzmh(3,ikk)
                     r1temp = SQRT(rxtemp*rxtemp + rytemp*rytemp + 
     &                    rztemp*rztemp)
                     dvxtemp = vxyzu(1,j) - vxyzu(1,ikk)
                     dvytemp = vxyzu(2,j) - vxyzu(2,ikk)
                     dvztemp = vxyzu(3,j) - vxyzu(3,ikk)
                     v2temp = dvxtemp*dvxtemp + dvytemp*dvytemp + 
     &                    dvztemp*dvztemp
                     vpottemp = -xyzmh(4,ikk)/r1temp + v2temp/2.
                     IF (vpottemp.LT.vpotdif) nboundto = kk
                  ENDIF
               END DO
               IF (nboundto.NE.iii) vpotdif = 1.0E+10

               IF(vpotdif.LE.0. .AND. specangmom2.LT.specangmomhacc2
     &              .OR. r2.LT.haccall2 .OR. iphase(i).EQ.1) THEN

                  IF (r2.LT.haccmin2 .OR. (iphase(i).EQ.3 .AND.
     &                 pmassjn.LT.pmassleast)) THEN
c
c--Accrete all of particle's mass, and set particle accreted (iphase=-1).
c
                     iphase(j) = -1
                     iremove(j) = 1
                     iaccr = 1
                     nlstacc = nlstacc + 1
                     IF (nlstacc.GT.nlstaccmax) THEN
                        WRITE (iprint,*) 'ERROR nlstacc'
                        CALL quit
                     ENDIF
                     listacc(nlstacc) = j
                     numberacc(iii) = numberacc(iii) + 1
                     pmassj = xyzmh(4,j)
                     pmassjn = xyzmh(4,j)
                     ptminner(iii) = ptminner(iii) + pmassj
                     IF (iphase(i).EQ.3) THEN
                        WRITE(iprint,*)'accreted ',iorig(j),r1,
     &                    pmassj/pmassleast, pmassleast
                     ENDIF
                     WRITE(iaccpr) iorig(j),realtime,
     &                    xyzmh(1,j),xyzmh(2,j),xyzmh(3,j),
     &                    vxyzu(1,j),vxyzu(2,j),vxyzu(3,j),
     &                    poten(j),dgrav(j)
                     CALL FLUSH (iaccpr)
                  ELSEIF (iphase(i).EQ.3) THEN
                     hacc2  = hacc*hacc
                     hacc11 = 2./hacc
                     hacc21 = 4./hacc2
                     hacc3  = hacc*hacc2
c
c--Point mass accretes part of a particle's mass at time
c
c--Use timestep for relevant accretion timestep
c
                     dtaccj = dt*isteps(j)/imaxstep

                     hj3 = xyzmh(5,j)**3
                     hratio = hacc3/hj3/8.
                     IF (hratio.GT.1.) hratio = 1.0
                     v2 = r2*hacc21

                     cnorm2 = 0.75 - haccall + 0.5*haccall**3 - 
     &                    0.1875*haccall**4
                     cnorm2 = 2. / cnorm2
c
c--Use smoothing kernel of point mass renormalised as accretion kernel
c
                     index = v2/dvtable
                     dxx = v2 - index*dvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dwdx = (wij(index1) - wij(index))/dvtable
                     wptj = 0.333333*(wij(index) + dwdx*dxx)*
     &                                                hacc11*cnorm2
                     pmassj = -wptj*xyzmh(4,j)*vrad*hratio*dtaccj/hacc
                     pmassjn = xyzmh(4,j) - pmassj
                  ELSE
                     CALL error(where,1)
                  ENDIF
c
c--accrete particle if it lies within h and above condtions okay
c      accrete mass
c      angular momentum
c      linear momentum
c
                  IF (notacc(j)) GOTO 888

                  pmiold = ptmsyn(iii) + ptmadd(iii)
                  ptmadd(iii) = ptmadd(iii) + pmassj
                  totalmass = ptmsyn(iii) + ptmadd(iii)
                  spinm = pmassj*pmiold/totalmass
                  spinxtemp = spinxtemp + spinm*(ry*dvz - dvy*rz)
                  spinytemp = spinytemp + spinm*(dvx*rz - rx*dvz)
                  spinztemp = spinztemp + spinm*(rx*dvy - dvx*ry)

                  xmomadd(iii) = xmomadd(iii) + pmassj*vxyzu(1,j)
                  ymomadd(iii) = ymomadd(iii) + pmassj*vxyzu(2,j)
                  zmomadd(iii) = zmomadd(iii) + pmassj*vxyzu(3,j)
                  angaddx(iii) = angaddx(iii) + pmassj*
     &                 (yi*vxyzu(3,j) - vxyzu(2,j)*zi)
                  angaddy(iii) = angaddy(iii) + pmassj*
     &                 (vxyzu(1,j)*zi - xi*vxyzu(3,j))
                  angaddz(iii) = angaddz(iii) + pmassj*
     &                 (xi*vxyzu(2,j) - vxyzu(1,j)*yi)

                  xtemp = xtemp + pmassj*xyzmh(1,j)
                  ytemp = ytemp + pmassj*xyzmh(2,j)
                  ztemp = ztemp + pmassj*xyzmh(3,j)

                  xi = (xiold*pmassi + xtemp)/totalmass
                  yi = (yiold*pmassi + ytemp)/totalmass
                  zi = (ziold*pmassi + ztemp)/totalmass
                  vxi = (xmomsyn(iii) + xmomadd(iii))/totalmass
                  vyi = (ymomsyn(iii) + ymomadd(iii))/totalmass
                  vzi = (zmomsyn(iii) + zmomadd(iii))/totalmass
                  f1vxtemp = f1vxtemp + pmassj*f1vxyzu(1,j)
                  f1vytemp = f1vytemp + pmassj*f1vxyzu(2,j)
                  f1vztemp = f1vztemp + pmassj*f1vxyzu(3,j)
                  utemp = utemp + pmassj*vxyzu(4,j)
                  xyzmh(4,j) = pmassjn

 888              CONTINUE
               ENDIF
            ENDIF
         END DO
         xyzmh(4,i) = ptmsyn(iii) + ptmadd(iii)
         spinx(iii) = spinx(iii) + spinxtemp
         spiny(iii) = spiny(iii) + spinytemp
         spinz(iii) = spinz(iii) + spinztemp
         spinadx(iii) = spinadx(iii) + spinxtemp
         spinady(iii) = spinady(iii) + spinytemp
         spinadz(iii) = spinadz(iii) + spinztemp

         pmassnew = xyzmh(4,i)
         vxyzu(1,i) = (xmomsyn(iii) + xmomadd(iii))/pmassnew
         vxyzu(2,i) = (ymomsyn(iii) + ymomadd(iii))/pmassnew
         vxyzu(3,i) = (zmomsyn(iii) + zmomadd(iii))/pmassnew
         vxyzu(4,i) = (ui*pmassi + utemp)/pmassnew

         xyzmh(1,i) = (xiold*pmassi + xtemp)/pmassnew
         xyzmh(2,i) = (yiold*pmassi + ytemp)/pmassnew
         xyzmh(3,i) = (ziold*pmassi + ztemp)/pmassnew
         xyzmh(5,i) = hi
         f1vxyzu(1,i) = (f1vxi*pmassi + f1vxtemp)/pmassnew
         f1vxyzu(2,i) = (f1vyi*pmassi + f1vytemp)/pmassnew
         f1vxyzu(3,i) = (f1vzi*pmassi + f1vztemp)/pmassnew

         ENDIF
      END DO
c
c--MERGER OF TWO POINT MASSES
c
      imerge = 0

      GOTO 1000

      DO iii = 1, nptmass
         iptm1 = listpm(iii) 
         x1 = xyzmh(1,iptm1)
         y1 = xyzmh(2,iptm1)
         z1 = xyzmh(3,iptm1)
         DO jjj = iii + 1, nptmass
            iptm2 = listpm(jjj)
            rx = x1 - xyzmh(1,iptm2)
            ry = y1 - xyzmh(2,iptm2)
            rz = z1 - xyzmh(3,iptm2)
            r2 = rx*rx + ry*ry + rz*rz
            IF (r2.LT.(MAX(xyzmh(5,iptm1),xyzmh(5,iptm2))**2.0)) THEN
               imerge = 1
               WRITE(iprint,*) 'POINT MASSES MERGED, TIME=',realtime
               WRITE(iprint,*) '   Radius**2 ',r2
               WRITE(iprint,*) '   Point masses ',iii,jjj,
     &              iorig(iptm1),iorig(iptm2)
               i1list = iii
               i2list = jjj
               IF (xyzmh(4,iptm2).GT.xyzmh(4,iptm1)) THEN
                  itemp = iptm1
                  iptm1 = iptm2
                  iptm2 = itemp
                  itemp = i1list
                  i1list = i2list
                  i2list = itemp
                  rx = -rx
                  ry = -ry
                  rz = -rz
                  x1 = xyzmh(1,iptm1)
                  y1 = xyzmh(2,iptm1)
                  z1 = xyzmh(3,iptm1)
               ENDIF
               iphase(iptm2) = -1

               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR nlstacc'
                  CALL quit
               ENDIF
               listacc(nlstacc) = iptm2

               numberacc(i1list) = numberacc(i1list) + 1
               pmass1 = xyzmh(4,iptm1)
               pmass2 = xyzmh(4,iptm2)
               WRITE(iprint,*) '   Point masses ',iii,jjj,
     &              iorig(iptm1),iorig(iptm2)
               WRITE(iprint,*) '   Mases ',pmass1, pmass2
               dvx = xyzmh(1,iptm1) - xyzmh(1,iptm2)
               dvy = xyzmh(2,iptm1) - xyzmh(2,iptm2)
               dvz = xyzmh(3,iptm1) - xyzmh(3,iptm2)
               vx1 = vxyzu(1,iptm1)
               vy1 = vxyzu(2,iptm1)
               vz1 = vxyzu(3,iptm1)

               totalmass = pmass1 + pmass2
               spinx(i1list) = spinx(i1list) + spinx(i2list)
               spiny(i1list) = spiny(i1list) + spiny(i2list)
               spinz(i1list) = spinz(i1list) + spinz(i2list)
               spinadx(i1list) = spinadx(i1list) + spinadx(i2list)
               spinady(i1list) = spinady(i1list) + spinady(i2list)
               spinadz(i1list) = spinadz(i1list) + spinadz(i2list)
               spinm = pmass2*pmass1/totalmass
               spinx(i1list) = spinx(i1list) + spinm*(ry*dvz - dvy*rz)
               spiny(i1list) = spiny(i1list) + spinm*(dvx*rz - rx*dvz)
               spinz(i1list) = spinz(i1list) + spinm*(rx*dvy - dvx*ry)
               spinadx(i1list) = spinadx(i1list) + 
     &              spinm*(ry*dvz - dvy*rz)
               spinady(i1list) = spinady(i1list) + 
     &              spinm*(dvx*rz - rx*dvz)
               spinadz(i1list) = spinadz(i1list) + 
     &              spinm*(rx*dvy - dvx*ry)
               vxyzu(1,iptm1) = (pmass1*vxyzu(1,iptm1) + 
     &              pmass2*vxyzu(1,iptm2))/totalmass
               vxyzu(2,iptm1) = (pmass1*vxyzu(2,iptm1) + 
     &              pmass2*vxyzu(2,iptm2))/totalmass
               vxyzu(3,iptm1) = (pmass1*vxyzu(3,iptm1) + 
     &              pmass2*vxyzu(3,iptm2))/totalmass
               vxyzu(4,iptm1) = (pmass1*vxyzu(4,iptm1) + 
     &              pmass2*vxyzu(4,iptm2))/totalmass
               xyzmh(1,iptm1) = (pmass1*xyzmh(1,iptm1) +
     &              pmass2*xyzmh(1,iptm2))/totalmass
               xyzmh(2,iptm1) = (pmass1*xyzmh(2,iptm1) +
     &              pmass2*xyzmh(2,iptm2))/totalmass
               xyzmh(3,iptm1) = (pmass1*xyzmh(3,iptm1) +
     &              pmass2*xyzmh(3,iptm2))/totalmass

               f1vxyzu(1,iptm1) = (pmass1*f1vxyzu(1,iptm1) + 
     &              pmass2*f1vxyzu(1,iptm2))/totalmass
               f1vxyzu(2,iptm1) = (pmass1*f1vxyzu(2,iptm1) + 
     &              pmass2*f1vxyzu(2,iptm2))/totalmass
               f1vxyzu(3,iptm1) = (pmass1*f1vxyzu(3,iptm1) + 
     &              pmass2*f1vxyzu(3,iptm2))/totalmass
               numberacc(i1list) = numberacc(i1list)+numberacc(i2list)
               ptminner(i1list) = ptminner(i1list) + ptminner(i2list)

               ptmsyn(i1list) = ptmsyn(i1list) + ptmsyn(i2list)
               ptmadd(i1list) = ptmadd(i1list) + ptmadd(i2list)
               xyzmh(4,iptm1) = ptmsyn(i1list) + ptmadd(i1list)
               xyzmh(5,iptm1) = MAX(xyzmh(5,iptm1),xyzmh(5,iptm2))

               xmomsyn(i1list) = xyzmh(4,iptm1)*vxyzu(1,iptm1)
               ymomsyn(i1list) = xyzmh(4,iptm1)*vxyzu(2,iptm1)
               zmomsyn(i1list) = xyzmh(4,iptm1)*vxyzu(3,iptm1)
               xmomadd(i1list) = 0.0
               ymomadd(i1list) = 0.0
               zmomadd(i1list) = 0.0

               GOTO 900
            ENDIF
         END DO
      END DO
c
c--Compactify list of point masses if merger
c
 900  IF (imerge.EQ.1) THEN
         icount = 0
         iaccr = 1
         DO iii = 1, nptmass
            iptm = listpm(iii)
            IF (iphase(iptm).GE.1) THEN
               icount = icount + 1
               listpm(icount) = listpm(iii)
               nactotal(icount) = nactotal(iii)
               ptmassinner(icount) = ptmassinner(iii)
               spinx(icount) = spinx(iii)
               spiny(icount) = spiny(iii)
               spinz(icount) = spinz(iii)
               spinadx(icount) = spinadx(iii)
               spinady(icount) = spinady(iii)
               spinadz(icount) = spinadz(iii)

               ptmsyn(icount) = ptmsyn(iii)
               ptmadd(icount) = ptmadd(iii)
               xmomsyn(icount) = xmomsyn(iii)
               ymomsyn(icount) = ymomsyn(iii)
               zmomsyn(icount) = zmomsyn(iii)
               xmomadd(icount) = xmomadd(iii)
               ymomadd(icount) = ymomadd(iii)
               zmomadd(icount) = zmomadd(iii)

               nptlist(icount) = nptlist(iii)
               numberacc(icount) = numberacc(iii)
               ptminner(icount) = ptminner(iii)
               DO jjj = 1, nptlist(icount)
                  nearpt(jjj,icount) = nearpt(jjj,iii)
               END DO
            ENDIF
         END DO
         nptmass = icount
      ENDIF
c
c--Consider accreted particles' effect on ghosts
c
 1000 DO i = 1, nghost
         j = ireal(i + npart)
         IF (iremove(j).EQ.1 .OR. 
     &        (icreate.EQ.1 .AND. j.EQ.irhonex)) THEN
            iphase(npart + i) = -1
         ENDIF
      END DO
c
c--Reset each type of particle count
c             
C$OMP PARALLEL default(none)
C$OMP& shared(nlst0,llist,iremove,nptmass,nactotal,numberacc)
C$OMP& shared(ptmassinner,ptminner)
C$OMP& private(i)
C$OMP& reduction(+:naccrete)
C$OMP& reduction(-:nactive)

C$OMP DO SCHEDULE(runtime)
      DO i = 1, nptmass
         nactotal(i) = nactotal(i) + numberacc(i)
         naccrete = naccrete + numberacc(i)
         ptmassinner(i) = ptmassinner(i) + ptminner(i)
         nactive = nactive - numberacc(i)
      END DO
C$OMP END DO
c
c--Reset iremove to -1
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst0
         iremove(llist(i)) = -1
      END DO
C$OMP END DO
C$OMP END PARALLEL
c
c--Dump point mass details to ptprint file
c
      IF (isave.EQ.1 .OR. icreate.EQ.1) THEN
         tcomp = SQRT((3 * pi) / (32 * rhozero))
         WRITE (iptprint) realtime/tcomp, realtime, nptmass

         DO i = 1, nptmass
            j = listpm(i)
            IF (ibound.EQ.8 .OR. ibound.GE.90) THEN
               WRITE (iptprint)iorig(j),xyzmh(1,j),xyzmh(2,j),
     &              xyzmh(3,j),vxyzu(1,j),vxyzu(2,j),vxyzu(3,j),
     &              xyzmh(4,j),rho(j),nactotal(i),ptmassinner(i),
     &              spinx(i),spiny(i),spinz(i),angaddx(i),angaddy(i),
     &              angaddz(i),spinadx(i),spinady(i),spinadz(i),
     &              naccrete,anglostx,anglosty,anglostz,nkill
            ELSE
               WRITE (iptprint)iorig(j),xyzmh(1,j),xyzmh(2,j),
     &              xyzmh(3,j),vxyzu(1,j),vxyzu(2,j),vxyzu(3,j),
     &              xyzmh(4,j),rho(j),nactotal(i),ptmassinner(i),
     &              spinx(i),spiny(i),spinz(i),angaddx(i),angaddy(i),
     &              angaddz(i),spinadx(i),spinady(i),spinadz(i),
     &              naccrete
            ENDIF
            nactotal(i) = 0
            ptmassinner(i) = 0.
         END DO
      ENDIF

      IF (icreate.EQ.1) THEN
         iaccr = 1
         icreate = 0
      ENDIF

      RETURN
      END
      

