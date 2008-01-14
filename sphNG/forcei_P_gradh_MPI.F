      SUBROUTINE forcei(nlst_in,nlst_end,nlst_tot,listp,dt,itime,npart,
     &      ntot,xyzmh,vxyzu,fxyzu,dha,trho,pr,vsound,alphaMMpass,
     &      ekcle,dedxyz,Bxyz,dBxyz,Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the forces on particle ipart    *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif

      DIMENSION xyzmh(5,mmax), vxyzu(4,idim)
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
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/cylinder'
      INCLUDE 'COMMONS/compact'
c
c--Used for listparents list to keep a list of iupdated
c
      INCLUDE 'COMMONS/treecom_P'

      CHARACTER*7 where
      LOGICAL inosink

      DATA where/'forcei'/
      DATA epsil/1.E-2/
      DATA epsil2/1.E-4/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine forcei')

c      print *,iproc,': Entered forcei'

      nlistupdated = 0
      realtime = dt*itime/imaxstep + gt
#ifdef MPI
c
c--Need to get all sink particle data from all MPI processes so can check
c     that viscous and magnetic forces are not applied across sinks.
c     Need xyzmh only.  This is also done in accrete.f
c
      IF (nptmasstot.GT.0) THEN
         CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
         CALL MPI_TYPE_COMMIT(i5REAL8,ierr)

         DO ii = 1, nptmass
            listpm(ii) = listpm(ii) - 1
         END DO

         CALL MPI_TYPE_INDEXED(nptmass, lblocklengths, listpm,
     &        i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

         DO ii = 1, nptmass
            listpm(ii) = listpm(ii) + 1
         END DO

         nptmasslocal = nptmass
         DO ii = 1, numproc - 1
            iahead = MOD(iproc+ii,numproc)
            ibehind = MOD(numproc+iproc-ii,numproc)
c
c--Send number of sink particles to be received
c
            CALL MPI_SENDRECV(nptmass,1,MPI_INTEGER,iahead, 110,
     &           nptmassrec,1,MPI_INTEGER,ibehind, 110,
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
            IF (inumber.NE.1) THEN
               WRITE (*,*) 'ERROR - nptmassrec not received (forcei)'
               CALL quit
            ENDIF
#ifdef MPIDEBUGS
            print *,iproc,': FORCEI send ',nptmass,' ptmasses to ',
     &           iahead,' recv ',nptmassrec,' from ',ibehind
#endif
            IF (nptmasslocal+nptmassrec.GT.iptdim) THEN
             WRITE (*,*) 'ERROR iptdim not large enough for all sinks',
     &              nptmasslocal,nptmassrec,iptdim
               CALL quit
            ENDIF
            IF (2*ntot+2+inumbertotal+nptmasslocal+nptmassrec-
     &           nptmass.GT.mmax) THEN
               WRITE (*,*) 'ERROR mmax not large enough for all sinks',
     &              nptmasslocal,nptmassrec,nptmass,ntot,mmax
               CALL quit
            ENDIF
c
c--Send sink particle information (xyzmh, vxyzu, iscurrent, iunique)
c
            CALL MPI_SENDRECV(xyzmh, 1, indexMPI5, iahead, 111,
     &           xyzmh(1,2*ntot+2+inumbertotal+nptmasslocal-nptmass+1), 
     &           idim, i5REAL8, ibehind, 111, 
     &           MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
            IF (inumber.NE.nptmassrec) THEN
               WRITE (*,*) 'ERROR - inumber.NE.nptmassrec ',inumber,
     &              nptmassrec
               CALL quit
            ENDIF

            DO iii = nptmasslocal + 1, nptmasslocal + nptmassrec
               listpm(iii) = 2*ntot + 2 + inumbertotal + iii - nptmass
            END DO

            nptmasslocal = nptmasslocal + nptmassrec
         END DO
         IF (nptmasslocal.NE.nptmasstot) THEN
            WRITE (*,*) 'ERROR - nptmasslocal.NE.nptmasstot ',
     &           nptmasslocal,nptmasstot
            CALL quit
         ENDIF
#ifdef MPIDEBUGS
         print *,iproc,': FORCEI: Transfer of all sink data complete ',
     &      nptmasslocal,nptmass,nptmassglobal
#endif
         CALL MPI_TYPE_FREE(indexMPI5,ierr)
         CALL MPI_TYPE_FREE(i5REAL8,ierr)
      ENDIF
#endif
c
c--For external magnetic pressure boundaries
c
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
c--Calculate forces
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,npart,ntot,listp,nneigh,neimin,neimax)
C$OMP& shared(nlst_tot,icall,dt,itime,imaxstep,isteps)
C$OMP& shared(xyzmh,vxyzu,dha,fxyzu,trho,pr,vsound,dq,gradhs)
C$OMP& shared(neighb,neighover,dvtable,ddvtable,psoft)
C$OMP& shared(fmass,fpoten,part2kernel,part1kernel,radkernel)
C$OMP& shared(part2potenkernel,part1potenkernel,grwij)
C$OMP& shared(divv,curlv,beta,alpha,poten,dgrav,nlst0)
C$OMP& shared(cnormk,epsil,epsil2,where,pext)
C$OMP& shared(iphase,listpm,iprint,nptmasstot,iorig)
C$OMP& shared(alphaMMpass,alphamin)
C$OMP& shared(igrp,igphi,ifsvi,iexf)
C$OMP& shared(ifcor,iexpan,iener,damp)
C$OMP& shared(realtime,ekcle,encal,dedxyz,acc)
C$OMP& shared(Bxyz,dBxyz,Bextxi,Bextyi,Bextzi,B2ext,divcurlB,Bevolxyz)
C$OMP& shared(gravxyzstore,potenstore,varmhd,iupdated,listparents)
C$OMP& shared(iscurrent,ireal,ibound,nlistupdated,iprocdebug)
C$OMP& private(n,ipart,iparttree,stepsi)
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
C$OMP& private(gradhi,gradsofti,gradp,gradpi,gradpj)
C$OMP& private(Bxi,Byi,Bzi,B2i,Bevolxi,Bevolyi,Bevolzi)
C$OMP& private(Bxj,Byj,Bzj,B2j,projBi,projBj,dBx,dBy,dBz)
C$OMP& private(dBxideali,dByideali,dBzideali,dBxdissi,dBydissi,dBzdissi)
C$OMP& private(Beqxi,Beqyi,Beqzi,dBxiequil,dByiequil,dBziequil)
C$OMP& private(divBi,curlBxi,curlByi,curlBzi,fanisoxi,fanisoyi,fanisozi)
C$OMP& private(dsoftxi,dsoftyi,dsoftzi,vsigi,vsigj,alphaB)
C$OMP& private(vsproji,vsprojj,termb,vs2i,vs2j,dB2,robar1,fmi)
C$OMP& private(grkerni,grpmi,phii,dsoftterm,dsofttermi,projdB,projbext)
C$OMP& private(grkernj,phij,dsofttermj,rhoij1,fmj,vsoundj)
C$OMP& private(grkern,grkerntable,deulalpha,deulbeta)
C$OMP& private(irealj,advBx,advBy,advBz,termeta)
C$OMP& private(fextx,fexty,fextz,projBeqi,inosink)
C$OMP& private(fanisox,fanisoy,fanisoz,forcex,forcey,forcez,xmasi)
C$OMP& private(t12jpmassi,t12jpmassj,r32grav,forcexi,forceyi,forcezi)
C$OMP& reduction(+:ioutmin,ioutsup,ioutinf)
C$OMP& reduction(MIN:inmin,inminsy)
C$OMP& reduction(MAX:inmax,inmaxsy)

C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_tot
c         ipart = listp(n)
         ipart = ivar(3,n)

c         IF (ipart.GT.ntot) print *,iprocdebug,': D ',n

c
c--Needed for MPI code
c
         IF (ipart.GT.npart) THEN
            iparttree = ipart + ntot + 2
         ELSE
            iparttree = ipart
         ENDIF
      
ccc         CALL treef(ipart,npart,ntot,xyzmh,acc,0,fxi,fyi,fzi,poteni)

         forcexi = 0.
         forceyi = 0.
         forcezi = 0.

         fxyzu(4,ipart) = 0.
         dha(1,ipart) = 0.0

         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit
         ELSEIF (iphase(ipart).GE.1) THEN
            GOTO 80
         ENDIF

         IF (icall.EQ.3 .AND. n.LE.nlst_end) THEN
ccc            numneigh = nneigh(ipart)
            numneigh = ivar(1,n)
            inmin = MIN(inmin,numneigh)
            inmax = MAX(inmax,numneigh)
            inminsy = MIN(inminsy,numneigh)
            inmaxsy = MAX(inmaxsy,numneigh)
            IF (xyzmh(5,iparttree).LT.hmin .AND. numneigh.GT.neimin)
     &           ioutmin = ioutmin + 1
            IF (numneigh.GT.neimax) ioutsup = ioutsup + 1
            IF (numneigh.LT.neimin) ioutinf = ioutinf + 1
         ENDIF
c
c--Compute forces on particle ipart
c
         xi = xyzmh(1,iparttree)
         yi = xyzmh(2,iparttree)
         zi = xyzmh(3,iparttree)
         pmassi = xyzmh(4,iparttree)
         hi = xyzmh(5,iparttree)
         hi1 = 1./hi
         hi21 = hi1*hi1
         hi41 = hi21*hi21
         
c         IF (ipart.GT.ntot) print *,iprocdebug,': D2 ',n

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
            IF (iBext.GT.0 .AND. iexf.NE.9) THEN
               Beqxi = BEXTERNAL(xi,yi,zi,1)
               Beqyi = BEXTERNAL(xi,yi,zi,2)
               Beqzi = BEXTERNAL(xi,yi,zi,3)
            ELSE
               Beqxi = 0.
               Beqyi = 0.
               Beqzi = 0.
            ENDIF
            B2i = Bxi**2 + Byi**2 + Bzi**2
            dBxideali = 0.
            dByideali = 0.
            dBzideali = 0.
            dBxdissi = 0.
            dBydissi = 0.
            dBzdissi = 0.
            dBxiequil = 0.
            dByiequil = 0.
            dBziequil = 0.
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
c
c--note that pressure term includes isotropic magnetic pressure         
c         pro2i = (pr(ipart) - pext)*rho21i
         pro2i = max((pr(ipart) - pext + 0.5*(B2i-B2ext)),0.0)*rho21i
         vsoundi = vsound(ipart)
         vs2i = vsoundi**2 + B2i*rho1i

         stepsi = dt*isteps(ipart)/imaxstep

c         IF (ipart.GT.ntot) print *,iprocdebug,': D3 ',n
c
c--Loop over neighbors
c
ccc         DO 70 k = 1, nneigh(ipart)
ccc            j = neighlist(k)
         ioffset = ivar(2,n)
c         IF (iproc.EQ.0) print *,': ipart ',ipart,iunique(iorig(ipart)),
c     &        ivar(1,n)
         DO 70 k = 1, ivar(1,n)
            j = ijvar(ioffset + k)
            
c            IF (iunique(iorig(ipart)).EQ.2757) THEN
c            IF (iproc.EQ.0) THEN
c               print *,iproc,' N2757 ',iunique(iorig(j)),k,ivar(1,n),
c     &              ioffset,j,xyzmh(1,j),xyzmh(2,j),xyzmh(3,j)
c            ENDIF

c         IF (ipart.GT.ntot) print *,iprocdebug,': D4 ',n,k

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
            IF (.FALSE.) THEN
ccc            IF (.NOT.iscurrent(j) .AND. .NOT.iupdated(j)) THEN
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
                        Bxyz(1,j) = Bevolxyz(1,j)*trho(j)
                        Bxyz(2,j) = Bevolxyz(2,j)*trho(j)
                        Bxyz(3,j) = Bevolxyz(3,j)*trho(j)
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
                        Bxyz(1,j) = Bevolxyz(1,j)*trho(j)
                        Bxyz(2,j) = Bevolxyz(2,j)*trho(j)
                        Bxyz(3,j) = Bevolxyz(3,j)*trho(j)
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
#ifdef PERIODIC_NO_GHOSTS
            IF (imhd.EQ.idim .AND. varmhd.EQ.'eulr') THEN
               deulalpha = Bevolxi - Bevolxyz(1,j)
               deulbeta = Bevolyi - Bevolxyz(2,j)
               CALL modboundeulr(dx,dy,dz,deulalpha,deulbeta)
            ELSE
               CALL modbound(dx,dy,dz)
            ENDIF
#endif
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
               projdB = dBx*runix + dBy*runiy + dBz*runiz
               
               projBeqi = Beqxi*runix + Beqyi*runiy + Beqzi*runiz
               projBext = Bextxi*runix + Bextyi*runiy + Bextzi*runiz
c               projBrho2i = (projBi - projBext)*rho21i
c               projBrho2j = (projBj - projBext)*rho1j*rho1j
            ENDIF
c
c--Determine whether or not there is a sink particle between i and j
c
c--No artificial viscosity, resistivity, pressure, or MHD between particles 
c     across a point mass
c
            inosink = .TRUE.
            DO ii = 1, nptmasstot
               iptcurv = listpm(ii)
               xii = xyzmh(1,iptcurv)
               yii = xyzmh(2,iptcurv)
               zii = xyzmh(3,iptcurv)
               vpos = (xii-xi)*(xii-xyzmh(1,j)) + 
     &              (yii-yi)*(yii-xyzmh(2,j)) +
     &              (zii-zi)*(zii-xyzmh(3,j))
               IF (vpos.LT.0.0) THEN
                  inosink = .FALSE.
                  GOTO 65
               ENDIF
            END DO
 65         CONTINUE
c
c--Using hi
c
            dsoftterm = 0.
            IF (vi.LT.radkernel) THEN
               index = v2i*ddvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dxx = v2i - index*dvtable               
               dgrwdx = (grwij(index1)-grwij(index))*ddvtable ! slope
c
c--(note that kernel gradient is multiplied by gradhi)
c
               grkerntable = (grwij(index)+ dgrwdx*dxx)
               grkerni = grkerntable*hi41*gradhi
               grpmi = grkerni*pmassj
c
c--i contribution to pressure gradient and pdv
c
               gradp = grkerni*pro2i
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
                  dsoftterm = 0.5*grkerni*gradsofti*cnormk
               ENDIF

               IF (imhd.EQ.idim) THEN
c
c--time derivative terms
c
                  dBxideali = dBxideali - grpmi*dvx*projBi*cnormk
                  dByideali = dByideali - grpmi*dvy*projBi*cnormk
                  dBzideali = dBzideali - grpmi*dvz*projBi*cnormk
c
c--time derivative of equilibrium B field, iBext.ne.0
c
                  IF (iBext.GT.0 .AND. iexf.NE.9) THEN
                     dBxiequil = dBxiequil - grpmi*dvx*projBeqi*cnormk
                     dByiequil = dByiequil - grpmi*dvy*projBeqi*cnormk
                     dBziequil = dBziequil - grpmi*dvz*projBeqi*cnormk
                  ENDIF
c
c--compute divB
c
c                  divBi = divBi - grpmi*projdB*cnormk
                  divBi = divBi - weight*grkerntable*hi1*projdB*cnormk
c
c--compute current
c
                  curlBxi = curlBxi +grpmi*(dBy*runiz -dBz*runiy)*cnormk
                  curlByi = curlByi +grpmi*(dBz*runix -dBx*runiz)*cnormk
                  curlBzi = curlBzi +grpmi*(dBx*runiy -dBy*runix)*cnormk
c                  curlBxi = curlBxi + weight*grkerntable*hi1*cnormk*
c     &                               (dBy*runiz - dBz*runiy)
c                  curlByi = curlByi + weight*grkerntable*hi1*cnormk*
c     &                               (dBz*runix - dBx*runiz)
c                  curlBzi = curlBzi + weight*grkerntable*hi1*cnormk*
c     &                               (dBx*runiy - dBy*runix)
               ENDIF

c
c--If i and j are not across a sink particle
c
               IF (inosink) THEN
c
c--Add terms to particle i
c-------------------------
c
c--Pressure force
c
                  gradpi = gradp*pmassj*cnormk
                  gradxi = gradxi - gradpi*runix
c                  IF (iproc.EQ.0 .AND. ipart.EQ.1164) 
c     &                 print *,j,- gradpi*runix
                  gradyi = gradyi - gradpi*runiy
                  gradzi = gradzi - gradpi*runiz

                  grkern = 0.5*grkerni*cnormk
c
c--Anisotropic magnetic force (Morris formalism)
c
                  IF (imhd.EQ.idim) THEN
                     rhoij1 = rho1i*rho1j
                     projBext = Bextxi*runix +Bextyi*runiy +Bextzi*runiz
                     
                     fanisox = grkern*((Bxj-Bextxi)*(projBj-projBext)
     &                    - (Bxi-Bextxi)*(projBi-projBext))*rhoij1
                     fanisoy = grkern*((Byj-Bextyi)*(projBj-projBext)
     &                    - (Byi-Bextyi)*(projBi-projBext))*rhoij1
                     fanisoz = grkern*((Bzj-Bextzi)*(projBj-projBext)
     &                    - (Bzi-Bextzi)*(projBi-projBext))*rhoij1

                     fanisoxi = fanisoxi + pmassj*fanisox
                     fanisoyi = fanisoyi + pmassj*fanisoy
                     fanisozi = fanisozi + pmassj*fanisoz
c
c--Signal velocity (MHD)
c
                     vsoundj = vsound(j)
                     vs2j = vsoundj**2 + B2j*rho1j
                     vsproji = 2.*vsoundi*projBi*sqrtrho1i
                     vsprojj = 2.*vsoundj*projBj*SQRT(rho1j)
                     vsigi = 0.5*(SQRT(vs2i - vsproji)
     &                    + SQRT(vs2i + vsproji))
                     vsigj = 0.5*(SQRT(vs2j - vsprojj)
     &                    + SQRT(vs2j + vsprojj))
                     vsbar = 0.5*(vsigi + vsigj)
c
c--Artificial resistivity
c
                     IF (j.LE.npart .OR. ibound.EQ.11.OR.ibound.EQ.12) 
     &                    THEN
                        alphaB = 0.5*(alphaMMpass(2,ipart)
     &                       + alphaMMpass(2,j))
c                        alphaB = 1.0
c                        termB = alphaB*pmassj*grkern*
c     &                       MAX(vsbar - projv,0.0)*robar1
                        termB = alphaB*grkern*
     &                       (vsbar + 2.0*abs(projv))*robar1
                        termBpmj = termB*pmassj

                        IF (varmhd.EQ.'eulr') THEN
#ifndef PERIODIC_NO_GHOSTS
                           deulalpha = Bevolxi - Bevolxyz(1,j)
                           deulbeta = Bevolyi - Bevolxyz(2,j)
#endif
                           dBxdissi = dBxdissi + termBpmj*deulalpha
                           dBydissi = dBydissi + termBpmj*deulbeta
                        ELSE
                           dBxdissi = dBxdissi + termBpmj*dBx*robar1
                           dBydissi = dBydissi + termBpmj*dBy*robar1
                           dBzdissi = dBzdissi + termBpmj*dBz*robar1
c        dBxdissi = dBxdissi + termBpmj*(dBx - runix*projdB)*robar1
c        dBydissi = dBydissi + termBpmj*(dBy - runiy*projdB)*robar1
c        dBzdissi = dBzdissi + termBpmj*(dBz - runiz*projdB)*robar1
c
c-- Real constant resistivity
c
c        !times rhoi below
                           termeta = etamhd*grkern*rho1j*rho1i
                           termetampj = termeta*pmassj
                           dBxdissi = dBxdissi + termetapmj*2.*dBx
                           dBydissi = dBydissi + termetapmj*2.*dBy
                           dBzdissi = dBzdissi + termetapmj*2.*dBz

c                        dBxdissi = dBxdissi + termetapmj*
c     &                             (5.*runix*projdB-3.*dBx)
c                        dBydissi = dBydissi +
c     &                             (5.*runiy*projdB-3.*dBy)
c                        dBzdissi = dBzdissi +
c     &                             (5.*runiz*projdB-3.*dBz)
                        ENDIF
c
c--add contribution to thermal energy
c
                        dB2 = dBx*dBx + dBy*dBy + dBz*dBz
c                        dqi = dqi - termBpmj/cnormk*(dB2 - projdB**2)*robar1
                        dqi = dqi - termBpmj/cnormk*dB2*robar1
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
     &               (j.LE.npart.OR.ibound.EQ.11.OR.ibound.EQ.12)) THEN
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
                        t12j = grkern*f*(beta*f - alpha*vsbar)*robar1
                     ELSEIF (ifsvi.EQ.6) THEN
                        alphamean = (alphaMMpass(1,ipart) +
     &                       alphaMMpass(1,j))/2.0
                        t12j = alphamean*grkern*f*(2.0*f - vsbar)*robar1
                     ELSE
                        t12j = grkern*f*(beta*f - alpha*vsbar)*robar1
                     ENDIF

                     t12jpmassj = t12j*pmassj

c                  IF (iproc.EQ.0 .AND. ipart.EQ.1164) 
c     &            print *,' j ',j,t12j,projv,alphamean,grkern,f,robar1


                     artxi = artxi - t12jpmassj*runix
                     artyi = artyi - t12jpmassj*runiy
                     artzi = artzi - t12jpmassj*runiz
                     dqi = dqi + t12jpmassj*projv/cnormk
                  ELSE
                     t12j = 0.
                  ENDIF


c
c--Add terms to particle j (if j is active)
c-------------------------
c
                  IF (iscurrent(j)) THEN
c
c--Pressure force
c
                     IF (igrp.NE.0) THEN
                        gradpj = gradp*pmassi*cnormk
                        forcex = gradpj*runix
                        forcey = gradpj*runiy
                        forcez = gradpj*runiz
                     ELSE
                        forcex = 0.
                        forcey = 0.
                        forcez = 0.
                     ENDIF
c
c--Anisotropic magnetic force (Morris formalism)
c
                     IF (imhd.EQ.idim) THEN
                        forcex = forcex + pmassi*fanisox
                        forcey = forcey + pmassi*fanisoy
                        forcez = forcez + pmassi*fanisoz
c
c--Artificial resistivity
c
                        termBpmi = termB*pmassi

                        dBxdissj = 0.
                        dBydissj = 0.
                        dBzdissj = 0.

                        IF (varmhd.EQ.'eulr') THEN
#ifndef PERIODIC_NO_GHOSTS
                           deulalpha = Bevolxi - Bevolxyz(1,j)
                           deulbeta = Bevolyi - Bevolxyz(2,j)
#endif
                           dBxdissj = dBxdissj + termBpmi*deulalpha
                           dBydissj = dBydissj + termBpmi*deulbeta
                        ELSE
                           dBxdissj = dBxdissj + termBpmi*dBx*robar1
                           dBydissj = dBydissj + termBpmi*dBy*robar1
                           dBzdissj = dBzdissj + termBpmi*dBz*robar1
c
c-- Real constant resistivity
c
c        !times rhoi below
                           termetapmi = termeta*pmassi*2.
                           dBxdissj = dBxdissj + termetapmi*dBx
                           dBydissj = dBydissj + termetapmi*dBy
                           dBzdissj = dBzdissj + termetapmi*dBz
                        ENDIF
C$OMP ATOMIC
                        dBxyz(1,j) = dBxyz(1,j) + dBxdissj
C$OMP ATOMIC
                        dBxyz(2,j) = dBxyz(2,j) + dBydissj
C$OMP ATOMIC
                        dBxyz(3,j) = dBxyz(3,j) + dBzdissj
c
c--add contribution to thermal energy
c
C$OMP ATOMIC
                        dq(j) = dq(j) - termBpmi/cnormk*dB2*robar1
                     ENDIF
c
c--Artificial viscosity
c
                     IF (ifsvi.NE.0 .AND. projv.LT.0.0) THEN
                        t12jpmassi = t12j*pmassi
                        forcex = forcex + t12jpmassi*runix
                        forcey = forcey + t12jpmassi*runiy
                        forcez = forcez + t12jpmassi*runiz
C$OMP ATOMIC
                        dq(j) = dq(j) + t12jpmassi*projv/cnormk
                     ENDIF
c
c--Add contributions to particle j
c
C$OMP ATOMIC
                     fxyzu(1,j) = fxyzu(1,j) + forcex
C$OMP ATOMIC
                     fxyzu(2,j) = fxyzu(2,j) + forcey
C$OMP ATOMIC
                     fxyzu(3,j) = fxyzu(3,j) + forcez
                  ENDIF
c
c--End if for inosink
c
               ENDIF
            ELSE
               grkerni = 0.
               fmi = 1.
ccc               fmi = 0.
               phii = -rij1
            ENDIF
c
c--Gravitational force softening
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
                  fm = 0.5*fmi
                  phi = 0.5*phii
               ELSE
                  CALL error(where,1)
               ENDIF
c
c--Gravitational force calculation
c
               r32grav = 1.0/(rij2grav*rijgrav)
               IF (j.LE.npart) THEN
                  xmasj = fm*pmassj*r32grav
                  gravxi = gravxi - xmasj*dx


c                  IF (iproc.EQ.0 .AND. ipart.EQ.1164)
c     &                 print *,j,- xmasj*dx,' ipart '


                  gravyi = gravyi - xmasj*dy
                  gravzi = gravzi - xmasj*dz
                  poteni = poteni + phi*pmassj
               ENDIF
c
c--Gravitational softening from kernel softening (only non-zero if isoft=0)
c
               dsofttermi = dsoftterm*pmassj
               dsoftxi = dsoftxi - dsofttermi*runix
               dsoftyi = dsoftyi - dsofttermi*runiy
               dsoftzi = dsoftzi - dsofttermi*runiz
c
c--Add terms to particle j (if j is active)
c-------------------------
c
               IF (iscurrent(j)) THEN
                  xmasi = fm*pmassi*r32grav

                  IF (rij.GT.hj*radkernel) THEN
                     xmasi = xmasi - 0.5*pmassi*r32grav
                     phi = phi + 0.5*rij1

                     xmasj = 0.5*pmassj*r32grav
                     gravxi = gravxi - xmasj*dx
                     gravyi = gravyi - xmasj*dy
                     gravzi = gravzi - xmasj*dz
                     poteni = poteni - 0.5*rij1*pmassj



c                     print *,hi,hj,pmassi,r32grav,rij1,xmasi,phi,dx,
c     &                    dy,dz,ipart,j
                  ENDIF


c                  IF (iproc.EQ.0 .AND. j.EQ.1164)
c     &                 print *,j, xmasi*dx,' j-act ',ipart,fxyzu(1,j)


                  dsofttermi = dsoftterm*pmassi
C$OMP ATOMIC
                  fxyzu(1,j) = fxyzu(1,j) + xmasi*dx + dsofttermi*runix
C$OMP ATOMIC
                  fxyzu(2,j) = fxyzu(2,j) + xmasi*dy + dsofttermi*runiy
C$OMP ATOMIC
                  fxyzu(3,j) = fxyzu(3,j) + xmasi*dz + dsofttermi*runiz
C$OMP ATOMIC
                  poten(j) = poten(j) + phi*pmassi


c                  IF (iproc.EQ.0 .AND. j.EQ.1164)
c     &                 print *,'cumm ',fxyzu(1,j),xmasi*dx

               ENDIF
            ENDIF


c
c--Using hj  (only need to be done if j is non-active since if it is active
c     j's contribution is added to i when j is done.
c
            IF (.NOT.iscurrent(j)) THEN
               dsofttermj = 0.
               IF (vj.LT.radkernel) THEN
                  index = v2j*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx = v2j - index*dvtable
                  dgrwdx = (grwij(index1)-grwij(index))*ddvtable ! slope
c
c--(note that kernel gradient is multiplied by gradhj)
c
                  grkernj = (grwij(index)+ dgrwdx*dxx)*hj41*gradhs(1,j)
c
c--j contribution to pressure gradient and isotropic mag force
c
c               pro2j = (pr(j) - pext)*rho1j*rho1j
                  pro2j = max((pr(j) - pext + 0.5*(B2j-B2ext)),0.0)
     &                 *rho1j*rho1j
                  gradp = grkernj*pro2j
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
                     dsofttermj = 0.5*grkernj*pmassj*gradhs(2,j)*cnormk
                  ENDIF
c
c--If i and j are not across a sink particle
c
                  IF (inosink) THEN
c
c--Add terms to particle i
c-------------------------
c
c--Pressure force
c
                     gradpj = gradp*pmassj*cnormk
                     gradxi = gradxi - gradpj*runix
                     gradyi = gradyi - gradpj*runiy
                     gradzi = gradzi - gradpj*runiz

                     grkern = 0.5*grkernj*cnormk
c
c--Anisotropic magnetic force (Morris formalism)
c
                     IF (imhd.EQ.idim) THEN
                        rhoij1 = rho1i*rho1j
                        projBext =Bextxi*runix +Bextyi*runiy +
     &                       Bextzi*runiz

                        fanisox = grkern*((Bxj-Bextxi)*(projBj-projBext)
     &                       - (Bxi-Bextxi)*(projBi-projBext))*rhoij1
                        fanisoy = grkern*((Byj-Bextyi)*(projBj-projBext)
     &                       - (Byi-Bextyi)*(projBi-projBext))*rhoij1
                        fanisoz = grkern*((Bzj-Bextzi)*(projBj-projBext)
     &                       - (Bzi-Bextzi)*(projBi-projBext))*rhoij1

                        fanisoxi = fanisoxi + pmassj*fanisox
                        fanisoyi = fanisoyi + pmassj*fanisoy
                        fanisozi = fanisozi + pmassj*fanisoz
c
c--Signal velocity (MHD)
c
                        vsoundj = vsound(j)
                        vs2j = vsoundj**2 + B2j*rho1j
                        vsproji = 2.*vsoundi*projBi*sqrtrho1i
                        vsprojj = 2.*vsoundj*projBj*SQRT(rho1j)
                        vsigi = 0.5*(SQRT(vs2i - vsproji)
     &                       +SQRT(vs2i + vsproji))
                        vsigj = 0.5*(SQRT(vs2j - vsprojj)
     &                       +SQRT(vs2j + vsprojj))
                        vsbar = 0.5*(vsigi + vsigj)
c
c--Artificial resistivity
c
                        IF (j.LE.npart .OR. ibound.EQ.11.OR.ibound.EQ.12
     &                                                            ) THEN
                           alphaB = 0.5*(alphaMMpass(2,ipart)
     &                          + alphaMMpass(2,j))
c                        alphaB = 1.0
c                        termB = alphaB*pmassj*grkern*
c     &                     MAX(vsbar - projv,0.0)*robar1
                           termB = alphaB*pmassj*grkern*
     &                          (vsbar + 2.0*abs(projv))*robar1

                           IF (varmhd.EQ.'eulr') THEN
#ifndef PERIODIC_NO_GHOSTS
                              deulalpha = Bevolxi - Bevolxyz(1,j)
                              deulbeta = Bevolyi - Bevolxyz(2,j)
#endif
                              dBxdissi = dBxdissi + termB*deulalpha
                              dBydissi = dBydissi + termB*deulbeta
                           ELSE
                              dBxdissi = dBxdissi + termB*dBx*robar1
                              dBydissi = dBydissi + termB*dBy*robar1
                              dBzdissi = dBzdissi + termB*dBz*robar1
c               dBxdissi = dBxdissi + termB*(dBx - runix*projdB)*robar1
c               dBydissi = dBydissi + termB*(dBy - runiy*projdB)*robar1
c               dBzdissi = dBzdissi + termB*(dBz - runiz*projdB)*robar1
c
c--Real constant resistivity
c
c        !times rhoi below
                              termeta = etamhd*pmassj*grkern*rho1j*rho1i
                              dBxdissi = dBxdissi + termeta*2.*dBx
                              dBydissi = dBydissi + termeta*2.*dBy
                              dBzdissi = dBzdissi + termeta*2.*dBz

c                           dBxdissi = dBxdissi + termeta*
c     &                             (5.*runix*projdB-3.*dBx)
c                           dBydissi = dBydissi + 
c     &                             (5.*runiy*projdB-3.*dBy)
c                           dBzdissi = dBzdissi +                    
c     &                             (5.*runiz*projdB-3.*dBz)
                           ENDIF
c
c--Add contribution to thermal energy
c
                           dB2 = dBx*dBx + dBy*dBy + dBz*dBz
c                        dqi = dqi - termB*(dB2 - projdB**2)*robar1
                           dqi = dqi - termB*dB2*robar1 
                        ENDIF
                     ELSE
c
c--Signal velocity (hydro)
c
                        vsbar = 0.5*(vsoundi + vsound(j))
                     ENDIF
c
c--Artificial viscosity and energy dissipation
c
                     IF (ifsvi.NE.0 .AND. projv.LT.0. .AND. 
     &               (j.LE.npart.OR.ibound.EQ.11.OR.ibound.EQ.12)) THEN
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
                           t12j = grkern*f*(beta*f - alpha*vsbar)*robar1
                        ELSEIF (ifsvi.EQ.6) THEN
                           alphamean = (alphaMMpass(1,ipart) + 
     &                          alphaMMpass(1,j))/2.0
                        t12j = alphamean*grkern*f*(2.0*f - vsbar)*robar1
                        ELSE
                           t12j = grkern*f*(beta*f - alpha*vsbar)*robar1
                        ENDIF

                        t12jpmassj = t12j*pmassj
                        artxi = artxi - t12jpmassj*runix
                        artyi = artyi - t12jpmassj*runiy
                        artzi = artzi - t12jpmassj*runiz
                        dqi = dqi + t12jpmassj*projv/cnormk
                     ENDIF
c
c--End if for inosink
c
                  ENDIF
               ELSE
                  grkernj = 0.
                  fmj = 1.
ccc                  fmj = 0.
                  phij = -rij1
               ENDIF
c
c--Gravitational force softening
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
                     fm = 0.5*fmj
                     phi = 0.5*phij
                  ELSE
                     CALL error(where,1)
                  ENDIF
c
c--Gravitational softening from kernel softening (only non-zero if isoft=0)
c
                  dsoftxi = dsoftxi - dsofttermj*runix
                  dsoftyi = dsoftyi - dsofttermj*runiy
                  dsoftzi = dsoftzi - dsofttermj*runiz
c
c--Gravitational force calculation
c
                  IF (j.LE.npart) THEN
                     xmasj = fm*pmassj/(rij2grav*rijgrav)


c                  IF (iproc.EQ.0 .AND. ipart.EQ.1164)
c     &                 print *,j, - xmasj*dx,' j-inactive '


                     gravxi = gravxi - xmasj*dx
                     gravyi = gravyi - xmasj*dy
                     gravzi = gravzi - xmasj*dz
                     poteni = poteni + phi*pmassj
                  ENDIF
               ENDIF
            ENDIF
            


 70      CONTINUE
c
c--Store quantities
c
c         IF (ipart.GT.ntot) print *,iprocdebug,': D5 ',n

         IF (igrape.EQ.0 .AND. igphi.NE.0) THEN
            forcexi = forcexi + gravxi + dsoftxi
            forceyi = forceyi + gravyi + dsoftyi
            forcezi = forcezi + gravzi + dsoftzi
ccc            forcexi = gravxi
ccc            forceyi = gravyi
ccc            forcezi = gravzi

c            IF (iproc.EQ.0 .AND. ipart.EQ.1164) print *,'gravxi ',
c     &           gravxi

C$OMP ATOMIC
            poten(ipart) = poten(ipart) + poteni 
c
c--Add self contribution to potential
c
            IF (ipart.LE.npart) THEN
C$OMP ATOMIC
               poten(ipart)=poten(ipart)+xyzmh(4,ipart)*fpoten(0)*hi1
            ENDIF
            dgrav(ipart) = 0.
         ENDIF
c
c--Pressure gradients
c
         IF (igrp.NE.0) THEN
            forcexi = forcexi + gradxi
            forceyi = forceyi + gradyi
            forcezi = forcezi + gradzi

c            IF (ipart.EQ.1164) THEN
c               print *,iproc,ipart,gradxi,xyzmh(1,ipart),
c     &           xyzmh(2,ipart),xyzmh(3,ipart)
c            ENDIF
         ENDIF
c
c--Anisotropic magnetic force, div/curl B
c
         IF (imhd.EQ.idim) THEN
            forcexi = forcexi + fanisoxi
            forceyi = forceyi + fanisoyi
            forcezi = forcezi + fanisozi
            divcurlB(1,ipart) = divBi !!*rho1i
            curlBxi = curlBxi*rho1i
            curlByi = curlByi*rho1i
            curlBzi = curlBzi*rho1i
            divcurlB(2,ipart) = curlBxi
            divcurlB(3,ipart) = curlByi
            divcurlB(4,ipart) = curlBzi
c
c--get contributions to dBxyz and force from an external B field
c  (get force only using iexf=9, so do not add another contribution here)
c
            IF (iBext.GT.0 .AND. iexf.NE.9) THEN
               CALL externBfield(xi,yi,zi,hi,vxi,vyi,vzi,rhoi,
     &          Bxi,Byi,Bzi,curlBxi,curlByi,curlBzi,Beqxi,Beqyi,Beqzi,
     &          fextx,fexty,fextz,advBx,advBy,advBz,'all')
               forcexi = forcexi + fextx
               forceyi = forceyi + fexty
               forcezi = forcezi + fextz
            ELSE
               advBx = 0.
               advBy = 0.
               advBz = 0.
            ENDIF
c
c--time derivative of B/rho, B (constructed from B/rho and rho derivatives) 
c  or the Euler potentials (dissipation only)
c
            IF (varmhd(1:1).EQ.'B') THEN
C$OMP ATOMIC
               dBxyz(1,ipart) = dBxyz(1,ipart) +
     &                        ((dBxideali+dBxiequil)*rho21i + dBxdissi)
C$OMP ATOMIC
               dBxyz(2,ipart) = dBxyz(2,ipart) + 
     &                        ((dByideali+dByiequil)*rho21i + dBydissi)
C$OMP ATOMIC
               dBxyz(3,ipart) = dBxyz(3,ipart) +
     &                        ((dBzideali+dBziequil)*rho21i + dBzdissi)
               
               !--add contribution from external B field advection
C$OMP ATOMIC
               dBxyz(1,ipart) = dBxyz(1,ipart) - advBx*rho1i
C$OMP ATOMIC
               dBxyz(2,ipart) = dBxyz(2,ipart) - advBy*rho1i
C$OMP ATOMIC
               dBxyz(3,ipart) = dBxyz(3,ipart) - advBz*rho1i
               
               !--translate to dB/dt and add divv term if B is evolved
               IF (varmhd.EQ.'Bvol') THEN
                  dBxyz(1,ipart) = rhoi*dBxyz(1,ipart) 
     &                         + (Bxyz(1,ipart)+Beqxi)*rho1i*divv(ipart)
                  dBxyz(2,ipart) = rhoi*dBxyz(2,ipart) 
     &                         + (Bxyz(2,ipart)+Beqyi)*rho1i*divv(ipart)
                  dBxyz(3,ipart) = rhoi*dBxyz(3,ipart) 
     &                         + (Bxyz(3,ipart)+Beqzi)*rho1i*divv(ipart)
               ENDIF
            ELSE
C$OMP ATOMIC
               dBxyz(1,ipart) = dBxyz(1,ipart) + dBxdissi
C$OMP ATOMIC
               dBxyz(2,ipart) = dBxyz(2,ipart) + dBydissi
C$OMP ATOMIC
               dBxyz(3,ipart) = dBxyz(3,ipart) + dBzdissi
            ENDIF
         ENDIF
c
c--Artificial viscosity
c
c         IF (ipart.GT.ntot) print *,iprocdebug,': D6 ',n

         IF (ifsvi.NE.0) THEN
            forcexi = forcexi + artxi
            forceyi = forceyi + artyi
            forcezi = forcezi + artzi
         ENDIF
         fxyzu(4,ipart) = pdvi
         dq(ipart) = dq(ipart) + dqi
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
cc     &           MIN(divv(ipart)*rho1i+0.5*vsoundi/hi,0.0)
            dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*tdecay1 -
     &             MIN(divv(ipart)*rho1i,0.0) !!*fi
            IF (imhd.EQ.idim) THEN
               dB2 = divcurlB(1,ipart)**2 + divcurlB(2,ipart)**2
     &             + divcurlB(3,ipart)**2 + divcurlB(4,ipart)**2
               IF (B2i.GT.tiny) THEN
                  source = SQRT(MAX(dB2*rho1i,
     &                           vs2i*divcurlB(1,ipart)**2/B2i))
               ELSE
                  source = SQRT(dB2*rho1i)
               ENDIF
               dha(3,ipart) = (alphamin(2)-alphaMMpass(2,ipart))*tdecay1
     &                        + source
            ENDIF
         ENDIF
c
c--Damp velocities if appropiate
c
         IF (damp.NE.0.) THEN
            forcexi = forcexi - damp*vxyzu(1,ipart)
            forceyi = forceyi - damp*vxyzu(2,ipart)
            forcezi = forcezi - damp*vxyzu(3,ipart)
         ENDIF
c         IF (ipart.GT.ntot) print *,iprocdebug,': D7 ',n

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
         IF (iexf.GE.1) 
     &        CALL externf(ipart,ntot,realtime,xyzmh,forcexi,forceyi,
     &        forcezi,trho,iexf)
c
c--Coriolis and centrifugal forces
c
         IF (ifcor.NE.0) CALL coriol(ipart,ntot,realtime,xyzmh,vxyzu,
     &        forcexi,forceyi,forcezi) 
c
c--Homologous expansion or contraction
c
         IF (iexpan.GT.0) CALL homexp(ipart,realtime,vxyzu,forcexi,
     &        forceyi,forcezi) 

c
c--Add local values to global force
c
C$OMP ATOMIC
         fxyzu(1,ipart) = fxyzu(1,ipart) + forcexi
C$OMP ATOMIC
         fxyzu(2,ipart) = fxyzu(2,ipart) + forceyi
C$OMP ATOMIC
         fxyzu(3,ipart) = fxyzu(3,ipart) + forcezi

      END DO
C$OMP END DO
C$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_tot
         ipart = listp(n)
         IF (iphase(ipart).EQ.0) THEN
c
c--Energy conservation
c
c            IF (n.LE.nlst_end)
c            IF (iunique(iorig(ipart)).EQ.2757) THEN
c               print *,iproc,': F ',ipart,fxyzu(4,ipart),dq(ipart),
c     &              nlst_in,nlst_end,nlst_tot
c            ENDIF
         CALL energ(ipart,npart,ntot,realtime,vxyzu,fxyzu,xyzmh,trho)
c            IF (iunique(iorig(ipart)).EQ.2757) THEN
c               print *,iproc,': F2 ',ipart,fxyzu(4,ipart)
c            ENDIF
c
c--Translate to dB/dt and add divv term if B is evolved
c
            IF (varmhd(1:1).EQ.'B') THEN
               IF (varmhd.EQ.'Bvol') THEN
                  rhoi = trho(ipart)
                  rho1i = 1./rhoi

                  DO k = 1, 3
                     dBxyz(k,ipart) = rhoi*dBxyz(k,ipart)
                  END DO

                  IF (n.LE.nlst_end) THEN
                     IF (ipart.GT.npart) THEN
                        WRITE (*,*) 'ERROR - ipart.GT.npart'
                        CALL quit
                     ELSE
                        iparttree = ipart
                     ENDIF
                     xi = xyzmh(1,iparttree)
                     yi = xyzmh(2,iparttree)
                     zi = xyzmh(3,iparttree)

                     IF (iBext.GT.0 .AND. iexf.NE.9) THEN
                        Beqxi = BEXTERNAL(xi,yi,zi,1)
                        Beqyi = BEXTERNAL(xi,yi,zi,2)
                        Beqzi = BEXTERNAL(xi,yi,zi,3)
                     ELSE
                        Beqxi = 0.
                        Beqyi = 0.
                        Beqzi = 0.
                     ENDIF

                     dBxyz(1,ipart) = dBxyz(1,ipart)
     &                    + (Bxyz(1,ipart)+Beqxi)*rho1i*divv(ipart)
                     dBxyz(2,ipart) = dBxyz(2,ipart)
     &                    + (Bxyz(2,ipart)+Beqyi)*rho1i*divv(ipart)
                     dBxyz(3,ipart) = dBxyz(3,ipart)
     &                    + (Bxyz(3,ipart)+Beqzi)*rho1i*divv(ipart)

                  ENDIF
               ENDIF
            ENDIF
         ENDIF
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
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine forcei')

c      print *,iproc,': Exited forcei'


c      CALL quit

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
      INCLUDE 'COMMONS/cgas'
      
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

      trho(j) = rho(j) + deltarho
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
c
c--For radiative transfer, need to update extrapolated cv and kappa
c
      IF (encal.EQ.'r') THEN
         ekcle(3,j) = GETCV(trho(j),vxyzu(4,j))
         ekcle(2,j) = GETKAPPA(vxyzu(4,j),ekcle(3,j),trho(j))
      ENDIF

      RETURN
      
      END
