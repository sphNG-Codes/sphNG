      SUBROUTINE forcei(nlst_in,nlst_end,nlst_tot,listp,dt,itime,npart,
     &      ntot,xyzmh,vxyzu,fxyzu,dha,trho,pr,vsound,alphaMMpass,
     &      ekcle,dedxyz,Bxyz,dBxyz,Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the forces on particle ipart    *
c                                                           *
c************************************************************

#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION xyzmh(5,mmax2), vxyzu(4,idim2)
      REAL*4 trho(idim2),pr(idim2),vsound(idim2)
      REAL*4 dha(1+isizealphaMM,idim2),alphaMMpass(isizealphaMM,idim2)
      DIMENSION listp(idim2)
      DIMENSION fxyzu(4,idim3)
      DIMENSION ekcle(5,iradtrans2)
      DIMENSION dedxyz(3,iradtrans2)
      DIMENSION Bxyz(3,imhd2),dBxyz(imhdevol,imhd2)
      DIMENSION Bevolxyz(imhdevol,imhd3) ! needed for prediction only
c
c--Parameters for switching between different implementations of SPMHD.
c     FALSE is normal for both iMorris and iPhantom_resistivity
c
      LOGICAL iMorris
      PARAMETER (iMorris = .FALSE.)
      LOGICAL iPhantom_resistivity
      PARAMETER (iPhantom_resistivity = .FALSE.)
c
c--Whether or not to use an average h formalism for two of the MHD 
c     equations: the anisotropic magnetic force, and the induction 
c     equation.  All other equations continue to use h_i and h_j.
c     See Lewis, Bate & Price (2015, MNRAS) for more details.
c     NOTE: ihaverage also requires imhd = idim (i.e. MHD is used).
c
      LOGICAL ihaverage_switch
      PARAMETER (ihaverage_switch = .FALSE.)
      LOGICAL ihaverage
      PARAMETER (ihaverage = (ihaverage_switch .AND. imhd.EQ.idim))
c
c--If ihaverage is TRUE, then four different forms of average can be
c     used.  Standard is to use '0'.  See Lewis, Bate, Monaghan & Price
c     (2015, SPHERIC conference proceedings for more details).
c
c     0 => arithmetic
c     1 => geometric
c     2 => harmonic
c     3 => quadratic
c
      INTEGER ihavtype
      PARAMETER (ihavtype = 0)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/numpa'
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
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/pxpy'
c
c--Used for listparents list to keep a list of iupdated
c
      INCLUDE 'COMMONS/treecom_P'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif

#ifdef MPICOPY
      INCLUDE 'COMMONS/mpicopy'
#endif

      CHARACTER*7 where
      LOGICAL inosink

#ifdef _OPENMP
      INTEGER*8 forcei_lock,revtree_lock
      COMMON /locks / forcei_lock(idim2/10+1),revtree_lock
#endif

#ifdef NONIDEAL
      REAL*8 dBnonideal(3),jcb(3),jcbcb(3)
#endif
      REAL*4 dBnonideali(3),dBnonidealj(3)

      LOGICAL ifirst
      DATA ifirst/.TRUE./

      DATA where/'forcei'/
      DATA epsil/1.E-2/
      DATA epsil2/1.E-4/
      
      EXTERNAL etafunc
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine forcei')

#ifdef MPIALL
#ifdef MPIDEBUG
      print *,iproc,': Entered forcei ',itime
#endif
#endif
c
c---------- FIRST TIME AROUND ----------
c
      IF (ifirst) THEN
#ifdef MPIALL
         IF (iproc.EQ.0) THEN
#endif
            IF (ihaverage) THEN
               WRITE (iprint,251) ihavtype
 251           FORMAT (' MHD using h_average formalism: ihavtype = ',I1)
            ENDIF
#ifdef MPIALL
         ENDIF
#endif
         ifirst = .FALSE.
      ENDIF
c
c---------- END FIRST TIME AROUND ----------
c

#ifdef DGRID
      IF (gasdrag) dragqchange(:,:,:) = 0.0
#endif

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
               CALL quit(1)
            ENDIF
#ifdef MPIDEBUGS
            print *,iproc,': FORCEI send ',nptmass,' ptmasses to ',
     &           iahead,' recv ',nptmassrec,' from ',ibehind
#endif
            IF (nptmasslocal+nptmassrec.GT.iptdim) THEN
             WRITE (*,*) 'ERROR iptdim not large enough for all sinks',
     &              nptmasslocal,nptmassrec,iptdim
               CALL quit(1)
            ENDIF
            IF (2*ntot+2+inumbertotal+nptmasslocal+nptmassrec-
     &           nptmass+1.GT.mmax2) THEN
               WRITE (*,*) 'ERROR mmax2 not large enough for all sinks',
     &              nptmasslocal,nptmassrec,nptmass,ntot,mmax2
               CALL quit(1)
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
               CALL quit(1)
            ENDIF

            DO iii = nptmasslocal + 1, nptmasslocal + nptmassrec
               listpm(iii) = 2*ntot + 2 + inumbertotal + iii - nptmass
            END DO

            nptmasslocal = nptmasslocal + nptmassrec
         END DO
         IF (nptmasslocal.NE.nptmasstot) THEN
            WRITE (*,*) 'ERROR - nptmasslocal.NE.nptmasstot ',
     &           nptmasslocal,nptmasstot
            CALL quit(1)
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
C$OMP& shared(part2potenkernel,part1potenkernel,wij,grwij)
C$OMP& shared(divv,curlv,beta,alpha,poten,dgrav,nlst0)
C$OMP& shared(cnormk,doublehumpnormk,epsil,epsil2,where,pext)
C$OMP& shared(iphase,listpm,iprint,nptmasstot,iorig)
C$OMP& shared(alphaMMpass,alphamin,iproc)
C$OMP& shared(h2ratio,abhpq,abeq,abHIq,abco)
C$OMP& shared(igrp,igphi,ifsvi,iexf)
C$OMP& shared(ifcor,iexpan,iener,damp)
C$OMP& shared(realtime,ekcle,encal,dedxyz,acc,it1,imax,gradB)
C$OMP& shared(Bxyz,dBxyz,Bextxi,Bextyi,Bextzi,B2ext,divcurlB,Bevolxyz)
C$OMP& shared(varmhd,iupdated,listparents)
C$OMP& shared(iscurrent,ireal,ibound,nlistupdated,iprocdebug,irotpot)
C$OMP& shared(ivar,ijvar,Bextx,Bexty,Bextz,etamhd,iresist)
C$OMP& shared(pi)
C$OMP& shared(forcei_lock)
C$OMP& shared(gasdrag,iregime)
C$OMP& shared(planetesimaltimestep,planetesimalnorm)
C$OMP& shared(gasdist,nearestgas,idragscheme)
C$OMP& shared(itot)
#ifdef DGRID
C$OMP& shared(dragenergy,dragqchange,px,py,pz)
#endif
#ifdef NONIDEAL
C$OMP& shared(jcurrent,Zg_on_ag,n_electronT,gamma)
C$OMP& private(Bi,Bi1,Bj,Bj1,dBnonideal,jcb,jcbcb,temperature,ierr)
C$OMP& private(eta_ohmi,eta_halli,eta_ambii)
C$OMP& private(eta_ohmj,eta_hallj,eta_ambij)
#endif
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend,numbertodoherekeep)
C$OMP& shared(realtransfer5to15,llisttrans,idisplacementsi,irecvcounti)
C$OMP& shared(real4transfer2)
#endif
C$OMP& private(n,ipart,iparttree,stepsi,dtime)
C$OMP& private(xi,yi,zi,vxi,vyi,vzi,pmassi,dhi,hi,gravxi,gravyi,gravzi)
C$OMP& private(fxi,fyi,fzi,numneigh,hmin)
C$OMP& private(poteni,dphiti,gradxi,gradyi,gradzi,artxi,artyi,artzi)
C$OMP& private(gasdragxi,gasdragyi,gasdragzi,ig,f_drag_co,dvij)
C$OMP& private(rsub,fgd)
C$OMP& private(rho_g,f_drag_i,ptsi,ptsj,ptsnormi,ptsnormj)
C$OMP& private(pdvi,dqi,rhoi,pro2i,pro2j,vsoundi,k,j,hj,dx,dy,dz)
#ifdef DGRID
C$OMP& private(dragr,dragp,dragt,idragr,idragp,idragt)
#endif
C$OMP& private(dq_drag,dq_dragi,dq_dragj)
C$OMP& private(rho1i,rho21i,sqrtrho1i,rho1j,rho21j)
C$OMP& private(rij2,rij,rij1,pmassj,runix,runiy,runiz,hmean,hmean21)
C$OMP& private(hi1,hi21,hi31,hi41,v2i,vi)
C$OMP& private(hj1,hj21,hj31,hj41,v2j,vj)
C$OMP& private(hav1,hav21,hav41,v2av,vav,grkernav,hav)
C$OMP& private(index,dxx,index1,rij2grav,rijgrav,fm)
C$OMP& private(phi,dphi,dfmassdx,dfptdx,dpotdh,xmasj,rhoj,robar)
C$OMP& private(wtij,dwdx,dgrwdx,grwtij,grpm,dvx,dvy,dvz,projv,vsbar)
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
C$OMP& private(dBnonideali,dBnonidealj)
C$OMP& private(Beqxi,Beqyi,Beqzi,dBxiequil,dByiequil,dBziequil)
C$OMP& private(divBi,curlBxi,curlByi,curlBzi,fanisoxi,fanisoyi,fanisozi)
C$OMP& private(gradBxxi,gradBxyi,gradBxzi)
C$OMP& private(gradByxi,gradByyi,gradByzi)
C$OMP& private(gradBzxi,gradBzyi,gradBzzi)
C$OMP& private(dsoftxi,dsoftyi,dsoftzi,vsigi,vsigj,alphaB)
C$OMP& private(vsproji,vsprojj,termb,vs2i,vsi,vs2j,vsj,dB2,robar1,fmi)
C$OMP& private(grkerni,grpmi,phii,dsoftterm,dsofttermi,projdB,projbext)
C$OMP& private(grkernj,phij,dsofttermj,rhoij1,fmj,vsoundj)
C$OMP& private(grkern,grkerntable,dBevolx,dBevoly,dBevolz)
C$OMP& private(irealj,advBx,advBy,advBz,termeta)
C$OMP& private(fextx,fexty,fextz,projBeqi,inosink,Borve_temp)
C$OMP& private(fanisox,fanisoy,fanisoz,forcex,forcey,forcez,xmasi)
C$OMP& private(dBpsix,dBpsiy,dBpsiz,dBpsixj,dBpsiyj,dBpsizj)
C$OMP& private(t12jpmassi,t12jpmassj,r32grav,forcexi,forceyi,forcezi)
C$OMP& private(Bevoldotdv,termBpmi,termBpmj,termetapmj,termetapmi)
C$OMP& private(dBxidealj,dByidealj,dBzidealj)
C$OMP& private(dBxdissj,dBydissj,dBzdissj,dBxyz_xi,dBxyz_yi,dBxyz_zi)
C$OMP& private(ioffset,etai,etaj,ui,uj,dqj)
#ifdef MPICOPY
C$OMP& private(i,ierr)
#endif
C$OMP& reduction(+:ioutmin,ioutsup,ioutinf)
C$OMP& reduction(MIN:inmin,inminsy)
C$OMP& reduction(MAX:inmax,inmaxsy)

C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
      DO n = numberstart, numberend
#else
      DO n = nlst_in, nlst_tot
#endif
c         ipart = listp(n)
         ipart = ivar(3,n)
         
         if (ipart.eq.0) then
            print*,' ERROR : ipart= ',ipart,'n = ',n,nlst_in
            CALL quit(1)
         endif

c         IF (ipart.GT.idim) print *,iproc,': F ',n,ipart

c
c--Needed for MPI code
c
         IF (ipart.GT.npart) THEN
            iparttree = ipart + ntot + 2
         ELSE
            iparttree = ipart
         ENDIF
      
         forcexi = 0.
         forceyi = 0.
         forcezi = 0.

         fxyzu(4,ipart) = 0.
c         dha(1,ipart) = 0.0

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

         gasdragxi = 0.
         gasdragyi = 0.
         gasdragzi = 0.

         ptsi = 0.
         ptsnormi = 0.

         pdvi = 0.
         dqi = 0.
         dq_dragi = 0.
         
         ddvxi = 0.
         ddvyi = 0.
         ddvzi = 0.
         ddvscalar = 0.

         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit(1)
c
c--For sink particles do not need pressure, viscous & gravity softening forces
c
         ELSEIF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
            GOTO 80
         ENDIF

         IF (iphase(ipart).EQ.0 .AND. icall.EQ.3 .AND. 
     &        n.LE.nlst_end) THEN
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
         hi31 = hi21*hi1
         hi41 = hi21*hi21

#ifdef DGRID
         IF (gasdrag) THEN
            CALL sphericalpos(xi, yi, zi, dragr, dragp, dragt)
            CALL dragcoordinates(dragr, dragp, dragt, idragr, idragp,
     &           idragt)
         ENDIF
#endif

c         IF (ipart.GT.ntot) print *,iprocdebug,': D2 ',n

         vxi = vxyzu(1,ipart)
         vyi = vxyzu(2,ipart)
         vzi = vxyzu(3,ipart)
         ui = vxyzu(4,ipart)
         gradhi = gradhs(1,ipart)
         gradsofti = gradhs(2,ipart)

         rhoi = trho(ipart) ! set rhoi here so we can define eta
         rho1i = 1./rhoi
         sqrtrho1i = SQRT(rho1i)
         rho21i = rho1i*rho1i         

         IF (imhd.EQ.idim .AND. iphase(ipart).EQ.0) THEN
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
            dBnonideali = 0.
            dBxiequil = 0.
            dByiequil = 0.
            dBziequil = 0.
            divBi = 0.
            curlBxi = 0.
            curlByi = 0.
            curlBzi = 0.
            gradBxxi = 0.
            gradBxyi = 0.
            gradBxzi = 0.
            gradByxi = 0.
            gradByyi = 0.
            gradByzi = 0.
            gradBzxi = 0.
            gradBzyi = 0.
            gradBzzi = 0.
            fanisoxi = 0.
            fanisoyi = 0.
            fanisozi = 0.
            dBpsix = 0.
            dBpsiy = 0.
            dBpsiz = 0.
            IF (iresist.GE.2) THEN
               etai = etafunc(rhoi,ui,etamhd,iresist)
            ELSEIF (iresist.EQ.1) THEN
               etai = etamhd
            ELSE
               etai = 0.
            ENDIF
         ELSE
            Bxi = 0.
            Byi = 0.
            Bzi = 0.
            Bevolxi = 0.
            Bevolyi = 0.
            Bevolzi = 0.
            B2i = 0.
            B2j = 0.
            etai = 0.
         ENDIF
c
c--Note that pressure term includes isotropic magnetic pressure.
c     It also allows pressure to be negativ due to an external pressure
c     boundary (allowing the pressure to be negative may be necessary
c     e.g. for confining a turbulent cloud
c
c         pro2i = (pr(ipart) - pext)*rho21i
c         pro2i = max((pr(ipart) - pext + 0.5*(B2i-B2ext)),0.0)*rho21i
         pro2i = (pr(ipart) - pext + 0.5*(B2i-B2ext))*rho21i
         vsoundi = vsound(ipart)
         vs2i = vsoundi**2 + B2i*rho1i
c Note: vsi != vsigi; the latter includes a `projection' term
         vsi = SQRT(vs2i)

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

c         IF (ipart.GT.ntot) print *,iprocdebug,': D4 ',n,k

            IF (iphase(j).GE.1 .AND. iphase(j).LT.10) GOTO 70

            IF (iphase(j).EQ.-1) THEN
               WRITE(iprint,*)'ERROR - Accreted particle as neighbour!'
               WRITE(iprint,*) j,iorig(j),xyzmh(1,j),xyzmh(2,j),
     &              xyzmh(3,j),vxyzu(1,j),vxyzu(2,j),vxyzu(3,j)
               WRITE(iprint,*) ipart,iorig(ipart),icall,xi,yi,zi
               CALL quit(1)
            ENDIF
c
c--Gravity and potential energy
c
            dx = xi - xyzmh(1,j)
            dy = yi - xyzmh(2,j)
            dz = zi - xyzmh(3,j)
#ifdef PERIODIC_NO_GHOSTS
            IF (imhd.EQ.idim .AND, iphase(j).EQ.0) THEN
               IF (varmhd.EQ.'eulr') THEN
                  dBevolx = Bevolxi - Bevolxyz(1,j)
                  dBevoly = Bevolyi - Bevolxyz(2,j)
                  CALL modboundeulr(dx,dy,dz,dBevolx,dBevoly)
               ELSEIF (varmhd.EQ.'vecp') THEN
                  WRITE(iprint,*) 'modbound not implemented for vecp'
                  CALL quit(1)
               ENDIF
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
            hj31 = hj21*hj1
            hj41 = hj21*hj21

            IF (ihaverage) THEN
               IF (ihavtype.EQ.0) THEN
                  hav = 0.5*(hi+hj)
               ELSEIF (ihavtype.EQ.1) THEN
                  hav = SQRT(hi*hj)
               ELSEIF (ihavtype.EQ.2) THEN
                  hav = (2.0*hi*hj)/(hi+hj)
               ELSEIF (ihavtype.EQ.3) THEN
                  hav = SQRT( 0.5*(hi*hi+hj*hj) )
               ELSE
                  WRITE (iprint,99008) ihavtype
                  WRITE (*,99008) ihavtype
99008             FORMAT(' ERROR - Unimplemented ihavtype = ',I2)
                  CALL quit(0)
               ENDIF
               hav1 = 1./hav
               hav21 = hav1*hav1
               hav41 = hav21*hav21
            ELSE
               hav   = 0.
               hav1  = 0.
               hav21 = 0. 
               hav41 = 0.
            ENDIF
c
c--Added for gas-dust drag
c
            vsoundj = vsound(j)
#ifdef DGRID
            IF (gasdrag .AND. iphase(j).EQ.0) THEN
               CALL sphericalpos(xi, yi, zi, dragr, dragp, dragt)
               CALL dragcoordinates(dragr, dragp, dragt, idragr, idragp,
     &              idragt)
            ENDIF
#endif
c
c--Unit vectors
c
            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1
            dvx = vxi - vxyzu(1,j)
            dvy = vyi - vxyzu(2,j)
            dvz = vzi - vxyzu(3,j)
            uj = vxyzu(4,j)
            projv = dvx*runix + dvy*runiy + dvz*runiz
            dvij = SQRT(dvx**2 + dvy**2 + dvz**2)

            v2i = rij2*hi21
            vi = rij*hi1
            v2j = rij2*hj21
            vj = rij*hj1
            IF (ihaverage) THEN
               v2av = rij2*hav21
               vav = rij*hav1
            ELSE
               vav = 1.0E+10
               v2av = vav*vav
            ENDIF
            rhoj = trho(j)
            rho1j = 1./rhoj
            rho21j = rho1j*rho1j
            robar1 = 2./(rhoi + rhoj)

            IF (imhd.EQ.idim) THEN
#ifdef NONIDEAL
               dBnonideal = 0.
#endif
               IF (iphase(ipart).EQ.0 .AND. iphase(j).EQ.0) THEN
                  Bxj = Bxyz(1,j)
                  Byj = Bxyz(2,j)
                  Bzj = Bxyz(3,j)
#ifndef PERIODIC_NO_GHOSTS 
c if not done above
                  IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
                     dBevolx = Bevolxi - Bevolxyz(1,j)
                     dBevoly = Bevolyi - Bevolxyz(2,j)
                     dBevolz = Bevolzi - Bevolxyz(3,j)
                  ENDIF
#endif
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
                  IF (iresist.GE.2) THEN
                     etaj = etafunc(rhoj,uj,etamhd,iresist)
                  ELSEIF (iresist.EQ.1) THEN
                     etaj = etamhd
                  ELSE
                     etaj = 0.
                  ENDIF
               ELSE
                  Bxj = 0.
                  Byj = 0.
                  Bzj = 0.
                  dBx = 0.
                  dBy = 0.
                  dBz = 0.
                  B2j = 0.
                  projBi = 0.
                  projBj = 0.
                  projdB = 0.
                  etaj = 0.
               ENDIF
            ENDIF
c
c--Determine whether or not there is a sink particle between i and j
c
c--No artificial viscosity, resistivity, pressure, or MHD between particles 
c     across a point mass or planet core
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
            IF (iexf.EQ.7) THEN
               xii = 1.
               yii = 0.
               zii = 0.
               vpos = (xii-xi)*(xii-xyzmh(1,j)) +
     &              (yii-yi)*(yii-xyzmh(2,j)) +
     &              (zii-zi)*(zii-xyzmh(3,j))
               IF (vpos.LT.0.0) inosink = .FALSE.
            ENDIF
c
c--Using hi
c
            dsoftterm = 0.
            wtij = 0.
            IF (vi.LT.radkernel .OR. vav.LT.radkernel) THEN
               IF (vi.LT.radkernel) THEN
                  index = v2i*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx = v2i - index*dvtable               
                  dgrwdx = (grwij(index1)-grwij(index))*ddvtable ! slope
c
c--Kernel required for gas drag calculation.
c
                  IF (gasdrag) THEN
                     dwdx = (wij(index1) - wij(index))*ddvtable
                     wtij = (wij(index) + dwdx*dxx)*hi31*v2i*
     &                    doublehumpnormk
                  ENDIF
c
c--(note that kernel gradient is multiplied by gradhi)
c
                  grkerntable = (grwij(index)+ dgrwdx*dxx)*cnormk
                  grkerni = grkerntable*hi41*gradhi
               ELSE
                  grkerni = 0.
               ENDIF

               IF (vav.LT.radkernel) THEN
                  index = v2av*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx = v2av - index*dvtable               
                  dgrwdx = (grwij(index1)-grwij(index))*ddvtable ! slope
c
c--(note that haverage kernel gradient is NOT multiplied by gradhi)
c
                  grkerntable = (grwij(index)+ dgrwdx*dxx)*cnormk
                  grkernav = grkerntable*hav41
               ELSE
                  grkernav = 0.
               ENDIF

               grpmi = grkerni*pmassj
c
c--i contribution to pressure gradient and pdv
c
               gradp = grkerni*pro2i
c
c--Only allow gas work between two gas particles
c
               IF (iphase(ipart).EQ.0 .AND. iphase(j).EQ.0) THEN
                  pdvi = pdvi + grpmi*projv
               ENDIF
c
c--i contribution to force softening (including pseudo-pressure term)
c
               IF (isoft.EQ.0) THEN
                  IF (vi.LT.radkernel) THEN
                     dfmassdx = (fmass(index1) - fmass(index))*ddvtable
                     fmi = (fmass(index) + dfmassdx*dxx)
                     dfptdx = (fpoten(index1) - fpoten(index))*ddvtable
                     phii = (fpoten(index) + dfptdx*dxx)*hi1
                     IF (vi.GT.part2kernel) THEN
                        phii = phii + rij1*part2potenkernel
                     ELSEIF (vi.GT.part1kernel) THEN
                        phii = phii + rij1*part1potenkernel
                     ENDIF
                     dsoftterm = 0.5*grkerni*gradsofti
                  ELSE
                     fmi = 1.
                     phii = -rij1
                  ENDIF
               ENDIF

               IF (imhd.EQ.idim .AND. iphase(ipart).EQ.0 .AND.
     &                                           iphase(j).EQ.0) THEN
c
c--Zero j terms for safety
c
                  dBxidealj = 0.
                  dByidealj = 0.
                  dBzidealj = 0.
c
c--Time derivative terms of `ideal' B field
c
                  IF (varmhd.EQ.'vecp') THEN
c--Axel Gauge dA_i /dt = -A^j grad_i v_j
                   Bevoldotdv = Bevolxi*dvx + Bevolyi*dvy + Bevolzi*dvz
                   dBxideali = dBxideali + Bevoldotdv*grpmi*runix
                   dByideali = dByideali + Bevoldotdv*grpmi*runiy
                   dBzideali = dBzideali + Bevoldotdv*grpmi*runiz
c--Alternative Gauge dA_i /dt = v^j grad_i A_j
c                   Bevoldotdv = vxi*dBevolx + vyi*dBevoly + vzi*dBevolz
c                   dBxideali = dBxideali - Bevoldotdv*grpmi*runix
c                   dByideali = dByideali - Bevoldotdv*grpmi*runiy
c                   dBzideali = dBzideali - Bevoldotdv*grpmi*runiz
                  ELSE
c
c--Set first IF to disable haverage on induction equation
c                     IF(.TRUE.)THEN
                     IF (.NOT.ihaverage) THEN
c
c--Standard Borve form
c
                        dBxideali = dBxideali - grpmi*dvx*projBi
                        dByideali = dByideali - grpmi*dvy*projBi
                        dBzideali = dBzideali - grpmi*dvz*projBi
                     ELSE
c
c--haverage Borve formalism 
c     For the derivation of the inequalties, think Iolanthe
c     "And noble statesmen do not itch, To interfere with matters which,
c     They do not understand"
c
c--Do Part 1 unconditionally
c
                        dBxideali= dBxideali -grkernav*pmassj*dvx*projBi
                        dByideali= dByideali -grkernav*pmassj*dvy*projBi
                        dBzideali= dBzideali -grkernav*pmassj*dvz*projBi
c
c--If we can see j, but j can't see us, do j's field too
c     Do this here, not with the force in P2 because of inosink 
c     protective IF
c
                        IF (iscurrent(j).AND.vj.GT.radkernel) THEN
                           dBxidealj =       -grkernav*pmassi*dvx*projBj
                           dByidealj =       -grkernav*pmassi*dvy*projBj
                           dBzidealj =       -grkernav*pmassi*dvz*projBj
                      ENDIF 
                   ENDIF
                  ENDIF
c
c--time derivative of equilibrium B field, iBext.ne.0
c
                  IF (iBext.GT.0 .AND. iexf.NE.9) THEN
                     dBxiequil = dBxiequil - grpmi*dvx*projBeqi
                     dByiequil = dByiequil - grpmi*dvy*projBeqi
                     dBziequil = dBziequil - grpmi*dvz*projBeqi
                  ENDIF
c
c--compute divB
c
                  divBi = divBi - grpmi*projdB
c
c--compute current
c
                  curlBxi = curlBxi +grpmi*(dBy*runiz -dBz*runiy)
                  curlByi = curlByi +grpmi*(dBz*runix -dBx*runiz)
                  curlBzi = curlBzi +grpmi*(dBx*runiy -dBy*runix)
c
c--compute magnetic gradients
c
                  gradBxxi = gradBxxi - grpmi*dBx*runix
                  gradBxyi = gradBxyi - grpmi*dBx*runiy
                  gradBxzi = gradBxzi - grpmi*dBx*runiz
                  gradByxi = gradByxi - grpmi*dBy*runix
                  gradByyi = gradByyi - grpmi*dBy*runiy
                  gradByzi = gradByzi - grpmi*dBy*runiz
                  gradBzxi = gradBzxi - grpmi*dBz*runix
                  gradBzyi = gradBzyi - grpmi*dBz*runiy
                  gradBzzi = gradBzzi - grpmi*dBz*runiz
               ENDIF

c
c--If i and j are not across a sink particle
c
               IF (inosink) THEN
c
c--Next region only for gas-gas interactions
c
                  IF (iphase(ipart).EQ.0 .AND. iphase(j).EQ.0) THEN
c
c--Add terms to particle i
c-------------------------
c
c--Pressure force
c
                  gradpi = gradp*pmassj
                  gradxi = gradxi - gradpi*runix
                  gradyi = gradyi - gradpi*runiy
                  gradzi = gradzi - gradpi*runiz

                  grkern = 0.5*grkerni

                  IF (imhd.EQ.idim) THEN
                     IF (iMorris) THEN
c
c--Anisotropic magnetic force (Morris formalism)
c
                        rhoij1 = rho1i*rho1j
                        projBext = Bextxi*runix + Bextyi*runiy +
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
c--Set first IF to disable haverage on anisotropic MHD force equation.
c     In the non-haverage formalism with the Borve correction 
c     (subtracting off the full B div B term) there is no anisotropic
c     force term that involves W(hi), and the W(hj) term is done by
c     particle j and added to the force on i (i.e. in Part 2 below).
c                     ELSEIF(.FALSE.)THEN
                     ELSEIF (ihaverage) THEN 
c
c--Anisotropic magnetic force (Borve formalism) - haverage formalism
c     Do Part 1 unconditionally
c     Part 2 is inside the icurrent(j) IF
c
                         Borve_temp = pmassj*grkernav*projBj*rho1j*rho1j
                         fanisoxi = fanisoxi + Borve_temp*(Bxj-Bxi)
                         fanisoyi = fanisoyi + Borve_temp*(Byj-Byi)
                         fanisozi = fanisozi + Borve_temp*(Bzj-Bzi)
                     ENDIF
c
c--Evolution of psi/c_h
c
                     IF (imhdevol.EQ.4) THEN
                        dBpsix = dBpsix - pmassj*grkerni*vsi*
     &                       Bevolxyz(4,ipart)*rho1i*rho1i*runix
                        dBpsiy = dBpsiy - pmassj*grkerni*vsi*
     &                       Bevolxyz(4,ipart)*rho1i*rho1i*runiy
                        dBpsiz = dBpsiz - pmassj*grkerni*vsi*
     &                       Bevolxyz(4,ipart)*rho1i*rho1i*runiz
                     ENDIF
c
c--Signal velocity (MHD)
c
                     vsoundj = vsound(j)
                     vs2j = vsoundj**2 + B2j*rho1j
c Note: vsj != vsigj; the latter includes a `projection' term
                     vsj = SQRT(vs2j)
                     IF (iPhantom_resistivity) THEN
                        vsigi = SQRT(vs2i)
                        vsigj = SQRT(vs2j)
                     ELSE
                        vsproji = 2.*vsoundi*projBi*sqrtrho1i
                        vsprojj = 2.*vsoundj*projBj*SQRT(rho1j)
                        vsigi = 0.5*(SQRT(vs2i - vsproji)
     &                       + SQRT(vs2i + vsproji))
                        vsigj = 0.5*(SQRT(vs2j - vsprojj)
     &                       + SQRT(vs2j + vsprojj))
                     ENDIF
                     vsbar = 0.5*(vsigi + vsigj)
c
c--Artificial and Non-artificial resistivity
c
                     IF (j.LE.npart .OR. ibound.EQ.11.OR.ibound.EQ.12
     &                    .OR. ibound.EQ.104) THEN
                        IF (ifsvi.GE.6) THEN
                           alphaB = 0.5*(alphaMMpass(2,ipart)
     &                          + alphaMMpass(2,j))
                        ELSE
                           alphaB = alphamin(2)
                        ENDIF
c                        alphaB = 1.0
                        IF (iPhantom_resistivity) THEN
c
c--Phantom version:  (the 1/robar1 is to cancel that term in the dBxdissi 
c     equation below)
c
                           termB = alphaMMpass(2,ipart)*vsigi*rho21i*
     &                          grkern/robar1
                        ELSE
c                           termB = alphaB*pmassj*grkern*
c     &                           MAX(vsbar - projv,0.0)*robar1
                           termB = alphaB*grkern*
     &                           (vsbar + 2.0*abs(projv))*robar1
                        ENDIF

                        termBpmj = termB*pmassj
c
c--Physical resistivity
c
                        termeta = (etai+etaj)*grkern*rij1
                        termetapmj = termeta*pmassj*rho1j

                        IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
c--artificial
                           dBxdissi = dBxdissi + termBpmj*dBevolx
                           dBydissi = dBydissi + termBpmj*dBevoly
                           dBzdissi = dBzdissi + termBpmj*dBevolz
c--physical
                           dBxdissi = dBxdissi + termetapmj*dBevolx
                           dBydissi = dBydissi + termetapmj*dBevoly
                           dBzdissi = dBzdissi + termetapmj*dBevolz
                        ELSE
c--artificial
                           dBxdissi = dBxdissi + termBpmj*dBx*robar1
                           dBydissi = dBydissi + termBpmj*dBy*robar1
                           dBzdissi = dBzdissi + termBpmj*dBz*robar1
c        dBxdissi = dBxdissi + termBpmj*(dBx - runix*projdB)*robar1
c        dBydissi = dBydissi + termBpmj*(dBy - runiy*projdB)*robar1
c        dBzdissi = dBzdissi + termBpmj*(dBz - runiz*projdB)*robar1

c--physical
                           dBxdissi = dBxdissi + termetapmj*dBx*rho1i
                           dBydissi = dBydissi + termetapmj*dBy*rho1i
                           dBzdissi = dBzdissi + termetapmj*dBz*rho1i

c                        dBxdissi = dBxdissi + termetapmj*
c     &                             (5.*runix*projdB-3.*dBx)
c                        dBydissi = dBydissi +
c     &                             (5.*runiy*projdB-3.*dBy)
c                        dBzdissi = dBzdissi +
c     &                             (5.*runiz*projdB-3.*dBz)
#ifdef NONIDEAL
c--physical (Ohmic resistivity + Hall effect + ambipolar diffusion)
                           Bi = SQRT(B2i)
                           Bi1 = 1.0/Bi
                           temperature = nimhd_gastemp(encal
     &                                  ,vxyzu(4,ipart),ekcle(3,ipart)
     &                                  ,rhoi,gamma)
                           CALL nicil_get_eta(eta_ohmi,eta_halli
     &                         ,eta_ambii,Bi,rhoi,temperature
     &                         ,Zg_on_ag(ipart),n_electronT(ipart),ierr)
                           IF (ierr/=0) THEN
                              CALL nicil_translate_error(ierr)
                              CALL quit(1)
                           ENDIF
                           CALL nimhd_get_jcbcb(jcbcb,jcb
     &                         ,jcurrent(:,ipart),Bxi,Byi,Bzi,Bi1)
                           CALL nimhd_get_dBdt(dBnonideal
     &                         ,eta_ohmi,eta_halli,eta_ambii
     &                         ,jcurrent(:,ipart),jcb,jcbcb
     &                         ,runix,runiy,runiz)
                           dBnonideali = dBnonideali
     &                                 - dBnonideal
     &                                 * (pmassj*rho1i*rho1i)*grkerni
#endif
                        ENDIF
c
c--add contribution to thermal energy
c
                        dB2 = dBx*dBx + dBy*dBy + dBz*dBz
c                        dqi = dqi - termBpmj*(dB2 - projdB**2)*robar1
                        dqi = dqi - termBpmj*dB2*robar1
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
     &                 (j.LE.npart.OR.ibound.EQ.11.OR.ibound.EQ.12
     &                 .OR. ibound.EQ.104)) THEN
c
c--Calculate artificial viscosity:
c     If ifsvi=1 then normal viscosity
c     If ifsvi=2 then divv/curl weighted viscosity
c     If ifsvi=6 then viscosity switch in Morris and Monaghan, 1997,
c                  J. Comp. Phys. 136, 41-50.  This formulation does not use
c                  beta - it sets beta to be 2*alpha automatically.
c     If ifsvi=7 then uses MM97 switch for the alpha, but beta is fixed at 2
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
                     ELSEIF (ifsvi.EQ.7) THEN
                        alphamean = (alphaMMpass(1,ipart) +
     &                       alphaMMpass(1,j))/2.0
                        t12j = grkern*f*(2.0*f - alphamean*vsbar)*robar1
                     ELSE
                        t12j = grkern*f*(beta*f - alpha*vsbar)*robar1
                     ENDIF

                     t12jpmassj = t12j*pmassj

                     artxi = artxi - t12jpmassj*runix
                     artyi = artyi - t12jpmassj*runiy
                     artzi = artzi - t12jpmassj*runiz

                     dqi = dqi + t12jpmassj*projv
                  ELSE
                     t12j = 0.
                  ENDIF
c
c--This is the end of gas-gas interactions (i.e. iphase(i) and iphase(j) both 0)
c
                  ENDIF
c
c--Dust/Planetesimal drag terms with modified kernel.
c
                  IF (gasdrag .AND. (iphase(ipart).EQ.0 .AND.
     &                 iphase(j).EQ.11)) THEN

                     f_drag_co = getcd(rhoi, dvij, vsoundi,
     &                    ipart, idragscheme)

                     f_drag_i = f_drag_co*projv*wtij*pmassj

                     ptsi = ptsi +  pmassj*wtij/(rhoi*rhoi*f_drag_co/3.)
                     ptsnormi = ptsnormi +  pmassj*wtij/rhoi

                     gasdragxi = gasdragxi - f_drag_i*runix
                     gasdragyi = gasdragyi - f_drag_i*runiy
                     gasdragzi = gasdragzi - f_drag_i*runiz

                     dq_drag = 2.*f_drag_i*projv
                     dq_dragi = dq_dragi + dq_drag
                     dqi = dqi + dq_drag
                  ELSE
                     f_drag_co = 0.0
                  ENDIF


                  IF (gasdrag .AND.  (iphase(ipart).EQ.11 .AND.
     &                 iphase(j).EQ.0)) THEN
                     IF (rij.LT.gasdist(ipart)) THEN
                        gasdist(ipart) = rij
                        nearestgas(ipart) = ipart
                     ENDIF
                  ENDIF
c
c--Add terms to particle j (if j is active)
c-------------------------
c
                  IF (iscurrent(j)) THEN

                     forcex = 0.
                     forcey = 0.
                     forcez = 0.
                     dqj = 0.
                     dq_dragj = 0.
                     ptsj = 0.
                     ptsnormj = 0.
c
c--Next region only for gas-gas interactions
c
                  IF (iphase(ipart).EQ.0 .AND. iphase(j).EQ.0) THEN
c
c--Pressure force
c
                     IF (igrp.NE.0) THEN
                        gradpj = gradp*pmassi
                        forcex = forcex + gradpj*runix
                        forcey = forcey + gradpj*runiy
                        forcez = forcez + gradpj*runiz
                     ENDIF

                     IF (imhd.EQ.idim) THEN
c
c--Anisotropic magnetic force (Morris formalism)
c
                        IF (iMorris) THEN
                           forcex = forcex + pmassi*fanisox
                           forcey = forcey + pmassi*fanisoy
                           forcez = forcez + pmassi*fanisoz
                        ELSE
c
c--Anisotropic magnetic force (Borve formalism)
c
c--Set first IF to disable haverage on anisotropic MHD force equation
c                           IF(.TRUE.)THEN
                           IF (.NOT.ihaverage) THEN
c Normal
                              Borve_temp =
     &                             pmassi*grkerni*projBi*rho1i*rho1i
                           ELSEIF (vj.GT.radkernel) THEN
c haverage
                              Borve_temp =
     &                             pmassi*grkernav*projBi*rho1i*rho1i
                           ELSE
c hav but done in j's part 1
                              Borve_temp = 0.
                           ENDIF
                           forcex = forcex - Borve_temp*(Bxi-Bxj)
                           forcey = forcey - Borve_temp*(Byi-Byj)
                           forcez = forcez - Borve_temp*(Bzi-Bzj)
                        ENDIF
c
c--Evolution of psi/c_h
c
                        IF (imhdevol.EQ.4) THEN
                           dBpsixj = pmassi*grkerni*vsi*
     &                        Bevolxyz(4,ipart)*rho1i*rho1i*runix
                           dBpsiyj = pmassi*grkerni*vsi*
     &                        Bevolxyz(4,ipart)*rho1i*rho1i*runiy
                           dBpsizj = pmassi*grkerni*vsi*
     &                        Bevolxyz(4,ipart)*rho1i*rho1i*runiz
                        ENDIF
c
c--Artificial resistivity
c
                        termBpmi = termB*pmassi

                        dBxdissj = 0.
                        dBydissj = 0.
                        dBzdissj = 0.
                        dBnonidealj = 0.

                        termeta = (etai+etaj)*grkern*rij1
                        termetapmi = termeta*pmassi*rho1i

                        IF (varmhd.EQ.'eulr' .OR. varmhd.EQ.'vecp') THEN
c--artificial
                           dBxdissj = dBxdissj + termBpmi*dBevolx
                           dBydissj = dBydissj + termBpmi*dBevoly
                           dBydissj = dBydissj + termBpmi*dBevolz
c--physical
                           dBxdissj = dBxdissj + termetapmi*dBevolx
                           dBydissj = dBydissj + termetapmi*dBevoly
                           dBzdissj = dBzdissj + termetapmi*dBevolz
                        ELSE
c--artificial
                           dBxdissj = dBxdissj + termBpmi*dBx*robar1
                           dBydissj = dBydissj + termBpmi*dBy*robar1
                           dBzdissj = dBzdissj + termBpmi*dBz*robar1
c--physical
                           dBxdissj = dBxdissj + termetapmi*dBx*rho1j
                           dBydissj = dBydissj + termetapmi*dBy*rho1j
                           dBzdissj = dBzdissj + termetapmi*dBz*rho1j
#ifdef NONIDEAL
                           dBnonidealj = -dBnonideal
     &                                   * (pmassi*rho1i*rho1i)*grkerni
#endif
                        ENDIF
c
c--add contribution to thermal energy
c
                        dqj = dqj - termBpmi*dB2*robar1
                     ENDIF
c
c--Artificial viscosity
c
                     IF (ifsvi.NE.0 .AND. projv.LT.0.0) THEN
                        t12jpmassi = t12j*pmassi

                        forcex = forcex + t12jpmassi*runix
                        forcey = forcey + t12jpmassi*runiy
                        forcez = forcez + t12jpmassi*runiz

c                     IF (j.EQ.0418) THEN
c                        print *,' P0418dqj ',j,t12jpmassi,projv
c                     ENDIF

                        dqj = dqj + t12jpmassi*projv
                     ENDIF
c
c--This is the end of gas-gas interactions (i.e. iphase(i) and iphase(j) both 0)
c
                  ENDIF
c
c--Dust/Planetesimal terms (contribution if j particle is active).
c
                  IF (gasdrag .AND. (iphase(ipart).EQ.0 .AND.
     &                 iphase(j).EQ.11)) THEN

                     f_drag_i = f_drag_co*pmassi*projv*wtij

                     ptsj = ptsj + pmassi*wtij/(rhoi*rhoi*f_drag_co/3.)
                     ptsnormj = ptsnormj + pmassi*wtij/rhoi

                     forcex = forcex + f_drag_i*runix
                     forcey = forcey + f_drag_i*runiy
                     forcez = forcez + f_drag_i*runiz

                     dq_drag = 2.*f_drag_i*projv
                     dq_dragj = dq_dragj + dq_drag
                     dqj = dqj + dq_drag
                  ENDIF
c
c--Add contributions to particle j
c
c--Note: C$OMP ATOMIC would be the best way to do these IF the compiler was
c     written properly, but ifort isn't
c
#ifdef _OPENMP
                     CALL OMP_SET_LOCK(forcei_lock(j/10+1))
#endif
                     fxyzu(1,j) = fxyzu(1,j) + forcex
                     fxyzu(2,j) = fxyzu(2,j) + forcey
                     fxyzu(3,j) = fxyzu(3,j) + forcez

                     dq(j) = dq(j) + dqj

#ifdef DGRID
                     IF (gasdrag) THEN
C$OMP CRITICAL
                        dragqchange(idragr, idragp, idragt) = 
     &                       dragqchange(idragr, idragp, idragt) +
     &                       dq_dragj*dt*isteps(j)/imaxstep
C$OMP END CRITICAL
                     ENDIF
#endif
                     IF (gasdrag) THEN
                        planetesimaltimestep(j) =
     &                    planetesimaltimestep(j) + ptsj
                        planetesimalnorm(j) = 
     &                       planetesimalnorm(j) + ptsnormj
                     ENDIF

                     IF (imhd.EQ.idim) THEN
c
c--NOTE: dBxidealj = 0. unless haverage code
c
                        dBxyz(1,j) = dBxyz(1,j) -dBxdissj-dBnonidealj(1)
     &                             + dBxidealj*rho21j
                        dBxyz(2,j) = dBxyz(2,j) -dBydissj-dBnonidealj(2)
     &                             + dByidealj*rho21j
                        dBxyz(3,j) = dBxyz(3,j) -dBzdissj-dBnonidealj(3)
     &                             + dBzidealj*rho21j
                        IF (imhdevol.EQ.4) THEN
                           dBxyz(1,j) = dBxyz(1,j) + dBpsixj
                           dBxyz(2,j) = dBxyz(2,j) + dBpsiyj
                           dBxyz(3,j) = dBxyz(3,j) + dBpsizj
                        ENDIF
                     ENDIF
#ifdef _OPENMP
                     CALL OMP_UNSET_LOCK(forcei_lock(j/10+1))
#endif
                  ENDIF
c
c--End if for inosink
c
               ENDIF
            ELSE
               grkerni = 0.
               fmi = 1.
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
c------------------------------------------
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
                  ENDIF

                  dsofttermi = dsoftterm*pmassi
c
c--Note: C$OMP ATOMIC would be the best way to do these IF the compiler was
c     written properly, but ifort isn't
c
#ifdef _OPENMP
                  CALL OMP_SET_LOCK(forcei_lock(j/10+1))
#endif
                  fxyzu(1,j) = fxyzu(1,j) + xmasi*dx + dsofttermi*runix
                  fxyzu(2,j) = fxyzu(2,j) + xmasi*dy + dsofttermi*runiy
                  fxyzu(3,j) = fxyzu(3,j) + xmasi*dz + dsofttermi*runiz

                  poten(j) = poten(j) + phi*pmassi
#ifdef _OPENMP
                  CALL OMP_UNSET_LOCK(forcei_lock(j/10+1))
#endif
               ENDIF
            ENDIF
c
c--Using hj  (only need to be done if j is non-active since if it is active
c     j's contribution is added to i when j is done.
c
            IF (.NOT.iscurrent(j)) THEN

               dsofttermj = 0.
               wtij = 0.
               IF (vj.LT.radkernel .OR. vav.LT.radkernel) THEN
                  IF (vj.LT.radkernel) THEN
                     index = v2j*ddvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dxx = v2j - index*dvtable
                     dgrwdx = (grwij(index1)-grwij(index))*ddvtable 

                     IF (gasdrag .AND. iphase(ipart).EQ.11) THEN
                        dwdx = (wij(index1) - wij(index))*ddvtable
                        wtij = (wij(index)+dwdx*dxx)*hj31*v2j*
     &                       doublehumpnormk
                     ENDIF
c
c--(note that kernel gradient is multiplied by gradhj)
c
                     grkernj = (grwij(index)+ dgrwdx*dxx)*
     &                    cnormk*hj41*gradhs(1,j)
                  ELSE
                     grkernj = 0.
                  ENDIF

                  IF (vav.LT.radkernel) THEN
                     index = v2av*ddvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dxx = v2av - index*dvtable
                     dgrwdx = (grwij(index1)-grwij(index))*ddvtable 
c
c--(note that haverage kernel gradient is NOT multiplied by gradhj)
c
                     grkernav = (grwij(index)+ dgrwdx*dxx)*
     &                    cnormk*hav41
                  ELSE
                     grkernav = 0.
                  ENDIF
c
c--j contribution to pressure gradient and isotropic mag force
c
c                  pro2j = (pr(j) - pext)*rho1j*rho1j
c                  pro2j = max((pr(j) - pext + 0.5*(B2j-B2ext)),0.0)
c     &                 *rho1j*rho1j
                  pro2j = (pr(j) - pext + 0.5*(B2j-B2ext))*rho1j*rho1j
                  gradp = grkernj*pro2j
c
c--j contribution to force softening (including pseudo-pressure term)
c
                  IF (isoft.EQ.0 .OR. isoft.EQ.2) THEN
                     IF (vj.LT.radkernel) THEN
                        dfmassdx= (fmass(index1) -fmass(index))*ddvtable
                        fmj = (fmass(index) + dfmassdx*dxx)
                        dfptdx= (fpoten(index1) -fpoten(index))*ddvtable
                        phij = (fpoten(index) + dfptdx*dxx)*hj1
                        IF (vj.GT.part2kernel) THEN
                           phij = phij + rij1*part2potenkernel
                        ELSEIF (vj.GT.part1kernel) THEN
                           phij = phij + rij1*part1potenkernel
                        ENDIF
                        dsofttermj = 0.5*grkernj*pmassj*gradhs(2,j)
                     ELSE
                        fmj = 1.
                        phij = -rij1
                     ENDIF
                  ENDIF
c
c--If i and j are not across a sink particle
c
                  IF (inosink) THEN
c
c--Next region only for gas-gas interactions
c
                  IF (iphase(ipart).EQ.0 .AND. iphase(j).EQ.0) THEN
c
c--Add terms to particle i
c-------------------------
c
c--Pressure force
c
                     gradpj = gradp*pmassj
                     gradxi = gradxi - gradpj*runix
                     gradyi = gradyi - gradpj*runiy
                     gradzi = gradzi - gradpj*runiz

                     grkern = 0.5*grkernj
c
c--Anisotropic magnetic force (Morris formalism)
c
                     IF (imhd.EQ.idim) THEN
                        IF (iMorris) THEN
                           rhoij1 = rho1i*rho1j
                           projBext =Bextxi*runix +Bextyi*runiy +
     &                          Bextzi*runiz

                        fanisox = grkern*((Bxj-Bextxi)*(projBj-projBext)
     &                          - (Bxi-Bextxi)*(projBi-projBext))*rhoij1
                        fanisoy = grkern*((Byj-Bextyi)*(projBj-projBext)
     &                          - (Byi-Bextyi)*(projBi-projBext))*rhoij1
                        fanisoz = grkern*((Bzj-Bextzi)*(projBj-projBext)
     &                          - (Bzi-Bextzi)*(projBi-projBext))*rhoij1

                           fanisoxi = fanisoxi + pmassj*fanisox
                           fanisoyi = fanisoyi + pmassj*fanisoy
                           fanisozi = fanisozi + pmassj*fanisoz
c
c--Anisotropic magnetic force (Borve formalism)
c
                        ELSE
c
c--Set first IF to disable haverage on anisotropic MHD force equation
c                           IF(.TRUE.)THEN
                           IF (.NOT.ihaverage) THEN
c Normal
                              Borve_temp =
     &                             pmassj*grkernj*projBj*rho1j*rho1j
                           ELSE
c haverage but done in Part 1 so a No-op
                              Borve_temp = 0.
                           ENDIF
                           fanisoxi = fanisoxi + Borve_temp*(Bxj-Bxi)
                           fanisoyi = fanisoyi + Borve_temp*(Byj-Byi)
                           fanisozi = fanisozi + Borve_temp*(Bzj-Bzi)
                        ENDIF
c
c--Evolution of psi/c_h
c
                        IF (imhdevol.EQ.4) THEN
                           dBpsix = dBpsix - pmassj*grkernj*vsj*
     &                          Bevolxyz(4,j)*rho1j*rho1j*runix
                           dBpsiy = dBpsiy - pmassj*grkernj*vsj*
     &                          Bevolxyz(4,j)*rho1j*rho1j*runiy
                           dBpsiz = dBpsiz - pmassj*grkernj*vsj*
     &                          Bevolxyz(4,j)*rho1j*rho1j*runiz
                        ENDIF
c
c--Signal velocity (MHD)
c
                        vsoundj = vsound(j)
                        vs2j = vsoundj**2 + B2j*rho1j
                        IF (iPhantom_resistivity) THEN
                           vsigi = SQRT(vs2i)
                           vsigj = SQRT(vs2j)
                        ELSE
                           vsproji = 2.*vsoundi*projBi*sqrtrho1i
                           vsprojj = 2.*vsoundj*projBj*SQRT(rho1j)
                           vsigi = 0.5*(SQRT(vs2i - vsproji)
     &                          +SQRT(vs2i + vsproji))
                           vsigj = 0.5*(SQRT(vs2j - vsprojj)
     &                          +SQRT(vs2j + vsprojj))
                        ENDIF
                        vsbar = 0.5*(vsigi + vsigj)
c
c--Artificial resistivity
c
                        IF (j.LE.npart .OR. ibound.EQ.11.OR.ibound.EQ.12
     &                       .OR. ibound.EQ.104) THEN
                           IF (ifsvi.GE.6) THEN
                              alphaB = 0.5*(alphaMMpass(2,ipart)
     &                                  + alphaMMpass(2,j))
                           ELSE
                              alphaB = alphamin(2)
                           ENDIF
c                           alphaB = 1.0

                           IF (iPhantom_resistivity) THEN
c
c--Phantom version:  (the 1/robar1 is to cancel that term in the dBxdissi 
c     equation below)
c
                              termB = alphaMMpass(2,j)*vsigj*
     &                             rho1j*rho1j*grkern/robar1
                           ELSE
c                              termB = alphaB*pmassj*grkern*
c     &                             MAX(vsbar - projv,0.0)*robar1
                              termB = alphaB*pmassj*grkern*
     &                             (vsbar + 2.0*abs(projv))*robar1
                           ENDIF
c
c--Real resistivity
c
                           termeta = (etai+etaj)*
     &                               pmassj*grkern*rho1j*rij1

                           IF(varmhd.EQ.'eulr'.OR.varmhd.EQ.'vecp') THEN
c--artificial
                              dBxdissi = dBxdissi + termB*dBevolx
                              dBydissi = dBydissi + termB*dBevoly
                              dBzdissi = dBzdissi + termB*dBevolz
c--physical
                              dBxdissi = dBxdissi + termeta*dBevolx
                              dBydissi = dBydissi + termeta*dBevoly
                              dBzdissi = dBzdissi + termeta*dBevolz
                           ELSE
c--artificial
                              dBxdissi = dBxdissi + termB*dBx*robar1
                              dBydissi = dBydissi + termB*dBy*robar1
                              dBzdissi = dBzdissi + termB*dBz*robar1
c               dBxdissi = dBxdissi + termB*(dBx - runix*projdB)*robar1
c               dBydissi = dBydissi + termB*(dBy - runiy*projdB)*robar1
c               dBzdissi = dBzdissi + termB*(dBz - runiz*projdB)*robar1

c--physical
                              dBxdissi = dBxdissi + termeta*dBx*rho1i
                              dBydissi = dBydissi + termeta*dBy*rho1i
                              dBzdissi = dBzdissi + termeta*dBz*rho1i

c                           dBxdissi = dBxdissi + termeta*
c     &                             (5.*runix*projdB-3.*dBx)
c                           dBydissi = dBydissi + 
c     &                             (5.*runiy*projdB-3.*dBy)
c                           dBzdissi = dBzdissi +                    
c     &                             (5.*runiz*projdB-3.*dBz)
#ifdef NONIDEAL
c--physical (Ohmic resistivity + Hall effect + ambipolar diffusion)
                              Bj  = SQRT(B2j)
                              Bj1 = 1.0/Bj
                              temperature = nimhd_gastemp(encal
     &                                     ,vxyzu(4,j),ekcle(3,j),rhoj
     &                                     ,gamma)
                              CALL nicil_get_eta(eta_ohmj,eta_hallj
     &                            ,eta_ambij,Bj,rhoj,temperature
     &                            ,Zg_on_ag(j),n_electronT(j),ierr)
                              IF (ierr/=0) THEN
                                 CALL nicil_translate_error(ierr)
                                 CALL quit(1)
                              ENDIF
                              CALL nimhd_get_jcbcb(jcbcb,jcb
     &                           ,jcurrent(:,j),Bxj,Byj,Bzj,Bj1)
                              CALL nimhd_get_dBdt(dBnonideal
     &                            ,eta_ohmj,eta_hallj,eta_ambij
     &                            ,jcurrent(:,j),jcb,jcbcb
     &                            ,runix,runiy,runiz)
                              dBnonideali = dBnonideali
     &                                    - dBnonideal
     &                                    * (pmassj*rho1j*rho1j)*grkernj

#endif
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
     &                    (j.LE.npart.OR.ibound.EQ.11.OR.ibound.EQ.12
     &                    .OR. ibound.EQ.104)) THEN
c
c--Calculate artificial viscosity:
c     If ifsvi=1 then normal viscosity
c     If ifsvi=2 then divv/curl weighted viscosity
c     If ifsvi=6 then viscosity switch in Morris and Monaghan, 1997, 
c                  J. Comp. Phys. 136, 41-50.  This formulation does not use
c                  beta - it sets beta to be 2*alpha automatically.
c     If ifsvi=7 then uses MM97 switch for the alpha, but beta is fixed at 2
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
                        ELSEIF (ifsvi.EQ.7) THEN
                           alphamean = (alphaMMpass(1,ipart) + 
     &                          alphaMMpass(1,j))/2.0
                        t12j = grkern*f*(2.0*f - alphamean*vsbar)*robar1
                        ELSE
                           t12j = grkern*f*(beta*f - alpha*vsbar)*robar1
                        ENDIF

                        t12jpmassj = t12j*pmassj
                        artxi = artxi - t12jpmassj*runix
                        artyi = artyi - t12jpmassj*runiy
                        artzi = artzi - t12jpmassj*runiz

c                     IF (ipart.EQ.0418) THEN
c                        print *,' P0418dqi ',j,t12jpmassj,projv
c                     ENDIF

                        dqi = dqi + t12jpmassj*projv
                     ENDIF
c
c--This is the end of gas-gas interactions (i.e. iphase(i) and iphase(j) both 0)
c
                  ENDIF
c
c--Dust/Planetesimal terms (contribution if j particle is inactive).
c
                  IF (gasdrag .AND. (iphase(ipart).EQ.11 .AND.
     &                 iphase(j).EQ.0)) THEN

                     f_drag_co = getcd(rhoj, dvij, vsoundj,
     &                    ipart, idragscheme)

                     f_drag_i = f_drag_co*projv*wtij*pmassj

                     ptsi = ptsi +  pmassj*wtij/
     &                      (rhoj*rhoj*f_drag_co/3.)
                     ptsnormi = ptsnormi +  pmassj*wtij/rhoj

                     gasdragxi = gasdragxi - f_drag_i*runix
                     gasdragyi = gasdragyi - f_drag_i*runiy
                     gasdragzi = gasdragzi - f_drag_i*runiz

                     dq_drag = 2.*f_drag_i*projv
                     dq_dragi = dq_dragi + dq_drag
                     dqi = dqi + dq_drag
                  ENDIF
c
c--End if for inosink
c
                  ENDIF
               ELSE
                  grkernj = 0.
                  fmj = 1.
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
c
c--Add self contribution to potential
c
            IF (ipart.LE.npart) THEN
               poteni = poteni+xyzmh(4,ipart)*fpoten(0)*hi1
               dgrav(ipart) = 0.
            ENDIF
         ENDIF
c
c--Add gas-dust drag interaction force
c
         IF (gasdrag) THEN
            forcexi = forcexi + gasdragxi
            forceyi = forceyi + gasdragyi
            forcezi = forcezi + gasdragzi
         ENDIF
c
c--The force contributions below only apply to gas particles
c
         IF (iphase(ipart).EQ.0) THEN
c
c--Pressure gradients
c
         IF (igrp.NE.0) THEN
            forcexi = forcexi + gradxi
            forceyi = forceyi + gradyi
            forcezi = forcezi + gradzi
         ENDIF
c
c--Anisotropic magnetic force, div/curl B
c
         IF (imhd.EQ.idim) THEN
            forcexi = forcexi + fanisoxi
            forceyi = forceyi + fanisoyi
            forcezi = forcezi + fanisozi

c            IF (ipart.LE.npart) THEN
               divcurlB(1,ipart) = divBi *rho1i
               curlBxi = curlBxi*rho1i
               curlByi = curlByi*rho1i
               curlBzi = curlBzi*rho1i
               divcurlB(2,ipart) = curlBxi
               divcurlB(3,ipart) = curlByi
               divcurlB(4,ipart) = curlBzi

               gradB(1,ipart) = gradBxxi*rho1i
               gradB(2,ipart) = gradBxyi*rho1i
               gradB(3,ipart) = gradBxzi*rho1i
               gradB(4,ipart) = gradByxi*rho1i
               gradB(5,ipart) = gradByyi*rho1i
               gradB(6,ipart) = gradByzi*rho1i
               gradB(7,ipart) = gradBzxi*rho1i
               gradB(8,ipart) = gradBzyi*rho1i
               gradB(9,ipart) = gradBzzi*rho1i
c            ENDIF
c
c--get contributions to dBxyz and force from an external B field
c  (get force only using iexf=9, so do not add another contribution here)
c
            IF (iBext.GT.0 .AND. iexf.NE.9 .AND. ipart.LE.npart) THEN
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
            dBxyz_xi = 0.
            dBxyz_yi = 0.
            dBxyz_zi = 0.
            IF (varmhd(1:1).EQ.'B') THEN
               dBxyz_xi = dBxyz_xi +
     &                        ((dBxideali+dBxiequil)*rho21i 
     &                        + dBxdissi + dBnonideali(1))
               dBxyz_yi = dBxyz_yi + 
     &                        ((dByideali+dByiequil)*rho21i 
     &                        + dBydissi + dBnonideali(2))
               dBxyz_zi = dBxyz_zi +
     &                        ((dBzideali+dBziequil)*rho21i 
     &                        + dBzdissi + dBnonideali(3))
               
               !--add contribution from external B field advection
               dBxyz_xi = dBxyz_xi - advBx*rho1i
               dBxyz_yi = dBxyz_yi - advBy*rho1i
               dBxyz_zi = dBxyz_zi - advBz*rho1i

            ELSEIF (varmhd.EQ.'vecp') THEN
               dBxyz_xi= dBxyz_xi +dBxideali*rho1i +dBxdissi
     &                       + (vyi*Bextz - vzi*Bexty)  ! add v x B_ext
               dBxyz_yi= dBxyz_yi +dByideali*rho1i +dBydissi
     &                       + (vzi*Bextx - vxi*Bextz)  ! add v x B_ext
               dBxyz_zi= dBxyz_zi +dBzideali*rho1i +dBzdissi
     &                       + (vxi*Bexty - vyi*Bextx)  ! add v x B_ext
            ELSE ! Euler potentials
               dBxyz_xi = dBxyz_xi + dBxdissi
               dBxyz_yi = dBxyz_yi + dBydissi
               dBxyz_zi = dBxyz_zi + dBzdissi
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
c
c--Morris & Monaghan switch source and decay terms for both artificial viscosity
c  and artificial resistivity (see Price & Monaghan 2005)
c
         IF (ifsvi.GE.6) THEN
            vsigi = SQRT(vs2i)
            tdecay1 = 0.2*vsigi/hi
c--balsara factor
c            adivi = ABS(divv(ipart)*rho1i)
c            acurlvi = ABS(curlv(ipart)*rho1i)
c            fi = adivi/(adivi+acurlvi+epsil2*vsigi/hi)

cc         dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*tdecay1 -
cc     &           MIN(divv(ipart)*rho1i+0.5*vsoundi/hi,0.0)
c            dha(2,ipart) = (alphamin(1)-alphaMMpass(1,ipart))*tdecay1 -
c     &             MIN(divv(ipart)*rho1i,0.0) !!*fi
            IF (imhd.EQ.idim .AND. ipart.LE.ntot) THEN
               dB2 = divcurlB(1,ipart)**2 + divcurlB(2,ipart)**2
     &             + divcurlB(3,ipart)**2 + divcurlB(4,ipart)**2
               IF (B2i.GT.tiny) THEN
                  source = SQRT(MAX(dB2*rho1i,
     &                           vs2i*divcurlB(1,ipart)**2/B2i))
               ELSE
                  source = SQRT(dB2*rho1i)
               ENDIF
c               dha(3,ipart) = (alphamin(2)-alphaMMpass(2,ipart))*tdecay1
c     &                        + source
            ENDIF
         ENDIF
c
c--Don't need these forces for foreign MPI particles
c
         IF (ipart.LE.ntot) THEN
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
               forcexi = forcexi -
     &              ekcle(4,ipart)/trho(ipart)*dedxyz(1,ipart)
               forceyi = forceyi -
     &              ekcle(4,ipart)/trho(ipart)*dedxyz(2,ipart)
               forcezi = forcezi -
     &              ekcle(4,ipart)/trho(ipart)*dedxyz(3,ipart)
            ENDIF
         ENDIF
c
c--END for: force contributions that only apply to gas particles
c
         ENDIF

c
c--Forces that operate on both sink particles and gas from here on:
c
 80      CONTINUE

c
c--External forces
c       
         IF (iexf.GE.1) 
     &        CALL externf(ipart,ntot,realtime,xyzmh,vxyzu,forcexi,
     &        forceyi,forcezi,trho,iexf,ibound,irotpot)
c
c--Coriolis and centrifugal forces
c
         IF (ifcor.NE.0) CALL coriol(ipart,ntot,realtime,xyzmh,vxyzu,
     &        forcexi,forceyi,forcezi) 
c
c--Homologous expansion or contraction
c
         IF (iexpan.GT.0) CALL homexp(ipart,ntot,realtime,vxyzu,forcexi,
     &        forceyi,forcezi) 

c
c--Add local values to global force
c
c--Note: C$OMP ATOMIC would be the best way to do these IF the compiler was
c     written properly, but ifort isn't
c
c         IF (ipart.GT.idim) print *,iproc,': Fx ',ipart
#ifdef _OPENMP
         CALL OMP_SET_LOCK(forcei_lock(ipart/10+1))
#endif
         fxyzu(1,ipart) = fxyzu(1,ipart) + forcexi
         fxyzu(2,ipart) = fxyzu(2,ipart) + forceyi
         fxyzu(3,ipart) = fxyzu(3,ipart) + forcezi

         poten(ipart) = poten(ipart) + poteni 

         IF (gasdrag) THEN
            planetesimaltimestep(ipart) =
     &           planetesimaltimestep(ipart) + ptsi
            planetesimalnorm(ipart) =
     &           planetesimalnorm(ipart) + ptsnormi
         ENDIF

         IF (iphase(ipart).EQ.0) THEN

c            IF (ipart.EQ.0418 .AND. itime.GE.134217728) THEN
c               print *,iproc,': P0418dqs ',dq(ipart),dqi
c            ENDIF

            dq(ipart) = dq(ipart) + dqi

#ifdef DGRID
            IF (gasdrag) THEN
C$OMP CRITICAL
               dragqchange(idragr, idragp, idragt) = 
     &              dragqchange(idragr, idragp, idragt) + dq_dragi*
     &              dt*isteps(ipart)/imaxstep
C$OMP END CRITICAL
            ENDIF
#endif

            IF (imhd.EQ.idim) THEN
               dBxyz(1,ipart) = dBxyz(1,ipart) + dBxyz_xi
               dBxyz(2,ipart) = dBxyz(2,ipart) + dBxyz_yi
               dBxyz(3,ipart) = dBxyz(3,ipart) + dBxyz_zi
               IF (imhdevol.EQ.4) THEN
                  dBxyz(1,ipart) = dBxyz(1,ipart) + dBpsix
                  dBxyz(2,ipart) = dBxyz(2,ipart) + dBpsiy
                  dBxyz(3,ipart) = dBxyz(3,ipart) + dBpsiz
               ENDIF
            ENDIF
         ENDIF

#ifdef _OPENMP
         CALL OMP_UNSET_LOCK(forcei_lock(ipart/10+1))
#endif
c         IF (ipart.GT.idim) print *,iproc,': Fx2 ',ipart
      END DO
C$OMP END DO

#ifdef MPICOPY
c
c--For MPICOPY job, unfortunately need to do MPI reduce because in the above
c     loop particle 'i' gives forces, dq, poten, and dBxyz to particle 'j'.
c
c--Transfer fxyzu(1-3,i) - do not need to transfer all 4 elements of this
c     arrary - only need to transfer 3 because change in 'u' is done
c     only using neighbours of 'i', not 'j'.
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst_end
         DO k = 1, 3
            realtransfer5to15((i-1)*3+k) = 
     &           fxyzu(k,llisttrans(i)+1)
            realtransfer5to15(3*idim+(i-1)*3+k) = 0.
         END DO
      END DO
C$OMP END DO
C$OMP SINGLE
c
c--All processes reduce the fxyzu(1-3) data
c
      CALL MPI_ALLREDUCE(realtransfer5to15,
     &     realtransfer5to15(3*idim+1),3*nlst_end,
     &     MPI_REAL8,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Now add reduced fxyz_reduce(1-3) on to existing fxyzu(1-3) data in correct 
c     places in fxyzu(1-4)
c
      j = iproc
C$OMP DO SCHEDULE(runtime)
      DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &     irecvcounti(j+1)
         DO k = 1, 3
            fxyzu(k,llisttrans(i) + 1) = 
     &           realtransfer5to15(3*idim+(i-1)*3 + k)
         END DO
      END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': fxyzu(1-3) reduce ',numbertodoherekeep,
     &     irecvcounti(j+1), nlst_end
#endif

c
c--Reduce poten(i)
c
C$OMP END SINGLE
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst_end
         real4transfer2(i) = poten(llisttrans(i)+1)
         real4transfer2(idim+i) = 0.
      END DO
C$OMP END DO
C$OMP SINGLE
c
c--All processes reduce the poten() data
c
      CALL MPI_ALLREDUCE(real4transfer2,real4transfer2(idim+1),
     &     nlst_end,MPI_REAL4,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
c
c--Now place reduced poten() data in correct places in poten()
c
C$OMP END SINGLE
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst_end
         poten(llisttrans(i) + 1) = real4transfer2(idim+i)
      END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': poten() reduce ',numbertodoherekeep,
     &     irecvcounti(j+1), nlst_end
#endif

c
c--Reduce dq(i)
c
C$OMP END SINGLE
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst_end
         real4transfer2(i) = dq(llisttrans(i)+1)
         real4transfer2(idim+i) = 0.
      END DO
C$OMP END DO
C$OMP SINGLE
c
c--All processes reduce the poten() data
c
      CALL MPI_ALLREDUCE(real4transfer2,real4transfer2(idim+1),
     &     nlst_end,MPI_REAL4,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Now transfer reduced real4_reduce() into correct places in dq().
c     Only need to keep the values of those particles being done on this
c     MPICOPY process.
c
      j = iproc
C$OMP DO SCHEDULE(runtime)
      DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &     irecvcounti(j+1)
         dq(llisttrans(i) + 1) = real4transfer2(idim+i)
      END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': dq() reduce ',numbertodoherekeep,
     &     irecvcounti(j+1), nlst_end
#endif
C$OMP END SINGLE

#ifdef DGRID
      IF (gasdrag) THEN
c--Sum dragqchange(i,j,k) across processes.
c     
         itot = idragresr*idragresp*idragrest
C$OMP DO SCHEDULE(runtime)
         DO i = 1, idragresr
            DO j = 1, idragresp
               DO k = 1, idragrest
                  ii = (i-1)*(idragresp*idragrest) + (j-1)*idragrest + k
                  realtransfer5to15(ii) = dragqchange(i,j,k)
                  realtransfer5to15(itot + ii) = 0.
               ENDDO
            ENDDO
         ENDDO
C$OMP END DO
C$OMP SINGLE
c
c--All processes reduce the dragqchange data
c
      CALL MPI_ALLREDUCE(realtransfer5to15,
     &     realtransfer5to15(itot+1),itot,
     &     MPI_REAL8,MPI_SUM,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Now return dragqchange to normal array in correct places
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, idragresr
         DO j = 1, idragresp
            DO k = 1, idragrest
               ii = (i-1)*(idragresp*idragrest) + (j-1)*idragrest + k
               dragenergy(i,j,k) = dragenergy(i,j,k) +
     &              realtransfer5to15(itot+ii)
            ENDDO
         ENDDO
      ENDDO
C$OMP END DO
C$OMP SINGLE

#ifdef MPIDEBUG
      print *,iproc,': dragenergy summed in force ',numbertodoherekeep,
     &     irecvcounti(j+1), nlst_end
#endif
C$OMP END SINGLE
      ENDIF
#endif

      IF (imhd.EQ.idim) THEN
c
c--Transfer dBxyz(1-3,i)
c
C$OMP DO SCHEDULE(runtime)
         DO i = 1, nlst_end
            DO k = 1, 3
               realtransfer5to15((i-1)*3+k) = 
     &              dBxyz(k,llisttrans(i)+1)
               realtransfer5to15(3*idim+(i-1)*3+k) = 0.
            END DO
         END DO
C$OMP END DO
C$OMP SINGLE
c
c--All processes reduce the dBxyz(1-3) data
c
         CALL MPI_ALLREDUCE(realtransfer5to15,
     &        realtransfer5to15(3*idim+1),3*nlst_end,
     &        MPI_REAL8,MPI_SUM,
     &        MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Now do reduction of dBxyzu(1-3) and put data in correct 
c     places in dBxyz(1-3)
c
C$OMP DO SCHEDULE(runtime)
         DO i = 1, nlst_end
            DO k = 1, 3
               dBxyz(k,llisttrans(i) + 1) = 
     &              realtransfer5to15(3*idim+(i-1)*3 + k)
            END DO
         END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
         print *,iproc,': dBxyz(1-3) reduce ',numbertodoherekeep,
     &        irecvcounti(j+1), nlst_end
#endif
C$OMP END SINGLE

      ENDIF
#else
#ifdef DGRID
c--Sum dragenergy in non-MPI case.
      IF (gasdrag) THEN
C$OMP DO SCHEDULE(runtime)
         DO i = 1, idragresr
            DO j = 1, idragresp
               DO k = 1, idragrest
                  dragenergy(i,j,k) = dragenergy(i,j,k) +
     &                 dragqchange(i,j,k)
               ENDDO
            ENDDO
         ENDDO
C$OMP END DO
      ENDIF
#endif
#endif

C$OMP DO SCHEDULE(runtime)
c
c--For MPICOPY, all processes do this locally because it is fast (no loop
c     over neighbours) and probably not worth the transfer time to copy
c     sections between MPI processes
c
      DO n = nlst_in, nlst_tot
         ipart = listp(n)
         IF (iphase(ipart).EQ.0) THEN
c
c--Energy conservation
c
         IF (idim_h2.EQ.1) THEN

           CALL energ(ipart,npart,ntot,realtime,vxyzu,fxyzu,xyzmh,trho)

         ELSE

           IF (it1(ipart).EQ.imax) THEN
             dtime = (dt*isteps(ipart)/2)/imaxstep
           ELSE
             dtime = (dt*isteps(ipart))/imaxstep
           ENDIF 

           CALL energ_cooling(ipart,npart,ntot,realtime,vxyzu,fxyzu,
     &           xyzmh,h2ratio,dtime,abhpq,abeq,abHIq,abco,trho)

         END IF 
c
c--Translate to dB/dt and add divv term if B is evolved
c
            IF (varmhd.EQ.'Bvol') THEN
               rhoi = trho(ipart)
               rho1i = 1./rhoi

               DO k = 1, 3
                  dBxyz(k,ipart) = rhoi*dBxyz(k,ipart)
               END DO

               IF (n.LE.nlst_end) THEN
                  IF (ipart.GT.npart) THEN
                     WRITE (*,*) 'ERROR - ipart.GT.npart'
                     CALL quit(1)
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
     &                 - (Bxyz(1,ipart)+Beqxi)*rho1i*divv(ipart)
                  dBxyz(2,ipart) = dBxyz(2,ipart)
     &                 - (Bxyz(2,ipart)+Beqyi)*rho1i*divv(ipart)
                  dBxyz(3,ipart) = dBxyz(3,ipart)
     &                 - (Bxyz(3,ipart)+Beqzi)*rho1i*divv(ipart)

               ENDIF
            ENDIF
         ENDIF
      END DO
C$OMP END DO
c
c--Set flag for whether a particle has been updated or not back to false
c
C$OMP DO SCHEDULE(runtime)
      DO ipart = 1, nlistupdated
         iupdated(listparents(ipart)) = .FALSE.
      END DO
C$OMP END DO 
C$OMP END PARALLEL

#ifdef MPICOPY
c
c--For MPICOPY job, need to transfer fxyzu(1-4,i), dha(1-1+isizealphaMM,i), 
c     poten(), dBxyz(1-3,i), divcurlB(1-4,i)
c
c     NOTE: dq() does not need to be transferred since it is only used for
c     each ipart in call to energ.f in this subroutine.
c
c     NOTE: Code assumes that dgrav() is NOT BEING USED (i.e. nlmax.EQ.1)
c
c--Transfer fxyzu(1-4,i)
c
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)

      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     i4REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their fxyzu(1-4) data
c
      CALL MPI_ALLGATHERV(fxyzu,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,i4REAL8,
     &     MPI_COMM_WORLD,ierr)
c
c--Put fxyzu(1-4) data into correct places in fxyzu(1-4)
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(idisplacementsi,irecvcounti,fxyzu,llisttrans)
C$OMP& shared(realtransfer5to15,numproc,iproc)
C$OMP& private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 4
                  fxyzu(k,llisttrans(i) + 1) =
     &                 realtransfer5to15((i-1)*4 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': fxyzu(1-4) move ',numbertodoherekeep
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
c
c--Transfer dha(1-1+isizealphaMM,i)
c
      CALL MPI_TYPE_CONTIGUOUS(1 + isizealphaMM, MPI_REAL4, iaREAL4, 
     &     ierr)
      CALL MPI_TYPE_COMMIT(iaREAL4,ierr)

      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     iaREAL4, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their dha(1-1+isizealphaMM) data
c
      CALL MPI_ALLGATHERV(dha,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,iaREAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put dha(1-1+isizealphaMM) data into correct places in dha()
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(idisplacementsi,irecvcounti,dha,llisttrans)
C$OMP& shared(real4transfer2,numproc,iproc)
C$OMP& private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 1+isizealphaMM
c                  dha(k,llisttrans(i) + 1) =
c     &                 real4transfer2((i-1)*(1+isizealphaMM) + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': dha(1-1+isizealphaMM) move ',numbertodoherekeep
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
      CALL MPI_TYPE_FREE(iaREAL4,ierr)

c
c--Unify nearest gas neighbour list for solid particles
c
      IF (gasdrag)
     &     CALL MPI_ALLGATHERV (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
     &     nearestgas, irecvcounti, idisplacementsi, MPI_INTEGER,
     &     MPI_COMM_WORLD, ierr)


      IF (imhd.EQ.idim) THEN
c
c--Transfer dBxyz(1-3,i)
c
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)

      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     i3REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their dBxyz(1-3) data
c
      CALL MPI_ALLGATHERV(dBxyz,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,i3REAL8,
     &     MPI_COMM_WORLD,ierr)
c
c--Put dBxyz(1-3) data into correct places in dBxyz
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(idisplacementsi,irecvcounti,dBxyz,llisttrans)
C$OMP& shared(realtransfer5to15,numproc,iproc)
C$OMP& private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 3
                  dBxyz(k,llisttrans(i) + 1) =
     &                 realtransfer5to15((i-1)*3 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': dBxyz(1-3) move ',numbertodoherekeep
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
c
c--Transfer divcurlB(1-4,i)
c
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL4, i4REAL4, ierr)
      CALL MPI_TYPE_COMMIT(i4REAL4,ierr)

      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     i4REAL4, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their divcurlB(1-4) data
c
      CALL MPI_ALLGATHERV(divcurlB,1,indexMPI_INT1,
     &     real4transfer2,irecvcounti,idisplacementsi,i4REAL4,
     &     MPI_COMM_WORLD,ierr)
c
c--Put divcurlB(1-4) data into correct places in divcurlB
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(idisplacementsi,irecvcounti,divcurlB,llisttrans)
C$OMP& shared(real4transfer2,numproc,iproc)
C$OMP& private(i,j,k)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 4
                  divcurlB(k,llisttrans(i) + 1) =
     &                 real4transfer2((i-1)*4 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPIDEBUG
      print *,iproc,': divcurlB(1-4) move ',numbertodoherekeep
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
      CALL MPI_TYPE_FREE(i4REAL4,ierr)

      ENDIF

#endif

 4567  CONTINUE
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine forcei')

#ifdef MPIALL
#ifdef MPIDEBUG
      print *,iproc,': Exited forcei'
#endif
#endif

      RETURN
      END


      SUBROUTINE extrapolate(j,dt,itime,trho,vxyzu,pr,vsound,
     &                       ekcle,Bxyz,Bevolxyz)
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
      INCLUDE 'COMMONS/phase'
      
      DIMENSION vxyzu(4,idim2)
      REAL*4 trho(idim2), pr(idim2), vsound(idim2)
      DIMENSION ekcle(5,iradtrans2)
      DIMENSION Bxyz(3,imhd2)
      DIMENSION Bevolxyz(imhdevol,imhd3)
      
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
      
      IF (iphase(j).EQ.0) THEN

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

         CALL eospg(j, vxyzu, trho, pr, vsound, ekcle)
      ENDIF

      RETURN

      END
