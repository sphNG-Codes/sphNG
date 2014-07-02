      SUBROUTINE step (dt)
c************************************************************
c                                                           *
c  This subroutine integrate the system of differential     *
c     equations for one timestep using                      *
c     a second-order Leap-Frog method.                      *
c                                                           *
c************************************************************
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/debpt'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/accurpt'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/binfile'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/hagain'
      INCLUDE 'COMMONS/crpart'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/out1'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/sync'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/radtran2'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/xforce'

c      DIMENSION nsteplist(30)

      LOGICAL ifirst
      LOGICAL icenter
      LOGICAL*1 irevise,iwillupdatetree,imaketree
      CHARACTER*7 where
      DATA ifirst/.true./
      DATA icenter/.false./
      DATA small/1.0E-04/

      DATA where/'step'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE(iprint,250)
  250 FORMAT(' entry subroutine step')
c
c--Set integrator to L-F
c
      integrator = 1
      IF (XSPH) THEN
         WRITE (iprint,*) 'ERROR - Leap-frog and XSPH incompatible'
         WRITE (*,*) 'ERROR - Leap-frog and XSPH incompatible'
         CALL quit(1)
      ENDIF
c
c--Compute next dump time and initialise variables
c
      ikilled = 0
      itnext = imaxstep
      iteighth = itnext/8
      iptout = imaxstep/iptoutnum
      xlog2 = 0.30103
      ncrit = INT(nactive/10.0)
      ilocal = 0
      naddedplanet = 0
      itreeupdate = .FALSE.
#ifdef MPI
      print *,iproc,': ENTERED STEP ',nptmasstot
#endif
c
c--Integration
c
c---------- FIRST TIME AROUND ----------
c
      IF (ifirst) THEN
         ifirst = .FALSE.

#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
         IF (integrator.EQ.1) THEN
            WRITE (iprint,*) 'Integrator: Leap-Frog'
            WRITE (iprint,*)
         ELSE
            WRITE (iprint,*) 'Integrator miss-match'
            CALL quit(1)
         ENDIF
#ifdef MPI
      ENDIF
#endif

         IF (gt.EQ.0.0) THEN
c
c--Set timesteps to zero before first call to derivi
c
            DO i = 1, npart
               isteps(i) = 0
            END DO
         ELSE
            istepmin = imax
            istepmax = -1
            DO i = 1, npart
               IF (iphase(i).NE.-1) THEN
                  istepmin = MIN(istepmin, isteps(i))
                  istepmax = MAX(istepmax, isteps(i))
               ENDIF
            END DO
         ENDIF
c
c--Set timestep of each bin and zero the number of particles in each bin
c
         nlst = 0
         DO i = 1, nbinmax
            nlstbins(i) = 0
            it2bin(i) = 2**i
         END DO
c
c--Set up particle timesteps and make dummy variables
c
         DO i = 1, npart
            IF (iphase(i).NE.-1) THEN
c
c--Set current integer time and integer time of next force evaluation
c
               it0(i) = 0
               it2(i) = isteps(i)
c
c--Set dummy variables at force evaluation time (need dummys in general
c     because need to predict properties of particles that are on different
c     timesteps
c
               DO k = 1, 5
                  dumxyzmh(k,i) = xyzmh(k,i)
               END DO
               DO k = 1, 4
                  dumvxyzu(k,i) = vxyzu(k,i)
               END DO
               IF (imhd.EQ.idim) THEN
                  DO k = 1, imhdevol
                     dumBevolxyz(k,i) = Bevolxyz(k,i)
                  END DO
               ENDIF

               IF (ifsvi.GE.6) THEN
                  DO k=1,isizealphaMM
                     dumalpha(k,i) = alphaMM(k,i)
                  ENDDO
               ENDIF

cccc               IF (ran1(1).LT.0.001) THEN
               nlst = nlst + 1
               llist(nlst) = i
cccc               ENDIF
               iscurrent(i) = .TRUE.
            ENDIF
         END DO
         nlst0 = nlst
c
c--Create ghost particles
c
         CALL ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)

         DO i = npart + 1, ntot
            DO k = 1, 5
               dumxyzmh(k,i) = xyzmh(k,i)
            END DO
            DO k = 1, 4
               dumvxyzu(k,i) = vxyzu(k,i)
            END DO
            IF (ifsvi.GE.6) THEN
               DO k=1,isizealphaMM
                  dumalpha(k,i) = alphaMM(k,ireal(i))
               ENDDO
            ENDIF
            IF (imhd.EQ.idim) THEN
               DO k = 1, imhdevol
                  dumBevolxyz(k,i) = Bevolxyz(k,i)
               END DO
            ENDIF
         END DO
c
c--Set integer time to zero
c
         itime = 0
c
c--Make the tree
c
         IF (igrape.EQ.0 .AND. (nactive.NE.nptmass .OR. 
     &                                      iptintree.GT.0)) THEN
            CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)
         ENDIF
c
c--Compute forces on all particles
c
         icall = 1

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , 1',
c     &        itime0,itime1
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , 1',
c     &        itime0,itime1
c            WRITE (*,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            WRITE (iprint,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            CALL quit(1)
c         ENDIF
c      END DO

         IF (itiming) CALL getused(ts101)

         CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,
     &        npart,ntot,ireal,dumalpha,ekcle,dumBevolxyz,f1Bxyz)
       
         IF (itiming) THEN
            CALL getused(ts102)
            ts10 = ts10 + (ts102 - ts101)
         ENDIF

         IF (gt.EQ.0.0) THEN
c
c--Set new timesteps - needs vsound, also needs default timestep to alter
c
            nlst = 0
            DO i = 1, npart
               IF (iphase(i).NE.-1) THEN
                  nlst = nlst + 1
                  llist(nlst) = i
                  isteps(i) = imaxstep
                  IF (iphase(i).EQ.0) THEN
                     CALL eospg(i,vxyzu,rho,pr,vsound,ekcle)
c
c--if rho not defined it is a problem
c
c                     vsound(i) = SQRT(2./3.*vxyzu(4,i))
                  ENDIF
               ENDIF
            END DO
            nlst0 = nlst

            CALL timestep(dt,itnext,nlst,llist,f1vxyzu,Bevolxyz)

            IF (individualtimesteps.EQ.0) THEN
c
c--All particles on smallest timestep OR NOT
c
               DO i = 1, npart
                  isteps(i) = istepmin
               END DO
            ELSEIF (individualtimesteps.EQ.1) THEN
c
c--Sinks ALL on smallest timestep OR NOT
c
               DO i = 1, nptmass
                  iptcur = listpm(i)
                  isteps(iptcur) = istepmin
               END DO
c
c--Additional constraint for planet with surface for migration
c
               IF (ibound.EQ.102 .AND. initialptm.EQ.5) THEN
                  DO j = 1, nptmass
                     IF (iphase(listpm(j)).EQ.5) THEN
                        radiushill2 = (xyzmh(4,listpm(j))/(3.*xmass))
     &                       **(2./3.)
                        DO i = 1, npart
                           IF ((iphase(i).LT.1 .OR. iphase(i).GT.10)
     &                          .AND. it0(i).EQ.itime) THEN
                              dx = xyzmh(1,i) - xyzmh(1,listpm(j))
                              dy = xyzmh(2,i) - xyzmh(2,listpm(j))
                              dz = xyzmh(3,i) - xyzmh(3,listpm(j))
                              radius = dx*dx + dy*dy + dz*dz
                              IF (radius.LT.uniformtslim2*radiushill2)
     &                             THEN
                                 isteps(i) = istepmin
                                 it2(i) = it0(i) + isteps(i)
                              ENDIF
                           ENDIF
                        END DO
                     ENDIF
                  ENDDO
               ENDIF

            ELSEIF (individualtimesteps.EQ.2) THEN
               DO i = 1, nptmass
                  iptcur = listpm(i)                  
                  IF (istepmingas/2.LT.isteps(iptcur)) THEN
                     isteps(iptcur) = MAX(istepmin,istepmingas/2)
                  ENDIF
               END DO
            ELSE
               WRITE (iprint,*) 'ERROR - individualtimesteps ',
     &              individualtimesteps
               CALL quit(1)
            ENDIF
         ENDIF
c
c--Set up particle timesteps and make dummy variables
c
         DO i = 1, npart
            IF (iphase(i).NE.-1) THEN
c
c--Set current integer time and integer time of next force evaluation
c
               it0(i) = 0
               it2(i) = isteps(i)
c
c--Set lists of all particles in a particular bin
c
               ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
               IF (ibin.GT.nbinmax) THEN
                  WRITE (*,*) 'ERROR - ibin.GT.nbinmax'
                  WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax'
                  CALL quit(1)
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               listbins(nlstbins(ibin),ibin) = i
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin'
                  WRITE (iprint,*) 'ERROR - it2bin'
                  CALL quit(1)
               ENDIF
            ENDIF
         END DO

c         write (*,*) f1vxyzu(1,1),f1vxyzu(1,2)

c         print *,'dha(1,1): ',dumxyzmh(1,1),dumvxyzu(1,1),
c     &        f1vxyzu(1,1),f1vxyzu(1,100),f1ha(1,1),f1ha(2,1)
c         DO i = idim-10, idim
c            print *,'f1ha(1,1): ',i,f1vxyzu(1,i),f1vxyzu(2,i),
c     &           f1vxyzu(3,i),f1vxyzu(4,i)
c         END DO
c         DO i = 1, 10
c            print *,'f1ha(1,1): ',f1ha(1,i),f1ha(2,i)
c         END DO
        
c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , 1b',
c     &        itime0,itime1
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , 1b',
c     &        itime0,itime1
c            WRITE (*,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            WRITE (iprint,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            CALL quit(1)
c         ENDIF
c      END DO

         rhomaxsync = 0.
         DO i = 1, npart
            iscurrent(i) = .FALSE.
            rhomaxsync = MAX(rhomaxsync, rho(i))
         END DO

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , 1c',
c     &        itime0,itime1
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , 1c',
c     &        itime0,itime1
c            WRITE (*,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            WRITE (iprint,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            CALL quit(1)
c         ENDIF
c      END DO
c
c--DIVERGENCE CLEANING
c
c         IF (imhd.EQ.idim) THEN
c            IF (varmhd.EQ.'Brho' .OR. varmhd.EQ.'Bvol') THEN
c               divBmax = 0.
c               DO i=1,npart
c                  divBmax = max(divcurlB(1,i),divBmax)
c               ENDDO
c               WRITE(iprint,*) 'div B max = ',divBmax
c               CALL divBclean(npart,ntot,xyzmh,rho,Bevolxyz)
c
c               DO i = 1, npart
c                  IF (iphase(i).NE.-1) THEN
c                     DO k = 1, 3
c                        dumBevolxyz(k,i) = Bevolxyz(k,i)
c                     END DO
c                  ENDIF
c               ENDDO
c               CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,
c     &              npart,ntot,ireal,dumalpha,ekcle,dumBevolxyz,f1Bxyz)
c            ENDIF
c         ENDIF
c         RETURN

      ELSE
         ntot = npart + nghost

c
c--DIVERGENCE CLEANING
c
c         IF (imhd.EQ.idim) THEN
c            IF (varmhd.EQ.'Brho' .OR. varmhd.EQ.'Bvol') THEN
c               divBmax = 0.
c               DO i=1,npart
c                  divBmax = max(divcurlB(1,i),divBmax)
c               ENDDO
c               WRITE(iprint,*) 'div B max = ',divBmax
c               CALL divBclean(npart,ntot,xyzmh,rho,Bevolxyz)

c               DO i = 1, npart
c                  IF (iphase(i).NE.-1) THEN
c                     DO k = 1, 3
c                        dumBevolxyz(k,i) = Bevolxyz(k,i)
c                     END DO
c                  ENDIF
c               ENDDO
c            ENDIF
c         ENDIF
         
      END IF
c
c---------- END FIRST TIME AROUND ----------
c
c--Find minimum time for force calculation
c
 100  itimelocal = imax
      IF (itiming) CALL getused(ts1p1)
      DO i = 1, nbinmax
         IF (nlstbins(i).GT.0) itimelocal = MIN(itimelocal, it2bin(i))
      END DO
#ifdef MPI
      CALL MPI_ALLREDUCE(itimelocal,itime,1,MPI_INTEGER,MPI_MIN,
     &     MPI_COMM_WORLD,ierr)
#ifdef MPIDEBUG
c      IF (gt.GT.63.7) 
          print *,iproc,': TIMES ',itimelocal, itime
#endif
#else
      itime = itimelocal
#endif
      IF (itiming) THEN
         CALL getused(ts1p2)
         ts1 = ts1 + (ts1p2 - ts1p1)
      ENDIF

c      IF (itime1.NE.itime1new .OR. itime0.NE.itime0new) THEN
c         WRITE (*,*) 'ERROR - times'
c         WRITE (iprint,*) 'ERROR - times'
c         WRITE (*,*) 'Choose:',itime1,itime0,itime1new,itime0new
c         WRITE (iprint,*) 'Choose:',itime1,itime0,itime1new,itime0new
c
c      DO j = 1, 30
c         nsteplist(j) = 0
c      END DO
c
c      DO j = 1, npart
c         IF (iphase(j).GE.0) THEN
c            ibin = INT(LOG10(REAL(isteps(j)))/xlog2+0.5)
c            IF ((ibin.GE.1) .AND. (ibin.LE.29)) THEN
c               nsteplist(ibin) = nsteplist(ibin) + 1
c            ELSE
c               nsteplist(30) = nsteplist(30) + 1
c            ENDIF
c         ENDIF
c      END DO
c         
c        DO i = 1, nbinmax 
c         WRITE (*,*) 2**i,it2bin(i),nlstbins(i),nsteplist(i)
c         WRITE (iprint,*) 2**i,it2bin(i),nlstbins(i),
c     &        nsteplist(i)
c        END DO
c        CALL quit(1)
c      ENDIF

      ioutinf = 0
      ioutsup = 0
      ioutmin = 0
      ioutmax = 0
      inmax = -1
      inmin = 10000

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , 2',
c     &        itime0,itime1
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , 2',
c     &        itime0,itime1
c            WRITE (*,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            WRITE (iprint,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            CALL quit(1)
c         ENDIF
c      END DO

c
c--Allow for tracing flow
c
      IF(itrace.EQ.'all') WRITE(iprint,251) dt*itime/imaxstep+gt, itime
 251  FORMAT(' step is doing time ',1F12.5, I10)
c
c---------- ADVANCING PARTICLES ----------
c
c--Identify particles to be advanced
c
      jj = 0
      imakeghost = 0
c
      IF (itiming) CALL getused(ts31)
      itbinupdate = 0
c
c--Loop over all timestep bins
c
      DO i = 1, nbinmax
c
c--For those that need to be updated at this time...
c
         IF (it2bin(i).EQ.itime) THEN
            itbinupdate = i
c
c--Loop over all particles in that bin, check validity of timestep info
c
            DO j = 1, nlstbins(i)
               ipart = listbins(j,i)
               IF (ipart.GT.idim) THEN
                  WRITE (*,*) 'ERROR - ipart.GT.idim 2'
                  WRITE (iprint,*) 'ERROR - ipart.GT.idim 2'
                  CALL quit(1)
               ELSEIF (ipart.LE.0) THEN
                  WRITE (*,*) 'ERROR - ipart.LE.0 2'
                  WRITE (iprint,*) 'ERROR - ipart.LE.0 2'
                  CALL quit(1)
               ENDIF
               IF (iphase(ipart).GE.0) THEN
                  IF (it2(ipart).NE.itime) THEN
                     WRITE (*,*) 'ERROR - it2(ipart).NE.itime 2'
                     WRITE (iprint,*) 'ERROR - it2(ipart).NE.itime 2'
                     WRITE (*,*) ipart,it2(ipart),itime,it2bin(i)
                     WRITE (*,*) it2(ipart),it0(ipart)
                     CALL quit(1)
                  ENDIF
                  IF (it0(ipart)+isteps(ipart).NE.itime) THEN
                     WRITE (*,*) 'ERROR - it0+istep.NE.itime 2'
                     WRITE (*,*) ipart,it2(ipart),itime,it2bin(i)
                     WRITE (*,*) isteps(ipart),it0(ipart),iphase(ipart)
                     WRITE (*,*) listpm(1),listpm(2),iproc
                     WRITE (iprint,*) 'ERROR - it0+istep.NE.itime 2'
                     CALL quit(1)
                  ENDIF
c
c--Make list
c
                  jj = jj + 1
                  IF (jj.GT.idim) THEN
                     WRITE (*,*) 'ERROR - jj.GT.idim 2'
                     WRITE (iprint,*) 'ERROR - jj.GT.idim 2'
                     CALL quit(1)
                  ENDIF
                  llist(jj) = ipart
                  IF (hasghost(ipart)) imakeghost = 1
               ENDIF
            END DO
         ENDIF
      END DO
      itbinupdate = MIN(itbinupdate + 7, nbinmax)

#ifdef MPIDEBUG
      print *,iproc,': itbinupdate ',itbinupdate
#endif

      nlst = jj
      nlstnneigh = nlst
      nlst0 = nlst

      IF (itiming) THEN
         CALL getused(ts32)
         ts3 = ts3 + (ts32 - ts31)
      ENDIF

      IF (itiming) CALL getused(ts71)
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst,llist,iscurrent,dt,isteps,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,rho)
C$OMP& shared(it0,itime,imaxdens,cnormk,iener)
C$OMP& shared(iprint,nneigh,hmaximum,iphase,nptmass,listpm)
C$OMP& shared(xmomsyn,ymomsyn,zmomsyn,pmass,listrealpm)
C$OMP& shared(ifsvi,alphaMM,alphamax,Bevolxyz,encal,gt)
C$OMP& shared(pdvvisc)
C$OMP& private(i,j,k,dtfull,dthalf,delvx,delvy,delvz)
C$OMP& private(iii,pmasspt)
C$OMP& reduction(+:ioutmax)
      DO j = 1, nlst
         i = llist(j)
         iscurrent(i) = .TRUE.
c
c--Synchronize the advanced particle times with current time
c
         it0(i) = itime

         dtfull = dt*isteps(i)/imaxstep
         dthalf = 0.5*dtfull
c
c--Update velocities (kick 1/2)
c
         delvx = dthalf*f1vxyzu(1,i)
         delvy = dthalf*f1vxyzu(2,i)
         delvz = dthalf*f1vxyzu(3,i)
c         IF (iphase(i).GT.0 .AND. nptmass.GE.2) 
c     &        print *,'f1 ',i,f1vxyzu(2,i),vxyzu(1,i),vxyzu(2,i)
         vxyzu(1,i) = vxyzu(1,i) + delvx
         vxyzu(2,i) = vxyzu(2,i) + delvy
         vxyzu(3,i) = vxyzu(3,i) + delvz
c
c--Update positions (drift full dt)
c
         DO k = 1, 3
            xyzmh(k,i) = xyzmh(k,i) + dtfull*vxyzu(k,i)
         END DO
c
c--Update momentum of sink particles
c
         IF (iphase(i).GE.1 .AND. iphase(i).LT.10) THEN
c            DO iii = 1, nptmass
c               IF (listpm(iii).EQ.i) GOTO 444
c            END DO
c 444        CONTINUE
            iii = listrealpm(i)
            pmasspt = xyzmh(4,i)
            xmomsyn(iii) = xmomsyn(iii) + delvx*pmasspt
            ymomsyn(iii) = ymomsyn(iii) + delvy*pmasspt
            zmomsyn(iii) = zmomsyn(iii) + delvz*pmasspt
         ELSE
c
c--Update u(i) and h(i)
c
            IF (encal.EQ.'r') THEN
               vxyzu(4,i) = dumvxyzu(4,i)
               pdvvisc(i) = dtfull*f1vxyzu(4,i)
            ELSEIF (encal.EQ.'m') THEN
               vxyzu(4,i) = dumvxyzu(4,i)
            ELSE
               vxyzu(4,i) = vxyzu(4,i) + dtfull*f1vxyzu(4,i)
            ENDIF

c correction of energy if using K & I cooling curve

            IF (iener.EQ.3) CALL coolcurve(vxyzu(4,i),rho(i),dtfull)

            IF (iener.EQ.2 .AND. vxyzu(4,i).LT.0.0) vxyzu(4,i)=0.15

            xyzmh(5,i) = xyzmh(5,i) + dtfull*f1ha(1,i)
            IF (xyzmh(5,i).LT.0) WRITE(iprint,*) 'error in h ', 
     &           xyzmh(5,i), f1ha(1,i), dtfull, nneigh(i)
            IF (hmaximum.GT.0.0) THEN
               IF (xyzmh(5,i).GT.hmaximum) THEN
                  xyzmh(5,i) = hmaximum
                  ioutmax = ioutmax + 1
               ENDIF
            ENDIF
c
c--Update viscosity switch
c
            IF (ifsvi.GE.6) THEN
               DO k=1,isizealphaMM
                  alphaMM(k,i) = MIN(alphamax(k), alphaMM(k,i) +
     &               dtfull*f1ha(k+1,i))
               ENDDO
            ENDIF
c
c--Update MHD
c
            IF (imhd.EQ.idim) THEN
               DO k = 1, imhdevol
                  Bevolxyz(k,i) = Bevolxyz(k,i) + dtfull*f1Bxyz(k,i)
               END DO
            ENDIF
         ENDIF
      END DO
C$OMP END PARALLEL DO
      IF (itiming) THEN
         CALL getused(ts72)
         ts7 = ts7 + (ts72 - ts71)
      ENDIF

c      PRINT *,'h(1)2: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)

c
c--Total check for particles that have moved outside boundary
c
      idonebound = 0
      IF (ibound.GT.0) THEN
         IF (ibound.LT.7 .OR. ibound.EQ.11) 
     &     CALL boundry(npart, llist, nlst, xyzmh, vxyzu, idonebound) 
c
c--Create ghost particles
c
         IF ((nlst.GT.ncrit .AND. imakeghost.EQ.1) 
     &                                    .OR. idonebound.EQ.1) THEN
c         IF (imakeghost.EQ.1 .OR. idonebound.EQ.1) THEN

         CALL ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)

         ENDIF
      ENDIF
c
c--Get ready to make or update tree
c
      iwillupdatetree = igrape.EQ.0 .AND. (nlst.GT.nptmass .OR. 
     &     nlstacc.GT.0 .OR. iptintree.GT.0 .OR. icreate.EQ.1)
      imaketree = iwillupdatetree .AND. (nlst.GT.ncrit .OR. 
     &     imakeghost.EQ.1 .OR. idonebound.EQ.1 .OR. icreate.EQ.1)
      irevise = iwillupdatetree .AND. (.NOT. imaketree) .AND.
     &     (.NOT.(iptintree.EQ.1 .AND. nlst.LE.nptmass .AND.
     &           nlstacc.EQ.0))
c      print *,iproc,': TREE ',iwillupdatetree,imaketree,irevise,
c     &     nlst,nptmass,nlstacc,iptintree
c
c--Predict variables at t=time
c
      IF (itiming) CALL getused(ts81)

      numberparents = 0

      IF (itbinupdate.GE.nbinmax-1 .OR. ipartialrevtree.EQ.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,dt,itime,it0,imaxstep,ireal)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,rho)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,dumBevolxyz)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,Bevolxyz)
C$OMP& shared(isibdaupar,iflagtree,imfac,numberparents,listparents)
C$OMP& private(j,k,deltat,deltat2,iparent)
      DO j = 1, npart
         IF (iphase(j).GE.0) THEN
c
c--Load predicted quantities into dummy arrays
c
            deltat = dt*(itime - it0(j))/imaxstep
            deltat2 = 0.5*deltat
            DO k = 1, 3
               dumvxyzu(k,j) = vxyzu(k,j) + deltat2*f1vxyzu(k,j)
            END DO
            DO k = 1, 3
               dumxyzmh(k,j) = xyzmh(k,j) + deltat*dumvxyzu(k,j)
            END DO
            dumxyzmh(4,j) = xyzmh(4,j)
            IF (iphase(j).EQ.0) THEN
               IF (encal.EQ.'r') THEN
                  dumvxyzu(4,j) = vxyzu(4,j)
               ELSEIF (encal.EQ.'m') THEN
                  dumvxyzu(4,j) = vxyzu(4,j)
               ELSE
                  dumvxyzu(4,j) = vxyzu(4,j) + deltat*f1vxyzu(4,j)
               ENDIF
              IF(iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
               dumxyzmh(5,j) = xyzmh(5,j) + deltat*f1ha(1,j)
               IF (ifsvi.GE.6) THEN
                  DO k=1,isizealphaMM
                    dumalpha(k,j) = MIN(alphamax(k),alphaMM(k,j)+
     &                 deltat*f1ha(1+k,j))
                  ENDDO
               ENDIF

           IF (iener.EQ.3) CALL coolcurve(dumvxyzu(4,j),rho(j),deltat)

               IF (imhd.EQ.idim) THEN
                  DO k = 1, imhdevol
                     dumBevolxyz(k,j) = Bevolxyz(k,j)+deltat*f1Bxyz(k,j)
                  END DO
               ENDIF
            ELSEIF (iphase(j).GE.1 .AND. iphase(j).LT.10) THEN
               dumxyzmh(5,j) = xyzmh(5,j)
            ENDIF
         ENDIF
      END DO
C$OMP END PARALLEL DO

      ELSEIF (nlst-nptmass.GT.10000 .OR. ipartialrevtree.EQ.1) THEN
         DO i = 1, itbinupdate
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(i,nlstbins,listbins,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,rho)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,dumBevolxyz)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,Bevolxyz)
C$OMP& shared(isibdaupar,iflagtree,imfac,npart)
C$OMP& shared(numberparents,listparents)
C$OMP& private(j,k,ipart,deltat,deltat2,iparent)
            DO j = 1, nlstbins(i)
               ipart = listbins(j,i)
               IF (iphase(ipart).GE.0) THEN
c
c--Load predicted quantities into dummy arrays
c
                  deltat = dt*(itime - it0(ipart))/imaxstep
                  deltat2 = 0.5*deltat
                  DO k = 1, 3
                     dumvxyzu(k,ipart) = vxyzu(k,ipart) + 
     &                    deltat2*f1vxyzu(k,ipart)
                  END DO
                  DO k = 1, 3
                     dumxyzmh(k,ipart) = xyzmh(k,ipart) + 
     &                    deltat*dumvxyzu(k,ipart)
                  END DO
                  dumxyzmh(4,ipart) = xyzmh(4,ipart)
                  IF (iphase(ipart).EQ.0) THEN
                     IF (encal.EQ.'r') THEN
                        dumvxyzu(4,ipart) = vxyzu(4,ipart)
                     ELSEIF (encal.EQ.'m') THEN
                        dumvxyzu(4,ipart) = vxyzu(4,ipart)
                     ELSE
                        dumvxyzu(4,ipart) = vxyzu(4,ipart) +
     &                       deltat*f1vxyzu(4,ipart)
                     ENDIF
                     IF (iener.EQ.2 .AND. dumvxyzu(4,ipart).LT.0.0) 
     &                    dumvxyzu(4,ipart)=0.15
                     dumxyzmh(5,ipart) = xyzmh(5,ipart) + 
     &                    deltat*f1ha(1,ipart)
                     IF (ifsvi.GE.6) THEN
                        DO k=1,isizealphaMM
                           dumalpha(k,ipart) = MIN(alphamax(k),
     &                        alphaMM(k,ipart) + deltat*f1ha(1+k,ipart))
                        ENDDO
                     ENDIF

           IF (iener.EQ.3) CALL coolcurve(dumvxyzu(4,ipart),rho(ipart),
     &      deltat)

                     IF (imhd.EQ.idim) THEN
                        DO k = 1, imhdevol
                           dumBevolxyz(k,ipart) = 
     &                         Bevolxyz(k,ipart)+deltat*f1Bxyz(k,ipart)
                        END DO
                     ENDIF
                  ELSEIF (iphase(ipart).GE.1 .AND. 
     &                    iphase(ipart).LT.10) THEN
                     dumxyzmh(5,ipart) = xyzmh(5,ipart)
                  ENDIF
               ENDIF
            END DO
C$OMP END PARALLEL DO
         END DO
c
c--Sinks ALL on smallest timestep OR NOT
c
c      ELSEIF (nlst.GT.nptmass) THEN
      ELSEIF (nlst.GT.0) THEN
         DO i = 1, itbinupdate
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(i,nlstbins,listbins,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,rho)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,dumBevolxyz)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,Bevolxyz)
C$OMP& shared(isibdaupar,iflagtree,imfac,npart)
C$OMP& shared(numberparents,listparents,irevise,iprint)
C$OMP& private(j,k,ipart,deltat,deltat2,iparent)
            DO j = 1, nlstbins(i)
               ipart = listbins(j,i)

               IF (iphase(ipart).GE.0) THEN
c
c--Set flag to state that parent node in tree needs to be recalculated
c                  
                  IF (irevise .AND.
     &                 (iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10 .OR.
     &                 iptintree.GE.1)) THEN
                     iparent = isibdaupar(3,ipart)
C$OMP CRITICAL (parentlist1)
                     IF (.NOT.iflagtree(iparent)) THEN
                        iflagtree(iparent) = .TRUE.
                        numberparents = numberparents + 1
                        IF (numberparents.GT.idim) THEN
                         WRITE (iprint,*) 'numberparents ',numberparents
                           CALL quit(1)
                        ENDIF
                        listparents(numberparents) = iparent
                     ENDIF
C$OMP END CRITICAL (parentlist1)
                     IF (ipart.GT.npart .OR.
     &                  (iphase(ipart).GE.1 .AND. 
     &                  iphase(ipart).LT.10 .AND. iptintree.EQ.1)) THEN
                        imfac(ipart) = 0
                     ELSE
                        imfac(ipart) = 1
                     ENDIF
                  ENDIF
c
c--Load predicted quantities into dummy arrays
c
                  deltat = dt*(itime - it0(ipart))/imaxstep
                  deltat2 = 0.5*deltat
                  DO k = 1, 3
                     dumvxyzu(k,ipart) = vxyzu(k,ipart) + 
     &                    deltat2*f1vxyzu(k,ipart)
                  END DO
                  DO k = 1, 3
                     dumxyzmh(k,ipart) = xyzmh(k,ipart) + 
     &                    deltat*dumvxyzu(k,ipart)
                  END DO
                  dumxyzmh(4,ipart) = xyzmh(4,ipart)
                  IF (iphase(ipart).EQ.0) THEN
                     IF (encal.EQ.'r') THEN
                        dumvxyzu(4,ipart) = vxyzu(4,ipart)
                     ELSEIF (encal.EQ.'m') THEN
                        dumvxyzu(4,ipart) = vxyzu(4,ipart)
                     ELSE
                        dumvxyzu(4,ipart) = vxyzu(4,ipart) +
     &                       deltat*f1vxyzu(4,ipart)
                     ENDIF
                     IF (iener.EQ.2 .AND. dumvxyzu(4,ipart).LT.0.0) 
     &                    dumvxyzu(4,ipart)=0.15
                     dumxyzmh(5,ipart) = xyzmh(5,ipart) + 
     &                    deltat*f1ha(1,ipart)
                     IF (ifsvi.GE.6) THEN
                        DO k=1,isizealphaMM
                           dumalpha(k,ipart) = MIN(alphamax(k),
     &                    alphaMM(k,ipart) + deltat*f1ha(k+1,ipart))
                        ENDDO
                     ENDIF


           IF (iener.EQ.3) CALL coolcurve(dumvxyzu(4,ipart),rho(ipart),
     &          deltat)

                     IF (imhd.EQ.idim) THEN
                        DO k = 1, imhdevol
                           dumBevolxyz(k,ipart) = 
     &                        Bevolxyz(k,ipart)+deltat*f1Bxyz(k,ipart)
                        END DO
                     ENDIF
                  ELSEIF (iphase(ipart).GE.1 .AND. 
     &                                 iphase(ipart).LT.10) THEN
                     dumxyzmh(5,ipart) = xyzmh(5,ipart)
                  ENDIF
               ENDIF
            END DO
C$OMP END PARALLEL DO
         END DO
      ELSE
c
c--Only update ptmasses, since nothing else is evolved
c     These should already be at the correct time, so don't need to extrapolate
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nptmass,listpm,itime,it0)
C$OMP& shared(xyzmh,vxyzu,dumxyzmh,dumvxyzu)
C$OMP& shared(isibdaupar,iflagtree,imfac)
C$OMP& shared(numberparents,listparents,irevise)
C$OMP& private(i,k,ipart,iparent)
         DO i = 1, nptmass
            ipart = listpm(i)
            IF (itime.NE.it0(ipart)) THEN
               WRITE (*,*) 'ERROR - ptmass not updated'
               CALL quit(1)
            ENDIF
c
c--Set flag to state that parent node in tree needs to be recalculated
c     but only if point masses in the tree and all gravity done by tree
c
            IF (irevise .AND. iptintree.EQ.2) THEN
               iparent = isibdaupar(3,ipart)
C$OMP CRITICAL (parentlist2)
               IF (.NOT.iflagtree(iparent)) THEN
                  iflagtree(iparent) = .TRUE.
                  numberparents = numberparents + 1
                  listparents(numberparents) = iparent
               ENDIF
C$OMP END CRITICAL (parentlist2)
               imfac(ipart) = 1
            ELSE
               imfac(ipart) = 0
            ENDIF
c
c--Load dummy arrays
c
            DO k = 1, 5
               dumxyzmh(k,ipart) = xyzmh(k,ipart)
            END DO
            DO k = 1, 3
               dumvxyzu(k,ipart) = vxyzu(k,ipart)
            END DO
         END DO
C$OMP END PARALLEL DO
      ENDIF

      IF (nghost.GT.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,ireal,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,dumxyzmh,dumvxyzu,encal)
C$OMP& shared(iphase,iener,ifsvi,dumalpha,alphaMM,alphamax)
C$OMP& shared(Bevolxyz,dumBevolxyz,rho)
C$OMP& shared(isibdaupar,iflagtree,imfac,irevise)
C$OMP& shared(numberparents,listparents)
C$OMP& private(j,k,l,deltat,iparent)
         DO j = npart + 1, npart + nghost
            IF (iphase(j).GE.0) THEN
c
c--Set flag to state that parent node in tree needs to be recalculated
c
               IF (irevise) THEN
                  iparent = isibdaupar(3,j)
C$OMP CRITICAL (parentlist3)
                  IF (.NOT.iflagtree(iparent)) THEN
                     iflagtree(iparent) = .TRUE.
                     numberparents = numberparents + 1
                     listparents(numberparents) = iparent
                  ENDIF
C$OMP END CRITICAL (parentlist3)
                  imfac(j) = 0
               ENDIF
c
c--Load predicted quantities into dummy arrays
c
               k = ireal(j)
               deltat = dt*(itime - it0(k))/imaxstep
               DO l = 1, 3
                  dumvxyzu(l,j) = vxyzu(l,j)
               END DO
               IF (encal.NE.'r' .AND. encal.NE.'m') THEN
                  dumvxyzu(4,j) = vxyzu(4,j) + deltat*f1vxyzu(4,k)
               ELSE
                  dumvxyzu(4,j) = vxyzu(4,j)
               ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
               DO l = 1, 3
                  dumxyzmh(l,j) = xyzmh(l,j) + deltat*dumvxyzu(l,j)
               END DO
               dumxyzmh(4,j) = xyzmh(4,j)
               dumxyzmh(5,j) = xyzmh(5,j) + deltat*f1ha(1,k)
               IF (ifsvi.GE.6) THEN
                  DO l=1,isizealphaMM
                     dumalpha(l,j) = MIN(alphamax(l),alphaMM(l,k)
     &              + deltat*f1ha(1+l,k))
                  ENDDO
               ENDIF

           IF (iener.EQ.3) CALL coolcurve(dumvxyzu(4,j),rho(j),deltat)

               IF (imhd.EQ.idim) THEN
                  DO l = 1, imhdevol
                     dumBevolxyz(l,j) = Bevolxyz(l,j)+deltat*f1Bxyz(l,k)
                  END DO
               ENDIF
            ENDIF
         END DO
C$OMP END PARALLEL DO
      ENDIF
      IF (itiming) THEN
         CALL getused(ts82)
         ts8 = ts8 + (ts82 - ts81)
      ENDIF

c      PRINT *,'h(1)3: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)

      IF (hmaximum.GT.0.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,hmaximum,dumxyzmh)
C$OMP& private(j)
         DO j = 1, npart + nghost
            IF (dumxyzmh(5,j).GT.hmaximum) dumxyzmh(5,j) = hmaximum
         END DO
C$OMP END PARALLEL DO
      ENDIF
c
c--Make or update the tree
c
      IF (iwillupdatetree) THEN
         IF (imaketree) THEN
            CALL insulate(1, ntot, npart, dumxyzmh, f1vxyzu)
         ELSEIF (irevise) THEN
            CALL insulate(2, ntot, npart, dumxyzmh, f1vxyzu)
         ELSEIF (numberparents.NE.0) THEN
            WRITE (*,*) 'ERROR: Setting iflagtree to ZERO'
            WRITE (iprint,*) 'ERROR: Setting iflagtree to ZERO'
            DO i = 1, numberparents
               iflagtree(listparents(i)) = .FALSE.
            END DO
         ENDIF
         iaccr = 0
         icreate = 0
         ikilled = 0
         nlstacc = 0
      ENDIF
c
c--Compute forces on list particles
c
 200  icall = 3

      IF (itiming) CALL getused(ts101)
      CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,npart,
     &     ntot,ireal,dumalpha,ekcle,dumBevolxyz,f1Bxyz)

c      IF (nptmasstot.GE.1) THEN
c      DO j = 1, nlst
c         i = llist(j)
c         IF (dumxyzmh(1,i)**2 + dumxyzmh(2,i)**2 + dumxyzmh(3,i)**2
c     &        .LT.2.0e-4) print *,iproc,': ',i,j,iunique(iorig(i)),
c     &     SQRT(dumxyzmh(1,i)**2 + dumxyzmh(2,i)**2 + dumxyzmh(3,i)**2),
c     &     f1vxyzu(1,i),poten(i)
c      END DO
c      ENDIF

      IF (itiming) THEN
         CALL getused(ts102)
         ts10 = ts10 + (ts102 - ts101)
      ENDIF

c      PRINT *,'h(1)4: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)

c
c--Update velocities (kick 1/2)
c
      IF (itiming) CALL getused(ts111)
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst,llist,dt,isteps,imaxstep)
C$OMP& shared(vxyzu,f1vxyzu,itime)
C$OMP& shared(iphase,nptmass,listpm,listrealpm)
C$OMP& shared(xmomsyn,ymomsyn,zmomsyn,xyzmh,gt)
C$OMP& private(i,j,dthalf,delvx,delvy,delvz,iii,pmasspt)
      DO j = 1, nlst
         i = llist(j)
         dthalf = 0.5*dt*isteps(i)/imaxstep
c
c--Update velocities (kick 1/2)
c
         delvx = dthalf*f1vxyzu(1,i)
         delvy = dthalf*f1vxyzu(2,i)
         delvz = dthalf*f1vxyzu(3,i)

c         IF (iphase(i).GT.0 .AND. nptmass.GE.2) 
c    &        print *,'f2 ',i,f1vxyzu(2,i),vxyzu(1,i),vxyzu(2,i)

         vxyzu(1,i) = vxyzu(1,i) + delvx
         vxyzu(2,i) = vxyzu(2,i) + delvy
         vxyzu(3,i) = vxyzu(3,i) + delvz

         IF (iphase(i).GE.1 .AND. iphase(i).LT.10) THEN
c            DO iii = 1, nptmass
c               IF (listpm(iii).EQ.i) GOTO 445
c            END DO
c 445        CONTINUE
            iii = listrealpm(i)
            pmasspt = xyzmh(4,i)
            xmomsyn(iii) = xmomsyn(iii) + delvx*pmasspt
            ymomsyn(iii) = ymomsyn(iii) + delvy*pmasspt
            zmomsyn(iii) = zmomsyn(iii) + delvz*pmasspt
         ENDIF
      END DO
C$OMP END PARALLEL DO
      IF (itiming) THEN
         CALL getused(ts112)
         ts11 = ts11 + (ts112 - ts111)
      ENDIF
c
c--Synchronization time
c
c      PRINT *,'h(1)5: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)
      idtsyn = itnext - itime
      IF (idtsyn.EQ.0) idtsyn = imaxstep
      ikilled = 0
      istepmin = imax
      istepmax = 0
      iptnum = 0
      time = dt*itime/imaxstep + gt
c
c--Set new timesteps
c
      IF (itiming) CALL getused(ts121)
      CALL timestep(dt,idtsyn,nlst,llist,f1vxyzu,Bevolxyz)
      IF (itiming) THEN
         CALL getused(ts122)
         ts12 = ts12 + (ts122 - ts121)
      ENDIF
c
c--Make point mass timesteps equal to the minimum time step used
c
c      PRINT *,'h(1)6: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)
      IF (itiming) CALL getused(ts131)

      IF (individualtimesteps.EQ.0) THEN
c
c--All particles on smallest timestep OR NOT
c
         DO i = 1, npart
            isteps(i) = istepmin
         END DO
      ELSEIF (individualtimesteps.EQ.1) THEN
c
c--Sinks ALL on smallest timestep OR NOT
c
         DO i = 1, nptmass
            iptcur = listpm(i)
            isteps(iptcur) = istepmin
         END DO
c
c--Additional constraint for planet with surface for migration
c
         IF (ibound.EQ.102 .AND. initialptm.EQ.5) THEN
            DO j = 1, nptmass
               IF (iphase(listpm(j)).EQ.5) THEN
                  radiushill2 = (xyzmh(4,listpm(j))/(3.*xmass))**(2./3.)
                  DO i = 1, npart
                     IF ((iphase(i).LT.1 .OR. iphase(i).GT.10)
     &                    .AND. it0(i).EQ.itime) THEN
                        dx = xyzmh(1,i) - xyzmh(1,listpm(j))
                        dy = xyzmh(2,i) - xyzmh(2,listpm(j))
                        dz = xyzmh(3,i) - xyzmh(3,listpm(j))
                        radius = dx*dx + dy*dy + dz*dz
                        IF (radius.LT.uniformtslim2*radiushill2) THEN
                           isteps(i) = istepmin
                           it2(i) = it0(i) + isteps(i)
                        ENDIF
                     ENDIF
                  END DO
               ENDIF
            ENDDO
         ENDIF

      ELSEIF (individualtimesteps.EQ.2) THEN
         DO i = 1, nptmass
            iptcur = listpm(i)
            IF (istepmingas/2.LT.isteps(iptcur)) THEN
               isteps(iptcur) = MAX(istepmin,istepmingas/2)
            ENDIF
         END DO
      ELSE
         WRITE (iprint,*) 'ERROR - individualtimesteps',
     &        individualtimesteps
         CALL quit(1)
      ENDIF

      IF (itiming) THEN
         CALL getused(ts132)
         ts13 = ts13 + (ts132 - ts131)
      ENDIF
c      PRINT *,'h(1)7: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)

c
c--Tidy up - kill particles, set u(i)=dumu(i), set new it2()
c
      IF (itiming) CALL getused(ts91)
C$OMP PARALLEL default(none)
C$OMP& shared(nlst,llist)
C$OMP& shared(it2,isteps,igphi)
C$OMP& shared(dgrav,dumvxyzu,vxyzu)
C$OMP& shared(it0,ibound,deadbound)
C$OMP& shared(iprint,xyzmh,dumxyzmh,poten,iphase,ikillpr)
C$OMP& shared(ikilled,nactive,nkill,time,iorig,nlstacc,listacc)
C$OMP& shared(anglostx,anglosty,anglostz)
C$OMP& shared(dumBevolxyz,Bevolxyz)
C$OMP& shared(variation,phibound,nlistinactive,listinactive)
C$OMP& shared(rcyl,rmind)
C$OMP& shared(pradfac,rplanet,px,py,irotpot,iplanetesimals)
C$OMP& private(i,j,k,ipkill,rlim)
C$OMP& private(r2,pmassi,phiparticle)
c
c--Compute new time steps for particles which have been advanced
c
C$OMP DO SCHEDULE(runtime)
      DO 855 j = 1, nlst
         i = llist(j)
c
c--Dead particle boundaries
c
         IF (ibound.EQ.8 .OR. ibound/10.EQ.9) THEN
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            IF (r2.GT.(1.05*deadbound)**2) THEN
               ipkill = iphase(i)
               iphase(i) = -1
               pmassi = xyzmh(4,i)
C$OMP CRITICAL (killparticle)
               ikilled = 1

               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR step nlstacc'
                  CALL quit(1)
               ENDIF
               listacc(nlstacc) = i

               nactive = nactive - 1
               nkill = nkill + 1
               WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              poten(i),dgrav(i),ipkill
               CALL FLUSH (ikillpr)
               anglostx = anglostx + pmassi*(xyzmh(2,i)*vxyzu(3,i) - 
     &              vxyzu(2,i)*xyzmh(3,i))
               anglosty = anglosty + pmassi*(vxyzu(1,i)*xyzmh(3,i) - 
     &              xyzmh(1,i)*vxyzu(3,i))
               anglostz = anglostz + pmassi*(xyzmh(1,i)*vxyzu(2,i) - 
     &              vxyzu(1,i)*xyzmh(2,i))
C$OMP END CRITICAL (killparticle)
               GOTO 855
            ENDIF
         ENDIF
c
c--Dead particle boundaries - planet embedded in disc
c
         IF (ibound.EQ.100) THEN
            phiparticle = ATAN2(xyzmh(2,i),xyzmh(1,i))
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            IF (r2.LT.(1.0-variation)**2 .OR. r2.GT.(1.0+variation)**2
     &                         .OR. ABS(phiparticle).GT.phibound) THEN
               ipkill = iphase(i)
               iphase(i) = -1
C$OMP CRITICAL (killparticle)
               ikilled = 1

               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR step nlstacc'
                  CALL quit(1)
               ENDIF
               listacc(nlstacc) = i

               nlistinactive = nlistinactive + 1
               IF (nlistinactive.GT.idim) THEN
                  WRITE (iprint,*) 'ERROR step nlistinactive'
                  CALL quit(1)
               ENDIF
               listinactive(nlistinactive) = i

               nactive = nactive - 1
               nkill = nkill + 1
               WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              poten(i),dgrav(i),ipkill
               CALL FLUSH(ikillpr)
C$OMP END CRITICAL (killparticle)
               GOTO 855
            ENDIF
         ENDIF

c
c--Dead particle boundaries - planet embedded in whole disc
c
         IF (ibound.EQ.102) THEN
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            IF ((r2.GT.rcyl**2 .AND. iphase(i).EQ.0) .OR.
     &           (r2.LE.rmind**2 .AND. iphase(i).GT.0) .OR.
     &           (r2.GT.(2.*rcyl)**2 .AND. iphase(i).GT.0)) THEN
               ipkill = iphase(i)
               iphase(i) = -1
C$OMP CRITICAL (killparticle)
               ikilled = 1

               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR step nlstacc'
                  CALL quit(1)
               ENDIF
               listacc(nlstacc) = i

               nactive = nactive - 1
               nkill = nkill + 1
               WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              poten(i),dgrav(i),ipkill
               CALL FLUSH(ikillpr)
C$OMP END CRITICAL (killparticle)
               GOTO 855
            ENDIF
         ENDIF
c
c--Kill planetesimals approaching planet modelled as potential.
c
         IF (irotpot.EQ.1 .AND. iplanetesimals.GT.0 .AND.
     &        iphase(i).EQ.11) THEN
            r2 = (xyzmh(1,i)-px)**2 + (xyzmh(2,i)-py)**2 +
     &           xyzmh(3,i)**2
            rlim = rplanet*pradfac(1)
            IF (r2.LT.rlim**2) THEN
               ipkill = iphase(i)
               iphase(i) = -1
#ifndef MPICOPY
C$OMP CRITICAL (killparticle)
#endif
               ikilled = 1

               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR step nlstacc'
                  CALL quit(1)
               ENDIF
               listacc(nlstacc) = i

               nactive = nactive - 1
               nkill = nkill + 1
               WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              poten(i),dgrav(i),ipkill
               CALL FLUSH(ikillpr)
#ifndef MPICOPY
C$OMP END CRITICAL (killparticle)
#endif
               GOTO 855
            ENDIF
         ENDIF
c
c--Set u to it's new value from DERIVI 
c     For polytropic equation of state, the du's are not used - u(i) 
c       is calculated directly from the density, in each derivi call and
c       put into dumu(i). Hence must be transferred from dumu(i) to u(i).
c     For adiabatic (or isothermal) u(i) is calculated via the du's
c       but this setting of u(i)=dumu(i) doesn't matter as the u(i)
c       updated at the full timestep above, then put into dumu(i)
c       but the derivi call doesn't alter them, so putting them back
c       into u(i) again changes nothing.
c
         vxyzu(4,i) = dumvxyzu(4,i)
c
c--also set Bevol to (possibly changed) value from derivi
c  (can be changed by divergence cleaning, B smoothing)
c
         IF (imhd.EQ.idim) THEN
            DO k=1,imhdevol
               Bevolxyz(k,i) = dumBevolxyz(k,i)
            ENDDO
         ENDIF
c
c--nlmax=1 is the sign that code is running using 'grad-h' so that h is set
c     inside derivi rather than being evolved.
c
         IF (nlmax.EQ.1) xyzmh(5,i) = dumxyzmh(5,i)
c
c--Set new timestep values
c
         it2(i) = it0(i) + isteps(i)
 855  CONTINUE
C$OMP END DO
C$OMP END PARALLEL
      IF (itiming) THEN
         CALL getused(ts92)
         ts9 = ts9 + (ts92 - ts91)
      ENDIF

c      PRINT *,'h(1)8: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)

c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint,252) dt*itime/imaxstep, 
     &     nlst, itime
 252  FORMAT (' step time ', 1F12.5,', particles moved ', I6,
     &     ' int time ', I10)
c
c--Update time
c
      IF (itime.GE.iteighth) THEN
#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
         WRITE(iprint,253) dt*itime/imaxstep + gt,
     &        dt*itime/imaxstep, nlst
 253     FORMAT ('Dynamic time = ', 1PE17.10,
     &        ' Step time = ', 1PE12.5,' Moved ', I10)

         CALL FLUSH (iprint)
         CALL FLUSH (iptprint)
#ifdef MPI
      ENDIF
#endif
         iteighth = itime + imaxstep/8
      ENDIF
c
c--Accrete particles near point mass, or create point mass
c
c      IF (nlst.GT.nptmass .AND. (nptmass.NE.0 .OR. icreatetot.EQ.1)) THEN
      IF (nptmasstot.NE.0 .OR. icreatetot.EQ.1) THEN
         isave = 0
         IF (itime.GE.iptout) THEN
            isave = 1
            iptout = itime + imaxstep/iptoutnum
         ENDIF
         realtime = dt*itime/imaxstep + gt
c
c--Accrete particles near point mass, or create point mass
c 
         IF (itiming) CALL getused(taccrete1)

         CALL accrete(dt, realtime, isave, ntot)

         IF (itiming) THEN
            CALL getused(taccrete2)
            taccrete = taccrete + (taccrete2 - taccrete1)
         ENDIF
      ENDIF
c
c--Update it2 for bins
c
      IF (itiming) CALL getused(ts141)
      DO i = 1, nbinmax
         IF (it2bin(i).EQ.itime) THEN
            nlstbins(i) = 0
            it2bin(i) = itime + 2**i
         ELSEIF (it2bin(i).LT.itime) THEN
            it2bin(i) = itime + 2**i
            IF (nlstbins(i).NE.0) THEN
               WRITE (*,*) 'ERROR - nlstbins(i).NE.0'
               WRITE (iprint,*) 'ERROR - nlstbins(i).NE.0'
               CALL quit(1)
            ENDIF
         ENDIF
      END DO
c
c--Update the lists of particles in each bin
c
      DO j = 1, nlst
         i = llist(j)
         iscurrent(i) = .FALSE.
         IF (i.GT.idim) THEN
            WRITE (*,*) 'ERROR - i.GT.idim 4'
            WRITE (iprint,*) 'ERROR - i.GT.idim 4'
            CALL quit(1)
         ELSEIF (i.LE.0) THEN
            WRITE (*,*) 'ERROR - i.LE.0 4'
            WRITE (iprint,*) 'ERROR - i.LE.0 4'
            WRITE (*,*) nlst,j,itime
            WRITE (iprint,*) nlst,j,itime
            DO i2 = 1, j
               WRITE (iprint,*) i2, llist(i2)
            END DO
            CALL quit(1)
         ENDIF
         IF (iphase(i).NE.-1) THEN
            ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
            IF (ibin.GT.nbinmax) THEN
               WRITE (*,*) 'ERROR - ibin.GT.nbinmax 4'
               WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax 4'
               CALL quit(1)
            ENDIF
            nlstbins(ibin) = nlstbins(ibin) + 1
            IF (nlstbins(ibin).GT.idim) THEN
               WRITE (*,*) 'ERROR - nlstbins(ibin).GT.idim 4'
               WRITE (iprint,*) 'ERROR - nlstbins(ibin).GT.idim 4'
               CALL quit(1)
            ENDIF
            listbins(nlstbins(ibin),ibin) = i
            IF (it2bin(ibin).NE.it2(i)) THEN
               WRITE (*,*) 'ERROR - it2bin 2',it2bin(ibin),it2(i),i,
     &              it0(i),itime,isteps(i)
               CALL quit(1)
            ENDIF
         ENDIF
      END DO
      IF (itiming) THEN
         CALL getused(ts142)
         ts14 = ts14 + (ts142 - ts141)
      ENDIF
c
c--Create NEW particles to keep number of particles within a shell
c     constant.  
c
      IF (nptmasstot.NE.0 .OR. icreatetot.EQ.1) THEN
         IF (ibound/10.EQ.9 .AND. nshell.GT.inshell) THEN
            IF (nlmax.EQ.1) THEN
               WRITE (*,*) 'Cannot use grad-h with shell boundaries'
               CALL quit(1)
            ENDIF
            iaccr = 0
            ikilled = 0
            nlstacc = 0
            nneightotsave = nneightot

            nlst = nshell-inshell
            nlst0 = nlst
            WRITE (iprint,*) 'Add ',nlst
            DO i = 1, nlst
               nnew = npart + i
               isort(nnew) = nnew
               iorig(nnew) = nnew
               llist(i) = nnew
               CALL phoenix(nnew, idtsyn, itime)
               nactive = nactive + 1
            END DO
            npart = npart + nlst
            n1 = n1 + nlst
            IF (npart.GT.idim) THEN
               CALL error(where,3)
            ENDIF

            DO j = 1, nlst
               iscurrent(llist(j)) = .TRUE.
            END DO

            nghost = 0
            CALL ghostp3(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
            ntot = npart + nghost

            DO j = npart - nlst + 1, npart
               IF (iphase(j).GE.0) THEN
                  DO k = 1, 5
                     dumxyzmh(k,j) = xyzmh(k,j)
                  END DO
                  DO k = 1, 4
                     dumvxyzu(k,j) = vxyzu(k,j)
                  END DO
                  IF (ifsvi.GE.6) THEN
                     DO k=1,isizealphaMM
                        dumalpha(k,j) = alphaMM(k,j)
                     ENDDO
                  ENDIF
               ENDIF
            END DO
            DO j = npart + 1, npart + nghost
               IF (iphase(j).GE.0) THEN
                  k = ireal(j)
                  deltat = dt*(itime - it0(k))/imaxstep
                  DO l = 1, 3
                     dumxyzmh(l,j) = xyzmh(l,j) + deltat*vxyzu(l,j)
                  END DO
                  dumxyzmh(4,j) = xyzmh(4,j)
                  dumxyzmh(5,j) = xyzmh(5,j) + deltat*f1ha(1,k)

                  DO l = 1, 3
                     dumvxyzu(l,j) = vxyzu(l,j)
                  END DO
                  IF (encal.NE.'r' .AND. encal.NE.'m') THEN
                     dumvxyzu(4,j) = vxyzu(4,j) + deltat*f1vxyzu(4,k)
                  ELSE
                     dumvxyzu(4,j) = vxyzu(4,j)
                  ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
                  IF (ifsvi.GE.6) THEN
                     DO l=1,isizealphaMM
                       dumalpha(l,j) = MIN(alphamax(l),
     &                               alphaMM(l,k) + deltat*f1ha(1+l,k))
                     ENDDO
                   ENDIF
               ENDIF
            END DO

            IF (iener.EQ.3) CALL coolcurve(dumvxyzu(4,j),rho(j),
     &      deltat)

            IF (igrape.EQ.0) THEN
               CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)
            ENDIF

            neighmean = (neimax + neimin)/2

            iokay = 1
            DO j = 1, nlst
               i = llist(j)

               ivalue = 0
               ichkloop = 0
 2000          ichkloop = ichkloop + 1

               IF (igrape.EQ.0) THEN
                  CALL insulate(3,ntot,npart,dumxyzmh,f1vxyzu)
                  numneigh = nneigh(i)
               ELSE
                  nlist = 0
                  CALL getneigh(i,ntot,xyzmh(5,i),dumxyzmh,mmax,nlist,
     &                 iptneigh,nearl)
                  numneigh = nlist
               ENDIF

               IF (numneigh.GT.1) THEN
                  xyzmh(5,i) = (xyzmh(5,i)/
     &                 (FLOAT(numneigh)/FLOAT(neighmean))**(1.0/3.0) +
     &                 ivalue*xyzmh(5,i))/(ivalue + 1)
               ELSE
                  xyzmh(5,i) = xyzmh(5,i)*2.0
               ENDIF
               dumxyzmh(5,i) = xyzmh(5,i)
               IF (ichkloop.GT.10) THEN
                  ivalue = 2
                  IF (ichkloop.GT.20) ivalue = 4
                  IF (ichkloop.GT.30) ivalue = 8
               ENDIF
               IF (ichkloop.GT.500) CALL error(where,2)

               IF (numneigh.GT.(neimax - nrange) .OR. 
     &              numneigh.LT.(neimin + nrange)) GOTO 2000
            END DO

            icall = 4
            CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,
     &           npart,ntot,ireal,dumalpha,ekcle,dumBevolxyz,f1Bxyz)

            time = dt*itime/imaxstep + gt
            DO j = 1, nlst
               i = llist(j)
               WRITE (ireasspr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),poten(i)
               CALL FLUSH (ireasspr)
               iscurrent(i) = .FALSE.
            END DO

            nneightot = nneightotsave

c
c--Update the lists of particles in each bin
c
            DO j = 1, nlst0
               i = llist(j)
               ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
               IF (ibin.GT.nbinmax) THEN
                  WRITE (*,*) 'ERROR - ibin.GT.nbinmax 5'
                  WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax 5'
                  CALL quit(1)
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               IF (nlstbins(ibin).GT.idim) THEN
                  WRITE (*,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  WRITE (iprint,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  CALL quit(1)
               ENDIF
               listbins(nlstbins(ibin),ibin) = i

c               IF (it1bin(ibin).NE.it1(i)) THEN
c                  WRITE (*,*) 'ERROR - it1bin 5',it1bin(ibin),it1(i),i,
c     &                 it0(i),itime,isteps(i),it2(i)
c                  CALL quit(1)
c               ENDIF
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin 5',it2bin(ibin),it2(i),i,
     &                 it0(i),itime,isteps(i),it1(i)
                  CALL quit(1)
               ENDIF
c
c--Check that particle is not in any other bin
c
               DO k = 1, nbinmax
                  IF (k.NE.ibin) THEN
                     DO l = 1, nlstbins(k)
                        IF (listbins(l,k).EQ.i) THEN
                           WRITE (*,*) 'Particle ',i,' in two bins!'
                           WRITE (*,*) 'Bins ',ibin,k,nlstbins(k)
                           CALL quit(1)
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
c
c--End creation of NEW particles
c
         ENDIF

         IF (nactive - nptmass.LT.50 .AND. nactive.NE.nptmass) 
     &                                      CALL error(where,1)

      ENDIF
c
c--Create NEW particles in box surrounding planet in disc (ibound=100)
c
      IF (ibound.EQ.100) THEN
         IF (nlmax.EQ.1) THEN
            WRITE (*,*) 'Cannot use grad-h with ibound=100'
            CALL quit(1)
         ENDIF
c
c--Work out how many particles to add
c
         timelocal = dt*itime/imaxstep
         nlst0 = INT(timelocal*flowrate) - naddedplanet
         IF (nlst0.GT.10) THEN
            nlst = nlst0
            WRITE (iprint,*) 'Add ',nlst0,' recyc ',nlistinactive,
     &           ' npart ',npart,' nkill ',nkill
            WRITE (*,*) 'Add ',nlst0,' recyc ',nlistinactive,
     &           ' npart ',npart,' nkill ',nkill
c
c--Set values
c
            iaccr = 0
            ikilled = 0
            nlstacc = 0
            nneightotsave = nneightot
c
c--Create new particle in free spot, add to list llist
c
            icountaddnpart = 0
            DO i = 1, nlst0
               IF (i.LE.nlistinactive) THEN
                  nnew = listinactive(nlistinactive-i+1)
               ELSE
                  icountaddnpart = icountaddnpart + 1
                  nnew = npart + icountaddnpart
                  IF (nnew.GT.idim) CALL error(where,3)
               ENDIF
               llist(i) = nnew
               isort(nnew) = nnew
               iorig(nnew) = nnew
               CALL phoenix2(nnew, idtsyn, itime)
               nactive = nactive + 1
            END DO
            IF (icountaddnpart.NE.0 .AND.
     &           icountaddnpart.NE.nlst0-nlistinactive) THEN
               WRITE (*,*) 'ERROR: icountaddnpart'
               CALL quit(1)
            ENDIF
            npart = npart + MAX(0,nlst0-nlistinactive)
            n1 = n1 + MAX(0,nlst0-nlistinactive)
            nlistinactive = MAX(0,nlistinactive-nlst0)
            WRITE (*,*) 'New nlistinactive = ',nlistinactive
c
c--Set is current of new particles
c
            DO i = 1, nlst0
               iscurrent(llist(i)) = .TRUE.
            END DO
c
c--Redo ghosts
c
            nghost = 0
            CALL ghostp100(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
            ntot = npart + nghost
c
c--Set dummy's for new particles
c
            DO i = 1, nlst0
               j = llist(i)
               IF (iphase(j).GE.0) THEN
                  DO k = 1, 5
                     dumxyzmh(k,j) = xyzmh(k,j)
                  END DO
                  DO k = 1, 4
                     dumvxyzu(k,j) = vxyzu(k,j)
                  END DO
                  IF (ifsvi.GE.6) THEN
                     DO k=1,isizealphaMM
                        dumalpha(k,j) = alphaMM(k,j)
                     ENDDO
                  ENDIF
               ENDIF
            END DO
c
c--Set dummy's for all ghosts
c
            DO j = npart + 1, npart + nghost
               IF (iphase(j).GE.0) THEN
                  k = ireal(j)
                  deltat = dt*(itime - it0(k))/imaxstep
                  DO l = 1, 3
                     dumxyzmh(l,j) = xyzmh(l,j) + deltat*vxyzu(l,j)
                  END DO
                  dumxyzmh(4,j) = xyzmh(4,j)
                  dumxyzmh(5,j) = xyzmh(5,j) + deltat*f1ha(1,k)
                  DO l = 1, 3
                     dumvxyzu(l,j) = vxyzu(l,j)
                  END DO
                  IF (encal.NE.'r' .AND. encal.NE.'m') THEN
                     dumvxyzu(4,j) = vxyzu(4,j) + deltat*f1vxyzu(4,k)
                  ELSE
                     dumvxyzu(4,j) = vxyzu(4,j)
                  ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
                  IF (ifsvi.GE.6) THEN
                     DO l=1,isizealphaMM
                        dumalpha(l,j) = MIN(alphamax(l),
     &                     alphaMM(l,k) + deltat*f1ha(1+l,k))
                     ENDDO
                  ENDIF
               ENDIF
            END DO
c
c--Rebuild the tree
c
            IF (igrape.EQ.0) THEN
               CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)
            ENDIF
c
c--Set smoothing lengths of new particles
c
            neighmean = (neimax + neimin)/2

            iokay = 1
            DO j = 1, nlst0
               i = llist(j)

               ivalue = 0
               ichkloop = 0
 2001          ichkloop = ichkloop + 1

               IF (igrape.EQ.0) THEN
                  CALL insulate(3,ntot,npart,dumxyzmh,f1vxyzu)
                  numneigh = nneigh(i)
               ELSE
                  nlist = 0
                  CALL getneigh(i,ntot,xyzmh(5,i),dumxyzmh,mmax,nlist,
     &                 iptneigh,nearl)
                  numneigh = nlist
               ENDIF

               IF (numneigh.GT.1) THEN
                  xyzmh(5,i) = (xyzmh(5,i)/
     &                 (FLOAT(numneigh)/FLOAT(neighmean))**(1.0/3.0) +
     &                 ivalue*xyzmh(5,i))/(ivalue + 1)
               ELSE
                  xyzmh(5,i) = xyzmh(5,i)*2.0
               ENDIF
               dumxyzmh(5,i) = xyzmh(5,i)
               IF (ichkloop.GT.100) THEN
                  ivalue = 2
                  IF (ichkloop.GT.200) ivalue = 4
                  IF (ichkloop.GT.300) ivalue = 8
               ENDIF
               IF (ichkloop.GT.500)  THEN
                  IF (numneigh.LE.neimax.AND.
     &              numneigh.GE.neimin) GOTO 2002

                  WRITE (iprint,*) 'i,h,nneigh ',i,xyzmh(5,i),numneigh
                  WRITE (iprint,*) 'x,y,z ',xyzmh(1,i),xyzmh(2,i),
     &                 xyzmh(3,i)
                  CALL error(where,2)
               ENDIF

               IF (numneigh.GT.(neimax - nrange) .OR.
     &              numneigh.LT.(neimin + nrange)) GOTO 2001
 2002          CONTINUE
            END DO
c
c--Set accelerations on new particles using call to derivi
c
            icall = 4
            PRINT *,"icall 4 triggered"
            CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,
     &         f1ha,npart,ntot,ireal,dumalpha,ekcle,dumBevolxyz,f1Bxyz)
c
c--Write new particles to file
c
            time = dt*itime/imaxstep + gt
            DO j = 1, nlst0
               i = llist(j)
               WRITE (ireasspr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),poten(i)
               iscurrent(i) = .FALSE.
            END DO
            CALL FLUSH (ireasspr)
c
c--Restore values for neighbours output
c
            nneightot = nneightotsave
            timeflowold = time
c
c--Check that nactive and new particles set up correctly
c
            IF (nactive - nptmass.LT.50 .AND. nactive.NE.nptmass)
     &           CALL error(where,4)
c
c--Update the lists of particles in each bin
c
            DO j = 1, nlst0
               i = llist(j)
               ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
               IF (ibin.GT.nbinmax) THEN
                  WRITE (*,*) 'ERROR - ibin.GT.nbinmax 5'
                  WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax 5'
                  CALL quit(1)
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               IF (nlstbins(ibin).GT.idim) THEN
                  WRITE (*,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  WRITE (iprint,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  CALL quit(1)
               ENDIF
               listbins(nlstbins(ibin),ibin) = i

c               IF (it1bin(ibin).NE.it1(i)) THEN
c                  WRITE (*,*) 'ERROR - it1bin 5',it1bin(ibin),it1(i),i,
c     &                 it0(i),itime,isteps(i),it2(i)
c                  CALL quit(1)
c               ENDIF
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin 5',it2bin(ibin),it2(i),i,
     &                 it0(i),itime,isteps(i),it1(i)
                  CALL quit(1)
               ENDIF
c
c--Check that particle is not in any other bin
c
               DO k = 1, nbinmax
                  IF (k.NE.ibin) THEN
                     DO l = 1, nlstbins(k)
                        IF (listbins(l,k).EQ.i) THEN
                           WRITE (*,*) 'Particle ',i,' in two bins!'
                           WRITE (*,*) 'Bins ',ibin,k,nlstbins(k)
                           CALL quit(1)
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO

            naddedplanet = naddedplanet + nlst0
         ENDIF
c
c--End creation of NEW particles for Embedded Planet (ibound=100)
c
      ENDIF
c
c--If accreted mass and angular momentum is large enough, then add on to
c     previous point mass's mass and momentum
c
      IF (itiming) CALL getused(ts151)
      DO ii = 1, nptmass
         i = listpm(ii)
         IF (ptmadd(ii)/ptmsyn(ii).GT.0.001 .OR. itime.GE.itnext) THEN
            ptmsyn(ii) = ptmsyn(ii) + ptmadd(ii)
            ptmadd(ii) = 0.0
            pmasspt = ptmsyn(ii)
            vxyzu(1,i) = (xmomsyn(ii) + xmomadd(ii))/pmasspt
            vxyzu(2,i) = (ymomsyn(ii) + ymomadd(ii))/pmasspt
            vxyzu(3,i) = (zmomsyn(ii) + zmomadd(ii))/pmasspt
         ENDIF
      END DO
      IF (itiming) THEN
         CALL getused(ts152)
         ts15 = ts15 + (ts152 - ts151)
      ENDIF
c
c--If binary with mean-Roche-lobe sized accretion radii under massive
c      accretion, then dynamically evolve the accretion radii.
c      Also, CONSTRAIN sink particles to be on a circular orbit.
c
      IF (nptmass.EQ.2 .AND. (iaccevol.EQ.'v' .OR. 
     &                                    iaccevol.EQ.'s')) THEN
         ipt1 = listpm(1)
         ipt2 = listpm(2)
         qratio = xyzmh(4,ipt2)/xyzmh(4,ipt1)
         IF (qratio.GT.1.0) THEN
            ipt1 = listpm(2)
            ipt2 = listpm(1)
            qratio = xyzmh(4,ipt2)/xyzmh(4,ipt1)
         ENDIF
         qratio1 = qratio + 1.0

         totmass = xyzmh(4,ipt1) + xyzmh(4,ipt2)
         dx = xyzmh(1,ipt1) - xyzmh(1,ipt2)
         dy = xyzmh(2,ipt1) - xyzmh(2,ipt2)
         dz = xyzmh(3,ipt1) - xyzmh(3,ipt2)
         dvx = vxyzu(1,ipt1) - vxyzu(1,ipt2)
         dvy = vxyzu(2,ipt1) - vxyzu(2,ipt2)
         dvz = vxyzu(3,ipt1) - vxyzu(3,ipt2)
         dr = SQRT(dx*dx + dy*dy + dz*dz)
c
c--  Problem with this is that when modify the orbit to be circular the
c      angular momentum of the binary is not conserved.  Rather, want to
c      base the separation on a=L^2/(GM (1-e^2) but where e=0 so that
c      when we change the orbit to be circular it is energy that is
c      dumped, NOT angular momentum!
c
ccccc         dv2 = dvx*dvx+ dvy*dvy+ dvz*dvz
ccccc         binenergy = -totmass/dr + dv2/2.0
ccccc         semiaxis = - totmass/2.0/binenergy

         binj = dvy*dx - dvx*dy
         semiaxis = binj*binj/totmass
c
c--Set accretion radii to be the roche lobe sizes 
c     (see Accretion Power in Astrophysics, Frank, King, & Raine)
c
         IF (iaccevol.EQ.'v') THEN
            IF (qratio.GE.0.05) THEN
               xyzmh(5,ipt1) = accfac*semiaxis*(0.38-0.20*LOG10(qratio))
            ELSE
               WRITE (iprint,*) 'ERROR: qratio < 0.05'
               CALL quit(1)
            ENDIF
            IF (qratio.LT.0.5) THEN
               xyzmh(5,ipt2) = accfac*semiaxis*
     &              (0.462*(qratio/qratio1)**(1.0/3.0))
            ELSE
               xyzmh(5,ipt2) = accfac*semiaxis*(0.38+0.20*LOG10(qratio))
            ENDIF
c
c--Set accretion radii to be some fraction of the separation 
c
         ELSEIF (iaccevol.EQ.'s') THEN
            xyzmh(5,ipt1) = accfac*semiaxis
            xyzmh(5,ipt2) = accfac*semiaxis
         ENDIF

         hacc = xyzmh(5,ipt2)
         haccall = xyzmh(5,ipt2)
c
c--Modify orbit to be circular
c
c         IF (totmass.GE.1.05) THEN
c            r2 = semiaxis/qratio1
c            r1 = qratio*r2
c            v2 = SQRT(totmass/semiaxis)/qratio1
c            v1 = qratio*v2
c            dxdr = dx/dr
c            dydr = dy/dr
c            xyzmh(1,ipt1) = r1*dxdr
c            xyzmh(2,ipt1) = r1*dydr
c            xyzmh(3,ipt1) = 0.
c            xyzmh(1,ipt2) = -r2*dxdr
c            xyzmh(2,ipt2) = -r2*dydr
c            xyzmh(3,ipt2) = 0.
c            vxyzu(1,ipt1) = -v1*dydr
c            vxyzu(2,ipt1) = v1*dxdr
c            vxyzu(3,ipt1) = 0.
c            vxyzu(1,ipt2) = v2*dydr
c            vxyzu(2,ipt2) = -v2*dxdr
c            vxyzu(3,ipt2) = 0.
c            xmomsyn(1) = xyzmh(4,ipt1)*vxyzu(1,ipt1)
c            xmomadd(1) = 0.
c            ymomsyn(1) = xyzmh(4,ipt1)*vxyzu(2,ipt1)
c            ymomadd(1) = 0.
c            zmomsyn(1) = 0.
c            zmomadd(1) = 0.
c            xmomsyn(2) = xyzmh(4,ipt2)*vxyzu(1,ipt2)
c            xmomadd(2) = 0.
c            ymomsyn(2) = xyzmh(4,ipt2)*vxyzu(2,ipt2)
c            ymomadd(2) = 0.
c            zmomsyn(2) = 0.
c            zmomadd(2) = 0.
c         ENDIF
c
c--Modify orbit to be circular and have a semimajor axis changing analytically
c     from 1.0 down to xxx and mass ratio kept at 0.6
c
         IF (totmass.LT.1.08) THEN
            icenter = .true.
ccc            qratio = 0.6 + (totmass-1.0)/0.08 * 0.007
            qratio = 0.6 + (totmass-1.0)/0.08 * 0.035
            qratio1 = 1.0 + qratio
ccc            semiaxis = 1.0 - (totmass-1.0)/0.08 * 0.020
            semiaxis = 1.0 + (totmass-1.0)/0.08 * 0.015
            xyzmh(4,ipt1) = totmass/qratio1
            xyzmh(4,ipt2) = qratio*xyzmh(4,ipt1)
            r2 = semiaxis/qratio1
            r1 = qratio*r2
            v2 = SQRT(totmass/semiaxis)/qratio1
            v1 = qratio*v2
            dxdr = dx/dr
            dydr = dy/dr
            xyzmh(1,ipt1) = r1*dxdr
            xyzmh(2,ipt1) = r1*dydr
            xyzmh(3,ipt1) = 0.
            xyzmh(1,ipt2) = -r2*dxdr
            xyzmh(2,ipt2) = -r2*dydr
            xyzmh(3,ipt2) = 0.
            vxyzu(1,ipt1) = -v1*dydr
            vxyzu(2,ipt1) = v1*dxdr
            vxyzu(3,ipt1) = 0.
            vxyzu(1,ipt2) = v2*dydr
            vxyzu(2,ipt2) = -v2*dxdr
            vxyzu(3,ipt2) = 0.
            ptmsyn(1) = xyzmh(4,ipt1)
            ptmsyn(2) = xyzmh(4,ipt2)
            ptmadd(1) = 0.
            ptmadd(2) = 0.
            xmomsyn(1) = xyzmh(4,ipt1)*vxyzu(1,ipt1)
            xmomadd(1) = 0.
            ymomsyn(1) = xyzmh(4,ipt1)*vxyzu(2,ipt1)
            ymomadd(1) = 0.
            zmomsyn(1) = 0.
            zmomadd(1) = 0.
            xmomsyn(2) = xyzmh(4,ipt2)*vxyzu(1,ipt2)
            xmomadd(2) = 0.
            ymomsyn(2) = xyzmh(4,ipt2)*vxyzu(2,ipt2)
            ymomadd(2) = 0.
            zmomsyn(2) = 0.
            zmomadd(2) = 0.
         ELSEIF (icenter) THEN
            icenter = .false.
            cmvx = 0.
            cmvy = 0.
            cmvz = 0.
            simmass = 0.
            DO i = 1, npart
               IF (iphase(i).GE.0) THEN
                  cmvx = cmvx + xyzmh(4,i)*vxyzu(1,i)
                  cmvy = cmvy + xyzmh(4,i)*vxyzu(2,i)
                  cmvz = cmvz + xyzmh(4,i)*vxyzu(3,i)
                  simmass = simmass + xyzmh(4,i)
               ENDIF
            END DO
            cmvx = cmvx/simmass
            cmvy = cmvy/simmass
            cmvz = cmvz/simmass
            DO i = 1, npart
               IF (iphase(i).GE.0) THEN
                  vxyzu(1,i) = vxyzu(1,i) - cmvx
                  vxyzu(2,i) = vxyzu(2,i) - cmvy
                  vxyzu(3,i) = vxyzu(3,i) - cmvz
               ENDIF
            END DO
            xmomsyn(1) = xmomsyn(1) - xyzmh(4,ipt1)*cmvx
            ymomsyn(1) = ymomsyn(1) - xyzmh(4,ipt1)*cmvy
            zmomsyn(1) = zmomsyn(1) - xyzmh(4,ipt1)*cmvz
            xmomsyn(2) = xmomsyn(2) - xyzmh(4,ipt2)*cmvx
            ymomsyn(2) = ymomsyn(2) - xyzmh(4,ipt2)*cmvy
            zmomsyn(2) = zmomsyn(2) - xyzmh(4,ipt2)*cmvz
            WRITE (iprint,*)
            WRITE (iprint,*) 'ZEROING CENTER OF MASS VELOCITY'
            WRITE (iprint,*) '   ', cmvx, cmvy, cmvz, simmass
            WRITE (iprint,*)
         ENDIF
c
c--Move in the outer boundary
c
c         realtimesync = dt*itime/imaxstep
c         IF (realtimesync.GT.dt) THEN
c            boundnew = dmax1 + hma1*0.05
c            rmax = boundnew
c            xmin = -boundnew
c            xmax = boundnew
c            ymin = -boundnew
c            ymax = boundnew
c            zmin = -boundnew
c            zmax = boundnew
c         ENDIF
      ENDIF

c      PRINT *,'h(1)9: ',xyzmh(5,1),itime,f1ha(1,1),dumxyzmh(5,1)

c
c--Return or loop again
c
      IF (itime.GE.itnext) THEN
         DO i = 1, nbinmax
            nlstbins(i) = 0
            it2bin(i) = 2**i
         END DO

         IF (itiming) CALL getused(ts161)
         DO i = 1, npart
            IF (iphase(i).NE.-1) THEN
               it0(i) = 0
               it2(i) = isteps(i)

               ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
               IF (ibin.GT.nbinmax) THEN
                  WRITE (*,*) 'ERROR - ibin.GT.nbinmax 5'
                  WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax 5'
                  CALL quit(1)
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               listbins(nlstbins(ibin),ibin) = i
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin 5'
                  CALL quit(1)
               ENDIF
            ENDIF
         END DO
         IF (itiming) THEN
            CALL getused(ts162)
            ts16 = ts16 + (ts162 - ts161)
         ENDIF

         IF (igrape.EQ.0) THEN
            nneightot = 0
            DO i = 1, nlstnneigh
               ipart = llist(i)
               nneightot = nneightot + nneigh(ipart)
            END DO
         ENDIF

#ifdef MPI
c         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         print *,iproc,': END OF STEP ',nptmasstot
c         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

         RETURN
      ENDIF

      ilocal = ilocal + 1
      IF (MOD(ilocal,50).EQ.0) THEN
         ilocal = 0
c         CALL secmes
      ENDIF

      GOTO  100

      END

