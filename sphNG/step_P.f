       SUBROUTINE step (dt)
c************************************************************
c                                                           *
c  This subroutine integrate the system of differential     *
c     equations for one timestep using                      *
c     a second-order Runge-Kutta-Fehlberg method.           *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/f2'
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
      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/out1'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/sync'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/perform'
c      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/timeextra'
c      INCLUDE 'COMMONS/steplocal'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'

      DIMENSION nsteplist(30)
      INTEGER*8 nneightotsave

      LOGICAL ifirst
      LOGICAL icenter
      CHARACTER*7 where
      CHARACTER*5 ptdebug
      CHARACTER*2 itemchar
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
c--Set integrator to R-K
c
      integrator = 0
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
c
c--Define coefficients for Runge-Kutta integrator
c
      f21 = 1./256.
      f22 = 255./256.
      e1 = 1./512.
c
c--Integration
c
c---------- FIRST TIME AROUND ----------
c
      IF (ifirst) THEN
         ifirst = .FALSE.

         IF (integrator.EQ.0) THEN
            WRITE (iprint,*) 'Integrator: Runga-Kutta-Fehlberg'
            WRITE (iprint,*)
         ELSE
            WRITE (iprint,*) 'Integrator miss-match'
            CALL quit
         ENDIF

         IF (gt.EQ.0.0) THEN
            ibin = INT(LOG10(dt/dtini)/xlog2) + 1
            IF (ibin.LT.0) THEN
               WRITE(iprint,*) 'Error with initial timesteps'
               CALL quit
            ENDIF
            idtini = imaxstep/2**ibin
            istepmin = idtini
            istepmax = idtini
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

         nlst = 0
         DO i = 1, nbinmax
            nlstbins(i) = 0
            it1bin(i) = 2**i/2
            it2bin(i) = 2**i
         END DO
         DO i = 1, npart
            IF (iphase(i).NE.-1) THEN
               it0(i) = 0
ccc               isteps(i) = idtini
               IF (gt.EQ.0.0) isteps(i) = idtini
               it1(i) = isteps(i)/2
               it2(i) = isteps(i)
               ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
               IF (ibin.GT.nbinmax) THEN
                  WRITE (*,*) 'ERROR - ibin.GT.nbinmax'
                  WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax'
                  CALL quit
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               listbins(nlstbins(ibin),ibin) = i
               IF (it1bin(ibin).NE.it1(i)) THEN 
                  WRITE (*,*) 'ERROR - it1bin'
                  WRITE (iprint,*) 'ERROR - it1bin'
                  CALL quit
               ENDIF
               IF (it2bin(ibin).NE.it2(i)) THEN 
                  WRITE (*,*) 'ERROR - it2bin'
                  WRITE (iprint,*) 'ERROR - it2bin'
                  CALL quit
               ENDIF
               DO k = 1, 5
                  dumxyzmh(k,i) = xyzmh(k,i) 
               END DO
               DO k = 1, 4
                  dumvxyzu(k,i) = vxyzu(k,i) 
               END DO
               IF (encal.EQ.'r') THEN
                  DO k = 1, 5
                     dumekcle(k,i) = ekcle(k,i)
                  END DO
               ENDIF
               IF (imhd.EQ.idim) THEN
                  DO k = 1, 3
                     dumBevolxyz(k,i) = Bevolxyz(k,i)
                  END DO
               ENDIF

               IF (ifsvi.EQ.6) dumalpha(i) = alphaMM(i)

c               IF (ran1(1).LT.0.001) THEN

                  nlst = nlst + 1
                  llist(nlst) = i

c               ENDIF

               iscurrent(i) = .TRUE.
            ENDIF
         END DO
         jj = nlst
c
c--Create ghost particles
c
         nghost = 0
         IF (ibound.EQ.1) CALL ghostp1(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.2) CALL ghostp2(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.3 .OR. ibound.EQ.8 .OR. ibound/10.EQ.9) 
     &                    CALL ghostp3(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.100) 
     &        CALL ghostp100(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.11) 
     &        CALL ghostp11(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         ntot = npart + nghost

         print *,npart,nghost

         DO i = npart + 1, ntot
            DO k = 1, 5
               dumxyzmh(k,i) = xyzmh(k,i) 
            END DO
            DO k = 1, 4
               dumvxyzu(k,i) = vxyzu(k,i)
            END DO
            IF (encal.EQ.'r') THEN
               DO k = 1, 5
                  dumekcle(k,i) = ekcle(k,i)
               END DO
            ENDIF
            IF (imhd.EQ.idim) THEN
               DO k = 1, 3
                  dumBevolxyz(k,i) = Bevolxyz(k,i)
               END DO
            ENDIF
            IF (ifsvi.EQ.6) dumalpha(i) = alphaMM(ireal(i))
         END DO
         itime = 0
c
c--Make the tree
c
c      WRITE (*,*) 'passed 1'
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
c            CALL quit
c         ENDIF
c      END DO
c      WRITE (*,*) 'passed 2'

 432     CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,npart,
     &        ntot,ireal,dumalpha,dumekcle,dumBevolxyz,f1Bxyz)

c         IF (icall.EQ.1) THEN
c             icall = 3
c          ELSE 
c             STOP
c          ENDIF

c          GOTO 432





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
c            CALL quit
c         ENDIF
c      END DO
c
c--BEGIN: For calculating instantaneous torques about density maximum
c
c         rhomaxtorq = 0.
c         DO i = 1, npart
c            IF (rho(i).GT.rhomaxtorq) THEN
c               rhomaxtorq = rho(i)
c               itorqmax = i
c            ENDIF
c         END DO
c         DO i = 1, npart
c            dx = xyzmh(1,i) - xyzmh(1,itorqmax)
c            dy = xyzmh(2,i) - xyzmh(2,itorqmax)
c            torqt(i)=(dx*f1vxyzu(2,i) - dy*f1vxyzu(1,i))
c            torqg(i)=(dx*gravy(i) - dy*gravx(i))
c            torqp(i)= - (dx*gradpy(i) - 
c     &           dy*gradpx(i))*cnormk
c            torqv(i)= - (dx*artviy(i) - 
c     &           dy*artvix(i))*cnormk
c         END DO
c         CALL wdump(idisk1)
c         STOP
c
c--END: For calculating instantaneous torques about density maximum
c
         rhomaxsync = 0.
         DO i = 1, npart
            iscurrent(i) = .FALSE.
            rhomaxsync = MAX(rhomaxsync, rho(i))

            IF (encal.EQ.'r') THEN
               DO k = 1, 5
                  ekcle(k,i) = dumekcle(k,i)
               END DO
            ENDIF

c            gravx1(i) = gravx(i)
c            gravy1(i) = gravy(i)
c            gradpx1(i) = gradpx(i)
c            gradpy1(i) = gradpy(i)
c            artvix1(i) = artvix(i)
c            artviy1(i) = artviy(i)
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
c            CALL quit
c         ENDIF
c      END DO


      ELSE
         ntot = npart + nghost
      END IF
c
c---------- END FIRST TIME AROUND ----------
c
c--Find minimum time for force calculation
c
 100  itime1 = imax
      itime0 = imax
      itime1new = imax
      itime0new = imax
      IF (itiming) CALL getused(ts1p1)
      DO i = 1, nbinmax
         IF (nlstbins(i).GT.0) THEN
c            itime1new = MIN(itime1new, it1bin(i))
c            itime0new = MIN(itime0new, it2bin(i))
            itime1 = MIN(itime1, it1bin(i))
            itime0 = MIN(itime0, it2bin(i))
         ENDIF
      END DO
cC$OMP PARALLEL DO SCHEDULE(runtime) default(none)
cC$OMP& shared(npart,iphase,it0,it1,isteps)
cC$OMP& private(j)
cC$OMP& reduction(MIN:itime0new,itime1new)
c      DO j = 1, npart
c         IF (iphase(j).GE.0) THEN
c            itime1new = MIN(itime1new, it1(j))
c            itime0new = MIN(itime0new, it0(j) + isteps(j))
cc            itime1 = MIN(itime1, it1(j))
cc            itime0 = MIN(itime0, it0(j) + isteps(j))
c         ENDIF
c      END DO
cC$OMP END PARALLEL DO
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
c         WRITE (*,*) 2**i,it1bin(i),it2bin(i),nlstbins(i),nsteplist(i)
c         WRITE (iprint,*) 2**i,it1bin(i),it2bin(i),nlstbins(i),
c     &        nsteplist(i)
c        END DO
c        CALL quit
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
c            CALL quit
c         ENDIF
c      END DO

c
c--Check if predictions at half time step are needed
c
      IF (itime1.GE.itime0) THEN
         itime = itime0
         GOTO 101
      ELSE
         itime = itime1
      ENDIF
c
c--Allow for tracing flow
c
      IF(itrace.EQ.'all') WRITE(iprint,251) dt*itime/imaxstep+gt, itime
 251  FORMAT(' step is doing time ',1F12.5, I10)
c
c---------- PREDICTIONS AT HALF TIME STEP ----------
c
c--Identify particles with matching times and make a list
c
      IF (itiming) CALL getused(ts21)
      jj = 0
      jjnew = 0
      itbinupdate = 0
      DO i = 1, nbinmax
         IF (it1bin(i).EQ.itime) THEN
            itbinupdate = i
            DO j = 1, nlstbins(i)
               ipart = listbins(j,i)
               IF (ipart.GT.idim) THEN
                  WRITE (*,*) 'ERROR - ipart.GT.idim 1'
                  WRITE (iprint,*) 'ERROR - ipart.GT.idim 1'
                  CALL quit
               ELSEIF (ipart.LE.0) THEN
                  WRITE (*,*) 'ERROR - ipart.LE.0 1'
                  WRITE (iprint,*) 'ERROR - ipart.LE.0 1'
                  CALL quit
               ENDIF
               IF (iphase(ipart).GE.0) THEN
                  IF (it1(ipart).NE.itime) THEN
                     WRITE (*,*) 'ERROR - it1(ipart).NE.itime 1'
                     WRITE (iprint,*) 'ERROR - it1(ipart).NE.itime 1'
                     CALL quit
                  ENDIF
                  jj = jj + 1
                  IF (jj.GT.idim) THEN
                     WRITE (*,*) 'ERROR - jj.GT.idim 1'
                     WRITE (iprint,*) 'ERROR - jj.GT.idim 1'
                     CALL quit
                  ENDIF
                  llist(jj) = ipart
                  iscurrent(ipart) = .TRUE.
               ENDIF
            END DO
         ENDIF
      END DO
      itbinupdate = MIN(itbinupdate + 7, nbinmax)

c      DO i = 1, npart
c         IF (it1(i).EQ.itime .AND. iphase(i).GE.0) THEN
c            jjnew = jjnew + 1
cc            llist(jj) = i
cc            iscurrent(i) = .TRUE.
c         END IF
c      END DO
      nlst = jj
      nlst0 = 0

c      IF (jj.NE.jjnew) THEN
c         WRITE (*,*) 'ERROR - Half-timestep'
c         WRITE (iprint,*) 'ERROR - Half-timestep'
c         WRITE (*,*) 'Half-timestep:',jj,jjnew
c         WRITE (iprint,*) 'Half-timestep:',jj,jjnew
c         CALL quit
c      ENDIF

      IF (itiming) THEN
         CALL getused(ts22)
         ts2 = ts2 + (ts22 - ts21)
      ENDIF
c
c--Predict variables at t=time
c
      IF (itiming) CALL getused(ts51)

      IF (itbinupdate.GE.nbinmax-1 .OR. (.NOT. ipartialrevtree)) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,dt,itime,it0,imaxstep,ireal)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,ekcle)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,Bevolxyz,dumekcle)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,dumBevolxyz)
C$OMP& private(j,k,deltat)
      DO j = 1, npart
         IF (iphase(j).GE.0) THEN
            deltat = dt*(itime - it0(j))/imaxstep
            DO k = 1, 3
               dumxyzmh(k,j) = xyzmh(k,j) + deltat*vxyzu(k,j)
            END DO
            dumxyzmh(4,j) = xyzmh(4,j)
            dumxyzmh(5,j) = xyzmh(5,j) + deltat*f1ha(1,j)
            DO k = 1, 3
               dumvxyzu(k,j) = vxyzu(k,j) + deltat*f1vxyzu(k,j)
            END DO
            IF (encal.NE.'r' .AND. encal.NE.'m') THEN 
               dumvxyzu(4,j) = vxyzu(4,j) + deltat*f1vxyzu(4,j)
            ELSE
               dumvxyzu(4,j) = vxyzu(4,j)
               DO k = 1, 5
                  dumekcle(k,j) = ekcle(k,j)
               END DO
            ENDIF
            IF (iener.EQ.2 .AND. dumvxyzu(4,j).LT.0.0) 
     &           dumvxyzu(4,j)=0.15
            IF (ifsvi.EQ.6) dumalpha(j) = MIN(alphamax,alphaMM(j)+
     &           deltat*f1ha(2,j))
            IF (imhd.EQ.idim) THEN
               DO k = 1, 3
                  dumBevolxyz(k,j) = Bevolxyz(k,j) + deltat*f1Bxyz(k,j)
               END DO
            ENDIF
         ENDIF
      END DO
C$OMP END PARALLEL DO

      ELSE

         IF (nlst.GT.nptmass) THEN
            DO i = 1, itbinupdate
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(i,nlstbins,listbins,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,ekcle)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,Bevolxyz,dumekcle)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,dumBevolxyz)
C$OMP& private(j,k,ipart,deltat)
               DO j = 1, nlstbins(i)
                  ipart = listbins(j,i)
                  IF (iphase(ipart).GE.0) THEN
                     deltat = dt*(itime - it0(ipart))/imaxstep
                     DO k = 1, 3
            dumxyzmh(k,ipart) = xyzmh(k,ipart) + deltat*vxyzu(k,ipart)
                     END DO
            dumxyzmh(4,ipart) = xyzmh(4,ipart)
            dumxyzmh(5,ipart) = xyzmh(5,ipart) + deltat*f1ha(1,ipart)
                     DO k = 1, 3
            dumvxyzu(k,ipart) = vxyzu(k,ipart) + deltat*f1vxyzu(k,ipart)
                     END DO
                     IF (encal.NE.'r' .AND. encal.NE.'m') THEN
            dumvxyzu(4,ipart) = vxyzu(4,ipart) + deltat*f1vxyzu(4,ipart)
                     ELSE
                        dumvxyzu(4,ipart) = vxyzu(4,ipart)
                        DO k = 1, 5
                           dumekcle(k,ipart) = ekcle(k,ipart)
                        END DO
                     ENDIF
            IF (iener.EQ.2 .AND. dumvxyzu(4,ipart).LT.0.0) 
     &           dumvxyzu(4,ipart)=0.15
            IF (ifsvi.EQ.6) dumalpha(ipart) = MIN(alphamax,
     &               alphaMM(ipart)+deltat*f1ha(2,ipart))
                     IF (imhd.EQ.idim) THEN
                        DO k = 1, 3
                           dumBevolxyz(k,ipart) = Bevolxyz(k,ipart) 
     &                                          + deltat*f1Bxyz(k,ipart)
                        END DO
                     ENDIF
                  ENDIF
               END DO
C$OMP END PARALLEL DO
            END DO
         ELSE
c
c--Only update ptmasses, since nothing else is evolved
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nptmass,listpm,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha)
C$OMP& shared(dumxyzmh,dumvxyzu)
C$OMP& private(i,k,ipart,deltat)
            DO i = 1, nptmass
               ipart = listpm(i)
               deltat = dt*(itime - it0(ipart))/imaxstep
               DO k = 1, 3
              dumxyzmh(k,ipart) = xyzmh(k,ipart) + deltat*vxyzu(k,ipart)
               END DO
               dumxyzmh(4,ipart) = xyzmh(4,ipart)
               DO k = 1, 3
            dumvxyzu(k,ipart) = vxyzu(k,ipart) + deltat*f1vxyzu(k,ipart)
               END DO
            END DO
C$OMP END PARALLEL DO
         ENDIF
      ENDIF

      IF (nghost.GT.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,ireal,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,f1Bxyz,ekcle)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,Bevolxyz,dumekcle)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,dumBevolxyz)
C$OMP& private(j,k,l,deltat)
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
               DO l = 1, 5
                  dumekcle(l,j) = ekcle(l,j)
               END DO
            ENDIF
         IF (iener.EQ.2 .AND. dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
            IF (ifsvi.EQ.6) dumalpha(j) = MIN(alphamax,alphaMM(k) +
     &                                    deltat*f1ha(2,k))
            IF (imhd.EQ.idim) THEN
               DO l = 1, 3
                  dumBevolxyz(l,j) = Bevolxyz(l,j) + deltat*f1Bxyz(l,j)
               END DO
            ENDIF
         ENDIF
      END DO
C$OMP END PARALLEL DO
      ENDIF
      IF (itiming) THEN
         CALL getused(ts52)
         ts5 = ts5 + (ts52 - ts51)
      ENDIF

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
c--Make the tree or update the tree
c
ccc       IF (nactive.NE.nptmass) THEN  
      IF (igrape.EQ.0 .AND. (nlst.GT.nptmass .OR. 
     &                                       iptintree.GT.0)) THEN    
ccc         IF (nlst.GE.ncrit .OR. iaccr.EQ.1 .OR. ikilled.EQ.1) THEN
         IF (nlst.GE.ncrit) THEN
ccc            iaccr = 0
ccc            ikilled = 0
            CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)
         ELSEIF (iptintree.GT.0 .OR. nlst.GT.nptmass) THEN
            CALL insulate(2,ntot,npart,dumxyzmh,f1vxyzu)
         ENDIF
         iaccr = 0
         ikilled = 0
         nlstacc = 0
      ENDIF
c
c--Compute forces on list particles
c
      icall = 2
      CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f2vxyzu,f2ha,npart,
     &     ntot,ireal,dumalpha,dumekcle,dumBevolxyz,f2Bxyz)
c
c--Save velocities at half time step
c
      IF (itiming) CALL getused(ts61)
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst,llist,iscurrent,dumvxyzu)
C$OMP& shared(dum2vxyz,it1,imax)
C$OMP& private(i,j)
      DO j = 1, nlst
         i = llist(j)
         iscurrent(i) = .FALSE.
         dum2vxyz(1,i) = dumvxyzu(1,i)
         dum2vxyz(2,i) = dumvxyzu(2,i)
         dum2vxyz(3,i) = dumvxyzu(3,i)
         it1(i) = imax
      END DO
C$OMP END PARALLEL DO
      DO i = 1, nbinmax
         IF (it1bin(i).EQ.itime) it1bin(i) = imax
      END DO
      IF (itiming) THEN
         CALL getused(ts62)
         ts6 = ts6 + (ts62 - ts61)
      ENDIF

      GOTO 100
c
c---------- ADVANCING PARTICLES ----------
c
c--Identify particles to be advanced
c
 101  jj = 0
      imakeghost = 0
c      write (*,*) '5 ',itime
c
c
c--BEGIN: Find maximum density for calculating torques
      imaxdens = 0
c      densmax = 0.
c      DO i = 1, npart
c         IF (iphase(i).GE.0 .AND. rho(i).GT.densmax) THEN
c            densmax = rho(i)
c            imaxdens = i
c         ENDIF
c      END DO
c--END:  Find maximum density for calculating torques
c
c
      IF (itiming) CALL getused(ts31)
      jjnew = 0
      itbinupdate = 0
      DO i = 1, nbinmax
         IF (it2bin(i).EQ.itime) THEN
            itbinupdate = i
            DO j = 1, nlstbins(i)
               ipart = listbins(j,i)
               IF (ipart.GT.idim) THEN
                  WRITE (*,*) 'ERROR - ipart.GT.idim 2'
                  WRITE (iprint,*) 'ERROR - ipart.GT.idim 2'
                  CALL quit
               ELSEIF (ipart.LE.0) THEN
                  WRITE (*,*) 'ERROR - ipart.LE.0 2'
                  WRITE (iprint,*) 'ERROR - ipart.LE.0 2'
                  CALL quit
               ENDIF
               IF (iphase(ipart).GE.0) THEN
                  IF (it2(ipart).NE.itime) THEN
                     WRITE (*,*) 'ERROR - it2(ipart).NE.itime 2'
                     WRITE (iprint,*) 'ERROR - it2(ipart).NE.itime 2'
                     WRITE (*,*) ipart,it2(ipart),itime,it2bin(i)
                     WRITE (*,*) it1(ipart),it0(ipart)
                     CALL quit
                  ENDIF
                  IF (it0(ipart)+isteps(ipart).NE.itime) THEN
                     WRITE (*,*) 'ERROR - it0+it1.NE.itime 2'
                     WRITE (iprint,*) 'ERROR - it0+it1.NE.itime 2'
                     CALL quit
                  ENDIF

                  jj = jj + 1
                  IF (jj.GT.idim) THEN
                     WRITE (*,*) 'ERROR - jj.GT.idim 2'
                     WRITE (iprint,*) 'ERROR - jj.GT.idim 2'
                     CALL quit
                  ENDIF
                  llist(jj) = ipart
                  IF (hasghost(ipart)) imakeghost = 1
               ENDIF
            END DO
         ENDIF
      END DO

c      DO i = 1, npart
c         IF (it0(i)+isteps(i).EQ.itime .AND. iphase(i).GE.0) THEN
c            jjnew = jjnew + 1
cc            llist(jj) = i
cc            IF (hasghost(i)) imakeghost = 1
c         ENDIF
c      END DO
      nlst0 = jj

c      IF (jj.NE.jjnew) THEN
c         WRITE (*,*) 'ERROR - Ad'
c         WRITE (iprint,*) 'ERROR - Ad'
c         WRITE (*,*) 'Advance:',jj,jjnew
c         WRITE (iprint,*) 'Advance:',jj,jjnew
c      ENDIF

      IF (itiming) THEN
         CALL getused(ts32)
         ts3 = ts3 + (ts32 - ts31)
      ENDIF

      IF (itiming) CALL getused(ts71)
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst0,llist,iscurrent,dt,isteps,imaxstep,f21,f22)
C$OMP& shared(xyzmh,vxyzu,dum2vxyz)
C$OMP& shared(f1vxyzu,f1ha,f2vxyzu,f2ha,f1Bxyz,f2Bxyz)
C$OMP& shared(gravx,gravy,gradpx,gradpy,artvix,artviy,gravx1,gravy1)
C$OMP& shared(gradpx1,gradpy1,artvix1,artviy1,torqt,torqg,torqp)
C$OMP& shared(torqv,torqc,it0,itime,imaxdens,cnormk,iener)
C$OMP& shared(iprint,nneigh,hmaximum,iphase,nptmass,listpm)
C$OMP& shared(xmomsyn,ymomsyn,zmomsyn,pmass,dumekcle)
C$OMP& shared(ifsvi,alphaMM,alphamax,encal,dumBevolxyz,Bevolxyz)
C$OMP& private(i,j,dtfull,dtf21,dtf22,xold,yold,delvx,delvy,delvz)
C$OMP& private(vxstore,vystore,dx,dy,delgx,delgy,delpx,delpy)
C$OMP& private(delax,delay,extraxt,extrayt,extraxg,extrayg)
C$OMP& private(extraxp,extrayp,extraxv,extrayv,dtf22dtfull)
C$OMP& private(iii,pmasspt)
C$OMP& reduction(+:ioutmax)
      DO j = 1, nlst0
         i = llist(j)
         IF (it0(i)+isteps(i).EQ.itime .AND. iphase(i).GE.0) THEN
            iscurrent(i) = .TRUE.
c
c--Update positions
c
            dtfull = (dt*isteps(i))/imaxstep
            dtf21 = dtfull*f21
            dtf22 = dtfull*f22

            xold = xyzmh(1,i)
            yold = xyzmh(2,i)

            xyzmh(1,i) = xyzmh(1,i) + dtf21*vxyzu(1,i) + 
     &           dtf22*dum2vxyz(1,i)
            xyzmh(2,i) = xyzmh(2,i) + dtf21*vxyzu(2,i) + 
     &           dtf22*dum2vxyz(2,i)
            xyzmh(3,i) = xyzmh(3,i) + dtf21*vxyzu(3,i) + 
     &           dtf22*dum2vxyz(3,i)
c
c--Keep old velocities for error calculation
c
            dum2vxyz(1,i) = vxyzu(1,i)
            dum2vxyz(2,i) = vxyzu(2,i)
            dum2vxyz(3,i) = vxyzu(3,i)
c
c--Update velocities
c
            delvx = dtf21*f1vxyzu(1,i) + dtf22*f2vxyzu(1,i)
            delvy = dtf21*f1vxyzu(2,i) + dtf22*f2vxyzu(2,i)
            delvz = dtf21*f1vxyzu(3,i) + dtf22*f2vxyzu(3,i)

            vxstore = vxyzu(1,i)
            vystore = vxyzu(2,i)

            vxyzu(1,i) = vxyzu(1,i) + delvx
            vxyzu(2,i) = vxyzu(2,i) + delvy
            vxyzu(3,i) = vxyzu(3,i) + delvz
            IF (iphase(i).GE.1) THEN
               DO iii = 1, nptmass
                  IF (listpm(iii).EQ.i) GOTO 444
               END DO
 444           CONTINUE
               pmasspt = xyzmh(4,i)
               xmomsyn(iii) = xmomsyn(iii) + delvx*pmasspt
               ymomsyn(iii) = ymomsyn(iii) + delvy*pmasspt
               zmomsyn(iii) = zmomsyn(iii) + delvz*pmasspt
            ENDIF
c
c--Save torques - derived by making the torques exactly equal to the change in 
c    angular momentum of the particle over the timestep.  The total torque
c    is also evaluated by a second method below as a check.
c
            IF (imaxdens.EQ.0) THEN
               dx = xyzmh(1,i)
               dy = xyzmh(2,i)
            ELSE
               dx = xyzmh(1,i) - xyzmh(1,imaxdens)
               dy = xyzmh(2,i) - xyzmh(2,imaxdens)
               xold = xold - xyzmh(1,imaxdens)
               yold = yold - xyzmh(2,imaxdens)
            ENDIF

c            delgx = dtf21*gravx1(i) + dtf22*gravx(i)
c            delgy = dtf21*gravy1(i) + dtf22*gravy(i)
c            delpx = dtf21*gradpx1(i) + dtf22*gradpx(i)
c            delpy = dtf21*gradpy1(i) + dtf22*gradpy(i)
c            delax = dtf21*artvix1(i) + dtf22*artvix(i)
c            delay = dtf21*artviy1(i) + dtf22*artviy(i)

c            dtf22dtfull = dtf22*dtfull/2.0

c            extraxt = dtf22dtfull*f1vxyzu(1,i)
c            extrayt = dtf22dtfull*f1vxyzu(2,i)
c            extraxg = dtf22dtfull*gravx1(i)
c            extrayg = dtf22dtfull*gravy1(i)
c            extraxp = dtf22dtfull*gradpx1(i)
c            extrayp = dtf22dtfull*gradpy1(i)
c            extraxv = dtf22dtfull*artvix1(i)
c            extrayv = dtf22dtfull*artviy1(i)

c            torqt(i)=torqt(i) + (dx*delvy - dy*delvx)
c     &           + (extraxt*vystore - extrayt*vxstore)
c            torqg(i)=torqg(i) + (dx*delgy - dy*delgx)
c     &           + (extraxg*vystore - extrayg*vxstore)
c            torqp(i)=torqp(i) - (dx*delpy - dy*delpx)*cnormk
c     &           - (extraxp*vystore - extrayp*vxstore)*cnormk
c            torqv(i)=torqv(i) - (dx*delay - dy*delax)*cnormk
c     &           - (extraxv*vystore - extrayv*vxstore)*cnormk
c
c--Check torques by using a less exact method - this method is the same
c    as above, but is derived using T = 1/256*T_1 + 255/256*T_2
c    (i.e. the same as for the intergration).  The above method
c    makes the torques *exactly* equal to the change in specific angular
c    momentum.  The two methods are identical in the limit that f1 and f2
c    are equal.  Hence this provides a check that the timestep is small
c    enough to give a good integration.
c
c            torqc(i)=torqc(i) + dtf21*(xold*f1vxyzu(2,i) - 
c     &            yold*f1vxyzu(1,i))
c     &           + dtf22*(xold*f2vxyzu(2,i) - yold*f2vxyzu(1,i))
c     &           + (vxstore*dtf22dtfull*f2vxyzu(2,i) - 
c     &           vystore*dtf22dtfull*f2vxyzu(1,i))
c
c--Update u(i) and h(i)
c
            IF (encal.NE.'r' .AND. encal.NE.'m') 
     &           vxyzu(4,i) = vxyzu(4,i) + 
     &           dtf21*f1vxyzu(4,i) + dtf22*f2vxyzu(4,i)

            IF (iener.EQ.2 .AND. vxyzu(4,i).LT.0.0) vxyzu(4,i)=0.15
            xyzmh(5,i) = xyzmh(5,i) + dtf21*f1ha(1,i) + dtf22*f2ha(1,i)
            IF (xyzmh(5,i).LT.0) WRITE(iprint,*) 'error in h ', 
     &           xyzmh(5,i),f1ha(1,i),dtf21,f2ha(1,i),dtf22,nneigh(i)
            IF (hmaximum.GT.0.0) THEN
               IF (xyzmh(5,i).GT.hmaximum) THEN
                  xyzmh(5,i) = hmaximum
                  ioutmax = ioutmax + 1
               ENDIF
            ENDIF
c
c--Update viscosity switch
c
            IF (ifsvi.EQ.6) alphaMM(i) = MIN(alphamax, alphaMM(i) +
     &           dtf21*f1ha(2,i) + dtf22*f2ha(2,i))
c
c-MHD
c
            IF (imhd.EQ.idim) THEN
               Bevolxyz(1,i) = Bevolxyz(1,i)
     &              + dtf21*f1Bxyz(1,i) + dtf22*f2Bxyz(1,i)
               Bevolxyz(2,i) = Bevolxyz(2,i)
     &              + dtf21*f1Bxyz(2,i) + dtf22*f2Bxyz(2,i)
               Bevolxyz(3,i) = Bevolxyz(3,i)
     &              + dtf21*f1Bxyz(3,i) + dtf22*f2Bxyz(3,i)
            ENDIF
c
c--Synchronize the advanced particle times with current time
c
            it0(i) = itime
         END IF
      END DO
C$OMP END PARALLEL DO
      IF (itiming) THEN
         CALL getused(ts72)
         ts7 = ts7 + (ts72 - ts71)
      ENDIF
c
c--Total number of advanced particles
c
      idonebound = 0
      IF (ibound.GT.0 .AND. ibound.LT.7 .OR. ibound.EQ.11) 
     &     CALL boundry(npart,
     &     llist, nlst0, xyzmh, vxyzu, idonebound)
c
c--Create ghost particles
c
      IF ((nlst0.GT.ncrit .AND. imakeghost.EQ.1) 
     &                                    .OR. idonebound.EQ.1) THEN
c      IF (imakeghost.EQ.1 .OR. idonebound.EQ.1) THEN
         nghost = 0
         IF (ibound.EQ.1) CALL ghostp1(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.2) CALL ghostp2(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.3 .OR. ibound.EQ.8 .OR. ibound/10.EQ.9)
     &                    CALL ghostp3(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.100)
     &        CALL ghostp100(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         IF (ibound.EQ.11)
     &        CALL ghostp11(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
         ntot = npart + nghost
      ENDIF
c
c--Identify particles to predict
c
      IF (itiming) CALL getused(ts41)
      DO i = 1, nbinmax
         IF (it1bin(i).EQ.itime) THEN
            itbinupdate = MAX(itbinupdate,i)
            DO j = 1, nlstbins(i)
               ipart = listbins(j,i)
               IF (ipart.GT.idim) THEN
                  WRITE (*,*) 'ERROR - ipart.GT.idim 3'
                  WRITE (iprint,*) 'ERROR - ipart.GT.idim 3'
                  CALL quit
               ELSEIF (ipart.LE.0) THEN
                  WRITE (*,*) 'ERROR - ipart.LE.0 3'
                  WRITE (iprint,*) 'ERROR - ipart.LE.0 3'
                  WRITE (*,*) ipart,i,j,itime,it1bin(i)
                  WRITE (iprint,*) ipart,i,j,itime,it1bin(i)
                  WRITE (*,*) nlst0,jj
                  WRITE (iprint,*) nlst0,jj
                  WRITE (*,*) nlstbins(i),listbins(j,i)
                  WRITE (iprint,*) nlstbins(i),listbins(j,i)

                     DO j2 = 1, 30
                        nsteplist(j2) = 0
                     END DO

                     DO j2 = 1, npart
                        IF (iphase(j2).GE.0) THEN
                        ibin = INT(LOG10(REAL(isteps(j2)))/xlog2+0.5)
                           IF ((ibin.GE.1) .AND. (ibin.LE.29)) THEN
                              nsteplist(ibin) = nsteplist(ibin) + 1
                           ELSE
                              nsteplist(30) = nsteplist(30) + 1
                           ENDIF
                        ENDIF
                     END DO

         DO i2 = 1, nbinmax
         WRITE (*,*) 2**i2,it1bin(i2),it2bin(i2),nlstbins(i2),
     &           nsteplist(i2)
         WRITE (iprint,*) 2**i2,it1bin(i2),it2bin(i2),nlstbins(i2),
     &        nsteplist(i2)
         END DO

         WRITE (*,*) 'list half-tstep bin',i,nlstbins(i)
         WRITE (iprint,*) 'list half-tstep bin',i,nlstbins(i)
                  DO i2 = 1,nlstbins(i)
         WRITE (iprint,*) i2,listbins(i2,i)
                  END DO

                  CALL quit
               ENDIF
               IF (iphase(ipart).GE.0) THEN
                  IF (it1(ipart).NE.itime) THEN
                     WRITE (*,*) 'ERROR - it1(ipart).NE.itime 3'
                     WRITE (iprint,*) 'ERROR - it1(ipart).NE.itime 3'
                     WRITE (*,*) i,ipart,iphase(ipart),iorig(ipart)
                WRITE (iprint,*) i,ipart,iphase(ipart),iorig(ipart)
                     WRITE (*,*) it1bin(i),it1(ipart),itime,j,nlst0
                WRITE (iprint,*) it1bin(i),it1(ipart),itime,j,nlst0
                     WRITE (*,*) it0(ipart),it2(ipart),nlstbins(i),jj
                WRITE (iprint,*) it0(ipart),it2(ipart),nlstbins(i),jj

                     DO j2 = 1, 30
                        nsteplist(j2) = 0
                     END DO

                     DO j2 = 1, npart
                        IF (iphase(j2).GE.0) THEN
                        ibin = INT(LOG10(REAL(isteps(j2)))/xlog2+0.5)
                           IF ((ibin.GE.1) .AND. (ibin.LE.29)) THEN
                              nsteplist(ibin) = nsteplist(ibin) + 1
                           ELSE
                              nsteplist(30) = nsteplist(30) + 1
                           ENDIF
                        ENDIF
                     END DO

         DO i2 = 1, nbinmax
         WRITE (*,*) 2**i2,it1bin(i2),it2bin(i2),nlstbins(i2),
     &           nsteplist(i2)
         WRITE (iprint,*) 2**i2,it1bin(i2),it2bin(i2),nlstbins(i2),
     &        nsteplist(i2)
         END DO


                     CALL quit
                  ENDIF
                  jj = jj + 1
                  IF (jj.GT.idim) THEN
                     WRITE (*,*) 'ERROR - jj.GT.idim 3'
                     WRITE (iprint,*) 'ERROR - jj.GT.idim 3'
                     CALL quit
                  ENDIF
                  llist(jj) = ipart
                  iscurrent(ipart) = .TRUE.
               ENDIF
            END DO
         ENDIF
      END DO
      itbinupdate = MIN(itbinupdate + 7, nbinmax)

c      DO j = 1, npart
c         IF (it1(j).EQ.itime) THEN
c            IF (iphase(j).GE.0) THEN
c               jjnew = jjnew + 1
cc               llist(jj) = j
cc               iscurrent(j) = .TRUE.
c            ENDIF
c         ENDIF
c      END DO

c      IF (jj.NE.jjnew) THEN
c         WRITE (*,*) 'ERROR - Half-timestep 2'
c         WRITE (iprint,*) 'ERROR - Half-timestep 2'
c         WRITE (*,*) 'Half-timestep 2:',jj,jjnew
c         WRITE (iprint,*) 'Half-timestep 2:',jj,jjnew
c         CALL quit
c      ENDIF

      IF (itiming) THEN
         CALL getused(ts42)
         ts4 = ts4 + (ts42 - ts41)
      ENDIF
c
c--Total number of particles to compute forces on
c
      nlst = jj
      nlstnneigh = nlst
c
c--Total number of particles on which forces are to be predicted
c
      nlst1 = nlst - nlst0
c
c--Predict variables at t=time
c
      IF (itiming) CALL getused(ts81)

      IF (itbinupdate.GE.nbinmax-1 .OR. (.NOT. ipartialrevtree)) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,dt,itime,it0,imaxstep,ireal)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,dumekcle,ekcle,f1Bxyz)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,dumBevolxyz)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,Bevolxyz)
C$OMP& private(j,k,deltat)
      DO j = 1, npart
         IF (iphase(j).GE.0) THEN
            deltat = dt*(itime - it0(j))/imaxstep
            DO k = 1, 3
               dumxyzmh(k,j) = xyzmh(k,j) + deltat*vxyzu(k,j)
            END DO
            dumxyzmh(4,j) = xyzmh(4,j)
            dumxyzmh(5,j) = xyzmh(5,j) + deltat*f1ha(1,j)
            DO k = 1, 3
               dumvxyzu(k,j) = vxyzu(k,j) + deltat*f1vxyzu(k,j)
            END DO
            IF (encal.NE.'r' .AND. encal.NE.'m') THEN
               dumvxyzu(4,j) = vxyzu(4,j) + deltat*f1vxyzu(4,j)
            ELSE
               dumvxyzu(4,j) = vxyzu(4,j)
               DO k = 1, 5
                  dumekcle(k,j) = ekcle(k,j)
               END DO
            ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
            IF (ifsvi.EQ.6) dumalpha(j) = MIN(alphamax,alphaMM(j)+
     &           deltat*f1ha(2,j))
            IF (imhd.EQ.idim) THEN
               DO k = 1, 3
                  dumBevolxyz(k,j) = Bevolxyz(k,j) + deltat*f1Bxyz(k,j)
               END DO
            ENDIF
         ENDIF
      END DO
C$OMP END PARALLEL DO

      ELSE

         IF (nlst.GT.nptmass) THEN
            DO i = 1, itbinupdate
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(i,nlstbins,listbins,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,dumekcle,ekcle,f1Bxyz)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,dumBevolxyz)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,Bevolxyz)
C$OMP& private(j,k,ipart,deltat)
               DO j = 1, nlstbins(i)
                  ipart = listbins(j,i)
                  IF (iphase(ipart).GE.0) THEN
                     deltat = dt*(itime - it0(ipart))/imaxstep
                     DO k = 1, 3
             dumxyzmh(k,ipart) = xyzmh(k,ipart) + deltat*vxyzu(k,ipart)
                     END DO
                     dumxyzmh(4,ipart) = xyzmh(4,ipart)
              dumxyzmh(5,ipart) = xyzmh(5,ipart) + deltat*f1ha(1,ipart)
                     DO k = 1, 3
           dumvxyzu(k,ipart) = vxyzu(k,ipart) + deltat*f1vxyzu(k,ipart)
                     END DO
                     IF (encal.NE.'r' .AND. encal.NE.'m') THEN
           dumvxyzu(4,ipart) = vxyzu(4,ipart) + deltat*f1vxyzu(4,ipart)
                     ELSE
                        dumvxyzu(4,ipart) = vxyzu(4,ipart)
                        DO k = 1, 5
                           dumekcle(k,ipart) = ekcle(k,ipart)
                        END DO
                     ENDIF
                     IF (iener.EQ.2 .AND. dumvxyzu(4,ipart).LT.0.0) 
     &                    dumvxyzu(4,ipart)=0.15
                     IF (ifsvi.EQ.6) dumalpha(ipart) = MIN(alphamax,
     &                    alphaMM(ipart) + deltat*f1ha(2,ipart))
                     IF (imhd.EQ.idim) THEN
                        DO k = 1, 3
                           dumBevolxyz(k,ipart) = Bevolxyz(k,ipart) 
     &                                          + deltat*f1Bxyz(k,ipart)
                        END DO
                     ENDIF
                  ENDIF
               END DO
C$OMP END PARALLEL DO
            END DO
         ELSE
c
c--Only update ptmasses, since nothing else is evolved
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nptmass,listpm,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha)
C$OMP& shared(dumxyzmh,dumvxyzu)
C$OMP& private(i,k,ipart,deltat)
            DO i = 1, nptmass
               ipart = listpm(i)
               deltat = dt*(itime - it0(ipart))/imaxstep
               DO k = 1, 3
            dumxyzmh(k,ipart) = xyzmh(k,ipart) + deltat*vxyzu(k,ipart)
               END DO
               dumxyzmh(4,ipart) = xyzmh(4,ipart)
               dumxyzmh(5,ipart) = xyzmh(5,ipart)
               DO k = 1, 3
            dumvxyzu(k,ipart) = vxyzu(k,ipart) + deltat*f1vxyzu(k,ipart)
               END DO
            END DO
C$OMP END PARALLEL DO
         ENDIF
      ENDIF

      IF (nghost.GT.0) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,nghost,ireal,dt,itime,it0,imaxstep)
C$OMP& shared(xyzmh,vxyzu,f1vxyzu,f1ha,dumekcle,ekcle,f1Bxyz)
C$OMP& shared(dumxyzmh,dumvxyzu,iphase,iener,dumBevolxyz)
C$OMP& shared(ifsvi,dumalpha,alphaMM,alphamax,encal,Bevolxyz)
C$OMP& private(j,k,l,deltat)
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
               DO l = 1, 5
                  dumekcle(l,j) = ekcle(l,j)
               END DO
            ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
            IF (ifsvi.EQ.6) dumalpha(j) = MIN(alphamax, alphaMM(k) 
     &           + deltat*f1ha(2,k))
            IF (imhd.EQ.idim) THEN
               DO l = 1, 3
                  dumBevolxyz(l,j) = Bevolxyz(l,j) + deltat*f1Bxyz(l,j)
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
ccc      IF (nactive.NE.nptmass) THEN    
      IF (igrape.EQ.0 .AND. (nlst.GT.nptmass .OR. 
     &                                       iptintree.GT.0)) THEN
         IF (nlst.GT.ncrit .OR.(nlst0.GT.ncrit.AND.imakeghost.EQ.1).OR.
ccc     &        idonebound.EQ.1 .OR. iaccr.EQ.1 .OR. ikilled.EQ.1) THEN
     &        idonebound.EQ.1) THEN
ccc            iaccr = 0
ccc            ikilled = 0
            CALL insulate(1, ntot, npart, dumxyzmh, f1vxyzu)
         ELSEIF (iptintree.GT.0 .OR. nlst.GT.nptmass) THEN
            CALL insulate(2, ntot, npart, dumxyzmh, f1vxyzu)
         ENDIF
         iaccr = 0
         ikilled = 0
         nlstacc = 0
      ENDIF
c
c--Compute forces on list particles - Note: May be full AND half time step
c     evaluations!!
c
 200  icall = 3
      CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f2vxyzu,f2ha,npart,
     &     ntot,ireal,dumalpha,dumekcle,dumBevolxyz,f2Bxyz)
c
c--Synchronization time
c
      idtsyn = itnext - itime
      IF (idtsyn.EQ.0) idtsyn = imaxstep
      ikilled = 0
      istminold = istepmin
      istepmin = imax
      istepmax = 0
      iptnum = 0
      time = dt*itime/imaxstep + gt
c
c--Keep velocities at half time step
c
      IF (itiming) CALL getused(ts91)
      DO i = 1, nbinmax
         IF (it1bin(i).EQ.itime) it1bin(i) = imax
      END DO

C$OMP PARALLEL default(none)
C$OMP& shared(nlst0,nlst1,llist,iscurrent,dumvxyzu)
C$OMP& shared(dum2vxyz,it1,it2,imax,igphi,encal)
C$OMP& shared(vxyzu,dgrav,dumxyzmh)
C$OMP& shared(f1vxyzu,f1ha,f2vxyzu,f2ha)
C$OMP& shared(tol,tolh,tolptm,dt,isteps,imaxstep,e1,small,divv,rho)
C$OMP& shared(alpha,beta,vsound,xlog2,gravx,gravy,gravz,gradpx,gradpy)
C$OMP& shared(gradpz,artvix,artviy,artviz,gravx1,gravy1,gradpx1,gradpy1)
C$OMP& shared(artvix1,artviy1,it0,idtsyn,ibound,deadbound)
C$OMP& shared(iprint,xyzmh,poten,nneigh,iphase,ikillpr)
C$OMP& shared(pmass,ikilled,nactive,nkill,time,iorig,nlstacc,listacc)
C$OMP& shared(anglostx,anglosty,anglostz,istminold)
C$OMP& shared(ifsvi,alphaMM,dumekcle,ekcle)
C$OMP& private(i,j,l,tolpart,errBx,errBy,errBz)
C$OMP& private(errx,erry,errz,errvx,errvy,errvz,erru,errh,errm)
C$OMP& private(errdivtol,rap,rmod1,divvi,aux1,aux2,aux3,denom)
C$OMP& private(crstepi,rmodcr,force21,force22,force2,rmodcr2,rmod)
C$OMP& private(stepi,ibin,factor,istep2,irat,r2,pmassi)
C$OMP& reduction(MAX:istepmax)
C$OMP& reduction(MIN:istepmin)

C$OMP DO SCHEDULE(runtime)
      DO j = nlst0 + 1, nlst0 + nlst1
         i = llist(j)
         iscurrent(i) = .FALSE.
         dum2vxyz(1,i) = dumvxyzu(1,i)
         dum2vxyz(2,i) = dumvxyzu(2,i)
         dum2vxyz(3,i) = dumvxyzu(3,i)
c
c--Label particles on which predicted forces have been computed
c
         it1(i) = imax
      END DO
C$OMP END DO
c
c--Compute new time steps for particles which have been advanced
c
C$OMP DO SCHEDULE(runtime)
      DO 855 j = 1, nlst0
         i = llist(j)
         iscurrent(i) = .FALSE.
c
c--Dead particle boundaries
c
         IF (ibound.EQ.8 .OR. ibound/10.EQ.9) THEN
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            IF (r2.GT.(1.05*deadbound)**2) THEN
               iphase(i) = -1
               pmassi = xyzmh(4,i)
C$OMP CRITICAL (killparticle)
               ikilled = 1

               nlstacc = nlstacc + 1
               IF (nlstacc.GT.nlstaccmax) THEN
                  WRITE (iprint,*) 'ERROR step nlstacc'
                  CALL quit
               ENDIF
               listacc(nlstacc) = i

               nactive = nactive - 1
               nkill = nkill + 1
               WRITE (ikillpr) iorig(i),time,xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              poten(i),dgrav(i)
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
c--Follow point masses with greater accuracy than the gas particles 
c
         IF (iphase(i).GE.1) THEN
            tolpart = tolptm
         ELSE
            tolpart = tol
         ENDIF
c
c--Set u to it's new value from DERIVI 
c     For polytropic equation of state, the du's are not used - u(i) 
c       is calculated directly from the density, in each derivi call and
c       put into dumu(i). Hence must be transferred from dumu(i) to u(i).
c     For radiative transfer, both u(i) and e(i) are set in DERIVI, and
c       hence u(i) and e(i) need to be copied back here.
c     For adiabatic (or isothermal) u(i) is calculated via the du's
c       but this setting of u(i)=dumu(i) doesn't matter as the u(i)
c       updated at the full timestep above, then put into dumu(i)
c       but the derivi call doesn't alter them, so putting them back
c       into u(i) again changes nothing.
c
         vxyzu(4,i) = dumvxyzu(4,i)
         IF (encal.EQ.'r') THEN
            DO l = 1, 5
               ekcle(l,i) = dumekcle(l,i)
            END DO
         ENDIF
c
c--nlmax=1 is the sign that code is running using 'grad-h' so that h is set
c     inside derivi rather than being evolved.
c
         IF (nlmax.EQ.1) xyzmh(5,i) = dumxyzmh(5,i)
c
c--Estimate maximum error
c
         errx = ABS(vxyzu(1,i) - dum2vxyz(1,i))
         erry = ABS(vxyzu(2,i) - dum2vxyz(2,i))
         errz = ABS(vxyzu(3,i) - dum2vxyz(3,i))
         errvx = ABS(f1vxyzu(1,i) - f2vxyzu(1,i))
         errvy = ABS(f1vxyzu(2,i) - f2vxyzu(2,i))
         errvz = ABS(f1vxyzu(3,i) - f2vxyzu(3,i))
         erru = ABS(f1vxyzu(4,i) - f2vxyzu(4,i))
         errh = 3.0*ABS(f1ha(1,i) - f2ha(1,i))
c
c--Find maximum error
c    
c         IF (imhd.EQ.idim) THEN
c            errBx = ABS(f1Bxyz(1,i) - f2Bxyz(1,i))
c            errBy = ABS(f1Bxyz(2,i) - f2Bxyz(2,i))
c            errBz = ABS(f1Bxyz(3,i) - f2Bxyz(3,i))
c            errm = MAX(errx,erry,errz,errvx,errvy,errvz,erru,
c     &                 errBx,errBy,errBz)
c         ELSE
            errm = MAX(errx,erry,errz,errvx,errvy,errvz,erru)
c         ENDIF
         errdivtol = MAX(errm/tolpart, errh/tolh)
c
c--Compute ratio of present to next timestep
c
         rap = dt*isteps(i)/imaxstep*errdivtol*e1 + small
         rmod1 = 1./SQRT(rap)
c
c--Compute the time step for the gas from the Courant condition
c
         IF (iphase(i).EQ.0) THEN
            divvi = divv(i)/rho(i)
            IF (ifsvi.EQ.6) THEN
               aux1 = alphaMM(i)*vsound(i)
               aux2 = xyzmh(5,i)*ABS(divvi)
               aux3 = aux2
               IF (divvi.LT.0.0) THEN
                  aux2 = 2.0*alphaMM(i)*aux2
               ELSE
                  aux2 = 0.0
               ENDIF
            ELSE
               aux1 = alpha*vsound(i)
               aux2 = xyzmh(5,i)*ABS(divvi)
               aux3 = aux2
               IF (divvi.LT.0.0) THEN
                  aux2 = beta*aux2
               ELSE
                  aux2 = 0.0
               ENDIF
            ENDIF
            denom = aux3 + vsound(i) + 1.2*(aux1 + aux2)
            crstepi = 0.3*xyzmh(5,i)/denom
            rmodcr = crstepi/(dt*isteps(i)/imaxstep)
            force21 = f1vxyzu(1,i)**2 + f1vxyzu(2,i)**2 + 
     &                                              f1vxyzu(3,i)**2
            force22 = f2vxyzu(1,i)**2 + f2vxyzu(2,i)**2 + 
     &                                              f2vxyzu(3,i)**2
            force2 = MAX(force21, force22)
            rmodcr2 = 0.3*SQRT(xyzmh(5,i)/SQRT(force2))
            rmodcr2 = rmodcr2/(dt*isteps(i)/imaxstep)
c
c--Decide which time step to take
c
            rmod = MIN(rmod1, rmodcr, rmodcr2)
         ELSE
            rmod = rmod1
         ENDIF
c
c--Reduce time step (no synchronization needed)
c
         IF (rmod.LT.0.97) THEN
            stepi = dt*isteps(i)/imaxstep*rmod
            ibin = INT(LOG10(dt/stepi)/xlog2) + 1

            isteps(i) = imaxstep/(2**ibin)

            IF (isteps(i).LT.2 .AND. rmod.EQ.rmod1 .AND. 
     &           errdivtol.GT.0.999*errh/tolh .AND.
     &           errdivtol.LT.1.001*errh/tolh) THEN
               factor = 2.0*(dt/imaxstep)/stepi
C$OMP CRITICAL (writeiprint)
               WRITE(iprint,99200)
C$OMP END CRITICAL (writeiprint)
99200 FORMAT('**** R-K TOLERANCE TEMPORALLY INCREASED FOR dh ****')
               IF (factor.GT.4.0) THEN
C$OMP CRITICAL (writeiprint)
                  WRITE (iprint,*) ' Tolerance increase too high'
C$OMP END CRITICAL (writeiprint)
                  GOTO 7777
               ENDIF
C$OMP CRITICAL (writeiprint)
               WRITE(iprint,*)' Increased by:',factor*factor
               WRITE(iprint,*) '  part. i x,y,z',i,iorig(i),
     &              xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)       
               WRITE(iprint,*) '    vx,vy,vz',vxyzu(1,i),vxyzu(2,i),
     &              vxyzu(3,i)
               WRITE(iprint,*) '    h,rho',xyzmh(5,i),rho(i)
               WRITE(iprint,*) '    rm1,rmc,rmc2',rmod1,rmodcr,rmodcr2
               WRITE(iprint,*) '    neigh ',nneigh(i)
               WRITE(iprint,*) '    divv ',divv(i)
               WRITE(iprint,*) '    f1h,f2h ',f1ha(1,i),f2ha(1,i)
               WRITE(iprint,*) '    tryst,omin ',isteps(i),istminold
C$OMP END CRITICAL (writeiprint)
               isteps(i) = 2
               GOTO 7777
            ENDIF
         ENDIF
c
c--Increase time step, if synchronization allows
c
         IF (rmod.GE.2.0) THEN
            istep2 = 2*isteps(i)
            irat = MOD(idtsyn, istep2)
            IF (irat.EQ.0) THEN
               isteps(i) = istep2
            ENDIF
         ENDIF
 7777    IF (isteps(i).GT.imaxstep) isteps(i) = imaxstep
c
c--Report if time steps are too small
c
         IF (isteps(i).LT.2) THEN
C$OMP CRITICAL (writeiprint)
            WRITE (iprint, 99300)
99300       FORMAT ('STEP : Step too small ! Nothing can help!')
            WRITE (iprint,*)i, iorig(i), iphase(i), xyzmh(1,i), 
     &           xyzmh(2,i), xyzmh(3,i)
            WRITE (iprint,*)vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           vxyzu(4,i), xyzmh(5,i)
            WRITE (iprint,*)f1vxyzu(1,i),f1vxyzu(2,i),f1vxyzu(3,i)
            WRITE (iprint,*)f2vxyzu(1,i),f2vxyzu(2,i),f2vxyzu(3,i)
            WRITE (iprint,*)rmod1,rmodcr,rmodcr2
            WRITE (iprint,*)nneigh(i),rho(i)
            WRITE (iprint,*)errx,erry,errz
            WRITE (iprint,*)errvx,errvy,errvz,erru,errh
            IF (imhd.EQ.idim) THEN
               WRITE (iprint,*)errBx,errBy,errBz
            ENDIF
c            WRITE (iprint,*)gravx(i),gravy(i),gravz(i)
c            WRITE (iprint,*)gradpx(i),gradpy(i),gradpz(i)
c            WRITE (iprint,*)artvix(i),artviy(i),artviz(i)
            CALL quit
C$OMP END CRITICAL (writeiprint)
         ENDIF
c
c--Update minimum time step
c
         istepmax = MAX(istepmax, isteps(i))
         istepmin = MIN(istepmin, isteps(i))
c
c--Reset variables
c
         f1vxyzu(1,i) = f2vxyzu(1,i)
         f1vxyzu(2,i) = f2vxyzu(2,i)
         f1vxyzu(3,i) = f2vxyzu(3,i)
         f1vxyzu(4,i) = f2vxyzu(4,i)
         f1ha(1,i) = f2ha(1,i)
         IF (ifsvi.EQ.6) f1ha(2,i) = f2ha(2,i)

c         gravx1(i) = gravx(i)
c         gravy1(i) = gravy(i)
c         gradpx1(i) = gradpx(i)
c         gradpy1(i) = gradpy(i)
c         artvix1(i) = artvix(i)
c         artviy1(i) = artviy(i)

         it1(i) = it0(i) + isteps(i)/2
         it2(i) = it0(i) + isteps(i)
 855  CONTINUE
C$OMP END DO
C$OMP END PARALLEL
      IF (itiming) THEN
         CALL getused(ts92)
         ts9 = ts9 + (ts92 - ts91)
      ENDIF
c      DO j = 1, nlst0
c         i = llist(j)
c         IF (istepmin.EQ.isteps(i)) ikeepstepmin = i
c      END DO
c      WRITE (iprint,*) 'Px ',istepmin,iorig(i),iphase(i)
c      CALL FLUSH (iprint)
c
c--Make point mass timesteps equal to the minimum time step used
c
      DO i = 1, nptmass
         iptcur = listpm(i)
         IF (istepmin.GT.isteps(iptcur)) THEN
            istep2 = 2*isteps(iptcur)
            irat = MOD(idtsyn, istep2)
            IF (irat.EQ.0) THEN
               isteps(iptcur) = istep2
               istepmin = isteps(iptcur)
               it1(iptcur) = it0(iptcur) + isteps(iptcur)/2
               it2(iptcur) = it0(iptcur) + isteps(iptcur)
            ENDIF
         ELSEIF (istepmin.LT.isteps(iptcur)) THEN
            isteps(iptcur) = istepmin  
            it1(iptcur) = it0(iptcur) + isteps(iptcur)/2        
            it2(iptcur) = it0(iptcur) + isteps(iptcur)
         ENDIF
      END DO
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint,252) dt*itime/imaxstep, 
     &     nlst0, itime
 252  FORMAT (' step time ', 1F12.5,', particles moved ', I6,
     &     ' int time ', I10)
c
c--Update time
c
      IF (itime.GE.iteighth) THEN
         WRITE(iprint,253) dt*itime/imaxstep + gt,
     &        dt*itime/imaxstep, nlst0
 253     FORMAT ('Dynamic time = ', 1PE17.10,
     &        ' Step time = ', 1PE12.5,' Moved ', I10)

c         IF (nptmass.EQ.0) GOTO 700

         GOTO 700

         iptdump = iptdump + 1
         WRITE (itemchar, 3333) iptdump
 3333    FORMAT (I2)
         ptdebug = 'PTD' // itemchar

         dytime = dt*itime/imaxstep + gt

         IF (dytime.LT.50.0) THEN
            WRITE(iprint,*)' Ptmass debug dump'
            
            OPEN (15, FILE = ptdebug)
         
            tcomp = SQRT((3*pi)/(32*rhozero))
            WRITE(15,*) '#', dytime, dytime/tcomp
         
            DO i = 1,npart
               IF (iphase(i).NE.-1) THEN
                  iipt = listpm(1)
                  rx = xyzmh(1,i) - xyzmh(1,iipt)
                  ry = xyzmh(2,i) - xyzmh(2,iipt)
                  rz = xyzmh(3,i) - xyzmh(3,iipt)
                  r2 = rx*rx + ry*ry + rz*rz
                  r = SQRT(r2)
                  IF (r.EQ.0.0) THEN
                     vr = 0.0
                     vth = 0.0
                  ELSE
                     vr = (vxyzu(1,i)*rx+vxyzu(2,i)*ry+vxyzu(3,i)*rz)/r
                     vv = vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2
                     vth = vv - vr*vr
                     IF (vth.LT.0.) THEN
                        vth = 0.
                     ELSE
                        vth = SQRT(vv - vr*vr)
                     ENDIF
                  ENDIF

                  IF (r.LT.hacc*5.0) WRITE(15,333) r,dumxyzmh(1,i),
     &                 dumxyzmh(2,i),dumxyzmh(3,i),rho(i),dumxyzmh(5,i),
     &                 vr,vth,dumvxyzu(4,i),divv(i),
     &                 nneigh(i),i,iorig(i),iphase(i)
 333              FORMAT(17(1PE13.6,1X),I6,1X,I6,1X,I6,1X,I6)
               ENDIF
            END DO

            CLOSE (15)
         ENDIF

 700     CONTINUE
         CALL FLUSH (iprint)
         CALL FLUSH (iptprint)
         iteighth = itime + imaxstep/8
      ENDIF
c
c--Accrete particles near point mass, or create point mass
c
      IF (nlst0.GT.nptmass .AND. (nptmass.NE.0 .OR. icreate.EQ.1)) THEN
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

         CALL accrete(dt, realtime, isave)

         IF (itiming) THEN
            CALL getused(taccrete2)
            taccrete = taccrete + (taccrete2 - taccrete1)
         ENDIF
      ENDIF
c
c--Update it1 and it2 for bins
c
c      WRITE (*,*) 'Zeroing '
c      WRITE (iprint,*) 'Zeroing '
      DO i = 1, nbinmax
         IF (it2bin(i).EQ.itime) THEN
c            WRITE (*,*) ' t-1 ',i,it2bin(i),itime
c            WRITE (iprint,*) ' t-1 ',i,it2bin(i),itime
            nlstbins(i) = 0
            it1bin(i) = itime + 2**i/2
            it2bin(i) = itime + 2**i
         ELSEIF (it2bin(i).LT.itime) THEN
c            WRITE (*,*) ' t-2 ',i,it2bin(i),itime
c            WRITE (iprint,*) ' t-2 ',i,it2bin(i),itime
            it1bin(i) = itime + 2**i/2
            it2bin(i) = itime + 2**i
            IF (nlstbins(i).NE.0) THEN
               WRITE (*,*) 'ERROR - nlstbins(i).NE.0'
               WRITE (iprint,*) 'ERROR - nlstbins(i).NE.0'
               CALL quit
            ENDIF
         ENDIF
      END DO
c
c--Update the lists of particles in each bin
c
c      WRITE (*,*) 'Updating ',nlst0
c      WRITE (iprint,*) 'Updating ',nlst0
      DO j = 1, nlst0
         i = llist(j)
         IF (i.GT.idim) THEN
            WRITE (*,*) 'ERROR - i.GT.idim 4'
            WRITE (iprint,*) 'ERROR - i.GT.idim 4'
            CALL quit
         ELSEIF (i.LE.0) THEN
            WRITE (*,*) 'ERROR - i.LE.0 4'
            WRITE (iprint,*) 'ERROR - i.LE.0 4'
            WRITE (*,*) nlst0,j,itime
            WRITE (iprint,*) nlst0,j,itime
            DO i2 = 1, j
               WRITE (iprint,*) i2, llist(i2)
            END DO
            CALL quit
         ENDIF
c         WRITE (*,*) '   part ',j,i,iphase(i)
c         WRITE (iprint,*) '   part ',j,i,iphase(i)
         IF (iphase(i).NE.-1) THEN
            ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
            IF (ibin.GT.nbinmax) THEN
               WRITE (*,*) 'ERROR - ibin.GT.nbinmax 4'
               WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax 4'
               CALL quit
            ENDIF
            nlstbins(ibin) = nlstbins(ibin) + 1
            IF (nlstbins(ibin).GT.idim) THEN
               WRITE (*,*) 'ERROR - nlstbins(ibin).GT.idim 4'
               WRITE (iprint,*) 'ERROR - nlstbins(ibin).GT.idim 4'
               CALL quit
            ENDIF
            listbins(nlstbins(ibin),ibin) = i
c            WRITE (*,*) '   step ',ibin,isteps(i),nlstbins(ibin)
c            WRITE (iprint,*) '   step ',ibin,isteps(i),nlstbins(ibin)
c            IF (i.EQ.1) THEN
c               WRITE (*,*) 'for i=1 upd, ',iphase(i),isteps(i),ibin,
c     &              nlstbins(ibin)
c              WRITE (iprint,*) 'for i=1 upd, ',iphase(i),isteps(i),
c     &              ibin,nlstbins(ibin)
c               WRITE (*,*) nlst0,j,it0(i),it1(i),it2(i),it1bin(ibin),
c     &              it2bin(ibin)
c          WRITE (iprint,*) nlst0,j,it0(i),it1(i),it2(i),it1bin(ibin),
c     &              it2bin(ibin)
c            ENDIF
            IF (it1bin(ibin).NE.it1(i)) THEN
               WRITE (*,*) 'ERROR - it1bin 2',it1bin(ibin),it1(i),i,
     &              it0(i),itime,isteps(i),it2(i)
               CALL quit
            ENDIF
            IF (it2bin(ibin).NE.it2(i)) THEN
               WRITE (*,*) 'ERROR - it2bin 2',it2bin(ibin),it2(i),i,
     &              it0(i),itime,isteps(i),it1(i)
               CALL quit
            ENDIF
         ENDIF
      END DO
c
c--Create NEW particles to keep number of particles within a shell
c     constant.  
c
      IF (nptmass.NE.0 .OR. icreate.EQ.1) THEN
         IF (ibound/10.EQ.9 .AND. nshell.GT.inshell) THEN
            iaccr = 0
            ikilled = 0
            nlstacc = 0
            nneightotsave = nneightot

            nlst0 = nshell-inshell
            nlst = nlst0
            WRITE (iprint,*) 'Add ',nlst0
            DO i = 1, nlst0
               nnew = npart + i
               isort(nnew) = nnew
               iorig(nnew) = nnew
               llist(i) = nnew
               CALL phoenix(nnew, idtsyn, itime)
               nactive = nactive + 1
            END DO
            npart = npart + nlst0
            n1 = n1 + nlst0
            IF (npart.GT.idim) THEN
               CALL error(where,3)
            ENDIF

            DO j = 1, nlst0
               iscurrent(llist(j)) = .TRUE.
            END DO

            nghost = 0
            CALL ghostp3(npart,xyzmh,vxyzu,ekcle,Bevolxyz)
            ntot = npart + nghost

            DO j = npart - nlst0 + 1, npart
               IF (iphase(j).GE.0) THEN
                  DO k = 1, 5
                     dumxyzmh(k,j) = xyzmh(k,j)
                  END DO
                  DO k = 1, 4
                     dumvxyzu(k,j) = vxyzu(k,j)
                  END DO
                  IF (ifsvi.EQ.6) dumalpha(j) = alphaMM(j)
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
                     DO l = 1, 5
                        dumekcle(l,j) = ekcle(l,j)
                     END DO
                  ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
                  IF (ifsvi.EQ.6) dumalpha(j) = MIN(alphamax,
     &                 alphaMM(k) + deltat*f1ha(2,k))
               ENDIF
            END DO

            IF (igrape.EQ.0) THEN
               CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)
            ENDIF

            neighmean = (neimax + neimin)/2

            iokay = 1
            DO j = 1, nlst0
               i = llist(j)

               ivalue = 0
               ichkloop = 0
 2000          ichkloop = ichkloop + 1

               IF (igrape.EQ.0) THEN
                  CALL insulate(3,ntot,npart,dumxyzmh,f1vxyzu)
                  numneigh = nneigh(i)
               ELSE
                  CALL getneigh(i,ntot,xyzmh(5,i),dumxyzmh,nlist,
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
            CALL derivi (dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,
     &      f1ha,npart,ntot,ireal,dumalpha,dumekcle,dumBevolxyz,f1Bxyz)

            time = dt*itime/imaxstep + gt
            DO j = 1, nlst0
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
                  CALL quit
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               IF (nlstbins(ibin).GT.idim) THEN
                  WRITE (*,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  WRITE (iprint,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  CALL quit
               ENDIF
               listbins(nlstbins(ibin),ibin) = i
            
               IF (it1bin(ibin).NE.it1(i)) THEN
                  WRITE (*,*) 'ERROR - it1bin 5',it1bin(ibin),it1(i),i,
     &                 it0(i),itime,isteps(i),it2(i)
                  CALL quit
               ENDIF
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin 5',it2bin(ibin),it2(i),i,
     &                 it0(i),itime,isteps(i),it1(i)
                  CALL quit
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
                           CALL quit
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
               CALL quit
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
                  IF (ifsvi.EQ.6) dumalpha(j) = alphaMM(j)
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
                     DO l = 1, 5
                        dumekcle(l,j) = ekcle(l,j)
                     END DO
                  ENDIF
            IF (iener.EQ.2.AND.dumvxyzu(4,j).LT.0.0) dumvxyzu(4,j)=0.15
                  IF (ifsvi.EQ.6) dumalpha(j) = MIN(alphamax,
     &                 alphaMM(k) + deltat*f1ha(2,k))
                  IF (imhd.EQ.idim) THEN
                     DO l = 1, 3
                        dumBevolxyz(l,j) = Bevolxyz(l,j) 
     &                                   + deltat*f1Bxyz(l,j)
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
                  CALL getneigh(i,ntot,xyzmh(5,i),dumxyzmh,nlist,
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
     &      f1ha,npart,ntot,ireal,dumalpha,dumekcle,dumBevolxyz,f1Bxyz)
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
                  CALL quit
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               IF (nlstbins(ibin).GT.idim) THEN
                  WRITE (*,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  WRITE (iprint,*) 'ERROR - nlstbins(ibin).GT.idim 5'
                  CALL quit
               ENDIF
               listbins(nlstbins(ibin),ibin) = i

               IF (it1bin(ibin).NE.it1(i)) THEN
                  WRITE (*,*) 'ERROR - it1bin 5',it1bin(ibin),it1(i),i,
     &                 it0(i),itime,isteps(i),it2(i)
                  CALL quit
               ENDIF
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin 5',it2bin(ibin),it2(i),i,
     &                 it0(i),itime,isteps(i),it1(i)
                  CALL quit
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
                           CALL quit
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
               xyzmh(5,ipt1) = accfac*semiaxis*
     &              (0.38 - 0.20*LOG10(qratio))
            ELSE
               WRITE (iprint,*) 'ERROR: qratio < 0.05'
               CALL quit
            ENDIF
            IF (qratio.LT.0.5) THEN
               xyzmh(5,ipt2) = accfac*semiaxis*
     &              (0.462*(qratio/qratio1)**(1.0/3.0))
            ELSE
               xyzmh(5,ipt2) = accfac*semiaxis*
     &              (0.38 + 0.20*LOG10(qratio))
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
c
c--Return or loop again
c
      IF (itime.GE.itnext) THEN

c         WRITE (*,*) 'itime.GE.itnext'
c         WRITE (iprint,*) 'itime.GE.itnext'

         DO i = 1, nbinmax
            nlstbins(i) = 0
            it1bin(i) = 2**i/2
            it2bin(i) = 2**i
         END DO

         DO i = 1, npart
            IF (iphase(i).NE.-1) THEN
               it0(i) = 0
               it1(i) = isteps(i)/2
               it2(i) = isteps(i)

               ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
               IF (ibin.GT.nbinmax) THEN
                  WRITE (*,*) 'ERROR - ibin.GT.nbinmax 5'
                  WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax 5'
                  CALL quit
               ENDIF
               nlstbins(ibin) = nlstbins(ibin) + 1
               listbins(nlstbins(ibin),ibin) = i
               IF (it1bin(ibin).NE.it1(i)) THEN
                  WRITE (*,*) 'ERROR - it1bin 5'
                  CALL quit
               ENDIF
               IF (it2bin(ibin).NE.it2(i)) THEN
                  WRITE (*,*) 'ERROR - it2bin 5'
                  CALL quit
               ENDIF
            ENDIF
         END DO

         IF (igrape.EQ.0) THEN
            nneightot = 0
            DO i = 1, nlstnneigh
               ipart = llist(i)
               nneightot = nneightot + nneigh(ipart)
            END DO
         ENDIF
         
         RETURN
      ENDIF

      ilocal = ilocal + 1
      IF (MOD(ilocal,50).EQ.0) THEN
         ilocal = 0
         CALL secmes
      ENDIF

      GOTO  100

      END

