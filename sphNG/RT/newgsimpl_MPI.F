      SUBROUTINE GSIMPLS(dtdo,dtmax,itime,npart,ntot,ekcle,xyzmh,vxyzu,
     &     dedxyz,rho,alphaMMpass,moresweep,nit,errorE,errorU,
     &     istepmindone)
 
#ifdef MPIALL
#include "mpi_sup.h"
#endif
 
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      DIMENSION ekcle(5,iradtrans2)
      DIMENSION xyzmh(5,mmax2)
      DIMENSION vxyzu(4,idim2)
      DIMENSION dedxyz(3,iradtrans2)
      REAL*4 rho(idim2),alphaMMpass(isizealphaMM,idim2)

      EXTERNAL boundaryparticle
      LOGICAL boundaryparticle

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/vsmooth'
c      INCLUDE 'COMMONS/rad2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/radsink'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/radstore'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/boundheight'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/andrea'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/radtran2'
      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/fullstep'
      INCLUDE 'COMMONS/rtstat'
      INCLUDE 'COMMONS/lyonfiles'

#ifdef MPICOPY
      INCLUDE 'COMMONS/mpicopy'
#endif

      PARAMETER (nswmax = 200)
      PARAMETER (nswmax2 = 100)

      COMMON /getkap/ iflag

      COMMON /rtvarr6/  EU0(2,iradtrans2)
      REAL U1i,E1i

      REAL*4 dvdx
      REAL*4 rhoreal4
      COMMON /rtvarr5/ dvdx(9,iradtrans2)

      INTEGER nosweep
      LOGICAL moresweep
      LOGICAL moresweep2
      REAL dx,dy,dz,dvz,dvx,dvy,lightspeed

      COMMON /rtvarr1/ transfer(iradtrans2),
     &     addenergyparticle(iradtrans2)
c      COMMON /rtvarr4/ varijnew(icompactmax)
      COMMON /rtvarr4/ varijnew(iradtrans2)

c      DIMENSION vari(2,iradtrans)
c      REAL*4 varij(4,icompactmax),varij2(3,icompactmax)
c      DIMENSION varinew(2,iradtrans)
cc      DIMENSION varijnew(icompactmax)
c      DIMENSION varijnew(iradtrans)
c      DIMENSION transfer(iradtrans)

      COMMON /rtvari/ neighblock(iradtrans)
c      DIMENSION neighblock(iradtrans), addenergyparticle(iradtrans)
      COMMON /rtvarr7/ fractionexposed(nmaxexposed)
c      DIMENSION fractionexposed(nmaxexposed)

      REAL*8 ptmasstime,dtstep,ptmassmodel,ptmassradius,ptmasstemp,
     &     accrate

#ifdef RTBAR
c Needed for LYON stellar evolution code only
      REAL*8 starradius,starluminosity,startime,starmass,evolveto,
     &lyonaccrate
      INTEGER nflag,starid
      COMMON /cepsilon/epsilon
#endif

      REAL maxerrE2,maxerrE2last,maxerrU2,maxerrU2last,maxerrU2old
#ifdef MPIALL
      REAL maxerrE2tot, maxerrU2tot
      LOGICAL moresweeptot
#endif

      LOGICAL ifirst
      DATA ifirst/.TRUE./

#ifdef MPIALL
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL8, i2REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i2REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i5REAL8, ierr)
#endif

c
c--Set up constants in Code units
c
      lightspeed = c / udist * utime
      uradconst = radconst / uergcc
      dtimax = dtdo/imaxstep

c      print *,'entry gsimpl ',rho(2550),rho(30788),rho(30789)
c     &  ,npart,ekcle(1,1),ekcle(2,1),xyzmh(5,1),vxyzu(1,1),
c     &     isteps(1)


c      boundtemp=12.037
c
c--Set errors to zero for iteration start
c
      numoscillations1 = 0
      numoscillations2 = 0
      numequal = 0
      numcomp = 0
      ipos = 1
      iphase1 = 0
c      print *,' ENTERED GSIMPLS ',nptmass, nptmasslast, itime

c
c--Sink particles radiate energy
c
      IF (ifirst) THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(iexposedold)
C$OMP& private(i)
         DO i = 1, iradtrans
            iexposedold(i) = .FALSE.
         END DO
C$OMP END PARALLEL DO
         ifirst = .FALSE.
      ENDIF
c
c--Find shortest timestep and initialise quantities to zero
c
      IF (nptmass.GE.1) THEN
c         IF (nptmasslast.EQ.0) CALL EGGSETUP(0)
      ENDIF

C$OMP PARALLEL default(none)
C$OMP& shared(xmaxerr1old,xmaxerr2old,ncompact,ivar)
C$OMP& shared(stellarenergy,addenergyparticle,transfer,varinew)
C$OMP& private(i,n,k)
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(runtime)
      DO i = 1, ntests
         xmaxerr1old(i) = 0.0
         xmaxerr2old(i) = 0.0
      END DO
C$OMP END DO
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)
         stellarenergy(i) = 0.
         addenergyparticle(i) = 0.

         transfer(i) = 0.

         DO k = 1, 2
            varinew(k,i) = 0.
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL

c      print *,iproc,' IN RT a'

c      rrmin = 1.0e+30
c      DO n = 1, ncompactlocal
c         i = ivar(3,n)
c         
c         r2 = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 +xyzmh(3,i)**2)
c         IF (r2.LT.rrmin) THEN
c            rrmin = r2
c            ikeepit = i
c            nkeepit = n
c         ENDIF
c      END DO

c      print *,iproc,': FOUND ',nkeepit,ivar(3,nkeepit),
c     &     iunique(iorig(ivar(3,nkeepit))),xyzmh(1,ivar(3,nkeepit)),
c     &     xyzmh(2,ivar(3,nkeepit)),xyzmh(3,ivar(3,nkeepit))

c      print *,'pre-stellar-feedback ',ncompact,ncompactlocal

      IF (istellarfeedback.EQ.0) GOTO 888
c
c--If new sink particle has been created, need to initialize object
c
      IF (nptmasslast.LT.nptmass) THEN
         DO iptn = nptmasslast + 1, nptmass

            nexposedold(iptn) = 0
            cummulativeenergy(iptn) = 0.
            actualcumenergy(iptn) = 0.
c
c--Stellar evolution model is initialised at time = 0
c     ptmassmodel is given in solar masses
c
            ptmasstime = 0.0
            ptmassform(iptn) = gt + (dtmax*itime)/imaxstep
            dtstep = dtimax*istepmindone*utime/31557600.0

            print *,'Starting ptmass ',iptn,ptmassform(iptn),dtstep,
     &           istepmindone

            IF (istellarfeedback.EQ.1) THEN
               print *,' '
               print *,' '
               print *,'CALLING LYONINIT ',ptmasstime,dtstep
               print *,' '
               print *,' '

c               CALL EGGINIT(iptn, ptmasstime, dtstep, ptmassmodel, 
c     &            ptmassradius, ptmasstemp)
               print *,'EGGINIT Values ',iptn,ptmasstime,
     &              ptmassform(iptn),dtstep,ptmassmodel,
     &              ptmassradius,ptmasstemp

            ELSEIF (istellarfeedback.EQ.2) THEN
#ifdef RTBAR
               nflag = 1
               accrate = 1.0E-4
               ptmassaccrate(iptn) = accrate

               evolveto = 10.0

               CALL evollyon(nflag,iptn,evolveto,accrate,starradius,
     &              starluminosity,startime,starmass)
               stellartime(iptn) = startime
               stellarmass(iptn) = starmass
               stellarradius(iptn) = starradius/solarr
               stellarluminosity(iptn) = starluminosity/solarl

               ptmassinitial(iptn)=-ABS(xyzmh(4,listpm(iptn))-starmass)
               ptmassinitial(iptn) = 0.0

               WRITE (53,*) 'First ',xyzmh(4,listpm(iptn)),starmass,
     &              ptmassinitial(iptn)
               CALL FLUSH(53)

               print *,'EvolLyon Values ',iptn,ptmasstime,
     &              ptmassform(iptn),accrate,
     &              starradius,starluminosity,startime,starmass
#else
               print *,'istellarfeedback = 2 but evollyon not compiled!'
               CALL quit(0)
#endif
            ELSEIF (istellarfeedback.EQ.3) THEN
               print *,'New star noted ',iptn,ptmasstime,
     &              ptmassform(iptn)
            ELSEIF (istellarfeedback.EQ.4) THEN
               ptmassinitial(iptn) = xyzmh(4,listpm(iptn))
               ptmassrecord(1,iptn) = ptmassinitial(iptn)
               ptaccrecord(1,iptn) = ptmassform(iptn)
               iptmassrecord(iptn) = 1
c
c--Note that for istellarfeedback=4, the stellar mass is precisely this - the
c     estimated mass of the stellar core, excluding the estimated mass of
c     the disc within the sink particle accretion radius
c     (i.e. stellarmass << xyzmh(4,listpm(iptn)) )
c
               stellarmass(iptn) = 0.006*solarm/umass

               print *,'New star noted ',iptn,ptmasstime,
     &              ptmassform(iptn),ptmassinitial(iptn)

            ELSEIF (istellarfeedback.EQ.5) THEN

               PRINT *, 'NEW STAR NOTED',
     &         iptn,nptmasstot,ptmassform(iptn)
               PRINT *, ' '
               PRINT *, 'DISK ACCRETION'
               PRINT *, ' '

               diskinit(iptn) = 0.0

               stellartime(iptn) = 0.0
               stellarmass(iptn) = 0.0 
               stellarradius(iptn) = 0.0
               stellarluminosity(iptn) = 0.0

            ELSEIF (istellarfeedback.EQ.6) THEN

               PRINT *, 'NEW STAR NOTED',
     &         iptn,nptmasstot,ptmassform(iptn)
               PRINT *, ' '
               PRINT *, 'DISK ACCRETION + STELLAR EVOLUTION'
               PRINT *, ' '

               diskinit(iptn) = 0.0

               stellartime(iptn) = 0.0
               stellarmass(iptn) = 0.0
               stellarradius(ipnt) = 0.0
               stellarluminosity(iptn) = 0.0

               lyontime(iptn) = 0.0
               lyonmass(iptn) = 0.0
               lyonradius(ipnt) = 0.0
               lyonluminosity(iptn) = 0.0
               listlyon(iptn) = nptmasstot

               file_stellarmodel_in(iptn) = file_stellarmodel_new
               file_stellarmodel_out(iptn) = file_stellarmodel_new

            ELSE
               WRITE (iprint,*) 'ERROR - istellarfeedback ',
     &              istellarfeedback
               CALL quit(1)
            ENDIF
         END DO
      ENDIF
      nptmasslast = nptmass
c
c--Set stellar luminosities
c
      DO iptn = 1, nptmass
         ptmasstime = gt + (dtmax*itime)/imaxstep - ptmassform(iptn)
         ptmasstime = ptmasstime*utime/31557600.0
         dtstep = dtimax*istepmindone*utime/31557600.0

         IF (istellarfeedback.EQ.1) THEN
c
c--Stellar evolution model begins with an 0.01 M_sol (10 M_Jup) object
c
            ptmasstime = gt + (dtmax*itime)/imaxstep - ptmassform(iptn)
            ptmasstime = ptmasstime*utime/31557600.0
         ptmassmodel=MAX(0.01*solarm,xyzmh(4,listpm(iptn))*umass)/solarm

c            ptmassmodel=0.01
            ptmassmodel= 0.01 + ptmasstime*1.0e-5/0.1/2.0*ptmasstime

            dtstep = dtimax*istepmindone*utime/31557600.0

            print *,' '
            print *,' '
            print *,'CALLING EGGLOOKUP ',iptn,ptmasstime,ptmassmodel,
     &           dtstep
            print *,' '
            print *,' '

c            CALL EGGLOOKUP(iptn,ptmasstime,ptmassmodel,dtstep,
c     &         ptmassradius, ptmasstemp)
c          ptmassluminosity(iptn) = (4.0*pi*((ptmassradius*solarr)**2.0)
c     &           *(stepboltz)*(ptmasstemp**4.0))/(uergg*umass/utime)

c            accluminosity = (((gg*ptmassmodel*solarm*1.0E-5*solarm)/
c     &        31557600.0)/(ptmassradius*solarr))/(uergg*umass/utime)

c        print *,'EGGLOOKUP Values ',iptn,ptmasstime,ptmassmodel,dtstep,
c     &        ptmassradius,ptmasstemp,ptmassluminosity(iptn),
c     &        accluminosity,xyzmh(4,listpm(iptn))

c            write (60,99333) iptn,ptmasstime,ptmassmodel,dtstep,
c     &        ptmassradius,ptmasstemp,ptmassluminosity(iptn),
c     &        accluminosity,xyzmh(4,listpm(iptn)),
c     &        gt + (dtmax*itime)/imaxstep
c            CALL FLUSH(60)
c99333        FORMAT(I4,8(1X,1PE16.9),1X,1PE16.9)

c            print *,' '
c            print *,' '

c         ptmassluminosity(iptn) =ptmassluminosity(iptn) + accluminosity

c         ptmassluminosity(iptn) = 1.0E+04*3.827E+33/(uergg*umass/utime)

            ptmassluminosity(iptn) = 0.0

         ELSEIF (istellarfeedback.EQ.2) THEN
#ifdef RTBAR
            iacclimit = 2
            IF (iacclimit.EQ.1) THEN
c
c--Time starts from when sink particle exceeds start mass
c
               actualmass = xyzmh(4,listpm(iptn))*umass/solarm
               ptmassmodelmass = stellarmass(iptn)

               diffmass = MAX(0.0,actualmass - ptmassmodelmass - 
     &              ptmassinitial(iptn))

               print *,'XXXX ',diffmass,actualmass,ptmassmodelmass,
     &              ptmassinitial(iptn)
c
c--Try to fix 90% of mass deficit during this timestep (one timescale 
c     to make sure stellar
c     evolution code does not break, another for the accretion rate)
c
               accrate = MAX(0.9*diffmass/dtstep,1.0E-8)
               accrate = MIN(accrate,1.0E-5)

               print *,'YYYY ',accrate,ifullstep(listpm(iptn))
c
c--Forgets about 90% of the mass deficit over some timescale
c
               timescalefix = 1.0

               IF (ifullstep(listpm(iptn))) THEN
                  ptmassinitial(iptn) = ptmassinitial(iptn) + 
     &              MAX(0.0,0.9*(diffmass-accrate*dtstep)/timescalefix)
               ENDIF
c
c--End section
c
            
            ELSEIF (iacclimit.EQ.2) THEN
               IF (0.99*xyzmh(4,listpm(iptn))*umass/solarm
     &              - stellarmass(iptn).GT.
     &              ptmassinitial(iptn)) THEN
                  ptmassinitial(iptn) = 0.99*xyzmh(4,listpm(iptn))*
     &                 umass/solarm - stellarmass(iptn)
               ENDIF

               diffmass = MAX(0.0,xyzmh(4,listpm(iptn)) - 
     &              stellarmass(iptn) -
     &              ptmassinitial(iptn))

               accrate = MAX(1.0*diffmass/dtstep,1.0E-8)
               accrate = MIN(accrate,1.0E-4)
c
c--Limit change in accrate
c
               IF (accrate/ptmassaccrate(iptn).GT.1.1) THEN
                  accrate = 1.5*ptmassaccrate(iptn)
               ELSEIF (accrate/ptmassaccrate(iptn).LT.0.9) THEN
                  accrate = 0.5*ptmassaccrate(iptn)
               ENDIF

            ELSE
c
c--Allow initial sink particle's mass to be that AFTER it has accreted all
c     particles within haccall
c
            IF (ptmassinitial(iptn).LE.0.0) THEN
               IF (icall.EQ.3) THEN
                  currentmass = xyzmh(4,listpm(iptn)) - 
     &                 stellarmass(iptn)

               WRITE (53,*) 'Next ',xyzmh(4,listpm(iptn)),
     &                 stellarmass(iptn),
     &                 ptmassinitial(iptn),
     &                 ABS((currentmass - ABS(ptmassinitial(iptn)))/
     &                 ptmassinitial(iptn)),
     &                 currentmass
               CALL FLUSH(53)

                  IF (ABS((currentmass - ABS(ptmassinitial(iptn)))/
     &                 ptmassinitial(iptn)).GT.0.05 .OR.
     &                 currentmass.LE.0.0) THEN
                     ptmassinitial(iptn) = - ABS(currentmass)
                  ELSE
                     ptmassinitial(iptn) = currentmass
                  ENDIF
               ENDIF
               ptmassinitiallocal = MAX(0.0,currentmass)
            ELSE
               ptmassinitiallocal = ptmassinitial(iptn)
            ENDIF

c
c--Time starts from when sink particle exceeds start mass
c
               actualmass = xyzmh(4,listpm(iptn))*umass/solarm
               ptmassmodelmass = stellarmass(iptn)

               diffmass = MAX(0.0,actualmass - ptmassmodelmass - 
     &              ptmassinitial(iptn))

               print *,'XXXX ',diffmass,actualmass,ptmassmodelmass,
     &              ptmassinitial(iptn)
c
c--Try to fix 90% of mass deficit during this timestep (one timescale 
c     to make sure stellar
c     evolution code does not break, another for the accretion rate)
c
               accrate = MAX(1.0*diffmass/dtstep,1.0E-8)
               accrate = MIN(accrate,1.0E-4)

               print *,'YYYY ',accrate,ifullstep(listpm(iptn))
c
c--Forgets about 90% of the mass deficit over some timescale
c
               timescalefix = 1.0

               IF (ifullstep(listpm(iptn))) THEN
c                  ptmassinitial(iptn) = ptmassinitial(iptn) + 
c     &                 0.9*(diffmass-accrate*dtstep)
c     &              MAX(0.0,0.7*(diffmass-accrate*dtstep)/timescalefix)
               ENDIF


            ENDIF
c
c--New time required for stellar model
c
            ptmasstime = ptmasstime + dtstep

c     &           (stellarmass(iptn)*
c     &           umass/solarm/accrate)/10.*0.0)

            IF (ifullstep(listpm(iptn)) .AND.
     &           ptmasstime.GT.stellartime(iptn)
c     &           + 0.016
c     &           *(1.0E-3/accrate)
     &           ) THEN
               evolveto = stellartime(iptn) + MAX(10.0,dtstep)

               nflag = 99
c               IF (dtstep.NE.0.0) THEN
c                  accrate = MAX(0.0,
c     &                 ((xyzmh(4,listpm(iptn))-ptmassinitiallocal)*
c     &                 umass/solarm - stellarmass(iptn))/dtstep)
c                  accrate = MIN(1.0E-3,accrate)
c               ELSE
c                  accrate = 1.0E-8
c               ENDIF
c               accrate = 7.0E-4

               print *,' '
               print *,' '
               print *,'CALLING EvolLyon ',iptn,ptmasstime,dtstep,
     &              stellartime(iptn),accrate
               print *,' '
               print *,' '

               CALL evollyon(nflag,iptn,evolveto,accrate,starradius,
     &              starluminosity,startime,starmass)

               stellartime(iptn) = startime
               stellarmass(iptn) = starmass
               stellarradius(iptn) = starradius/solarr
               stellarluminosity(iptn) = starluminosity/solarl

               WRITE (54,78779) iptn,ptmasstime,stellartime(iptn),
     &              accrate,ptmassaccrate(iptn)
78779          FORMAT(I3,1X,4(1PE12.5,1X))
               ptmassaccrate(iptn) = accrate

               ptmassluminosity(iptn) = starluminosity/
     &              (uergg*umass/utime)
               accluminosity = (((gg*starmass*solarm*accrate*solarm)/
     &              31557600.0)/(starradius))/(uergg*umass/utime)

               print *,'LyonEvol Values ',iptn,ptmasstime,dtstep,
     &              accrate,xyzmh(4,listpm(iptn)),starmass,
     &              starradius/solarr,starluminosity/solarl,
     &              startime,ptmassluminosity(iptn),
     &              accluminosity

               WRITE(52,78780) iptn,ptmasstime,dtstep,
     &              accrate,xyzmh(4,listpm(iptn)),starmass,
     &              starradius/solarr,starluminosity/solarl,
     &              startime,ptmassluminosity(iptn),
     &              accluminosity,ptmassinitial(iptn),
     &              ptmassinitial(iptn)+starmass
78780          FORMAT(I3,1X,12(1PE12.5))
               CALL FLUSH(52)
               CALL FLUSH(20)

             ptmassluminosity(iptn)=ptmassluminosity(iptn)+accluminosity
c            ptmassluminosity(iptn) = 0.0
            ELSE
               print *,' '
               print *,' '
               print *,'RE-USING EvolLyon ',iptn,ptmasstime,dtstep,
     &              stellartime(iptn)
               print *,' '
               print *,' '

               starluminosity = stellarluminosity(iptn)
               starmass = stellarmass(iptn)
               starradius = stellarradius(iptn)

               ptmassluminosity(iptn) = starluminosity*solarl/
     &              (uergg*umass/utime)
               accluminosity = (((gg*starmass*solarm*accrate*solarm)/
     &              31557600.0)/(starradius*solarr))/(uergg*umass/utime)

             ptmassluminosity(iptn)=ptmassluminosity(iptn)+accluminosity
c            ptmassluminosity(iptn) = 0.0
            ENDIF
#else
              print *,'istellarfeedback = 2 but evollyon not compiled!'
              CALL quit(0)
#endif
         ELSEIF (istellarfeedback.EQ.3) THEN
            ptmassluminosity(iptn)=1.0E+01*solarl/(uergg*umass/utime)
c            ptmassluminosity(iptn) = 0.0
         ELSEIF (istellarfeedback.EQ.4) THEN
            realtime = gt + (dtmax*itime)/imaxstep
c
c--Accretion rate since last recorded time in M_sun/yr
c
            accrate = (xyzmh(4,listpm(iptn)) - ptmassrecord(1,iptn))*
     &           umass/solarm/
     &           ((realtime - ptaccrecord(1,iptn))*utime/31557600.0)
c
c--Save time is minimum of 0.1 years and max of 100 years, but is set by
c     accretion rate (lower accretion rate gives longer period between
c     saves)
c
            savetime = MIN(100.,MAX(0.1,(xyzmh(4,listpm(iptn))*
     &           umass/solarm)/accrate/1000.0))

c
c--If sink particle is doing a full step
c
            IF (ifullstep(listpm(iptn))) THEN
c
c--If there is only one recorded value and the sink particle mass is still
c     changing rapidly, then over-write the 'first' record.
c
               IF (iptmassrecord(iptn).EQ.1 .AND. 
     &              ABS(xyzmh(4,listpm(iptn))/ptmassrecord(1,iptn))
     &              .GT.1.1) THEN
                  idonevalue = 0
                  ptmassrecord(1,iptn) = xyzmh(4,listpm(iptn))
                  ptaccrecord(1,iptn) = realtime
ccc                  stellarmass(iptn) = xyzmh(4,listpm(iptn))
c
c--Else see whether the current time is more than savetime from the last
c     saved point mass data
c
               ELSEIF ((realtime - ptaccrecord(iptmassrecord(iptn),
     &                 iptn))*utime/31557600.0.GT.savetime) THEN
c
c--Only iaccrec values are saved (e.g. 10 past values).  If there are less
c     values than this currently, simply add new record of the point mass's
c     mass and the time it had that mass.
c
                  IF (iptmassrecord(iptn).LT.iaccrec) THEN
                     idonevalue = 1
                     iptmassrecord(iptn) = iptmassrecord(iptn) + 1

                     ptmassrecord(iptmassrecord(iptn),iptn) = 
     &                    xyzmh(4,listpm(iptn))
                     ptaccrecord(iptmassrecord(iptn),iptn) = realtime
c
c--If already have iaccrec records, then delete the oldest record and store
c     the new one.
c
                  ELSE
                     idonevalue = 2
                     DO iloop = 1, iaccrec-1
                     ptmassrecord(iloop,iptn)=ptmassrecord(iloop+1,iptn)
                     ptaccrecord(iloop,iptn) = ptaccrecord(iloop+1,iptn)
                     END DO
                     ptmassrecord(iaccrec,iptn) = xyzmh(4,listpm(iptn))
                     ptaccrecord(iaccrec,iptn) = realtime
                  ENDIF
               ENDIF
            ENDIF
c
c--Set accretion rate based on the current time and mass and the oldest
c     recorded value of the point mass's mass.  Accretion rate is in
c     Msun/yr.
c
            IF (realtime - ptaccrecord(1,iptn).GT.0.0) THEN
               accrate = (xyzmh(4,listpm(iptn)) - ptmassrecord(1,iptn))*
     &              umass/solarm/
     &              ((realtime - ptaccrecord(1,iptn))*utime/31557600.0)
            ELSE
               accrate = 0.0
            ENDIF
c
c--Place limits on accrate
c
            IF (accrate.LT.0.0) accrate = 0.0
            IF (accrate.GT.1.0E-5) accrate = 1.0E-5
c
c--Ignore intrinsic protostellar luminosity
c
            ptmassluminosity(iptn) = 0.0
c
c--Set accretion luminosity using the total mass of the point mass and
c     a fixed radius.
c
            accluminosity = (((gg*stellarmass(iptn)*umass*
     &           accrate*solarm)/
     &           31557600.0)/(2.0*solarr))/(uergg*umass/utime)

            ptmassluminosity(iptn)=ptmassluminosity(iptn)+accluminosity
c            ptmassluminosity(iptn) = 0.0

c
c--Record the stellar mass as given based on the accretion rate
c     Output to a file.
c
            IF (ifullstep(listpm(iptn))) THEN
               stellarmass(iptn) = stellarmass(iptn) + accrate*dtstep

               WRITE(52,78789) iptn,ptmasstime,dtstep,
     &              accrate,xyzmh(4,listpm(iptn)),stellarmass(iptn),
     &              ptmassluminosity(iptn),
     &              accluminosity,ptaccrecord(1,iptn),
     &              ptmassrecord(1,iptn),idonevalue
78789          FORMAT(I3,1X,9(1PE12.5),1X,I2)
               CALL FLUSH(52)
            ENDIF
         
         ELSEIF (istellarfeedback.EQ.5) THEN

            accluminosity = 0.

            IF (ptmasstime.LT.10.0) THEN

c               PRINT *, 'PTMASS TIME LESS THAN 10 YEARS'

               stellartime(iptn) = 0.
               stellarmass(iptn) = 0.
               stellarradius(iptn) = 0.
               stellarluminosity(iptn) = 0.

c               PRINT *, 'DISK MODEL DISABLED'

               diskinit(iptn) = 0.

            ELSEIF (ptmasstime.GE.10.0) THEN

c               PRINT *, 'PTMASS TIME GREATER THAN 10 YEARS'

               IF (stellartime(iptn).EQ.0.0) THEN
  
                  stellartime(iptn) = 0.1
                  stellarmass(iptn) = 1.0E-03*solarm
                  stellarradius(iptn) = 2.0*solarr
                  stellarluminosity(iptn) = 0.0
 
c                  PRINT *, 'DISK MODEL INITIALISING'

                  diskinit(iptn) = (xyzmh(4,listpm(iptn))*umass)
     &            -stellarmass(iptn)
                 
               ENDIF

c               PRINT *, 'DISK MODEL RUNNING'

               accrate = (1.0E-05*solarm/31557600.0)*
     &         ((xyzmh(4,listpm(iptn))*umass)-
     &         stellarmass(iptn))/diskinit(iptn)

               accluminosity = 1.0*((gg*stellarmass(iptn)*
     &         accrate)/stellarradius(iptn))

c               PRINT *, 'accrate (Msun/yr) = ',
c     &         accrate*31557600.0/solarm
c               PRINT *, 'accluminosity (Lsun) = ',
c     &         accluminosity/solarl

               stellartime(iptn) = stellartime(iptn)+dtstep
               stellarmass(iptn) = stellarmass(iptn)+
     &         (accrate*dtstep*31557600.0)
               stellarradius(iptn) = stellarradius(iptn)
               stellarluminosity(iptn) = stellarluminosity(iptn)

c               PRINT *, 'stellartime (yr) = ',
c     &         stellartime(iptn)
c               PRINT *, 'stellarmass (Msun) = ',
c     &         stellarmass(iptn)/solarm
c               PRINT *, 'stellarradius (Rsun) = ',
c     &         stellarradius(iptn)/solarr
c               PRINT *, 'stellarluminosity (Lsun) = ',
c     &         stellarluminosity(iptn)/solarl

            ENDIF

            ptmassluminosity(iptn) = accluminosity/
     &      (uergg*umass/utime)

         ELSEIF (istellarfeedback.EQ.6) THEN

            accluminosity = 0.

            IF (ptmasstime.LT.10.0) THEN
            
c               PRINT *, 'PTMASS TIME LESS THAN 10 YEARS'

               stellartime(iptn) = 0.
               stellarmass(iptn) = 0.
               stellarradius(iptn) = 0.
               stellarluminosity(iptn) = 0.

c               PRINT *, 'DISK MODEL DISABLED'

               diskinit(iptn) = 0.
                               
            ELSEIF (ptmasstime.GE.10.0) THEN

               PRINT *, 'PTMASS TIME GREATER THAN 10 YEARS'

               IF (stellartime(iptn).EQ.0.0) THEN
                 
                  PRINT *, 'STELLAR MODEL INITIALISING'          

                  nflag = 1
                  evolveto = 0.1
                  accrate = 1.0E-07
                  epsilon = 0.0
                  iptnlyon = listlyon(iptn)

                  CALL evollyon(nflag,iptnlyon,evolveto,accrate,
     &                 lyonradius(iptn),lyonluminosity(iptn),
     &                 lyontime(iptn),lyonmass(iptn),
     &                 file_stellarmodel_in(iptn),
     &                 file_stellarmodel_out(iptn))

                  file_stellarmodel_in(iptn) = 
     &                 file_stellarmodel_out(iptn)

                  stellartime(iptn) = lyontime(iptn)
                  stellarmass(iptn) = lyonmass(iptn)*solarm
                  stellarradius(iptn) = lyonradius(iptn)
                  stellarluminosity(iptn) = lyonluminosity(iptn)

                  PRINT *, 'DISK MODEL INITIALISING' 

                  diskinit(iptn) = (xyzmh(4,listpm(iptn))*umass)
     &            -stellarmass(iptn)
                
               ENDIF

               PRINT *, 'DISK MODEL RUNNING'

               accrate = (1.0E-05*solarm/31557600.0)*
     &         ((xyzmh(4,listpm(iptn))*umass)-
     &         stellarmass(iptn))/diskinit(iptn)

               accluminosity = 1.0*((gg*stellarmass(iptn)*
     &         accrate)/stellarradius(iptn))

               PRINT *, 'accrate (Msun/yr) = ',
     &         accrate/solarm*31557600.0
               PRINT *, 'accluminosity (Lsun) = ',
     &         accluminosity/solarl

               IF (stellartime(iptn).GE.(lyontime(iptn)+1.0)) THEN
 
                  PRINT *, 'STELLAR MODEL RUNNING'

                  nflag = 99 
                  evolveto = stellartime(iptn)+dtstep
                  accrate = accrate/solarm*31557600.0
                  epsilon = 0.02
                  iptnlyon = listlyon(iptn)

                  CALL evollyon(nflag,iptnlyon,evolveto,accrate,
     &                 lyonradius(iptn),lyonluminosity(iptn),
     &                 lyontime(iptn),lyonmass(iptn),
     &                 file_stellarmodel_in(iptn),
     &                 file_stellarmodel_out(iptn))

                  file_stellarmodel_in(iptn) = 
     &                 file_stellarmodel_out(iptn)

                  stellartime(iptn) = lyontime(iptn)
                  stellarmass(iptn) = lyonmass(iptn)*solarm
                  stellarradius(iptn) = lyonradius(iptn)
                  stellarluminosity(iptn) = lyonluminosity(iptn)

                  PRINT *, 'lyontime (yr) = ',
     &            lyontime(iptn)
                  PRINT *, 'lyonmass (Msun) = ',
     &            lyonmass(iptn)
                  PRINT *, 'lyonradius (Rsun) = ',
     &            lyonradius(iptn)/solarr
                  PRINT *, 'lyonluminosity (Lsun) = ',
     &            lyonluminosity(iptn)/solarl

               ELSE

                  PRINT *, 'STELLAR MODEL IDLING'

                  stellartime(iptn) = stellartime(iptn)+dtstep
                  stellarmass(iptn) = stellarmass(iptn)+
     &            (accrate*dtstep*31557600.0)
                  stellarradius(iptn) = stellarradius(iptn)
                  stellarluminosity(iptn) = stellarluminosity(iptn)

               ENDIF

c               PRINT *, 'stellartime (yr) = ',
c     &         stellartime(iptn)
c               PRINT *, 'stellarmass (Msun) = ',
c     &         stellarmass(iptn)/solarm
c               PRINT *, 'stellarradius (Rsun) = ',
c     &         stellarradius(iptn)/solarr
c               PRINT *, 'stellarluminosity (Lsun) = ',
c     &         stellarluminosity(iptn)/solarl

            ENDIF
c
c--Total sink luminosity is sum of accretion luminosity and stellar
c  model luminosity.            
c
            ptmassluminosity(iptn) =
     &      (stellarluminosity(iptn)+accluminosity)/
     &      (uergg*umass/utime)

         ELSE
            WRITE (iprint,*) 'ERROR - istellarfeedback ',
     &           istellarfeedback
            CALL quit(1)
         ENDIF

      END DO
c
c--For each star make list of particles that receive radiation from that star
c     Particles that receive the energy should not be blocked by another
c     particle (i.e. a closer particle whose smoothed region blocks the
c     particle in question.
c
      DO iptn = 1, nptmass
         iptcur = listpm(iptn)
c         print *,'Doing sink ',iptn,iptcur
         xipt = xyzmh(1,iptcur)
         yipt = xyzmh(2,iptcur)
         zipt = xyzmh(3,iptcur)
c
c--Find whether active particle is blocked or not
c
         DO n = 1, nexposedold(iptn)
            iexposedold(listexposedold(n,iptn)) = .TRUE.
         END DO

         nexposed(iptn) = 0
         nexposedolddoing = 0
         
c
c--Only expose particles doing full timestep
c
         DO n = 1, ncompact
            i = ivar(3,n)
            IF (ifullstep(i)) THEN
               IF (iexposedold(i)) nexposedolddoing = nexposedolddoing+1
               IF (iphase(i).EQ.0 .AND. rho(i).GT.bounddens) THEN
                  CALL getblock(i,xipt,yipt,zipt,numblock,neighblock,
     &                 xyzmh)
                  IF (numblock.EQ.0) THEN
                     nexposed(iptn) = nexposed(iptn) + 1
                     IF (nexposed(iptn).GT.nmaxexposed) THEN
                        WRITE (*,*) 'ERROR - nmaxexposed exceeded'
                        CALL quit(1)
                     ENDIF
                     listexposed(nexposed(iptn),iptn) = i

c            print *,xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),xyzmh(5,i),
c     &           SQRT(xyzmh(1,i)**2+xyzmh(2,i)**2+xyzmh(3,i)**2),
c     &           rho(i),i,' done'

                  ENDIF
               ENDIF
            ENDIF
         END DO

         DO n = 1, nexposedold(iptn)
            iexposedold(listexposedold(n,iptn)) = .FALSE.
         END DO
c
c--Need to assign fraction of stellar luminosity that goes to exposed particle
c     Use CURRENT values of kappa to calculate
c     Relative exposure goes as Luminosity*dt/(particle mass) * solid angle
c     where solid angle is 2*pi*(1.0-COS(ASIN(MIN(2*h,r)/r)))
c
c     To consider extinction, could also add a exp(-tau) term
c
         total = 0.0
         totaladd = 0.0
         DO n = 1, nexposed(iptn)
            i = listexposed(n,iptn)
            rad = SQRT((xyzmh(1,i)-xipt)**2 + (xyzmh(2,i)-yipt)**2 +
     &           (xyzmh(3,i)-zipt)**2)
c            fractionexposed(n) = (1.0-COS(ASIN(MIN(1.0,MIN(
c     &           2*xyzmh(5,i),rad)/rad))))/rho(i)
            fractionexposed(n) = (1.0-COS(ASIN(MIN(1.0,MIN(
     &           2*xyzmh(5,i),rad)/rad))))
c            fractionexposed(n) = (1.0-COS(ASIN(MIN(1.0,MIN(
c     &           2*xyzmh(5,i),rad)/rad))))/xyzmh(4,i)

            total = total + fractionexposed(n)
c
c
c--Works for both R-K and L-F Integrators
c     Particles on half-timestep receive energy from the star but this energy
c     isn't kept
c     'totaladd' is for keeping the amount of energy actually imparted to
c     the gas particles (not for half-timestep particles)
c
            IF (ifullstep(i)) totaladd = totaladd + 
     &           fractionexposed(n)*dtimax*isteps(i)

            print *,i,rad,rho(i),fractionexposed(n),total
         END DO
c
c--Check that energy injection is correct
c
         print *,' '
         print *,'NEXP ',nexposedolddoing, nexposedold(iptn),
     &        nexposed(iptn)
         print *,' '
         IF (icall.EQ.1 .OR. nexposedolddoing.EQ.nexposedold(iptn)) THEN
            IF (icall.NE.1) THEN
c
c--All exposed particles have just been calculated
c     Make sure that energy injected since last time is correct
c
c
c--This is the extra energy that needs to be added to be consistent with the
c     star's luminosity
c
               addenergy = cummulativeenergy(iptn) -
     &              actualcumenergy(iptn)
               print *,' Test ',cummulativeenergy(iptn),
     &              actualcumenergy(iptn)
c
c--Need to divide this extra energy up between particles
c     Depends on their exposure and on the timestep they are taking
c
               DO n = 1, nexposed(iptn)
                  i = listexposed(n,iptn)
c
c--Works for both R-K and L-F Integrators
c     Particles on half-timestep receive energy from the star but this energy
c     isn't kept
c
                  IF (ifullstep(i)) THEN
                     dti = dtimax*isteps(i)
                     addenergyparticle(i) = addenergyparticle(i) +
     &                    addenergy*fractionexposed(n)*dti/totaladd
                  ENDIF
               END DO
            ENDIF

            actualcumenergy(iptn) = 0.
            cummulativeenergy(iptn) = 0.
c
c--Save list of exposed particles for each sink particle
c
            DO n = 1, nexposed(iptn)
               listexposedold(n,iptn) = listexposed(n,iptn)
            END DO
            nexposedold(iptn) = nexposed(iptn)
         ENDIF
c
c--Set quantities to go into radiation hydrodynamics calculation as a source 
c     term
c
         DO n = 1, nexposed(iptn)
            i = listexposed(n,iptn)

            amount = ptmassluminosity(iptn)*fractionexposed(n)/total
            stellarenergy(i) = stellarenergy(i) + amount
c
c--Works for both R-K and L-F Integrators
c     Particles on half-timestep receive energy from the star but this energy
c     isn't kept
c
            IF (ifullstep(i)) THEN
               dti = dtimax*isteps(i)
               actualcumenergy(iptn) =actualcumenergy(iptn) +amount*dti
            ENDIF

c            print *,'Frac ',i,ptmassluminosity(iptn),
c     &           fractionexposed(n)/total,amount,amount*dti

         END DO
c
c--Record the actual amount of energy that should have gone into simulation
c     from each star
c
         cummulativeenergy(iptn) = cummulativeenergy(iptn) + 
     &        ptmassluminosity(iptn)*dtimax*istepmindone
      END DO
c
c--End sink particle radiation
c

 888  CONTINUE

      ihasghostcount = 0
      ihasghost = 0

      IF (ibound.EQ.100 .AND. iandrea.EQ.1)
     &     hmass_e = hillmass(planetmass(1), 0, xmass, 0) -planetmass(1)

C$OMP PARALLEL default(none)
C$OMP& shared(EU0,uergg,ekcle)
C$OMP& shared(vxyzu,lightspeed,uradconst,iprint,ddvtable)
C$OMP& shared(icall,dtimax,dtdo,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij,ijvar2,varij2,bounddens)
C$OMP& shared(xyzmh,dvtable,grwij,cnormk)
C$OMP& shared(dvdx,ihasghost,iflag,origEU,vsmooth,iphase,ifullstep)
C$OMP& shared(ntot,ncompactlocal,radkernel,gradhs,ncompact)
C$OMP& shared(ibound,hmass_e,e_planetesimal,iandrea)
C$OMP& shared(prmin,princ,pmmin,pminc,ptable,udist)
C$OMP& shared(dust_tk)
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend)
#endif
C$OMP& private(radius,index_r,index_m,y1,y2,y3,y4,w,v)
C$OMP& private(enp)
C$OMP& private(n,i,j,k,rhoi,icompact,pmi)
C$OMP& private(dvxdxi,dvxdyi,dvxdzi,dvydxi,dvydyi,dvydzi)
C$OMP& private(dvzdxi,dvzdyi,dvzdzi,dti,dx,dy,dz)
C$OMP& private(rij2,rij,rij1,dr,pmj,rhoj,hi,hj,hi21,hj21,hi41,hj41)
C$OMP& private(v2i,vi,v2j,vj,index,dxx,index1,dgrwdx,grwtij)
C$OMP& private(dWi,dWj,dvx,dvy,dvz)
C$OMP& private(rhomean,dvdotdr,dv,vmu,dvdWimj,dvdWimi,dvdWjmj)
C$OMP& private(dWidrlightrhorhom,pmjdWrijrhoi,dvdotdrs,dvs)
C$OMP& private(dWjdrlightrhorhom,dWiidrlightrhorhom,cv_effective)
C$OMP& private(pmjdWrunix,pmjdWruniy,pmjdWruniz,dvxs,dvys,dvzs)
C$OMP& private(iparttree)
C$OMP& private(dust_kappai,dust_cooling,heatingISRi,dust_gas)
C$OMP& reduction(+:ihasghostcount)
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
      DO n = numberstart, numberend
#else
      DO n = 1, ncompact
#endif
         i = ivar(3,n)
c
c--Needed for MPI code
c
         IF (i.GT.ntot) THEN
            iparttree = i + ntot + 2
         ELSE
            iparttree = i
         ENDIF

         IF (iphase(i).EQ.0) THEN
            EU0(1,i) = ekcle(1,i)
            EU0(2,i) = vxyzu(4,i)

            ekcle(3,i) = GETCV(rho(i),vxyzu(4,i))
            ekcle(2,i) = GETKAPPA(vxyzu(4,i),ekcle(3,i),rho(i))
c
c--Diffuse ISM: Set dust temperature and opacity
c
            IF (idustRT.GT.0 .AND. n.LE.ncompactlocal) THEN
               rhoi = rho(i)
               dust_tk(1,i) = dust_temperature(i,ntot,uradconst,
     &              ekcle(1,i),vxyzu(4,i),ekcle,rhoi,
     &              dust_kappai,dust_cooling,heatingISRi,dust_gas)
               dust_tk(2,i) = dust_kappai
            ENDIF
c
c--Note that CV and Kappa have already been done in ASS
c
            cv_effective = ekcle(3,i)/GET1OVERMU(rho(i),vxyzu(4,i))

            IF (n.LE.ncompactlocal) THEN
               IF(hasghost(i)) THEN
                  ihasghostcount = ihasghostcount + 1
               ENDIF

               IF (icall.EQ.1) THEN
                  IF (isteps(i).EQ.0) THEN
                     dti = dtdo*1.0d-12
                  ELSE
                     dti = dtimax*isteps(i)*1.0d-12
                  ENDIF
               ELSEIF (ifullstep(i)) THEN
                  dti = dtimax*isteps(i)
               ELSE
cTEST
                  dti = dtimax*isteps(i)/2.0
               ENDIF
            ENDIF

c
c--Addition of energy due to planetesimal collisions with planetary
c  atmosphere. Tables provided by Andrea Fortier.
c
            IF (ibound.EQ.100 .AND. iandrea.EQ.1) THEN
               radius = sqrt((xyzmh(1,i)-1.0)**2 + xyzmh(2,i)**2 +
     &                 xyzmh(3,i)**2)
               radius = log10(radius*udist)
               index_r = INT((radius - prmin)/princ) + 1
               index_m = INT((hmass_e - pmmin)/pminc) + 1
               IF (index_r.GE.rbins) index_r = rbins-1
               IF (index_m.GE.mbins) index_m = mbins-1
               IF (radius.LT.prmin) THEN
                  print *, 'Radius too small for Fortier tables'
                  CALL quit(1) 
               ENDIF
               
               y1=ptable(index_m,index_r)
               y2=ptable(index_m+1,index_r)
               y3=ptable(index_m+1,index_r+1)
               y4=ptable(index_m,index_r+1)

               IF (y1.LT.tiny .OR. y2.LT.tiny .OR. y3.LT.tiny .OR.
     &              y4.LT.tiny) THEN
                  e_planetesimal(i) = 0.0
                  goto 250
               ENDIF

               y1=log10(y1)
               y2=log10(y2)
               y3=log10(y3)
               y4=log10(y4)
               
               w = (log10(radius) - log10(prmin+(index_r-1)*princ))/
     &              (log10(prmin+(index_r*princ)) - 
     &              log10(prmin+((index_r-1)*princ)))
               v = (log10(hmass_e) - log10(pmmin+(index_m-1)*pminc))/
     &              (log10(pmmin+(index_m*pminc)) - 
     &              log10(pmmin+((index_m-1)*pminc)))
                  
               enp = (1.0-v)*(1.0-w)*y1 + v*(1.0-w)*y2 +
     &              w*v*y3 + (1.0-v)*w*y4
               e_planetesimal(i) = (10.0**enp)

 250           continue
            ENDIF

            dvxdxi = 0.
            dvxdyi = 0.
            dvxdzi = 0.
            dvydxi = 0.
            dvydyi = 0.
            dvydzi = 0.
            dvzdxi = 0.
            dvzdyi = 0.
            dvzdzi = 0.

            pmi = xyzmh(4,iparttree)
            hi = xyzmh(5,iparttree)
            hi21 = 1./(hi*hi)
            hi41 = hi21*hi21
            rhoi = rho(i)

            DO k = 1, ivar(1,n)
               icompact = ivar(2,n) + k
               j = ijvar(icompact)

               IF (iphase(j).EQ.0) THEN
c
c--Need to make sure that E and U values are loaded for non-active neighbours
c
                  EU0(1,j) = ekcle(1,j)
                  EU0(2,j) = vxyzu(4,j)
c
c--Note that CV and Kappa have already been done in ASS
c
                  cv_effective =ekcle(3,j)/GET1OVERMU(rho(j),vxyzu(4,j))
c
c--Calculate other quantities
c
                  dx = xyzmh(1,iparttree) - xyzmh(1,j)
                  dy = xyzmh(2,iparttree) - xyzmh(2,j)
                  dz = xyzmh(3,iparttree) - xyzmh(3,j)
                  rij2 = dx*dx + dy*dy + dz*dz + tiny
                  rij = SQRT(rij2)
                  rij1 = 1./rij
                  dr = rij

                  pmj = xyzmh(4,j)
                  rhoj = rho(j)

                  hj = xyzmh(5,j)
                  hj21 = 1./(hj*hj)
                  hj41 = hj21*hj21

                  v2i = rij2*hi21
                  vi = rij/hi

                  v2j = rij2*hj21
                  vj = rij/hj

                  IF (vi.LT.radkernel) THEN
                     index = v2i*ddvtable
                     dxx = v2i - index*dvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                     grwtij = (grwij(index) + dgrwdx*dxx)*hi41
                     dWi = grwtij * cnormk
                  ELSE
                     dWi = 0.
                  ENDIF

                  IF (vj.LT.radkernel) THEN
                     index = v2j*ddvtable
                     dxx = v2j - index*dvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                     grwtij = (grwij(index) + dgrwdx*dxx)*hj41
                     dWj = grwtij * cnormk
                  ELSE
                     dWj = 0.
                  ENDIF

                  dvx = vxyzu(1,i) - vxyzu(1,j)
                  dvy = vxyzu(2,i) - vxyzu(2,j)
                  dvz = vxyzu(3,i) - vxyzu(3,j)
                  IF (XSPH) THEN
                     dvxs = vsmooth(1,i) - vsmooth(1,j)
                     dvys = vsmooth(2,i) - vsmooth(2,j)
                     dvzs = vsmooth(3,i) - vsmooth(3,j)
                  ELSE
                     dvxs = dvx
                     dvys = dvy
                     dvzs = dvz
                  ENDIF

                  rhomean = 0.5*(rhoi+rhoj)

                  dvdotdr = dvx*dx + dvy*dy + dvz*dz
                  dv = dvdotdr/dr

                  dvdotdrs = dvxs*dx + dvys*dy + dvzs*dz
                  dvs = dvdotdr/dr

                  IF(dvdotdr.GT.0.0) THEN
cc                  IF(dvdotdrs.GT.0.0) THEN
                     vmu = 0.0
                  ELSE
                     vmu = dv
cc                     vmu = dvs
                  END IF

                  dvdWimj = pmj*dvs*dWi*gradhs(1,i)
                  dvdWimi = pmi*dvs*dWi*gradhs(1,i)
                  dvdWjmj = pmj*dvs*dWj*gradhs(1,j)

                  dWidrlightrhorhom = lightspeed*dWi/dr*pmj/(rhoi*rhoj)
                  dWiidrlightrhorhom = lightspeed*dWi/dr*pmi/(rhoi*rhoj)
                  dWjdrlightrhorhom = lightspeed*dWj/dr*pmj/(rhoi*rhoj)

                  pmjdWrijrhoi = pmj*dWi*rij1/rhoi
                  pmjdWrunix = pmjdWrijrhoi*dx
                  pmjdWruniy = pmjdWrijrhoi*dy
                  pmjdWruniz = pmjdWrijrhoi*dz
c
c--Calculates density(i) times the gradient of velocity
c
                  dvxdxi = dvxdxi - dvxs*pmjdWrunix
                  dvxdyi = dvxdyi - dvxs*pmjdWruniy
                  dvxdzi = dvxdzi - dvxs*pmjdWruniz
                  dvydxi = dvydxi - dvys*pmjdWrunix
                  dvydyi = dvydyi - dvys*pmjdWruniy
                  dvydzi = dvydzi - dvys*pmjdWruniz
                  dvzdxi = dvzdxi - dvzs*pmjdWrunix
                  dvzdyi = dvzdyi - dvzs*pmjdWruniy
                  dvzdzi = dvzdzi - dvzs*pmjdWruniz

                  varij(1,icompact) = rhoj
                  varij(2,icompact) = dWiidrlightrhorhom
                  varij(3,icompact) = dWidrlightrhorhom
                  varij(4,icompact) = dWjdrlightrhorhom

                  varij2(1,icompact) = pmjdWrunix
                  varij2(2,icompact) = pmjdWruniy
                  varij2(3,icompact) = pmjdWruniz
               ENDIF   ! iphase(j).EQ.0
            END DO

            dvdx(1,i) = dvxdxi
            dvdx(2,i) = dvxdyi
            dvdx(3,i) = dvxdzi
            dvdx(4,i) = dvydxi
            dvdx(5,i) = dvydyi
            dvdx(6,i) = dvydzi
            dvdx(7,i) = dvzdxi
            dvdx(8,i) = dvzdyi
            dvdx(9,i) = dvzdzi

            vari(1,n) = dti
            vari(2,n) = rhoi
         ENDIF   ! iphase(i).EQ.0
      END DO
C$OMP END DO
C$OMP END PARALLEL

#ifdef MPI
c
c--Need to add in non-local contributions to dvdx from particles on other 
c     MPI processes (note - only needed on local process, so don't need to
c     transfer total values back again).  Also, don't need dti except on 
c     local process.
c      
      istart = ntot + 1
      istartrec = istart + inumbertotal
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back forces, du, potential energy, dBevolxyz to be added on
c     to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUGRT
               print *,iproc,': expecting to get ',inumofsends,
     &              ' dvdx ',j
#endif
               CALL MPI_RECV(dvdx(1,istartrec),9*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,150,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.9*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.9*nnsentback dvdx'
                  CALL quit(1)
               ENDIF

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  DO k = 1, 9
                     dvdx(k,ipos) = dvdx(k,ipos) + dvdx(k,jpos)
                  END DO
               END DO

#ifdef MPIDEBUGRT
               print *,iproc,': got dvdx from ',iprocrec,' put into ',
     &              llistsentback(1,iprocrec+1)+1
#endif

            END DO
c
c--Other processes send the forces back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUGRT
                  print *,iproc,': sending dvdx to ',iprocsend,
     &                 ' istart ',istart,' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(dvdx(1,istart),9*inumber,MPI_REAL8,
     &                 iprocsend,150,MPI_COMM_WORLD,ierr)
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                  GOTO 7701
               ENDIF
            END DO
 7701       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO

c     
c--When using idustRT.GT.0, also need to transfer values of dust_tk
c
      IF (idustRT.GT.0) THEN
         inumbertotallocal = 0
         inumberreturned = 0
         DO i = 0, numproc - 1
            IF (iproc.EQ.i) THEN
               IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
c     is because here it is assumed that the total number of remote particles
c     on which forces need to be calculated by this MPI process is less
c     than idim.
c
                  DO j = 0, numproc - 1
                     IF (j.NE.iproc) THEN
c
c--Only sends anything at all to another MPI process if it needs to (the
c     other MPI process knows when it doesn't need to expect anything)
c
                        IF (nneighsentany(j+1)) THEN
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i2REAL8,indexDust,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i2REAL8,indexDust,ierr)
                CALL MPI_TYPE_COMMIT(indexDust,ierr)

                           CALL MPI_SEND(dust_tk,1,indexDust,j,161,
     &                        MPI_COMM_WORLD, ierr)

                           CALL MPI_TYPE_FREE(indexDust,ierr)
#ifdef MPIDEBUGRT
                           print *,iproc,' sent dust_tk (1)'
#endif
                        ENDIF
                     ENDIF
                  END DO
               ENDIF
c
c--Other processes receive the particles being sent
c
            ELSE
               IF (inumberreturned.LT.inumofreturns .AND.
     &              nneightogetback(i+1).GT.0) THEN
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUGRT
                  print *,iproc,': receiving dust_tk (1) data ',
     &                 inumberreturned,inumofreturns
#endif
                  istart = ntot + inumbertotallocal + 1
                  CALL MPI_RECV(dust_tk(1,istart), idim, MPI_REAL8,
     &                 i, 161, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
                  iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUGRT
                  print *,iproc,': got dust_tk (1) from ',iprocrec,
     &                 inumber
#endif
                  inumbertotallocal = inumbertotallocal + inumber/2
                  IF (inumberindiv(inumberreturned)*2.NE.inumber) THEN
                     WRITE (*,*) iproc,': ERR - inumberindiv ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                     CALL quit
                  ENDIF
                  IF (inumbercumm(inumberreturned).NE.istart) THEN
                     WRITE (*,*) iproc,': ERR - inumbercumm ',istart,
     &                    inumbercumm(inumberreturned),inumberreturned
                     CALL quit
                  ENDIF
                  IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                     WRITE (*,*) iproc,': ERR - inumberproc ',iprocrec,
     &                    inumberproc(inumberreturned),inumberreturned
                     CALL quit
                  ENDIF
               ENDIF
            ENDIF
         END DO
         IF (inumbertotallocal.NE.inumbertotal) THEN
            WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &           inumbertotallocal,inumbertotal
            CALL quit
         ENDIF
      ENDIF

#endif

      IF (ihasghostcount.GE.1) ihasghost = 1

C$OMP PARALLEL DO SCHEDULE(runtime)
C$OMP& shared(npart,nghost,ihasghost,ireal,dvdx)
C$OMP& private(i,j,k)
      DO i = npart + 1, npart + nghost*ihasghost
         j = ireal(i)
         DO k = 1, 9
            dvdx(k,i) = dvdx(k,j)
         END DO
      END DO
C$OMP END PARALLEL DO
c
c--Begin iterating
c
c      print *,iproc, 'begin iterations '
      maxerrE2last = 1.0E+30
      maxerrU2last = 1.0E+30

      DO nosweep = 1, nswmax
c         print *,iproc, 'it ',nosweep
c
c--Set error to zero for this iteration    
c
         maxerrE2 = 0.0
         maxerrU2 = 0.0
c
c--Calculate fluxlimiter values without using separate subroutine
c
C$OMP PARALLEL default(none)
C$OMP& shared(vari,ivar,varij2,ijvar,varij,uergg)
C$OMP& shared(dedxyz,ekcle,EU0,uradconst,lightspeed,rho,itime)
C$OMP& shared(nosweep,boundtemp,bounddens)
C$OMP& shared(npart,nghost,ihasghost,ireal,dvdx,ncompact)
C$OMP& shared(origEU,moresweep,xyzmh,gamma,alpha,beta)
C$OMP& shared(iflag,icall,ifsvi,alphaMMpass)
C$OMP& shared(iphase,maxerrU2last)
C$OMP& shared(stellarenergy,addenergyparticle,maxerrE2last,varinew)
C$OMP& shared(ncompactlocal,iscurrent)
C$OMP& shared(iiimax,xxmaxerrU2,xorigEU,xdUcomb,xdti,xU1i)
C$OMP& shared(xpres_numerator,xpres_denominator,xtfour)
C$OMP& shared(pdvvisc)
C$OMP& shared(ibound,iprint)
C$OMP& shared(rmind,deltar,zoverh,tprof,use_tprof)
C$OMP& shared(e_planetesimal,pcolen)
C$OMP& shared(iunique,iorig,n1,n2,ntot)
C$OMP& shared(udens,uergcc,utime,dust_tk,metallicity)
#ifdef MPI
C$OMP& shared(istart,istartrec,inumbertotal,inumberreturned)
C$OMP& shared(numproc,iproc,inumofsends,maxnneighsentback,istatus)
C$OMP& shared(ierr,ireturned,iprocrec,nneighsentback,ipos,jpos,inumber)
C$OMP& shared(llistsentback,inumofreturns,jjj,inumberproc,iprocsend)
C$OMP& shared(inumbertotallocal,nneighsentanyatall,nneighsentany)
C$OMP& shared(transfer,nneightogetback,inumberindiv,inumbercumm)
#endif
#ifdef MPICOPY
C$OMP& shared(numberstart,numberend,numbertodoherekeep,llisttrans)
C$OMP& shared(numproc,iproc,lblocklengths,indexMPI_INT1,i5REAL8)
C$OMP& shared(irecvcounti,idisplacementsi,realtransfer5to15,i2REAL8)
#endif
C$OMP& private(boundtest,radu)
C$OMP& private(bounddenslocal,izoverh)
C$OMP& private(i,j,k,n,dedxi,dedyi,dedzi,rhoi,pmj1,rhoj)
C$OMP& private(icompact,dti,pres_denominator,diffusion_numerator)
C$OMP& private(diffusion_denominator,pres_numerator,tfour)
C$OMP& private(radpresdenom,rhomean,dvdWimj,dvdWimi,dvdWjmj,vmu)
C$OMP& private(dWidrlightrhorhom,dWjdrlightrhorhom,dWiidrlightrhorhom)
C$OMP& private(cs,vpi,bi,bj,b1,gradEi2,rpdiag,rpall,gradvPi,vsig)
C$OMP& private(pmjdWrunix,pmjdWruniy,pmjdWruniz,Eij1)
C$OMP& private(betaval,chival,gammaval,u4term,u1term,u0term)
C$OMP& private(gradE1i,tsr1i,moresweep2,U1i,E1i,dUcomb,dEcomb)
C$OMP& private(presioverrhoi2,presjoverrhoj2,pmjdW,convection)
C$OMP& private(cvcold,ucold,ecold,residualE,residualU,alphamean)
C$OMP& private(stellarradiation,xchange,maxerrU2old,tempval1,tempval2)
C$OMP& private(gas_temp,xnH2,gas_dust_val,dustgammaval)
C$OMP& private(gas_dust_cooling,dust_tempi,dust_kappai)
C$OMP& private(cooling_line,cosmic_ray,photoelectric,dust_cooling)
C$OMP& private(heatingISRi,dust_gas,dust_term,betaval_d)
C$OMP& private(dust_heating,gas_dust_dustT,opacity)
C$OMP& private(ieqtype,rhoreal4,func_plus,func_old)
C$OMP& private(u_found,u_last,u_plus,cv1,cooling_line1,itry)
C$OMP& private(t_found,t_last,t_orig,t_plus,tdiff)
C$OMP& private(cooling_line2,dline_du,func,derivative,iterationloop)
C$OMP& private(photoelectric1,photoelectric2,dphoto_du)
C$OMP& private(errU2,errorEi,h2form,h2form1,h2form2,dh2form_dt)
#ifdef MPI
C$OMP& private(l)
#endif
#ifdef MPICOPY
C$OMP& private(ierr)
#endif
C$OMP& reduction(MAX:maxerrU2)
C$OMP& reduction(MAX:maxerrE2)

#ifdef MPICOPY
c
c--Need to zero all of varinew because calculated using reduction across
c     MPICOPY processes
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            varinew(1,i) = 0.
            varinew(2,i) = 0.
         END DO
C$OMP END DO
#endif

C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, ncompact
#endif
            i = ivar(3,n)

#ifndef MPICOPY
c
c--For MPICOPY code varinew is zeroed for entire list of particles above
c
            varinew(1,i) = 0.
            varinew(2,i) = 0.
#endif
c            varinew(3,i) = 0.
            IF (iphase(i).EQ.0 .AND. 
     &           (n2.EQ.0 .OR. iunique(iorig(i)).LE.n1)) THEN

               dedxi = 0.
               dedyi = 0.
               dedzi = 0.

               rhoi = vari(2,n)

               DO k = 1, ivar(1,n)
                  icompact = ivar(2,n) + k
                  j = ijvar(icompact)

                  IF (iphase(j).EQ.0) THEN
                     rhoj = varij(1,icompact)
                     pmjdWrunix = varij2(1,icompact)
                     pmjdWruniy = varij2(2,icompact)
                     pmjdWruniz = varij2(3,icompact)
c
c--Calculates the gradient of E (where E=rho*e)
c
                     Eij1 = rhoi*EU0(1,i) - rhoj*EU0(1,j)

                     dedxi = dedxi - Eij1*pmjdWrunix
                     dedyi = dedyi - Eij1*pmjdWruniy
                     dedzi = dedzi - Eij1*pmjdWruniz
                  ENDIF
               END DO

               dedxyz(1,i) = dedxi
               dedxyz(2,i) = dedyi
               dedxyz(3,i) = dedzi
            ELSE
               dedxyz(1,i) = 0.
               dedxyz(2,i) = 0.
               dedxyz(3,i) = 0.
            ENDIF
         END DO
C$OMP END DO

#ifdef MPI
C$OMP SINGLE
c
c--Need to add in non-local contributions to dedxyz from particles on other
c     MPI processes (note - only needed on local process, so don't need to 
c     transfer total values back again)
c
      istart = ntot + 1
      istartrec = istart + inumbertotal
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back forces, du, potential energy, dBevolxyz to be added on
c     to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUGRT
               print *,iproc,': expecting to get ',inumofsends,
     &              ' dedxyz ',j
#endif
               CALL MPI_RECV(dedxyz(1,istartrec),3*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,151,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.3*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.3*nnsentback dedxyz'
                  CALL quit(1)
               ENDIF

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  DO k = 1, 3
                     dedxyz(k,ipos) = dedxyz(k,ipos) + dedxyz(k,jpos)
                  END DO
               END DO

#ifdef MPIDEBUGRT
               print *,iproc,': got dedxyz from ',iprocrec,' put into ',
     &              llistsentback(1,iprocrec+1)+1
#endif

            END DO
c
c--Other processes send the forces back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUGRT
                  print *,iproc,': sending dedxyz to ',iprocsend,
     &                 ' istart ',istart,' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(dedxyz(1,istart),3*inumber,MPI_REAL8,
     &                 iprocsend,151,MPI_COMM_WORLD,ierr)
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                  GOTO 7702
               ENDIF
            END DO
 7702       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO

C$OMP END SINGLE
#endif

C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, ncompactlocal
#endif
            i = ivar(3,n)
            IF (iphase(i).EQ.0 .AND. 
     &           (n2.EQ.0 .OR. iunique(iorig(i)).LE.n1)) THEN

               gradE1i = SQRT(dedxyz(1,i)**2 + dedxyz(2,i)**2 + 
     &              dedxyz(3,i)**2)

               rhoi = vari(2,n)
c
c--If using diffuse ISM, use Rosseland mean opacity from the frequency
c     dependent opacity when dust temperatures are cold (T_d<100 K).
c     Otherwise use the tabulated grey dust opacities.
c
               IF (idustRT.GT.0) THEN
                  IF (dust_tk(1,i).LT.100.) THEN
                     opacity = dust_tk(2,i)
                  ELSE
                     opacity = ekcle(2,i)
                  ENDIF
               ELSE
                  opacity = ekcle(2,i)
               ENDIF

               IF (opacity.LE.0.) THEN
                  WRITE (*,*) 'ERROR - opacity1 ',opacity,
     &              dust_tk(1,i),ekcle(2,i),dust_tk(2,i),i
                  CALL quit(1)
               ENDIF

               tsr1i = ABS(gradE1i)/(EU0(1,i)*(rhoi**2)*opacity)

               ekcle(4,i) = (2. + tsr1i ) / (6. + 3.0*tsr1i + tsr1i**2) 
               ekcle(5,i) = ekcle(4,i) + ekcle(4,i)**2 * tsr1i**2
c
c--NOTE: ***** Forcing lambda and eddington to be 1/3 *****
c
c               IF (rhoi.GT.10000.0) THEN
c                  ekcle(4,i) = 1.0/3.0
c                  ekcle(5,i) = 1.0/3.0
c               ENDIF
            ELSE
               ekcle(4,i) = 1.0/3.0
               ekcle(5,i) = 1.0/3.0
            ENDIF
         END DO
C$OMP END DO

#ifdef MPICOPY
c
c--For MPICOPY job, need to transfer ekcle(4,i) values (only).
c
c--Transfer ekcle(4,i)
c
C$OMP DO SCHEDULE(runtime)
      DO i = numberstart, numberstart + numbertodoherekeep - 1
          llisttrans(i) = llisttrans(i)*5 + 3
      END DO
C$OMP END DO

C$OMP SINGLE
      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     MPI_REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
C$OMP END SINGLE

C$OMP DO SCHEDULE(runtime)
      DO i = numberstart, numberstart + numbertodoherekeep - 1
         llisttrans(i) = llisttrans(i)/5
      END DO
C$OMP END DO
c
c--All processes transfer their ekcle(4,i) data
c
C$OMP SINGLE
      CALL MPI_ALLGATHERV(ekcle,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,MPI_REAL8,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Put ekcle(4,i) data into correct places in ekcle
c
C$OMP DO SCHEDULE(runtime)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               ekcle(4,llisttrans(i) + 1) = realtransfer5to15(i)
            END DO
         ENDIF
      END DO
C$OMP END DO

C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': ekcle(4,i) move ',numbertodoherekeep
#endif

      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
C$OMP END SINGLE
#endif

#ifdef MPI
C$OMP SINGLE
c
c--Need to send back overall values of ekcle(4,i) only
c
      inumbertotallocal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
c     is because here it is assumed that the total number of remote particles
c     on which forces need to be calculated by this MPI process is less
c     than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
c
c--Only sends anything at all to another MPI process if it needs to (the
c     other MPI process knows when it doesn't need to expect anything)
c
                     IF (nneighsentany(j+1)) THEN
#ifdef MPIDEBUGRT
                        print *,iproc,': sending ekcle(4) to ',j
#endif
                        DO kkk = 1, nneighsentback(j+1)
                           transfer(kkk) = ekcle(4,
     &                          llistsentback(kkk,j+1)+1)
                        END DO

                        CALL MPI_SEND(transfer,nneighsentback(j+1),
     &                       MPI_REAL8,j,152,MPI_COMM_WORLD, ierr)
#ifdef MPIDEBUGRT
                        print *,iproc,' sent ekcle(4)'
#endif
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned = inumberreturned + 1
#ifdef MPIDEBUGRT
               print *,iproc,': receiving ekcle(4) data ',
     &              inumberreturned,inumofreturns
#endif
               CALL MPI_RECV(transfer, idim, MPI_REAL8,
     &              i, 152, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUGRT
               print *,iproc,': got ekcle(4) data ',iprocrec,inumber
#endif

               istart = ntot + inumbertotallocal
               DO jjj = 1, inumber
                  ekcle(4,istart+jjj) = transfer(jjj)
               END DO

#ifdef MPIDEBUGRT
               print *,iproc,': got ekcle(4) from ',iprocrec,inumber
#endif
               inumbertotallocal = inumbertotallocal + inumber
               IF (inumberindiv(inumberreturned).NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumbercumm(inumberreturned).NE.istart+1) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm ',istart+1,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberproc ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit(1)
      ENDIF

C$OMP END SINGLE
#endif

C$OMP DO SCHEDULE(runtime)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)

            ekcle(4,i) = ekcle(4,j)
            ekcle(5,i) = ekcle(5,j)
            dedxyz(1,i) = dedxyz(1,j)
            dedxyz(2,i) = dedxyz(2,j)
            dedxyz(3,i) = dedxyz(3,j)

            IF (idustRT.GT.0) THEN
               dust_tk(1,i) = dust_tk(1,j)
               dust_tk(2,i) = dust_tk(2,j)
            ENDIF
         END DO
C$OMP END DO

c
c--Particle I loop for calculating I-J quantities
c

C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, ncompact
#endif
            i = ivar(3,n)

            IF (iphase(i).EQ.0) THEN

               rhoi = vari(2,n)
c
c--NOTE: Needs to do this loop even for boundaryparticles because active
c     boundary particles will need to contribute to the varinew() 
c     quantities (i.e. diffusion terms) of particle j due to the way that
c     particle j only finds neighbours inside h_j or non-active particles
c     inside h_i.  The varinew() quantities of a boundaryparticle are
c     not used, but its contributions to j are.
c
c--Initialising counters to zero for this particle
c
               diffusion_numerator = 0.0
               diffusion_denominator = 0.0
c
c--All the neighbours loop
c
               DO k = 1, ivar(1,n)
                  icompact = ivar(2,n) + k
                  j = ijvar(icompact)
                  IF (iphase(j).EQ.0) THEN
                     rhoj = varij(1,icompact)
                     dWiidrlightrhorhom = varij(2,icompact)
                     dWidrlightrhorhom = varij(3,icompact)
                     dWjdrlightrhorhom = varij(4,icompact)
c
c--Set c*lambda/kappa*rho term for current quantities
c
                     IF (idustRT.GT.0) THEN
                        IF (dust_tk(1,i).LT.100.) THEN
                           opacity = dust_tk(2,i)
                        ELSE
                           opacity = ekcle(2,i)
                        ENDIF
                     ELSE
                        opacity = ekcle(2,i)
                     ENDIF
                     IF (opacity.LE.0.) THEN
                        WRITE (*,*) 'ERROR - opacity2 ',opacity,
     &                       dust_tk(1,i),ekcle(2,i),dust_tk(2,i),i
                        CALL quit(1)
                     ENDIF
                     bi = ekcle(4,i)/(opacity*rhoi)

                     IF (idustRT.GT.0) THEN
                        IF (dust_tk(1,j).LT.100.) THEN
                           opacity = dust_tk(2,j)
                        ELSE
                           opacity = ekcle(2,j)
                        ENDIF
                     ELSE
                        opacity = ekcle(2,j)
                     ENDIF
c
c--Check for zero or negative opacity, but only for particles in object
c     "n1", or ghosts if there is no object "n2" (typically particles in
c     "n2" (and their ghosts) have no opacity).
c
                     IF (opacity.LE.0.) THEN
                        IF (n2.EQ.0) THEN
                           WRITE (*,*) 'ERROR - opacity3 ',opacity,
     &                          dust_tk(1,j),ekcle(2,j),dust_tk(2,j),j,
     &                          npart
                           CALL quit(1)
                        ELSEIF (j.LT.npart) THEN
                           IF (iunique(iorig(j)).LE.n1) THEN
                              WRITE (*,*) 'ERROR - opacity4 ',opacity,
     &                             dust_tk(1,j),ekcle(2,j),dust_tk(2,j),
     &                             j,npart
                              CALL quit(1)
                           ENDIF
                        ENDIF
                     ENDIF
                     bj = ekcle(4,j)/(opacity*rhoj)
c
c--Choose the 'average' diffusion value.  The (bi+bj) quantity biased in
c     favour of the particle with the lowest opacity.  The other average
c     is that original recommended in Cleary & Monaghan for heat diffusion.
c
c               b1 = (4.0*bi*bj)/(bi+bj)
                     b1 = bi + bj
c	    
c--Diffusion numerator and denominator
c
                     diffusion_numerator = diffusion_numerator -
     &                    0.5*dWidrlightrhorhom*b1*EU0(1,J)*rhoj

                     diffusion_denominator = diffusion_denominator +
     &                    0.5*dWidrlightrhorhom*b1*rhoi
c
c--For current(j) need to add contribution due to i for h_j
c
                     IF (iscurrent(j)) THEN
                        tempval1 = 0.5*dWiidrlightrhorhom*b1
                        tempval2 = tempval1*rhoj
                        tempval1 = tempval1*EU0(1,i)*rhoi

C$OMP ATOMIC
                        varinew(1,j) = varinew(1,j) - tempval1
C$OMP ATOMIC
                        varinew(2,j) = varinew(2,j) + tempval2

                     ELSE
                        diffusion_numerator = diffusion_numerator -
     &                       0.5*dWjdrlightrhorhom*b1*EU0(1,J)*rhoj

                        diffusion_denominator = diffusion_denominator +
     &                       0.5*dWjdrlightrhorhom*b1*rhoi

                     ENDIF

                  ENDIF         !--iphase(j).EQ.0
               END DO           !J-loop

C$OMP ATOMIC
            varinew(1,i) = varinew(1,i) + diffusion_numerator
C$OMP ATOMIC
            varinew(2,i) = varinew(2,i) + diffusion_denominator

            ENDIF
         END DO   !--iphase(i).EQ.0
C$OMP END DO

#ifdef MPICOPY
c
c--For MPICOPY job, unfortunately need to do MPI reduce because in the above
c     loop particle 'i' gives varinew contributions to particle 'j'.
c
c--Transfer varinew(1-2,i)
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO k = 1, 2
               realtransfer5to15((n-1)*2+k) = varinew(k,i)
               realtransfer5to15(3*idim+(n-1)*2+k) = 0.
            END DO
         END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
         print *,iproc,': about to do varinew(1-2) reduce ',ncompact
#endif
c
c--All processes reduce the fxyzu(1-3) data
c
         CALL MPI_ALLREDUCE(realtransfer5to15,
     &        realtransfer5to15(3*idim+1),2*ncompact,
     &        MPI_REAL8,MPI_SUM,
     &        MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Now do reduction of varinew(1-2)
c
C$OMP DO SCHEDULE(runtime)
         DO n = 1, ncompact
            i = ivar(3,n)
            DO k = 1, 2
               varinew(k,i) = realtransfer5to15(3*idim+(n-1)*2 + k)
            END DO
         END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
         print *,iproc,': varinew(1-2) reduce ',ncompact
#endif

C$OMP END SINGLE
#endif

#ifdef MPI
c
c--Need to add contributions to varinew(1-4) from other MPI processes
c
C$OMP SINGLE
      istart = ntot + 1
      istartrec = istart + inumbertotal
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back forces, du, potential energy, dBevolxyz to be added on
c     to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUGRT
               print *,iproc,': expecting to get ',inumofsends,
     &              ' varinew ',j
#endif
               CALL MPI_RECV(varinew(1,istartrec),
     &              nelements_varinew*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,153,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.
     &              nelements_varinew*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nelements',
     &                 '*nnsentback varinw'
                  CALL quit(1)
               ENDIF

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1

                  DO k = 1, 2
                     varinew(k,ipos) = varinew(k,ipos) + varinew(k,jpos)
                  END DO
               END DO

#ifdef MPIDEBUGRT
               print *,iproc,': got varinew from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif

            END DO
c
c--Other processes send the values back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUGRT
                  print *,iproc,': sending varinew to ',iprocsend,
     &                 ' istart ',istart,' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(varinew(1,istart),
     &                 nelements_varinew*inumber,MPI_REAL8,
     &                 iprocsend,153,MPI_COMM_WORLD,ierr)
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                  GOTO 7703
               ENDIF
            END DO
 7703       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO
C$OMP END SINGLE
#endif

c
c--Particle I loop for calculating I only quantities
c
C$OMP DO SCHEDULE(runtime)
#ifdef MPICOPY
         DO n = numberstart, numberend
#else
         DO n = 1, ncompactlocal
#endif
            i = ivar(3,n)

            IF (iphase(i).EQ.0) THEN

            dti = vari(1,n)
            rhoi = vari(2,n)

            IF (.NOT.boundaryparticle(i,xyzmh,rhoi)) THEN

            diffusion_numerator = varinew(1,i)
            diffusion_denominator = varinew(2,i)
c            pres_numerator = varinew(3,i)
c            pres_numerator = 0.

            pres_numerator = pdvvisc(i)/dti

c            pres_denominator = varinew(4,i)
            pres_denominator = 0.
c
c--Radiation pressure...
c
            gradEi2 = (dedxyz(1,i)**2+dedxyz(2,i)**2+dedxyz(3,i)**2)
          
            IF (gradEi2.EQ.0.0) THEN
               gradvPi = 0.0
            ELSE
               rpdiag=0.5*(1.0-ekcle(5,i))
               rpall=0.5*(3.0*ekcle(5,i)-1.0)/gradEi2
               gradvPi=(((rpdiag+rpall*dedxyz(1,i)**2)*dvdx(1,i))+
     $              ((rpall*dedxyz(1,i)*dedxyz(2,i))*dvdx(2,i))+
     $              ((rpall*dedxyz(1,i)*dedxyz(3,i))*dvdx(3,i))+
     $              ((rpall*dedxyz(2,i)*dedxyz(1,i))*dvdx(4,i))+
     $              ((rpdiag+rpall*dedxyz(2,i)**2)*dvdx(5,i))+
     $              ((rpall*dedxyz(2,i)*dedxyz(3,i))*dvdx(6,i))+
     $              ((rpall*dedxyz(3,i)*dedxyz(1,i))*dvdx(7,i))+
     $              ((rpall*dedxyz(3,i)*dedxyz(2,i))*dvdx(8,i))+
     $              ((rpdiag+rpall*dedxyz(3,i)**2)*dvdx(9,i)))
            ENDIF ! gradEi2.EQ.0.0

            radpresdenom = gradvPi * EU0(1,i)

c         radpresdenom=0.0

            stellarradiation = (dti*stellarenergy(i) + 
     &           addenergyparticle(i))/xyzmh(4,i)
c            stellarradiation = 0.
c
c--For idustRT>0, replace the gas-dust coupling term in the u equation
c     by the dust-radiation term, but keep the T_d rather than T_g
c
            ieqtype = 0

 456        IF (idustRT.GT.0) THEN
               rhoreal4 = rhoi
               ekcle(2,i) = GETKAPPA(EU0(2,i),ekcle(3,i),rhoreal4)
               dust_tempi = dust_temperature(i,ntot,uradconst,
     &              EU0(1,i),EU0(2,i),ekcle,rhoi,dust_kappai,
     &              dust_cooling,heatingISRi,dust_gas)
               gas_temp = EU0(2,i)/ekcle(3,i)
               xnH2 = rhoi*udens/(gmw*mH)
            ENDIF ! idustRT.GT.0
c
c--For low densities and temperatures, use the form of the equations that
c     includes the gas-dust coupling term.  Also explicitly set the gas
c     opacity to zero (i.e. betaval=gammaval=tfour=0).  This works well
c     until the gas-dust coupling term gets very large, where upon the
c     convergence fails because even a small difference (<1 K) in the
c     gas and dust temperatures results in a large term.  The gas-dust
c     coupling term potentially becomes large for high densities, but 
c     a density criterion alone is not sufficient.  What we really want
c     to see is that the matter is well coupled with the radiation field
c     already.  So if we have high densities AND the difference between
c     the temperatures of the dust and local radiation field is small
c     (<1 K), then we abandon the gas-dust coupling term.
c
            IF (idustRT.GT.0 .AND. ((
     &           ABS(dust_tempi-(rhoi*EU0(1,i)/uradconst)**0.25).GT.1.0
     &           .AND.
     &           xnH2.LT.1.0E+11/metallicity)
     &           .OR. ieqtype.EQ.2)
     &           ) THEN
               ieqtype = ieqtype + 3

               betaval = 0.
               gammaval = 0.
               tfour = 0.

               betaval_d = 0.

               chival = dti*(diffusion_denominator-
     &              radpresdenom/EU0(1,I))-betaval - betaval_d

               gas_dust_val = 0.

               gas_dust_val = gas_dust_collisional_term(
     &              xnH2,metallicity,gas_temp) / rhoi/uergcc*utime
c
c--Implement crude effect of dust sublimation on gas-dust thermal coupling
c

               IF (gas_temp.GT.1000.0) gas_dust_val = MAX(gas_dust_val*
     &              (1500.0-gas_temp)/500.0, 0.0)
               gas_dust_cooling = gas_dust_val*(gas_temp - dust_tempi)
               gas_dust_dustT = gas_dust_val*dust_tempi

               cosmic_ray = cosmic_ray_heating(xnH2) /rhoi/uergcc*utime
c
c--This adds the dust heating term into the equations.  
c
               dust_heating = heatingISRi/uergg*utime

               cooling_line = 
     &              cooling_line_rate(i,gas_temp,xnH2,metallicity)/
     &              rhoi/uergcc*utime
c
c--H_2 formation heating
c
               IF (iH2formation_heating) THEN
                  h2form = h2_formation(i,gas_temp,dust_tempi,xnH2) *
     &                 4.48*eleccharge/rhoi/uergcc*utime
c
c--Potentially add heating from UV destruction and pumping of H_2
c     (see Bate 2015) for effect of this.
c
c     &                 + h2_destruction(i,xnH2,.FALSE.) *
c     &            (0.4 + 2.0/(1.0+criticaln(i,gas_temp)/(2.0*xnH2)))*
c     &                 eleccharge/rhoi/uergcc*utime
               ELSE
                  h2form = 0.
               ENDIF ! iH2formation_heating

               dustgammaval = 0.
               dust_term = 0.

               photoelectric =
     &              photoelectric_heating(i,gas_temp,xnH2,metallicity)/
     &              rhoi/uergcc*utime
c
c--In the low-density regime, the line cooling term (and sometimes 
c     other terms) can be very large, so need to solve for U1i using
c     Newton-Raphson so that it doesn't go negative
c
               IF (.TRUE.) THEN
                  itry = 1
 231              u_found = EU0(2,i)
                  rhoreal4 = rhoi
                  t_found = EU0(2,i)/GETCV(rhoreal4,u_found)
                  t_orig = t_found
                  DO iterationloop = 1, 100
                     u_last = u_found
                     t_last = t_found

                     u_found = getu(rhoreal4,t_found)
                     cv1 = GETCV(rhoreal4,u_found)
                     t_found = u_found/cv1
c
c--For calculating numerical derivative with gas temperature, 
c     set t_plus to be 1 K higher, or 1% higher if temperature is small
c
                     IF (t_found.GT.10.) THEN
                        t_plus = t_found + 1.0
                     ELSE
                        t_plus = t_found * 1.01
                     ENDIF ! t_found.GT.10.
                     tdiff = t_plus - t_found
                     u_plus = getu(rhoreal4,t_plus)
c
c--Molecular line cooling.  NOTE that because cooling_line_rate() sets the
c     abundances of carbon (in the chemistry() array) it it important that
c     the t_found call is done after the t_plus call otherwise the carbon
c     abundance that is stored in the chemistry() array is produced using
c     the wrong temperature!
c
                     cooling_line2 =
     &                    cooling_line_rate(i,t_plus,xnH2,
     &                    metallicity)/
     &                    rhoi/uergcc*utime
                     cooling_line1 =
     &                    cooling_line_rate(i,t_found,xnH2,
     &                    metallicity)/
     &                    rhoi/uergcc*utime

                     dline_du = (cooling_line2-cooling_line1)/tdiff
c
c--H_2 formation heating
c
                     IF (iH2formation_heating) THEN
                        h2form1 = h2_formation(i,t_found,dust_tempi,
     &                       xnH2)*4.48*eleccharge/
     &                       rhoi/uergcc*utime
c
c--Potentially add heating from UV destruction and pumping of H_2
c     (see Bate 2015) for effect of this.
c
c     &                          + h2_destruction(i,xnH2,.FALSE.) *
c     &               (0.4 + 2.0/(1.0+criticaln(i,t_found)/(2.0*xnH2)))*
c     &                          eleccharge/rhoi/uergcc*utime
                        h2form2 = h2_formation(i,t_plus,dust_tempi,
     &                       xnH2)*4.48*eleccharge/
     &                       rhoi/uergcc*utime
c
c--Potentially add heating from UV destruction and pumping of H_2
c     (see Bate 2015) for effect of this.
c
c     &                          + h2_destruction(i,xnH2,.FALSE.) *
c     &               (0.4 + 2.0/(1.0+criticaln(i,t_plus)/(2.0*xnH2)))*
c     &                          eleccharge/rhoi/uergcc*utime

                        dh2form_dt = (h2form2-h2form1)/tdiff
                     ELSE
                        h2form1 = 0.
                        h2form2 = 0.
                        dh2form_dt = 0.
                     ENDIF ! iH2formation_heating
c
c--Photoelectric heating
c
                     photoelectric1 =
     &                    photoelectric_heating(i,t_found,xnH2,
     &                    metallicity)/
     &                    rhoi/uergcc*utime
                     photoelectric2 =
     &                    photoelectric_heating(i,t_plus,xnH2,
     &                    metallicity)/
     &                    rhoi/uergcc*utime

                     dphoto_du = (photoelectric2-photoelectric1)/
     &                    tdiff
c
c--Now perform Newton-Raphson iteration
c
                     func = u_found + dti*gas_dust_val*
     &                    (u_found/cv1 - dust_tempi)
     &                    - origEU(2,i) - dti*pres_numerator
     &                    - dti*e_planetesimal(i) - dti*cosmic_ray
     &                    + dti*cooling_line1 - dti*photoelectric1
     &                    + dti*dust_term
     &                    - dti*h2form1

                     derivative = (u_plus - u_found)/tdiff +
     &                    dti*gas_dust_val
     &                    + dti*dline_du - dti*dphoto_du
     &                    - dti*dh2form_dt
c
c--If failed, do it again, but this time print out diagnostics
c
                     IF (itry.EQ.2) THEN
C$OMP CRITICAL(nrprint)
                        print *,'N-R',iterationloop,i,t_orig,rhoi
                        print *,'t_f,t_p',
     &                       t_found,t_plus
                        print *,'u_f,u_l,u_p',
     &                       u_found,u_last,u_plus,cv1
                        print *,
     &                     dti*gas_dust_val*(u_found/cv1 - dust_tempi),
     &                       origEU(2,i),-dti*cosmic_ray,
     &                       dti*cooling_line1,-dti*photoelectric1,
     &                       dust_tempi,-dti*pres_numerator,
     &                       dust_term,func,derivative
                        print *,u_plus,u_found,tdiff,
     &                       (u_plus - u_found)/tdiff,
     &                       dti*gas_dust_val,dti*dline_du,
     &                       - dti*dphoto_du,photoelectric2,
     &                       photoelectric1

C$OMP END CRITICAL(nrprint)
                     ENDIF ! itry.EQ.2
c
c--Limit the change to be no larger than a certain fraction each iteration
c
                     func_old = func
                     IF (func/derivative / t_found .GT. 0.3) THEN
                        func = 0.3*derivative * t_found
                     ELSEIF (func/derivative / t_found.LT.-0.3) THEN
                        func = -0.3*derivative * t_found
                     ENDIF ! func/derivative / t_found .GT. 0.3
                     t_found = t_found - func/derivative

                     IF (t_found.LT.1.) THEN
cC$OMP CRITICAL(nrprint)
c                           WRITE (81,*) 'T<1: ',i,t_found,
c     &                          t_found+func/derivative,func,derivative,
c     &                          func_old,- dti*h2form1,- dti*dh2form_dt
cC$OMP END CRITICAL(nrprint)
                        t_found = 1.
                     ENDIF ! t_found.LT.1

                     u_found = getu(rhoreal4,t_found)
c
c--Test for success
c
                     IF (ABS((t_found - t_last)/t_orig).LT.1E-3)
     &                    THEN
                        U1i = getu(rhoreal4,t_found)
                        ekcle(3,i) = cv1
                        photoelectric = photoelectric1
                        cooling_line = cooling_line1
                        h2form = h2form1
                        GOTO 233
                     ENDIF ! ABS((t_found - t_last)/t_orig).LT.1E-3
                  END DO
C$OMP CRITICAL(quart)
                  print *,"N-R failed for ieqtype = 3, try again"
            print *,"ngs ",u4term,u1term,u0term,betaval,chival,gammaval
            print *,"    ",ekcle(2,i),rhoi,dti
            print *,"    ",diffusion_denominator,diffusion_numerator
            print *,"    ",pres_denominator,pres_numerator,uradconst
            print *,"    ",radpresdenom,EU0(1,I),EU0(2,I),ekcle(3,i)
            print *,"    ",lightspeed,origEU(1,i),origEU(2,i)
            print *,"    ",dti*gas_dust_val*dust_tempi
            print *,"    ",dti*dust_heating,dti*gas_dust_val*dust_tempi
            print *,"    ",dti*cosmic_ray,dti*cooling_line
            print *,"    ",dti*photoelectric,dti*dust_term
            print *,"    ",ieqtype,u_found,u_last,u_plus,cv1
            print *,"    ",cooling_line1,cooling_line2,dline_du,func
            print *,"    ",derivative,h2form1,h2form2,dh2form_dt
C$OMP END CRITICAL(quart)
                  IF (itry.EQ.1) THEN
                     itry = 2
                     GOTO 231
                  ELSE
C$OMP CRITICAL (moresweepset)
                     moresweep=.TRUE.
C$OMP END CRITICAL (moresweepset)
                     GOTO 200
                  ENDIF ! itry.EQ.1
               ENDIF  ! Turn on/off N-R solve
c
c--Replaces the gas-dust coupling term in the u equation by the dust-radiation
c     term and assumes that T_g=T_d so that the dust-radiation term is
c     actually the gas-radiation term (i.e. uses kappa from opacity tables
c     which is based on the gas temperature, and uses (u/cv) rather than T_d.)
c
            ELSEIF (idustRT.GT.0) THEN
 777           ieqtype = ieqtype + 2

               gas_temp = EU0(2,i)/ekcle(3,i)
               xnH2 = rhoi*udens/(gmw*mH)
c
c--Use existing values of betaval, gammaval, chival, tfour
c     because the values from getkappa include both
c     the gas and dust opacities already.
c
               gas_dust_val = 0.
               dustgammaval = 0.
               gas_dust_cooling = 0.
               dust_heating = 0.
               dust_term = 0.
c
               cosmic_ray = cosmic_ray_heating(xnH2) /rhoi/uergcc*utime
c
c--This adds the dust heating term into the equations.  However, because the
c     above assumption is that T_g=T_d, this makes the low-density gas
c     as hot as the dust whereas in fact it should be a lot cooler.
c
               cosmic_ray = cosmic_ray + heatingISRi/uergg*utime

               cooling_line = 
     &              cooling_line_rate(i,gas_temp,xnH2,metallicity)/
     &              rhoi/uergcc*utime
c
c--H_2 formation heating
c
               IF (iH2formation_heating) THEN
                  h2form = h2_formation(i,gas_temp,dust_tempi,xnH2)*
     &                 4.48*eleccharge/rhoi/uergcc*utime
c
c--Potentially add heating from UV destruction and pumping of H_2
c     (see Bate 2015) for effect of this.
c
c     &                 + h2_destruction(i,xnH2,.FALSE.) *
c     &              (0.4 + 2.0/(1.0+criticaln(i,gas_temp)/(2.0*xnH2)))*
c     &                 eleccharge/rhoi/uergcc*utime
               ELSE
                  h2form = 0.
               ENDIF ! iH2formation_heating

               photoelectric =
     &              photoelectric_heating(i,gas_temp,xnH2,metallicity)/
     &              rhoi/uergcc*utime
c
c--Else, this is the original version of radiative transfer
c     (Whitehouse & Bate 2006), which does not include the 
c     diffuse ISM model of Bate (2015).
c
            ELSE
               ieqtype = 0

               dust_tempi = 0.
               gas_dust_val = 0.
               dustgammaval = 0.
               gas_dust_cooling = 0.

               cosmic_ray = 0.
               cooling_line = 0.
               photoelectric = 0.
               h2form = 0.

               dust_heating = 0.
               dust_term = 0.
            ENDIF ! idustRT.GT.0 .AND. ((
                  ! ABS(dust_tempi-(rhoi*EU0(1,i)/uradconst)**0.25).GT.1.0
                  ! .AND.
                  ! xnH2.LT.1.0E+11/metallicity)
                  ! .OR. ieqtype.EQ.2)
                
c
c--Now solve those equations...
c
            betaval = lightspeed*ekcle(2,i)*rhoi*dti
            chival = dti*(diffusion_denominator-radpresdenom/EU0(1,I))-
     &         betaval

            gammaval = uradconst*lightspeed*ekcle(2,i)/(ekcle(3,i))**4

            tfour=uradconst*lightspeed*ekcle(2,i)* 
     &         ((rhoi*EU0(1,i)/uradconst) - ((EU0(2,i)/ekcle(3,i))**4))

            u4term = gammaval*dti*(dti*(diffusion_denominator-
     &           radpresdenom/EU0(1,I)) - 1.0)
            u1term = (chival-1.0)*(1.0-dti*pres_denominator
     &           + dti*gas_dust_val/ekcle(3,i)
     &           ) 
     &           - betaval*dti*gas_dust_val/ekcle(3,i)
            u0term = betaval*(origEU(1,i) -
     &           dti*gas_dust_val*dust_tempi + dti*dust_heating) + 
     &           (chival-1.0)*(-origEU(2,i)-
     &           dti*pres_numerator - dti*e_planetesimal(i)
     &           - dti*gas_dust_val*dust_tempi - dti*cosmic_ray
     &           + dti*cooling_line - dti*photoelectric
     &           - dti*h2form
     &           + dti*dust_term) + 
     &           dti*diffusion_numerator*betaval +
     &           stellarradiation*betaval - (chival-1.0)*pcolen(i)

            IF (u1term.GT.0.0 .AND. u0term.GT.0.0 .OR. u1term.LT.0.0 
     &           .AND. u0term.LT.0.0) THEN
C$OMP CRITICAL(quart)
            print *,"ngs-2 ",u4term,u1term,u0term,betaval,chival
            print *,"    ",ekcle(2,i),rhoi,dti,gammaval
            print *,"    ",diffusion_denominator,diffusion_numerator
            print *,"    ",pres_denominator,pres_numerator,uradconst
            print *,"    ",radpresdenom,EU0(1,I),EU0(2,I),ekcle(3,i)
            print *,"    ",lightspeed,origEU(1,i),origEU(2,i)
            print *,"    ",dti*gas_dust_val*dust_tempi
            print *,"    ",dti*dust_heating,dti*gas_dust_val*dust_tempi
            print *,"    ",dti*cosmic_ray,dti*cooling_line
            print *,"    ",dti*photoelectric,dti*dust_term
            print *,"    ",ieqtype,dti*h2form,h2form,cooling_line
            print *,"    ",dust_tempi,origEU(2,i)/ekcle(3,i)
            print *,"    ",dti*diffusion_numerator*betaval,
     &           stellarradiation*betaval,gas_dust_val
C$OMP END CRITICAL(quart)
               IF (ieqtype.EQ.0 .OR. ieqtype.GE.3) THEN
C$OMP CRITICAL (moresweepset)
                  moresweep=.TRUE.
C$OMP END CRITICAL (moresweepset)
                  GOTO 200
               ELSE
                  GOTO 456
               ENDIF ! ieqtype.EQ.0 .OR. ieqtype.GE.3
            ENDIF ! u1term.GT.0.0 .AND. u0term.GT.0.0 .OR. 
                  ! u1term.LT.0.0 .AND. u0term.LT.0.0

            u1term = u1term/u4term
            u0term = u0term/u4term

            moresweep2 = .FALSE.

            CALL quartic_gs1t(u1term,u0term,EU0(2,I),U1i,
     &           moresweep2,i)

            IF (moresweep2) THEN
               IF (ieqtype.EQ.0 .OR. ieqtype.GE.3) THEN
C$OMP CRITICAL (moresweepset)
                  moresweep=.TRUE.
                  PRINT *,"Info: ",EU0(2,I)/ekcle(3,i),xyzmh(1,i),
     &                 xyzmh(2,i),xyzmh(3,i)
                  PRINT *,"Info2: ",u0term,u1term,u4term,gammaval,
     &                 ekcle(2,i),ekcle(3,i)
                  PRINT *,"Info3: ",chival,betaval,dti
                  PRINT *,"Info4: ",pres_denominator,
     &                 origEU(1,i),pres_numerator
                  PRINT *,"Info5: ",diffusion_numerator,
     &                 stellarradiation,diffusion_denominator
                  PRINT *,"Info6: ",radpresdenom,EU0(1,I)
                  PRINT *,"Info7: ",cosmic_ray,heatingISRi/uergg*utime
                  PRINT *,"Info8: ",cooling_line,photoelectric,h2form
C$OMP END CRITICAL (moresweepset)
                  GOTO 200
               ELSE
                  GOTO 456
               ENDIF ! ieqtype.EQ.0 .OR. ieqtype.GE.3
            ENDIF ! moresweep2

 233        CONTINUE
c
c - Use solution of quartic to solve for E
c
            E1i = (origEU(1,i)+dti*diffusion_numerator
     &           + gammaval*dti*U1i**4
     &           + dustgammaval*dti
     &           + dti*gas_dust_val*(U1i/ekcle(3,i)-dust_tempi)
     &           + dti*dust_heating
     &           + stellarradiation)/(1.0-chival)

            dUcomb = pres_numerator + pres_denominator*EU0(2,I) + tfour
     &           - gas_dust_cooling + cosmic_ray - cooling_line
     &           + photoelectric + h2form
     &           + e_planetesimal(i) + pcolen(i)
            dEcomb = diffusion_numerator + diffusion_denominator*
     &           EU0(1,I) - tfour - radpresdenom + stellarradiation
     &           + dust_heating 
     &           + gas_dust_cooling
c
c--Tests for negativity
c
            IF(U1i.LE.0.0) THEN
C$OMP CRITICAL (moresweepset)
               moresweep=.TRUE.
               PRINT *,"GSIMPL: U has gone negative ",i,xyzmh(1,i),
     &              xyzmh(2,i),xyzmh(3,i)
C$OMP END CRITICAL (moresweepset)
            ENDIF ! U1i.LE.0.0

            IF(E1i.LE.0.0) THEN
               IF (ieqtype.EQ.0 .OR. ieqtype.GT.3) THEN
C$OMP CRITICAL (moresweepset)
                  moresweep=.TRUE.
                  PRINT *,"GSIMPL: E has gone negative ",i,xyzmh(1,i),
     &                 xyzmh(2,i),xyzmh(3,i),dust_tempi,gas_temp,rhoi
C$OMP END CRITICAL (moresweepset)
               ELSEIF (ieqtype.EQ.3) THEN
                  GOTO 777
               ELSEIF (ieqtype.EQ.2) THEN
                  GOTO 456
               ELSE
                  WRITE (iprint,*) 'ERROR - ieqtype = ',ieqtype
                  CALL quit(1)
               ENDIF ! ieqtype.EQ.0 .OR. ieqtype.GT.3
            END IF ! E1i.LE.0.0
c
c--And the error is...
c
c            IF (EU0(2,I)/ekcle(3,i).GT.1.E+30) THEN
            IF (EU0(2,I)/ekcle(3,i).GT.0.0) THEN
               maxerrE2 = MAX(maxerrE2,1.0*ABS((EU0(1,I) - E1i) /E1i) 
     &              )

               residualE = 0.0
            ELSE
               xchange = ABS((origEU(1,i) + (dEcomb)*
     &              dti - E1i) /E1i)
               maxerrE2 = MAX(maxerrE2,xchange)

               residualE = origEU(1,i) + (dEcomb)*dti - E1i
            ENDIF ! EU0(2,I)/ekcle(3,i).GT.0.0

c            maxerrE2 = MAX(maxerrE2,ABS((origEU(1,i) + (dEcomb)*
c     &           dti - EU0(1,I)) /EU0(1,I)))

c            maxerrU2 = MAX(maxerrU2,ABS((origEU(2,i)+(dUcomb)* 
c     &           dti - EU0(2,I))/EU0(2,I)))


            IF (EU0(2,I)/ekcle(3,i).GT.2000.0) THEN
c            IF (EU0(2,I)/ekcle(3,i).GT.1.0E+30) THEN
               maxerrU2 = MAX(maxerrU2,1.0*ABS((EU0(2,I) - U1i) /U1i)
     &              )

               residualU = 0.0
            ELSE
               maxerrU2old = maxerrU2
               maxerrU2 = MAX(maxerrU2,ABS((origEU(2,i)+(dUcomb)* 
     &              dti - U1i)/U1i))

               IF (maxerrU2.NE.maxerrU2old) THEN
C$OMP CRITICAL (storeUerror)
                  iiimax = i
                  xxmaxerrU2 = maxerrU2
                  xorigEU = origEU(2,i)
                  xdUcomb = dUcomb
                  xdti = dti
                  xU1i = U1i
                  xpres_numerator = pres_numerator
                  xpres_denominator = pres_denominator*EU0(2,I)
                  xtfour = tfour
C$OMP END CRITICAL (storeUerror)
               ENDIF ! maxerrU2.NE.maxerrU2old
               residualU = origEU(2,i)+(dUcomb)*dti - U1i
            ENDIF ! EU0(2,I)/ekcle(3,i).GT.2000.0)
c
c--Copy values
c
            IF (E1i.GT.0.) EU0(1,I) = E1i
            EU0(2,I) = U1i

            ekcle(3,i) = GETCV(rho(i),EU0(2,i))
            ekcle(2,i) = GETKAPPA(EU0(2,i),ekcle(3,i),rho(i))

            IF (idustRT.GT.0) THEN
               rhoi = rho(i)
               dust_tk(1,i) = dust_temperature(i,ntot,uradconst,
     &              ekcle(1,i),EU0(2,i),ekcle,rhoi,
     &              dust_kappai,dust_cooling,heatingISRi,dust_gas)
               dust_tk(2,i) = dust_kappai
            ENDIF ! idustRT.GT.0      

 200        CONTINUE
            ENDIF ! .NOT.boundaryparticle(i,xyzmh,rhoi)

            ENDIF ! iphase(i).EQ.0
         END DO ! I-loop
C$OMP END DO

#ifdef MPICOPY
C$OMP SINGLE
c
c--Need to transfer ekcle(2,i) around MPI processes.  Need to do it
c         here to make sure that the right ekcle(2,i) values are stored
c         on exit, and transferred to ghosts.
c
c         Also need to transfer ekcle(3,i) if another iteration is to be
c         done, so do it here at the same time.
c
c         And need to transfer EU0(1-2,i) so that when exits, can load
c         into ekcle and vxyzu, or, if another iteration is do, so that
c         the latest values for all particles are used.
c
      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     i5REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their ekcle(1-5) data
c
      CALL MPI_ALLGATHERV(ekcle,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,i5REAL8,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Put ekcle(1-5) data into correct places in ekcle
c
C$OMP DO SCHEDULE(runtime)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 5
                  ekcle(k,llisttrans(i) + 1) =
     &                 realtransfer5to15((i-1)*5 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': ekcle(1-5) move ',numbertodoherekeep
#endif
      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
c
c--Transfer EU0(1-2,i)
c
      CALL MPI_TYPE_INDEXED(numbertodoherekeep,
     &     lblocklengths, llisttrans(numberstart),
     &     i2REAL8, indexMPI_INT1, ierr)
      CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
c
c--All processes transfer their EU0(1-2) data
c
      CALL MPI_ALLGATHERV(EU0,1,indexMPI_INT1,
     &     realtransfer5to15,irecvcounti,idisplacementsi,i2REAL8,
     &     MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Put EU0(1-2) data into correct places in EU0
c
C$OMP DO SCHEDULE(runtime)
      DO j = 0, numproc - 1
         IF (j.NE.iproc) THEN
            DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &           irecvcounti(j+1)
               DO k = 1, 2
                 EU0(k,llisttrans(i) + 1) =
     &                 realtransfer5to15((i-1)*2 + k)
               END DO
            END DO
         ENDIF
      END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
      print *,iproc,': EU0(1-2) move ',numbertodoherekeep
#endif
C$OMP END SINGLE
      IF (idustRT.GT.0) THEN
c
c--Transfer dust_tk(1-2,i)
c
c--All processes transfer their dust_tk(1-2) data
c
C$OMP SINGLE
         CALL MPI_ALLGATHERV(dust_tk,1,indexMPI_INT1,
     &        realtransfer5to15,irecvcounti,idisplacementsi,i2REAL8,
     &        MPI_COMM_WORLD,ierr)
C$OMP END SINGLE
c
c--Put dust_tk(1-2) data into correct places in dust_tk
c
C$OMP DO SCHEDULE(runtime)
         DO j = 0, numproc - 1
            IF (j.NE.iproc) THEN
               DO i = idisplacementsi(j+1)+1,idisplacementsi(j+1) +
     &              irecvcounti(j+1)
                  DO k = 1, 2
                     dust_tk(k,llisttrans(i) + 1) =
     &                    realtransfer5to15((i-1)*2 + k)
                  END DO
               END DO
            ENDIF
         END DO
C$OMP END DO
C$OMP SINGLE
#ifdef MPIDEBUG
         print *,iproc,': dust_tk(1-2) move ',numbertodoherekeep
#endif
C$OMP END SINGLE
      ENDIF
C$OMP SINGLE
      CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
C$OMP END SINGLE
#endif

C$OMP DO SCHEDULE(runtime)
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            ekcle(2,i) = ekcle(2,j)
            IF (idustRT.GT.0) THEN
               dust_tk(1,i) = dust_tk(1,j)
               dust_tk(2,i) = dust_tk(2,j)
            ENDIF
         END DO
C$OMP END DO
C$OMP END PARALLEL

#ifdef MPIALL
c
c--Need to do global check for moresweep, and for maxerrE2, maxerrU2 to check
c     for errors and convergence across all MPI processes
c
c--Applies to BOTH MPI and MPICOPY
c
         CALL MPI_ALLREDUCE(moresweep,moresweeptot,1,MPI_LOGICAL,
     &        MPI_LOR,MPI_COMM_WORLD,ierr)
         moresweep = moresweeptot
         CALL MPI_ALLREDUCE(maxerrE2,maxerrE2tot,1,MPI_REAL8,MPI_MAX,
     &        MPI_COMM_WORLD,ierr)
         maxerrE2 = maxerrE2tot
         CALL MPI_ALLREDUCE(maxerrU2,maxerrU2tot,1,MPI_REAL8,MPI_MAX,
     &        MPI_COMM_WORLD,ierr)
         maxerrU2 = maxerrU2tot
#endif


#ifdef MPIALL
         IF (iproc.EQ.0) THEN
#endif
c         IF (MOD(nosweep,1).EQ.0) THEN
c           PRINT 97001,nosweep,nswmax2,maxerrE2,maxerrU2
c97001      FORMAT('GSIMPL: Finished iteration ',I4,' of ',I4,' (',
c     &        1PE12.5,1PE12.5,')')
c         ENDIF
#ifdef MPIALL
         ENDIF
#endif

         IF (moresweep) GOTO 999
c
c--The actual test
c
         IF(maxerrE2.LE.tolerance_rt.AND.maxerrU2.LE.tolerance_rt) THEN
c            PRINT *,"Complete with ",nosweep," iterations"
            GOTO 150
         ENDIF
         maxerrU2last = maxerrU2
c
c--Test for convergence to non-zero value (incorrect minimum)
c     Must have equal value at least twice to stop detecting up and down
c     as non-convergence
c
c         GOTO 333
            
         DO itest = ntests,1,-1

c            GOTO 332

            IF (nosweep.GT.10) THEN
               IF (maxerrE2.GT.0.99999*xmaxerr1old(itest).AND.
     &              maxerrE2.LT.1.00001*xmaxerr1old(itest)) THEN
                  DO iii = 1, numcomp
                     IF (maxerrE2.GT.0.99999*xmaxerrcomp1(iii).AND.
     &                    maxerrE2.LT.1.00001*xmaxerrcomp1(iii)) THEN
                        IF (irtnoisy) THEN
                         PRINT *,"GSIMPL: Non-convergence1 ",maxerrE2,
     &                   xmaxerrcomp1(iii),(xmaxerr1old(ii),ii=1,ntests)
                        ENDIF
                        nrtnc1 = nrtnc1 + 1
                        moresweep = .TRUE.
                        GOTO 998
                     ENDIF
                  END DO
                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
                  xmaxerrcomp1(numcomp) = maxerrE2
                  IF (irtnoisy) THEN
                    PRINT *,"GSIMPL: Almost Non-convergence ",
     &                   maxerrE2,(xmaxerr1old(ii),ii=1,ntests)
                  ENDIF
                  nrtnca = nrtnca + 1
               ENDIF
               IF (maxerrU2.GT.0.99999*xmaxerr2old(itest).AND.
     &              maxerrU2.LT.1.00001*xmaxerr2old(itest)) THEN
                  DO iii = 1, numcomp
                     IF (maxerrU2.GT.0.99999*xmaxerrcomp2(iii).AND.
     &                    maxerrU2.LT.1.00001*xmaxerrcomp2(iii)) THEN
                        IF (irtnoisy) THEN
                         PRINT *,"GSIMPL: Non-convergence1-2 ",maxerrU2,
     &                  xmaxerrcomp2(iii),(xmaxerr2old(ii),ii=1,ntests)
                        ENDIF
                        nrtnc2 = nrtnc2 + 1 
                        moresweep = .TRUE.
                        GOTO 998
                     ENDIF
                  END DO
                  numcomp = MAX(1,MOD(numcomp + 1,ntests + 1))
                  xmaxerrcomp2(numcomp) = maxerrU2
                  IF (irtnoisy) THEN
                    PRINT *,"GSIMPL: Almost Non-convergence ",
     &                   maxerrU2,(xmaxerr2old(ii),ii=1,ntests)
                  ENDIF
                  nrtnca = nrtnca + 1
               ENDIF
            ENDIF
c 332        CONTINUE
            IF (itest.NE.1) THEN
               xmaxerr1old(itest) = xmaxerr1old(itest-1)
               xmaxerr2old(itest) = xmaxerr2old(itest-1)
            ELSE
               xmaxerr1old(itest) = maxerrE2
               xmaxerr2old(itest) = maxerrU2
            ENDIF
         END DO
c 333     CONTINUE

#ifdef MPI
c
c--If another iteration is required, transfer new values of EU0, ekcle(2-3)
c     to other MPI processes
c
      inumbertotallocal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
c     is because here it is assumed that the total number of remote particles
c     on which forces need to be calculated by this MPI process is less
c     than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
c
c--Only sends anything at all to another MPI process if it needs to (the
c     other MPI process knows when it doesn't need to expect anything)
c
                     IF (nneighsentany(j+1)) THEN
#ifdef MPIDEBUGRT
                        print *,iproc,': sending EU0 to ',j,
     &                     ' starting at ',llistsentback(1,j+1),'+1 ',
     &                       nneighsentback(j+1),' list '
c     &             ,(llistsentback(kkk,j+1),kkk=1,nneighsentback(j+1))
#endif
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i2REAL8,indexEU,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i2REAL8,indexEU,ierr)
                CALL MPI_TYPE_COMMIT(indexEU,ierr)

                        CALL MPI_SEND(EU0,1,indexEU,j,154,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexEU,ierr)
#ifdef MPIDEBUGRT
                        print *,iproc,' sent EU0'
#endif
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i5REAL8,indexekcle,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i5REAL8,indexekcle,ierr)
                        CALL MPI_TYPE_COMMIT(indexekcle,ierr)

                        CALL MPI_SEND(ekcle,1,indexekcle,j,155,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexekcle,ierr)
#ifdef MPIDEBUGRT
                        print *,iproc,' sent ekcle'
#endif
                        IF (idustRT.GT.0) THEN
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i2REAL8,indexDust,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i2REAL8,indexDust,ierr)
                CALL MPI_TYPE_COMMIT(indexDust,ierr)

                        CALL MPI_SEND(dust_tk,1,indexDust,j,159,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexDust,ierr)
#ifdef MPIDEBUGRT
                        print *,iproc,' sent dust_tk'
#endif
                        ENDIF
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned = inumberreturned + 1
#ifdef MPIDEBUGRT
               print *,iproc,': receiving EU0 data ',
     &              inumberreturned,inumofreturns
#endif
               istart = ntot + inumbertotallocal + 1
               CALL MPI_RECV(EU0(1,istart), idim, MPI_REAL8,
     &              i, 154, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUGRT
               print *,iproc,': got EU0 from ',iprocrec,inumber
#endif
               CALL MPI_RECV(ekcle(1,istart), idim, MPI_REAL8,
     &              i, 155, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUGRT
               print *,iproc,': got ekcle from ',iprocrec,inumber
#endif
               inumbertotallocal = inumbertotallocal + inumber/5
               IF (inumberindiv(inumberreturned)*5.NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF

               IF (idustRT.GT.0) THEN
                  CALL MPI_RECV(dust_tk(1,istart), idim, MPI_REAL8,
     &                 i, 159, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
                  iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUGRT
                  print *,iproc,': got dust_tk from ',iprocrec,inumber
#endif
               ENDIF

               IF (inumbercumm(inumberreturned).NE.istart) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm ',istart,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberproc ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit(1)
      ENDIF

#endif

      END DO ! Iterations loop
c
c--Maximum number of iterations reached
c
      IF (irtnoisy) THEN
        PRINT *,"GSIMPL: Warning. Maximum iterations reached"
      ENDIF
      nrtmaxit = nrtmaxit + 1
      IF (iphase1.GE.10) THEN
         moresweep = .TRUE.
         GOTO 998
         CALL quit(1) 
      ELSE
         iphase1 = iphase1 + 1
         moresweep = .TRUE.
         GOTO 998
      ENDIF
c
c--Output success
c
 150  CONTINUE
      IF (irtnoisy) THEN
#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
        WRITE (*,38001) nosweep,maxerrE2,maxerrU2
#ifdef MPI
      ENDIF
#endif
      ENDIF
      nintegsuc = nintegsuc + 1
      nintegtot = nintegtot + nosweep
38001 FORMAT('Succeeded with ',I4,' iterations E U ',2(1PE12.5,1X))

 998  CONTINUE

      nit = nosweep
      errorE = maxerrE2
      errorU = maxerrU2
c
c--And that done, return everything to ASS
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompactlocal,ivar,vxyzu,EU0,iflag,vari)
C$OMP& shared(uergg,boundtemp,uradconst,rho,ekcle)
C$OMP& shared(iphase,bounddens)
C$OMP& shared(hoverr,centralmass,ibound,xyzmh)
C$OMP& shared(rmind,deltar,zoverh,tprof,use_tprof)
C$OMP& shared(gamma,dedxyz)
C$OMP& shared(iunique,iorig,n1,n2,uzero_n2)
C$OMP& shared(ifullstep,pcolen)
C$OMP& shared(dust_tk,ntot)
C$OMP& private(temp, rho0)
C$OMP& private(boundtest,radu)
C$OMP& private(bounddenslocal,izoverh)
C$OMP& private(n,i)
C$OMP& private(cvcold,ucold,ecold,cv_effective)
C$OMP& private(rhoi,dust_kappai,dust_cooling,heatingISRi,dust_gas)
      DO n = 1, ncompactlocal
         i = ivar(3,n)
         IF (iphase(i).EQ.0) THEN

            IF (.NOT.boundaryparticle(i,xyzmh,vari(2,n))) THEN
               ekcle(1,i) = EU0(1,i)
               vxyzu(4,i) = EU0(2,i)
            ELSEIF (n2.GT.0) THEN
               vxyzu(4,i) = uzero_n2
               ekcle(1,i) = uradconst*(boundtemp)**4/rho(i)
               ekcle(4,i) = 0.
               ekcle(5,i) = 0.
               dedxyz(1,i) = 0.
               dedxyz(2,i) = 0.
               dedxyz(3,i) = 0.
               IF (idustRT.GT.0) THEN
                  rhoi = rho(i)
                  dust_tk(1,i) = dust_temperature(i,ntot,uradconst,
     &                 ekcle(1,i),vxyzu(4,i),ekcle,rhoi,
     &                 dust_kappai,dust_cooling,heatingISRi,dust_gas)
                  dust_tk(2,i) = dust_kappai
               ENDIF
            ELSE
               cvcold = 1.5*Rg/(gmw*uergg)

               IF (ibound/10.EQ.10) THEN
                  radu = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 + 
     &                 xyzmh(3,i)**2)
                  IF (use_tprof) THEN
                     temp = gmw*centralmass*hoverr**2*radu**
     &                    (tprof+1)/((Rg/uergg)*gamma*radu)
                     ucold = getu(rho(i), temp)
                  ELSE
                     ucold = centralmass*hoverr**2/
     &                    (radu*gamma*(gamma-1.0))
                  ENDIF
                  ecold = uradconst*(ucold/GETCV(rho(i),ucold))**4/
     &                 rho(i)
               ELSE
                  ucold = boundtemp*cvcold
                  ecold = uradconst*(boundtemp)**4/rho(i)
               ENDIF

               ekcle(1,i) = ecold

               rho0 = 0.05
               IF(ibound/10.EQ.10 .AND. radu.LT.0.385 .AND.
     &              rho(i).GT.rho0 .AND. use_tprof) THEN
                  vxyzu(4,i) = ucold*(rho(i)/rho0)**0.8
               ELSE
                  vxyzu(4,i) = ucold
               ENDIF
               ekcle(3,i) = GETCV(rho(i),vxyzu(4,i))
               ekcle(2,i) = GETKAPPA(vxyzu(4,i),ekcle(3,i),rho(i))
               IF (idustRT.GT.0) THEN
                  rhoi = rho(i)
                  dust_tk(1,i) = dust_temperature(i,ntot,uradconst,
     &                 ekcle(1,i),vxyzu(4,i),ekcle,rhoi,
     &                 dust_kappai,dust_cooling,heatingISRi,dust_gas)
                  dust_tk(2,i) = dust_kappai
               ENDIF
            ENDIF
c
c--Zero collision energy array after its addition.
c
            IF (ifullstep(i)) THEN
               pcolen(i) = 0.0
            ENDIF

            cv_effective = ekcle(3,i)/GET1OVERMU(rho(i),vxyzu(4,i))

         ENDIF
      END DO
C$OMP END PARALLEL DO

 999  CONTINUE

#ifdef MPIALL
      CALL MPI_TYPE_FREE(i2REAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)
#endif

      RETURN

      END !SUBROUTINE GSIMPL
