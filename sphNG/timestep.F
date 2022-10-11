       SUBROUTINE timestep (dt,idtsyn,nlst,llist,fxyzu,Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine calculates the new timestep for each     *
c     particle in the list.                                 *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      PARAMETER (Courant=0.30)

      DIMENSION llist(idim),fxyzu(4,idim)
      DIMENSION Bevolxyz(imhdevol,imhd)

      REAL omega, omega1

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/cylinder'
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/timecrit'

      REAL rmodall(15)

      xlog2 = 0.30103
      istepmin = imax
      istepmingasnew = imax
      istepmax = -1
      istepeta = 0
      etafac = 0.

C$OMP PARALLEL default(none)
C$OMP& shared(nlst,llist)
C$OMP& shared(xyzmh,vxyzu)
C$OMP& shared(fxyzu,Bxyz,iexf)
C$OMP& shared(dt,isteps,imaxstep,divv,rho)
C$OMP& shared(alpha,beta,vsound,xlog2)
C$OMP& shared(idtsyn,nearpt,nptlist)
C$OMP& shared(iprint,nneigh,iphase)
C$OMP& shared(iorig,nptmass,listpm,listrealpm,encal)
C$OMP& shared(ifsvi,alphaMM,iptsoft,ptsoft,iresist,etamhd)
C$OMP& shared(iptmass,hacc,haccall,ptmcrit,radcrit)
C$OMP& shared(iplanetesimals,gasdrag,planetesimaltimestep)
C$OMP& shared(timestepcriterion)
C$OMP& private(rmodgd)
C$OMP& private(i,j,l,ll,xmindist,iptcur,rad2)
C$OMP& private(divvi,aux1,aux2,aux3,denom,valfven2,rmodcool)
C$OMP& private(crstepi,rmodcr,rmodvel,force2,rmod,vel2,coolstepi)
C$OMP& private(stepi,ibin,istep2,irat,softrad,rad,omega,omega1,tcool)
C$OMP& private(Bxi,Byi,Bzi,xi,yi,zi,dti,dteta,rmodeta,etai)
C$OMP& private(rmodall,rmod_du)
C$OMP& reduction(MAX:istepmax)
C$OMP& reduction(MIN:istepmin)
C$OMP& reduction(MIN:istepmingasnew)
C$OMP& reduction(+:istepeta,etafac)
c
c--Compute new time steps for particles which have been advanced
c
C$OMP DO SCHEDULE(runtime)
      DO j = 1, nlst
         i = llist(j)
         rmodall(:) = 1.0E+30
c
c--Compute timestep from force condition
c
         force2 = fxyzu(1,i)**2 + fxyzu(2,i)**2 + fxyzu(3,i)**2
         IF (force2.NE.0.0) THEN
            IF (iphase(i).GT.0 .AND. iphase(i).LT.10 .AND. 
     &           iptsoft.EQ.1) THEN
               softrad = ptsoft
            ELSE
               softrad = xyzmh(5,i)
            ENDIF
            rmod = Courant*SQRT(softrad/SQRT(force2))
c
c--If gas is involved, sink timesteps divided by 1000 (pre-existing),
c  regardless of whether planetesimals exist. If the case is just
c  planetesimals (iplanetesimals=1), then the sink can be on a
c  longer timestep, but all particle timesteps reduced according to
c  the Courant number. The planetesimals must always be reduced by the
c  Courant factor.
c
            IF (iplanetesimals.NE.1) THEN
               IF (iphase(i).GT.0 .AND. iphase(i).LT.10) rmod =
     &              rmod/10.0
               IF (iphase(i).EQ.11) rmod = rmod*Courant
            ELSEIF (iplanetesimals.EQ.1) THEN
               rmod = rmod*Courant
            ENDIF
            rmod = rmod/(dt*isteps(i)/imaxstep)
         ELSE
            rmod = 1.0E+30
         ENDIF
         rmodall(5) = rmod
c
c--For sink particles, also use velocity criterion
c
c         IF (iphase(i).GT.0 .AND. iphase(i).LT.10) THEN
c            vel2 = vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2
c            IF (vel2.NE.0.0) THEN
c               rmodvel = 0.1*xyzmh(5,i)/SQRT(vel2)
cc               write (*,*) 'h ',xyzmh(5,i),rmodvel
c               rmodvel = rmodvel/(dt*isteps(i)/imaxstep)
c            ELSE
c               rmodvel = 1.0E+30
c            ENDIF
c            rmod = MIN(rmod, rmodvel)
c         ENDIF
c
c--Compute the time step for the gas from the Courant condition
c
         IF (iphase(i).EQ.0) THEN
            IF (rho(i).GT.0.0) THEN
               divvi = divv(i)/rho(i)
            ELSE
               divvi = 0.0
            ENDIF
            IF (ifsvi.GE.6) THEN
               IF (imhd.EQ.idim) THEN
                  valfven2 = (Bxyz(1,i)*Bxyz(1,i) +  Bxyz(2,i)*Bxyz(2,i)
     &                      + Bxyz(3,i)*Bxyz(3,i))/rho(i)
                  aux1 = max(alphaMM(1,i),alphaMM(isizealphaMM,i))*
     &                   sqrt(vsound(i)**2 + valfven2)
               ELSE
                  valfven2 = 0.
                  aux1 = alphaMM(1,i)*vsound(i)
               ENDIF
               aux2 = xyzmh(5,i)*ABS(divvi)
               aux3 = aux2
               IF (divvi.LT.0.0) THEN
                  aux2 = 2.0*alphaMM(1,i)*aux2
               ELSE
                  aux2 = 0.0
               ENDIF
            ELSE
               IF (imhd.EQ.idim) THEN
                   xi = xyzmh(1,i)
                   yi = xyzmh(2,i)
                   zi = xyzmh(3,i)
                   Bxi = Bxyz(1,i)
                   Byi = Bxyz(2,i)
                   Bzi = Bxyz(3,i)
                   IF (iBext.GT.0 .AND. iexf.NE.9) THEN
                      Bxi = Bxi + BEXTERNAL(xi,yi,zi,1)
                      Byi = Byi + BEXTERNAL(xi,yi,zi,2)
                      Bzi = Bzi + BEXTERNAL(xi,yi,zi,3)
                   ENDIF
                   valfven2 = (Bxi*Bxi +  Byi*Byi + Bzi*Bzi)/rho(i)
                   
c                  valfven2 = (Bxyz(1,i)*Bxyz(1,i) +  Bxyz(2,i)*Bxyz(2,i)
c     &                      + Bxyz(3,i)*Bxyz(3,i))/rho(i)
                  aux1 = alpha*sqrt(vsound(i)**2 + valfven2)
               ELSE
                  valfven2 = 0.
                  aux1 = alpha*vsound(i)
               ENDIF
               aux2 = xyzmh(5,i)*ABS(divvi)
               aux3 = aux2
               IF (divvi.LT.0.0) THEN
                  aux2 = beta*aux2
               ELSE
                  aux2 = 0.0
               ENDIF
            ENDIF
            denom = aux3 + sqrt(vsound(i)**2 + valfven2)
     &                   + 1.2*(aux1 + aux2)
            crstepi = Courant*xyzmh(5,i)/denom
            rmodcr = crstepi/(dt*isteps(i)/imaxstep)
            rmod = MIN(rmod, rmodcr)
            rmodall(6) = rmodcr
         ELSEIF (iphase(i).GE.10) THEN
            IF (iphase(i).EQ.11 .AND. gasdrag) THEN
               rmodgd = 0.1*planetesimaltimestep(i)/
     &              (dt*isteps(i)/imaxstep)
            ELSE
               rmodgd = 1.0E+30
            ENDIF
            rmod = MIN(rmod, rmodgd)
            rmodall(8) = rmodgd
         ENDIF
c
c--Timestep from change in thermal energy
c
         IF (iradtrans.EQ.idim .AND. encal.EQ.'r') THEN
            IF (fxyzu(4,i).GT.tiny) THEN
               rmod_du = 0.1*ABS(vxyzu(4,i)/fxyzu(4,i))/
     &              (dt*isteps(i)/imaxstep)
               rmod = MIN(rmod,rmod_du)
               rmodall(9) = rmod_du
            ENDIF
         ENDIF
c
c--New timestep for cooling equation of state
c
         IF(iphase(i).EQ.0 .AND. encal.EQ.'c') THEN
             rad = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 + 
     &               xyzmh(3,i)**2 )
             omega = sqrt(xyzmh(4,listpm(1))/ rad**3)
             omega1 = 1./omega
             tcool = 5.0 * omega1
             coolstepi = 0.01 * tcool
             rmodcool = coolstepi/(dt*isteps(i)/imaxstep)
             rmod = MIN(rmod,rmodcool)
             rmodall(13) = rmodcool
         END IF
c
c--New timestep if particles cross boundaries
c
         IF(iphase(i).EQ.0 .AND. iexf.GT.0 .OR. iBext.GT. 0) THEN
             dti = rmod*(dt*isteps(i)/imaxstep)
             rmod = MIN(rmod,dtcrossbound(xyzmh(1,i),xyzmh(2,i),
     &              xyzmh(3,i),vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              fxyzu(1,i),fxyzu(2,i),fxyzu(3,i),dti)/
     &              (dt*isteps(i)/imaxstep))   
         END IF
c
c--Check resistive diffusion time
c
         IF (imhd.EQ.idim) THEN
            IF (iphase(i).EQ.0) THEN
               IF (iresist.GE.2) THEN
c           resistivity which is a function of density
                  etai = etafunc(real(rho(i)),vxyzu(4,i),etamhd,iresist)
                  IF (etai.GT.tiny) THEN
                     dteta = 0.15*xyzmh(5,i)**2/etai
                     rmodeta = dteta/(dt*isteps(i)/imaxstep)
                     IF (rmodeta.LT.rmod) THEN
                        istepeta = istepeta + 1
                        etafac = etafac 
     &                       + rmod*(dt*isteps(i)/imaxstep)/dteta
                        rmod = MIN(rmod, rmodeta)
                        rmodall(14) = rmodeta
                     ENDIF
                  ENDIF
               ELSEIF (iresist.EQ.1) THEN
c           constant resistivity
                  IF (etamhd.GT.tiny) THEN
                     dteta = 0.15*xyzmh(5,i)**2/etamhd
                     rmodeta = dteta/(dt*isteps(i)/imaxstep)
                     IF (rmodeta.LT.rmod) THEN
                        istepeta = istepeta + 1
                        etafac = etafac 
     &                       + rmod*(dt*isteps(i)/imaxstep)/dteta
                        rmod = MIN(rmod, rmodeta)
                        rmodall(14) = rmodeta
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
c
c--Record limiting criterion for timestep
c
         timestepcriterion(i) = MINLOC(rmodall,DIM=1)
c
c--Reduce time step (no synchronization needed)
c
         IF (rmod.LT.0.97) THEN
            stepi = dt*isteps(i)/imaxstep*rmod
            ibin = INT(LOG10(dt/stepi)/xlog2) + 1

            IF (ibin.LT.0 .OR. ibin.GT.30) THEN
               print *,'ERROR - ibin.LT.0 .OR. ibin.GT.30 ',ibin,stepi,
     &              dt,isteps(i),i,rmod,rmodcr,rmodgd,xyzmh(5,i),
     &              vsound(i),alphaMM(1,i),denom,crstepi,iphase(i),
     &              force2,timestepcriterion(i)
               IF (ibin.LT.0) THEN
                  CALL quit(1)
               ELSE
                  ibin = 30
               ENDIF
c
c--Allow to break time step criterion based on change of thermal energy
c
               IF (timestepcriterion(i).NE.9) CALL quit(1)
            ENDIF

            isteps(i) = imaxstep/(2**ibin)
c
c--Increase time step, if synchronization allows
c
         ELSEIF (rmod.GE.2.0) THEN
            istep2 = 2*isteps(i)
            irat = MOD(idtsyn, istep2)
            IF (irat.EQ.0) THEN
               isteps(i) = istep2
            ENDIF
         ENDIF
         IF (isteps(i).GT.imaxstep) isteps(i) = imaxstep
c
c--Report if time steps are too small
c
         IF (isteps(i).LT.2) THEN
            IF (iphase(i).GE.1 .AND. iphase(i).LT.10) THEN
C$OMP CRITICAL (writeiprint)
               WRITE (iprint, 99301)
99301          FORMAT ('STEP : Sink allowed to break TIme step')

               isteps(i) = 2

               WRITE (iprint,*)i, iorig(i), listrealpm(i), stepi, ibin
               WRITE (iprint,*)rmod, nlst

               xmindist = 1.0e+30
               DO l = 1, nptmass
                  iptcur = listpm(l)
                  IF (iptcur.NE.i) THEN
                     rad2 = (xyzmh(1,i) - xyzmh(1,iptcur))**2 + 
     &                    (xyzmh(2,i) - xyzmh(2,iptcur))**2 + 
     &                    (xyzmh(3,i) - xyzmh(3,iptcur))**2
                     xmindist = MIN(xmindist,rad2)
                  ENDIF
c                  WRITE (iprint,*)'sink ',l,listpm(l),nptlist(l),rad2
c                  WRITE (iprint,*)iptmass,hacc,haccall,ptmcrit,radcrit,
c     &                 xyzmh(5,iptcur)
c                  WRITE (iprint,*)
               END DO
               WRITE (iprint,*)'rmin ',xmindist
C$OMP END CRITICAL (writeiprint)
            ELSE
c                        write (13,*) isteps(i), 'print4'
C$OMP CRITICAL (writeiprint)
               WRITE (iprint, 99300)
99300          FORMAT ('STEP : Step too small ! Nothing can help!')
               WRITE (iprint,*)i, iorig(i), iphase(i), xyzmh(1,i), 
     &              xyzmh(2,i), xyzmh(3,i)
               WRITE (iprint,*)vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              vxyzu(4,i),xyzmh(5,i)
               WRITE (iprint,*)fxyzu(1,i),fxyzu(2,i),fxyzu(3,i)
               IF (imhd.EQ.idim) THEN
                  WRITE (iprint,*) 'B=',Bxyz(1,i),Bxyz(2,i),Bxyz(3,i)
                  WRITE (iprint,*) 'valfven=',valfven2,aux1
               ENDIF
               WRITE (iprint,*)rmod,rmodcr,timestepcriterion(i)
               WRITE (iprint,*)(rmodall(l),l=1,15)
               WRITE (iprint,*)nneigh(i),rho(i),nlst

               xmindist = 1.0e+30
               DO l = 1, nptmass
                  iptcur = listpm(l)
                  IF (iptcur.NE.i) THEN
                     rad2 = (xyzmh(1,i) - xyzmh(1,iptcur))**2 + 
     &                    (xyzmh(2,i) - xyzmh(2,iptcur))**2 + 
     &                    (xyzmh(3,i) - xyzmh(3,iptcur))**2
                     xmindist = MIN(xmindist,rad2)
                  WRITE (iprint,*)'sink ',l,listpm(l),nptlist(l),rad2
                  WRITE (iprint,*)iptmass,hacc,haccall,ptmcrit,radcrit,
     &                    xyzmh(5,iptcur)
                     DO ll = 1, nptlist(l)
                        WRITE (iprint,*) '  ',ll,nearpt(ll,l)
                     END DO
                     WRITE (iprint,*)
                  ENDIF
               END DO
               WRITE (iprint,*)'rmin ',xmindist

               IF (timestepcriterion(i).EQ.9) THEN
c--Limit due to criterion for dudt ignored
                  isteps(i) = 2
               ELSE
                  CALL quit(1)
               ENDIF
C$OMP END CRITICAL (writeiprint)
            ENDIF
         ENDIF
c
c--Update minimum time step
c
         istepmax = MAX(istepmax, isteps(i))
         istepmin = MIN(istepmin, isteps(i))
         IF(iphase(i).EQ.0) istepmingasnew=MIN(istepmingasnew,isteps(i))
      END DO
C$OMP END DO
C$OMP END PARALLEL

#ifdef MPI
      CALL MPI_ALLREDUCE(istepmax,istepmaxtot,1,MPI_INTEGER,MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(istepmin,istepmintot,1,MPI_INTEGER,MPI_MIN,
     &     MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(istepmingasnew,istepmingastot,1,MPI_INTEGER,
     &     MPI_MIN,MPI_COMM_WORLD, ierr)
      istepmax = istepmaxtot
      istepmin = istepmintot
      istepmingasnew = istepmingastot
#endif

      IF (imhd.eq.idim .and. istepeta.GT.1) THEN
         WRITE(*,*) 'WARNING: resistivity constraining dt on ',
     &            istepeta,' parts by ratio ',etafac/real(istepeta)
      ENDIF

      IF (istepmingasnew.NE.imax) THEN
         istepmingas = istepmingasnew
         IF (istepmingas.LT.2) THEN
            WRITE (iprint,*) 'ERROR - istepmingas.LT.2'
            CALL quit(1)
         ENDIF
      ENDIF

      RETURN

      END
