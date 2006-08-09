       SUBROUTINE timestep (dt,idtsyn,nlst,llist,fxyzu,Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine calculates the new timestep for each     *
c     particle in the list.                                 *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      PARAMETER (Courant=0.30)

      DIMENSION llist(idim),fxyzu(4,idim)
      DIMENSION Bevolxyz(3,imhd)

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

      xlog2 = 0.30103
      istepmin = imax
      istepmingasnew = imax
      istepmax = -1

C$OMP PARALLEL default(none)
C$OMP& shared(nlst,llist)
C$OMP& shared(xyzmh,vxyzu)
C$OMP& shared(fxyzu,Bxyz)
C$OMP& shared(dt,isteps,imaxstep,divv,rho)
C$OMP& shared(alpha,beta,vsound,xlog2)
C$OMP& shared(idtsyn,nearpt,nptlist)
C$OMP& shared(iprint,nneigh,iphase)
C$OMP& shared(iorig,nptmass,listpm,listrealpm,encal)
C$OMP& shared(ifsvi,alphaMM,iptsoft,ptsoft)
C$OMP& shared(iptmass,hacc,haccall,ptmcrit,radcrit)
C$OMP& private(i,j,l,ll,xmindist,iptcur,rad2)
C$OMP& private(divvi,aux1,aux2,aux3,denom,valfven2,rmodcool)
C$OMP& private(crstepi,rmodcr,rmodvel,force2,rmod,vel2,coolstepi)
C$OMP& private(stepi,ibin,istep2,irat,softrad,rad,omega,omega1,tcool)
C$OMP& reduction(MAX:istepmax)
C$OMP& reduction(MIN:istepmin)
C$OMP& reduction(MIN:istepmingasnew)
c
c--Compute new time steps for particles which have been advanced
c
C$OMP DO SCHEDULE(runtime)
      DO j = 1, nlst
         i = llist(j)
c
c--Compute timestep from force condition
c
         force2 = fxyzu(1,i)**2 + fxyzu(2,i)**2 + fxyzu(3,i)**2
         IF (force2.NE.0.0) THEN
            IF (iphase(i).GT.0 .AND. iptsoft.EQ.1) THEN
               softrad = ptsoft
            ELSE
               softrad = xyzmh(5,i)
            ENDIF
            rmod = Courant*SQRT(softrad/SQRT(force2))
            IF (iphase(i).GT.0) rmod = rmod/1000.0
            rmod = rmod/(dt*isteps(i)/imaxstep)
         ELSE
            rmod = 1.0E+30
         ENDIF
c
c--For sink particles, also use velocity criterion
c
c         IF (iphase(i).GT.0) THEN
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
            IF (ifsvi.EQ.6) THEN
               IF (imhd.EQ.idim) THEN
                  valfven2 = (Bxyz(1,i)*Bxyz(1,i) +  Bxyz(2,i)*Bxyz(2,i)
     &                      + Bxyz(3,i)*Bxyz(3,i))/rho(i)
                  aux1 = max(alphaMM(1,i),alphaMM(2,i))*
     &                   sqrt(vsound(i)**2 + valfven2)
               ELSE
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
                  valfven2 = (Bxyz(1,i)*Bxyz(1,i) +  Bxyz(2,i)*Bxyz(2,i)
     &                      + Bxyz(3,i)*Bxyz(3,i))/rho(i)
                  aux1 = alpha*sqrt(vsound(i)**2 + valfven2)
               ELSE
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
            denom = aux3 + vsound(i) + 1.2*(aux1 + aux2)
            crstepi = Courant*xyzmh(5,i)/denom
            rmodcr = crstepi/(dt*isteps(i)/imaxstep)
            rmod = MIN(rmod, rmodcr)
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
c             PRINT *,rmod,rmodcool
             rmod = MIN(rmod,rmodcool)
         END IF
c
c--Reduce time step (no synchronization needed)
c
         IF (rmod.LT.0.97) THEN
            stepi = dt*isteps(i)/imaxstep*rmod
            ibin = INT(LOG10(dt/stepi)/xlog2) + 1
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
            IF (iphase(i).GE.1) THEN
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
C$OMP CRITICAL (writeiprint)
               WRITE (iprint, 99300)
99300          FORMAT ('STEP : Step too small ! Nothing can help!')
               WRITE (iprint,*)i, iorig(i), iphase(i), xyzmh(1,i), 
     &              xyzmh(2,i), xyzmh(3,i)
               WRITE (iprint,*)vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &              vxyzu(4,i),xyzmh(5,i)
               WRITE (iprint,*)fxyzu(1,i),fxyzu(2,i),fxyzu(3,i)
               IF (imhd.EQ.idim) THEN
                  WRITE (iprint,*) Bxyz(1,i),Bxyz(2,i),Bxyz(3,i)
               ENDIF
               WRITE (iprint,*)rmod,rmodcr
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

               CALL quit
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

      istepmingas = istepmingasnew
      IF (istepmingas.LE.2) THEN
         WRITE (iprint,*) 'ERROR - istepmingas.LE.2'
         CALL quit
      ENDIF

      RETURN

      END
