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

      xlog2 = 0.30103
      istepmin = imax
      istepmingasnew = imax
      istepmax = -1

C$OMP PARALLEL default(none)
C$OMP& shared(nlst,llist)
C$OMP& shared(xyzmh,vxyzu)
C$OMP& shared(fxyzu)
C$OMP& shared(dt,isteps,imaxstep,divv,rho)
C$OMP& shared(alpha,beta,vsound,xlog2)
C$OMP& shared(idtsyn)
C$OMP& shared(iprint,nneigh,iphase)
C$OMP& shared(iorig)
C$OMP& shared(ifsvi,alphaMM,iptsoft,ptsoft)
C$OMP& private(i,j)
C$OMP& private(divvi,aux1,aux2,aux3,denom)
C$OMP& private(crstepi,rmodcr,rmodvel,force2,rmod,vel2)
C$OMP& private(stepi,ibin,istep2,irat,softrad)
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
            crstepi = Courant*xyzmh(5,i)/denom
            rmodcr = crstepi/(dt*isteps(i)/imaxstep)
            rmod = MIN(rmod, rmodcr)
         ENDIF
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
C$OMP CRITICAL (writeiprint)
            WRITE (iprint, 99300)
99300       FORMAT ('STEP : Step too small ! Nothing can help!')
            WRITE (iprint,*)i, iorig(i), iphase(i), xyzmh(1,i), 
     &           xyzmh(2,i), xyzmh(3,i)
            WRITE (iprint,*)vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),
     &           vxyzu(4,i),xyzmh(5,i)
            WRITE (iprint,*)fxyzu(1,i),fxyzu(2,i),fxyzu(3,i)
            WRITE (iprint,*)rmod,rmodcr
            WRITE (iprint,*)nneigh(i),rho(i)
            CALL quit
C$OMP END CRITICAL (writeiprint)
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

      RETURN

      END
