      SUBROUTINE ASS(nlst_in,nlst_end,nlstall,list,dt,itime,npart,
     $     xyzmh,vxyzu,ekcle,trho,dedxyz)

      INCLUDE 'idim'
      
      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(5,iradtrans)
      REAL*4 trho(idim)
      DIMENSION list(idim)
      DIMENSION dedxyz(3,iradtrans)

      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/logun'

      DIMENSION saveuekcdedxyz(7,idim)
      
      CHARACTER*7 where
      LOGICAL*1 moresweep

      DATA where/'ass'/
      WRITE (*,*) 'ASS, icall ',icall,nlstall,nlst_end
c     &     ,xyzmh(5,1),xyzmh(5,2),xyzmh(5,3),xyzmh(5,4)
c
c--Keep original versions of U and E and set cv and kappa
c
C$OMP PARALLEL DO SCHEDULE (static) default(none)
C$OMP& shared(list,ekcle,saveuekcdedxyz,vxyzu,trho,nlst_in,nlstall)
C$OMP& shared(dedxyz)
C$OMP& private(n,i)
      DO n = nlst_in, nlstall
         i = list(n)

         ekcle(3,i) = GETCV(trho(i),vxyzu(4,i))
         ekcle(2,i) = GETKAPPA(vxyzu(4,i),ekcle(3,i),trho(i))

         saveuekcdedxyz(1,i) = vxyzu(4,i)
         saveuekcdedxyz(2,i) = ekcle(1,i)
         saveuekcdedxyz(3,i) = ekcle(2,i)
         saveuekcdedxyz(4,i) = ekcle(3,i)
         saveuekcdedxyz(5,i) = dedxyz(1,i)
         saveuekcdedxyz(6,i) = dedxyz(2,i)
         saveuekcdedxyz(7,i) = dedxyz(3,i)
      END DO
C$OMP END PARALLEL DO
c
c--Loop until has arrived at a solution to the implicit radiative transfer
c      
      Nstep = 1
      DO
         moresweep = .FALSE.
c
c--Split required timestep up into Nstep steps and integrate energies
c
         increment = 1
         ifrac = 0
         DO k = 1, Nstep
            ifrac = ifrac + increment

            PRINT *,"DOSEG: ",K,Nstep,ifrac
c     &           ,increment*dt/Nstep,nlst_in,nlst_end,nlstall,npart,
c     &     list(1),ekcle(1,1),ekcle(2,1),xyzmh(5,1),vxyzu(1,1)
            CALL GSIMPL(increment*dt/Nstep,nlst_in,nlst_end,
     &           nlstall,npart,list,ekcle,xyzmh,vxyzu,dedxyz,trho,
     &           moresweep,nit,error)
            
            IF (moresweep) GOTO 100
            IF (ifrac.EQ.Nstep) GOTO 150
            IF (nit.EQ.1 .AND. MOD(ifrac,2*increment).EQ.0
     &         .AND.error.LT.0.5*tolerance/increment**2
     &         .AND.increment.EQ.1) THEN
               increment = increment * 2
            ENDIF
         END DO
c
c--Integration failed, increase number of sub-steps and restore initial state
c                  
 100     Nstep = Nstep * 2
         IF (Nstep.GT.1000000) THEN
            WRITE (iprint,*) 'ERROR in ASS: Nsteps > 10^6'
            CALL quit
         ENDIF
C$OMP PARALLEL DO SCHEDULE (static) default(none)
C$OMP& shared(list,ekcle,saveuekcdedxyz,vxyzu,trho,nlst_in,nlstall)
C$OMP& shared(dedxyz)
C$OMP& private(n,i)
         DO n = nlst_in, nlstall
            i = list(n)
            vxyzu(4,i) = saveuekcdedxyz(1,i)
            ekcle(1,i) = saveuekcdedxyz(2,i)
            ekcle(2,i) = saveuekcdedxyz(3,i)
            ekcle(3,i) = saveuekcdedxyz(4,i)
            dedxyz(1,i) = saveuekcdedxyz(5,i)
            dedxyz(2,i) = saveuekcdedxyz(6,i)
            dedxyz(3,i) = saveuekcdedxyz(7,i)
         END DO
C$OMP END PARALLEL DO
      END DO
      
 150  CONTINUE

c      PRINT *,"SEGMENT: Call complete with ",Nstep," subdivisions"

      RETURN

      END
