      SUBROUTINE ASS(nlst_in,nlst_end,list,dt,itime,
     $     npart,ntot,xyzmh,vxyzu,ekcle,trho,dedxyz,alphaMM)

      INCLUDE 'idim'
      
      DIMENSION xyzmh(5,mmax2)
      DIMENSION vxyzu(4,idim2)
      DIMENSION ekcle(5,iradtrans2)
      REAL*4 trho(idim2),alphaMM(isizealphaMM,idim2)
      DIMENSION list(idim2)
      DIMENSION dedxyz(3,iradtrans2)

      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/radstore'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/implicit'

      DIMENSION EUsave(2,iradtrans2)

      LOGICAL moresweep

c      WRITE (*,*) 'ASS, icall ',icall,nlst_end
c
c--Keep original versions of U and E and set cv and kappa
c
      istepmindone = imax
C$OMP PARALLEL DO SCHEDULE (static) default(none)
C$OMP& shared(list,ekcle,origEU,vxyzu,trho,nlst_in,nlst_end)
C$OMP& shared(icall,iphase,ifullstep,nlst0,EUsave,isteps)
C$OMP& private(n,i)
C$OMP& reduction(MIN:istepmindone)
      DO n = nlst_in, nlst_end
         i = list(n)
         IF (n.LE.nlst0) THEN
            istepmindone = MIN(istepmindone,isteps(i))
            ifullstep(i) = .TRUE.
         ELSE
            ifullstep(i) = .FALSE.
         ENDIF

         IF (iphase(i).EQ.0) THEN
            origEU(1,i) = ekcle(1,i)
            origEU(2,i) = vxyzu(4,i)

            EUsave(1,i) = ekcle(1,i)
            EUsave(2,i) = vxyzu(4,i)
         ENDIF
      END DO
C$OMP END PARALLEL DO
      IF (istepmindone.EQ.imax) istepmindone = 0
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

c            PRINT *,"DOSEG: ",K,Nstep,ifrac
c     &           ,increment*dt/Nstep,nlst_in,nlst_end,npart,
c     &     list(1),ekcle(1,1),ekcle(2,1),xyzmh(5,1),vxyzu(1,1)

c            IF (icall.EQ.1) THEN
c            CALL GSIMPL(increment*dt/Nstep,dt,itime,
c     &           npart,ntot,ekcle,xyzmh,vxyzu,dedxyz,trho,
c     &           alphaMM,moresweep,nit,error,istepmindone)
c            ELSE
            CALL GSIMPLS(increment*dt/Nstep,dt,itime,
     &           npart,ntot,ekcle,xyzmh,vxyzu,dedxyz,trho,
     &           alphaMM,moresweep,nit,error,istepmindone)
c            ENDIF

c            IF (moresweep) GOTO 100

            IF (moresweep) THEN
         print *,'Integration failed - using U & E values anyway ',
     &              nlst_end
       WRITE(iprint,*) 'Integration failed - using U & E values anyway',
     &        nlst_end
               GOTO 150
            ENDIF


            IF (ifrac.EQ.Nstep) GOTO 150
c            IF (nit.EQ.1 .AND. MOD(ifrac,2*increment).EQ.0
c     &         .AND.error.LT.0.5*tolerance/increment**2
c     &         .AND.increment.EQ.1) THEN
c               increment = increment * 2
c            ENDIF

C$OMP PARALLEL DO SCHEDULE (static) default(none)
C$OMP& shared(list,ekcle,origEU,vxyzu,trho,nlst_in,nlst_end)
C$OMP& shared(iphase)
C$OMP& private(n,i)
            DO n = nlst_in, nlst_end
               i = list(n)
               IF (iphase(i).EQ.0) THEN
                  origEU(1,i) = ekcle(1,i)
                  origEU(2,i) = vxyzu(4,i)

                  ekcle(3,i) = GETCV(trho(i),vxyzu(4,i))
                  ekcle(2,i) = GETKAPPA(vxyzu(4,i),ekcle(3,i),trho(i))
               ENDIF
            END DO
C$OMP END PARALLEL DO
            

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
C$OMP& shared(list,ekcle,origEU,vxyzu,trho,nlst_in,nlst_end)
C$OMP& shared(iphase,EUsave)
C$OMP& private(n,i)
         DO n = nlst_in, nlst_end
            i = list(n)
            IF (iphase(i).EQ.0) THEN
               ekcle(1,i) = EUsave(1,i)
               vxyzu(4,i) = EUsave(2,i)

               ekcle(3,i) = GETCV(trho(i),vxyzu(4,i))
               ekcle(2,i) = GETKAPPA(vxyzu(4,i),ekcle(3,i),trho(i))
            ENDIF
         END DO
C$OMP END PARALLEL DO

c         print *,'Integration failed - kept old U & E values',nlst_end
c         WRITE (iprint,*) 'Integration failed - kept old U & E values',
c     &        nlst_end
c         GOTO 150
      END DO
      
 150  CONTINUE

c      PRINT *,"SEGMENT: Call complete with ",Nstep," subdivisions"

      RETURN

      END
