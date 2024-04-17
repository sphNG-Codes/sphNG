      SUBROUTINE ASS_DUST1(dt,itime,npart,ntot,xyzmh,dustvar,rho)
    
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/fullstep'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/dustimplstat'

      REAL, DIMENSION(ndusttypes,idim_dustFluid2) :: dustvar
      REAL, DIMENSION(ndusttypes,idim_dustFluid2) :: dustvarsave
      REAL, DIMENSION(5,mmax2) :: xyzmh
      REAL, DIMENSION(4,idim2) :: vxyzu
      
      REAL*4, DIMENSION(idim2) :: rho

      REAL dt, dtdo
      REAL errorS

      INTEGER n,i,k
      INTEGER increment
      INTEGER ifrac
      INTEGER Nstep
      INTEGER npart, ntot
      INTEGER iprint
      INTEGER nit
      INTEGER itime

      LOGICAL moresweep
c
c--Need a coefficient to use gas sound speed to calculate the thermal
c     velocity of the gas, i.e. vthermal = SQRT(8 kB T / pi mu mH).
c     Do that here to avoid doing SQRT() A LOT of times later on
c
      IF (encal.EQ.'i'.OR.encal.EQ.'t') THEN
         vthermal_coef = SQRT(8.0/pi)
      ELSE
         vthermal_coef = SQRT(8.0/pi/gamma)
      ENDIF
c
c--Keep original version of dustvar
c
C$OMP PARALLEL default(none)      
C$OMP& shared(ivar, ncompact)
C$OMP& shared(dustvar,dustvarsave)
C$OMP& shared(iphase)
C$OMP& private(i,n)
C$OMP DO SCHEDULE(static)
      DO n = 1, ncompact
         i = ivar(3,n)

         IF (iphase(i).EQ.0) THEN
            dustvarsave(:,i) = dustvar(:,i)
         ENDIF ! iphase(i).EQ.0
      END DO ! n = 1, ncompact
C$OMP END DO
C$OMP END PARALLEL
c
c--Loop until has arrived at a solution to the implicit 1fluid dust
c
      Nstep = 1
      DO
         moresweep = .FALSE.
c
c--Split required timestep up into Nstep steps and integrate dustvar
c
         increment = 1
         ifrac = 0
         DO k = 1, Nstep
            ifrac = ifrac + increment
            dtdo = increment*dt/Nstep
            CALL dust1_implicit(dtdo,itime,npart,ntot,
     &            xyzmh,rho,dustvar,nit,errorS,moresweep)

c            print *,'After dust1_implicit ',itime,ifrac,Nstep,moresweep
            IF (moresweep) THEN
               WRITE (iprint,38000) Nstep, ncompact, errorS
38000          FORMAT('Integration failed - increasing Nstep',2(I8,1X),
     &         1PE10.3)
               ndustimpl_fail = ndustimpl_fail + 1
               GOTO 300
            ENDIF ! moresweep

            IF (ifrac.EQ.Nstep) GOTO 150
         END DO ! k = 1, Nstep 
c
c--Integration failed, increase number of sub-steps and restore
c     initial state
c
300      CONTINUE
         Nstep = Nstep * 2
         PRINT *, 'ASS_DUST1 Nstep doubling = ',Nstep,' at ',itime
         IF (Nstep.GT.1000000) THEN
            WRITE (iprint,*) 'ERROR in 1fluid subscycle: Nsteps > 10^6'
            CALL quit(1)
         ENDIF ! Nstep.GT.1000000
C$OMP PARALLEL default(none)
C$OMP& shared(ivar,ncompact)
C$OMP& shared(dustvar,dustvarsave)
C$OMP& shared(iphase)
C$OMP& private(i,n)
C$OMP DO SCHEDULE (static)
            DO n = 1, ncompact
               i = ivar(3,n)

               IF (iphase(i).EQ.0) THEN
                  dustvar(:,i) = dustvarsave(:,i)
               ENDIF ! (iphase(i).EQ.0)
            ENDDO ! n = 1, ncompact
C$OMP END DO
C$OMP END PARALLEL
      END DO ! This is a do while loop

150   CONTINUE

c      PRINT *, 'ASS_DUST1 Nstep = ', Nstep

      RETURN

      END
