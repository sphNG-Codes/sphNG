      SUBROUTINE ASS_DUST1(nlst_in,nlst_end,list,dt,itime,npart,ntot,
     &       xyzmh,dustvar,rho)
    
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/eosq'
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

      REAL, DIMENSION(idim_dustFluid) :: dustvar
      REAL, DIMENSION(idim_dustFluid) :: origdustvar
      REAL, DIMENSION(idim_dustFluid) :: dustvarsave
      REAL, DIMENSION(5,mmax2) :: xyzmh
      REAL, DIMENSION(4,idim2) :: vxyzu
      
      REAL*4, DIMENSION(idim2) :: rho

      INTEGER, DIMENSION(idim2) :: list

      REAL dt
      REAL itime
      REAL nit
      REAL errorS

      INTEGER nlst_in, nlst_end
      INTEGER n,i,k
      INTEGER nlst0
      INTEGER increment
      INTEGER ifrac
      INTEGER Nstep
      INTEGER npart, ntot
      INTEGER iprint

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
C$OMP& shared(list,nlst_in,nlst_end)
C$OMP& shared(dustvar,origdustvar,dustvarsave,nlst0)
C$OMP& shared(iphase)
C$OMP& private(i,n)
C$OMP DO SCHEDULE(static)
      DO n = nlst_in, nlst_end
         i = list(n)

         IF (iphase(i).EQ.0) THEN
            origdustvar(i) = dustvar(i)
            dustvarsave(i) = dustvar(i)
         ENDIF ! iphase(i).EQ.0
      END DO ! n = nlst_in, nlst_end
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
            CALL dust1_implicit(increment*dt/Nstep,itime,npart,ntot,
     &            xyzmh,rho,dustvar,nit,errorS,moresweep)
            IF (moresweep) THEN
               WRITE (iprint,38000) nlst_end, errorS
38000          FORMAT('Integration failed - increasing Nstep',I8,1X,
     &         1PE11.4)
               GOTO 300
            ENDIF ! moresweep

            IF (ifrac.EQ.Nstep) GOTO 150
C$OMP PARALLEL default(none)
C$OMP& shared(list,nlst_in,nlst_end)
C$OMP& shared(dustvar,origdustvar)
C$OMP& shared(iphase)
C$OMP& private(i,n)
C$OMP DO SCHEDULE (static)
            DO n = nlst_in, nlst_end
               i = list(n)
               IF (iphase(i).EQ.0) THEN
                  origdustvar(i) = dustvar(i)
               ENDIF ! (iphase(i).EQ.0)
            ENDDO ! n = nlst_in, nlst_end
C$OMP END DO
C$OMP END PARALLEL
         END DO ! k = 1, Nstep 
c
c--Integration failed, increase number of sub-steps and restore
c     initial state
c
300      CONTINUE
         Nstep = Nstep * 2
c         PRINT *, 'ASS_DUST1 Nstep doubling = ',Nstep
         IF (Nstep.GT.1000000) THEN
            WRITE (iprint,*) 'ERROR in 1fluid subscycle: Nsteps > 10^6'
            CALL quit(1)
         ENDIF ! Nstep.GT.1000000
C$OMP PARALLEL default(none)
C$OMP& shared(list,nlst_in,nlst_end)
C$OMP& shared(dustvar,dustvarsave)
C$OMP& shared(iphase)
C$OMP& private(i,n)
C$OMP DO SCHEDULE (static)
            DO n = nlst_in, nlst_end
               i = list(n)
               IF (iphase(i).EQ.0) THEN
                  dustvar(i) = dustvarsave(i)
               ENDIF ! (iphase(i).EQ.0)
            ENDDO ! n = nlst_in, nlst_end
C$OMP END DO
C$OMP END PARALLEL
      END DO ! This is a do while loop

150   CONTINUE

c      PRINT *, 'ASS_DUST1 Nstep = ', Nstep

      RETURN

      END
