      SUBROUTINE divBdiffuse(dt,nlst_in,nlst_end,npart,list,
     &     xyzmh,rho,Bxyz)

      INCLUDE 'idim'
      
      DIMENSION xyzmh(5,idim)
      DIMENSION Bxyz(3,imhd)
      REAL*4 rho(idim)
      DIMENSION list(idim)
      PARAMETER (tolerance=1.e-3)

      INCLUDE 'COMMONS/logun'
      
      CHARACTER*7 where
      LOGICAL*1 moresweep

      DATA where/'divBdiffuse'/
      WRITE (*,*) 'ASS, ',nlst_end
c     &     ,xyzmh(5,1),xyzmh(5,2),xyzmh(5,3),xyzmh(5,4)

c
c--Loop until has arrived at a solution to the implicit diffusion equation
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
            CALL divBiterate(increment*dt/Nstep,nlst_in,nlst_end,
     &           npart,list,xyzmh,rho,Bxyz,
     &           moresweep,nit,error)
      
            IF (moresweep) GOTO 100
            IF (ifrac.EQ.Nstep) GOTO 150
c
c--sometimes needs a small dt to start with 
c  but can then double dt if it works fine
c            
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
            WRITE (iprint,*) 'ERROR in divBdiffuse: Nsteps > 10^6'
            CALL quit
         ENDIF

      END DO
      
 150  CONTINUE

c      PRINT *,"SEGMENT: Call complete with ",Nstep," subdivisions"

      RETURN

      END
