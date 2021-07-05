      SUBROUTINE dust_norm_rhoG(nlst_in,nlst_end,npart,ntot,
     &     xyzmh,ilenx,vxyzu,ilenv,rgrain,ilenr)
c***********************************************************************
c                                                                      *
c     This subroutine calculates the dustnorm() which essentially      *
c     is a measure of how well the dust particle kernel is normalised  *
c     (the values should be unity if the kernel is well sampled).      *
c     It also allows calculation of the gas density in the vicinity of *
c     each dust particle (though this is not currently used anywhere   *
c     else in the code).                                               *
c     The quantity dustnorm() is used in the implementation of         *
c     the method of Loren-Aguilar & Bate (2014) which treats drag      *
c     forces between gas and dust imposing momentum conservation +     *
c     energy dissipation.                                              *
c                                                                      *
c***********************************************************************
c

#ifdef MPIALL
#define IMPLICITNONE
#include "mpi_sup.h"
#else
      IMPLICIT NONE
#endif

c
c--Common variables
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/dustimplicit'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/eosq'
c 
c--MPI headers
c
#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif
c
c--I/O variables
c
      INTEGER   nlst_in, nlst_end, npart, ntot, ilenx, ilenv, ilenr
      REAL      xyzmh(5,ilenx), vxyzu(4,ilenv),rgrain(ilenr)
c
c--Functions
c
      REAL thermal_velocity

c--Local parameters
c
      REAL alpha_turb
      PARAMETER (alpha_turb=1.0E-02)
c
c--Local variables
c
      INTEGER   ipart, j, k, n, ioffset, index, index1, iparttot2
      REAL      m_g,
     &          dr2, drx, dry, drz, h_g, h2_g, h3_g,
     &          wk_sym, wk_sym2, dwdx, dxx, u2_g, vthermal, 
     &          r_planetesimal, rho_planetesimal, mD_Ks, ts, St, vt, Sc,
     &          turb, r, Omega, ddvx, ddvy, ddvz, scalar

c      print *,'Entered dust_norm_rhoG ',nlst_end

c
c--Want to calculate dustnorm using gas particles that are within the
c     *gas* smoothing length of dust.  Cannot simply loop over dust
c     particles to look at its neighbours, because this will not include
c     *active* gas particles for which h_g > h_d.
c     However, can't loop over only active gas particles either, because
c     even with global timesteps, this will miss interactions between
c     dust and ghost gas particles (because ghosts are not active).
c     Therefore, need to loop over both, but make sure interactions
c     are not included twice!
c     With individual timesteps, if gas is inactive, then it will appear
c     in the dust list and so will be done.  If dust is inactive, then
c     it will appear in the gas list, but will need to be done 
c     regardless of whether h_d < h_g or not.
c

C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,ivar,iphase,ijvar,xyzmh,vxyzu)
C$OMP& shared(dumrho,radkernel,ddvtable,dvtable)
C$OMP& shared(wij,dustnorm,cnormk,rhoG)
C$OMP& shared(doublehumpnormk,npart,ntot,iscurrent)
C$OMP& shared(vturb,rgrain)
C$OMP& private(j,k,n,ipart,ioffset,drx,dry,drz,dr2)
C$OMP& private(ddvx,ddvy,ddvz,scalar)
C$OMP& private(m_g,h_g,h2_g,h3_g,u2_g)
C$OMP& private(wk_sym,wk_sym2)
C$OMP& private(dwdx,dxx,index,index1)
C$OMP& private(iparttot2,vthermal)
C$OMP& private(r_planetesimal,rho_planetesimal)
C$OMP& private(mD_Ks,ts,St,vt,Sc,turb,r,Omega)
c
c--Normalization calculation from looping over dust particles
c     (and zeroes the values initially too)
c
C$OMP DO SCHEDULE(runtime)
      DO 20 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
         iparttot2 = ipart
         dustnorm(ipart) = 0.0
         rhoG(ipart) = 0.0
         IF (idim_grow.EQ.idim) vturb(ipart) = 0.0
c
c--Check for non-MPI code that the particle number makes sense.
c     For MPI code, need to define correct index to the xyzmh() array,
c     which is different to most arrays because the tree structure is
c     stored in the array.
c
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#else
         IF (ipart.GT.npart) THEN
            WRITE (*,*) 'ERROR - ipart.GT.npart'
            CALL quit(0)
         ENDIF
#endif
c
c--Loop only through dust particles
c
         IF (iphase(ipart).LT.11) GOTO 20
c
c--Loop over neighbouring particles
c
         ioffset = ivar(2,n)
         DO k = 1, ivar(1,n)
c
c--Get neighbouring particle index from the list
c
            j = ijvar(ioffset + k)
c
c--Loop only through gas particles
c
            IF (iphase(j).EQ.0) THEN

               IF (dumrho(j).EQ.0.0) THEN
                  WRITE (*,*) 'ERROR - dumrho(j).EQ.0.0 ',j
                  CALL quit(0)
               ENDIF
c
c--Dust-gas pair distances
c
               drx = xyzmh(1,iparttot2) - xyzmh(1,j)
               dry = xyzmh(2,iparttot2) - xyzmh(2,j)
               drz = xyzmh(3,iparttot2) - xyzmh(3,j)
               dr2  = drx*drx + dry*dry + drz*drz
c
c--Relative velocities needed for dust growth
c
               IF (idim_grow.EQ.idim) THEN
                  IF (dr2.NE.0.0) THEN
                     ddvx = vxyzu(1,ipart) - vxyzu(1,j)
                     ddvy = vxyzu(2,ipart) - vxyzu(2,j)
                     ddvz = vxyzu(3,ipart) - vxyzu(3,j)
                     scalar = (ddvx*drx + ddvy*dry + ddvz*drz)/SQRT(dr2)
                  ELSE
                     scalar = 0.
                  ENDIF
               ENDIF
c
c--Mass coefficients
c
               m_g     = xyzmh(4,j)
c
c--Kernel (double hump) calculations
c
               h_g      = xyzmh(5,j)
               h2_g     = h_g*h_g
               h3_g     = h2_g*h_g
               u2_g     = dr2/h2_g
               IF (u2_g.LT.radkernel**2) THEN
                  index  = u2_g*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx    = u2_g - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  wk_sym   = doublehumpnormk*(wij(index) + dwdx*dxx)*
     &                 u2_g/h3_g
                  wk_sym2  = cnormk*(wij(index) + dwdx*dxx)/h3_g
c
c--SPH normalization equation
c
                  dustnorm(ipart) = dustnorm(ipart) +
     &                 (m_g/dumrho(j))*wk_sym
c
c--Gas density at dust particle location
c
                  rhoG(ipart) = rhoG(ipart) + m_g*wk_sym2
c
c--Gas turbulent velocity at dust particle location.
c     Only needed in grain growth simulations.
c
                  IF (idim_grow.EQ.idim) THEN
                    r =SQRT(xyzmh(1,iparttot2)**2+xyzmh(2,iparttot2)**2)
                    Omega = SQRT(1./r**3)

                    r_planetesimal = rgrain(ipart)
                    CALL rho_grain(rho_planetesimal)
                    IF (rgrain(ipart).EQ.0.0) THEN
                      WRITE(*,*)'ERROR - rgrain is zero in dust_norm'
                      CALL quit(0)
                    ENDIF

                    vthermal = thermal_velocity(j)
                    mD_Ks    = r_planetesimal*rho_planetesimal/vthermal
                    ts       = mD_Ks/(dumrho(ipart)+dumrho(j))
                    St       = Omega*ts
                    !turb    = sqrt(4.5*alpha_turb/(St+1./St))*vthermal
c--------------------------------SV970----------------------------
                    vt       = SQRT(2.0**0.5*3.0*alpha_turb)*vthermal
                    Sc       = (1.0+St)*SQRT(1.0+(scalar/vt)**2)
                    turb     = 2.0**0.5*vt*(SQRT(Sc-1.0)/Sc)
c--------------------------------SV970----------------------------

                    vturb(ipart) = vturb(ipart) +
     &                     (m_g/dumrho(j))*turb*wk_sym
                  ENDIF

               ENDIF   ! EndIf kernel
            ENDIF   ! EndIf iphase = gas
         END DO   ! End neighbour loop

        IF (dustnorm(ipart).GT.1.5) THEN
           print *,'ERROR -- dustnorm seems too big ',ipart,
     &       dustnorm(ipart),iphase(ipart)
           CALL quit(1)
        ENDIF

 20   END DO
C$OMP END DO
c
c--Now, add contributions from looping over gas particles
c
C$OMP DO SCHEDULE(runtime)
      DO 30 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
         iparttot2 = ipart
c
c--Check for non-MPI code that the particle number makes sense.
c     For MPI code, need to define correct index to the xyzmh() array,
c     which is different to most arrays because the tree structure is
c     stored in the array.
c
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#else
         IF (ipart.GT.npart) THEN
            WRITE (*,*) 'ERROR - ipart.GT.npart'
            CALL quit(0)
         ENDIF
#endif
c
c--Loop only through gas particles
c
         IF (iphase(ipart).NE.0) GOTO 30
         m_g = xyzmh(4,iparttot2)
         h_g = xyzmh(5,iparttot2)
c
c--Loop over neighbouring particles
c
         ioffset = ivar(2,n)
         DO k = 1, ivar(1,n)
c
c--Get neighbouring particle index from the list
c
            j = ijvar(ioffset + k)
c
c--Dust-gas pair distances
c
            drx = xyzmh(1,iparttot2) - xyzmh(1,j)
            dry = xyzmh(2,iparttot2) - xyzmh(2,j)
            drz = xyzmh(3,iparttot2) - xyzmh(3,j)
            dr2  = drx*drx + dry*dry + drz*drz
c
c--Relative velocities needed for dust growth
c
            IF (idim_grow.EQ.idim) THEN
               IF (dr2.NE.0.0) THEN
                  ddvx = vxyzu(1,ipart) - vxyzu(1,j)
                  ddvy = vxyzu(2,ipart) - vxyzu(2,j)
                  ddvz = vxyzu(3,ipart) - vxyzu(3,j)
                  scalar = (ddvx*drx + ddvy*dry + ddvz*drz)/SQRT(dr2)
               ELSE
                  scalar = 0.
               ENDIF
            ENDIF
c
c--If the dust particle is active, only include dust particles for 
c     which h_gas > h_dust (otherwise their contribution has already 
c     been accounted for).  Alternately, if the dust particle is 
c     inactive, need to add the contribution here regardless of how
c     h_gas and h_dust compare.
c
            IF (iphase(j).GE.11 .AND. 
     &           (dr2.GT.(radkernel*xyzmh(5,j))**2 .OR. 
     &           .NOT.iscurrent(j))) THEN
c
c--Kernel (double hump) calculations
c
               h2_g     = h_g*h_g
               h3_g     = h2_g*h_g
               u2_g     = dr2/h2_g
               IF (u2_g.LT.radkernel**2) THEN
                  index  = u2_g*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx    = u2_g - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  wk_sym   = doublehumpnormk*(wij(index) + dwdx*dxx)*
     &                 u2_g/h3_g
                  wk_sym2  = cnormk*(wij(index) + dwdx*dxx)/h3_g
c
c--SPH normalization equation
c

                  IF (dumrho(ipart).EQ.0.0) THEN
                     WRITE (*,*) 'ERROR - dumrho(ipart).EQ.0.0 ',ipart
                     CALL quit(0)
                  ENDIF
C$OMP ATOMIC
                  dustnorm(j) = dustnorm(j) + (m_g/dumrho(ipart))*wk_sym
c
c--Gas density at dust particle location
c
C$OMP ATOMIC
                  rhoG(j) = rhoG(j) + m_g*wk_sym2
c
c--Gas turbulent velocity at dust particle location.
c     Only needed in grain growth simulations.
c
                  IF (idim_grow.EQ.idim) THEN
                    r = SQRT(xyzmh(1,j)**2+xyzmh(2,j)**2)
                    Omega = SQRT(1./r**3)

                    r_planetesimal = rgrain(j)
                    CALL rho_grain(rho_planetesimal)
                    IF (rgrain(j).EQ.0.0) THEN
                      WRITE(*,*)'ERROR - rgrain is zero in dust_norm'
                      CALL quit(0)
                    ENDIF

                    vthermal = thermal_velocity(ipart)
                    mD_Ks    = r_planetesimal*rho_planetesimal/vthermal
                    ts       = mD_Ks/(dumrho(ipart)+dumrho(j))
                    St       = Omega*ts
                    !turb    = sqrt(4.5*alpha_turb/(St+1./St))*vthermal
c--------------------------------SV970----------------------------
                    vt       = SQRT(2.0**0.5*3.0*alpha_turb)*vthermal
                    Sc       = (1.0+St)*SQRT(1.0+(scalar/vt)**2)
                    turb     = 2.0**0.5*vt*(SQRT(Sc-1.0)/Sc)
c--------------------------------SV970----------------------------
C$OMP ATOMIC
                    vturb(j) = vturb(j) +
     &                   (m_g/dumrho(ipart))*turb*wk_sym
                  ENDIF

               ENDIF
            ENDIF
         END DO   ! End neighbour loop
 30   END DO
C$OMP END DO
C$OMP END PARALLEL

      END SUBROUTINE dust_norm_rhoG
