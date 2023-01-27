      SUBROUTINE dust_drag(nlst_in,nlst_end,npart,ntot,
     &     xyzmh,ilenx,vxyzu,ilenv,rgrain,ilenr,dt)
c***********************************************************************
c                                                                      *
c     This subroutine calculates the drag forces between gas and dust  *
c     using the method of Loren-Aguilar & Bate (2014), which           *
c     imposes momentum conservation + energy dissipation               *
c                                                                      *
c***********************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif

c
c--Functions
c
      REAL rho_grain
c
c--Common variables
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/dustimplicit'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/fullstep'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/f2'
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
c--Uncomment this line if you want to check momentum conservation 
c
c#define debug_drag
c
c--I/O variables
c
      INTEGER   nlst_in, nlst_end, npart, ntot, ilenx, ilenv
      REAL      xyzmh(5,ilenx), vxyzu(4,ilenv), rgrain(ilenr), dt
c
c--Functions
c
      REAL thermal_velocity
c
c--Local variables
c
      INTEGER   ipart, j, k, n, ioffset, index, index1, idust, igas
      REAL      m_ipart, m_j, scalar, ddvx, ddvy, ddvz, dr, 
     &          dr2, drx, dry, drz, chi, h_g, h2_g, h3_g, u2_g,
     &          wk_sym, nu, lambda,
     &          mD_Ks, dwdx, dxx, vthermal, ts
      PARAMETER (nu=3.0)
#ifdef debug_drag
      REAL      px1, py1, pz1, px0, py0, pz0
      REAL      Lx1, Ly1, Lz1, Lx0, Ly0, Lz0, L0, L1
#endif

c      print *,'Entered dust_drag ',nlst,nlst0
c
c--Want to calculate the change in velocity due to drag between dust and
c     gas particles that are within the *gas* smoothing length of dust.
c     Cannot simply loop over dust particles to look at its neighbours, 
c     because this will not include *active* gas particles for which 
c     the separation is less than radkernel*h_g, but greater than
c     radkernel*h_d.
c     However, can't loop over only active gas particles either, because
c     even with global timesteps, this will miss interactions between
c     dust and ghost gas particles (because ghosts are not active).
c     Therefore, need to loop over both, but make sure interactions
c     are not included twice!
c     With individual timesteps, if gas is inactive, then it will appear
c     in the dust list if the separation less than radkernel*h_g 
c     regardless of the value of h_d, and so will be done.  If dust is 
c     inactive, then it will appear in the gas list, and the interaction
c     will need to be done regardless of whether h_d < h_g or not.
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,ivar,iphase,ijvar,xyzmh,vxyzu)
C$OMP& shared(dumrho,dt,radkernel,ddvtable,dvtable,nlst0)
C$OMP& shared(wij,dustnorm,r_planetesimals,rho_planetesimals)
C$OMP& shared(utime,udist,cnormk)
C$OMP& shared(doublehumpnormk,npart,ntot,imaxstep,isteps)
C$OMP& shared(n_planetesimal_types,vxyzu_change)
C$OMP& shared(vsound,ifullstep,iscurrent)
C$OMP& shared(f1vxyzu,f2vxyzu)
C$OMP& shared(rgrain)
C$OMP& private(j,k,n,ipart,ioffset,drx,dry,drz,dr,ddvx,ddvy,ddvz)
C$OMP& private(m_ipart,m_j,dr2,igas,idust)
C$OMP& private(h_g,h2_g,h3_g,chi,u2_g,lambda)
C$OMP& private(wk_sym,ts)
C$OMP& private(dwdx,dxx,index,index1,beta,scalar,deltat)
C$OMP& private(r_planetesimal,rho_planetesimal)
C$OMP& private(vthermal,mD_Ks,radius,r2,iparttot2)
C$OMP& private(ax_i,ay_i,az_i,ax_j,ay_j,az_j,scalar_accel,term)
c
c--Drag calculation
c
C$OMP DO SCHEDULE(runtime)
      DO 40 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
         iparttot2 = ipart
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#endif
c
c--Loop through dust *and* gas particles
c
         IF (iphase(ipart).GT.0 .AND. iphase(ipart).LT.11) GOTO 40
c
c--Check that dumrho() is defined (to avoid divide by zero)
c
         IF (dumrho(ipart).EQ.0.0) THEN
            WRITE (*,*) 'ERROR - idumrho(ipart).EQ.0.0 ',ipart,
     &           iphase(ipart)
            CALL quit(0)
         ENDIF
c
c--Loop over neighbouring particles
c
         ioffset = ivar(2,n)
         DO 30 k = 1, ivar(1,n)
c
c--Get neighbouring particle index from the list
c
            j = ijvar(ioffset + k)  
c
c--Loop only through particles of opposite type to ipart
c
            IF (iphase(ipart).EQ.0 .AND. iphase(j).GE.11 .OR.
     &           iphase(ipart).GE.11 .AND. iphase(j).EQ.0) THEN
c  
c--Dust-gas pair distances (order doesn't matter, because the sign of
c     the velocity change depends both on the unit vectors, and whether
c     it is gas-dust or dust-gas below.
c
               drx = xyzmh(1,iparttot2)-xyzmh(1,j)
               dry = xyzmh(2,iparttot2)-xyzmh(2,j)
               drz = xyzmh(3,iparttot2)-xyzmh(3,j)
               dr2 = drx*drx + dry*dry + drz*drz

               IF (dumrho(j).EQ.0.0) THEN
                  WRITE (*,*) 'ERROR - dumrho(j).EQ.0.0 ',j,iphase(j)
                  CALL quit(0)
               ENDIF
c
c--Define gas or dust quantities, depending on the phase of ipart
c
               IF (iphase(ipart).EQ.0) THEN
                  igas = ipart
                  idust = j
                  h_g = xyzmh(5,iparttot2)
c
c--If we are looping over dust neighbours of a gas particle, then if the
c     dust particle is active, we only need to compute drag contributions
c     if the separation is greater than radkernel*h_dust, otherwise the 
c     interaction will be done by the loop over the active dust particle.
c     If the dust particle is not active (shouldn't happen if dust 
c     always has a smaller timestep than the gas), then the interaction 
c     needs to be done when looping over the gas regardless of the
c     values of h_gas and h_dust.
c
                  IF (dr2.LE.(radkernel*xyzmh(5,j))**2 .AND. 
     &                 iscurrent(idust)) GOTO 30
               ELSE
                  igas = j
                  idust = ipart
                  h_g = xyzmh(5,j)
               ENDIF
c
c--Dust particle time-step
c
               deltat = dt*isteps(idust)/imaxstep
c
c--Drag coefficient -- in code units (Epstein regime)
c     NOTE: r_planetesimal & rho_planetesimal specified in code units
c     Example: 1km is r_planetesimal = 1.0E5/udist
c
               IF (idim_grow.EQ.idim) THEN
                  r_planetesimal   = rgrain(idust)
                  IF (rgrain(idust).EQ.0.0) THEN
                     WRITE(*,*)'ERROR - rgrain is zero in dust drag'
                     CALL quit(0)
                  ENDIF
               ELSE
                  IF (n_planetesimal_types.EQ.0) THEN
                     WRITE(*,*) 'ERROR - you are using fixed-size ',        
     &                    'implicit dust-drag but have not defined ',
     &                    'any dust grain size in your ASCII file'
                     CALL quit(0)
                  ENDIF
  
                  IF (iphase(idust)-10.LE.n_planetesimal_types) THEN
                     r_planetesimal = r_planetesimals(iphase(idust)-10)
                     rho_planetesimal = rho_grain(iphase(idust))
                  ELSE
                     WRITE(*,*)'ERROR - Not enough n_planetesimal_types'
                     CALL quit(0)
                  ENDIF
               ENDIF
c
c--Calculate unit vectors
c
               dr = SQRT(dr2)
               IF (dr.NE.0.0) THEN
                  drx = drx/dr
                  dry = dry/dr
                  drz = drz/dr
               ELSE
                  drx = 0.0
                  dry = 0.0
                  drz = 0.0
               ENDIF
c
c--Velocity differences (order doesn't matter, since scalar involves
c     dot product).
c
               ddvx = vxyzu(1,ipart)-vxyzu(1,j)
               ddvy = vxyzu(2,ipart)-vxyzu(2,j)
               ddvz = vxyzu(3,ipart)-vxyzu(3,j)
               scalar = ddvx*drx + ddvy*dry + ddvz*drz
c
c--Acceleration differences (excluding gas-dust drag)
c
c
c--Gas and dust particle accelerations (excluding drag).
c     Assumes dust timestep is always smaller than or equal to gas.
c
               IF (ifullstep(ipart)) THEN
                  ax_i = f2vxyzu(1,ipart)
                  ay_i = f2vxyzu(2,ipart)
                  az_i = f2vxyzu(3,ipart)
               ELSE
                  ax_i = f1vxyzu(1,ipart)
                  ay_i = f1vxyzu(2,ipart)
                  az_i = f1vxyzu(3,ipart)
               ENDIF
               IF (ifullstep(j)) THEN
                  ax_j = f2vxyzu(1,j)
                  ay_j = f2vxyzu(2,j)
                  az_j = f2vxyzu(3,j)
               ELSE
                  ax_j = f1vxyzu(1,j)
                  ay_j = f1vxyzu(2,j)
                  az_j = f1vxyzu(3,j)
               ENDIF
               IF (igas.EQ.ipart) THEN
                  IF (.NOT.iscurrent(ipart)) THEN
                     ax_i = f1vxyzu(1,ipart)
                     ay_i = f1vxyzu(2,ipart)
                     az_i = f1vxyzu(3,ipart)
                  ENDIF
               ELSE
                  IF (.NOT.iscurrent(j)) THEN
                     ax_j = f1vxyzu(1,j)
                     ay_j = f1vxyzu(2,j)
                     az_j = f1vxyzu(3,j)
                  ENDIF
               ENDIF
               scalar_accel = (ax_i-ax_j)*drx + (ay_i-ay_j)*dry +
     &              (az_i-az_j)*drz               
c
c--Mass coefficients
c
               m_ipart = xyzmh(4,iparttot2)
               m_j     = xyzmh(4,j)
c
c--Kernel (double hump) calculations
c
               dr2      = dr*dr
               h2_g     = h_g*h_g
               h3_g     = h2_g*h_g
               u2_g     = dr2/h2_g

               IF (u2_g.LT.radkernel**2) THEN
                  index  = u2_g*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx  = u2_g - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  wk_sym = doublehumpnormk*(wij(index) + dwdx*dxx)*u2_g/
     &                 h3_g
c
c--Calculate drag coefficient
c
c               vthermal = 1.0E5*utime/udist ! Fixed value in code units
                  vthermal = thermal_velocity(igas)
                  mD_Ks = r_planetesimal*rho_planetesimal/vthermal
c
c--Constant drag
c
c              ts   = mD_Ks/(1.0+m_ipart/mG(ipart))
c
c--Epstein profile
c
                  ts = mD_Ks/(dumrho(ipart)+dumrho(j))
c
c--Non-linear profile
c
c              ts   = mD_Ks/((dumrho(ipart)+dumrho(j))*vDG(ipart))
c
c--Chi calculation
c
                  IF (ts.GT.0.0) THEN
                     beta = 1.0/(dumrho(ipart)+ dumrho(j))
                     IF (ifullstep(idust)) THEN
                        chi = beta*(1.0-dexp(-deltat/ts))
                     ELSE
                        chi = beta*(1.0-dexp(-0.5*deltat/ts))
                     ENDIF
                     lambda = (deltat+ts)*chi - deltat*beta
                  ELSE
                     WRITE (*,*) 'ERROR - ts.LE.0.0 ',ts,idust,igas
                     CALL quit(1)
                  ENDIF
c
c--Careful with norm!!!
c
                  IF (dustnorm(idust).EQ.0.0) THEN
                     WRITE (*,*) 'ERROR - dustnorm(idust).EQ.0.0 ',
     &                    idust,dustnorm(idust),ipart,ivar(1,n),j,
     &                    u2_g,wk_sym
                     CALL quit(1)
                  ENDIF
c
c--Change dust and gas velocities
c  If lambda=0 & scalar_accel=0 reverts to Loren-Aguilar & Bate (2014)
c                  lambda = 0.
c                  scalar_accel = 0.

                  term = nu*m_j*(chi*scalar-lambda*scalar_accel)*
     &                 wk_sym/dustnorm(idust)
C$OMP ATOMIC
                  vxyzu_change(1,ipart) = vxyzu_change(1,ipart) - 
     &                 term*drx
C$OMP ATOMIC
                  vxyzu_change(2,ipart) = vxyzu_change(2,ipart) - 
     &                 term*dry
C$OMP ATOMIC
                  vxyzu_change(3,ipart) = vxyzu_change(3,ipart) - 
     &                 term*drz
c
c--Other particle:
c
                  term = term*m_ipart/m_j
C$OMP ATOMIC
                  vxyzu_change(1,j) = vxyzu_change(1,j) + term*drx
C$OMP ATOMIC
                  vxyzu_change(2,j) = vxyzu_change(2,j) + term*dry
C$OMP ATOMIC
                  vxyzu_change(3,j) = vxyzu_change(3,j) + term*drz
c
c--Change gas internal energy
cC$OMP ATOMIC
c               vxyzu_change(4,igas)=vxyzu_change(4,igas) + nu*m_ipart*
c     &              chi*scalar*scalar*wk_sym/dustnorm(idust)
               ENDIF
            ENDIF
 30      END DO                 ! End of neighbours loop
 40   END DO                    ! End of active particles loop
C$OMP END DO 
C$OMP END PARALLEL
c
c--Check momentum after dust-gas drag interaction (NOTE: Only works for
c     global timesteps, with no ghost boundaries).
c
#ifdef debug_drag
      px1 = 0.0
      py1 = 0.0
      pz1 = 0.0
      Lx1 = 0.0
      Ly1 = 0.0
      Lz1 = 0.0
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,npart,xyzmh,vxyzu_change,ivar)
C$OMP& private(n,ipart,iparttot2)
C$OMP& reduction(+:px1)
C$OMP& reduction(+:py1)
C$OMP& reduction(+:pz1)
C$OMP& reduction(+:Lx1)
C$OMP& reduction(+:Ly1)
C$OMP& reduction(+:Lz1)
      DO 60 n = nlst_in, nlst_end
         ipart = ivar(3,n)
         iparttot2 = ipart
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#endif
         px1 = px1 + xyzmh(4,iparttot2)*vxyzu_change(1,ipart)
         py1 = py1 + xyzmh(4,iparttot2)*vxyzu_change(2,ipart)
         pz1 = pz1 + xyzmh(4,iparttot2)*vxyzu_change(3,ipart)
         Lx1 = Lx1 + xyzmh(4,iparttot2)*
     &        (xyzmh(2,iparttot2)*vxyzu_change(3,ipart) -   
     &        vxyzu_change(2,ipart)*xyzmh(3,iparttot2))
         Ly1 = Ly1 - xyzmh(4,iparttot2)*
     &        (xyzmh(1,iparttot2)*vxyzu_change(3,ipart) -   
     &        vxyzu_change(1,ipart)*xyzmh(3,iparttot2))
         Lz1 = Lz1 + xyzmh(4,iparttot2)*
     &        (xyzmh(1,iparttot2)*vxyzu_change(2,ipart) -   
     &        vxyzu_change(1,ipart)*xyzmh(2,iparttot2))
 65      CONTINUE
 60   END DO
C$OMP END PARALLEL DO
      L1 = SQRT(Lx1*Lx1 + Ly1*Ly1 + Lz1*Lz1)
      PRINT *, 'Momentum loss in dust-gas drag =',px1,py1,pz1
      PRINT *, 'Angular momentum loss in dust-gas drag =',L1
#endif

      END SUBROUTINE dust_drag
