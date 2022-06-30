      SUBROUTINE hy09dustevol(nlst,llist,npart,ntot,xyzmh,
     &     vxyzu,ekcle,HY09_bin_rho,rho,HY09_drhodt,HY09vrel_h)
c
c--Set evolve dust densities for Hirashita & Yan (2009) and 
c     Hirashita & Omukai (2009) method of dust coagulation.
c
c     Matthew R. Bate, 21 July 2021
c
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/HY09accel'

      INTEGER nlst,npart,ntot
      INTEGER llist(idim)
      REAL xyzmh(5,idim), vxyzu(4,idim), ekcle(5,iradtrans)
      REAL*8 HY09_bin_rho(HY09_ndust_bins,idimHY09)
      REAL*4 rho(idim)
      REAL*8 HY09_drhodt(HY09_ndust_bins,idimHY09)
      REAL*8 HY09vrel_h(idimHY09)
c
c--Functions
c
      REAL HY09_alpha, HY09_vreldust
c
c--Local variables
c
      INTEGER n, ipart, i, k, j, istart, iend, ilocation
      REAL sum, sum_mass, rhoi, radius, temperature, vthermal
      REAL vrel_press_coeff, vrel_rad_drift_coeff, omega, vg2
      REAL vrel_settling_coeff, freqJeanssoundcrossing
      REAL Stokes_k, Stokes_j, size_k
      REAL rgrav_accel_code_units, rgrav_accel_cgs, dlnPdlnr, eta
      REAL vel_dot_grav_accel, vel2, vel_tangential
      REAL ratio_vtan_to_circular_orbit, ratio_vrad_to_vtan
      REAL ratio_vrad_to_cs
      LOGICAL istartset,iDisc
c
c--Loop over active gas particles
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst,llist,iphase,rho,encal,xyzmh,vxyzu,ekcle)
C$OMP& shared(HY09_bin_rho,HY09_drhodt,HY09binmass,HY09binmass_max)
C$OMP& shared(udens,utime,udist,HY09vrel_h,HY09_dust_density)
C$OMP& shared(gas_accel,grav_accel,HY09discflag,vsound)
C$OMP& private(n,ipart,rhoi,radius,temperature,vthermal)
C$OMP& private(vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2)
C$OMP& private(vrel_settling_coeff,freqJeanssoundcrossing)
C$OMP& private(rgrav_accel_code_units,rgrav_accel_cgs,dlnPdlnr,eta)
C$OMP& private(i,j,k,sum,sum_mass,istart,iend,istartset)
C$OMP& private(ilocation,Stokes_k,Stokes_j,size_k)
C$OMP& private(vel_dot_grav_accel,vel2,vel_tangential,ratio_vrad_to_cs)
C$OMP& private(ratio_vtan_to_circular_orbit,ratio_vrad_to_vtan)
C$OMP& private(iDisc)

      DO n = 1, nlst
         ipart = llist(n)

         IF (ipart.EQ.0) THEN
            PRINT *,'ERROR - ipart= ',ipart,'n= ',n,nlst
            CALL quit(0)
         ENDIF

         IF (iphase(ipart).NE.0) CYCLE
c
c--Pre-compute values that don't depend on grain size
c
         rhoi = rho(ipart)
         radius = SQRT(xyzmh(1,ipart)**2 + xyzmh(2,ipart)**2 + 
     &        xyzmh(3,ipart)**2)
         IF (encal.EQ.'r') THEN
            temperature = vxyzu(4,ipart)/ekcle(3,ipart)
         ELSE
            temperature = 10.0
         ENDIF
         vthermal = SQRT(8.0*boltzmannk*temperature/(pi*gmw*mH))
c
c--Pre-compute coefficient for vrel from pressure gradient
c
         vrel_press_coeff = HY09_dust_density/(rhoi*udens*vthermal)*
     &        SQRT(gas_accel(1,ipart)**2 + gas_accel(2,ipart)**2 + 
     &        gas_accel(3,ipart)**2)*udist/(utime**2)
c
c--Pre-compute coefficient for radial drift.  The rgrav_accel_cgs
c     is used because this is equal to Omega*v_Keplerian, with the 
c     Omega coming from the definition of the Stokes number.
c
c     Compute the magnitude of the gravitational acceleration, GM/r^2
c
         rgrav_accel_code_units = SQRT(grav_accel(1,ipart)**2 +
     &           grav_accel(2,ipart)**2 +
     &           grav_accel(3,ipart)**2)
         IF (rgrav_accel_code_units.EQ.0.) THEN
            WRITE (*,*) 'ERROR - rgrav_accel.EQ.0.'
            CALL quit(0)
         ENDIF
         rgrav_accel_cgs = rgrav_accel_code_units*udist/(utime)**2
c
c--The computation of "eta", determining how the gas pressure gradient
c     leads to the gas circular speed being different to Keplerian.
c     The value below assumes H/r=0.05, and the dlnPdlnr is a guess.
c
c         dlnPdlnr = -1.0
c         eta = - 0.5*(0.05)**2*dlnPdlnr
c
c--This is the way to determine "eta" based on the pressure acceleration
c     and the gravitational acceleration directly, without assumptions
c     about H/r and having to measure dlnPdlnr.  "eta" is dimensionless
c     so using quantities in code units is fine.
c
         eta = - 0.5*(gas_accel(1,ipart)*grav_accel(1,ipart) +
     &        gas_accel(2,ipart)*grav_accel(2,ipart) +
     &        gas_accel(3,ipart)*grav_accel(3,ipart))/
     &        rgrav_accel_code_units**2

         vrel_rad_drift_coeff = -2.0*HY09_dust_density/
     &      (rhoi*udens*vthermal)*eta*rgrav_accel_cgs
c
c--Pre-compute quantities for turbulence relative velocities
c
         omega = SQRT(rgrav_accel_cgs/(radius*udist+solarr))
         vg2 = pi/8.0*vthermal**2
         freqJeanssoundcrossing = 2.0*SQRT(gg*rhoi*udens/pi)
c
c--Pre-compute coefficient for vertical settling
c
         vrel_settling_coeff = - xyzmh(3,ipart)*udist*omega
c
c--Determine whether or not particle is in a disc (determines whether
c     or not the relative velocities for dust grain growth include:
c       pressure gradient driven terminal velocity (envelope only)
c       turbulent relative velocities (disc only)
c       radial drift (disc only)
c     Based on ratio of radial velocity to tangential velocity (with
c       the direction being defined by the local gravitational 
c       acceleration), and the ratio of the tangential velocity to
c       the circular orbital speed (defined by local gravitational 
c       acceleration and radius).
c
         vel_dot_grav_accel = (vxyzu(1,ipart)*grav_accel(1,ipart) +
     &        vxyzu(2,ipart)*grav_accel(2,ipart) +
     &        vxyzu(3,ipart)*grav_accel(3,ipart)) /
     &        rgrav_accel_code_units
         vel2 = vxyzu(1,ipart)**2 +vxyzu(2,ipart)**2 +vxyzu(3,ipart)**2
         vel_tangential = SQRT( vel2 - vel_dot_grav_accel**2 )
         ratio_vrad_to_vtan = vel_dot_grav_accel/vel_tangential
         ratio_vrad_to_cs = vel_dot_grav_accel/vsound(ipart)
c         ratio_vtan_to_circular_orbit = vel_tangential/
c     &        SQRT(rgrav_accel_code_units*(radius+solarr/udist))
c         HY09discflag(ipart) = ratio_vtan_to_circular_orbit
c
c--Disc flag is that ABS(ratio_vrad_to_vtan).LT.0.5
c     I also tried: ratio_vtan_to_circular_orbit.GT.0.6 and
c     ABS(ratio_vrad_to_vtan).LT.0.5
c     but the additional criterion only helps at the outer boundary
c     of the cloud where dust coagulation is negligible anyway.
c     This seems to provide good disc extraction, for example, 
c     excluding particles in the high "atmosphere" of a disc.
c
c         IF (ratio_vtan_to_circular_orbit.GT.0.6) THEN
         IF (ABS(ratio_vrad_to_vtan).LT.0.5) THEN
            HY09discflag(ipart) = ratio_vrad_to_vtan
            iDisc = .TRUE.
c
c     As pointed out in Bate (2022), when including envelope turbulence
c     in non-rotating cloud calculations it is unreasonable to have the 
c     first hydrostatic core have envelope-type turbulence rather than
c     more quiescent disc-type turbulence, so the additional criteria
c     below was added to use disc-type turbulence if the radial
c     velocity is smaller than 10% of the sound speed.  This selects
c     the first hydrostatic core in non-rotating calculations as being
c     disc-type turbulence.
c
         ELSEIF (ABS(ratio_vrad_to_cs).LT.0.1) THEN
            HY09discflag(ipart) = ratio_vrad_to_cs
            iDisc = .TRUE.
         ELSE
            HY09discflag(ipart) = 1.0
            iDisc = .FALSE.
         ENDIF
c
c--Find empty bins to skip would optimise speed (not fully tested,
c     hence commented out)
c
         istart = 1
         iend = HY09_ndust_bins
c         istartset = .FALSE.
c         DO i = 1, HY09_ndust_bins
c            IF (istartset .AND. HY09_bin_rho(i,ipart).LT.1.0E-30)
c     &           THEN
c               iend = i
c               EXIT
c            ENDIF
c            IF (.NOT.istartset .AND. HY09_bin_rho(i,ipart).LT.1.0E-30)
c     &           THEN
c               istart = i
c            ELSE
c               istartset = .TRUE.
c            ENDIF
c         END DO
c
c--Find unused bins (optimises speed)
c
         DO i = HY09_ndust_bins, 1, -1
            IF (HY09_bin_rho(i,ipart).GT.0.) THEN
               iend = i+1
               EXIT
            ENDIF
         END DO
c
c--Loop over dust bins
c
         DO i = 1, HY09_ndust_bins
c
c--Compute first terms -- NOTE: rhoi only appears once because the
c     result (drhodt) is actually (drho/dt)/rho_G.  To get everything
c     in cgs units, the "rhoi" which is in code units would need to
c     be multiplied by udens, but this is done once at the end when we
c     also multiply by utime to get (drhodt) to use code units of time
c     for when bin_rho is integrated in step_P.F
c
            sum = 0.
            DO k = istart, iend
               sum = sum + 
     &            HY09_alpha(ipart,k,i,temperature,radius,rhoi,vthermal,
     &              vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &              vrel_settling_coeff,freqJeanssoundcrossing,iDisc)*
     &              HY09_bin_rho(k,ipart)
            END DO
            HY09_drhodt(i,ipart) = - HY09binmass(i)*
     &           HY09_bin_rho(i,ipart)*sum*rhoi
c
c--Compute second terms -- NOTE: rhoi only appears once because the
c     result (drhodt) is actually (drho/dt)/rho_G
c
            sum = 0.
            DO j = istart, iend
               DO k = istart, iend
                  sum_mass = HY09binmass(k) + HY09binmass(j)
                  IF (sum_mass.GE.HY09binmass_max(i-1) .AND.
     &                 sum_mass.LT.HY09binmass_max(i)) THEN
                     sum = sum + 
     &         (HY09_alpha(ipart,k,j,temperature,radius,rhoi,vthermal,
     &          vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &          vrel_settling_coeff,freqJeanssoundcrossing,iDisc)*
     &                    HY09_bin_rho(k,ipart))*
     &                    HY09_bin_rho(j,ipart)*sum_mass
c
c--Original algorithm from HY09 below, doesn't conserve dust mass
c
c     &                    HY09_bin_rho(j,ipart)*HY09binmass(i)
                  ENDIF
               END DO
            END DO
            HY09_drhodt(i,ipart) = HY09_drhodt(i,ipart) + 
     &           0.5*sum*rhoi

         END DO ! End loop over HY09_ndust_bins
c
c--Convert result to code units. Most quantities above are in cgs units,
c     but HY09_bin_rho is \tilda{rho}/gas_density and when we multiply
c     by the gas density we use values in code units (i.e. we miss a 
c     factor of udens).  The RHS of the full cgs equation has two 
c     \tidla{rho} terms (i.e. dust density squared), but since our
c     definition of bin_rho is actually \tidla{rho}/gas_density and we 
c     want to evaluate (drho/dt)/gas_density on the LHS, then only the
c     single factor of udens is missing (not two).
c     Finally, we need to multiply by utime to get (drho/dt)/rho_G into
c     code units for when it is used to integrate bin_rho in step_P.F.
c
         HY09_drhodt(:,ipart) = HY09_drhodt(:,ipart)*(udens*utime)
c
c--Determine relative velocity for peak grain size (w.r.t. gas, assuming
c     smallest grain size is tied to gas), and divide by the smoothing
c     length.  When muliplied by the timestep, this then gives the 
c     distance that the grain has moved as a fraction of h.
c     NOTE: This quantity is computed in code units, not cgs.
c
         ilocation = MAXLOC(HY09_bin_rho(:,ipart), DIM=1 )
c
c--Turn off Brownian motion for computing relative speed between gas
c     and the most common dust grains (because Brownian motion is
c     a maximum for the smallest dust grains, whereas for all other
c     sources of dust velocity the smallest grains are tied to the gas.
c
         HY09vrel_h(ipart) = HY09_vreldust(ipart,ilocation,1,
     &        temperature,radius,rhoi,vthermal,
     &        vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &        vrel_settling_coeff,freqJeanssoundcrossing,
     &        Stokes_k,Stokes_j,size_k,.FALSE.,iDisc)*utime/udist / 
     &        xyzmh(5,ipart)

      END DO ! End loop over gas particles      
C$OMP END PARALLEL DO

      RETURN
      END
