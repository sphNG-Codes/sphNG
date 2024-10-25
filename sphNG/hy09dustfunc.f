c
c--Contains various function definitions for Hirashita & Yan (2009)
c     method of dust coagulation.
c
c     Matthew R. Bate: 19 July 2021
c
c-----------------------------------------------------------

      FUNCTION HY09_alpha(ipart,k,j,temperature,radius,rhoi,vthermal,
     &     vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &     vrel_settling_coeff,freqJeanssoundcrossing,iDisc)
c
c--Growth function, alpha
c
      IMPLICIT NONE

      INCLUDE 'COMMONS/HY09dustprops'
c
c--Function types
c
      REAL HY09_vreldust,HY09_vcoag,HY09_cross_sec

      INTEGER ipart,k,j
      REAL HY09_alpha,temperature,radius,rhoi,vthermal
      REAL vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2
      REAL vrel_settling_coeff,freqJeanssoundcrossing
      LOGICAL iDisc

      REAL vrel

      REAL Stokes_k,Stokes_j,size_k

      vrel = HY09_vreldust(ipart,k,j,temperature,radius,rhoi,vthermal,
     &     vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &     vrel_settling_coeff,freqJeanssoundcrossing,
     &     Stokes_k,Stokes_j,size_k,.TRUE.,iDisc)
      IF (vrel.LT.HY09_vcoag(k,j)) THEN
         HY09_alpha = HY09_cross_sec(k,j)*vrel/
     &        (HY09binmass(k)*HY09binmass(j))
      ELSE
         HY09_alpha = 0.
      ENDIF

      RETURN
      END FUNCTION HY09_alpha

c-----------------------------------------------------------

      FUNCTION HY09_vreldust(ipart,k,j,temperature,radius,rhoi,vthermal,
     &     vrel_press_coeff,vrel_rad_drift_coeff,omega,vg2,
     &     vrel_settling_coeff,freqJeanssoundcrossing,
     &     Stokes_k,Stokes_j,size_k,iBrownian,iDisc)
c
c--Relative dust velocity function.  Computed in cgs units.
c
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09accel'
      INCLUDE 'COMMONS/dustfluidvelu'
      INCLUDE 'COMMONS/typef'

      LOGICAL iBrownian, iDisc
      REAL HY09_vreldust,temperature,radius,rhoi,vthermal
      INTEGER ipart,k,j
      REAL vrel_brownian2, vrel_radial_disc, vrel_turb2, vrel_terminal
      REAL vrel_settling
      REAL eta, dlnPdlnr,Stokes_max,Stokes_min,radius_use, rgrav_accel
      REAL stoppingtime_k, stoppingtime_j, Stokes_k, Stokes_j
      REAL sqrt_inv_Re, omega, size_k
      REAL vrel_press_coeff, vrel_rad_drift_coeff, vg2
      REAL vrel_settling_coeff, freqJeanssoundcrossing
      REAL HY09discflag_i
      REAL epsilon, vg2_coeff_disc_or_env
c
c--Brownian motion (vrel^2 to avoid lots of square roots)
c
      vrel_brownian2 = 8.0*boltzmannk*temperature/
     &     (pi*HY09binmass(k)*HY09binmass(j)/
     &     (HY09binmass(k) + HY09binmass(j)))
c
c--Velocity difference due to terminal velocity from pressure gradient
c
      IF (.TRUE. .AND. .NOT.iDisc) THEN
c      IF (.TRUE. .AND. iDisc) THEN
c      IF (.TRUE.) THEN
         vrel_terminal = (HY09binsizes(k)-HY09binsizes(j))*
     &        vrel_press_coeff
c     &        HY09_dust_density/(rhoi*udens*vthermal)*
c     &        SQRT(gas_accel(1,ipart)**2 + gas_accel(2,ipart)**2 + 
c     &        gas_accel(3,ipart)**2)*udist/(utime**2)
      ELSE
         vrel_terminal = 0.
      ENDIF
c
c--Radial drift speed (assumes Stokes number <1)
c
      IF (.TRUE. .AND. iDisc .AND. iexf.NE.11) THEN
c      IF (.TRUE. .AND. .NOT.iDisc) THEN
c      IF (.FALSE.) THEN
         dlnPdlnr = -1.0
         eta = - 0.5*(0.05)**2*dlnPdlnr
         IF (.FALSE.) THEN
         radius_use = MAX(radius,0.001)
         vrel_radial_disc = -2.0*(HY09binsizes(k)-HY09binsizes(j))*
     &        HY09_dust_density/(rhoi*udens*vthermal)*
     &        eta* gg*(solarm*0.1)/
c     &        eta* gg*(solarm*(radius_use*udist/(0.1*pc))**3)/
     &        (radius_use*udist)**2
         ELSE
c            rgrav_accel = SQRT(grav_accel(1,ipart)**2 + 
c     &           grav_accel(2,ipart)**2 +
c     &           grav_accel(3,ipart)**2)*udist/(utime)**2
c            vrel_radial_disc = -2.0*(HY09binsizes(k)-HY09binsizes(j))*
c     &           HY09_dust_density/(rhoi*udens*vthermal)*eta*
c     &           rgrav_accel
            vrel_radial_disc = vrel_rad_drift_coeff*
     &           (HY09binsizes(k)-HY09binsizes(j))
         ENDIF
      ELSE
         vrel_radial_disc = 0.
      ENDIF
c
c--Set stopping times and Stokes numbers
c
      stoppingtime_j = HY09binsizes(j)*
     &     HY09_dust_density/(rhoi*udens*vthermal) ! Epstein regime
      stoppingtime_k = HY09binsizes(k)*
     &     HY09_dust_density/(rhoi*udens*vthermal) ! Epstein regime

      IF (.FALSE.) THEN
         Stokes_max = MAX(stoppingtime_j,stoppingtime_k)*
     &        SQRT(gg*(solarm*0.1)/
c     &        SQRT(gg*(solarm*(radius_use*udist/(0.1*pc))**3)/
     &        (radius_use*udist)**3)
      ELSE
         Stokes_j =stoppingtime_j*omega
         Stokes_k =stoppingtime_k*omega
         Stokes_max = MAX(Stokes_j,Stokes_k)
         size_k = HY09binsizes(k)
      ENDIF
c
c--Settling relative speed
c
      IF (iexf.EQ.11) THEN
         vrel_settling = udist/utime*
     &        (dustfluidvxyz(3,k,ipart) - dustfluidvxyz(3,j,ipart))
      ELSEIF (.TRUE. .AND. iDisc) THEN
         vrel_settling = vrel_settling_coeff*(Stokes_k-Stokes_j)
      ELSE
         vrel_settling = 0.
      ENDIF
c
c--Turbulent relative speed
c
      IF (.FALSE.) THEN
c--From Adriens literature review -- equal mass particles
         vrel_turb2 = 3.0*1.0E-03/(Stokes_max+1.0/Stokes_max)*
     &        pi/8.0*vthermal**2
      ELSEIF (.TRUE.) THEN
c      ELSEIF (.TRUE. .AND. iDisc) THEN
c      ELSEIF (.TRUE. .AND. .NOT.iDisc) THEN
c      ELSEIF (.FALSE.) THEN
c
c--Inverse square root of Reynolds numnber (assumed constant)
c
         sqrt_inv_Re = 1.0E-4

         IF (iDisc) THEN
            vg2_coeff_disc_or_env = 0.001
c            vg2_coeff_disc_or_env = 1.5
         ELSE
c
c--For envelope turbulence, redefine Vg2 and Stokes numbers from
c     those used for disc turbulence.  Means that from this point on
c     in this subroutine, the Stokes numbers will be either disc
c     or envelope values (prior to this point all numbers are for
c     disc turbulence).
c
            vg2_coeff_disc_or_env = 1.5
            Stokes_max = Stokes_max*freqJeanssoundcrossing/omega
            Stokes_j = Stokes_j*freqJeanssoundcrossing/omega
            Stokes_k = Stokes_k*freqJeanssoundcrossing/omega
         ENDIF
c         vg2 = 1.0E-03*pi/8.0*vthermal**2

c         vrel_turb2 = 1.0E-03*pi/8.0*vthermal**2*
c     &        ABS( (Stokes-sqrt_inv_Re) + 
c     & Stokes_k**2*(1.0/(Stokes_k+Stokes) -1.0/(Stokes_k+sqrt_inv_Re))+
c     & Stokes_j**2*(1.0/(Stokes_j+Stokes) -1.0/(Stokes_j+sqrt_inv_Re)) )

         epsilon = Stokes_j/Stokes_k
         IF (epsilon.GT.1.) epsilon = 1./epsilon

c--Ormel & Cuzzi, equation 26
c         vrel_turb2 = 1.0E-03*pi/8.0*vthermal**2*
c     &        ABS( (Stokes_k-Stokes_j)*
c     &   (Stokes_k**2/(Stokes_k+sqrt_inv_Re)-
c     &   Stokes_j**2/(Stokes_j+sqrt_inv_Re))/
c     &        (Stokes_k+Stokes_j) )
c--Ormel & Cuzzi, equation 26 (tightly coupled particles)
         Stokes_min = MIN(0.9,epsilon)*Stokes_max
         vrel_turb2 = vg2*vg2_coeff_disc_or_env*
     &        ABS( (Stokes_max-Stokes_min)*
     &        (Stokes_max**2*(1.0/(Stokes_max+sqrt_inv_Re) - 
     &        1.0/(1.0+Stokes_max))-
     &        Stokes_min**2*(1.0/(Stokes_min+sqrt_inv_Re) - 
     &        1.0/(1.0+Stokes_min)))/
     &        (Stokes_max+Stokes_min) )

c
c--Ormel & Cuzzi, equation 28 (fully intermediate regime, needed for
c     1.6*Stokes>Re^(-1/2) )
c
         IF (1.6*Stokes_max.GT.sqrt_inv_Re) THEN
c            vrel_turb2 = vrel_turb2 + vg2*vg2_coeff_disc_or_env*
            vrel_turb2 = vg2*vg2_coeff_disc_or_env*
     &           ( 2.2 - epsilon + 2.0/(1.0+epsilon)*
     &           (1.0/2.6 + epsilon**3/(1.6+ epsilon)))*Stokes_max
         ENDIF

      ELSE
         vrel_turb2 = 0.
      ENDIF
c
c--Turn off Brownian when only want to determine grain speed relative
c     to gas (for estimating drift as a fraction of h).  Brownian
c     velocity is largest for smallest particles, which are otherwise
c     stuck to the gas.
c
      IF (.NOT.iBrownian) vrel_brownian2 = 0.
c
c--Allow to turn off different components
c
c      vrel_brownian2 = 0.
c      vrel_terminal = 0.
c      vrel_radial_disc = 0.
c      vrel_settling = 0.
c      vrel_turb2 = 0.

      HY09_vreldust = SQRT(vrel_brownian2 + vrel_terminal**2 +
     &     vrel_radial_disc**2 + vrel_settling**2 + vrel_turb2)

      RETURN
      END FUNCTION HY09_vreldust

c----------------------------------------------------------- 

      FUNCTION HY09_vcoag(k,j)
c
c--Critical relative velocity for coagulation.  In cgs units.
c
      IMPLICIT NONE

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/HY09dustprops'

      REAL*8 HY09_vcoag
      INTEGER k,j

      HY09_vcoag =  HY09_vcoag_coeff*
     &     SQRT((HY09binsizes(k)**3 + HY09binsizes(j)**3)/
     &     (HY09binsizes(k) + HY09binsizes(j))**3) *
     &     (HY09binsizes(k)*HY09binsizes(j)/
     &     (HY09binsizes(k) + HY09binsizes(j)))**(-5.0/6.0)

c      HY09_vcoag = 1.0E+10

      RETURN
      END FUNCTION HY09_vcoag

c-----------------------------------------------------------

      FUNCTION HY09_cross_sec(k,j)
c
c--Collisional cross section
c
      IMPLICIT NONE

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/HY09dustprops'

      REAL HY09_cross_sec
      INTEGER k,j

      HY09_cross_sec = pi*(HY09binsizes(k) + HY09binsizes(j))**2

      RETURN
      END FUNCTION HY09_cross_sec

