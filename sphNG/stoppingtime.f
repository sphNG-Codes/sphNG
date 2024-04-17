      SUBROUTINE get_ts(ts,sgrain,densgrain,rhogas,rhodust,
     &                 cs_gas,abs_dv,iregime)

c***********************************************************************
c                                                                      *
c    This subroutine calculates the stopping time for gas + dust       *
c    unifying all of the previous drag prescriptions:                  *
c                 (1) Ayliffe & Bate (2012) two-fluid drag             *
c                 (2) Loren-Aguilar & Bate (2014) implicit drag        *
c                 (3) Price & Laibe (2015) one-fluid drag              *
c                                                                      *
c***********************************************************************
 
      IMPLICIT NONE 

      INCLUDE 'idim'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/tstop'
c
c--I/O variables
c
      INTEGER, INTENT(OUT) :: iregime
      REAL,    INTENT(OUT) :: ts
      REAL,    INTENT(IN)  :: sgrain,densgrain,abs_dv
      REAL,    INTENT(IN)  :: rhogas,rhodust,cs_gas
c
c--Local variables
c
      REAL :: tol_super
      REAL :: rhosum,kwok
      REAL :: lambda,kn_eff,viscmol_nu,Re_dust
      REAL :: mach,knudsen
      REAL :: Cdrag,dragcoeff,f,ts1
      REAL :: dv2,vratio,vratio2
      
      ! initialise variables
      tol_super  = 0.1
      dragcoeff  = 0.
      f          = 0.
      ts1        = 0.
      ts         = 0.
      rhosum     = rhogas + rhodust

!*******************************************************
      ! compute quantities specific to the drag regime
      SELECT CASE(idragscheme)
      CASE(1)
!*******************************************************
!  Laibe & Price (2012b)  &  Price & Laibe (2015)
!       (two-fluid)              (one-fluid)
!  physical drag (Epstein and Stokes regime)
!*******************************************************
         lambda = seff/rhogas
         IF (sgrain > 0.) THEN
            kn_eff = 9.*lambda/(4.*sgrain)
         ELSE
            kn_eff = HUGE(kn_eff)
         ENDIF
         IF (kn_eff >= 1.) THEN
            !
            ! Epstein regime
            !
            IF (densgrain > 0.) THEN
               dragcoeff = coeff_gei_1*cs_gas/(densgrain*sgrain)
            ELSE
               dragcoeff = HUGE(dragcoeff) ! so get ts=0 in this case
            ENDIF

            dv2 = abs_dv*abs_dv
            IF (cs_gas > 0. .and. dv2 > 0.) THEN
               kwok = 9.*pi/128.*dv2/(cs_gas*cs_gas)
               f = SQRT(1.+kwok)
            ELSE
               kwok = 0.
               f = 1. ! not important
            ENDIF
            iregime   = 1
            ! count where Kwok (1975) correction for supersonic drag 
            ! is important
            IF (kwok > tol_super) iregime = 2
         ELSE
            !
            ! Stokes regime
            !
            viscmol_nu = cste_mu*lambda*cs_gas  ! kinematic viscosity
            !--compute the local Stokes number
            Re_dust = 2.*sgrain*abs_dv/viscmol_nu
            IF (Re_dust  <=  1.) THEN
               dragcoeff = 4.5*viscmol_nu/(densgrain*sgrain*sgrain)
               f         = 1.
               iregime   = 3
            ELSEIF (Re_dust  <=  800.) THEN
               dragcoeff = 9./(densgrain*sgrain*Re_dust**0.6)
               f         = abs_dv
               iregime   = 4
            ELSE
               ! coeff is (3/8)*24/800**0.6
               dragcoeff = 0.163075/(densgrain*sgrain)
               f         = abs_dv
               iregime   = 5
            ENDIF
         ENDIF

      CASE(2)
!*******************************************************
!  Perets & Murray-Cley 2011
!  physical drag (Stokes regime)
!*******************************************************
         
         ! lambda = 1.619E-15/rho_g
         ! Constant using lambda defined in Perets & Murray-Clay 2011
         lambda     = 6.0E-15/rhogas
         Re_dust    = 2.*sgrain*abs_dv/(0.5*coeff_gei_1*cs_gas*lambda)
         Cdrag      = (24./Re_dust)*(1.+0.27*Re_dust)**0.43 + 
     &                0.47*(1.-EXP(-0.04*Re_dust**0.38))
         f          = abs_dv
         dragcoeff  = 3.*Cdrag/(8.*sgrain*densgrain)
         IF (Re_dust <= 1.) THEN
            iregime = 3
         ELSEIF (Re_dust <= 800.) THEN
            iregime = 4
         ELSE
            iregime = 5
         ENDIF

      CASE(3)
!*******************************************************
!  Baines et al. 1968 + Brasser et al. 2007
!  physical drag (Epstein and Stokes regime)
!*******************************************************

         !--Constant from Stepinski & Valageas 1996 (I).
         lambda = 1.619E-15/rhogas
         IF (sgrain > 0.) THEN
            kn_eff = 9.*lambda/(4.*sgrain)
         ELSE
            kn_eff = HUGE(kn_eff)
         ENDIF

         IF (kn_eff >= 1.) THEN
            !
            ! Epstein regime
            !
            iregime = 1
            vratio  = sqrt(gamma/2.)*abs_dv/cs_gas
            vratio2 = vratio**2

            !--Factor of 8 to match transition to Stokes...I think
            !  other coefficient is 1./(2*sqrt(pi))
            Cdrag = 8.*2.820948E-1*((1./vratio+1./(2.*vratio2*vratio))*
     &              exp(-vratio2) + (1. + 1./vratio2 - 1./
     &              (4.*vratio2*vratio2))*erf(vratio)*sqrt(pi))
 
         ELSE
            !
            ! Stokes regime
            !
            
            !
            !--Equation 19, Brasser et al. 2007
            !  Coefficient divided by 10 since Brasser uses SI units.
            !
            knudsen = (1.67E-9*udist**2/umass)/(rhogas*sgrain)
            !
            !--Equation 14, Brasser et al. 2007
            !
            mach = abs_dv/cs_gas
            !
            !--Equation 18, Brasser et al. 2007
            !  Using gamma = 5/3, rather than the 7/5 assumed by Brasser
            !
            Re_dust = 4.85*mach/knudsen
            !
            !--Set CD depending upon regime defined by above quantities.
            !
            IF (mach.GE.1.0) THEN
               Cdrag   = 2.0
               iregime = 6
            ELSEIF (mach.LT.1.0 .AND. Re_dust.GE.1.0E3) THEN
               Cdrag   = 0.44 + 1.56*mach**2
               iregime = 5
            ELSEIF (Re_dust.LT.1.0E3) THEN
               Cdrag   = 2.*mach**2 + (1. - mach**2)*(24.*(1. +
     &                   0.15*Re_dust**0.687)/Re_dust)
               iregime = 4
            ELSEIF (knudsen.GE.1.0) THEN
               !--If knudsen > 1, Brasser's formulation gives rho_g = 0.
               !  This means the gas drag should be zero, hence CD = 0.
               Cdrag   = 0.0
               iregime = 0
            ENDIF

         ENDIF

         f         = abs_dv
         dragcoeff = 3.*Cdrag/(8.*sgrain*densgrain)

      CASE(4)
!*******************************************************
! constant drag coefficient
!*******************************************************
         IF (K_code > 0. .and. rhogas > 0. .and. rhodust > 0.) THEN
            ! WARNING! When ndusttypes > 1, K_code ONLY makes sense
            ! if all of the grains are identical to one another.
            dragcoeff = K_code/(rhogas*rhodust)
         ELSE
            dragcoeff = tiny
         ENDIF
         f       = 1.
         iregime = 0

      CASE(5)
!*******************************************************
! constant ts
!*******************************************************
         dragcoeff = 1./(K_code*rhosum)
         f         = 1.
         iregime   = 0
      
      CASE DEFAULT
         !ts = 0.
         !iregime = 0 ! unknown
         WRITE (iprint,*) 'Invalid idrag selected ',idragscheme
         CALL quit(0) 
      END SELECT
!*******************************************************

      !--Calculate the stopping time
      IF (dragcoeff >= HUGE(dragcoeff)) THEN
         ts1 = HUGE(ts1)
      ELSE
         ts1 = dragcoeff*f*rhosum
      ENDIF
      IF (ts1 > 0.) THEN
         ts  = 1./ts1
      ELSE
         ts = HUGE(ts)
      ENDIF

      RETURN
      END SUBROUTINE get_ts
