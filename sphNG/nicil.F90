!----------------------------------------------------------------------!
!                               N I C I L                              !                                  
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!                                                                      !
! This is a stand-alone library that will calculate ionisation values  !
! and the coefficients for the non-ideal MHD terms: Ohmic resistivity, !
! Hall Effect and Ambipolar diffusion.                                 !          
!                                                                      !
!                   Copyright (c) 2015 James Wurster                   !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
!  MODULE: nicil
!          
!  DESCRIPTION:
!  Contains routines to calculation the ionisation rate, and related
!  useful quantities.
!  Copyright (c) 2015 James Wurster
!  See LICENCE file for usage and distribution conditions
!
!
!  REFERENCES: Asplund et al. (2009)
!              Fujii et al. (2011)
!              Keith & Wardle (2014)
!              Liu et al. (2003)
!              Nakano, Nishi & Umebayashi (2002)
!              Pandey & Wardle (2008)
!              Pinto & Galli (2008)
!              Pollack et al. (1994)
!              Umebayashi & Nakano (2009)
!              Umebayashi & Nakano (1990)
!              Wardle (2007)
!              Wardle & Ng (1999)
!
!  AUTHOR: James Wurster
!
!  PRIMARY RUNTIME PARAMETERS:
!    use_ohm        -- Calcualate and use the coefficient for Ohmic resistivity
!    use_hall       -- Calcualate and use the coefficient for the Hall effect
!    use_ambi       -- Calcualate and use the coefficient for ambipolar diffusion
!    g_cnst         -- Set constant grain size (true) or approximate an MRN grain distribution (false)
!    ion_rays       -- Include ionisation from cosmic (or x-) rays
!    ion_thermal    -- Include thermal ionisation
!    use_metalXY    -- set metallicity of H and He (true) or set abundances of 5 elements (false)
!    eta_constant   -- use a constant resistivity
!    fdg            -- Grain Parameter: gas to dust mass ratio
!    a0_grain        -- Grain Parameter if g_cnst=true : grain radius for constant grain size
!    an_grain       -- Grain Parameter if g_cnst=false: minimum grain radius for ~MRN distribution
!    ax_grain       -- Grain Parameter if g_cnst=false: maximum grain radius for ~MRN distribution
!    rho_bulk       -- Grain Parameter: bulk grain density
!    zeta_cgs       -- Ionisation rate
!    mass_MionR_mp  -- Mass of a metallic ion (for cosmic ray ionisation)
!    metal_X        -- If use_metalXY=true: Mass fraction of Hydrogen
!    metal_Y        -- If use_metalXY=true: Mass fraction of Helium
!    eta_const_calc -- If eta_constant=true: calculate eta using fixed environmental parameters
!    C_OR           -- If eta_constant=true & eta_constant_calc=false: eta_OR = C_OR
!    C_HE           -- If eta_constant=true & eta_constant_calc=false: eta_HE = C_HE*B
!    C_AD           -- If eta_constant=true & eta_constant_calc=false: eta_AD = C_AD*v_A^2
!    n_e_cnst       -- If eta_constant=true & eta_constant_calc=true: Constant electron number density
!    rho_i_cnst     -- If eta_constant=true & eta_constant_calc=true: Density of ionised gas
!    gamma_AD       -- If eta_constant=true & eta_constant_calc=true: Collisional coupling coefficient
!    hall_lt_zero   -- If eta_constant=true & eta_constant_calc=true: The sign of the Hall coefficient
!
!----------------------------------------------------------------------!
module nicil
 implicit none
 !--INPUT PAREMETERS
 !--Turn on/off individual non-ideal MHD coefficients
 logical,       public  :: use_ohm           = .true.           ! Calculate the coefficient for Ohmic resistivity
 logical,       public  :: use_hall          = .true.           ! Calculate the coefficient for the Hall effect
 logical,       public  :: use_ambi          = .true.           ! Calculate the coefficient for ambipolar diffusion
 !--Set constant grain size (g_cnst=true) or approximate an MRN grain distribution (g_cnst=false)
 logical,       public  :: g_cnst            = .true.
 !--Set ionisation type (can have one or both = true)
 logical,       public  :: ion_rays          = .true.           ! Include ionisation from cosmic (or x-) rays (=true)
 logical,       public  :: ion_thermal       = .true.           ! Include thermal ionisation (=true)
  !--Use metallicity (.true.) or abundances (.false.)
 logical,       public  :: use_metalXY       = .false.
 !--Assume that all the Hydrogen is molecular Hydrogen (used only for calculating mean molecular mass)
 logical,       public  :: H_is_H2           = .true.
 !--Use constant resistivity coefficients for all three resistivity terms
 logical,       public  :: eta_constant      = .false.
 !--Use the modified Hall parameters
 logical,       public  :: mod_beta          = .true.
 !
 !--Grain properties
 real,          public  :: fdg               = 0.01             ! gas to dust mass ratio
 real,          public  :: a0_grain          = 1.0d-5           ! grain radius for constant grain size [cm]
 real,          public  :: an_grain          = 5.0d-7           ! minimum grain radius for ~MRN distribution [cm]
 real,          public  :: ax_grain          = 2.5d-5           ! maximum grain radius for ~MRN distribution [cm]
 real,          public  :: rho_bulk          = 3.0              ! bulk grain density [g/cm^3]
 integer,    parameter  :: na_max            = 40               ! number of bins of grain size for MRN distribution

 !--Cosmic ray ionisation
 integer,public,parameter :: nimass          =  2               ! Number of ion masses for cosmic ray ionisation
 real,          public  :: zeta_cgs          =  1.0d-17         ! ionisation rate [s^-1]
 real,          public  :: mass_MionR_mp     = 24.3             ! mass of ion (default is mass of magnesium) [m_proton]
 ! 
 real,          public  :: delta_gn          = 1.3              ! Ionisation parameter: multiplicative factor for sigmavgnbyT 
 real,          public  :: sigmav_eH2        = 3.16d-11         ! Ionisation parameter: momentum transfer coefficient for e-H2 [cm^3/s]
 real,          public  :: sigmav_eHe        = 7.08d-11         ! Ionisation parameter: momentum transfer coefficient for e-He [cm^3/s]
 real,          public  :: pnH2              = 0.804            ! Ionisation parameter: polarizability for H2 [angstroms^3] 
 real,          public  :: pnHe              = 0.207            ! Ionisation parameter: polarizability for He [angstroms^3]
 !
 !--Metallicities (used if use_metalXY=.true.)
 real,        public    :: metal_X           =     0.70         ! Mass fraction of Hydrogen
 real,        public    :: metal_Y           =     0.28         ! Mass fraction of Helium
 !--Thermal ionisation
 integer,public,parameter :: nelements_max   =     5            ! Maximum number of elements
 integer,public,parameter :: nlevels         =     2            ! Number of calculated ionisation levels
 integer,     parameter :: n_nuin            =  1000            ! Number of ni_inT_coef values to pre-caculate
 real,        parameter :: se                =    1.0           ! Electron ticking coefficient: se \in (1.0d-3, 1.0)
 real,        parameter :: m_min             =    1.0           ! minimum ion mass (units of m_proton) in the table
 real,        parameter :: m_max             =  101.0           ! maximum ion mass (units of m_proton) in the table
 !--Resistivity coefficients (if fixed as constants)
 logical,     public    :: eta_const_calc    = .false.          ! Calculate constant coefficients using physical parameters (F = user set)
 real,        public    :: C_OR              =  0.1             ! eta_OR = C_OR       if eta_const_calc = .false.
 real,        public    :: C_HE              = -0.5             ! eta_HE = C_HE*B     if eta_const_calc = .false.
 real,        public    :: C_AD              =  0.01            ! eta_AD = C_AD*v_A^2 if eta_const_calc = .false.
 real,        public    :: n_e_cnst          = 1.0d19           ! Constant electron number density  [cm^-3]      if eta_const_calc = .true.
 real,        public    :: rho_i_cnst        = 3.8d-11          ! Density of ionised gas            [g/cm^3]     if eta_const_calc = .true.
 real,        public    :: gamma_AD          = 2.6e13           ! Collisional coupling coefficient  [cm^3/(s g)] if eta_const_calc = .true.
 logical,     public    :: hall_lt_zero      = .true.           ! The sign of the Hall coefficient               if eta_const_calc = .true.
 !
 !--Threshholds
 real,          public  :: T_thresh          = 1.0d9            ! Temperature [K] above which non-ideal MHD will turn off
 real,          public  :: ne_thresh_coef    = 1.0d-30          ! Minimum allowed ratio of n_e/n_n
 real,          public  :: Texp_thresh0      = 0.005            ! Will set exp(-chi/kT) = 0 is T is too low
 !
 !--Additional parameters
 real,        parameter :: epsilon_coef      = 100.0            ! for subtractions, will assume 0 if abs(a-b)<epsilon_coef*epsilon(a)
 real,        parameter :: NRtol             = 1.0e-12          ! tolerance on Newton–Raphson iterations
 integer,     parameter :: NRctrmax          = 2000             ! maximum number of Newton–Raphson before quiting
 real,        public    :: C_nimhd           = 0.1591549        ! Coefficient to control the timestep (==1/2pi)
 !--END OF INPUT PARAMETERS
 !
 !--Physical Constants (CGS)
 real,        parameter :: pi                =  3.1415926536d0  !  pi
 real,        parameter :: twopi             =  6.2831853072d0  ! 2pi
 real,        parameter :: fourpi            = 12.5663706144d0  ! 4pi
 real,        parameter :: c                 =  2.997924d10     ! Speed of light [cm/s]
 real,        parameter :: qe                =  4.8032068d-10   ! charge on electron [esu == statC]
 real,        parameter :: mass_electron_cgs =  9.10938291d-28  ! Electron mass [g]
 real,        parameter :: eV                =  1.60217657d-12  ! Electron volts [ergs]
 real,        parameter :: mass_proton_cgs   =  1.67262158d-24  ! Proton mass [g]
 real,        parameter :: kboltz            =  1.38066d-16     ! Boltzmann constant  [erg/K]
 real,        parameter :: planckh           =  6.6260755d-27   ! Planck's Constant [erg s]
 real,        parameter :: Amrn              =  1.5d-25         ! Constant for MRN grain distribution [cm^2.5]
 !
 !--Indicies for the various species; will allow for cleaner bookkeeping arrays
 integer,     parameter :: nspecies_max      =  6+na_max-1      ! The maximum number of charged species (do not modify)
 integer,     parameter :: ine               =  1               ! index for electron number density
 integer,     parameter :: iniHR             =  2               ! index for ion density of light elements (for CR ionisation)
 integer,     parameter :: iniMR             =  3               ! index for ion density of metallic elements (for CR ionisation)
 integer,     parameter :: inisT             =  4               ! index for ion density of singly ionised atoms (for thermal ionisation)
 integer,     parameter :: inidT             =  5               ! index for ion density of doubly ionised atoms (for thermal ionisation)
 integer,     parameter :: ing               =  6               ! index for the number density of the first grain
 !
 !--Local Variables to the NICIL module that may be required elsewhere in the user's code
 integer,     public    :: nelements
 real,        public    :: meanmolmass,unit_eta
 !--Local Variables to the NICIL module
 integer,     private   :: iprint,nspecies,na
 real,        private   :: csqbyfourpi
 real,        private   :: mass_proton,mump1,mass_neutral_cgs,mass_neutral1
 real,        private   :: nu_ei_coef,Saha_coef,n_total_grain_coef,ne_thresh
 real,        private   :: sigma_coef,sigmavenXbyT,sigmavenYbyT,sigmavenX_LR,sigmavenY_LR
 real,        private   :: eta_ohm_cnst,eta_hall_cnst,eta_ambi_cnst
 real,        private   :: zeta,Za_iter_coef,k_fac_coef,dm_ion,umass0,unit_density0
 real,        private   :: a_grain(na_max),a_grain_cgs(na_max),n_grain_coef(na_max)
 real,        private   :: kn_ig_coef(nimass,na_max),kn_eg_coef(na_max)
 real,        private   :: dkn_ig_coef(nimass,na_max),dkn_eg_coef(na_max)
 real,        private   :: n_grain_coef_cgs(na_max)
 real,        private   :: log_abunj(nelements_max),mj(nelements_max),abundancej(nelements_max)
 real,        private   :: mass_frac(nelements_max),sqrtmj1(nelements_max)
 real,        private   :: chij(nlevels,nelements_max),gej(nlevels,nelements_max),Texp_thresh(nlevels,nelements_max)
 real,        private   :: Zj0(nspecies_max),aZj0(nspecies_max),sZj0(nspecies_max)
 real,        private   :: massj(nspecies_max),beta_coef(nspecies_max),nu_jn_coef(nspecies_max)
 real,        private   :: nu_inT_coef(n_nuin)
 character(2),private   :: symj(nelements_max)
 !
 !--Subroutines
 public  :: nicil_initialise,nicil_get_ion_nZ,nicil_get_eta
 public  :: nimhd_get_jcbcb,nimhd_get_dBdt,nimhd_get_dudt,nimhd_get_dt
 public  :: nicil_translate_error
 private :: nicil_initialise_species,nicil_nuinTcoef,nicil_print_summary,nicil_ic_error,nicil_version
 private :: nicil_ion_get_sigma,nicil_ionR_get_Zgrain,nicil_ionR_get_n
 private :: nicil_ionT_get_ne,nicil_ionT_get_n,nicil_ionT_get_nj_Kjk,nicil_ionT_get_nion
 private :: nicil_nimhd_get_eta,nicil_nimhd_get_eta_cnst,nimhd_get_DcrossR
 !
 private
!
contains
!+
!----------------------------------------------------------------------!
!+
! Define the properties of the species that will be used for thermal 
! ionisation
!+
!----------------------------------------------------------------------!
subroutine nicil_initialise_species
 !  
 !--Initialise/Zero the values
 log_abunj     = 0.0
 chij          = huge(chij)         ! Default the ionisation potential to huge
 mj            = 0.0
 gej           = 2.0                ! Statistical weight of the electron
 if (use_metalXY) then
   nelements     = 2                ! Total number of elements
   !--Hydrogen
   symj(1)       = "H"              ! Chemical Symbol
   mj(1)         =  1.00            ! Mass [m_proton]
   chij(1,1)     = 13.60            ! First  ionisation potential [eV]
   mass_frac(1)  = metal_X
   !--Helium 
   symj(2)       = "He"             ! Chemical Symbol
   mj(2)         =  4.00            ! Mass [m_proton]
   chij(1,2)     = 24.59            ! First  ionisation potential [eV]
   chij(2,2)     = 54.42            ! Second ionisation potential [eV]
   mass_frac(2)  = metal_Y
 else
   nelements     = 5                ! Total number of elements
   !--Hydrogen
   symj(1)       = "H"              ! Chemical Symbol
   log_abunj(1)  = 12.0             ! Log abundance
   mj(1)         =  1.01            ! Mass [m_proton]
   chij(1,1)     = 13.60            ! First  ionisation potential [eV]
   !--Helium 
   symj(2)       = "He"             ! Chemical Symbol
   log_abunj(2)  = 10.93            ! Log abundance
   mj(2)         =  4.00            ! Mass [m_proton]
   chij(1,2)     = 24.59            ! First  ionisation potential [eV]
   chij(2,2)     = 54.42            ! Second ionisation potential [eV]
   !--Sodium 
   symj(3)       = "Na"             ! Chemical Symbol
   log_abunj(3)  =  6.24            ! Log abundance
   mj(3)         = 22.98            ! Mass [m_proton]
   chij(1,3)     =  5.14            ! First  ionisation potential [eV]
   chij(2,3)     = 47.29            ! Second ionisation potential [eV]
   !--Magnesium 
   symj(4)       = "Mg"             ! Chemical Symbol
   log_abunj(4)  =  7.60            ! Log abundance
   mj(4)         = 24.31            ! Mass [m_proton]
   chij(1,4)     =  7.65            ! First  ionisation potential [eV]
   chij(2,4)     = 15.03            ! Second ionisation potential [eV]
   !--Potassium 
   symj(5)       = "K"              ! Chemical Symbol
   log_abunj(5)  =  5.03            ! Log abundance
   mj(5)         = 39.10            ! Mass [m_proton]
   chij(1,5)     =  4.34            ! First  ionisation potential [eV]
   chij(2,5)     = 31.62            ! Second ionisation potential [eV]
 end if
 !
end subroutine nicil_initialise_species
!======================================================================!
! INITIALISATION & CONTROL ROUTINES                                    !
!======================================================================!
!----------------------------------------------------------------------!
!+
! Initialisation subroutine
! This will initialise all the variable required for NICIL, including 
! the frequently used coefficients.  
! All coefficients will be converted to code units.
!+
!----------------------------------------------------------------------!
subroutine nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint_in)
 implicit none
 real,              intent(in)  :: utime,umass,udist,unit_Bfield
 integer,           intent(out) :: ierr
 integer, optional, intent(in)  :: iprint_in
 integer                        :: j,k
 real                           :: unit_velocity,unit_density,unit_erg
 real                           :: sigmaviRn(nimass),sigmavgnbyT_coef
 real                           :: mass_proton_mp,mump
 real                           :: mass_neutral,massj_cgs(nspecies_max),massj_mp(nspecies_max)
 real                           :: a01_grain,dloga,n_denom_cgs,n_denom
 real                           :: vrms_mks,mu_iX,mu_iY,mu_eX,mu_eY
 real                           :: mass_total,n_rel,nj_rel(nelements_max)
 !
 !--Initialise species properties for thermal ionisation & Calculate abundances & mass fractions
 !  By construction, Sum (abundancej), Sum(mass_frac) == 1 
 call nicil_initialise_species
 !
 !--Verify input parameters are realistic; print error messages for each invalid error
 ierr = 0
 if (fdg  >= 1.0 .or. fdg < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid dust-to-gas fraction','fdg',fdg)
 if ( g_cnst ) then
   if (a0_grain  < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid fixed grain radius:','a_grain',a0_grain)
 else
   if (an_grain < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid minimum MRN grain radius:','an_grain',an_grain)
   if (ax_grain < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid maximum MRN grain radius:','ax_grain',ax_grain)
   if (ax_grain < an_grain) call nicil_ic_error(ierr,'nicil_initialise','Invalid radii ordering of MRN grain radii')
 end if
 if (rho_bulk < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid bulk grain density:','rho_bulk',rho_bulk)
 if (zeta_cgs < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid ionisation rate:'   ,'zeta'    ,zeta_cgs)
 if (mass_MionR_mp < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid ion mass:'     ,'m_Mion'  ,mass_MionR_mp)
 if (use_metalXY) then
    if (metal_X < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid Hydrogen mass fraction:','X',metal_X)
    if (metal_X < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid helium mass fraction:'  ,'Y',metal_Y)
    if (metal_X+metal_Y>1.0) &
                       call nicil_ic_error(ierr,'nicil_initialise','Invalid metalicities:'  ,'metal_X + metal_Y',metal_X+metal_Y)
 end if
 if (epsilon_coef <= 1.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid subtraction constraint','epsilon_coef',epsilon_coef)
 if (NRtol < 1.0e-15 .or. NRtol > 0.1) & 
                     call nicil_ic_error(ierr,'nicil_initialise','Poor constraint on Newton–Raphson tolerance','NRtol',NRtol)
 if (NRctrmax <= 10) call nicil_ic_error(ierr,'nicil_initialise','Too few maximum permitted Newton–Raphson iterations' &
 ,'NRctrmax',real(NRctrmax))
 if (nelements > nelements_max)  call nicil_ic_error(ierr,'nicil_initialise','nelements > nelements_max')
 if (iniHR > iniMR)  call nicil_ic_error(ierr,'nicil_initialise','indicies for light & heavy ion densities out of order')
 if (inisT > inidT)  call nicil_ic_error(ierr,'nicil_initialise','indicies for singly & doubly ionised densities out of order')
 if (trim(symj(1))/="H" ) & 
    call nicil_ic_error(ierr,'nicil_initialise_species','Hygrogen does not have the correct array number of 1')
 if (trim(symj(2))/="He") &
    call nicil_ic_error(ierr,'nicil_initialise_species','Helium does not have the correct array number of 2')
 !--Abort setup if errors fatal exist
 if (ierr/=0) return   
 !--If no ionisation or no non-ideal mhd coefficients, then set all to false
 !  (this is not a fatal error, since the user may require ideal mhd while not removing the nicil algorithm)
 if ((.not.ion_rays .and. .not.ion_thermal).or.(.not.use_ohm .and. .not.use_hall .and. .not.use_ambi)) then
   ion_rays    = .false.
   ion_thermal = .false.
   use_ohm     = .false.
   use_hall    = .false.
   use_ambi    = .false.
 end if
 !
 !--Calculate metallicities or abundances, as required
 if (use_metalXY) then
   ! Note: Abundance is approximate since we do not explicity accound for Z
   mass_total    = 1.0/(metal_X/mj(1) + metal_Y/mj(2))
   abundancej(1) = metal_X*(mass_total)/(mj(1))
   abundancej(2) = metal_Y*(mass_total)/(mj(2))
   massj_mp(iniHR)  = mass_total                            ! Light element ion mass
 else
   n_rel       = 0.0
   mass_total  = 0.0
   do j = 1,nelements
     nj_rel(j) = 10**(log_abunj(j) - 12.0)                  ! Number density relative to H
     n_rel     = n_rel + nj_rel(j)
   end do
   do j = 1,nelements
     abundancej(j) = nj_rel(j)/n_rel                        ! Abundance
     mass_total    = mass_total + mj(j)*abundancej(j)
   end do
   do j = 1,nelements
     mass_frac(j)  = mj(j)*abundancej(j)/mass_total
   end do
 end if
 !--Calculate the mean molar mass
 meanmolmass = 0.0
 sqrtmj1     = 0.0
 do j = 1,nelements
   if (H_is_H2 .and.  trim(symj(j))=="H") then 
     meanmolmass = meanmolmass + mass_frac(j)/(2.0*mj(j))
   else
     meanmolmass = meanmolmass + mass_frac(j)/mj(j)
   end if
   if (j==2) massj_mp(iniHR) = 1.0/meanmolmass              ! Light element ion mass
   sqrtmj1(j) = 1.0/sqrt(mj(j))
 end do
 meanmolmass = 1.0/meanmolmass
 !
 !--Unit conversions from cgs-> code
 unit_velocity     = udist/utime
 unit_density      = umass/udist**3
 unit_erg          = umass*unit_velocity**2
 unit_eta          = udist**2/utime 
 umass0            = umass
 unit_density0     = unit_density
 !
 !--Initialise masses (m_proton units)
 mass_proton_mp    = 1.0
 massj_mp(ine)     = mass_electron_cgs/mass_proton_cgs
 massj_mp(iniMR)   = mass_MionR_mp
 !
 !--Initialise variables (CGS variables)
 massj_cgs(ine)         = mass_electron_cgs                                  ! Electron mass
 massj_cgs(iniHR:iniMR) = massj_mp(iniHR:iniMR)*mass_proton_cgs              ! Ion mass for cosmic ray ionisation
 mass_neutral_cgs       = meanmolmass  *mass_proton_cgs                      ! Neutral mass
 mump                   = meanmolmass  *mass_proton_cgs                      ! mean particle mass assuming composition is H2 & He
 if (.not. use_metalXY) then
   metal_X         = mass_frac(1)                                            ! Metallicity of Hydrogen
   metal_Y         = mass_frac(2)                                            ! Metallicity of Helium
 end if
 !
 !--Determine grain distribution and initialise
 a_grain_cgs      = 0.0
 n_grain_coef_cgs = 0.0
 nspecies         = 5
 if (ion_rays) then
   if (g_cnst) then
     na = 1
     a_grain_cgs(1)      = a0_grain
     massj_cgs(ing)      = fourpi/3.0*a_grain_cgs(1)**3*rho_bulk
     n_grain_coef_cgs(1) = mass_neutral_cgs*fdg/massj_cgs(ing)
     n_denom_cgs         = n_grain_coef_cgs(1)*a_grain_cgs(1)
   else
     na          = na_max
     dloga       = (log10(ax_grain) - log10(an_grain))/float(na)
     n_denom_cgs = 0.0
     do j = 1,na
       a_grain_cgs(j)      = 10**(log10(an_grain) +(j-0.5)*dloga )
       massj_cgs(ing+j-1)  = fourpi/3.0*a_grain_cgs(j)**3*rho_bulk
       n_grain_coef_cgs(j) = 0.4*Amrn*(( 10**(log10(an_grain) +(j-1)*dloga ) )**(-2.5) &
                           - (10**(log10(an_grain) +(j)*dloga ))**(-2.5))
       n_denom_cgs         = n_denom_cgs + n_grain_coef_cgs(j)*a_grain_cgs(j)
     end do
   end if
   massj_mp(ing:ing+na-1) = massj_cgs(ing:ing+na-1)/mass_proton_cgs
   nspecies = nspecies + na
 end if
 !
 !--Initialise variables (Code unit variables)
 csqbyfourpi       = c**2/fourpi              / unit_eta
 zeta              = zeta_cgs                 * utime
 mass_proton       = mass_proton_cgs          / umass
 mass_neutral      = mass_neutral_cgs         / umass
 massj             = massj_cgs                / umass
 a_grain           = a_grain_cgs              / udist
 n_grain_coef      = n_grain_coef_cgs
 !
 !--Charge capture Rates coefficients (Code units)
 k_fac_coef         = qe**2/kboltz  / udist
 do j = 1,na
   do k = 1,nimass
     kn_ig_coef(k,j)  = a_grain_cgs(j)**2*sqrt(8.0*pi*kboltz/massj_cgs(k-1+iniHR) )*n_grain_coef(j) * utime/udist**3
     dkn_ig_coef(k,j) = kn_ig_coef(k,j)*k_fac_coef
   end do
   kn_eg_coef(j)      = a_grain_cgs(j)**2*sqrt(8.0*pi*kboltz/mass_electron_cgs)*n_grain_coef(j)     * utime/udist**3
   dkn_eg_coef(j)     = kn_eg_coef(j)*k_fac_coef
 end do
 !
 !--Coefficients for grain charge calculation (Code units)
 n_denom            = n_denom_cgs   / udist
 n_total_grain_coef = 0.0
 do j = 1,na
   n_total_grain_coef = n_total_grain_coef + n_grain_coef(j)
 end do
 !
 !--Rate coefficientes (CGS; sigma units are cm^3/s)
 !  *_LR are Langevine rates, and the maximum sigmaenX is used at any given time
 vrms_mks      = sqrt(8.0*kboltz/(pi*mass_electron_cgs)) * 1.0d-5  ! km/s
 mu_eX         = mass_proton_cgs*(2.0*mj(1) * massj_mp(ine))        /(2.0*mj(1) + massj_mp(ine))
 mu_eY         = mass_proton_cgs*(    mj(2) * massj_mp(ine))        /(    mj(2) + massj_mp(ine))
 do j = 1,nimass
   mu_iX        = mass_proton_cgs*(2.0*mj(1) * massj_mp(j-1+iniHR)) /(2.0*mj(1) + massj_mp(j-1+iniHR))
   mu_iY        = mass_proton_cgs*(    mj(2) * massj_mp(j-1+iniHR)) /(    mj(2) + massj_mp(j-1+iniHR))
   sigmaviRn(j) = 2.81d-9*(metal_X*sqrt(pnH2 * mass_proton_cgs/mu_iX) + metal_Y*sqrt(pnHe * mass_proton_cgs/mu_iY) )
 end do
 sigmavenXbyT  = metal_X*sigmav_eH2*vrms_mks**1.3
 sigmavenYbyT  = metal_Y*sigmav_eHe*vrms_mks
 sigmavenX_LR  = metal_X*2.81d-9*sqrt(pnH2 * mass_proton_cgs/mu_eX)                
 sigmavenY_LR  = metal_Y*2.81d-9*sqrt(pnHe * mass_proton_cgs/mu_eY)                              
 sigmavgnbyT_coef = delta_gn * sqrt(128.0*pi*kboltz/(9.0*mass_neutral_cgs))
 !
 !--Collisional Frequencies (CGS = s^-1)
 nu_ei_coef              = 51.0                                                                    / udist**3        ! nu_ei
 nu_jn_coef(ine)         = 1.0       /(mass_neutral_cgs+massj_cgs(ine)  )                          * unit_density    ! nu_en
 nu_jn_coef(iniHR:iniMR) = sigmaviRn /(mass_neutral_cgs+massj_cgs(iniHR:iniMR))                    * unit_density    ! nu_in
 do j = 1,na
   nu_jn_coef(ing+j-1) = sigmavgnbyT_coef*a_grain_cgs(j)**2/(mass_neutral_cgs+massj_cgs(ing+j-1) ) * unit_density    ! nu_gn
 end do
 if (ion_thermal) then
   call nicil_nuinTcoef
 else
   nu_inT_coef = 0.0
 end if
 !
 !--Fill charge arrays (electric charge,Z; absolute value of Z; sign of Z)
 Zj0(ine)   = -1.0; aZj0(ine)   = 1.0; sZj0(ine)   = -1.0        ! electrons
 Zj0(iniHR) =  1.0; aZj0(iniHR) = 1.0; sZj0(iniHR) =  1.0        ! ions of light elements (for CR ionisation)
 Zj0(iniMR) =  1.0; aZj0(iniMR) = 1.0; sZj0(iniMR) =  1.0        ! ions of metallic elements (for CR ionisation)
 Zj0(inisT) =  1.0; aZj0(inisT) = 1.0; sZj0(inisT) =  1.0        ! ions of singly ionised atoms (for thermal ionisation)
 Zj0(inidT) =  2.0; aZj0(inidT) = 2.0; sZj0(inidT) =  1.0        ! ions of doubly ionised atoms (for thermal ionisation)
 !
 !--Hall parameters (dimensionless after multipiled by B_code (and mass_ionT_code))
 beta_coef(ine)          = qe/(massj_cgs(ine)          *c) * unit_Bfield
 beta_coef(iniHR:iniMR)  = qe/(massj_cgs(iniHR:iniMR)  *c) * unit_Bfield
 beta_coef(inisT:inidT)  = qe/(                         c) * unit_Bfield    /umass
 beta_coef(ing:ing+na-1) = qe/(massj_cgs(ing:ing+na-1) *c) * unit_Bfield
 !
 !--Conductivity coefficient (code units)
 sigma_coef   = qe*c                                       / (udist**3 * unit_Bfield)
 !
 !--Particle Masses (Code units)
 mass_neutral1= 1.0/mass_neutral
 mump1        = 1.0/mump                                   * umass
 ne_thresh    = ne_thresh_coef*mump1
 !
 !--Grain charge coefficient (Code units)
 if (ion_rays) then
   Za_iter_coef  = zeta/(n_denom*mump1)
 end if
 !
 !--Coefficient for the Saha equation (Code units)
 Saha_coef   = (2.0*pi*massj(ine)*kboltz/planckh**2 * utime**2*unit_erg )**1.5
 chij        = chij*eV/kboltz
 Texp_thresh = chij*Texp_thresh0
 !
 !--Set constant coefficients
 if ( eta_const_calc ) then
   eta_ohm_cnst  = mass_electron_cgs*c**2/(fourpi*qe**2*n_e_cnst)
   eta_hall_cnst = c/(fourpi*qe*n_e_cnst)
   eta_ambi_cnst = 1.0/(fourpi*gamma_AD*rho_i_cnst)
   if (hall_lt_zero) eta_hall_cnst = -eta_hall_cnst
 else
   eta_ohm_cnst  = C_OR
   eta_hall_cnst = C_HE / sqrt(fourpi)
   eta_ambi_cnst = C_AD /      fourpi
 end if
 !  Convert units as required
 eta_ohm_cnst  = eta_ohm_cnst                                  / unit_eta
 eta_hall_cnst = eta_hall_cnst *  unit_Bfield                  / unit_eta
 eta_ambi_cnst = eta_ambi_cnst * (unit_Bfield**2/unit_density) / unit_eta
 !
 !--Set the location of the printouts
 if (present(iprint_in)) then
   iprint = iprint_in
 else
   iprint = 6
 end if
 !
 !--Print Statements to summarise conditions used
 call nicil_print_summary(a01_grain,massj_mp(iniHR),massj_mp(ing),meanmolmass)
 !
 return
end subroutine nicil_initialise
!----------------------------------------------------------------------!
!+
! Routines to calculate nu_inT_coef tables, and extract the relevant 
! value
!+
!----------------------------------------------------------------------!
!--Control subroutine to make the table
subroutine nicil_nuinTcoef
 implicit none
 integer            :: i
 real               :: mion
 !
 dm_ion  = (m_max-m_min)/float(n_nuin)
 do i = 1,n_nuin
   mion = ( m_min + float(i-1)*dm_ion)*mass_proton
   call nicil_calc_nuinTcoef(nu_inT_coef(i),mion)
 end do
!
end subroutine nicil_nuinTcoef
!--Subroutine to make the individual entries
pure subroutine nicil_calc_nuinTcoef(nu_inT_coef0,mion)
 implicit none
 real,    intent(in)  :: mion
 real,    intent(out) :: nu_inT_coef0
 real                 :: mion_mp,mu_iX,mu_iY,sigmaviTn
 mion_mp      = mion/mass_proton
 mu_iX        = mass_proton_cgs*(2.0*mj(1) * mion_mp)/(2.0*mj(1) + mion_mp)
 mu_iY        = mass_proton_cgs*(    mj(2) * mion_mp)/(    mj(2) + mion_mp)
 sigmaviTn    = 2.81d-9*(metal_X*sqrt(pnH2 * mass_proton_cgs/mu_iX) + metal_Y*sqrt(pnHe * mass_proton_cgs/mu_iY) )
 nu_inT_coef0 = sigmaviTn   /(mass_neutral_cgs+mion*umass0) * unit_density0
end subroutine nicil_calc_nuinTcoef
!
!--Subroutine to extract a coefficient from the table
pure subroutine nicil_extract_nuinTcoef(nu_inT_coef0,mion)
 implicit none
 real,    intent(in)  :: mion
 real,    intent(out) :: nu_inT_coef0
 integer              :: ibin
 real                 :: rbin,fbin,mion_mp
 !
 mion_mp = mion/mass_proton
 rbin    = (mion_mp - m_min)/dm_ion + 1.0
 ibin    = int(rbin)
 if (ibin < 1 .or. ibin > n_nuin-1) then
  call nicil_calc_nuinTcoef(nu_inT_coef0,mion)
 else
   fbin         = rbin - ibin
   nu_inT_coef0 = nu_inT_coef(ibin) + fbin*(nu_inT_coef(ibin+1) - nu_inT_coef(ibin))
 end if
end subroutine nicil_extract_nuinTcoef
!----------------------------------------------------------------------!
!+
! Print Statements to summarise what is being used here
!+
!----------------------------------------------------------------------!
subroutine nicil_print_summary(a01_grain,mass_HionR_mp,mass_grain_mp,meanmolmass)
 implicit none
 real, intent(in) :: a01_grain,mass_HionR_mp,mass_grain_mp,meanmolmass
 integer          :: j
 character(200)   :: ni_terms,version
 !
 call nicil_version(version)
 write(iprint,'(2a)') "NICIL: ",trim(version)
 write(iprint,'(a)' ) "NICIL: Copyright (c) 2015 James Wurster"
 write(iprint,'(a)' ) "NICIL: See LICENCE file for usage and distribution conditions"

 !
 if(.not.use_ohm.and. .not.use_hall .and. .not.use_ambi) ni_terms=""
 if(     use_ohm.and. .not.use_hall .and. .not.use_ambi) ni_terms="Ohmic Resistivity"
 if(.not.use_ohm.and.      use_hall .and. .not.use_ambi) ni_terms="Hall Effect"
 if(.not.use_ohm.and. .not.use_hall .and.      use_ambi) ni_terms="Ambipolar Diffusion"
 if(     use_ohm.and.      use_hall .and. .not.use_ambi) ni_terms="Ohmic Resistivity, Hall Effect"
 if(     use_ohm.and. .not.use_hall .and.      use_ambi) ni_terms="Ohmic Resistivity, Ambipolar Diffusion"
 if(.not.use_ohm.and.      use_hall .and.      use_ambi) ni_terms="Hall Effect, Ambipolar Diffusion"
 if(     use_ohm.and.      use_hall .and.      use_ambi) ni_terms="Ohmic Resistivity, Hall Effect, Ambipolar Diffusion"
 !
 if (eta_constant) then
   write(iprint,'(2a)')                       "NICIL: Non-ideal terms used: ",trim(ni_terms)
   write(iprint,'(a)' )                       "NICIL: All resistivity coefficients are constant."
   if (use_ohm ) write(iprint,'(a,Es10.3  )') "NICIL: Eta_ohm  = ", eta_ohm_cnst
   if (use_hall) write(iprint,'(a,Es10.3,a)') "NICIL: Eta_Hall = ", eta_hall_cnst,"*|B|"
   if (use_ambi) write(iprint,'(a,Es10.3,a)') "NICIL: Eta_ambi = ", eta_ambi_cnst,"*|B|^2/rho"
 else
   if (ion_rays) then
     if (g_cnst) then 
       write(iprint,'(a)')           "NICIL: Using constant grain size."
     else
       write(iprint,'(a)')           "NICIL: Approximating grain distribution with MRN distribution."
     end if
   end if
   write(iprint,'(a)')             "NICIL: Species  abundance   mass fraction"
   do j = 1,nelements
     if (mass_frac(j) > 0.01) then 
       write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),"     ",abundancej(j),mass_frac(j)
     else
       write(iprint,'(3a,Es10.3,4x,Es10.3)') "NICIL: ",symj(j),"     ",abundancej(j),mass_frac(j)
     end if
   end do
   if (H_is_H2) then 
     write(iprint,'(a,Es10.3,a)')  "NICIL: Mean molecular mass:             ",meanmolmass," m_proton  (assuming molecular hydrogen)"
   else
     write(iprint,'(a,Es10.3,a)')  "NICIL: Mean molecular mass:             ",meanmolmass," m_proton  (assuming atomic hydrogen)"
   end if
   if (ion_rays) then
     write(iprint,'(a,Es10.3,a)')    "NICIL: Light element ion mass for CR's: ",mass_HionR_mp," m_proton"
     write(iprint,'(a,Es10.3,a)')    "NICIL: Metallic ion mass for CR's:      ",mass_MionR_mp," m_proton"
     if (g_cnst) then 
       write(iprint,'(a,Es10.3,a)')  "NICIL: Grain mass:                      ",mass_grain_mp," m_proton"
       write(iprint,'(a,Es10.3,a)')  "NICIL: Grain radius:                    ",a0_grain," cm"
     else
       write(iprint,'(a,Es10.3,a)')  "NICIL: Average grain mass:              ",mass_grain_mp/mass_proton_cgs," m_proton"
       write(iprint,'(a,Es10.3,a)')  "NICIL: Average grain radius:            ",a01_grain," cm"
     end if
     write(iprint,'(a,Es10.3,a)')  "NICIL: Cosmic ray ionisation rate:      ",zeta_cgs," s^{-1}"
     write(iprint,'(a)')           "NICIL: Including ionisation from Cosmic rays"
   end if
   if (ion_thermal) then
     write(iprint,'(a)')           "NICIL: Including thermal ionisation"
   end if
   if(use_ohm.or.use_hall .or.use_ambi) then
     write(iprint,'(2a)')          "NICIL: Non-ideal terms used: ",trim(ni_terms)
   end if
   if (.not.ion_thermal .and. .not.ion_rays) then
     write(iprint,'(a)')           "NICIL: WARNING! No ionisation sources included!"
     write(iprint,'(a)')           "NICIL: WARNING! No non-ideal MHD coefficients will be calculated!"
     write(iprint,'(a)')           "NICIL: WARNING! This is Ideal MHD!"
   end if
 end if
 !
end subroutine nicil_print_summary

!----------------------------------------------------------------------!
!+
! Theres will print an error message and increment a counter for each
! error found in NICIL. 
! This will NOT end the main programme, but pass the error code such
! that the code can be properly and cleanly terminated.
! The first programme is specifically for initialising NICIL, and the
! second is during runtime.
!+
!----------------------------------------------------------------------!
subroutine nicil_ic_error(num_errors,wherefrom,reason,var,val)
 implicit none
 integer,                    intent(inout) :: num_errors
 character(len=*),           intent(in)    :: wherefrom,reason
 character(len=*), optional, intent(in)    :: var
 real,             optional, intent(in)    :: val
 !
 !--Print cause of  error
 if (present(var) .and. present(val)) then 
   write(iprint,'(7a,Es16.4)') 'NICIL: ERROR: ',wherefrom,': ',reason,'. ',var,'=',val
 else
   write(iprint,'(4a)') 'NICIL: ERROR: ',wherefrom,': ',reason
 end if
 num_errors = num_errors + 1
 !
end subroutine nicil_ic_error
!
subroutine nicil_translate_error(ierr)
 implicit none
 integer,                    intent(in)    :: ierr
 !
 if (ierr== 1.or.ierr== 3.or.ierr== 5.or.ierr== 7.or.ierr== 9.or.ierr==11.or.ierr==13.or.ierr==15) &
   write(iprint,'(a)') 'NICIL: ERROR: nicil_ion_get_Zgrain: Zg_on_ag did not converge'
 if (ierr== 2.or.ierr== 3.or.ierr== 6.or.ierr== 7.or.ierr==10.or.ierr==11.or.ierr==14.or.ierr==15) &
   write(iprint,'(a)') 'NICIL: ERROR: nicil_ion_get_Zgrain: Zg_on_ag > 0 is not permitted with this physics'
 if (ierr== 4.or.ierr== 5.or.ierr== 6.or.ierr== 7.or.ierr==12.or.ierr==13.or.ierr==14.or.ierr==15) &
   write(iprint,'(a)') 'NICIL: ERROR: nicil_ionT_get_ne: n_electronT did not converge'
 if (ierr== 8.or.ierr== 9.or.ierr==10.or.ierr==11.or.ierr==12.or.ierr==13.or.ierr==14.or.ierr==15) &
   write(iprint,'(a)') 'NICIL: ERROR: nicil_ionT_get_ne: n_electronT < 0'
 !
 if (ierr==10.or.ierr==30) &
   write(iprint,'(a)') 'NICIL: ERROR: nicil_ionR_get_n: unrealistic n_elctronR value (i.e. n_e < 0)'
 if (ierr==20.or.ierr==30) &
   write(iprint,'(a)') 'NICIL: ERROR: nicil_ionR_get_n: unrealistic n_ionR value (i.e. n_i < 0)'
 !
end subroutine nicil_translate_error
!----------------------------------------------------------------------!
!+
! Internal Version Control
!+
!----------------------------------------------------------------------!
pure subroutine nicil_version(version)
 implicit none
 character(200), intent(out) :: version 
 !
 version = "Version 1.0: 8 Dec 2015: Initial Version"
 !
end subroutine nicil_version
!----------------------------------------------------------------------!
!+
! This is the primary control routine for NICIL.  It will call the 
! required ionisation and non-ideal MHD routines to determine the 
! coefficients of the non-ideal MHD terms
!+
!----------------------------------------------------------------------!
pure subroutine nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,T,Zg_on_ag,n_electronT,ierr,data_out)
 implicit none
 integer,        intent(out)   :: ierr
 real,           intent(out)   :: eta_ohm,eta_hall,eta_ambi
 real,           intent(in)    :: Bfield,rho,T,Zg_on_ag,n_electronT
 real, optional, intent(out)   :: data_out(13+nelements_max*nlevels-1)
 real                          :: n_electron,n_electronR
 real                          :: n_ionR(nimass),n_ionT(nlevels),mass_ionT
 real                          :: sigmaOHPpa(5),n_densities(9)
 real                          :: njk(nlevels,nelements)
 logical                       :: get_data_out
 !
 ierr = 0                                 ! initialise error code
 !--Determine if we will be returning bookkeeping data
 if (present(data_out)) then
    get_data_out = .true.
 else
    get_data_out = .false.
 end if
 !
 !--Return constant coefficient version and exit
 if (eta_constant) then
   call nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
   if (present(data_out)) data_out = 0.0  ! bookkeeping has no meaning here
   return
 end if
 !
 !--Calcualte the grain charge and electron number density from cosmic rays
 if (ion_rays) then
   call nicil_ionR_get_n(n_electronR,n_ionR,T,Zg_on_ag,ierr)
 else
   n_electronR = 0.0
   n_ionR      = 0.0
 end if
 !
 !--Calculate the grain charge from thermal ionisaion and update electron number density
 if (ion_thermal) then
   call nicil_ionT_get_n(n_ionT,mass_ionT,njk,n_electronT,T,rho)
 else
   n_ionT      = 0.0
   mass_ionT   = 0.0
   njk         = 0.0
 end if
 n_electron    = n_electronR + n_electronT
 !
 !--Calculate the conductivities
 call nicil_ion_get_sigma(Bfield,rho,n_electron,n_ionR,Zg_on_ag,n_ionT,mass_ionT,T, &
                         sigmaOHPpa,get_data_out,n_densities)
 !
 !--Calculate the coefficients
 call nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,T,sigmaOHPpa)
 !
 !--Copy optional arrays to output, if requested
 if (present(data_out)) then
   data_out( 1: 3) = sigmaOHPpa(1:3)  ! Ohmic conductivity, Hall conductivity, Pedersen conductivity
   data_out(    4) = n_densities(  1) ! Average grain charge
   data_out( 5: 6) = n_densities(2:3) ! rho_neutral, rho_ion (total)
   data_out(    7) = n_electron
   data_out( 8:13) = n_densities(4:9) ! n_neutral, n_ionR (light,metallic), n_ionT(singly,doubly), n_grainR
   data_out(14:14+nelements-1) = njk(1,:)
   data_out(14+nelements:14+nelements*nlevels-2) = njk(2,2:nelements)
 end if
 !
 return
end subroutine nicil_get_eta
!----------------------------------------------------------------------!
!+
!  calculates the condictivities
!  Terms have been modified such that they require input of sound speed 
!  rather than temperature.
!  sigmaOHPpa(1): Ohmic conductivity
!  sigmaOHPpa(2): Hall conductivity
!  sigmaOHPpa(3): Pedersen conductivity
!  sigmaOHPpa(4): inverse square of the perpendicular conductivity:
!                 1.0/(sigmaOHPpa(2)**2+sigmaOHPpa(3)**2)
!  sigmaOHPpa(5): rearrangement for the coefficient to Ambipolar
!                 diffusion (see notes in nimhd_get_eta)
!  The n_densities arrays is for bookkeeping only, and has no real 
!  computational value. 
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ion_get_sigma(Bmag,rho,n_electron,n_ionR,Zg_on_agi,n_ionT,mass_ionT,T, &
                                    sigmaOHPpa,get_data_out,n_densities)
 implicit none
 real,   intent(out) :: sigmaOHPpa(5),n_densities(9)
 real,   intent(in)  :: Bmag,rho,T,n_electron,Zg_on_agi,mass_ionT
 real,   intent(in)  :: n_ionR(:),n_ionT(:)
 logical,intent(in)  :: get_data_out
 integer             :: j,k
 real                :: sqrtT,sigmavenX,sigmavenY,sigma_coef_onB
 real                :: rho_n,nu_ei,nu_inT_coef0,sigmaHp,sigmaHn
 real                :: Zj(nspecies_max),aZj(nspecies_max),sZj(nspecies_max)
 real                :: nu_jn(nspecies_max),ns(nspecies_max),rho_j(nspecies_max)
 real                :: betaj(nspecies_max),beta2p11(nspecies_max)
 !
 !--Initialise values
 sqrtT          = sqrt(T)
 sigma_coef_onB = sigma_coef/Bmag
 sigmaOHPpa     = 1.0
 !
 !--Number densities
 ns(ine)          = n_electron
 ns(iniHR:iniMR)  = n_ionR
 ns(inisT:inidT)  = n_ionT
 ns(ing:ing+na-1) = n_grain_coef(1:na)*(rho*mump1)
 !
 !--Densities 
 rho_j(ine)         = n_electron*massj(ine)
 rho_j(iniHR:iniMR) = n_ionR *massj(iniHR:iniMR)
 rho_j(inisT:inidT) = n_ionT *mass_ionT
 rho_n              = rho - rho_j(ine)
 do j = 1,nimass
   rho_n            = rho - rho_j(j-1+iniHR)
 end do
 do k = 1,nlevels
   rho_n            = rho - rho_j(k-1+inisT)
 end do
 !
 if (rho_n < tiny(rho_n) .or. n_electron < rho_n*ne_thresh) then
   ! This is the Ideal MHD regieme.  Turn off non-ideal terms.  
   sigmaOHPpa  = -1.0
   if (rho_n < tiny(rho_n)) rho_n =  0.0
 else
   !
   !--Update the grain charges in the array
   Zj       = Zj0
   aZj      = aZj0
   sZj      = sZj0
   Zj( ing:ing+na-1) = Zg_on_agi*a_grain(1:na)
   aZj(ing:ing+na-1) = abs(      Zj(ing:ing+na-1) )
   sZj(ing:ing+na-1) = sign(1.0, Zj(ing:ing+na-1) )
   !
   !--Rate Coefficients; use Langevine rate if high temperature
   sigmavenX = sigmavenXbyT*sqrtT**1.3
   sigmavenY = sigmavenYbyT*sqrtT
   if (sigmavenX > sigmavenX_LR) sigmavenX = sigmavenX_LR
   if (sigmavenY > sigmavenX_LR) sigmavenY = sigmavenX_LR
   !
   !--Collisional Frequencies
   if (ion_thermal .and.  mass_ionT > 0.0) then
     call nicil_extract_nuinTcoef(nu_inT_coef0,mass_ionT)
   else
     nu_inT_coef0 = 0.0
   end if
   nu_ei               = nu_ei_coef                *n_electron/sqrtT**3
   nu_jn(ine)          = nu_jn_coef(ine)           *rho_n*(sigmavenX + sigmavenY)
   nu_jn(iniHR:iniMR)  = nu_jn_coef(iniHR:iniMR)   *rho_n
   nu_jn(inisT:inidT)  = nu_inT_coef0              *rho_n
   nu_jn(ing:ing+na-1) = nu_jn_coef(ing:ing+na-1)  *rho_n*sqrtT
   !
   !--Hall parameters
   betaj        = 0.0
   beta2p11     = 0.0
   betaj(ine)   = beta_coef(ine) *         Bmag/(nu_jn(ine)  + nu_ei)
   betaj(iniHR) = betaIj(beta_coef(iniHR),rho_j(iniHR),rho_j(ine),nu_jn(iniHR),nu_ei,aZj(iniHR),Bmag,ion_rays)
   betaj(iniMR) = betaIj(beta_coef(iniMR),rho_j(iniMR),rho_j(ine),nu_jn(iniMR),nu_ei,aZj(iniMR),Bmag,ion_rays)
   betaj(inisT) = betaIj(beta_coef(inisT),rho_j(inisT),rho_j(ine),nu_jn(inisT),nu_ei,aZj(inisT),Bmag,ion_thermal,mass_ionT)
   betaj(inidT) = betaIj(beta_coef(inidT),rho_j(inidT),rho_j(ine),nu_jn(inidT),nu_ei,aZj(inidT),Bmag,ion_thermal,mass_ionT)
   betaj(ing:ing+na-1) = beta_coef(ing:ing+na-1) *aZj(ing:ing+na-1)*Bmag/ nu_jn(ing:ing+na-1)
   do j = 1,nspecies
     if (betaj(j) > 0.0) beta2p11(j) = 1.0/(1.0+betaj(j)**2 )
   end do
   !
   !--Conductivities
   sigmaOHPpa = 0.0
   sigmaHp    = 0.0
   sigmaHn    = 0.0
   do j = 1,nspecies
     sigmaOHPpa(1) = sigmaOHPpa(1) + ns(j)*aZj(j)*betaj(j)
     sigmaOHPpa(3) = sigmaOHPpa(3) + ns(j)*aZj(j)*betaj(j)*beta2p11(j)
     if (Zj(j) > 0.0) then
       sigmaHp     = sigmaHp       + ns(j)* Zj(j)         *beta2p11(j)
     else
       sigmaHn     = sigmaHn       + ns(j)* Zj(j)         *beta2p11(j)
     end if
     do k = j+1,nspecies
       sigmaOHPpa(5) = sigmaOHPpa(5) + ns(j)*aZj(j)*betaj(j)*beta2p11(j) &
                                     * ns(k)*aZj(k)*betaj(k)*beta2p11(k) &
                                     * (sZj(j)*betaj(j) - sZj(k)*betaj(k))**2
     end do
   end do
   if ( abs(sigmaHp+sigmaHn) > epsilon_coef*epsilon(sigmaHp)*max(sigmaHp,-sigmaHn)) then
     sigmaOHPpa(2) = sigmaHp + sigmaHn
   end if
   sigmaOHPpa(1:3) = sigmaOHPpa(1:3)*sigma_coef_onB
   sigmaOHPpa(5)   = sigmaOHPpa(  5)*sigma_coef_onB**2
   if ( abs(sigmaOHPpa(2)) > 0.0 .or. sigmaOHPpa(3) > 0.0) then
     sigmaOHPpa(4)   = 1.0/(sigmaOHPpa(2)**2 + sigmaOHPpa(3)**2)
   end if
 end if
 !
 !--Number densities (for bookkeeping)
 n_densities = 0.0
 if ( get_data_out ) then
   do j = 1,na
     n_densities(1) = n_densities(1) + Zg_on_agi*a_grain(j)*n_grain_coef(j) ! average grain charge
   end do
   n_densities(1) = n_densities(1)/n_total_grain_coef
   n_densities(2) = rho_n                          ! neutral density
   n_densities(3) = rho_j(iniHR)+rho_j(iniMR) &
                  + rho_j(inisT)+rho_j(inidT)      ! total ion density
   n_densities(4) = rho_n*mass_neutral1            ! neutral number density
   n_densities(5) = n_ionR(1)                      ! light element ion number density from cosmic rays
   n_densities(6) = n_ionR(2)                      ! metallic ion number density from cosmic rays
   n_densities(7) = n_ionT(1)                      ! singly ionised ion number density from thermal radiation
   n_densities(8) = n_ionT(2)                      ! doubly ion number density from thermal radiation
   n_densities(9) = n_total_grain_coef*rho_n*mump1 ! grain number density
 end if
 !
end subroutine nicil_ion_get_sigma
!
pure real function betaIj(betaj_coef,rhoj,rhoe,nu_jn,nu_ei,Bmag,aZ,calc_beta,mass)
 ! Function to calculate the Hall parameter for ions
 ! Recall: nu_ie = nu_ei*rho_e/rho_i
 real,    intent(in)           :: betaj_coef,rhoj,rhoe,nu_jn,nu_ei,Bmag,aZ
 real,    intent(in), optional :: mass
 logical, intent(in)           :: calc_beta
 !
 if ( calc_beta .and. rhoj > 0.0 ) then
   if (mod_beta) then
     betaIj = betaj_coef*aZ*Bmag/(nu_jn + nu_ei*rhoe/rhoj)
   else
     betaIj = betaj_coef*aZ*Bmag/ nu_jn
   end if
   if (present(mass)) betaIj = betaIj/mass
 else
   betaIj = 0.0
 end if
 !
end function betaIj
!----------------------------------------------------------------------!
!+
! This is a control routine for NICIL.  It will update
! Zg_on_ag and n_electronT.  These are calculated iteratively, thus
! will be done prior to calculating resistivities.  
!+
!----------------------------------------------------------------------!
pure subroutine nicil_get_ion_nZ(rho,T,Zg_on_ag,n_electronT,ierr)
 implicit none
 real,    intent(in)    :: rho,T
 real,    intent(inout) :: Zg_on_ag,n_electronT
 integer, intent(out)   :: ierr
 logical                :: use_new_guess
 !
 ierr = 0
 !--Exit if using constant resistivities
 if (eta_constant) return
 !
 !--Calcualte the grain charge and electron number density from cosmic rays
 if (ion_rays) then
   use_new_guess = .false.
   if (Zg_on_ag > -epsilon(Zg_on_ag)) Zg_on_ag = -epsilon(Zg_on_ag)                 ! special consideration for Zg_on_ag ~ 0
   call nicil_ionR_get_Zgrain(Zg_on_ag,rho,T,ierr,use_new_guess)
   if (use_new_guess) call nicil_ionR_get_Zgrain(Zg_on_ag,rho,T,ierr,use_new_guess) ! Try with new guess
 end if
 !
 !--Calculate the grain charge from thermal ionisaion and update electron number density
 if (ion_thermal) then
   use_new_guess = .false.
   if (n_electronT < epsilon(n_electronT)) n_electronT = epsilon(n_electronT)      ! special consideration for n_electronT ~ 0
   call nicil_ionT_get_ne(n_electronT,rho,T,ierr,use_new_guess)
   if (use_new_guess) call nicil_ionT_get_ne(n_electronT,rho,T,ierr,use_new_guess) ! Try with new guess
 end if
 !
 return
end subroutine nicil_get_ion_nZ
!======================================================================!
! COSMIC RAY IONISATION-RELATED SUBROUTINES                            !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  This will calculate Zg_on_ag iteratively, using the Newton–Raphson
!  method. 
!  NOTE: By construction, n_grain = (n_e - n_i)/Z_g = m_n/m_g*f_dg*n
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionR_get_Zgrain(Zg_on_ag,rho,T,ierr,use_new_guess)
 implicit none
 integer               :: iter
 integer,intent(inout) :: ierr
 real,   intent(inout) :: Zg_on_ag
 real,   intent(in)    :: rho,T
 logical,intent(inout) :: use_new_guess
 integer               :: j,k
 real                  :: T1,sqrtT1,sqrtT
 real                  :: k_fac,exp_kfac
 real                  :: Za_coef,Za_rat,Za_old,Za_new,fatZa,fatZadZa
 real                  :: sum_k_ig(nimass),sum_k_ig1(nimass),sum_dk_ig(nimass)
 real                  :: sum_k_eg,sum_k_eg1,sum_dk_eg
 real                  :: sum_k1,sum_dk
 logical               :: iterate
 !
 !--Initialise values
 iter    = 0
 Za_rat  = NRtol*2.0
 Za_old  = Zg_on_ag
 T1      = 1.0/T
 sqrtT1  = sqrt(T1)
 sqrtT   = sqrt(T)
 Za_coef = Za_iter_coef/rho
 iterate = .true.
 !
 !--Perform the iterations
 do while ( iterate )
   k_fac     = k_fac_coef*Za_old*T1
   exp_kfac  = exp(k_fac)
   !
   sum_k_ig  = 0.0
   sum_k_eg  = 0.0
   sum_k_ig1 = 0.0
   sum_k_eg1 = 0.0
   sum_dk_ig = 0.0
   sum_dk_eg = 0.0
   do j = 1,na
     sum_k_ig  = sum_k_ig  +  kn_ig_coef(:,j)*sqrtT*(1.0-k_fac)
     sum_k_eg  = sum_k_eg  +  kn_eg_coef(  j)*sqrtT*exp_kfac
     sum_dk_ig = sum_dk_ig + dkn_ig_coef(:,j)*sqrtT1
     sum_dk_eg = sum_dk_eg - dkn_eg_coef(  j)*sqrtT1*exp_kfac
   end do
   do k = 1,nimass
     if (sum_k_ig(k) > 0.0) sum_k_ig1(k) = 1.0/sum_k_ig(k)
   end do
   if (sum_k_eg > 0.0) sum_k_eg1 = 1.0/sum_k_eg
   sum_k1    =  sum_k_eg1
   sum_dk    =  sum_dk_eg*sum_k_eg1**2
   do k = 1,nimass
     sum_k1 = sum_k1 - sum_k_ig1(k)
     sum_dk = sum_dk - sum_dk_ig(k)*sum_k_ig1(k)**2
   end do
   fatZa    = Za_old - Za_coef*sum_k1
   fatZadZa = 1.0    - Za_coef*sum_dk
   !
   Za_new   = Za_old - fatZa/fatZadZa
   Za_rat   = abs( 1.0 - Za_new/Za_old )
   !
   Za_old   = Za_new
   iter     = iter + 1
   if (iter >= NRctrmax .or. Za_rat < NRtol .or. Za_new > 0.0) iterate = .false.
 end do
 !
 ! Set the new value of Z/a
 Zg_on_ag = Za_new
 !
 ! Actions if errors occurred
 if (iter >= NRctrmax .or. Za_new > 0.0) then
   if (use_new_guess) then
     ! New guess failed; trigger warnings     
     if (iter >= NRctrmax) ierr = ierr + 1  ! Zg_on_ag did not converge
     if (Zg_on_ag > 0.0)   ierr = ierr + 2  ! Zg_on_ag > 0
   else
     ! Try again with default guess rather than value from previous iteration
     Zg_on_ag      = -epsilon(Zg_on_ag)
     use_new_guess = .true.
   end if
 end if
 !
end subroutine nicil_ionR_get_Zgrain
!----------------------------------------------------------------------!
!+
! Calculate the electron and ion number densities for
! cosmic ray ionisation
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionR_get_n(n_eleR,n_ionR,T,Zg_on_agi,ierr)
 implicit none
 integer, intent(out) :: ierr
 real,    intent(out) :: n_eleR,n_ionR(nimass)
 real,    intent(in)  :: T,Zg_on_agi
 integer              :: j
 real                 :: sqrtT,k_fac,sum_k_ig(nimass),sum_k_eg
 !
 ierr      = 0
 sqrtT     = sqrt(T)
 k_fac     = k_fac_coef*Zg_on_agi/T
 sum_k_ig  = 0.0
 sum_k_eg  = 0.0
 do j = 1,na
   sum_k_ig  = sum_k_ig  + kn_ig_coef(:,j)*sqrtT*(1.0-k_fac)
   sum_k_eg  = sum_k_eg  + kn_eg_coef(  j)*sqrtT*exp(k_fac)
 end do
 n_ionR = zeta/sum_k_ig
 n_eleR = zeta/sum_k_eg
 !
 !--Verifications against absurb number densities
 if (n_eleR    < 0.0 )                      ierr = ierr + 10 ! Trigger warning if n_electronR < 0
 if (n_ionR(1) < 0.0 .or. n_ionR(2) < 0.0 ) ierr = ierr + 20 ! Trigger warning if n_ionR < 0
 !
 end subroutine nicil_ionR_get_n
!======================================================================!
! THERMAL IONISATION-RELATED SUBROUTINES                               !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  This will solve the Saha equation to determine the electron number
!  density.  The species of interest and thier properties are given
!  at the beginnig of this module.  
!  This assumes for each element, there are only neutral and singly-
!  ionised species.  
!  We will include the electron number density calculated by the cosmic
!  ray ionisaiton as well.
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_ne(n_electronT,rho,T,ierr,use_new_guess)
 implicit none
 integer               :: iter
 integer,intent(inout) :: ierr
 real,   intent(inout) :: n_electronT
 real,   intent(in)    :: rho,T
 logical,intent(inout) :: use_new_guess
 real                  :: mass_ionT_mp,n_ionT(nlevels)
 real                  :: fatn,fatndn,nerat,neold,nenew
 real                  :: neS,dneSdne
 real                  :: Kjk(nlevels,nelements),nj(nelements),njk(nlevels,nelements),ne0
 logical               :: iterate
 !
 !--Initialise values
 iter    = 0
 nerat   = NRtol*2.0
 neold   = n_electronT
 ne0     = neold
 iterate = .true.
 call nicil_ionT_get_nj_Kjk(nj,Kjk,T,rho)
 !
 !--Calcualte n_e
 do while ( iterate )
   call nicil_ionT_get_nion(n_ionT,neS,dneSdne,mass_ionT_mp,neold,njk,nj,Kjk,iter)
   fatn   = neold - neS
   fatndn = 1.0   - dneSdne
   if ( iterate ) then 
     nenew = neold - fatn/fatndn
     if (nenew > 0.0) then
       nerat = abs( 1.0 - neold/nenew )
     else
       nerat = 0.0
     end if
     if (iter >= NRctrmax .or. nerat < NRtol .or. nenew < 0.0) iterate = .false.
   end if
   iter  = iter + 1
   neold = nenew
 end do
 !
 ! Set new electron number density
 n_electronT = nenew
 !
 ! Actions if errors occurred
 if (iter >= NRctrmax .or. nenew < 0.0) then
   if (use_new_guess) then
     ! New guess failed; trigger warnings     
     ierr = ierr + 4  ! n_electronT did not converge
     ierr = ierr + 8  ! n_electronT < 0
   else
     ! Try again with default guess rather than value from previous iteration
     n_electronT   = epsilon(n_electronT)
     use_new_guess = .true.
   end if
 end if
 !
end subroutine nicil_ionT_get_ne
!----------------------------------------------------------------------!
!+
!  This is a stand-alone routine to calculate n_ionT, mass_ionT & njk
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_n(n_ionT,mass_ionT,njk,n_electronT,T,rho)
 implicit none
 real,    intent(out) :: mass_ionT,n_ionT(:),njk(:,:)
 real,    intent(in)  :: n_electronT,T,rho
 integer              :: iter
 real                 :: Kjk(nlevels,nelements),nj(nelements),neS,dneSdne,mass_ionT_mp
 !
 iter = 0
 call nicil_ionT_get_nj_Kjk(nj,Kjk,T,rho)
 call nicil_ionT_get_nion(n_ionT,neS,dneSdne,mass_ionT_mp,n_electronT,njk,nj,Kjk,iter)
 mass_ionT = mass_ionT_mp*mass_proton
 !
end subroutine nicil_ionT_get_n
!----------------------------------------------------------------------!
!+
!  This is a stand-alone routine to calculate nj & Kjk
!  (done here due to restriction on Kjk)
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_nj_Kjk(nj,Kjk,T,rho)
 implicit none 
 real,    intent(out) :: nj(:),Kjk(:,:)
 real,    intent(in)  :: T,rho
 integer              :: j,k,isize
 !
 isize = size(nj)
 nj  = abundancej(1:isize)*rho*mump1
 Kjk = gej(:,1:isize)*Saha_coef*sqrt(T)**3*exp(-chij(:,1:isize)/T)
 !
 do k = 1,nelements
   do j = 1,nlevels
     if (T < Texp_thresh(j,k) ) Kjk(j,k) = 0.0
   end do
 end do
 !
end subroutine nicil_ionT_get_nj_Kjk
!----------------------------------------------------------------------!
!+
!  This calculates the total ion number density, and average ion mass
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_nion(ni,neS,dneSdne,mi_mp,ne,njk,nj,Kjk,iter)
 implicit none
 integer, intent(in)  :: iter
 real,    intent(out) :: ni(:),neS,dneSdne,mi_mp,njk(:,:)
 real,    intent(in)  :: ne,nj(:),Kjk(:,:)
 integer              :: k
 real                 :: term,term1,ne1,KKonne,mD
 !
 njk     = 0.0
 ni      = 0.0
 neS     = 0.0
 dneSdne = 0.0
 mi_mp   = 0.0
 mD      = 0.0
 if ( ne > 0.0 .and. (iter==0 .or. (iter>0 .and. abs(ne-epsilon(ne)) > epsilon(ne)**2 ) ) ) then 
   ne1 = 1.0/ne
   do k = 1,nelements
     KKonne = Kjk(1,k)*Kjk(2,k)*ne1
     term = ne + Kjk(1,k) + KKonne
     if (term > 0.0) then
       term1 = 1.0/term
     else
       term1 = 0.0
     end if
     if (Kjk(1,k) > 0.0) njk(1,k) = nj(k)*Kjk(1,k)*term1
     if (Kjk(2,k) > 0.0) njk(2,k) = njk(1,k)*Kjk(2,k)*ne1
     ni      = ni      + njk(:,k)
     neS     = neS     + njk(1,k) + 2.0*njk(2,k)
     dneSdne = dneSdne - njk(1,k)*term1*(1.0-KKonne*ne1)   &
                       - 2.0*njk(2,k)*term1*ne1*(2.0*ne+Kjk(1,k))
     mD      = mD      + (njk(1,k)+njk(2,k))*sqrtmj1(k)
   end do
   if (mD > 0.0 .and. ni(1)+ni(2) > 0.0) mi_mp  = ((ni(1)+ni(2))/mD)**2
 end if
 !
end subroutine nicil_ionT_get_nion
!======================================================================!
! NON-IDEAL MHD-RELATED SUBROUTINES                                    !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  Calculates the coefficients for the non-ideal MHD terms:
!    Ohmic Resistivity
!    Hall effect 
!    Ambipolar diffusion
!  Note: eta_ambi = csqbyfourpi * sigmaP * sigmaPerp21 - eta_ohm 
!                 = eta_ohm*(sigmaO*sigmaP - sigmaPerp2)*sigmaPerp21
!                 = eta_ohm*sigma2A*sigmaPerp21
!+
!----------------------------------------------------------------------!
pure subroutine nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,T,sigmaOHPpa)
 implicit none
 real, intent(out) :: eta_ohm,eta_hall,eta_ambi
 real, intent(in)  :: T,sigmaOHPpa(5)
 !
 eta_ohm  = 0.0
 eta_hall = 0.0
 eta_ambi = 0.0
 if (T < T_thresh .and. sigmaOHPpa(1) > 0.0) then
                 eta_ohm  = csqbyfourpi / sigmaOHPpa(1)
   if (use_hall) eta_hall = csqbyfourpi * sigmaOHPpa(2) * sigmaOHPpa(4)
   if (use_ambi) eta_ambi = eta_ohm     * sigmaOHPpa(5) * sigmaOHPpa(4)
   if (.not. use_ohm ) eta_ohm = 0.0
 end if
 !
 return
end subroutine nicil_nimhd_get_eta
!-----------------------------------------------------------------------
pure subroutine nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
 implicit none
 real, intent(out) :: eta_ohm,eta_hall,eta_ambi
 real, intent(in)  :: Bfield,rho
 !
 eta_ohm  = 0.0
 eta_hall = 0.0
 eta_ambi = 0.0
 if ( use_ohm  ) eta_ohm  = eta_ohm_cnst
 if ( use_hall ) eta_hall = eta_hall_cnst*Bfield
 if ( use_ambi ) eta_ambi = eta_ambi_cnst*Bfield**2/rho
 !
 return
end subroutine nicil_nimhd_get_eta_cnst
!-----------------------------------------------------------------------
!+
!  Calculates JxB and (JxB)xB
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_jcbcb(jcbcb,jcb,jcurrent,Bx,By,Bz,B1)
 implicit none
 real, intent(out) :: jcb(3), jcbcb(3)
 real, intent(in)  :: jcurrent(3)
 real, intent(in)  :: Bx, By, Bz, B1
 !
 jcb(1)   = ( jcurrent(2)*Bz - jcurrent(3)*By )*B1
 jcb(2)   = ( jcurrent(3)*Bx - jcurrent(1)*Bz )*B1
 jcb(3)   = ( jcurrent(1)*By - jcurrent(2)*Bx )*B1
 !
 jcbcb(1) = ( jcb(2)*Bz      - jcb(3)*By      )*B1
 jcbcb(2) = ( jcb(3)*Bx      - jcb(1)*Bz      )*B1
 jcbcb(3) = ( jcb(1)*By      - jcb(2)*Bx      )*B1
 !
 return
end subroutine nimhd_get_jcbcb
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD contributions to the magnetic field in SPH
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dBdt(dBnonideal,etaohm,etahall,etaambi &
                              ,jcurrent,jcb,jcbcb,dxr1,dyr1,dzr1)
 implicit none
 real, intent(out) :: dBnonideal(3)
 real, intent(in)  :: jcurrent(3),jcb(3),jcbcb(3)
 real, intent(in)  :: etaohm,etahall,etaambi,dxr1,dyr1,dzr1
 real              :: dBohm(3),dBhall(3),dBambi(3)
 !
 dBohm  = 0.0
 dBhall = 0.0
 dBambi = 0.0
 !
 if (use_ohm ) call nimhd_get_DcrossR(dBohm ,jcurrent,dxr1,dyr1,dzr1,etaohm )
 if (use_hall) call nimhd_get_DcrossR(dBhall,jcb     ,dxr1,dyr1,dzr1,etahall)
 if (use_ambi) call nimhd_get_DcrossR(dBambi,jcbcb   ,dxr1,dyr1,dzr1,etaambi)
 dBnonideal = dBambi - dBhall - dBohm 
 !
 return
end subroutine nimhd_get_dBdt
!-----------------------------------------------------------------------
!+
!  performs simple cross product 
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_DcrossR(DcrossR,D_in,dx,dy,dz,eta)
 implicit none
 real, intent (in)  :: D_in(3)
 real, intent (out) :: DcrossR(3)
 real, intent (in)  :: dx, dy, dz, eta
 !
 DcrossR(1) = (D_in(2)*dz - D_in(3)*dy)*eta
 DcrossR(2) = (D_in(3)*dx - D_in(1)*dz)*eta
 DcrossR(3) = (D_in(1)*dy - D_in(2)*dx)*eta
 !
 return
end subroutine nimhd_get_DcrossR
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD contributions to energy
!  Note: dudthall==0
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dudt(dudtnonideal,etaohm,etaambi,rho,J,B)
 implicit none
 real, intent(out) :: dudtnonideal
 real, intent(in)  :: J(3),B(3)
 real, intent(in)  :: etaohm,etaambi,rho
 !
 dudtnonideal = 0.0
 dudtnonideal = dudtnonideal + etaohm /rho*dot_product(J,J)
 dudtnonideal = dudtnonideal + etaambi/rho &
                             * (dot_product(J,J)*dot_product(B,B) &
                               -dot_product(B,J)*dot_product(B,J))
 !
 return
end subroutine nimhd_get_dudt
!-----------------------------------------------------------------------
!+
!  Calculates the timesteps for non-ideal MHD
!+
!-----------------------------------------------------------------------
subroutine nimhd_get_dt(dtohm,dthall,dtambi,h,etaohm,etahall,etaambi)
 implicit none
 real, intent (in)  :: h,etaohm,etahall,etaambi
 real, intent (out) :: dtohm,dthall,dtambi
 real               :: h2
 !
 dtohm  = huge(dtohm)
 dthall = huge(dthall)
 dtambi = huge(dtambi)
 h2     = h*h
 !
 if (use_ohm  .and.     etaohm   > tiny(etaohm ) ) dtohm  =      C_nimhd*h2/etaohm
 if (use_hall .and. abs(etahall) > tiny(etahall) ) dthall = abs( C_nimhd*h2/etahall )
 if (use_ambi .and.     etaambi  > tiny(etaambi) ) dtambi =      C_nimhd*h2/etaambi       
 !
 return
end subroutine nimhd_get_dt
!----------------------------------------------------------------------!
end module nicil
