!----------------------------------------------------------------------!
! This is a stand-alone library that will calculate ionisation values
! and the coefficients for the non-ideal MHD terms (Ohmic resistivity,
! Hall Effect and Ambipolar diffusion)
!----------------------------------------------------------------------!
!+
!  MODULE: nicil
!          
!  DESCRIPTION:
!  Contains routines to calculation the ionisation rate, and related
!  useful quantities.
!
!  REFERENCES: WN99 = Wardle & Ng (1999)
!              PG08 = Pinto & Galli (2008)
!              PW08 = Pandey & Wardle (2008)
!              KW14 = Keith & Wardle (2014)
!              ****THERE ARE MORE!
!
!  AUTHOR: James Wurster
!
!  RUNTIME PARAMETERS:
!    a_grain     -- grain radius (cm)
!    fdg         -- gas to dust mass ratio
!    mass_ion_mp -- ion mass (m_proton) 
!    metal_X     -- metalicity: Hydrogen fraction
!    metal_Y     -- metalicity: Helium fraction
!    rho_bulk    -- bulk grain density (g/cm^3) 
!    zeta        -- cosmic ray ionisation rate (s^-1) 
!    ****THERE ARE MORE!
!----------------------------------------------------------------------!
module nicil
 implicit none
 !--INPUT PAREMETERS
 !--Turn on/off individual non-ideal MHD coefficients
 logical, private :: use_ohm           = .false.           ! Calcualate coefficient for Ohmic resistivity
 logical, private :: use_hall          = .false.           ! Calcualate coefficient for the Hall effect
 logical, private :: use_ambi          = .true.           ! Calcualate coefficient for ambipolar diffusion
 !--Set constant grain size (g_cnst=true) or approximate an MRN grain distribution (g_cnst=false)
 logical, private :: g_cnst            = .true.           ! Calcualate coefficient for ambipolar diffusion
 !--Set ionisation type (can have one or both = true)
 logical, private :: ion_rays          = .true.           ! Include ionisation from cosmic (or x-) rays (=true)
 logical, private :: ion_thermal       = .true.           ! Include thermal ionisation (=true)
 !--Assume that all the Hydrogen is molecular Hydrogen (used only for calculating mean molecular mass)
 logical, private :: H_is_H2           = .true.
 !--Use constant resistivity coefficients for all three resistivity terms
 logical, private :: eta_constant      = .true. 
 !
 real,    private :: T_thresh          = 1.0e3            ! Temperature [K] above which non-ideal MHD will turn off
 !
 real,    private :: fdg               = 0.01             ! Grain Parameter: gas to dust mass ratio 
 real,    private :: a_grain           = 1.0d-5           ! Grain Parameter: grain radius for constant grain size [cm]
 real,    private :: a0_grain          = 5.0d-7           ! Grain Parameter: minimum grain radius for ~MRN distribution [cm]
 real,    private :: a1_grain          = 2.5d-5           ! Grain Parameter: maximum grain radius for ~MRN distribution [cm]
 real,    private :: rho_bulk          = 3.0              ! Grain Parameter: bulk grain density [g/cm^3]
 real,    private :: zeta_cgs          = 1.0d-17          ! ionisation rate [s^-1]
 real,    private :: mass_ion_mp       = 24.3             ! mass of ion (default is mass of magnesium) [m_proton]
 ! 
 real,    private :: delta_gn          = 1.3              ! Ionisation parameter: multiplicative factor for sigmavgnbyT 
 real,    private :: sigmav_eH2        = 3.16d-11         ! Ionisation parameter: momentum transfer coefficient for e-H2 [cm^3/s]
 real,    private :: sigmav_eHe        = 7.08d-11         ! Ionisation parameter: momentum transfer coefficient for e-He [cm^3/s]
 real,    private :: pnH2              = 0.804            ! Ionisation parameter: polarizability for H2 [angstroms^3] 
 real,    private :: pnHe              = 0.207            ! Ionisation parameter: polarizability for He [angstroms^3]
 !real,   private :: Z_electron        = -1.0             ! Ionisation parameter: Charge on electron (definition)
 !real,   private :: Z_ion             =  1.0             ! Ionisation parameter: Charge on ion (assumption)
 !
 !--Thermal ionisation
 integer, parameter :: nelements       =     5            ! Number of elements
 integer, parameter :: n_nuin          =  1000            ! Number of ni_inT_coef values to pre-caculate
 real,    parameter :: se              =    1.0           ! Electron ticking coefficient: se \in (1.0d-3, 1.0)
 real,    parameter :: m_min           =    1.0           ! minimum ion mass (units of m_proton) in the table
 real,    parameter :: m_max           =  101.0           ! maximum ion mass (units of m_proton) in the table
 !--The number of species
 integer, parameter :: nspecies        =    5             ! The number of charged species (do not modify)
 !--Resistivity coefficients (if fixed as constants)
 real,    private   :: n_e_cnst        = 1.0e7            ! Constant electron number density
 real,    private   :: rho_i_cnst      = 0.1              ! Density of ionised gas
 real,    private   :: gamma_AD        = 1000.0           ! Collisional coupling coefficient
 logical, private   :: hall_lt_zero    = .false.          ! The sign of the Hall coefficient
 !
 !--Physical Constants (CGS)
 real, parameter  :: pi                =  3.1415926536d0  !  pi
 real, parameter  :: twopi             =  6.2831853072d0  ! 2pi
 real, parameter  :: fourpi            = 12.5663706144d0  ! 4pi
 real, parameter  :: c                 =  2.997924d10     ! Speed of light [cm/s]
 real, parameter  :: qe                =  4.8032068d-10   ! charge on electron [esu == statC]
 real, parameter  :: mass_electron_cgs =  9.10938291d-28  ! Electron mass [g]
 real, parameter  :: eV                =  1.60217657d-12  ! Electron volts [ergs]
 real, parameter  :: mass_proton_cgs   =  1.67262158d-24  ! Proton mass [g]
 real, parameter  :: kboltz            =  1.38066d-16     ! Boltzmann constant  [erg/K]
 real, parameter  :: planckh           =  6.6260755d-27   ! Planck's Constant [erg s]
 real, parameter  :: Amrn              =  1.5d-25         ! Constant for MRN grain distribution [cm^2.5]
 !
 !--Local Variables to the NICIL module
 integer, private :: iprint
 real,    private :: C_nimhd, cs_thresh,csqbyfourpi
 real,    private :: mass_electron,mass_proton,mass_ion,mump1,mass_neutral_cgs,mass_neutral1
 real,    private :: n_ele_coef,n_ion_coef,n_grain_coef
 real,    private :: Zgrain0,Z_iter_coef,nelectron0
 real,    private :: k_ig_coef,k_eg_coef,k_ig_coef1,k_eg_coef1,k_fac_coef
 real,    private :: beta_e_coef,beta_iR_coef,beta_iT_coef,beta_g_coef
 real,    private :: nu_ei_coef,nu_en_coef,nu_inR_coef,nu_gn_coef
 real,    private :: sigma_coef,sigmavenXbyT,sigmavenYbyT,sigmavenX_LR,sigmavenY_LR
 real,    private :: eta_ohm_cnst,eta_hall_cnst,eta_ambi_cnst
 real,    private :: Saha_coef,y_ionT_coef,tau0_coef,tau_coef
 real,    private :: log_abunj(nelements),mj(nelements),abundancej(nelements),mass_frac(nelements)
 real,    private :: depletionj(nelements),chij(nelements),gej(nelements)
 real,    private :: dm_ion,nu_inT_coef(n_nuin)
 real,    private :: metal_X,metal_Y,umass0,unit_density0
 character(2), private :: symj(nelements)
 !
 !--Subroutines
 public  :: nicil_initial,nicil_get_eta
 public  :: nimhd_get_jcbcb,nimhd_get_sphdBdt,nimhd_get_dudt,nimhd_get_dt
 public  :: nicil_ionR_get_Zgrain
 private :: nicil_ion_get_sigma
 private :: nicil_ionT_get_ne,nicil_ionT_get_ng
 private :: nicil_nimhd_get_eta,nimhd_get_DcrossR
 private :: nicil_fatal, nicil_test_eta
 !
 private
!
contains
!
!----------------------------------------------------------------------!
! Define the properties of the species that will be used for thermal 
! ionisation
subroutine nicil_initial_species
 !  
 gej           =  2.0             ! Statistical weight of the electron
 !--Hydrogen
 symj(1)       = "H"              ! Chemical Symbol
 log_abunj(1)  = 12.0             ! Log abundance
 mj(1)         =  1.01            ! Mass [m_proton]
 chij(1)       = 13.60            ! Ionisation potential eV]
 depletionj(1) =  0.0             ! Depletion factor
 !--Helium 
 symj(2)       = "He"             ! Chemical Symbol
 log_abunj(2)  = 10.93            ! Log abundance
 mj(2)         =  4.00            ! Mass [m_proton]
 chij(2)       = 24.59            ! Ionisation potential [eV]
 depletionj(2) =  0.0             ! Depletion factor
 !--Sodium 
 symj(3)       = "Na"             ! Chemical Symbol
 log_abunj(3)  =  6.24            ! Log abundance
 mj(3)         = 22.98            ! Mass [m_proton]
 chij(3)       =  5.14            ! Ionisation potential [eV]
 depletionj(3) = -0.92            ! Depletion factor
 !--Magnesium 
 symj(4)       = "Mg"             ! Chemical Symbol
 log_abunj(4)  =  7.60            ! Log abundance
 mj(4)         = 24.31            ! Mass [m_proton]
 chij(4)       =  7.65            ! Ionisation potential [eV]
 depletionj(4) = -0.92            ! Depletion factor
 !--Potassium 
 symj(5)       = "K"              ! Chemical Symbol
 log_abunj(5)  =  5.03            ! Log abundance
 mj(5)         = 39.10            ! Mass [m_proton]
 chij(5)       =  4.34            ! Ionisation potential [eV]
 depletionj(5) = -0.92            ! Depletion factor
 !
end subroutine nicil_initial_species
!======================================================================!
! INITIALISATION & CONTROL ROUTINES                                    !
!======================================================================!
!----------------------------------------------------------------------!
!
! Initialisation subroutine
! This will initialise all the variable required for NICIL, including 
! the frequently used coefficients.  
! All coefficients will be converted to code units.
!
!----------------------------------------------------------------------!
subroutine nicil_initial(utime,udist,umass,unit_Bfield,Z_grainR,n_electron,iprint_in,test_type)
 implicit none
 real,              intent(in)    :: utime,umass,udist,unit_Bfield
 real,    optional, intent(inout) :: Z_grainR(:),n_electron(:)
 integer, optional, intent(in)    :: iprint_in
 integer, optional, intent(out)   :: test_type
 integer                          :: j
 real                             :: unit_velocity,unit_density,unit_erg,unit_eta
 real                             :: sigmavin,sigmavgnbyT
 real                             :: meanmolmass,mass_proton_mp,mass_electron_mp,mump
 real                             :: mass_grain, mass_neutral, mass_grain_cgs, mass_ion_cgs
 real                             :: a01_grain
 real                             :: zeta
 real                             :: vrms_cgs,vrms_mks, mu_iX,mu_iY,mu_eX,mu_eY
 real                             :: mass_total,n_rel,nj_rel(nelements)
 logical                          :: test_eta
 !
 !--Initialise species properties for thermal ionisation & Calculate abundances & mass fractions
 !  By construction, Sum (abundancej), Sum(mass_frac) == 1 
 call nicil_initial_species
 n_rel       = 0.0
 mass_total  = 0.0
 meanmolmass = 0.0
 do j = 1,nelements
   nj_rel(j) = 10**(log_abunj(j) - 12.0 + depletionj(j))  ! Number density relative to H
   n_rel     = n_rel + nj_rel(j)
 end do
 do j = 1,nelements
   abundancej(j) = nj_rel(j)/n_rel                     ! Abundance
   mass_total    = mass_total + mj(j)*abundancej(j)
 end do
 do j = 1,nelements
   mass_frac(j)  = mj(j)*abundancej(j)/mass_total
   if (H_is_H2 .and.  trim(symj(j))=="H") then 
     meanmolmass = meanmolmass + mass_frac(j)/(2.0*mj(j))
   else
     meanmolmass = meanmolmass + mass_frac(j)/mj(j)
   end if
 end do
 meanmolmass = 1.0/meanmolmass
 !
 !--Additional unit conversions from cgs-> code
 unit_velocity     = udist/utime
 unit_density      = umass/udist**3
 unit_erg          = umass*unit_velocity**2
 unit_eta          = udist**2/utime 
 umass0            = umass
 unit_density0     = unit_density
 !
 !--Initialise variables (CGS variables)
 C_nimhd           = 1.0/twopi                                               ! 'Courant' factor for the non-ideal MHD timestep
 mass_ion_cgs      = mass_ion_mp*mass_proton_cgs                             ! Ion mass 
 mass_neutral_cgs  = meanmolmass*mass_proton_cgs                             ! Neutral mass 
 mump              = meanmolmass*mass_proton_cgs                             ! mean particle mass assuming composition is H2 & He 
 if (g_cnst) then 
   mass_grain_cgs  = fourpi/3.0*a_grain**3*rho_bulk                          ! Grain particle mass (constant)
   a01_grain       = a_grain
 else
   mass_grain_cgs  = 5.0*fourpi/3.0*rho_bulk*(a0_grain*a1_grain)**2.5 &      ! Average grain particle mass (approximately MRN)
                   * (a1_grain**0.5-a0_grain**0.5)/(a1_grain**2.5-a0_grain**2.5)
   a01_grain       =  5.0/3.0*a0_grain*a1_grain                       &      ! Average grain radius (approximately MRN)
                   * (a1_grain**1.5-a0_grain**1.5)/(a1_grain**2.5-a0_grain**2.5)
 end if
 if (trim(symj(1))=="H") then
   metal_X         = mass_frac(1)                                          ! Metallicity of Hydrogen
 else
   call nicil_fatal('nicil_test_eta','Hygrogen does not have the correct array number of 1')
 end if
 if (trim(symj(2))=="He") then
   metal_Y         = mass_frac(2)                                          ! Metallicity of Helium
 else
   call nicil_fatal('nicil_test_eta','Helium does not have the correct array number of 1')
 end if
 !
 !--Initialise variables (Code unit variables)
 csqbyfourpi       = c**2/fourpi              / unit_eta                     ! coefficient for sigma
 cs_thresh         = sqrt(kboltz*T_thresh/mass_neutral_cgs) /unit_velocity   ! Temperature threshhold in soundspeed 
 zeta              = zeta_cgs * utime                                        ! Ionisation rate
 mass_proton       = mass_proton_cgs          / umass
 mass_neutral      = mass_neutral_cgs         / umass
 mass_electron     = mass_electron_cgs        / umass
 mass_ion          = mass_ion_cgs             / umass
 mass_grain        = mass_grain_cgs           / umass
 !
 !--Initialise variables (m_proton units)
 mass_proton_mp    = 1.0
 mass_electron_mp  = mass_electron_cgs/mass_proton_cgs
 !
 !--Charge capture Rates coefficients (Code units)
 if (g_cnst) then 
   k_ig_coef   = a_grain**2*sqrt(8.0*pi*mump/mass_ion_cgs     ) /udist**2         
   k_eg_coef   = a_grain**2*sqrt(8.0*pi*mump/mass_electron_cgs) /udist**2         
   k_fac_coef  = qe**2/(a_grain*mump)                           /unit_velocity**2 
 else
   k_ig_coef   = 5.0*sqrt(8.0*pi*mump/mass_ion_cgs     )*(a0_grain*a1_grain)**2 &
               * (a1_grain**0.5-a0_grain**0.5)/(a1_grain**2.5-a0_grain**2.5)    /udist**2         
   k_eg_coef   = 5.0*sqrt(8.0*pi*mump/mass_electron_cgs)*(a0_grain*a1_grain)**2 &
               * (a1_grain**0.5-a0_grain**0.5)/(a1_grain**2.5-a0_grain**2.5)    /udist**2    
   k_fac_coef  = 5.0*qe**2/(7.0*mump*a0_grain*a1_grain)                         &
               * (a1_grain**3.5-a0_grain**3.5)/(a1_grain**2.5-a0_grain**2.5)    /unit_velocity**2 
 end if
 k_ig_coef1  = 1.0/k_ig_coef
 k_eg_coef1  = 1.0/k_eg_coef
 !
 !--Rate coefficientes (CGS; sigma units are cm^3/s)
 vrms_cgs      = sqrt(8.0*mump/(pi*mass_electron_cgs))                            *unit_velocity 
 vrms_mks      = vrms_cgs                                                         *1.0e-5
 mu_iX         = mass_proton_cgs*(2.0*mj(1) * mass_ion_mp)     /(2.0*mj(1) + mass_ion_mp)
 mu_iY         = mass_proton_cgs*(    mj(2) * mass_ion_mp)     /(    mj(2) + mass_ion_mp)
 mu_eX         = mass_proton_cgs*(2.0*mj(1) * mass_electron_mp)/(2.0*mj(1) + mass_electron_mp)
 mu_eY         = mass_proton_cgs*(    mj(2) * mass_electron_mp)/(    mj(2) + mass_electron_mp)
 sigmavin      = 2.81d-9*(metal_X*sqrt(pnH2 * mass_proton_cgs/mu_iX) + metal_Y*sqrt(pnHe * mass_proton_cgs/mu_iY) )
 sigmavenXbyT  = metal_X*sigmav_eH2*vrms_mks**1.3                                               
 sigmavenYbyT  = metal_Y*sigmav_eHe*vrms_mks                                                     
 sigmavenX_LR  = metal_X*2.81d-9*sqrt(pnH2 * mass_proton_cgs/mu_eX)                
 sigmavenY_LR  = metal_Y*2.81d-9*sqrt(pnHe * mass_proton_cgs/mu_eY)                              
 if (g_cnst) then
   sigmavgnbyT = a_grain**2*delta_gn * sqrt(128.0*pi*mump/(9.0*mass_neutral_cgs)) *unit_velocity 
 else
   sigmavgnbyT = 5.0*delta_gn * sqrt(128.0*pi*mump/(9.0*mass_neutral_cgs))*(a0_grain*a1_grain)**2 &
               * (a1_grain**0.5-a0_grain**0.5)/(a1_grain**2.5-a0_grain**2.5)       *unit_velocity 
 end if
 !
 !--Collisional Frequencies (CGS = s^-1)
 nu_ei_coef  = 51.0 * (kboltz/mump)**1.5                        /(udist**3*unit_velocity**3)     
 nu_inR_coef = sigmavin   /(mass_neutral_cgs+mass_ion_cgs     ) * unit_density                  
 nu_en_coef  = 1.0        /(mass_neutral_cgs+mass_electron_cgs) * unit_density                   
 nu_gn_coef  = sigmavgnbyT/(mass_neutral_cgs+mass_grain_cgs   ) * unit_density
 if (ion_thermal) then
   call nicil_nuinTcoef
 else
   nu_inT_coef = 0.0
 end if
 !
 !--Number density coefficients (Code units)
 if (g_cnst) then
   n_grain_coef = mass_neutral_cgs*fdg/mass_grain_cgs
 else
   n_grain_coef = 0.4*Amrn*(1.0/a0_grain**2.5 - 1.0/a1_grain**2.5)
 end if
 n_ele_coef   = zeta/n_grain_coef*k_eg_coef1
 n_ion_coef   = zeta/n_grain_coef*k_ig_coef1
 !
 !--Hall parameters (dimensionless after multipiled by B_code (and m_ionT_code))
 beta_e_coef  = qe/(mass_electron_cgs*c) * unit_Bfield                             
 beta_iR_coef = qe/(mass_ion_cgs     *c) * unit_Bfield                             
 beta_iT_coef = qe/(                  c) * unit_Bfield    /umass
 beta_g_coef  = qe/(mass_grain_cgs   *c) * unit_Bfield                            
 !
 !--Conductivity coefficient (code units)
 sigma_coef   = qe*c                     / (udist**3 * unit_Bfield)                
 !
 !--Particle Masses (Code units)
 mass_neutral1= 1.0/mass_neutral
 mump1        = 1.0/mump                 * umass 
 !
 !--Grain charge coefficient and initial guess (Code units)
 Z_iter_coef  = zeta/n_grain_coef**2
 Zgrain0      = -epsilon(Zgrain0)
 if (present(Z_grainR)) Z_grainR = Zgrain0
 !--Initial guess of n_electron (required only for thermal ionisation)
 nelectron0   = epsilon(nelectron0)
 if (present(n_electron)) n_electron = nelectron0
 !
 !--Coefficient for the Saha equation (Code units)
 Saha_coef = (2.0*pi*(mump/umass)*mass_electron/(planckh/(unit_erg*utime))**3)**1.5
 chij      = chij*eV/unit_erg*mump1
 !--Coefficients to calculate grain properties from thermal ionisation
 !  Tau is dimensionless after being multiples by (mu c_s^2)_code
 y_ionT_coef = exp(1.0)*sqrt(mass_proton/mass_electron)
 tau0_coef   = 8.0*mass_electron/(pi*mass_proton)
 if (g_cnst) then
   tau_coef  = a_grain*mass_proton_cgs/qe**2                      * unit_velocity**2
 else
   tau_coef  = 5.0*mass_proton_cgs/(3.0*qe**2)*a0_grain*a1_grain  * unit_velocity**2 &
             * (a1_grain**1.5-a0_grain**1.5)/(a1_grain**2.5-a0_grain**2.5)
 end if
 !
 !--Set constant coefficients
 eta_ohm_cnst  = mass_electron_cgs*c**2/(fourpi*qe**2*n_e_cnst)                                 / unit_eta
 eta_hall_cnst = c/(fourpi*qe*n_e_cnst)                         *  unit_Bfield                  / unit_eta
 eta_ambi_cnst = 1.0/(fourpi*gamma_AD*rho_i_cnst)              * (unit_Bfield**2/unit_density) / unit_eta
 if (hall_lt_zero) eta_hall_cnst = -eta_hall_cnst
 !
 !--Set the location of the printouts
 if (present(iprint_in)) then
   iprint = iprint_in
 else
   iprint = 6
 end if
 !
 !--Determine what we are testing, if requested
 test_eta = .false.
 if (present(test_type)) then
   test_type = 0
   if (use_ohm .and. use_hall .and. use_ambi) then
     if      (g_cnst      .and.      ion_rays .and. .not.ion_thermal) then 
       test_type = 1
       test_eta  = .true.
     else if (.not.g_cnst .and.      ion_rays .and. .not.ion_thermal) then 
       test_type = 2
     else if (g_cnst      .and. .not.ion_rays .and.      ion_thermal) then 
       test_type = 3
     else if (.not.g_cnst .and. .not.ion_rays .and.      ion_thermal) then 
       test_type = 4
     else if (g_cnst      .and.      ion_rays .and.      ion_thermal) then 
       test_type = 5
     else if (.not.g_cnst .and.      ion_rays .and.      ion_thermal) then 
       test_type = 6
     end if       
   end if
 end if
 !
 !--Print Statements to summarise conditions used
 call nicil_print_summary(a01_grain,mass_grain_cgs,meanmolmass)
 !
 !--Test calculations of eta
 test_eta = .false.
 call nicil_test_eta(unit_density,unit_Bfield,unit_eta,sqrt(kboltz*30.0/mass_neutral_cgs) /unit_velocity,test_eta)
 return
end subroutine nicil_initial
!----------------------------------------------------------------------!
!
! Routines to calculate nu_inT_coef tables, and extract the relevant 
! value
!
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
subroutine nicil_calc_nuinTcoef(nu_inT_coef0,mion)
 real,    intent(in)  :: mion
 real,    intent(out) :: nu_inT_coef0
 real                 :: mion_mp,mu_iX,mu_iY,sigmavin
 mion_mp      = mion/mass_proton
 mu_iX        = mass_proton_cgs*(2.0*mj(1) * mion_mp)/(2.0*mj(1) + mion_mp)
 mu_iY        = mass_proton_cgs*(    mj(2) * mion_mp)/(    mj(2) + mion_mp)
 sigmavin     = 2.81d-9*(metal_X*sqrt(pnH2 * mass_proton_cgs/mu_iX) + metal_Y*sqrt(pnHe * mass_proton_cgs/mu_iY) )
 nu_inT_coef0 = sigmavin   /(mass_neutral_cgs+mion*umass0) * unit_density0
end subroutine nicil_calc_nuinTcoef
!
!--Subroutine to extract a coefficient from the table
subroutine nicil_extract_nuinTcoef(nu_inT_coef0,mion)
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
!
! Print Statements to summarise what is being used here
!
!----------------------------------------------------------------------!
subroutine nicil_print_summary(a01_grain,mass_grain_cgs,meanmolmass)
 real, intent(in)                 :: a01_grain,mass_grain_cgs,meanmolmass
 integer                          :: j
 character(2)                     :: comma
 character(200)                   :: ni_terms
 !
 ni_terms = ""
 comma    = ""
 if (use_ohm ) then
   write(ni_terms,'(2a)') trim(ni_terms),"Ohmic resistivity"
   comma = ", "
 end if
 if (use_hall) then
   write(ni_terms,'(3a)') trim(ni_terms),comma,"Hall effect"
   comma = ", "
 end if
 if (use_ambi) write(ni_terms,'(3a)') trim(ni_terms),comma,"Ambipolar Diffusion"
 write(iprint,'(a)')             "NICIL: Initialisation complete."
 if (eta_constant) then
   write(iprint,'(2a)')            "NICIL: Non-ideal terms used: ",trim(ni_terms)
   write(iprint,'(a)' )           "NICIL: All resistivity coefficients are constant."
   if (use_ohm ) write(iprint,'(a,Es10.3)')           "NICIL: Eta_ohm          = ", eta_ohm_cnst
   if (use_hall) write(iprint,'(a,Es10.3)')           "NICIL: Eta_Hall/B       = ", eta_hall_cnst
   if (use_ambi) write(iprint,'(a,Es10.3)')           "NICIL: Eta_ambi*rho/B^2 = ", eta_ambi_cnst
 else
   if (g_cnst) then 
     write(iprint,'(a)')           "NICIL: Using constant grain size."
   else
     write(iprint,'(a)')           "NICIL: Approximating grain distribution with MRN distribution."
   end if
   write(iprint,'(a)')             "NICIL: Species  abundance   mass fraction"
   do j = 1,nelements
     if (mass_frac(j) > 0.01) then 
       write(iprint,'(3a,Es10.2,a,F5.3)')   "NICIL: ",symj(j),"     ",abundancej(j),"     ",mass_frac(j)
     else
       write(iprint,'(3a,Es10.2,a,Es10.2)') "NICIL: ",symj(j),"     ",abundancej(j),"   ",mass_frac(j)
     end if
   end do
   if (H_is_H2) then 
     write(iprint,'(a,Es10.3,a)')  "NICIL: Mean molecular mass: ",meanmolmass," m_proton  (assuming molecular hydrogen)"
   else
     write(iprint,'(a,Es10.3,a)')  "NICIL: Mean molecular mass: ",meanmolmass," m_proton  (assuming atomic hydrogen)"
   end if
   write(iprint,'(a,Es10.3,a)')    "NICIL: Ion mass for CR's:   ",mass_ion_mp," m_proton"
   if (g_cnst) then 
     write(iprint,'(a,Es10.3,a)')  "NICIL: Grain mass:          ",mass_grain_cgs/mass_proton_cgs," m_proton"
     write(iprint,'(a,Es10.3,a)')  "NICIL: Grain radius:        ",a_grain," cm"
   else
     write(iprint,'(a,Es10.3,a)')  "NICIL: Average grain mass:  ",mass_grain_cgs/mass_proton_cgs," m_proton"
     write(iprint,'(a,Es10.3,a)')  "NICIL: Average grain radius:",a01_grain," cm"
   end if
   if (ion_rays) then 
     write(iprint,'(a,Es10.3,a)')  "NICIL: Cosmic ray ionisation rate: ",zeta_cgs," s^{-1}"
     write(iprint,'(a)')           "NICIL: Including ionisation from Cosmic rays"
   end if
   if (ion_thermal) then
     write(iprint,'(a)')           "NICIL: Including thermal ionisation"
   end if
   if (.not.ion_thermal .and. .not.ion_rays) then
     write(iprint,'(a)')           "NICIL: WARNING! No ionisation selectioned.  This will be ideal MHD!"
   end if
   write(iprint,'(2a)')            "NICIL: Non-ideal terms used: ",trim(ni_terms)
 end if
 !
end subroutine nicil_print_summary
!----------------------------------------------------------------------!
!
! If the default values are used, test to ensure that the correct eta's
! are obtained for five preset densities & magnetic fields
!
!----------------------------------------------------------------------!
subroutine nicil_test_eta(unit_density,unit_Bfield,unit_eta,cs,test_eta)
 real,   intent(in) :: unit_density,unit_Bfield,unit_eta
 real,   intent(in) :: cs
 logical,intent(in) :: test_eta
 integer            :: i, ierr
 real               :: rho, Bfield,Z_grain,eta_ohm,eta_hall,eta_ambi
 real               :: test_vals(7,6)
 real, parameter    :: tol = 1.0d-10
 !
 if (  test_eta                                         .and. &
      (abs(fdg-0.01)            < epsilon(fdg)        ) .and. &
      (abs(a_grain-1.0d-5)      < epsilon(a_grain)    ) .and. &
      (abs(rho_bulk-3.0)        < epsilon(rho_bulk)   ) .and. &
      (abs(zeta_cgs-1.0d-17)    < epsilon(zeta_cgs)   ) .and. &
      (abs(mass_ion_mp-24.3)    < epsilon(mass_ion_mp)) .and. & 
      (abs(delta_gn-1.3)        < epsilon(delta_gn)   ) .and. &
      (abs(sigmav_eH2-3.16d-11) < epsilon(sigmav_eH2) ) .and. &
      (abs(sigmav_eHe-7.08d-11) < epsilon(sigmav_eHe) ) .and. &
      (abs(pnH2-0.804)          < epsilon(pnH2)       ) .and. & 
      (abs(pnHe-0.207)          < epsilon(pnHe)       ) ) then
   write(iprint,'(a)') "NICIL: Testing eta-calculations"
   !--The array values to test
   !  test_vals (in CGS): rho,B,z_grain,eta_ohm,eta_hall,eta_ambi
   test_vals(1,:)= (/4.0d-16, 3.165d-3, dble(-0.67965359128367264) &
                   , 7.0710765231099207d11, -1.5901423978163131d18, 1.1026782373994166d19/) 
   test_vals(2,:)= (/4.0d-15, 5.623d-3, dble(-0.67742896292003218) & 
                   , 6.3602523224303701d12, -7.2940546547500122d17, 3.4013040910039188d18/) 
   test_vals(3,:)= (/4.0d-14, 1.000d-2, dble(-0.65488470475132110) &
                   , 7.1201363533821625d13, -2.0345421822999002d17, 1.1821377537178232d18/) 
   test_vals(4,:)= (/4.0d-13, 1.779d-2, dble(-0.44798087497272804) &
                   , 2.2474258865389345d15, -4.4830850928372320d16, 1.1992110013820787d18/) 
   test_vals(5,:)= (/4.0d-12, 3.165d-2, dble(-0.11686824645557491) &
                   , 1.3804136628980962d17,  3.0604552538848179d17, 4.3336523410589404d18/) 
   test_vals(6,:)= (/4.0d-11, 5.630d-2, dble(-0.01772546056530152) &
                   , 2.2840406304506004d18,  8.6492943686429164d18, 2.4280715535144219d19/) 
   test_vals(7,:)= (/4.0d-10, 1.000d-1, dble(-0.00192860503181799) &
                   , 2.4587096107352924d19,  1.0540931084049741d20, 5.8670864377822757d19/) 
   do i = 1,7
     rho      = test_vals(i,1)/unit_density
     Bfield   = test_vals(i,2)/unit_Bfield
     Z_grain  = Zgrain0
     call nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,cs,Z_grain)
     eta_ohm  = eta_ohm *unit_eta
     eta_hall = eta_hall*unit_eta
     eta_ambi = eta_ambi*unit_eta
     ierr     = 0
     if (abs(1.0-test_vals(i,3)/Z_grain ) > tol) ierr = ierr + 1
     if (abs(1.0-test_vals(i,4)/eta_ohm ) > tol) ierr = ierr + 1
     if (abs(1.0-test_vals(i,5)/eta_hall) > tol) ierr = ierr + 1
     if (abs(1.0-test_vals(i,6)/eta_ambi) > tol) ierr = ierr + 1
     if (ierr > 0) then
       write(iprint,'(a,I2,a,I8)') "NICIL: Testing of eta-calculations failed on ",ierr," of 4 values on set ",i
       call nicil_fatal('nicil_test_eta','eta calculation test failed','Z_grain p.d. ',abs(1.0-test_vals(i,3)/Z_grain ),.false.)
       call nicil_fatal('nicil_test_eta','eta calculation test failed','eta_ohm p.d. ',abs(1.0-test_vals(i,4)/eta_ohm ),.false.)
       call nicil_fatal('nicil_test_eta','eta calculation test failed','eta_hall p.d.',abs(1.0-test_vals(i,5)/eta_hall),.false.)
       call nicil_fatal('nicil_test_eta','eta calculation test failed','eta_ambi p.d.',abs(1.0-test_vals(i,6)/eta_ambi),.true. )
     end if
   end do
   write(iprint,'(a)') "NICIL: Testing of eta-calculations complete.  Passed all tests."
 else
   write(iprint,'(a)') "NICIL: Unable to test eta since default values are not used."
 end if
 !
end subroutine nicil_test_eta
!----------------------------------------------------------------------!
!
! This is the primary control routine for NICIL.  It will call the 
! required ionisation and non-ideal MHD routines to determine the 
! coefficients of the non-ideal MHD terms
!
!----------------------------------------------------------------------!
subroutine nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,cs,Z_grainR,n_electron,Z_grainT,densities_out,sigmas_out)
 real,           intent(out)   :: eta_ohm,eta_hall,eta_ambi
 real,           intent(in)    :: Bfield,rho,cs
 real, optional, intent(inout) :: Z_grainR,n_electron
 real, optional, intent(out)   :: Z_grainT,densities_out(5),sigmas_out(5)
 real                          :: ZgrainR,ZgrainT,nelectron,nelectronR,nionR,nionT,mionT,ngrainT,mionT_mp,ngrainR
 real                          :: sigmaOHPpa(5),n_densities(5)
 !
 !--Return constant coefficient version and exit
 if (eta_constant) then
   call nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
   if (present(Z_grainR)     ) Z_grainR      = 0.0
   if (present(Z_grainT)     ) Z_grainT      = 0.0
   if (present(densities_out)) densities_out = 0.0
   if (present(n_electron)   ) n_electron    = 0.0
   if (present(sigmas_out)   ) sigmas_out    = 0.0
   return
 end if
 !
 !--Copy grain charge and electron number density to proper values, if present
 if (present(Z_grainR)) then
   ZgrainR = Z_grainR
 else
   ZgrainR = Zgrain0
 end if
 if (present(n_electron)) then
   nelectron = n_electron
 else
   nelectron = nelectron0
 end if
 !
 !--Calcualte the grain charge and electron number density from cosmic rays
 if (ion_rays) then
   call nicil_ionR_get_Zgrain(ZgrainR,nelectronR,nionR,rho,cs)
 else
   ZgrainR    = 0.0
   nelectronR = 0.0
   nionR      = 0.0
 end if
 !
 !--Calculate the grain charge from thermal ionisaion and update electron number density
 if (ion_thermal) then
  if (.false.) then
   call nicil_ionT_get_ne(nelectron,nelectronR,nionT,mionT_mp,rho,cs)
   call nicil_ionT_get_ng(nelectron-nelectronR,nionT,mionT_mp,cs,ZgrainT,ngrainT)
  end if
  if (.false.) then
   call nicil_ionT_get_ne(nelectron,0.0,nionT,mionT_mp,rho,cs)
   call nicil_ionT_get_ng(nelectron,nionT,mionT_mp,cs,ZgrainT,ngrainT)
   nelectron = nelectron + nelectronR
  end if
  if (.true.) then
   ngrainR     = n_grain_coef*(rho*mump1)
   call nicil_ionT_get_ne(nelectron,nelectronR,nionT,mionT_mp,rho,cs)
   call nicil_ionT_get_ng(nelectron,nionT,mionT_mp,cs,ZgrainT,ngrainT,nionR+ZgrainR*ngrainR)
   !print*, nelectron-nelectronR,nionT,nionR+ZgrainR*ngrainR
  end if
   !print*, nelectron,nelectronR,ZgrainT,ngrainT,mionT_mp,mass_proton
   mionT = mionT_mp*mass_proton
 else
   nelectron = nelectronR
   nionT     = 0.0
   mionT     = 0.0
   ngrainT   = 0.0
   ZgrainT   = 0.0
 end if
 !
 !--Calculate the conductivities
 call nicil_ion_get_sigma(Bfield,rho,nelectron,nionR,ZgrainR,nionT,mionT,ngrainT,ZgrainT,cs,sigmaOHPpa,n_densities)
 !
 !--Calculate the coefficients
 call nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,cs,sigmaOHPpa)
 !
 !--Copy optional arrays to output, if requested
 if (present(Z_grainR)     ) Z_grainR      = ZgrainR
 if (present(Z_grainT)     ) Z_grainT      = ZgrainT
 if (present(densities_out)) densities_out = n_densities
 if (present(n_electron)   ) n_electron    = nelectron
 if (present(sigmas_out)   ) sigmas_out    = sigmaOHPpa
 !
 return
end subroutine nicil_get_eta
!----------------------------------------------------------------------!
!
! This will terminate the primary programme if there is an error in  
! NICIL and unreasonable values are initialised or calculated.  
!
!----------------------------------------------------------------------!
subroutine nicil_fatal(wherefrom,reason,var,val,kill_passed)
 character(len=*),           intent(in) :: wherefrom,reason
 character(len=*), optional, intent(in) :: var
 real,             optional, intent(in) :: val
 logical,          optional, intent(in) :: kill_passed
 logical                                :: kill_local
 !
 !--Print cause of fatal error
 if (present(var) .and. present(val)) then 
   write(iprint,'(7a,Es16.4)') 'NICIL: ',wherefrom,': ',reason,'. ',var,'=',val
 else
   write(iprint,'(4a)') 'NICIL: ',wherefrom,': ',reason
 end if
 !--Terminate programme if requested
 kill_local = .true.
 if (present(kill_passed)) kill_local = kill_passed
 if (kill_local) then 
   write(iprint,'( a)') 'NICIL: Aborting.'
   stop
 end if
 !
end subroutine nicil_fatal
!======================================================================!
! COSMIC RAY IONISATION-RELATED SUBROUTINES                            !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  This will calculate Z_grain iteratively, using the Newton–Raphson 
!  method. 
!  NOTE: By construction, n_grain = (n_e - n_i)/Z_g = m_n/m_g*f_dg*n
!+
!----------------------------------------------------------------------!
subroutine nicil_ionR_get_Zgrain(Z_grain,n_electron,n_ion,rho,cs)
 implicit none
 integer, parameter   :: ctrmax = 200
 real,    parameter   :: tol    = 1.0e-8
 integer              :: ctr
 real, intent(inout)  :: Z_grain
 real, intent(out)    :: n_electron,n_ion
 real, intent(in)     :: rho,cs
 real                 :: cs1, n_total1
 real                 :: Zrat,Zold,Znew
 real                 :: k_fac,expkZ,termkZ1,fatZ,fatZdZ
 logical              :: iterate
 !
 !--Initialise values
 ctr      = 0
 Zrat     = tol*2.0
 Zold     = Z_grain
 cs1      = 1.0/cs
 k_fac    = k_fac_coef*cs1**2
 n_total1 = 1.0/(rho*mump1)
 iterate  = .true. 
 !
 ! Perform the iterations
 do while ( iterate )
   if ( Zold > epsilon(Zold)**2 ) then 
     expkZ   = exp(k_fac*Zold)     
     termkZ1 = 1.0 / ( 1.0 + k_fac*Zold )
     fatZ    = Zold - Z_iter_coef*n_total1          &
             *( k_eg_coef1*cs1*termkZ1    - k_ig_coef1*cs1*expkZ )
     fatZdZ  = 1.0  + Z_iter_coef*n_total1*k_fac    &
             *( k_eg_coef1*cs1*termkZ1**2 + k_ig_coef1*cs1*expkZ )
   else if (Zold < -epsilon(Zold)**2 ) then
     expkZ   = exp(-k_fac*Zold)     
     termkZ1 = 1.0 / ( 1.0 - k_fac*Zold )
     fatZ    = Zold - Z_iter_coef*n_total1       &
             *( k_eg_coef1*cs1*expkZ   - k_ig_coef1*cs1*termkZ1    )
     fatZdZ  = 1.0  + Z_iter_coef*n_total1*k_fac &
             *( k_eg_coef1*cs1*expkZ   + k_ig_coef1*cs1*termkZ1**2 )
   else
     expkZ   = 0.0
     termkZ1 = 0.0
     Znew    = 0.0
     iterate = .false.
   end if
   if ( iterate ) then 
     Znew    = Zold - fatZ/fatZdZ
     Zrat    = abs( 1.0 - Znew/Zold )
     Zold    = Znew
     ctr     = ctr + 1
     if (ctr >= ctrmax .or. Zrat < tol) iterate = .false.
   end if
 end do
 if (ctr >= ctrmax) call nicil_fatal('nicil_ion_get_Zgrain','Z_grain did not converge')
 !
 !--The new value of Z & n_electron
 Z_grain    = Znew 
 if (Z_grain >= 0) then 
   n_electron = n_ele_coef*cs1*termkZ1
   n_ion      = n_ion_coef*cs1*expkZ
 else
   n_electron = n_ele_coef*cs1*expkZ
   n_ion      = n_ion_coef*cs1*termkZ1
 end if
 !
 !--Verifications against absurb number densities
 if (n_electron < 0.0 ) then
   call nicil_fatal('nicil_ion_get_Zgrain','unrealistic n_elctron','n_electron',n_electron)
 end if
 if (n_ion < 0.0 ) then
   call nicil_fatal('nicil_ion_get_Zgrain','unrealistic n_ion','n_ion',n_ion)
 end if
 !
end subroutine nicil_ionR_get_Zgrain
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
subroutine nicil_ion_get_sigma(Bmag,rho,n_electron,n_ionR,Z_grainR,n_ionT,m_ionT,n_grainT,Z_grainT,cs,sigmaOHPpa,n_densities)
 real,              intent(out) :: sigmaOHPpa(5),n_densities(5)
 real,              intent(in)  :: Bmag,rho,cs,n_electron,n_ionR,Z_grainR,n_ionT,m_ionT,n_grainT,Z_grainT
 integer                        :: j,k
 real                           :: n_grainR,rho_n,rho_iR,rho_iT,rho_e,cs1
 real                           :: sigmavenX,sigmavenY,sigma_coef_onB
 real                           :: nu_ei,nu_en,nu_inR,nu_inT,nu_gn,nu_inT_coef0
 real                           :: beta_e,beta_iR,beta_iT,beta_gR,beta_gT
 real                           :: betae2p11,betaiR2p11,betaiT2p11,betagR2p11,betagT2p11
 real                           :: nj(nspecies),Zj(nspecies),aZj(nspecies),sZj(nspecies)
 real                           :: betaj(nspecies), beta2p11(nspecies)
 !
 !--Initialise values
 cs1            = 1.0/cs
 sigma_coef_onB = sigma_coef/Bmag
 !
 !--Densities 
 n_grainR     = n_grain_coef*(rho*mump1)
 rho_iR       = n_ionR*mass_ion
 rho_iT       = n_ionT*m_ionT
 rho_e        = n_electron*mass_electron
 rho_n        = rho - (rho_iR + rho_iT + rho_e) 
 !rho_n        = rho - (rho_iR + rho_e) 
 !print*, rho, rho_n,rho_iR,rho_iT,m_ionT
 if (rho_n < tiny(rho_n)) then
   ! This is the Ideal MHD regieme.  Turn off non-ideal terms.  
   rho_n      =  0.0
   sigmaOHPpa = -1.0
 else
   !
   !--Rate Coefficients
   sigmavenX = sigmavenXbyT*cs**1.3
   sigmavenY = sigmavenYbyT*cs
   if (sigmavenX > sigmavenX_LR) sigmavenX = sigmavenX_LR
   if (sigmavenY > sigmavenX_LR) sigmavenY = sigmavenX_LR
   !
   !--Collisional Frequencies
   call nicil_extract_nuinTcoef(nu_inT_coef0,m_ionT)
   nu_ei      = nu_ei_coef  *n_electron*cs1**3 
   nu_inR     = nu_inR_coef *rho_n   
   nu_inT     = nu_inT_coef0*rho_n   
   nu_en      = nu_en_coef  *rho_n*(sigmavenX + sigmavenY)
   nu_gn      = nu_gn_coef  *rho_n*cs   
   !
   !--Hall parameters
   if (ion_rays .and. rho_iR > 0.0) then
     beta_iR    = beta_iR_coef*          Bmag/(nu_inR + nu_ei*rho_e/rho_iR)  !recall nu_ie = nu_ei*rho_e/rho_i
     betaiR2p11 = 1.0/(1.0+beta_iR**2)
   else
     beta_iR    = 0.0
     betaiR2p11 = 0.0
   end if
   if (ion_thermal .and. rho_iT > 0.0) then
     beta_iT    = beta_iT_coef*          Bmag/(m_ionT*(nu_inT + nu_ei*rho_e/rho_iT))  !recall nu_ie = nu_ei*rho_e/rho_i
     betaiT2p11 = 1.0/(1.0+beta_iT**2)
   else
     beta_iT    = 0.0
     betaiT2p11 = 0.0
   end if
   beta_e       = beta_e_coef *              Bmag/(nu_en  + nu_ei)
   beta_gR      = beta_g_coef *abs(Z_grainR)*Bmag/nu_gn
   beta_gT      = beta_g_coef *abs(Z_grainT)*Bmag/nu_gn
   betae2p11    = 1.0/(1.0+beta_e**2 )
   betagR2p11   = 1.0/(1.0+beta_gR**2)
   betagT2p11   = 1.0/(1.0+beta_gT**2)
   !--fill in arrays for simplicity when calculating conductivities
   nj       = (/n_electron, n_ionR,    n_ionT,    n_grainR,      n_grainT     /)
   Zj       = (/-1.0,       1.0,       1.0,           Z_grainR,      Z_grainT /)
   aZj      = (/ 1.0,       1.0,       1.0,       abs(Z_grainR), abs(Z_grainT)/)
   sZj      = (/ 1.0,       1.0,       1.0,       sign(1.0,Z_grainR), sign(1.0,Z_grainT)/)
   betaj    = (/beta_e   ,  beta_iR,   beta_iT,   beta_gR,       beta_gT      /)
   beta2p11 = (/betae2p11,  betaiR2p11,betaiT2p11,betagR2p11,    betagT2p11   /)
   !print*, beta_e, beta_iR,beta_iT,beta_gR,beta_gT,betae2p11,betagR2p11,betagT2p11,aZ_grainT
   !
   !--Conductivities
   sigmaOHPpa    = 0.0
   do j = 1,nspecies
     sigmaOHPpa(1) = sigmaOHPpa(1) + nj(j)*aZj(j)*betaj(j)
     sigmaOHPpa(2) = sigmaOHPpa(2) + nj(j)* Zj(j)         *beta2p11(j)
     sigmaOHPpa(3) = sigmaOHPpa(3) + nj(j)*aZj(j)*betaj(j)*beta2p11(j)
     do k = j+1,nspecies
       sigmaOHPpa(5) = sigmaOHPpa(5) + nj(j)*aZj(j)*betaj(j)*beta2p11(j) &
                                     * nj(k)*aZj(k)*betaj(k)*beta2p11(k) &
                                     * (sZj(j)*betaj(j) - sZj(k)*betaj(k))**2
     end do
   end do
   sigmaOHPpa(1:3) = sigmaOHPpa(1:3)*sigma_coef_onB
   sigmaOHPpa(5)   = sigmaOHPpa(  5)*sigma_coef_onB**2
   sigmaOHPpa(4)   = 1.0/(sigmaOHPpa(2)**2 + sigmaOHPpa(3)**2)
 end if
 !print*, sigmaOHPpa, n_electron,n_ionR,n_ionT,n_grainR,n_grainT
 !
 !--Number densities (for bookkeeping)
 n_densities(1) = rho_n*mass_neutral1 ! total number density
 n_densities(2) = n_ionR              ! ion number density
 n_densities(3) = n_ionT              ! ion number density
 n_densities(4) = n_grainR            ! grain number density
 n_densities(5) = n_grainT            ! grain number density
 !
end subroutine nicil_ion_get_sigma
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
!  NOTE: We cannot calculate this at the beginning of the step since
!        in this form, it requires ne_zeta!
!+
!----------------------------------------------------------------------!
subroutine nicil_ionT_get_ne(ne_saha,ne_zeta,nion,mion_mp,rho,cs)
 implicit none
 integer, parameter   :: ctrmax = 200
 real,    parameter   :: tol    = 1.0e-8
 integer              :: j,ctr
 real, intent(inout)  :: ne_saha
 real, intent(out)    :: mion_mp,nion
 real, intent(in)     :: ne_zeta,rho,cs
 real                 :: fatn,fatndn,nerat,neold,nenew,termj
 real                 :: Kj(nelements),nj(nelements),nij(nelements)
 logical              :: iterate
 !
 !--Initialise values
 ctr        = 0
 nerat      = tol*2.0
 neold      = ne_saha
 iterate    = .true. 
 !  The 'constant' part of the Saha equation & the species number densities
 !  note that chij has been modified to include mump1
 do j = 1,nelements
   Kj(j) = gej(j)*Saha_coef*cs**3*exp(-chij(j)/cs**2)
   nj(j) = abundancej(j)*rho*mump1
 end do
 !
 !--Calcualte n_e
 do while ( iterate )
   fatn   = neold - ne_zeta 
   fatndn = 1.0  
   do j = 1,nelements
     termj  = nj(j)*Kj(j)/(neold + Kj(j))
     fatn   = fatn   - termj
     fatndn = fatndn + termj/(neold + Kj(j))
   end do
   if ( iterate ) then 
     nenew    = neold - fatn/fatndn
     nerat    = abs( 1.0 - nenew/neold )
     neold    = nenew
     ctr      = ctr + 1
     if (ctr >= ctrmax .or. nerat < tol) iterate = .false.
     !print*, neold, ne_zeta, nenew,ne_saha
   end if
 end do
 if (ctr >= ctrmax) call nicil_fatal('nicil_thermal_get_ne','n_e did not converge')
 !
 !--Get the ion number density (total and per species)
 nion    = 0.0
 mion_mp = 0.0
 do j = 1,nelements
   nij(j)   = nj(j)*Kj(j)/(nenew + Kj(j))
   nion    = nion    + nij(j)
   mion_mp = mion_mp + nij(j)/sqrt(mj(j))
 end do
 if (mion_mp > 0.0) mion_mp = (nion/mion_mp)**2
! write(89,*) rho*unit_density0,mion_mp, nion*unit_density0/umass0,nenew*unit_density0/umass0
 ne_saha = nenew
! write(88,*) ne_saha, nenew,ne_zeta,cs,cs**3,Saha_coef!,exp(-chij/cs**2),chij
 !
end subroutine nicil_ionT_get_ne
!----------------------------------------------------------------------!
! Calculate grain number density and charge
subroutine nicil_ionT_get_ng(ne_saha,nion,mion_mp,cs,Z_g_ionT,n_g_ionT,nZR)
 implicit none
 real, intent(in)     :: ne_saha,mion_mp,nion,cs
 real, intent(in), optional :: nZR
 real, intent(out)    :: Z_g_ionT,n_g_ionT
 real                 :: mu_ionT,y_ionT,lnyp1,Psi,tau0,tau,dn
 !
 if (nion > 0.0) then 
   mu_ionT  = (se*ne_saha/nion)**2*mion_mp
   y_ionT   = y_ionT_coef*sqrt(mu_ionT) 
   lnyp1    = log(y_ionT + 1.0)
   Psi      = 1.0 - lnyp1 + lnyp1/(1.0+lnyp1)*log( (1.0+1.0/y_ionT)*lnyp1 )
   tau0     = tau0_coef/mu_ionT
   tau      = tau_coef*cs**2*mu_ionT
   Z_g_ionT = Psi*tau - 1.0/(1.0 + sqrt(tau0/tau))
   dn       = ne_saha - nion
   if (present(nZR)) dn = ne_saha - nion - nZR
!   dn = 1.0 - nion/ne_saha
!   if ( dn > 2.0*max(ne_saha,nion)*epsilon(ne_saha - nion)) then 
   if ( dn > epsilon(dn)) then 
 !     n_g_ionT = ne_saha*dn/Z_g_ionT
     n_g_ionT = dn/Z_g_ionT
   else
     Z_g_ionT = 0.0
     n_g_ionT = 0.0
     dn = 0.0
   end if
  ! write(86,*) ne_saha, nion, dn, dn*nion,n_g_ionT,Z_g_ionT
   !print*, mu_ionT,y_ionT,lnyp1,Psi,tau0,tau,Z_g_ionT,n_g_ionT
 else
   Z_g_ionT = 0.0
   n_g_ionT = 0.0
 end if
! write(87,*) Z_g_ionT,n_g_ionT

 !
end subroutine nicil_ionT_get_ng
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
pure subroutine nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,cs,sigmaOHPpa)
 implicit none
 real,    intent(out) :: eta_ohm,eta_hall,eta_ambi
 real,    intent(in)  :: cs,sigmaOHPpa(5)
 !
 eta_ohm  = 0.0
 eta_hall = 0.0
 eta_ambi = 0.0
 if (cs < cs_thresh .and. sigmaOHPpa(1) > 0.0) then
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
 real,    intent(out) :: eta_ohm,eta_hall,eta_ambi
 real,    intent(in)  :: Bfield,rho
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
pure subroutine nimhd_get_sphdBdt(dBnonideal,etaohm,etahall,etaambi &
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
end subroutine nimhd_get_sphdBdt
!----------------------------------------------------------------
!+
!  performs simple cross product 
!+
!----------------------------------------------------------------
pure subroutine nimhd_get_DcrossR(DcrossR,D_in,dx,dy,dz,eta)
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
!----------------------------------------------------------------
!+
!  Calculates the timesteps for non-ideal MHD
!+
!----------------------------------------------------------------
pure subroutine nimhd_get_dt(dtohm,dthall,dtambi,h,etaohm,etahall,etaambi)
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
