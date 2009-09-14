c
c  Based on work by G. Suttner (Univ. Wuerzburg, 1995) and  
c  M. D. Smith (Armagh Observatory, 2000-2001).
c  Substantially modified by S. Glover (AMNH, 2002-2005, AIP 2006)
c

c On some machines we may want our reals to be real*4
c In practice, I have never had cause to use anything 
c other than real*8
c
#define REAL real*8

c ABORT is called if something is broken beyond repair
c in the cooling code. 99% of the time, this is because
c the temperature has become NaN, at which point we 
c may as well just quit...
c
c If your code requires cleanup when terminating, then
c replace the simple 'stop' here with a call to the
c appropriate function
#define ABORT(x) stop

c Boltzmann constant
      REAL kboltz
      parameter (kboltz = 1.38066e-16)

c One electron volt, in ergs
      REAL eV
      parameter (eV = 1.60219e-12)
c
c He:H ratio by number (=> ratio by mass is 4*abhe)
c
      REAL abhe
      parameter(abhe = 0.1d0)

c Number of entries in cooling table
      integer nmd
      parameter(nmd = 10000)

c Number of cooling / heating rates computed in cooling fn.
      integer nrates
      parameter(nrates = 11)

c Size of abundance array that is passed to cool_func as input
      integer nabn
      parameter(nabn = 9)

c Number of different quantities stored in cooling look-up table
      integer ncltab
      parameter (ncltab = 54) 

c These variables are initialized in coolinmo
      REAL temptab(nmd)
      REAL cltab(ncltab,nmd), dtcltab(ncltab,nmd)
      REAL dtlog, tmax, tmin
c
c These variables must be initialized during problem setup
c
c Total abundances of C, O, Si:
      REAL abundc
      parameter (abundc = 2d-4)      

      REAL abundo
      parameter (abundo = 4.5d-4)

      REAL abundsi
      parameter (abundsi = 3d-5)
c
c Strength of UV field (in Habing units):
      REAL G0
      parameter (G0 = 1d0)
c
c  Dust temperature (in K)
      REAL tdust
      parameter (tdust = 1d1)
c
c Dust:gas ratio in units of standard Galactic value.
c
c [At lower metallicities, I generally assume for 
c  simplicity that this scales linearly with metallicity,
c   but it need not do so]
c
      REAL  dust_to_gas_ratio
      parameter (dust_to_gas_ratio  = 1d0)
c
c Visual extinction (A_V) per unit column density (in cm^-2)
      REAL AV_conversion_factor
      parameter (AV_conversion_factor = 5.348d-22)
c
c Cosmic ray ionization rate of HI (in s^-1)
      REAL cosmic_ray_ion_rate
      parameter (cosmic_ray_ion_rate  = 1d-17)

c
c Flag controlling treatment of photoelectric heating
c iphoto = 0 ==> optically thin gas assumed
c iphoto = 1 ==> approximate treatment of effects of extinction
c [See cool_func.F for details]
c
      integer iphoto
      parameter (iphoto = 1)
c
c Flag controlling which atomic cooling function is used.
c iflag_atom = 1 is the appropriate choice for the Galactic ISM
c [iflag_atom = 2 is used for Z=0 gas]
      integer iflag_atom
      parameter (iflag_atom = 1)
c
      common /coolr/ temptab, cltab, dtcltab, dtlog, tmax, 
     $               tmin

