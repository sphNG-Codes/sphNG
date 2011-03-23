c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            C O O L _ F U N C              \\\\\\\\\\
c
c=======================================================================
c
      subroutine cool_func(temp, yn, dl, abundances, ylam,rates)
c                                                                       
c    Based on cool_h, written for ZEUS-3D by Michael D. Smith and
c    Georgi Pavlovski (Armagh Observatory, 2003) and substantially
c    modified by S. Glover (AMNH, 2003-2005, AIP 2006).
c
c    PURPOSE:  Compute the cooling function for the gas. Note that 
c    the convention used here is that Lambda (the net _cooling_ rate)
c    is positive and has units of erg s^-1 cm^-3, so that for gas at
c    rest the equation for the internal energy density, e, can be 
c    written as:   
c                 de/dt = - Lambda
c    This means that any heating terms which we include in Lambda
c    (e.g. photoelectric heating) must be NEGATIVE; i.e. heating is
c    treated as negative cooling.
c
c    N.B. When used for cosmological simulations, it is necessary to
c         account for the effects of the CMB on the various heating & 
c         cooling rates. In the case of the fine structure coolants,
c         this is done explicitly (and exactly) by including the 
c         appropriate stimulated emission & absorption terms in the 
c         level population equations. For other species (HD, CO, H2O),
c         this effect is not treated explicitly -- instead, we 
c         approximate it by first calling cool_func with the gas 
c         temperature and then with the CMB temperature, and then 
c         subtracting the cooling rate given by the second call from 
c         that given by the first. This approximation is of questionable
c         accuracy near T_cmb -- for instance, it can lead to ~50% errors
c         in the HD cooling rate -- but is probably sufficient for our 
c         present purposes, given the other uncertainties
c                  
c  INPUT VARIABLES: temp       -- temperature
c                   yn         -- density
c                   divv       -- velocity divergence
c                   abundances -- chemical abundances
c
c  OUTPUT VARIABLES: Computed values of ylam & rates
c
c  LOCAL VARIABLES:  rates  -- parts of the cooling function
c                       cl[1-*] -- cooling table parameters
c                     dtcl[1-*] -- their change with T = T + dT
c
c  NB: the dt... terms are TEMPERATURE derivatives, not time derivatives
c
c  EXTERNALS:
c
c-----------------------------------------------------------------------
c
      implicit none
#include "cool.h"
#include "fs_data.h"

      integer    itemp
      REAL       temp      , dtemp   , yn      , dl

      REAL       abundances(nabn)
      REAL       ynh2      , ynh     , yne     , ynhp
      REAL       rates(nrates)
c
      REAL       ylam
c
      REAL       abh2      , abo       ,
     $           abcI      , abcII     , absiI   , absiII   , abe  ,
     $           abhp      , abhI
c
      REAL       h2var0    , h2var1
c
      REAL       oxc01, oxc02,
     $           oxc10, oxc20, oxn1  , oxn2  , oxn0  , oxr01 ,
     $           oxr02, oxr12, oxr10 , oxr20 , oxr21 , oxa   ,
     $           oxb  , oxc  , oxc21 , oxc12
c
      REAL       cIc01, cIc02, cIc10 ,
     $           cIc20, cIn1 , cIn2  , cIn0  , cIr01 ,
     $           cIr02, cIr12, cIr10 , cIr20 , cIr21 , cIa   ,
     $           cIb  , cIc  , cIc21 , cIc12
c
      REAL       siIc01, siIc02,  siIc10,
     $           siIc20, siIn1 , siIn2  , siIn0 , siIr01 ,
     $           siIr02, siIr12, siIr10 , siIr20 , siIr21 , siIa,
     $           siIb  , siIc  , siIc21 , siIc12
c
      REAL       cIIc10,   cIIc01, cIIn1
c
      REAL       siIIc10, siIIc01, siIIn1
c
      REAL       eps       , PEvar0   , PEvar1   , PEvar2,
     $           G_dust    , AV       , fdust    , NH_tot
c
      REAL       cl1       , cl2     , cl3     , cl6     , cl14
     $         , cl16      , cl17    , cl18    , cl19    , cl20
     $         , cl21      , cl22    , cl23    , cl24    , cl25
     $         , cl26      , cl27    , cl28    , cl29    , cl30
     $         , cl31      , cl32    , cl33    , cl34    , cl35
     $         , cl36      , cl37    , cl38    , cl39    , cl40
     $         , cl41      , cl42    , cl43    , cl44    , cl45
     $         , cl46      , cl47    , cl48    , cl49    , cl50
     $         , cl51      , cl52    , cl53    , cl54
c
      REAL       dtcl1     , dtcl2   , dtcl3   , dtcl6   , dtcl14
     $         , dtcl16    , dtcl17  , dtcl18  , dtcl19  , dtcl20
     $         , dtcl21    , dtcl22  , dtcl23  , dtcl24  , dtcl25
     $         , dtcl26    , dtcl27  , dtcl28  , dtcl29  , dtcl30
     $         , dtcl31    , dtcl32  , dtcl33  , dtcl34  , dtcl35
     $         , dtcl36    , dtcl37  , dtcl38  , dtcl39  , dtcl40
     $         , dtcl41    , dtcl42  , dtcl43  , dtcl44  , dtcl45
     $         , dtcl46    , dtcl47  , dtcl48  , dtcl49  , dtcl50
     $         , dtcl51    , dtcl52  , dtcl53  , dtcl54
c
c ---------------------------------------------------------------------
c
c Read out tables.
c
      if (temp .le. tmin) then
        itemp = 1
        dtemp = 0d0
      elseif (temp .gt. tmax) then
        itemp = nmd
        dtemp = temp - temptab(itemp)
      else
        itemp = int(dlog10(temp) / dtlog) + 1
        if (itemp .le. 0 .or. itemp .gt. nmd) then
          print*, 'Fatal error in cool_func.F', itemp, temp
          ABORT(0)
        endif
        dtemp = temp - temptab(itemp)
      endif
c
      dtcl1  = dtcltab(1,itemp)
      dtcl2  = dtcltab(2,itemp)
      dtcl3  = dtcltab(3,itemp)
c
      dtcl6  = dtcltab(6,itemp)
c
      dtcl14 = dtcltab(14,itemp)             
c
      dtcl16 = dtcltab(16,itemp)
      dtcl17 = dtcltab(17,itemp)
      dtcl18 = dtcltab(18,itemp)
      dtcl19 = dtcltab(19,itemp)
      dtcl20 = dtcltab(20,itemp)
      dtcl21 = dtcltab(21,itemp)
      dtcl22 = dtcltab(22,itemp)
      dtcl23 = dtcltab(23,itemp)
      dtcl24 = dtcltab(24,itemp)
      dtcl25 = dtcltab(25,itemp)
      dtcl26 = dtcltab(26,itemp)
      dtcl27 = dtcltab(27,itemp)
      dtcl28 = dtcltab(28,itemp)
      dtcl29 = dtcltab(29,itemp)
      dtcl30 = dtcltab(30,itemp)
      dtcl31 = dtcltab(31,itemp)
      dtcl32 = dtcltab(32,itemp)
      dtcl33 = dtcltab(33,itemp)
      dtcl34 = dtcltab(34,itemp)
      dtcl35 = dtcltab(35,itemp)
      dtcl36 = dtcltab(36,itemp)
      dtcl37 = dtcltab(37,itemp)
      dtcl38 = dtcltab(38,itemp)
      dtcl39 = dtcltab(39,itemp)
      dtcl40 = dtcltab(40,itemp)
      dtcl41 = dtcltab(41,itemp)
      dtcl42 = dtcltab(42,itemp)
      dtcl43 = dtcltab(43,itemp)
      dtcl44 = dtcltab(44,itemp)
      dtcl45 = dtcltab(45,itemp)
      dtcl46 = dtcltab(46,itemp)
      dtcl47 = dtcltab(47,itemp)
      dtcl48 = dtcltab(48,itemp)
      dtcl49 = dtcltab(49,itemp)
      dtcl50 = dtcltab(50,itemp)
      dtcl51 = dtcltab(51,itemp)
      dtcl52 = dtcltab(52,itemp)
      dtcl53 = dtcltab(53,itemp)
      dtcl54 = dtcltab(54,itemp)

      cl1  = cltab(1,itemp) + dtemp * dtcl1
      cl2  = cltab(2,itemp) + dtemp * dtcl2
      cl3  = cltab(3,itemp) + dtemp * dtcl3
c
      cl6  = cltab(6,itemp) + dtemp * dtcl6
c
      cl14 = cltab(14,itemp) + dtemp * dtcl14
c
      cl16 = cltab(16,itemp) + dtemp * dtcl16
      cl17 = cltab(17,itemp) + dtemp * dtcl17
      cl18 = cltab(18,itemp) + dtemp * dtcl18
      cl19 = cltab(19,itemp) + dtemp * dtcl19
      cl20 = cltab(20,itemp) + dtemp * dtcl20
      cl21 = cltab(21,itemp) + dtemp * dtcl21
      cl22 = cltab(22,itemp) + dtemp * dtcl22
      cl23 = cltab(23,itemp) + dtemp * dtcl23
      cl24 = cltab(24,itemp) + dtemp * dtcl24
      cl25 = cltab(25,itemp) + dtemp * dtcl25
      cl26 = cltab(26,itemp) + dtemp * dtcl26
      cl27 = cltab(27,itemp) + dtemp * dtcl27
      cl28 = cltab(28,itemp) + dtemp * dtcl28
      cl29 = cltab(29,itemp) + dtemp * dtcl29
      cl30 = cltab(30,itemp) + dtemp * dtcl30
      cl31 = cltab(31,itemp) + dtemp * dtcl31
      cl32 = cltab(32,itemp) + dtemp * dtcl32
      cl33 = cltab(33,itemp) + dtemp * dtcl33
      cl34 = cltab(34,itemp) + dtemp * dtcl34
      cl35 = cltab(35,itemp) + dtemp * dtcl35
      cl36 = cltab(36,itemp) + dtemp * dtcl36
      cl37 = cltab(37,itemp) + dtemp * dtcl37
      cl38 = cltab(38,itemp) + dtemp * dtcl38
      cl39 = cltab(39,itemp) + dtemp * dtcl39
      cl40 = cltab(40,itemp) + dtemp * dtcl40
      cl41 = cltab(41,itemp) + dtemp * dtcl41
      cl42 = cltab(42,itemp) + dtemp * dtcl42
      cl43 = cltab(43,itemp) + dtemp * dtcl43
      cl44 = cltab(44,itemp) + dtemp * dtcl44
      cl45 = cltab(45,itemp) + dtemp * dtcl45
      cl46 = cltab(46,itemp) + dtemp * dtcl46
      cl47 = cltab(47,itemp) + dtemp * dtcl47
      cl48 = cltab(48,itemp) + dtemp * dtcl48
      cl49 = cltab(49,itemp) + dtemp * dtcl49
      cl50 = cltab(50,itemp) + dtemp * dtcl50
      cl51 = cltab(51,itemp) + dtemp * dtcl51
      cl52 = cltab(52,itemp) + dtemp * dtcl52
      cl53 = cltab(53,itemp) + dtemp * dtcl53
      cl54 = cltab(54,itemp) + dtemp * dtcl54
c
c Set abundances
c
      abh2   = abundances(1)
      abhI   = abundances(2)
      abe    = abundances(3) 
      abhp   = abundances(4)
      abo    = abundances(5)
      abcI   = abundances(6)
      abcII  = abundances(7)
      absiI  = abundances(8)
      absiII = abundances(9)
c
c Compute useful auxiliary variables

      ynh2 = abh2 * yn
      ynh  = abHI * yn
      yne  = abe  * yn
      ynhp = abhp * yn
c
c (R1) -- gas-grain cooling-heating -- dust:gas ratio already incorporated
c         into rate coefficient in coolinmo
c
      rates(1) = cl14 * (temp - tdust) * yn**2
c
c
c (R2) -- H2 (vr) cooling
c
      if (ynh2 .eq. 0d0) then
        rates(2) = 0d0
      else
        h2var0   = ynh * cl2 + ynh2 * cl3
        h2var1   = cl1 + h2var0
        rates(2) = ynh2 * cl1 * h2var0 / h2var1
      endif
c
c (R3) -- atomic cooling
c
      rates(3) = cl6 * yn**2 + cl52 * ynh * yne
c
c (R4) -- cosmic ray heating; independent of gas temperature
c       -- following Goldsmith & Langer (1978), we assume that 
c          each ionization deposits 20eV as heat
c
      rates(4) =  -3.2d-11 * (ynh2 + ynh) * cosmic_ray_ion_rate
c
c (R5-6) -- photoelectric heating, PAH recombination cooling
c
c If there's no UV field, or if the electron density is very low (in which
c case the photoheating efficiency will also be very low), then we set the
c rates to zero. Otherwise, we compute the heating rate using the 
c Bakes & Tielens (1994) formula (as modified by Wolfire et al, 2003),
c and the recombination cooling using the formula from Wolfire et al 2003.
c
c Treatment of dust attenuation follows Bergin et al (2004, ApJ,  612, 921)
c Calculation of A_V uses the conversion factor that was initialized during
c problem setup. This treatment neglects any variation of the photoelectric 
c heating efficiency with increasing A_V.
c
      if (G0 .eq. 0d0 .or. dust_to_gas_ratio .eq. 0d0) then
        rates(5)  = 0d0
        rates(6) = 0d0
      elseif (yne .lt. 1d-9 * G0 * cl16) then
        rates(5)  = 0d0
        rates(6) = 0d0
      else
        if (iphoto .eq. 0) then
          G_dust = G0
        else
          NH_tot = 0.5d0 * dl * yn * (2d0 * abh2 + abhI + abhp)
          AV     = AV_conversion_factor * dust_to_gas_ratio * NH_tot
          fdust  = dexp(-2.5d0 * AV)
          G_dust = G0 * fdust
        endif

        PEvar0 = G_dust * cl16 / yne
        PEvar1 = (1d0 + 4d-3 * PEvar0**0.73)
        PEvar2 = (1d0 + 2d-4 * PEvar0)

        eps = (4.9d-2 / PEvar1) + (cl17 / PEvar2)

c Photoelectric heating:
        rates(5) = -1.3d-24 * eps * G_dust * yn * dust_to_gas_ratio
c Recombination cooling:
        rates(6) = cl54 * PEvar0**cl53 * yne * yn * dust_to_gas_ratio
      endif
c
c (R7) -- OI fine-structure cooling
c
c Total collisional rates:
c
      oxc10  = cl18 * ynh + cl21 * ynh2 + cl24 * yne + cl27 * ynhp
      oxc20  = cl19 * ynh + cl22 * ynh2 + cl25 * yne + cl28 * ynhp
      oxc21  = cl20 * ynh + cl23 * ynh2 + cl26 * yne + cl29 * ynhp

      if (abo .le. 1d-5 * abundo) then
        rates(7) = 0d0
      elseif (oxc10 .eq. 0d0 .and. oxc20 .eq. 0d0 .and.
     $    oxc21 .eq. 0d0) then
        rates(7) = 0d0
      else
        oxa =  dexp(-oxe10 / (kboltz * temp))
        oxb =  dexp(-oxe20 / (kboltz * temp))
        oxc =  dexp(-oxe21 / (kboltz * temp))

        oxc01  = 0.6d0 * oxc10 * oxa
        oxc02  = 0.2d0 * oxc20 * oxb
        oxc12  = (1d0 / 3d0) * oxc21 * oxc
c
c Total transition rates:
c
        oxR01  = oxc01
        oxR02  = oxc02
        oxR12  = oxc12
        oxR10  = oxc10 + oxa10
        oxR20  = oxc20 + oxa20
        oxR21  = oxc21 + oxa21
c
        call three_level_pops(oxR01, oxR02, oxR12, oxR10, oxR20,
     $                        oxR21, oxn0, oxn1, oxn2)
c
c Total emitted energy:
c
        rates(7) = (oxa10 * oxe10 * oxn1 + (oxa20 * oxe20 + 
     $               oxa21 * oxe21) * oxn2) * abo * yn
      endif
c
c (R8) --  CI fine-structure cooling 
c
c Collisional rates:
c
      cIc10  = cl30 * ynh + cl33 * ynh2 + cl36 * yne + cl39 * ynhp
      cIc20  = cl31 * ynh + cl34 * ynh2 + cl37 * yne + cl40 * ynhp
      cIc21  = cl32 * ynh + cl35 * ynh2 + cl38 * yne + cl41 * ynhp

      if (abcI .le. 1d-5 * abundc) then
        rates(8) = 0d0
      elseif (cIc10 .eq. 0d0 .and. cIc20 .eq. 0d0 .and.
     $    cIc21 .eq. 0d0) then
        rates(8) = 0d0
      else
        cIa =  dexp(-cIe10 / (kboltz * temp))
        cIb =  dexp(-cIe20 / (kboltz * temp))
        cIc =  dexp(-cIe21 / (kboltz * temp))

        cIc01  = 3d0 * cIc10 * cIa
        cIc02  = 5d0 * cIc20 * cIb
        cIc12  = (5d0 / 3d0) * cIc21 * cIc
c
c Total transition rates:
c
        cIR01  = cIc01
        cIR02  = cIc02
        cIR12  = cIc12
        cIR10  = cIc10 + cIa10
        cIR20  = cIc20 + cIa20
        cIR21  = cIc21 + cIa21
c
        call three_level_pops(cIR01, cIR02, cIR12, cIR10, cIR20,
     $                        cIR21, cIn0, cIn1, cIn2)
c
c Total emitted energy:
c
        rates(8) = (cIa10 * cIe10 * cIn1 + (cIa20 * cIe20 +
     $               cIa21 * cIe21) * cIn2) * abcI * yn
      endif
c
c (R9) --  SiI fine-structure cooling 
c
c Proton rates (from HM89) are constant and so there's no point
c tabulating them in coolinmo.
c
      siIc10 = cl42 * ynh + 7.2d-9 * ynhp
      siIc20 = cl43 * ynh + 7.2d-9 * ynhp
      siIc21 = cl44 * ynh + 2.2d-8 * ynhp

      if (absiI .le. 1d-5 * abundsi) then
        rates(9) = 0d0
      elseif (siIc10 .eq. 0d0 .and. siIc20 .eq. 0d0 .and.
     $    siIc21 .eq. 0d0) then
        rates(9) = 0d0
      else
        siIa =  dexp(-siIe10 / (kboltz * temp))
        siIb =  dexp(-siIe20 / (kboltz * temp))
        siIc =  dexp(-siIe21 / (kboltz * temp))
c
        siIc01  = 3d0 * siIc10 * siIa
        siIc02  = 5d0 * siIc20 * siIb
        siIc12  = (5d0 / 3d0) * siIc21 * siIc
c
c Total transition rates:
c
        siIR01  = siIc01
        siIR02  = siIc02
        siIR12  = siIc12
        siIR10  = siIc10 + siIa10
        siIR20  = siIc20 + siIa20
        siIR21  = siIc21 + siIa21
c
        call three_level_pops(siIR01, siIR02, siIR12, siIR10,
     $                        siIR20, siIR21, siIn0, siIn1, siIn2)
c
c Total emitted energy:
c
        rates(9) = (siIa10 * siIe10 * siIn1 + (siIa20 * siIe20 + 
     $               siIa21 * siIe21) * siIn2) * absiI * yn
      endif
c
c (R10) -- CII fine-structure cooling
c
      cIIc10 = cl45 * ynh + cl46 * ynh2 + cl47 * yne
      cIIc01 = cl48 * cIIc10

      if (cIIc10 .eq. 0d0 .or. abcII .le. 1d-5 * abundc) then
        rates(10) = 0d0
      else
c
        cIIn1 = cIIc01 / (cIIc01 + cIIc10 + cIIa10)
c
        rates(10) = cIIa10 * cIIe10 * cIIn1 * abcII * yn
      endif
c
c (R11) -- SiII fine-structure cooling
c
      siIIc10 = cl49 * ynh + cl50 * yne
      siIIc01 = cl51 * siIIc10

      if (siIIc10 .eq. 0d0 .or. absiII .le. 1d-5 * abundsi) then
        rates(11) = 0d0
      else
c
        siIIn1 = siIIc01 / (siIIc10 + siIIa10 + siIIc01)
c
        rates(11) = siIIa10 * siIIe10 * siIIn1 * absiII * yn
      endif
c
c Benchmarking suggests that writing this out explicitly is more efficient
c than using a loop (although this is probably only true if the compiler
c optimization is poor).
c
      ylam = rates(1) + rates(2) + rates(3) + rates(4) + rates(5) + 
     $       rates(6) + rates(7) + rates(8) + rates(9) + rates(10) +
     $       rates(11)
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             C O O L _ F U N C             \\\\\\\\\\
c
c=======================================================================

c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
c
c=======================================================================
c
      subroutine three_level_pops(r01, r02, r12, r10, r20, r21, 
     $                            n0, n1, n2)
      implicit none
#include "cool.h"
      REAL r01, r02, r12, r10, r20, r21
      REAL n0 , n1 , n2
      REAL a1 , a2 , a3 , b1 , b2 , b3
c
c If excitation rates are negligibly small, then we assume that all
c of the atoms are in level 0:
c
      if (r01 .eq. 0d0 .and. r02 .eq. 0d0) then
        n0 = 1d0
        n1 = 0d0
        n2 = 0d0
        return
      endif

      a1 = r01 + r02
      a2 = -r10
      a3 = -r20
      b1 = r01
      b2 = -(r10 + r12)
      b3 = r21
c
      n2 = -a1 * (a1 * b2 - b1 * a2) / ((a1 - a2) * 
     $     (a1 * b3 - b1 * a3) - (a1 - a3) * 
     $     (a1 * b2 - b1 * a2))
c
      n1 = (a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2
c
      n0 = 1d0 - n1 - n2
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
c
c=======================================================================
