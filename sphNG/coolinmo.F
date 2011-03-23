c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              C O O L I N M O              \\\\\\\\\\
c
c=======================================================================
c
       subroutine coolinmo
c
c    Based on an original routine by G. Suttner (University Wuerzburg, 1995)
c    OI, CI cooling added by M. D. Smith (Armagh Observatory, 2000-2001)
c    Extensively rewritten by S. Glover (AMNH, 2003-2005)
c    -- fine structure cooling from CII, SiI, SiII added
c    -- fine structure cooling treatment is now exact
c    -- new H2 cooling function used (taken from Le Bourlot et al 1999)
c    -- high T metal-free cooling from Sutherland & Dopita (1993)
c    -- many other minor changes
c
c  PURPOSE:  Tabulate cooling function for various ions, atoms and molecules
c            before starting run
c
c  REFERENCES:
c       FL77  -- Flower & Launay, 1977, J. Phys. B, 10, 3673
c       B81   -- Black, 1981, MNRAS, 197, 553 
c       K86   -- Keenan et al, 1986, MNRAS, 220, 571
c       JBK87 -- Johnson, Burke & Kingston, 1987, J. Phys B, 20, 2553
c       HM89  -- Hollenbach & McKee, 1989, ApJ, 342, 306
c       RLB90 -- Roueff & Le Bourlot, 1990, A&A, 236, 515
c       P90   -- Pequignot, 1990, A&A, 231, 499
c       R90   -- Roueff, 1990, A&A, 234, 567
c       DK91  -- Dufton & Kingston, 1991, MNRAS, 248, 827
c       S91   -- Schroeder et al, 1991, J. Phys. B, 24, 2487
c       C92   -- Cen, 1992, ApJS, 78, 341 
c       SD93  -- Sutherland & Dopita, 1993, ApJS, 88, 253
c       BT94  -- Bakes & Tielens, 1994, ApJ, 427, 822
c       W95   -- Wolfire et al, 1995, ApJ,  443, 152
c       WBV96 -- Warin, Benayoun & Viala, 1996, A&A, 308, 535
c       BBT98 -- Bell, Berrington & Thomas, 1998, MNRAS, 293, L83
c       LPF99 -- Le Bourlot, Pineau des Forets & Flower, 1999, MNRAS, 305, 802
c       WB02  -- Wilson & Bell, 2002, MNRAS, 337, 1027
c       W03   -- Wolfire et al, 2003, ApJ, 587, 278
c       AKD07 -- Abrahamsson, Krems & Dalgarno, 2007, ApJ, 654, 1171
c
      implicit NONE
#include "cool.h"
#include "mol_data.h"
#include "fs_data.h"
c
       integer  i, j, natom, natom2, itemp
c
       parameter(natom = 81)
       REAL coolatom(natom), coolatom_temp(natom)
       parameter(natom2 = 76)
       REAL ca2(natom2), ca2_temp(natom2)
c
       REAL rate0(nmd), rate1(nmd), rate2(nmd), rate3(nmd)
c
       REAL   temp    , temp2   , f       , gg       , hh      
     $      , dtemp   , tinv    , tau     , tsqrt    , opratio 
     $      , fortho  , brem    , fpara   , atomic   , tloge
     $      , h2e20   , h2e31   , h2n2    , h2n3     , h2q02
     $      , h2q13   , phi_pah , tinth   , tinq     , tintq
     $      , tisqt   , tfix    , tfintq  , tfinq    , tfinth
c
c atomic cooling:  table of SD93 for zero radiation field, non-equil., 
c Fe = -0.5 (table 10)
c
      DATA coolatom /-68.,-67.,-66.,-65.,-64.,-63.,-62.,-61.,-60.,-59.,
     *     -58.,-57.,-56.,-55.,-54.,-53.,-52.,-51.,-50.,-49.,-48.,-47.,
     *     -46.,-45.,-44.,-43.,-42.,-41.,-40.,-39.,-38.,-37.,-36.,-35.,
     *     -34.,-33.,-32.,-31.,-30.,-26.67,-23.16,-22.83,-22.53,-22.30,
     *     -22.20,-22.11,-22.00,-21.88,-21.77,-21.69,-21.60,-21.54,
     *     -21.53,-21.58,-21.73,-21.91,-22.00,-22.07,-22.23,-22.32,
     *     -22.35,-22.34,-22.43,-22.62,-22.75,-22.80,-22.83,-22.86,
     *     -22.85,-22.84,-22.85,-22.87,-22.87,-22.86,-22.85,-22.83,
     *     -22.79,-22.75,-22.70,-22.65,-22.60/
c
c atomic cooling #2: from SD93 -- zero radiation field, zero metals,
c collisional ionization equlibrium assumed. Does not include effects
c of HI collisional excitation. [Based on table 6 of SD93]
c
c NB. Data runs from log(T) = 4.3 - 8.0 in increments of 0.05.
c     First value corresponds to T = T_min, but is utterly arbitrary 
c     -- this is here purely to make the spline fitting somewhat 
c     easier, but as the actual entries for T < 10**4.3 K are handled
c     specially below, it should not significantly affect the computed
c     rates. 
c
      DATA ca2 /-60.00, -23.38, -23.17, -23.06, -22.95, -22.89, -22.89, 
     *          -22.93, -22.93, -22.86, -22.69, -22.47, -22.27, -22.18, 
     *          -22.19, -22.27, -22.36, -22.46, -22.56, -22.66, -22.74, 
     *          -22.82, -22.89, -22.95, -23.01, -23.05, -23.09, -23.13, 
     *          -23.16, -23.18, -23.20, -23.22, -23.23, -23.24, -23.25, 
     *          -23.25, -23.25, -23.25, -23.25, -23.25, -23.24, -23.23, 
     *          -23.23, -23.23, -23.21, -23.20, -23.19, -23.18, -23.16, 
     *          -23.15, -23.13, -23.11, -23.10, -23.08, -23.06, -23.05, 
     *          -23.03, -23.01, -22.99, -22.97, -22.95, -22.93, -22.91, 
     *          -22.89, -22.86, -22.84, -22.82, -22.80, -22.78, -22.76, 
     *          -22.74, -22.71, -22.69, -22.67, -22.65, -22.62/
c
c H2 ortho-para ratio (N.B. Only used for OI fine structure cooling)
c
      opratio = 2.4d0
      fortho  = opratio / (1.d0 + opratio)
      fpara   = 1.d0 / (1.d0 + opratio)
c
c establish temperature table     
c
      tmin   = 1.d0
      tmax   = 1.d8
      dtlog  = dlog10(tmax) / (nmd - 1)
c
      do 100 itemp = 1, nmd
        temptab(itemp) = 1.D1**( (itemp - 1) * dtlog )
100   continue
c
c tabulate cooling functions for each temp
c
      do 200 itemp = 1, nmd
         temp  = temptab(itemp)
         temp2 = temp * 1d-2
         tau   = temp * 1d-3 + 1d0
         tinv  = 1.d0 / temp
         tinth = tinv**(1d0/3d0)
         tinq  = tinv**0.25d0
         tintq = tinv**0.75d0
         tsqrt = dsqrt(temp)
         tisqt = 1.d0 / tsqrt
         tloge = log(temp)
c
c cl[1-3]: H2 cooling (LPF99)
c
c   Values in range 100 -- 10000K are handled at end. 
c
c   Rate at T > 10000K is arbitrarily fixed to be the same as at 10000K
c   (this is incorrect, but simple, and in any case Lyman-alpha will 
c    dominate at these temperatures)
c
c   Below 100K, we assume that all cooling comes from the J=2-0, 3-1 
c   transitions in the vibrational ground state.
c
c   All rates assume an ortho:para ratio of 3:1 
c
c (cl1) -- LTE rate
c
         h2e20 = 510.06d0
         h2e31 = 844.94d0
c
         if (temp .lt. 5d0) then
           cltab(1,itemp) = 1d-60
         else if (temp .lt. 1d2) then
           h2n2 = 0.25d0 * 5d0 * dexp(-h2e20 / temp)
           h2n3 = 0.75d0 * (7d0 / 3d0) * dexp(-h2e31 / temp)
           f    = 2.94d-11 * h2e20 * h2n2 * kboltz +
     $            4.76d-10 * h2e31 * h2n3 * kboltz
           cltab(1,itemp) = max(f, 1d-60)
         else if (temp .gt. 1d4) then
           cltab(1,itemp) = 1d1**(-18.253d0)
         endif
c
c (cl2) -- H rate -- 0->2 rate tweaked slightly to ensure that we match up
c                    properly with the tabulated rates at 100K
c
         h2q02 = 1d1**(-9.121d0 - (3.983d0 / tau) - (1.844d-1 / tau**2))
     $           * 5d0 * dexp(-h2e20 / temp)
c
         h2q13 = 1d1**(-9.493d0 - (4.435d0 / tau) + (2.707d-1 / tau**2))
     $           * (7d0 / 3d0) * dexp(-h2e31 / temp)
c
         if (temp .lt. 5d0) then
           cltab(2,itemp) = 1d-60
         else if (temp .lt. 1d2) then
           cltab(2,itemp) = (fpara  * h2q02 * h2e20 + 
     $                       fortho * h2q13 * h2e31) * kboltz
         else if (temp .gt. 1d4) then
           cltab(2,itemp) = 1d1**(-21.943d0)
         endif
c
c (cl3) -- H2 rate
c
         h2q02 = 1d1**(-9.946d0 - (2.688d0 / tau) + (2.020d-1 / tau**2))
     $           * 5d0 * dexp(-h2e20 / temp)
c
         h2q13 = 1d1**(-9.785d0 - (3.024d0 / tau) + (2.930d-1 / tau**2))
     $           * (7d0 / 3d0) * dexp(-h2e31 / temp)
c
         if (temp .lt. 5d0) then
           cltab(3,itemp) = 1d-60
         else if (temp .lt. 1d2) then
           cltab(3,itemp) = (fpara  * h2q02 * h2e20  + 
     $                       fortho * h2q13 * h2e31) * kboltz
         else if (temp .gt. 1d4) then
           cltab(3,itemp) = 1d1**(-22.758d0)
         endif
c 
c (cl4, cl5) -- currently unused
c     
c (cl6) --  the atomic cooling function - this is computed by fitting
c           a cubic spline to the data specified in coolatom above, and
c           so is calculated after the main loop is done
c
c (cl7 -- cl13) -- unused in this version of the code
c
c (cl14) -- gas-grain cooling (HM89, eqn 2.15)
c
         gg              = 1d0 - 0.8 * dexp(-75d0 * tinv)
c
         cltab(14,itemp) = 3.8d-33 * tsqrt * gg * dust_to_gas_ratio
c
c (cl15) --  unused in this version of the code
c
c (cl[16-17]) -- Photoelectric heating
c
c [phi_pah is an adjustable parameter introduced by W03. They quote
c  a value of 0.5, which we also adopt here. Note that the 
c  photoelectric heating rate given by BT94 and W95 corresponds to
c  the case phi_pah = 1.0]
c
         phi_pah = 0.5d0

         cltab(16,itemp) = dsqrt(temp) / phi_pah
c
         cltab(17,itemp) = 3.7e-2 * (temp / 1d4)**0.7d0
c
c Fine structure cooling -- three-level atoms
c
c We use the convention that the ground state is level 0, the first
c excited state is level 1, and the second excited state is level 2.
c
c We tabulate here the collisional de-excitation rates for all relevent
c collision partners. Everything else is taken care of in the cooling
c function itself at run time.
c
c (cl[18-29]): OI fine-structure lines 
c
c Collisional de-excitation rates:
c 
c HI: taken from AKD07 below 1000K, extended to higher temperatures
c     with a simple power-law extrapolation
c
c 1 -> 0
         if (temp .lt. 5d0) then
           tfix = 5d0
           tfintq = 1d0 / tfix**0.75d0 
           cltab(18, itemp) = (5d-11 / 3d0) * exp(4.581
     $                      - 156.118       * tfintq
     $                      + 2679.979      * tfintq**2
     $                      - 78996.962     * tfintq**3
     $                      + 1308323.468   * tfintq**4
     $                      - 13011761.861  * tfintq**5
     $                      + 71010784.971  * tfintq**6
     $                      - 162826621.855 * tfintq**7)
     $                      * exp(oxe10 / (kboltz * tfix))
         elseif (temp .lt. 1d3) then
           cltab(18, itemp) = (5d-11 / 3d0) * exp(4.581
     $                      - 156.118       * tintq
     $                      + 2679.979      * tintq**2
     $                      - 78996.962     * tintq**3
     $                      + 1308323.468   * tintq**4
     $                      - 13011761.861  * tintq**5
     $                      + 71010784.971  * tintq**6
     $                      - 162826621.855 * tintq**7)
     $                      * exp(oxe10 / (kboltz * temp))
         else
           cltab(18, itemp) = 6.81d-11 * temp**0.376d0
         endif
c 2 -> 0
         if (temp .lt. 5d0) then
           tfix = 5d0
           tfintq = 1d0 / tfix**0.75d0
           cltab(19, itemp) = 5d-11 * exp(3.297
     $                      - 168.382       * tfintq
     $                      + 1844.099      * tfintq**2
     $                      - 68362.889     * tfintq**3
     $                      + 1376864.737   * tfintq**4
     $                      - 17964610.169  * tfintq**5
     $                      + 134374927.808 * tfintq**6
     $                      - 430107587.886 * tfintq**7)
     $                      * exp(oxe20 / (kboltz * tfix)) 
         elseif (temp .lt. 1d3) then
           cltab(19, itemp) = 5d-11 * exp(3.297
     $                      - 168.382       * tintq
     $                      + 1844.099      * tintq**2
     $                      - 68362.889     * tintq**3
     $                      + 1376864.737   * tintq**4
     $                      - 17964610.169  * tintq**5
     $                      + 134374927.808 * tintq**6
     $                      - 430107587.886 * tintq**7)
     $                      * exp(oxe20 / (kboltz * temp))
         else
           cltab(19, itemp) = 6.34d-11 * temp**0.36d0
         endif
c
c 2 -> 1:  Low T extrapolation here is necessary because the AKD07 
c          fitting function blows up
c
         if (temp .lt. 5d1) then
           cltab(20, itemp) = 2.62d-12 * temp**0.74d0
         elseif (temp .lt. 1d3) then
           cltab(20, itemp) = 3d-11 * exp(3.437
     $                      + 17.443    * tisqt
     $                      - 618.761   * tisqt**2
     $                      + 3757.156  * tisqt**3
     $                      - 12736.468 * tisqt**4
     $                      + 22785.266 * tisqt**5
     $                      - 22759.228 * tisqt**6
     $                      + 12668.261 * tisqt**7)
     $                      * exp(oxe21 / (kboltz * temp))
         else
           cltab(20,itemp) = 3.61d-10 * temp**0.158d0
         endif
c
c H2 rates supplied by Flower (priv. comm.)
c (NB. If the H2 rate is based on the data from J92, then strictly 
c speaking it is only applicable for T < 1500K; however, the rate 
c doesn't misbehave too badly at higher T).
c
c H2 - ortho and para states must be accounted for separately
c
c 1 -> 0
         f               = fortho * 2.70d-11 * (temp**0.362)
         hh              = fpara  * 3.46d-11 * (temp**0.316)
         cltab(21,itemp) = f + hh
c 2 -> 0
         f               = fortho * 5.49d-11 * (temp**0.317)
         hh              = fpara  * 7.07d-11 * (temp**0.268)
         cltab(22,itemp) = f + hh
c 2 -> 1
         f               = fortho * 2.74d-14 * (temp**1.060)
         hh              = fpara  * 3.33d-15 * (temp**1.360)
         cltab(23,itemp) = f + hh
c
c Electron rate -- from my fits to BBT98.
c
c 1 -> 0
         cltab(24,itemp) = 5.12d-10  * (temp**(-0.075))
c 2 -> 0
         cltab(25,itemp) = 4.863d-10 * (temp**(-0.026))
c 2 -> 1
         cltab(26,itemp) = 1.082d-14 * (temp**(0.926))
c
c Proton rate -- from P90
c
c 1 -> 0
         if (temp .lt. 194) then
           cltab(27,itemp) = 6.38d-11 * (temp**0.40)
         else if (temp .lt. 3686) then
           cltab(27,itemp) = 7.75d-12 * (temp**0.80)
         else
           cltab(27,itemp) = 2.65d-10 * (temp**0.37)
         endif
c 2 -> 0
         if (temp .lt. 511) then
           cltab(28,itemp) = 6.10d-13 * (temp**1.10)
         else if (temp .lt. 7510) then
           cltab(28,itemp) = 2.12d-12 * (temp**0.90)
         else
           cltab(28,itemp) = 4.49d-10 * (temp**0.30)
         endif
c 2 -> 1
         if (temp .lt. 2090) then
           cltab(29,itemp) = 2.029d-11 * (temp**0.56)
         else
           cltab(29,itemp) = 3.434d-10 * (temp**0.19)
         endif
c
c (cl[30-41]): CI fine-structure lines 
c
c Collisional de-excitation rates:  
c
c HI: taken from AKD07 below 1000K, extended to higher temperatures
c     with a simple power-law extrapolation
c
c 1 -> 0
         if (temp .lt. 5d0) then
           tfix = 5d0
           tfinq = 1d0 / tfix**0.25d0
           cltab(30, itemp) =  (1d-11 / 3d0) * dexp(3.6593
     $                      + 56.6023    * tfinq
     $                      - 802.9765   * tfinq**2
     $                      + 5025.1882  * tfinq**3 
     $                      - 17874.4255 * tfinq**4 
     $                      + 38343.6655 * tfinq**5 
     $                      - 49249.4895 * tfinq**6 
     $                      + 34789.3941 * tfinq**7 
     $                      - 10390.9809 * tfinq**8)
     $                      * dexp(cIe10 / (kboltz * tfix))
         elseif (temp .lt. 1d3) then
           cltab(30, itemp) = (1d-11 / 3d0) * dexp(3.6593
     $                      + 56.6023    * tinq
     $                      - 802.9765   * tinq**2
     $                      + 5025.1882  * tinq**3 
     $                      - 17874.4255 * tinq**4 
     $                      + 38343.6655 * tinq**5 
     $                      - 49249.4895 * tinq**6 
     $                      + 34789.3941 * tinq**7 
     $                      - 10390.9809 * tinq**8)
     $                      * dexp(cIe10 / (kboltz * temp))
         else
           cltab(30, itemp) = 2.57d-11 * temp**0.31d0
         endif
c 2 -> 0

         if (temp .lt. 5d0) then
           tfix = 5d0
           tfinth = 1d0 / tfix**(1d0/3d0)
           cltab(31, itemp) = 2.0d-12 * dexp(10.8377 
     $                      - 173.4153    * tfinth 
     $                      + 2024.0272   * tfinth**2 
     $                      - 13391.6549  * tfinth**3
     $                      + 52198.5522  * tfinth**4 
     $                      - 124518.3586 * tfinth**5 
     $                      + 178182.5823 * tfinth**6
     $                      - 140970.6106 * tfinth**7 
     $                      + 47504.5861  * tfinth**8)
     $                      * dexp(cIe20 / (kboltz * tfix))
         elseif (temp .lt. 1d3) then
           cltab(31, itemp) = 2.0d-12 * dexp(10.8377 
     $                      - 173.4153    * tinth 
     $                      + 2024.0272   * tinth**2 
     $                      - 13391.6549  * tinth**3
     $                      + 52198.5522  * tinth**4 
     $                      - 124518.3586 * tinth**5 
     $                      + 178182.5823 * tinth**6
     $                      - 140970.6106 * tinth**7 
     $                      + 47504.5861  * tinth**8)
     $                      * dexp(cIe20 / (kboltz * temp))
         else
           cltab(31, itemp) = 1.69d-11 * temp**0.35d0
         endif
c 2 -> 1
         if (temp .lt. 5d0) then
           tfix = 5d0
           tfinq = 1d0 / tfix**0.25d0
           cltab(32, itemp) = 6.0d-12 * dexp(15.8996 
     $                      - 201.3030   * tfinq
     $                      + 1533.6164  * tfinq**2 
     $                      - 6491.0083  * tfinq**3
     $                      + 15921.9239 * tfinq**4 
     $                      - 22691.1632 * tfinq**5 
     $                      + 17334.7529 * tfinq**6 
     $                      - 5517.9360  * tfinq**7) 
     $                      * dexp(cIe21 / (kboltz * tfix))
         elseif (temp .lt. 1d3) then
           cltab(32, itemp) = 6.0d-12 * dexp(15.8996 
     $                      - 201.3030   * tinq
     $                      + 1533.6164  * tinq**2 
     $                      - 6491.0083  * tinq**3
     $                      + 15921.9239 * tinq**4 
     $                      - 22691.1632 * tinq**5 
     $                      + 17334.7529 * tinq**6 
     $                      - 5517.9360  * tinq**7) 
     $                      * dexp(cIe21 / (kboltz * temp))
         elseif (temp .lt. 1d3) then
           cltab(32, itemp) = 6.0d-12 * dexp(15.8996 
     $                      - 201.3030   * tinq
     $                      + 1533.6164  * tinq**2 
     $                      - 6491.0083  * tinq**3
     $                      + 15921.9239 * tinq**4 
     $                      - 22691.1632 * tinq**5 
     $                      + 17334.7529 * tinq**6 
     $                      - 5517.9360  * tinq**7) 
     $                      * dexp(cIe21 / (kboltz * temp))
         else
           cltab(32,itemp) = 4.95d-11 * temp**0.35
         endif
c
c H2 -- ortho and para states must be accounted for separately
c    -- rates from S91 using the fit from WBV96
c  
c 1 -> 0
         f  = (8.7d-11 - 6.6d-11 * dexp(-temp / 218.3d0) 
     $         + 6.6d-11 * exp(-2d0 * temp / 218.3d0) ) * fortho
c
         hh = (7.9d-11 - 8.7d-11 * dexp(-temp / 126.4d0) 
     $         + 1.3d-10 * exp(-2d0 * temp / 126.4d0) ) * fpara
c
         cltab(33,itemp) = f + hh
c 2 -> 0
         f  = (1.2d-10 - 6.1d-11 * dexp(-temp / 387.3d0)) * fortho
c
         hh = (1.1d-10 - 8.6d-11 * dexp(-temp / 223.0d0) 
     $         + 8.7d-11 * exp(-2d0 * temp / 223.0d0) ) * fpara
c
         cltab(34,itemp) = f + hh
c 2 -> 1
         f  = (2.9d-10 - 1.9d-10 * dexp(-temp / 348.9d0)) * fortho
c
         hh = (2.7d-10 - 2.6d-10 * dexp(-temp / 250.7d0) 
     $        + 1.8d-10 * exp(-2d0 * temp / 250.7d0)) * fpara
c
         cltab(35,itemp) = f + hh
c
c Electrons -- from JBK87. Note that the 'electron' rate given in HM89 is 
c              actually the proton rate, as footnote f of table 8 makes 
c              clear. Note that fits are only valid for T < 10^4 K --
c              at higher temperatures, we assume that the collision 
c              strength is constant (and so the rate scales as T^-1/2)
c 1 -> 0
         if (temp .lt. 1d3) then
             f = -9.25141d0 - 0.773782d0 * tloge + 0.361184d0 * tloge**2
     $           - 0.150892d-1 * tloge**3 - 0.656325d-3 * tloge**4
         elseif (temp .le. 1d4) then
             f = 4.44600d2 - 2.27913d2 * tloge + 42.5952 * tloge**2 
     $           - 3.47620d0 * tloge**3 + 0.105085d0 * tloge**4
         else
             f = -0.990634367510893d0
         endif
c
         cltab(36,itemp) = (1d0 / 3d0) * 8.629d-6 * dexp(f) / 
     $                     dsqrt(temp)
c
c 2 -> 0  
c
         if (temp .lt. 1d3) then
             f = -7.69735d0 - 1.30743d0 * tloge + 0.697638d0 * tloge**2 
     $           - 0.111338 * tloge**3 + 0.705277d-2 * tloge**4
         elseif (temp .le. 1d4) then
             f = 3.50609d2 - 1.87474d2 * tloge + 3.61803d1 * tloge**2 
     $           - 3.03283d0 * tloge**3 + 0.938138d-1 * tloge**4
         else
             f = -1.40040241654697d0
         endif
c
         cltab(37,itemp) = 0.2d0 * 8.629d-6 * dexp(f) / dsqrt(temp)
c
c 2 -> 1
         if (temp .lt. 1d3) then
             f = -7.4387d0 - 0.57443d0 * tloge + 0.358264d0 * 
     $           tloge**2  - 0.418166d-1 * tloge**3 + 0.235272d-2 * 
     $           tloge**4
         elseif (temp .le. 1d4) then
             f = 3.86186d2 - 2.02192d2 * tloge + 3.85049d1 * tloge**2 
     $           - 3.19268d0 * tloge**3 + 0.978573d-1 * tloge**4
         else
             f = 0.0198312027880547d0
         endif
c
         cltab(38,itemp) = 0.2d0 * 8.629d-6 * dexp(f) / dsqrt(temp)
c
c Protons -- from RLB90, fits by WBV96 (T < 5000K), SCOG (T > 5000K).
c Assume rate becomes constant above limit of tabulated data.
c
c 1 -> 0
         if (temp .lt. 5d3) then
             f = (9.6d-11 - 1.8d-14 * temp + 1.9d-18 * temp**2)
     $           * temp**0.45
         elseif (temp .lt. 2d4) then
             f = 8.873d-10 * temp**0.117
         else
             f = 2.8268d-9
         endif
c
         cltab(39,itemp) = f
c
c 2 -> 0
         if (temp .lt. 5d3) then
             f = (3.1d-12 - 6.0d-16 * temp + 3.9d-20 * temp**2)
     $           * temp
         elseif (temp .lt. 2d4) then
             f = 2.314d-9 * temp**0.0965
         else
             f = 6.0175e-09
         endif
c
         cltab(40,itemp) = f
c
c 2 -> 1
         if (temp .lt. 5d3) then
           f = (1.0d-10 - 2.2d-14 * temp + 1.7d-18 * temp**2) * 
     $         temp**0.70
         elseif (temp .lt. 2d4) then
           f = 9.198d-9 * temp**0.0535
         else
           f = 1.5624e-08
         endif
c
         cltab(41,itemp) = f
c
c (cl[42-43]) -- SiI fine structure lines
c
c Collisional rates -- HI, protons from HM89
c
c HI:
c 1 -> 0
         cltab(42,itemp) = 3.5d-10 * temp2**(-0.03)
c 2 -> 0
         cltab(43,itemp) = 1.7d-11 * temp2**0.17
c 2 -> 1
         cltab(44,itemp) = 5.0d-10 * temp2**0.17
c
c Proton rates are independent of T, so we don't bother to tabulate
c them here -- instead, they're listed in cool_func.F
c
c No data for H2, electrons.
c
c Fine structure cooling -- two-level atoms
c
c We use the convention that the ground state is level 0 and 
c the excited state is level 1.
c
c The cooling rate is simply:
c
c Lambda = A(1->0) * E(1->0) * f(1,LTE) * C(1->0) / ( A(1->0) + 
c          C(1->0) * (1 + f(1,LTE)) )
c
c where f(1,LTE) is the level population of level 1 in LTE.
c We tabulate f(1,LTE) and the temperature-dependent bits of C(1->0),
c and compute everything else in the cooling function.
c
c (cl[45-48]) -- CII fine structure lines
c
c Collisional rates:
c
c (cl45) CII - HI (HM89 below 2000K; K86 above 2000K)
c
c (Note that the high T coefficient has been tweaked slightly to ensure that
c  the rate is continuous at 2000K -- the adjustment is well within the 
c  errors).
c
         if (temp .le. 2d3) then
            cltab(45,itemp) = 8d-10 * temp2**0.07d0
         else
            cltab(45,itemp) = 3.113619d-10 * temp2**0.385d0
         endif
c
c (cl46) CII - H2 (Below 250K, we use the fit from WBV96 to the data from
c                  FL77. Above 250K, we assume that the rate scales the
c                  same as the low-temp. HI rate)
c
         if (temp .lt. 250d0) then
             f  = (4.7d-10 + 4.6d-13 * temp) * fortho
             hh = (2.5d-10 * temp**0.12) * fpara
         else
             f  = (5.85d-10 * temp**0.07) * fortho
             hh = (4.85d-10 * temp**0.07) * fpara
         endif
c
         cltab(46,itemp) = f + hh
c
c (cl47) CII - electron (WB02). 
c
c (Note that the high T coefficient has been tweaked slightly to ensure that
c  the rate is continuous at 2000K -- the adjustment is well within the 
c  errors).
c
         if (temp .le. 2d3) then
            cltab(47,itemp) = 3.86d-7 / dsqrt(temp2)
         else
            cltab(47,itemp) = 2.426206d-7 / temp2**0.345d0
         endif
c
c Proton rate negligible below 10^4K, ignorable below 10^5K.
c
c Finally, cl48 holds f(1,LTE)
c
         cltab(48,itemp) = 2d0 * dexp(-91.25d0 / temp)
c
c (cl[49-51]) -- SiII fine structure lines
c
c Collisional rates:
c
c (cl49) SiII - HI (fit by SCOG to data from R90).
c
         cltab(49,itemp) = 4.95d-10 * temp2**0.24
c
c (cl50) SiII - electron (DK91 -- extrapolated to T < 4000K, 
c                         assuming constant collision strength)
c
         cltab(50,itemp) = 1.2d-6 / dsqrt(temp2)
c
c cl51 holds f(1,LTE)
c
         cltab(51, itemp) = 2d0 * dexp(-412.24d0 / temp)
c
c cl52 -- HI excitation cooling (aka Lyman-alpha cooling) - for
c         iflag_atom = 1, this is included in the atomic cooling 
c         term computed below
c
c cl53, cl54 -- these values are used in the calculation of the
c cooling rate due to electron recombination with PAHs
c
         cltab(53, itemp) = 0.74d0 / temp**0.068
         cltab(54, itemp) = 4.65d-30 * phi_pah * temp**0.94d0
c
 200  continue
c
c  Rates requiring spline fits to tabulated data:
c
c H2 rovibrational cooling -- data in mol_data.h
c
      do I = 1, nh2data
        h2_temp(i) = 1d1**(2d0 + 5d-2 * (I - 1))
      enddo
c
      call spline_eval(nh2data, h2_temp, h2_lte, nmd, temptab, rate1)
      call spline_eval(nh2data, h2_temp, h2_h_rate, nmd, temptab, rate2)
      call spline_eval(nh2data, h2_temp, h2_h2_rate, nmd, temptab,rate3)
c
      do itemp = 1, nmd
        if (temptab(itemp) .ge. 1d2 .and. temptab(itemp) .le. 1d4) then 
          cltab(1,itemp) = 1d1**rate1(itemp)
          cltab(2,itemp) = 1d1**rate2(itemp)
          cltab(3,itemp) = 1d1**rate3(itemp)
        endif
      enddo
c
c (cl6) --  the atomic cooling function
c
c Calculate the temperatures corresponding to the data in coolatom
c
      if (iflag_atom .eq. 1) then
        do 210 i = 1, natom
          coolatom_temp(i) = 10**(8d0 * (i-1) / 79d0)
210     continue
c
c Compute and evaluate a spline fit to the data in coolatom
c
        call spline_eval(natom, coolatom_temp, coolatom, nmd, temptab,
     $                   rate0)

        do 220 itemp = 1, nmd
          temp   = temptab(itemp)
          atomic = 1d1**(rate0(itemp))
c
c For temp > 10^4 K, include thermal bremsstrahlung component
c
          if (temp .gt. 1d4)  then
            brem = 1.42d-27 * sqrt(temp)
          else
            brem = 0d0
          endif

          cltab(6,itemp)  = atomic + brem
c
c HI excitation cooling is included in this rate, so we set its entry to zero
          cltab(52,itemp) = 0d0
220     continue
c
      elseif (iflag_atom .eq. 2) then
c
        ca2_temp(1) = tmin
        do 230 i = 2, natom2
          ca2_temp(i) = 10**(4.25d0 + 0.05d0 * i)
230     continue
c
        call spline_eval(natom2, ca2_temp, ca2, nmd, temptab, rate0)
c
        do 240 itemp = 1, nmd
          temp   = temptab(itemp)
          if (temp .lt. 10**4.3d0) then
            cltab(6,itemp) = 0d0
          else
            cltab(6,itemp) = 1d1**(rate0(itemp)) 
          endif
c
c Tabulate HI excitation cooling rate separately. 
c This expression is from C92, based on B81.
c
          if (temp .lt. 1d3) then
            cltab(52,itemp) = 0d0
          else
            cltab(52,itemp) = 7.5d-19 * dexp(-1.18348d5 / temp) /
     $                        (1d0 + dsqrt(temp / 1d5))
          endif
c
240     continue
      endif
c
c Approximate temperature derivatives
c
      do 300 itemp = 1,nmd-1
         dtemp           = 1d0 / ( temptab(itemp + 1) - temptab(itemp) )
c
         do 250 j = 1, ncltab 
            dtcltab(j,itemp) = ( cltab(j, itemp+1)
     $                         - cltab(j, itemp  ) ) * dtemp
c
 250     continue
 300  continue
c
c Manually add on the final tabulated value at itemp = ntabtemp
c
      do 310 j = 1, ncltab
         dtcltab(j,nmd) = dtcltab(j,nmd-1)
 310  continue
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             C O O L I N M O               \\\\\\\\\\
c
c=======================================================================
