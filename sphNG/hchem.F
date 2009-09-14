c
c Written by S. Glover, AIP, September 2006.
c
c This routine computes the creation term C and 
c destruction term D appearing in the rate equation
c for the number density of atomic hydrogen, n_H, 
c i.e.
c
c  dn_H / dt = C - D n_H
c
c The inputs required are the temperature and density
c of the gas, the total HI column density (which
c controls the charge state of the dust and hence
c the recombination rate on grain surfaces) and the
c fractional abundances of electrons and protons
c (which will not, in general, be equal, as elements
c other than hydrogen may also contribute electrons).
c
      subroutine hchem(temp, yn, NH, abe, abhp, C, D)
      implicit none
#include "cool.h"
      REAL temp, yn, NH, abe, abhp, C, D
 
c Temperature in Kelvin such that kb * temp_1ev = 1eV
      REAL temp_1eV
      parameter (temp_1ev = eV / kboltz)

      REAL lnte, tinv, var1, var2, var3, var4
      REAL AV, G_dust, phi
      REAL yne, ynhp

      REAL k_ci, k_rec, k_gr

c Compute rate coefficients
      lnte = log(temp / temp_1eV)
      tinv = 1d0 / temp
c
c Collisional ionization of HI by electrons
c From A97; based on data from J87
c
      k_ci = dexp(-32.71396786d0 
     $            + 13.5365560d0  * lnte 
     $            - 5.73932875d0  * (lnte**2) 
     $            + 1.56315498d0  * (lnte**3) 
     $            - 0.28770560d0  * (lnte**4) 
     $            + 3.48255977d-2 * (lnte**5)
     $            - 2.63197617d-3 * (lnte**6) 
     $            + 1.11954395d-4 * (lnte**7) 
     $            - 2.03914985d-6 * (lnte**8))
c
c Case B gas-phase recombination
c From F92.
c
      k_rec = 2.753d-14 * (315614d0 * tinv)**1.500d0 / 
     $        ((1d0 + (115188d0 * tinv)**0.407d0)**2.242d0)
c Recombination on grain surfaces. Rate from WD01.
c
      if (abe .lt. 1d-20) then
c We do this to avoid numerical problems at v. small abe
        phi = 1d20
      elseif (iphoto .eq. 0) then
        phi = G0 * dsqrt(temp) / (yn * abe)
      else
        AV     = AV_conversion_factor * dust_to_gas_ratio * NH
        G_dust = G0 * dexp(-2.5d0 * AV)
        phi    = G_dust * dsqrt(temp) / (yn * abe)
      endif

      var1 = 5.087d2 * temp**1.586d-2
      var2 = - 0.4723d0 - 1.102d-5 * log(temp)

      if (phi .eq. 0d0) then
        k_gr = 1.225d-13 * dust_to_gas_ratio
      else
        var3  = 8.074d-6 * phi**1.378d0
        var4  = (1d0 + var1 * phi**var2)
        k_gr     = 1.225d-13 * dust_to_gas_ratio /
     $             (1d0 + var3 * var4)
      endif
c
c Compute creation and destruction rates for hydrogen
c
      yne  = abe  * yn
      ynhp = abhp * yn
c
      C = k_rec * yne * ynhp + k_gr * ynhp * yn
c
      D = cosmic_ray_ion_rate + k_ci * yne
c
      return
      end
c
c REFERENCES:
c
c      J87    -- Janev et al, 1987, 'Elementary Processes in 
c                Hydrogen-Helium Plasmas', Springer
c      F92    -- Ferland et al, 1992, ApJ, 387, 95
c      A97    -- Abel et al, 1997, New Astron, 2, 181 
c      WD01   -- Weingartner & Draine, 2001, ApJ, 563, 842
c
