      SUBROUTINE hy09dustset(npart)
c
c--Set initial bins and dust densities for Hirashita & Yan (2009)
c     method of dust coagulation.
c
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09rho'
      INCLUDE 'COMMONS/HY09check'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'

      INTEGER npart

      REAL*8 xnorm, xintegral, sum, gasdensity
      REAL*8 bin_numdens, bin_rho(HY09_ndust_bins)
      INTEGER ipart, i
c
c--Set min, max sizes and slope (in cm, or XX*1E-04 for microns)
c
      HY09_size_min = 0.005 * 1.0E-04
      HY09_size_max = 1000.0 * 1.0E-04
c      HY09_size_max = 100.0 * 1.0E-04
      HY09_cutoff   = 0.25 * 1.0E-04
      HY09_slope    = -3.5
      HY09_dust_density = 2.26
      HY09_dustgas_ratio = 0.01
      HY09_vcoag_coeff = 21.4 * (12.0)**(5./6.)/(3.4E+10)**(1./3.)/
     &     SQRT(HY09_dust_density)
c
c--Set up bins (bounds and mid-points)
c
      CALL hy09dustsetbins
c
c--Set up initial bin mass densities
c
      DO i = 1, HY09_ndust_bins
         IF (HY09binsize_bounds(i).LE.HY09_cutoff) THEN
            bin_numdens = (HY09binsize_bounds(i)**(HY09_slope+4.0) - 
     &           HY09binsize_bounds(i-1)**(HY09_slope+4.0)) /
     &           (HY09_slope+4.0)
            bin_rho(i) = bin_numdens
         ELSEIF (HY09binsize_bounds(i).GT.HY09_cutoff .AND.
     &           HY09binsize_bounds(i-1).LE.HY09_cutoff) THEN
            bin_numdens = (HY09_cutoff**(HY09_slope+4.0) -
     &           HY09binsize_bounds(i-1)**(HY09_slope+4.0)) /
     &           (HY09_slope+4.0)
            bin_rho(i) = bin_numdens
         ELSE
            bin_rho(i) = 0.
         ENDIF
c         WRITE (70,*) HY09binsizes(i)*1.0E+04,HY09bin_rho(i)
      END DO
c
c--Normalise.  Following Hirashita & Omukai (2009) = HO09.
c     N(a)*da is the number of grains in interval da, so N(a)
c     has units of number per cm (per unit mass of dust).
c     The integral of 4*pi/3*grain_density*a^3*N(a) da = mtot
c     the total mass of dust.  Now da = a*d(ln a), and rearranging we
c     get: 1/mtot *integral{ a^4 N(a) d(ln a) = 3/(4*pi*grain_density) }
c
c     Note:  I think HO09's Fig.3 is slightly wrong.  The normalisation
c     only makes sense if you use dust density of 2.26, but they have
c     a mixture of 2.26 and 3.3 g/cm^3 dust.
c
      gasdensity = 1.0
c      gasdensity = 1.67E-24*(1.0+4*0.083)*num_dens_H
      xnorm = HY09_dustgas_ratio*gasdensity
      xintegral = (HY09_cutoff**(4.0+HY09_slope) - 
     &     HY09_size_min**(4.0+HY09_slope)) / (4.0+HY09_slope)
      print *,'xnorm, xint ',xnorm, xintegral
      xnorm = xnorm/xintegral
      bin_rho(:) = bin_rho(:)*xnorm
c
c--Check total grain mass (g/cm^3)
c
      sum = 0.
      DO i = 1, HY09_ndust_bins
         sum = sum + bin_rho(i)
      END DO
      PRINT *,'Total mass of grains (per cm^3): ',sum
      PRINT *,'Should be                      : ',
     &     HY09_dustgas_ratio*gasdensity
c
c--Set all gas particles to have the same dust to mass ratio
c     NOTE: Units of HY09bin_rho are cgs/(gas density), in other words
c     to get \tilda{rho} from HY09 you need to multiply by the gas
c     density (see the xnorm and gasdensity definitions above).  
c     If you multiply by the gas density in cgs units, then
c     you will get \tilda{rho} in full cgs units.  This is because we
c     want the dust density to of an SPH particle to follow the gas
c     density (e.g. during compression).
c
      DO ipart = 1, npart
         HY09bin_rho(:,ipart) = bin_rho(:)
      END DO
      HY09frac_h_moved(:) = 0.

      RETURN
      END

c--------------------------------------------------------------

      SUBROUTINE hy09dustsetbins

      IMPLICIT NONE

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/HY09dustprops'

      INTEGER i
      REAL delta

      delta = EXP(LOG(HY09_size_max/HY09_size_min)/HY09_ndust_bins)
      HY09binsize_bounds(0) = HY09_size_min
      HY09binmass_max(0) = 0.
      DO i = 1, HY09_ndust_bins
         HY09binsize_bounds(i) = HY09binsize_bounds(i-1)*delta
         HY09binsizes(i) = 0.5*(HY09binsize_bounds(i-1) + 
     &        HY09binsize_bounds(i))
         HY09binmass(i) = 4.0*pi/3.0*HY09binsizes(i)**3*
     &        HY09_dust_density
         HY09binmass_max(i) = 4.0*pi/3.0*HY09binsize_bounds(i)**3*
     &        HY09_dust_density
      END DO
      PRINT *,'Bin sizes (microns) are: '
      PRINT *,(HY09binsize_bounds(i)*1.0E+04, i=0, HY09_ndust_bins)
      PRINT *,'Bin centres (microns) are: '
      PRINT *,(HY09binsizes(i)*1.0E+04, i=1, HY09_ndust_bins)

      RETURN
      END
