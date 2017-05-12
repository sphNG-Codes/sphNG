      FUNCTION Qv(v)
c
c--This is the absorption efficiency of the dust with frequency
c
c--This is based on Zucconi et al. (2001)
c
c--NOTE: The values are in cgs units, not code units.
c     The units are cm^2/H_2 and it assumes a standard 100:1 gas to dust
c     ratio.  Therefore, to allow different metallicities, we multiply
c     the result by metallicity which is in solar units.
c
      IMPLICIT NONE

      REAL*8 v,Qv

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rbnd'

      IF (c/v .LT. 1e-3) THEN
         Qv = 3.9E-22 * (1.0E-4/(c/v))**1.4 * metallicity
      ELSEIF (c/v .GT. 0.04) THEN
         Qv = 3.3E-26 * (1.06E-1/(c/v))**2.0 * metallicity
      ELSE
         Qv = 1.5E-24 * (1.4E-2/(c/v))**1.6 * metallicity
      ENDIF

      RETURN
      END

c-----------------------------------------------------------

