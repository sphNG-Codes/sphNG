      FUNCTION rho_grain(iphase)
c************************************************************
c                                                           *
c     This function returns the dust grain intrinsic        *
c     density.                                              *
c                                                           *
c     NOTE: rho_planetesimal specified in code units        *
c     Example:                                              *
c     Density of 1.0 g/cm^3 is rho_planetesimal = 1.0/udens *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/planetesimal'

      REAL rho_grain
      INTEGER*1 iphase

c      rho_grain = 1.0/udens
      rho_grain = rho_planetesimals(iphase - 10)

      END FUNCTION rho_grain
