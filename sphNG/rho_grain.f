      SUBROUTINE rho_grain(rho_planetesimal)
c************************************************************
c                                                           *
c     This subroutine returns the dust grain intrinsic      *
c     density.                                              *
c     As it is now, it is just a fixed constant but         *
c     in the future we could make this a dynamic quantity.  *
c                                                           *
c     NOTE: rho_planetesimal specified in code units        *
c     Example:                                              *
c     Density of 1.0 g/cm^3 is rho_planetesimal = 1.0/udens *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMONS/units'

      REAL rho_planetesimal

      rho_planetesimal = 1.0/udens

      END SUBROUTINE rho_grain
