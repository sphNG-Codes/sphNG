      REAL FUNCTION GETCD(rho_g, vr, cs, r_planetesimal)

c************************************************************
c                                                           *
c  Routine to find drag coefficient for a planetesimal in   *
c  gas, based upon Perets & Murray-Clay 2011                *
c                                                           *
c  Ben Ayliffe, 8th July 2011.                              *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'

      REAL rho_g, vr, cs
      REAL reynolds, lambda_mfp
      REAL vthermal, r_planetesimal
      
c      lambda_mfp = 1.619E-15/rho_g
c--Constant using lambda as defined in Perets & Murray-Clay 2011
      lambda_mfp = 6.0E-15/rho_g
      vthermal = sqrt(8./pi)*cs
      reynolds = 2.*r_planetesimal*vr/(0.5*vthermal*lambda_mfp)
         
      getcd = (24./reynolds)*(1.+0.27*reynolds)**0.43 + 
     &     0.47*(1.-exp(-0.04*reynolds**0.38))

      RETURN

      END FUNCTION GETCD
