      REAL FUNCTION GETCD(rho_g, vr, cs)

c************************************************************
c                                                           *
c  Routine to find drag coefficient for a planetesimal in   *
c  gas, based upon Brasser et al. 2007.                     *
c                                                           *
c  Ben Ayliffe, 12th April 2011.                            *
c                                                           *
c************************************************************

      IMPLICIT NONE

      REAL rho_g
      REAL r_p, vr, cs
      REAL reynolds, knudsen, mach

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/planetesimal'
c
c--Equation 19, Brasser et al. 2007
c  Coefficient divided by 10 to convert from Brasser's SI units.
c
      knudsen = (1.67E-9*udist**2/umass)/(rho_g*r_planetesimal)
c
c--Equation 14, Brasser et al. 2007
c
      mach = vr/cs
c
c--Equation 18, Brasser et al. 2007
c  Using gamma = 5/3, rather than the 7/5 assumed by Brasser et al.
c
      reynolds = 4.85*mach/knudsen
c
c--Set CD depending upon regime defined by above quantities.
c
      IF (mach.GE.1.0) THEN
         getcd = 2.0
      ELSEIF (mach.LT.1.0 .AND. reynolds.GE.1.0E3) THEN
         getcd = 0.44 + 1.56*mach**2
      ELSEIF (reynolds.LT.1.0E3) THEN
         getcd = 2.*mach**2 + (1. - mach**2)*(24.*(1. +
     &        0.15*reynolds**0.687)/reynolds)
      ELSEIF (knudsen.GE.1.0) THEN
         getcd = 0.0
c--If knudsen > 1, Brasser's formulation gives rho_g = 0.
c  This means the gas drag should be zero, hence CD = 0.
      ENDIF

      RETURN
      END FUNCTION GETCD
