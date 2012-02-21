      REAL FUNCTION GETCD(rho_g, vr, cs, r_planetesimal, iregime, ipart,
     &     dragscheme)

      IMPLICIT NONE
      
      INCLUDE 'COMMONS/logun'

      REAL rho_g, vr, cs, r_planetesimal
      REAL getcd_pmc, getcd_bb
      INTEGER iregime, ipart, dragscheme

      IF (dragscheme.EQ.0) THEN
         getcd = getcd_pmc(rho_g, vr, cs, r_planetesimal)
      ELSEIF (dragscheme.EQ.1) THEN
         getcd = getcd_bb(rho_g, vr, cs, iregime, ipart)
      ELSE
         WRITE (iprint,*) 'Invalid dragscheme selected'
         STOP
      ENDIF

      RETURN

      END FUNCTION GETCD


      REAL FUNCTION GETCD_PMC(rho_g, vr, cs, r_planetesimal)

c************************************************************
c                                                           *
c                     dragscheme = 0                        *
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
         
      getcd_pmc = (24./reynolds)*(1.+0.27*reynolds)**0.43 + 
     &     0.47*(1.-exp(-0.04*reynolds**0.38))

      RETURN

      END FUNCTION GETCD_PMC


      REAL FUNCTION GETCD_BB(rho_g, vr, cs, iregime, ipart)

c************************************************************
c                                                           *
c                     dragscheme = 1                        *
c                                                           *
c  Routine to find drag coefficient for a planetesimal in   *
c  gas, based upon Brasser et al. 2007.                     *
c                                                           *
c  Ben Ayliffe, 12th April 2011.                            *
c                                                           *
c                                                           *
c  Added iregime switch to calculate coefficient for  both  *
c  Stokes drag (0) and Epstein drag (1).                    *
c                                                           *
c  Epstein drag stuff from Baines et al. 1965, Eq. 4.5      *
c                                                           *
c  Regime switch in forcei. Epstein if gas mean free path   *
c  lambda_mfp.GT.4.*r_planetesimal/9.                       *
c                                                           *
c  Ben Ayliffe, 8th June 2011                               *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'

      REAL rho_g, vratio, vratio2, gamma_l
      REAL r_p, vr, cs, coeff
      REAL getcv, get1overmu
      REAL reynolds, knudsen, mach
      INTEGER iregime, ipart

      IF (iregime.EQ.0) THEN
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
            getcd_bb = 2.0
         ELSEIF (mach.LT.1.0 .AND. reynolds.GE.1.0E3) THEN
            getcd_bb = 0.44 + 1.56*mach**2
         ELSEIF (reynolds.LT.1.0E3) THEN
            getcd_bb = 2.*mach**2 + (1. - mach**2)*(24.*(1. +
     &           0.15*reynolds**0.687)/reynolds)
         ELSEIF (knudsen.GE.1.0) THEN
            getcd_bb = 0.0
c--If knudsen > 1, Brasser's formulation gives rho_g = 0.
c  This means the gas drag should be zero, hence CD = 0.
         ENDIF

      ELSEIF (iregime.EQ.1) THEN
         IF (encal.EQ.'r' .OR. (encal.EQ.'i' .AND. use_tprof)) THEN
            gamma_l = Rg/uergg/getcv(rho(ipart),vxyzu(4,ipart))*
     &        get1overmu(rho(ipart),vxyzu(4,ipart)) + 1.
         ELSEIF (encal.EQ.'i') THEN
            gamma_l = gamma
         ENDIF

         vratio = sqrt(gamma_l/2.)*vr/cs
         vratio2 = vratio**2
c--coeff = 1./(2*sqrt(pi))
         coeff = 2.820948E-01

c--Factor of 8 introduced to match transition to Stokes I think. Not sure.
         getcd_bb = 8.*coeff*( (1./vratio + 1./(2.*vratio2*vratio))*
     &        exp(-vratio2) + (1. + 1./vratio2 - 1./
     &        (4.*vratio2*vratio2))*erf(vratio)*sqrt(pi))
      ENDIF

      RETURN
      END FUNCTION GETCD_BB
