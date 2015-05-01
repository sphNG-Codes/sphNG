      REAL FUNCTION GETCD(rho_g, dvij, cs, ipart, idragscheme)

c************************************************************
c                                                           *
c  Wrapper for various drag_coefficient schemes that can    *
c  be selected using the idragscheme parameter.             *
c                                                           *
c  Ben Ayliffe, 20th March 2012.                            *
c                                                           *
c************************************************************

      IMPLICIT NONE
      
      INCLUDE 'idim'

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/planetesimal'

      REAL rho_g, dvij, cs
      REAL getcd_pmc, getcd_bb, getcd_lp
      REAL r_planetesimal, rho_planetesimal
      INTEGER iregime, ipart, idragscheme

      r_planetesimal = r_planetesimals(1)
      rho_planetesimal = rho_planetesimals(1)

      IF (idragscheme.EQ.0) THEN
c
c--Perret-Murray-Cley
c
         getcd = getcd_pmc(rho_g, dvij, cs, r_planetesimal)
         getcd = getcd*9.*dvij/(8.*r_planetesimal*rho_planetesimal)

      ELSEIF (idragscheme.EQ.1) THEN
c
c--Brasser et al. 2007 + Baines et al. 1968
c
         getcd = getcd_bb(rho_g, dvij, cs, iregime, ipart)
         getcd = getcd*9.*dvij/(8.*r_planetesimal*rho_planetesimal)

      ELSEIF (idragscheme.EQ.2) THEN
c
c--Laibe & Price 2012b.
c
         getcd = getcd_lp(ipart, rho_g, dvij, cs)

      ELSE
         WRITE (iprint,*) 'Invalid idragscheme selected'
         STOP
      ENDIF

      RETURN

      END FUNCTION GETCD


      REAL FUNCTION GETCD_PMC(rho_g, dvij, cs, r_planetesimal)

c************************************************************
c                                                           *
c                     idragscheme = 0                       *
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

      REAL rho_g, dvij, cs
      REAL reynolds, lambda_mfp
      REAL vthermal, r_planetesimal

c      lambda_mfp = 1.619E-15/rho_g
c--Constant using lambda as defined in Perets & Murray-Clay 2011
      lambda_mfp = 6.0E-15/rho_g
      vthermal = sqrt(8./pi)*cs
      reynolds = 2.*r_planetesimal*dvij/(0.5*vthermal*lambda_mfp)
         
      getcd_pmc = (24./reynolds)*(1.+0.27*reynolds)**0.43 + 
     &     0.47*(1.-exp(-0.04*reynolds**0.38))

      RETURN

      END FUNCTION GETCD_PMC


      REAL FUNCTION GETCD_BB(rho_g, dvij, cs, iregime, ipart)

c************************************************************
c                                                           *
c                     idragscheme = 1                       *
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
      REAL dvij, cs, coeff
      REAL getcv, get1overmu
      REAL reynolds, knudsen, mach
      REAL r_planetesimal, rho_planetesimal
      INTEGER iregime, ipart

      r_planetesimal = r_planetesimals(1)
      rho_planetesimal = rho_planetesimals(1)

c--Constant from Stepinski & Valageas 1996 (I).
      iregime = 0               ! Stokes
      IF (1.619E-15/rho_g.GT.(4.*r_planetesimal/9.))
     &     iregime = 1          ! Epstein 

      
      IF (iregime.EQ.0) THEN
c
c--Equation 19, Brasser et al. 2007
c  Coefficient divided by 10 to convert from Brasser's SI units.
c
         knudsen = (1.67E-9*udist**2/umass)/(rho_g*r_planetesimal)
c
c--Equation 14, Brasser et al. 2007
c
         mach = dvij/cs
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

         vratio = sqrt(gamma_l/2.)*dvij/cs
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



      REAL FUNCTION getcd_lp (ipart, rho_g, dvij, cs)

c************************************************************
c                                                           *
c                     idragscheme = 2                       *
c                                                           *
c  Routine to find drag coefficient for a planetesimal in   *
c  gas, based upon Laibe & Price 2012b.                     *
c                                                           *
c  Ben Ayliffe, 20th March 2012.                            *
c                                                           *
c                                                           *
c  Different iregime switch between Stokes drag (0) and     *
c  Epstein drag (1) when compared with Brasser switch       *
c  above.                                                   *
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
      INCLUDE 'COMMONS/logun'

      REAL rho_g, dvij, cs
      INTEGER ipart

      REAL reynolds
      REAL gamma_l, muvisc, muconst, lambda_mfp
      REAL getcv, get1overmu
      REAL r_planetesimal, rho_planetesimal
      INTEGER idrag_structure, iregime

      REAL dragcoeff, volfrac_gas

      r_planetesimal = r_planetesimals(1)
      rho_planetesimal = rho_planetesimals(1)

      dragcoeff = 0.
      volfrac_gas = 1.0
c
c--Calculate value of gamma for RT calculations
c
      IF (encal.EQ.'r' .OR. (encal.EQ.'i' .AND. use_tprof)) THEN
         gamma_l = Rg/uergg/getcv(rho(ipart),vxyzu(4,ipart))*
     &        get1overmu(rho(ipart),vxyzu(4,ipart)) + 1.
      ELSEIF (encal.EQ.'i') THEN
         gamma_l = gamma
      ENDIF

c--H2 collision cross section sigmah2 = 2.367E−15 cm^2
c--H2 mass mh2 = 3.3474472E-24 g
c--muconst = 5*mh2/(64*sigmah2)

      muconst = 3.35818117E-16
      muvisc = muconst*sqrt(pi/gamma_l)*cs
      lambda_mfp = sqrt(pi*gamma_l/2.)*muvisc/(rho_g*cs)

      iregime = 0                                           ! Stokes
      IF (lambda_mfp.GT.(4.*r_planetesimal/9.)) iregime = 1 ! Epstein

      idrag_structure = 0
      IF (iregime.EQ.0) THEN    !--Stokes
         reynolds = 2.*r_planetesimal*rho_g*dvij/muvisc
      
         IF (reynolds.LE.1.) THEN
            idrag_structure = 1
         ELSEIF (reynolds.LE.800.) THEN
            idrag_structure = 2
         ELSE
            idrag_structure = 3
         ENDIF


         IF (idrag_structure.EQ.1) THEN !--linear
            dragcoeff = 9.*6.*muvisc/
     &           (rho_g*4.*r_planetesimal**2*rho_planetesimal)

         ELSEIF (idrag_structure.EQ.2) THEN !--power-law
            dragcoeff = 27.*(muvisc/(2.*rho_g))**0.6/
     &           (volfrac_gas**0.4*r_planetesimal**1.6*
     &           rho_planetesimal)*dvij**0.4
            
         ELSEIF (idrag_structure.EQ.3) THEN !--quadratic
            dragcoeff = 9.*0.22*dvij/
     &           (4.*r_planetesimal*rho_planetesimal*volfrac_gas)
         ELSE
            WRITE (iprint,*) 'getcd'
            WRITE (iprint,*) 'invalid idrag_structure for Stokes'
            STOP
         ENDIF
      ELSEIF (iregime.EQ.1) THEN !--Epstein
c
c--Using Paardekooper & Mellema 2006 scheme that interpolates between
c  asymptotic schemes that apply at low and high mach numbers in the
c  Epstein regime.
c
         dragcoeff = 3.*sqrt(8.*pi/gamma_l)*cs/
     &        (pi*r_planetesimal*rho_planetesimal*volfrac_gas)*
     &        sqrt(1. + (9.*pi/128.)*dvij**2/(cs*cs))

      ELSE
         WRITE (iprint,*) 'Invalid value for iregime in getcd'
         STOP
      ENDIF

      getcd_lp = dragcoeff

      END FUNCTION getcd_lp
