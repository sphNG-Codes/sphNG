       SUBROUTINE constan
c************************************************************
c                                                           *
c  This subroutine initializes the mathematical, physical   *
c     and astronomical constants (refs: Handbook of         *
c     Chemistry and Physics, 1975-1976 edition;             *
c     Astrophysical Concepts)                               *
c                                                           *
c************************************************************

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/cgas'
c
c--Mathematical constants
c
      pi = 3.141592653589
c
c--Physical constants (in cgs units)
c
      c = 2.997924d10   
      gg = 6.672041d-8
      Rg = 8.314d7
      radconst = 7.5646d-15
      cgsmu0 = 4.d0*pi
c
c--charge on electron (in esu)
c
      qe = 4.8032068d-10
c
c--ratio of proton to electron mass
c
      ratiomptome = 1836.152755656d0
c
c--collection of constants in front
c  of conductivity term
c
      sigmaterm = 4.d0*pi/(c*c)*qe*qe
c
c--term required in rate coefficient
c  for physical resistivity
c     sqrt(8*k_b/(pi*me))
c
      ratecoeffconst = dsqrt(8.*boltzmannk/(pi*elecm))
c
c--Astronomical constants (in cgs units)
c
c--Solar mass and radius
c
      solarm = 1.991d33 
      solarr = 6.959500d10
      solarl = 3.9d33
c
c--Earth mass and radius
c
      earthm = 5.979d27
      earthr = 6.371315d8
c
c--Distance scale
c
      au = 1.496d13  
      pc = 3.086d18
c
c--Gas molecular weight, mu
c
c      gmw = 2.0
c      gmw = 2.46
      gmw = 4.0/(2*0.7+0.28)

      RETURN
      END
