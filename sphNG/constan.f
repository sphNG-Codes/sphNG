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
      pi = 3.141592654
c
c--Physical constants (in cgs units)
c
      c = 2.997924e10   
      gg = 6.672041e-8
      Rg = 8.314e7
      radconst = 7.5646e-15
      cgsmu0 = 4.*pi
c
c--Astronomical constants (in cgs units)
c
c--Solar mass and radius
c
      solarm = 1.991e33 
      solarr = 6.959500e10
      solarl = 3.9e33
c
c--Earth mass and radius
c
      earthm = 5.979e27
      earthr = 6.371315e8
c
c--Distance scale
c
      au = 1.496e13  
      pc = 3.086e18
c
c--Gas molecular weight, mu
c
      gmw = 2.0
c      gmw = 2.46
      IF (encal.EQ.'r') gmw = 4.0/(2*0.7+0.28)

      RETURN
      END
