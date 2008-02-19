      FUNCTION etafunc(rhoi,ui)
c************************************************************
c                                                           *
c  This function specifies the (physical) resistivity as a  *
c  function of density and/or temperature                   *
c                                                           *
c  takes in density in code units and returns resistivity   *
c  in code units                                            *
c                                                           *
c************************************************************
      IMPLICIT NONE ! music to my ears
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/units'
      REAL*4 rhoi
      REAL ui
      REAL etafunc
      
      REAL*8 ratecoeff,densne
      REAL*8 ratecoefft,Te,sigmaterme,sigmae,sigma,rhoreal,etareal
      REAL*8 fac
c
c--formula from Wardle & Ng (1999) 
c  ratecoefft = 1.d-15*dsqrt(128.d0*kb/(9.d0*pi*me))*dsqrt(Te)
c  parameter is the pre-calculated non-temperature dependent term
c
      PARAMETER (ratecoeff = 8.2833728D-10)
c
c--Number density of electrons (cm^-3)
c  this is taken from Figure 1 of Wardle & Ng (1999)
c  which is an adaption of Figure 2 in Umebayashi & Nakano (1990)
c
      PARAMETER (densne = 1.d-2)
c
c--get temperature in K
c
      Te = 2./3.*ui/(Rg/gmw/uergg)
c
c--introduce rapid shutoff for temperatures above 2000K
c
      IF (Te.gt.2000) THEN
         fac = exp(-(Te-2000.)/10.)
      ELSEIF (Te.le.0.) THEN
c--handle T=0 case - prevent floating exceptions         
         etafunc = 0.
         RETURN
      ELSE
         fac = 1.
      ENDIF
c
c--calculate rate coefficient for electron-neutral scattering
c
c   formula from Wardle & Ng (1999) 
c   ratecoefft = 1.d-15*dsqrt(128.d0*kb/(9.d0*pi*me))*dsqrt(Te)
c   where we have pre-calculated the non-temperature dependent term
c
      ratecoefft = ratecoeff*dsqrt(Te)  ! units cm^3/s
c
c     conductivity sigma (multiplied by mu0 and needs to be divided by rho)
      sigmaterme = sigmaterm*densne
      sigmae = sigmaterme*(1.d0 + gmw*ratiomptome)/ratecoefft

c     contribution from ion-neutral collisions (negligible)
c      sigmai = qe*qe*ndensi*(1.d0 + 34.*ratiomptome)/(ratecoeffi)

c
c--total conductivity is sum of electron-neutral ion-neutral and grain contributions
c
      sigma = sigmae !+ sigmai
      
      rhoreal = rhoi*umass/udist**3
      etareal = rhoreal/sigma
c
c--multiply by rho in code units and sqrt(Te) in K
c
      etafunc = etareal*utime/udist**2
      
    !print*,' rho = ',rhoi,rhoreal,' ui = ',ui,' T = ',Te,' gmw = ',gmw
      !print*,' eta (physical) = ',etareal,' eta (code) = ',etafunc
      
      END FUNCTION etafunc


      SUBROUTINE test_etafunc
      IMPLICIT NONE
      INCLUDE 'COMMONS/units'
      INTEGER npts,i
      PARAMETER (npts = 100)
      REAL rho,rhologstart,deltarholog,eta,etareal,ucode
      REAL etafunc
      REAL*4 rhocode

      PRINT*,' TESTING RESISTIVITY FUNCTION '
      
      rhologstart = -20.
      deltarholog = 20./dble(npts-1)
      
      OPEN(UNIT=1,FILE='etatest.out',STATUS='REPLACE',FORM='FORMATTED')
      DO i=1,npts
         rho = 10.d0**(rhologstart + (i-1)*deltarholog)
         rhocode = rho*udist**3/umass
         ucode = 5.e-2
         eta = etafunc(rhocode,ucode)
         etareal = eta*udist**2/utime
         WRITE(1,*) rho,eta,etareal
      ENDDO
      CLOSE(UNIT=1)
      
      END SUBROUTINE test_etafunc
