      FUNCTION etafunc(rhoi,ui,fac,iresist)
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
      INCLUDE 'COMMONS/dissi'
      REAL,    INTENT(IN) :: rhoi
      REAL,    INTENT(IN) :: ui,fac
      INTEGER, INTENT(IN) :: iresist
      REAL etafunc
      
      REAL*8 ratecoefftH,ratecoefftHe,densne,densnei,densnh,vrms
      REAL*8 ratecoefft,Te,sigmaterme,sigmae,sigma,rhoreal,etareal
      REAL*8 X,Y,Z,specific,uoverT,gmwi,gasmw
      REAL*4 rhoreal4
      REAL, EXTERNAL :: getcv
c
c--get temperature in K
c
c      Te = 2./3.*ui/(Rg/gmw/uergg)
      rhoreal4 = rhoi
      Te = ui/getcv(rhoreal4,ui)
c
c--handle T=0 case - prevent floating exceptions         
c
      IF (Te.le.0.) THEN
         etafunc = 0.
         RETURN
      ENDIF
c
c--assumed composition
c
      X = 0.7
      Y = 0.28
      Z = 0.02
      gasmw = 4./(2.*X + Y)
      densnh = X*rhoi/mH
c
c--Firstly, work out the rate coefficient for electron-neutral scattering
c  -this is a function of the temperature and composition.
c  We take the rates from Galli & Pinto (2008), combined with an assumed composition
c  (their results are parametrised between 0 and 10,000K)
c
      IF (Te.lt.10000) THEN
         !vrms = dsqrt(8.d0*kb*Te/(pi*me))*sqrt(Te)
         vrms = ratecoeffconst*sqrt(Te)
         ratecoefftH  = 3.6d-11*(vrms/1.d5)**1.3
         ratecoefftHe = 0.428d-9*sqrt(Te)
      ELSE
c
c--above 10,000K we adopt the Langevin rates
c  this doesn't really matter anyway though, as we are
c  almost certainly fully ionised here anyway
c
         ratecoefftH  = 1.1d-7
         ratecoefftHe = 4.5d-8
      ENDIF
      ratecoefft = X*ratecoefftH + Y*ratecoefftHe
c
c--Secondly, get the electron number density from a model for
c  the ionisation fraction. We do this based on the models of Nakano et al. (2002)
c  as described in Figs 1 and 3 of Wardle and Ng.
c  (really we should do this self-consistently ourselves
c  by solving for the ionisation/dissociation equilibrium from cosmic rays)
c  The ionisation depends most strongly on assumptions about the 
c  dust composition and in particular the grain size.
c
c      densne = 1.d-2
      
      IF (iresist.EQ.3) THEN
c
c--parametrisation of ne in Fig. 1 of Wardle & Ng (1999)
c  (this is assuming a single grain size of 0.1 micron)
c  << gives LOWER resistivity >>
c
         IF (densnh < 2.e8) THEN
            densne = 0.15*(densnh/2.e8)**(0.4) !0.7d-4*ndens**0.4
         ELSEIF (densnh < 3.e10) THEN
            densne = 0.15
         ELSEIF (densnh < 3.e12) THEN
            densne = 0.15*(densnh/3.e10)**(-0.6)
         ELSE
            densne = 0.15*(densnh/3.e10)**(-0.6)*(densnh/3.e12)**0.5
         ENDIF
      ELSE
c
c--Parametrisation of Fig. 3 of Wardle & Ng (1999)
c  (this is for an MRN grain distribution)
c  << this is the default option, gives HIGH resistivity >>
c     
         IF (densnh < 1.e6) THEN
            densne = 2.d-3*(densnh/1.e6)**(0.4)
         ELSEIF (densnh < 6.e6) THEN
            densne = 2.d-3
         ELSEIF (densnh <  1.e9) THEN
            densne = 2.d-3*(densnh/6.e6)**(-0.4)
         ELSE
            densne = 2.d-3*(1.e9/6.e6)**(-0.4)
         ENDIF
      ENDIF
c
c--At high temperatures, the number density of electrons starts
c  to rise due to thermal ionisation. We do this self-consistently
c  by solving the Saha equation
c
      IF (Te.gt.200.) THEN
         call ionisation(rhoi,Te,X,Y,Z,specific,uoverT,gmwi,densnei)
c       print*,' Te = ',Te, ' ne = ',ndense,'mu = ',mu,' ionisation fraction = ',ndense/densnh
         densne = densne + densnei
      ENDIF
c
c--finally, construct the conductivity term from e-n collisions
c  sigmaterm is the term 4.*pi/(c^2)*qe*qe - this is precomputed in constan.f
c     
      sigmaterme = sigmaterm*densne
      sigmae = sigmaterme*(1.d0 + gasmw*ratiomptome)/(ratecoefft)

c     contribution from ion-neutral collisions (negligible)
c      sigmai = qe*qe*ndensi*(1.d0 + 34.*ratiomptome)/(ratecoeffi)

c
c--total conductivity is sum of electron-neutral, ion-neutral and grain-neutral contributions
c  We consider only the electron-neutral contribution
c
      sigma = sigmae !+ sigmai
c
c--finally, get rho in physical units and multiply eta by it
c
      rhoreal = rhoi*umass/udist**3
      etareal = etamhd*rhoreal/sigma
c
c--get eta in code units instead of physical units
c
      etafunc = fac*etareal*utime/udist**2
      
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
      REAL rhocode

      PRINT*,' TESTING RESISTIVITY FUNCTION '
      
      rhologstart = -20.
      deltarholog = 20./dble(npts-1)
      
      OPEN(UNIT=1,FILE='etatest.out',STATUS='REPLACE',FORM='FORMATTED')
      DO i=1,npts
         rho = 10.d0**(rhologstart + (i-1)*deltarholog)
         rhocode = rho*udist**3/umass
         ucode = 5.e-2
         eta = etafunc(rhocode,ucode,1.0,2)
         etareal = eta*udist**2/utime
         WRITE(1,*) rho,eta,etareal
      ENDDO
      CLOSE(UNIT=1)
      
      END SUBROUTINE test_etafunc
