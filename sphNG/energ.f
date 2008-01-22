      SUBROUTINE energ(ipart, npart, ntot, ti, vxyzu, dvxyzu, xyzmh,
     &     trho)
c************************************************************
c                                                           *
c  This routine computes the change in internal energy or   *
c     in entropy.                                           *
c                                                           *
c    varsta = entropy : compute change in specific entropy  *
c    varsta = intener : compute change in specific internal *
c                          energy                           *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION vxyzu(4,idim), dvxyzu(4,idim), xyzmh(5,mmax)
      REAL*4 trho(idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/tcooling'

      REAL omega1
c
c--For MPI code do not need to add energy changes that are added only added
c     by the process that actually holds the particle
c
      IF (ipart.GT.ntot) THEN
         xfac = 0.0
         ipartntot = ipart + ntot + 2
      ELSE
         xfac = 1.0
         ipartntot = ipart
      ENDIF
c
c--Initialisation
c
      gamma1 = gamma - 1.
c
c--pdv term is stored in dvxyzu(4,ipart) in forcei to use memory efficiently
c
      pdv = dvxyzu(4,ipart)
      dvxyzu(4,ipart) = 0.0

      IF (iener.EQ.3) THEN
c
c--Parameters for cooling curve (proton mass and heating rate) 
c     initiated in thermeq
c     Heating rate from Vazquez-Semadeni 2006 is divided by hydrogen mass
c     as required for energy equation (Koyama & Inutsuka 2002).
c     Koyama & Inutsuka 2002 use d(rho*u)/dt rather than du/dt which
c     is used here for SPH. Everything is done in g and cm for now.
c
c--Find temperature in K

         tempiso = 2./3.*vxyzu(4,ipart)/(Rg/gmw/uergg)
c
c--Cooling rate from Vazquez-Semadeni
c
         coolingr=(1.e7*exp(-1.148e5/(tempiso+1000.))+
     &        0.014*sqrt(tempiso)*exp(-92./tempiso))*heatingr
c
c--Find change in energy - in the energy equation (Koyama & Inutsuka 2002)
c     the cooling term is multiplied by rho/mass of hydrogen

         dvxyzu(4,ipart)=heatingr-trho(ipart)*udens*coolingr/mpenerg
c
c--Convert from cm^2 s^-3 to code units
c
         dvxyzu(4,ipart)=dvxyzu(4,ipart)*utime**3./udist**2.
      END IF
c
c--End of section for cooling curve
c

c
c--dq from expansion
c
      chi = 4.0 - 3.0*gamma
      CALL scaling(ti, rscale, drdt, dlnrdt)
      dqexp = chi*dlnrdt

      IF (varsta.NE.'entropy') THEN
         IF (iener.EQ.2) THEN
c
c--Cooling instability
c
c            critlambda = 0.1
c            rhoratio = 100.0
c            rho_zero = 0.001
c            rhocrit = 3.0*rho_zero

            critlambda = 0.1
            rhoratio = 10.0
            rho_zero = 0.001
            rhocrit = 3.0*rho_zero

            dvxyzu(4,ipart) = 1.5/critlambda*(1.0 - trho(ipart)/rhocrit)
            IF (vxyzu(4,ipart).GE.1.5 .AND. dvxyzu(4,ipart).GT.0.0) 
     &           dvxyzu(4,ipart) = 0.0
            IF (vxyzu(4,ipart).LE.1.5/rhoratio .AND. 
     &           dvxyzu(4,ipart).LT.0.) dvxyzu(4,ipart) = 0.0
         ELSEIF ((iener.EQ.1).OR.(iener.EQ.3)) THEN

c
c--Compute change in specific internal energy
c
c  a) pdv term first
c
            IF (iexpan.EQ.0) THEN
               dvxyzu(4,ipart) = dvxyzu(4,ipart)+pdv*pr(ipart)/
     &              trho(ipart)**2

            ELSE
               dvxyzu(4,ipart) = dvxyzu(4,ipart)+pdv*pr(ipart)/
     &              trho(ipart)**2 - vxyzu(4,ipart)*dqexp*xfac
            ENDIF
c
c  b) Shock dissipation if appropriate
c
            IF (ichoc.NE.0) THEN
               dvxyzu(4,ipart) = dvxyzu(4,ipart) + 0.5*dq(ipart)
            ENDIF

c
c  c) Cooling for adiabatic-with-cooling eos
c
            IF(encal.EQ.'c' .AND. xfac.EQ.1.0) THEN
               r = SQRT(xyzmh(1,ipartntot)**2 + xyzmh(2,ipartntot)**2 + 
     &               xyzmh(3,ipartntot)**2 )

c Omega1 is 1/omega
c !! I think listpm(1) should be the array index of the 
c !! central object. If not, this is wrong.
               omega = sqrt(xyzmh(4,listpm(1))/ r**3)
               omega1 = 1./omega

c Cooling times used by Rice et al 2003 are 5/omega and 3/omega
c              tcool = 2.0 * omega1
c              tcool = 3.0 * omega1
               tcool = 5.0 * omega1

c               PRINT *, dvxyzu(4,ipart),vxyzu(4,ipart)/tcool,omega
               dvxyzu(4,ipart) = dvxyzu(4,ipart) - vxyzu(4,ipart) 
     &              / tcool
            ENDIF
         ENDIF
c
c--Compute change in specific entropy
c
c  a) Shock dissipation
c
      ELSEIF (ichoc.NE.0) THEN
          dvxyzu(4,ipart) = vxyzu(4,ipart)*gamma1*dq(ipart)*
     &        trho(ipart)*0.5/pr(ipart) - pr(ipart)*dqexp*xfac
      ENDIF

      RETURN
      END
