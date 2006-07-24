      SUBROUTINE energ(ipart, ti, vxyzu, dvxyzu, xyzmh)
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

      DIMENSION vxyzu(4,idim), dvxyzu(4,idim), xyzmh(5,idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/ptmass'

      REAL omega1
c
c--Initialisation
c
      cnormk05 = cnormk*0.5
      gamma1 = gamma - 1.
c
c--pdv term is stored in dvxyzu(4,ipart) in forcei to use memory efficiently
c
      pdv = dvxyzu(4,ipart)
      dvxyzu(4,ipart) = 0.0
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

            dvxyzu(4,ipart) = 1.5/critlambda*(1.0 - rho(ipart)/rhocrit)
            IF (vxyzu(4,ipart).GE.1.5 .AND. dvxyzu(4,ipart).GT.0.0) 
     &           dvxyzu(4,ipart) = 0.0
            IF (vxyzu(4,ipart).LE.1.5/rhoratio .AND. 
     &           dvxyzu(4,ipart).LT.0.) dvxyzu(4,ipart) = 0.0
         ELSEIF (iener.EQ.1) THEN

c
c--Compute change in specific internal energy
c
c  a) pdv term first
c
            IF (iexpan.EQ.0) THEN
               dvxyzu(4,ipart) = cnormk*pdv*pr(ipart)/
     &              rho(ipart)**2
            ELSE
               dvxyzu(4,ipart) = cnormk*pdv*pr(ipart)/
     &              rho(ipart)**2 - vxyzu(4,ipart)*dqexp
            ENDIF
c
c  b) Shock dissipation if appropriate
c
            IF (ichoc.NE.0) THEN
               dvxyzu(4,ipart) = dvxyzu(4,ipart) + cnormk05*dq(ipart)
            ENDIF
c
c  c) Cooling for adiabatic-with-cooling eos
c
            IF(encal.EQ.'c') THEN
               r = SQRT(xyzmh(1,ipart)**2 + xyzmh(2,ipart)**2 + 
     &               xyzmh(3,ipart)**2 )

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
          dvxyzu(4,ipart) = vxyzu(4,ipart)*gamma1*dq(ipart)*rho(ipart)* 
     &                 cnormk05/pr(ipart) - pr(ipart)*dqexp
      ENDIF

      RETURN
      END
