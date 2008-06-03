      SUBROUTINE eospg(ipart,vxyzu,rho,pr,vsound,ekcle)
c************************************************************
c                                                           *
c  This routine computes the pressure and sound speed       *
c     according to a perfect gas equation of state          *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      REAL*4 rho,pr,vsound
      DIMENSION vxyzu(4,idim), rho(idim), pr(idim), vsound(idim)
      DIMENSION ekcle(5,iradtrans)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/vargam'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/physeos'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/bodys'

      CHARACTER*7 where

      DATA where/'eospg'/

      gama1 = gamma - 1.

      IF (vxyzu(4,ipart).LE.0) THEN
         WRITE (*,*) 'ERROR - vxyzu(4,ipart).LE.0 ',vxyzu(4,ipart),ipart
         CALL quit
      ENDIF
c
c--Variable is internal energy
c
c--Variable gamma equation of state
c
      IF (encal.EQ.'v') THEN
          gam1 = gam - 1.
          gamdh1 = gamdh - 1.
          gamah1 = gamah - 1.
          IF (iorig(ipart).GT.n1 .AND. iorig(ipart).LE.n1+n2) THEN
             uzero = uzero_n2
          ELSE
             uzero = RK2 * (rhozero ** gama1)
          ENDIF
          RK3 = uzero / (rhocrit ** gam1)
          RK4 = RK3 * (rhocrit2 ** gam1) / (rhocrit2 ** gamdh1)
          RK5 = RK4 * (rhocrit3 ** gamdh1) / (rhocrit3 ** gamah1)
          IF (rho(ipart).LT.rhocrit) THEN
             vxyzu(4,ipart) = uzero
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(pr(ipart)/rho(ipart))
          ELSE IF (rho(ipart).LT.rhocrit2) THEN
             vxyzu(4,ipart) = RK3 * (rho(ipart) ** gam1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(gam * pr(ipart) / rho(ipart))
          ELSE IF (rho(ipart).LT. rhocrit3) THEN
             vxyzu(4,ipart) = RK4 * (rho(ipart) ** gamdh1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamdh * pr(ipart) / rho(ipart))
          ELSE 
             vxyzu(4,ipart) = RK5 * (rho(ipart) ** gamah1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamah * pr(ipart) / rho(ipart))
          ENDIF
c
c--Physical equation of state
c
      ELSEIF (encal.EQ.'x') THEN
          gam1 = gamphys1 - 1.
          gam2 = gamphys2 - 1.
          gam3 = gamphys3 - 1.
          IF (iorig(ipart).GT.n1+1 .AND. iorig(ipart).LE.n1+n2) THEN
             uzero = uzero_n2
          ELSE
             uzero = RK2 * (rhozero ** gama1)
          ENDIF
          IF (rho(ipart).LT.rhochange1) THEN
             vxyzu(4,ipart) = uzero*(1.0 + (rho(ipart)/rhoref1)**gam1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamphys1*pr(ipart)/rho(ipart))
          ELSE IF (rho(ipart).LT.rhochange2) THEN
             vxyzu(4,ipart) = uzero*(1.0 + (rho(ipart)/rhoref2)**gam2)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamphys2* pr(ipart) / rho(ipart))
          ELSE 
             vxyzu(4,ipart) = uzero*(1.0 + (rho(ipart)/rhoref3)**gam3)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamphys3* pr(ipart) / rho(ipart))
          ENDIF
c
c--Isothermal equation of state
c
      ELSEIF (encal.EQ.'i') THEN
            pr(ipart) = (2./3.) * vxyzu(4,ipart) * rho(ipart)
            vsound(ipart) = SQRT(pr(ipart)/rho(ipart))
c
c--Thermal instability equation of state
c
      ELSEIF (encal.EQ.'t') THEN
            pr(ipart) = (2./3.) * vxyzu(4,ipart) * rho(ipart)
            vsound(ipart) = SQRT(pr(ipart)/rho(ipart))
c
c--Adiabatic equation of state
c
      ELSEIF (encal.EQ.'a' .OR. encal.EQ.'c') THEN
            pr(ipart) =  gama1 * vxyzu(4,ipart) * rho(ipart)
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
c
c--Radiative transfer equation of state
c
      ELSEIF (encal.EQ.'r') THEN
c           pr(ipart)= Rg*rho(ipart)*(vxyzu(4,ipart)/ekcle(3,ipart))*
c     &        get1overmu(rho(ipart),vxyzu(4,ipart))/uergg
         pr(ipart)= Rg*rho(ipart)*(vxyzu(4,ipart)/
     &        getcv(rho(ipart),vxyzu(4,ipart)))*
     &        get1overmu(rho(ipart),vxyzu(4,ipart))/uergg
         vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
c
c--Radiative transfer equation of state using Monte Carlo code
c
      ELSEIF (encal.EQ.'m') THEN
           pr(ipart)= Rg*rho(ipart)*(vxyzu(4,ipart)/
     &        getcv(rho(ipart),vxyzu(4,ipart)))*
     &        get1overmu(rho(ipart),vxyzu(4,ipart))/uergg
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
c
c--Polytropic equation of state
c
      ELSEIF (encal.EQ.'p') THEN
            vxyzu(4,ipart) = RK2 * (rho(ipart) ** gama1)
            pr(ipart) = (2./3.) * vxyzu(4,ipart) * rho(ipart)
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
c
c--Variable of state is entropy
c
      ELSEIF (varsta.EQ.'entropy') THEN
            pr(ipart) = vxyzu(4,ipart)*rho(ipart)**gamma
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
      ELSE
         WRITE (iprint,99500) encal
99500    FORMAT ('encal = ',A1)
         CALL error(where,1)
      ENDIF

      RETURN
      END
