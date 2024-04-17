      SUBROUTINE eospg(ipart,vxyzu,rho,pr,vsound,ekcle,dustvar)
c************************************************************
c                                                           *
c  This routine computes the pressure and sound speed       *
c     according to a perfect gas equation of state          *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      REAL*4 rho,pr,vsound
      DIMENSION vxyzu(4,idim), rho(idim), pr(idim), vsound(idim)
      DIMENSION ekcle(5,iradtrans),dustvar(ndusttypes,idim_dustFluid)

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

      DIMENSION dustfraci(ndusttypes)
      SAVE dustfraci
C$OMP THREADPRIVATE(dustfraci)

      CHARACTER*7 where
      REAL*4 rho_gas_r4
      DATA where/'eospg'/

      gama1 = gamma - 1.

      IF (vxyzu(4,ipart).LE.0) THEN
         WRITE (*,*) 'ERROR - vxyzu(4,ipart).LE.0 ',vxyzu(4,ipart),ipart
         CALL quit(1)
      ENDIF
c
c--For one fluid dust, need to reduce the gas density based on the
c     dust fraction since the density is the total density, not the
c     density of the gas
c
      IF (idustFluid.NE.0) THEN
         CALL dustvar2dustfrac(1.0/rho(ipart),dustvar(:,ipart),
     &        dustfraci(:),dustfracisum)
         rho_gas = rho(ipart)*(1. - dustfracisum)
      ELSE
         rho_gas = rho(ipart)
      ENDIF
      rho_gas_r4 = rho_gas
c
c--Variable is internal energy
c
c--Variable gamma equation of state
c
      IF (encal.EQ.'v') THEN
          gam1 = gam - 1.
          gamdh1 = gamdh - 1.
          gamah1 = gamah - 1.
          IF (iunique(iorig(ipart)).GT.n1 .AND. 
     &         iunique(iorig(ipart)).LE.n1+n2) THEN
             uzero = uzero_n2
          ELSE
             uzero = RK2 * (rhozero ** gama1)
          ENDIF
          RK3 = uzero / (rhocrit ** gam1)
          RK4 = RK3 * (rhocrit2 ** gam1) / (rhocrit2 ** gamdh1)
          RK5 = RK4 * (rhocrit3 ** gamdh1) / (rhocrit3 ** gamah1)
          IF (rho_gas.LT.rhocrit) THEN
             vxyzu(4,ipart) = uzero
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(pr(ipart)/rho_gas)
          ELSE IF (rho_gas.LT.rhocrit2) THEN
             vxyzu(4,ipart) = RK3 * (rho_gas ** gam1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(gam * pr(ipart) / rho_gas)
          ELSE IF (rho_gas.LT. rhocrit3) THEN
             vxyzu(4,ipart) = RK4 * (rho_gas ** gamdh1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(gamdh * pr(ipart) / rho_gas)
          ELSE 
             vxyzu(4,ipart) = RK5 * (rho_gas ** gamah1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(gamah * pr(ipart) / rho_gas)
          ENDIF
c
c--Physical equation of state
c
      ELSEIF (encal.EQ.'x') THEN
          gam1 = gamphys1 - 1.
          gam2 = gamphys2 - 1.
          gam3 = gamphys3 - 1.
          IF (iunique(iorig(ipart)).GT.n1 .AND. 
     &         iunique(iorig(ipart)).LE.n1+n2) THEN
             uzero = uzero_n2
          ELSE
             uzero = RK2 * (rhozero ** gama1)
          ENDIF
          IF (rho_gas.LT.rhochange1) THEN
             vxyzu(4,ipart) = uzero*(1.0 + (rho_gas/rhoref1)**gam1)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(gamphys1*pr(ipart)/rho_gas)
          ELSE IF (rho_gas.LT.rhochange2) THEN
             vxyzu(4,ipart) = uzero*(1.0 + (rho_gas/rhoref2)**gam2)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(gamphys2* pr(ipart) / rho_gas)
          ELSE 
             vxyzu(4,ipart) = uzero*(1.0 + (rho_gas/rhoref3)**gam3)
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(gamphys3* pr(ipart) / rho_gas)
          ENDIF
c
c--Isothermal equation of state
c
      ELSEIF (encal.EQ.'i') THEN
            pr(ipart) = (2./3.) * vxyzu(4,ipart) * rho_gas
            vsound(ipart) = SQRT(pr(ipart)/rho_gas)
c
c--Thermal instability equation of state
c
      ELSEIF (encal.EQ.'t') THEN
            pr(ipart) = (2./3.) * vxyzu(4,ipart) * rho_gas
            vsound(ipart) = SQRT(pr(ipart)/rho_gas)
c
c--Adiabatic equation of state
c
      ELSEIF (encal.EQ.'a' .OR. encal.EQ.'c') THEN
            pr(ipart) =  gama1 * vxyzu(4,ipart) * rho_gas
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho_gas)
c
c--Radiative transfer equation of state
c
      ELSEIF (encal.EQ.'r') THEN
          IF (iunique(iorig(ipart)).GT.n1 .AND. 
     &         iunique(iorig(ipart)).LE.n1+n2) THEN
             vxyzu(4,ipart) = uzero_n2
             pr(ipart) = 2./3. * vxyzu(4,ipart) * rho_gas
             vsound(ipart) = SQRT(pr(ipart)/rho_gas)
          ELSE
             pr(ipart)= Rg*rho_gas*(vxyzu(4,ipart)/
     &            getcv(rho_gas_r4,vxyzu(4,ipart)))*
     &            get1overmu(rho_gas_r4,vxyzu(4,ipart))/uergg
             vsound(ipart) = SQRT(gamma*pr(ipart)/rho_gas)
          ENDIF
c
c--Radiative transfer equation of state using Monte Carlo code
c
      ELSEIF (encal.EQ.'m') THEN
           pr(ipart)= Rg*rho_gas*(vxyzu(4,ipart)/
     &        getcv(rho_gas_r4,vxyzu(4,ipart)))*
     &        get1overmu(rho_gas_r4,vxyzu(4,ipart))/uergg
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho_gas)
c
c--Polytropic equation of state
c
      ELSEIF (encal.EQ.'p') THEN
            vxyzu(4,ipart) = RK2 * (rho_gas ** gama1)
            pr(ipart) = (2./3.) * vxyzu(4,ipart) * rho_gas
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho_gas)
c
c--Variable of state is entropy
c
      ELSEIF (varsta.EQ.'entropy') THEN
            pr(ipart) = vxyzu(4,ipart)*rho_gas**gamma
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho_gas)
      ELSE
         WRITE (iprint,99500) encal
99500    FORMAT ('encal = ',A1)
         CALL error(where,1)
      ENDIF

      RETURN
      END
