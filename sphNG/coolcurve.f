      SUBROUTINE coolcurve(ucool,rhocool,dtcool)
c************************************************************
c                                                           *
c     This subroutine compares the temperature change       *
c     as found in energ with the equilibrium temperature    *
c     from the corresponding cooling curve. This subroutine *
c     follows the implementation of Vazquez-Semadeni et.al. *
c     to correct the value of u if du is very large.        *
c************************************************************

      INCLUDE 'COMMONS/tcooling'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'

      REAL tempiso, teq, ueq, coolingtime, coolingr, rhon
      REAL*4 rhocool
      INTEGER ic

c Find temperature in K

      tempiso = 2./3.*ucool/(Rg/gmw/uergg)

c Cooling rate from Vazquez-Semadeni

      coolingr=(1.e7*exp(-1.148e5/(tempiso+1000.))+
     &   0.014*sqrt(tempiso)*exp(-92./tempiso))*heatingr

 
c Find equilibrium temperature from previously defined values

      rhon=rhocool*udens/mpenerg
      DO ic=1,35000
      IF (rhon.GE.cooling(ic)) THEN
       teq=(rhon-cooling(ic))/(tguess(ic+1)-tguess(ic))+tguess(ic)
       GOTO 446
      END IF
      END DO

      ic=35000
      teq=(rhon-cooling(ic))/(tguess(ic+1)-tguess(ic))+tguess(ic) 
  
 446  ueq=1.5*teq*Rg/gmw/uergg
         
      coolingtime=ABS((ucool-ueq)/
     &   ((coolingr*rhocool*udens/mpenerg-heatingr)*utime**3/udist**2.))

      IF (ucool.LE.0.) THEN
           ucool=ueq
      ELSE
           ucool=ueq+(ucool-ueq)*exp(-dtcool/coolingtime)
      ENDIF           

      IF (coolingtime.EQ.0.) ucool=ueq

      RETURN
      END
