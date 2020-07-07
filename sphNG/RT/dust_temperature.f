      FUNCTION dust_temperature(ipart,ntot,uradconst,Erad,Ugas,
     &     ekcle,rhoi,dust_kappai,dust_cooling,heatingISRi,dust_gas)
c
c--This function calculates the dust temperature by balancing the
c     dust heating from the interstellar radiation (ISR) field
c     (including extinction), the dust cooling/heating via radiation 
c     (taking into account the (grey) local radiation field -- i.e. not 
c     including the ISR), and the dust cooling/heating via dust-gas 
c     collisions.
c
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'idim'
      INCLUDE 'igrape'      

      PARAMETER (sboltzmann = 5.67E-5)

      REAL*8 gas_dust_collisional_term
      REAL dust_temperature,dust_cooling
      REAL uradconst,Erad,Ugas,rhoi,dust_kappai
      REAL ekcle(5,iradtrans2)
      REAL*4 rho_real4

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/rbnd'

c      COMMON /tracenp/ ntrace

      rho_real4 = rhoi
c
c--The value of heatingISR is the value of the heating rate of the dust
c     due to the ISR.  However, it only gives the extincted value of 
c     integral of Q_v*J_v, where Q_v is in units of cm^2/H_2 (from 
c     Zucconi et al. 2001).
c     For the dust emission term, the dust_kappa value which is used 
c     is the Planck mean opacity in units of cm^2/g.  So we need to 
c     divide the heatingISR value by mu*mH change Q_v into cm^2/g.  
c     The 4*pi comes from the fact that a*c = 4*sigma_B and the Planck 
c     mean opacity is normalised by the integral over the Planck 
c     function which is sigma_B*T^4/pi.  So a*c times one over the 
c     normalisation is 4*pi.
c
c     The heatingISR must be balanced by the net cooling from the 
c     interaction with the non-ISR radiation field and collisional 
c     cooling with the gas.
c
c--heatingISR is in erg/s/H_2 so must be divided by mu*mH to get 
c     into erg/g/s
c
      heatingISRi = G0*heatingISR(1,ipart) * 4.0*pi /(gmw*mH)
c
c--We assume that the net cooling from the interaction with the non-ISR
c     radiation can be expressed as
c
c     a*c*kappa*(T_r^4-T_d^4)
c     which takes the Planck mean opacities and has units of erg/g/s
c
c     while the latter is (from Zucconi et al. 2001 or Goldsmith 2001)
c     10^-33 * n_H2^2 * T_gas^0.5 *(T_dust - T_gas) / rho
c
c     where we have divided by rho to convert from erg/cm^3/s to erg/g/s
c
      constant1 = radconst*c
      xnH2 = rhoi*udens/(gmw*mH)
      gas_temp = Ugas/ekcle(3,ipart)

c      IF (heatingISRi.LE.0.0) WRITE (*,*) 'ERROR - heatingISRi ',
c     &     ipart,ntrace,gas_temp,xnH2,heatingISRi
c
c--If the gas is dense and hot, then the gas and dust will be well 
c     coupled thermally, so just take the dust temperature to be equal
c     to the gas temperature
c
      IF (xnH2.GT.1.0E+12/metallicity .AND. gas_temp.GE.100.) THEN
         t_found = gas_temp
         index = INT(t_found*10.0)
         IF (index.GT.nPlanckPoints) THEN
c            WRITE (*,*) 'ERROR - off end of Planck Opacity table'
            index = nPlanckPoints-1
         ELSEIF (index.LT.1) THEN
            index = 1
         ENDIF
         dust_kappai = Planck_table(index)
         dust_gas = 0.
         dust_cooling = constant1*dust_kappai*t_found**4
         GOTO 222
      ENDIF
c
c--Set constant2 for gas-dust cooling term (from Zucconi et al. 2001,
c     but modified by the metallicity, assuming that the dust number 
c     density scales linearly with the metallicity).
c
      constant2 = gas_dust_collisional_term(xnH2,metallicity,gas_temp)/
     &     (rhoi*udens)
      rad_temp4 = rhoi*Erad/uradconst
c      rad_temp4 = 0.
c
c--If both radiation and gas temperatures are >1000, then set dust
c     temperature directly
c
      IF (rad_temp4.LT.0.) THEN
         WRITE (*,*) 'ERROR - radiation temp. in dust_temperature',
     &        ' is < 0 ',rad_temp4,ipart,gas_temp
c         CALL quit(1)
      ELSE
         rad_temp = rad_temp4**(0.25)
         IF (rad_temp.GT.1000. .AND. gas_temp.GT.1000.) THEN
            t_found = rad_temp
            index = INT(t_found*10.0)
            IF (index.GT.nPlanckPoints) THEN
               index = nPlanckPoints-1
            ELSEIF (index.LT.1) THEN
               index = 1
            ENDIF
            dust_kappai = Planck_table(index)
            dust_gas = 0.
            dust_cooling = constant1*dust_kappai*t_found**4
            GOTO 222
         ENDIF
      ENDIF
c
c--Find initial dust temperature by assuming that the dust is in thermal
c     balance with the radiation field (i.e. ignore gas-dust collisional
c     cooling).  This will be good if the density is low.
c
      iattempt = 1
 200  t_found = MIN(gas_temp,1000.0)
      DO i = 1, 100
         t_last = t_found
c
c--Dust opacity as Planck mean in cm^2/g.
c
         index = INT(t_found*10.0)
         IF (index.GT.nPlanckPoints) THEN
c            WRITE (*,*) 'ERROR - off end of Planck Opacity table'
            index = nPlanckPoints-1
         ELSEIF (index.LT.1) THEN
            index = 1
         ENDIF
         dust_kappai = Planck_table(index)
c
c--T_dust assuming balance with ISR and FLD radiation fields.  
c     The equation being solved, essentially is 
c     Integral(Q_v*B_v) = Integral(Q_v*J_ISR).
c     Remember, that a=4*sigma/c.
c
         t_found = (heatingISRi / (radconst*c*dust_kappai) + 
     &        rad_temp4)**0.25

         IF (ABS(t_found-t_last).LT.0.01) GOTO 221
c         IF (ABS(t_found-t_last).LT.0.01) GOTO 222
      END DO
 221  CONTINUE
      IF (iattempt.EQ.2) GOTO 222

      t_found = MIN(t_found,1000.0)
c
c--Now find actual solution (including gas-dust collisional cooling)
c     using Newton-Raphson method
c
      DO i = 1, 1000
         t_last = t_found

         index = INT(t_found*10.0)
         IF (index.GT.nPlanckPoints) THEN
c            WRITE (*,*) 'ERROR - off end of Planck Opacity table'
            index = nPlanckPoints-1
         ELSEIF (index.LT.1) THEN
            index = 1
         ENDIF
         dust_kappai = Planck_table(index)
c
c--Need to determine dkappa/dT for the derivative of the function.
c     The "10x" is because the Planck mean opacity table is stored
c     in increments of 1/10 K.
c
         dkappa_dust_dT = 10.0*(Planck_table(index+1) - dust_kappai)

         dust_cooling = constant1*dust_kappai*t_found**4

         dust_gas = constant2*(t_found - gas_temp)

         func = heatingISRi + 
     &        constant1*dust_kappai*(rad_temp4 - t_found**4) -
     &        constant2*(t_found - gas_temp)
         derivative = constant1*dkappa_dust_dT*(rad_temp4 - t_found**4)-
     &        4.0*constant1*dust_kappai*t_found**3 - constant2

         IF (func/derivative / t_found .GT. 0.1) THEN
            func = 0.1*derivative*t_found 
         ELSEIF (func/derivative / t_found .LT. -0.1) THEN
            func = -0.1*derivative*t_found 
         ENDIF
         t_found = t_found - func/derivative

         IF (ABS(t_found-t_last).LT.0.05
     &        .AND. 
     &        ABS(func).LT.5.0E-2*(ABS(heatingISRi)+
     &        ABS(constant1*dust_kappai*(rad_temp4 - t_found**4))+
     &        ABS(constant2*(t_found - gas_temp)))
     &        ) GOTO 222
      END DO
      WRITE (*,*) 'ERROR - dust_temperature failed ',t_found,t_last,
     &     func,(ABS(heatingISRi)+
     &     ABS(constant1*dust_kappai*(rad_temp4 - t_found**4))+
     &     ABS(constant2*(t_found - gas_temp))),
     &     heatingISRi,constant1*dust_kappai*(rad_temp4 - t_found**4),
     &     constant2*(t_found - gas_temp),constant1,dust_kappai,
     &     rad_temp4,constant2,gas_temp,Ugas/ekcle(3,ipart),ipart
c      CALL quit
c      t_found = 20.
      iattempt = 2
      GOTO 200

 222  dust_temperature = t_found
      IF (t_found.GT.1500.0) THEN
c         WRITE (*,*) 'ERROR - dust temperature > 1500 K'
c         CALL quit
      ENDIF
c
c--Convert dust_kappa to code units on exit
c
      dust_kappai = dust_kappai/udist**2*umass

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION gas_temperature_equilibrium(ipart,dust_temp,rhoi,
     &     photoelectric,cooling,gas_dust,cosmic_ray,func,derivative)
c
c--This function calculates the INITIAL gas temperature by assuming a 
c     balance between cosmic ray heating, line cooling, and transfer of 
c     energy between the dust and gas via collisions.
c     Assumes that pdV is zero, and that the continuum gas opacity is
c     negligible (i.e. T < ~1000 K).
c     It needs a pre-calculated (input) dust temperature.
c     The quantities in this subroutine are in cgs units, NOT code units
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'      

      REAL*8 gas_dust_collisional_term,cosmic_ray_heating,
     &     photoelectric_heating
      REAL gas_temperature_equilibrium
      REAL dust_temp, rhoi, photoelectric,cosmic_ray
      REAL t_found,xNH2,value1,t_last,cooling,sqrt_t,gas_dust
      REAL func,dcooling_dT,derivative,cooling_line_rate
      REAL dphotoelectric_dT,func_use
      INTEGER i,ipart

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/rbnd'
c
c--Initial guess temperature
c
      t_found = dust_temp
      xnH2 = rhoi*udens/(gmw*mH)

      cosmic_ray = cosmic_ray_heating(xnH2)
c
c--Now find actual solution using Newton-Raphson method
c
      DO i = 1, 10000
         t_last = t_found
         cooling = cooling_line_rate(ipart,t_found,xnH2,metallicity)
         photoelectric = photoelectric_heating(ipart,t_found,xnH2,
     &        metallicity)

         sqrt_t = SQRT(t_found)

         value1 = gas_dust_collisional_term(xnH2,metallicity,t_found)

         IF (t_found.GT.1500.) value1 = 0.



c         value1 = 0.





         func = value1*(dust_temp-t_found)           ! collisions
     &        - cooling                              ! line cooling
     &        + cosmic_ray                           ! cosmic rays
     &        + photoelectric                        ! photoelectric
c
c--Use numerical derivative to calculation dcooling/dT
c
         dcooling_dT = cooling_line_rate(ipart,t_found+1.0,xnH2,
     &        metallicity) - cooling
         dphotoelectric_dT = photoelectric_heating(ipart,t_found+1.0,
     &        xnH2,metallicity) - photoelectric

         derivative = 0.5*value1/t_found*dust_temp - 
     &        1.5*value1 - dcooling_dT + dphotoelectric_dT

c         IF (ipart.EQ.181) THEN
c      WRITE (*,*) 'iterating ',ipart,i,
c     &     t_found,t_last,(ABS(func/MAX(photoelectric,cosmic_ray,
c     &        value1*(dust_temp-t_found)))),
c     &     ABS(t_found-t_last)/MAX(t_found,t_last),func,
c     &     MAX(photoelectric,cosmic_ray,
c     &        value1*(dust_temp-t_found)),photoelectric,cosmic_ray,
c     &      value1*(dust_temp-t_found),cooling,dcooling_dT,
c     &           dphotoelectric_dT,derivative
c      ENDIF
      
         IF (func/derivative / t_found .GT. 0.1) THEN
            func_use = 0.1*derivative*t_found
         ELSE
            func_use = func
         ENDIF
         t_found = t_found - func_use/derivative
c         IF (t_found.GT.1500.0) THEN
c            WRITE (*,*) 'ERROR - Initial dust temperature >1500 K'
c            CALL quit
c         ENDIF
         IF (ABS(func/MAX(ABS(photoelectric),ABS(cooling),
     &        ABS(cosmic_ray),
     &        ABS(value1*(dust_temp-t_found)))).LT.1.0E-3 .AND. 
     &        (ABS(t_found-t_last).LT.0.01 .OR. t_found.GT.100.0 .AND. 
     &        ABS(t_found-t_last)/MAX(t_found,t_last).LT.1E-3)) GOTO 223
      END DO
      WRITE (*,*) 'ERROR - gas_temperature_equilibrium failed ',ipart,
     &     t_found,t_last,(ABS(func/MAX(photoelectric,cosmic_ray,
     &        value1*(dust_temp-t_found)))),
     &     ABS(t_found-t_last)/MAX(t_found,t_last),func,
     &     MAX(photoelectric,cosmic_ray,
     &        value1*(dust_temp-t_found)),photoelectric,cosmic_ray,
     &      value1*(dust_temp-t_found),dcooling_dT,derivative
      CALL quit(1)


 223  gas_temperature_equilibrium = t_found
      gas_dust = value1*(dust_temp-t_found)

c      WRITE (88,991) xnH2,cooling,photoelectric,dcooling_dT,
c     &     t_found,func
c 991  FORMAT(6(1PE12.5,1X))
      RETURN
      END

c-----------------------------------------------------------

      SUBROUTINE build_interstellar_radiation_table
c
c--Calculate integral of Q_v * J_v with extinction for absorption rate from
c     ISM radiation field.  Builds a lookup table with extinction values
c     from A_V = 10^-3 to A_V = ~1200 in log-space.
c
c--This is based on Zucconi et al. (2001)
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'idim'
      INCLUDE 'igrape'      

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'

      xISRfac = EXP(xISRlogfac)
      dlogv = LOG(1.1)

      A_V_init = 1.0E-3
      DO j = 1, nISRpoints
         A_V = A_V_init*xISRfac**(j-1)

         x_integral = 0.
         DO i = 1, 1000
            v = (1.1**i)
            v = v / 1.0E+10
            xJ = xJ_is(v)
            IF (xJ.GT.1e-40 .AND. v.GT.0.1) THEN
               x_integral = x_integral + EXP(-A_V/1.086*Qv(v)/
     &              Qv(c/0.0000550) )*Qv(v)*xJ*v*dlogv
            ENDIF
c            WRITE (23,*) v, xJ, Qv(v)
c            WRITE (24,*) v, exp(-10.0/1.086*Qv(v)/
c     &           Qv(c/0.0000550) )*Qv(v)*xJ*v*dlogv
         END DO
         ISR_table(j) = x_integral
c         print *,'IS value ',j,x_integral
      END DO

      RETURN
      END

c-----------------------------------------------------------

      SUBROUTINE build_photoelectric_table
c
c--Calculate Ge, which is the integral of exp(-tau) * J_v, but only above 6 eV 
c     (for photons that produce photoelectic heating of the gas) with 
c     extinction for absorption rate from ISM radiation field.  
c     Builds a lookup table with extinction values
c     from A_V = 10^-3 to A_V = ~1200 in log-space.
c
c--This is based on Keto & Casselli (2008) and Young et al. (2004).
c
c     The table value is normalised (i.e. with no extinction it is unity).
c     The x_integral value is per steradian (i.e. it is not multiplied by 4pi
c     and should has a value of about 2.3x10^-4 erg/s/cm^2/steradian (for
c     the provided xJ) whereas the Habing (1968) value is about 1.3x10^-4
c     and the Drain (1978) value is about 2.2x10^-4.
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'idim'
      INCLUDE 'igrape'      

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'

      xISRfac = EXP(xISRlogfac)
      dlogv = LOG(1.1)

      A_V_init = 1.0E-3
      DO j = 1, nISRpoints
         A_V = A_V_init*xISRfac**(j-1)

         x_integral = 0.
         x_norm = 0.
         DO i = 1, 1000
            v = (1.1**i)
            v = v / 1.0E+10
            xJ = xJ_is(v)
c
c--Only calculate the integral from 6 to 13.6 eV.
c
            IF (xJ.GT.1e-40 .AND. v.GT.6.0*eleccharge/planckconst 
     &          .AND.  v.LT.13.6*eleccharge/planckconst
     &           ) THEN
               x_integral = x_integral + EXP(-A_V/1.086*Qv(v)/
     &              Qv(c/0.0000550) )*xJ*v*dlogv
               x_norm = x_norm + xJ*v*dlogv
            ENDIF
         END DO
         Ge_ISR_table(j) = x_integral/x_norm
c         print *,'Ge IS value ',j,x_integral,Ge_ISR_table(j)
      END DO

      RETURN
      END

c-----------------------------------------------------------

      SUBROUTINE build_planck_mean_table
c
c--Calculate integral of Q_v * B_v (absorption times Planck function),
c     normalised by the integral of the Planck function (sigma*T^4/pi)
c     and divided by gmw*mH because the absorption is given per H_2
c     molecule and we want opacity in cm^2/g (gmw*mH is the mass in g
c     per H_2 molecule).
c
c     Builds a lookup table for temperature, T.
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'idim'
      INCLUDE 'igrape'      

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'

      dlogv = LOG(1.1)

      T_init = 0.1
      DO j = 1, nPlanckPoints
         T = T_init*j

         x_integral = 0.
         DO i = 1, 1000
            v = (1.1**i)
            v = v / 1.0E+10
            xJ = planck(v,T)
            IF (xJ.GT.1e-40 .AND. v.GT.0.1) THEN
               x_integral = x_integral + Qv(v)*xJ*v*dlogv
            ENDIF
         END DO
         x_integral = x_integral*pi/(T**4*stefanboltz)/(gmw*mH)
         Planck_table(j) = x_integral

c         IF (j.LE.nPlanckPoints) print *,'Planck value ',j,
c     &        Planck_table(j)
c         WRITE (21,*) T,x_integral
      END DO

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION xJ_is(v)
c
c--This is the Interstellar Radiation (ISR) field
c--From Zucconi et al. (2001), following Black (1994)
c     But adds Draine (1978) UV field (equation 11) which is similar to
c     Black (1994) but doesn't seem to be included in Zucconi.
c
c--NOTE: The values are in cgs units, not code units.
c
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/interstellar'

      REAL*8 xJ_is, v
      REAL*8 xlambda(5), power(5), weight(5), temperature(5)
      REAL*8 hoverk, val, vp, E_eV, val2, val3
      INTEGER i

      xlambda(1) = 0.4E-4
      xlambda(2) = 0.75E-4
      xlambda(3) = 1.0E-4
      xlambda(4) = 140.0E-4
      xlambda(5) = 1.06E-1 * (1.0 + redshift)

      power(1) = 0.
      power(2) = 0.
      power(3) = 0.
      power(4) = 1.65
      power(5) = 0.

      weight(1) = 1.0E-14
      weight(2) = 1.0E-13
      weight(3) = 4.0E-13
      weight(4) = 2.0E-4
      weight(5) = 1.0

      temperature(1) = 7500.0
      temperature(2) = 4000.0
      temperature(3) = 3000.0
      temperature(4) = 23.3
      temperature(5) = 2.728 * (1.0 + redshift)

      hoverk = planckconst/boltzmannk

      val = 0.
      DO i = 1, 5
         val3 = hoverk*v/temperature(i)
         IF (val3.LT.600.) THEN
            val2 = EXP( hoverk*v/temperature(i)) - 1.0
            IF (val2.GT.0.) THEN
               val = val + (xlambda(i)*v/c)**power(i) * weight(i)
     &              / val2
            ENDIF
         ENDIF
      END DO
      val = val * 2.0*((planckconst*v)*v)*v/c**2
c
c--Add mid-infrared which has a cut-off longer than 100 microns
c
      IF (iISR_MIR .AND. c/v.LT.100.0E-4) THEN
         vp = c/100.0E-4
         val = val + 5.0E-7* 2.0*planckconst*(vp/c)**2*vp * 
     &        (v/vp)**(-1.8)
      ENDIF


c      val = val/10.

c
c--Add Draine 1978 UV
c
      IF (iISR_Draine .AND. 
     &     v.GT.5.0*eleccharge/planckconst .AND. 
     &     v.LT.13.6*eleccharge/planckconst) THEN
         E_eV = v*planckconst/eleccharge
         val = val + (1.658E+06*E_eV - 2.152E+05*E_eV**2 + 
     &        6.919E+03*E_eV**3)*planckconst*E_eV
      ENDIF

      xJ_is = val

      RETURN
      END

c-----------------------------------------------------------


      FUNCTION cooling_line_rate(ipart,gas_temp,xnH2,metallicity)
c
c--Molecular line cooling rate from Goldsmith (2001) Table 2 (also used by
c     Keto & Field 2005).  In cgs units, provides cooling: erg/s/cm^3.
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/interstellar'

      REAL cooling_line_rate, gas_temp, xnH2, logT, metallicity
      REAL logxnH2, fac, cooling_log1, cooling_log2
      REAL diff1, diff2, diff3, diff4, xdiff12, xdiff13, xdiff23
      REAL ydiff1, ydiff2, ydiff12, ypos1, ypos2, tenlog
      REAL log_depletion
      REAL xdiff24, xdiff34, cooling_log3, cooling_log4
      REAL xpos1,xpos2,xpos3,xpos4
      REAL xnH, xne, brackets, Gr, recombination, phiPAH
      REAL oxygen
      REAL electron_fraction, x_Cplus
      REAL carbon_chemistry, depletion, cooling_line_rate_dep, x_CO
      INTEGER ipos,ipos_depletion,ipart

      REAL Goldsmith2001_table2(3,7), Goldsmith2001_depletion(2,5,4)

      cooling_line_rate = 0.

      IF (iDRT_line) THEN

         Goldsmith2001_table2(1,1) = LOG(6.3E-26)
         Goldsmith2001_table2(1,2) = LOG(3.2E-25)
         Goldsmith2001_table2(1,3) = LOG(1.1E-24)
         Goldsmith2001_table2(1,4) = LOG(5.6E-24)
         Goldsmith2001_table2(1,5) = LOG(2.3E-23)
         Goldsmith2001_table2(1,6) = LOG(4.9E-23)
         Goldsmith2001_table2(1,7) = LOG(7.4E-23)

         Goldsmith2001_table2(2,1) = 1.4
         Goldsmith2001_table2(2,2) = 1.8
         Goldsmith2001_table2(2,3) = 2.4
         Goldsmith2001_table2(2,4) = 2.7
         Goldsmith2001_table2(2,5) = 3.0
         Goldsmith2001_table2(2,6) = 3.4
         Goldsmith2001_table2(2,7) = 3.8

         Goldsmith2001_table2(3,1) = 2.0
         Goldsmith2001_table2(3,2) = 2.5
         Goldsmith2001_table2(3,3) = 3.0
         Goldsmith2001_table2(3,4) = 4.0
         Goldsmith2001_table2(3,5) = 5.0
         Goldsmith2001_table2(3,6) = 6.0
         Goldsmith2001_table2(3,7) = 7.0

         logxnH2 = LOG10(xnH2)
         logT = LOG(gas_temp/10.0)
c
c--Work out cooling rates at the three points surrounding the requested point
c     and then use three-point polynomial interpolation from those rates 
c     (in log space).  For an explanation of polynomial interpolation see 
c     Numerical Recipes 2, Section 3.1
c
c     Straight 3-point polynomial interpolation works pretty well, but
c     gives a "bulge" between log10(n_H2)=3-4.  So to fix this I calculate
c     two polynomial interpolants using the points at (2.5, 3, 4) and
c     (3,4,5) and then use a linear fraction (still in log space) of 
c     these between 3 and 4 (i.e. smoothly going from one to the other).
c     This works really nicely.
c
c     Using linear interpolation of the alpha and beta values or even of the
c     actually cooling values doesn't give a smooth resulting temperature 
c     profile because the slope of beta is has a strong discontinuity 
c     at log(n_H2)=3.
c
         IF (logxnH2.LT.3.0) THEN
            IF (logxnH2.LT.2.5) THEN
               ipos = 1
            ELSEIF (logxnH2.GE.2.5) THEN
               ipos = 2
            ENDIF
            xpos1 = Goldsmith2001_table2(3,ipos)
            xpos2 = Goldsmith2001_table2(3,ipos+1)
            xpos3 = Goldsmith2001_table2(3,ipos+2)
            fac = (logxnH2-xpos1)/0.5
            diff1 = logxnH2-xpos1
            diff2 = logxnH2-xpos2
            diff3 = logxnH2-xpos3
            xdiff12 = xpos1 - xpos2
            xdiff13 = xpos1 - xpos3
            xdiff23 = xpos2 - xpos3
c
c--Altered because of trouble with overcooling of gas at low metallicities
c     around nH2 = 10^7 - 10^8
c
c         ELSEIF (logxnH2.LT.7.3) THEN
         ELSEIF (logxnH2.LT.8.0) THEN
            IF (logxnH2.GT.7.0) THEN
               ipos = 5
            ELSE
               ipos = MAX(3.,logxnH2)-1
            ENDIF
            xpos1 = Goldsmith2001_table2(3,ipos)
            xpos2 = Goldsmith2001_table2(3,ipos+1)
            xpos3 = Goldsmith2001_table2(3,ipos+2)
            fac = logxnH2-xpos1
            diff1 = logxnH2-xpos1
            diff2 = logxnH2-xpos2
            diff3 = logxnH2-xpos3
            xdiff12 = xpos1 - xpos2
            xdiff13 = xpos1 - xpos3
            xdiff23 = xpos2 - xpos3
            
         ELSE
            cooling_line_rate = 0.
            GOTO 1000
         ENDIF
         cooling_log1 = Goldsmith2001_table2(1,ipos) + 
     &        Goldsmith2001_table2(2,ipos)*logT
         cooling_log2 = Goldsmith2001_table2(1,ipos+1) + 
     &        Goldsmith2001_table2(2,ipos+1)*logT
         cooling_log3 = Goldsmith2001_table2(1,ipos+2) + 
     &        Goldsmith2001_table2(2,ipos+2)*logT

         cooling_line_rate =diff2*diff3/(xdiff12*xdiff13)*cooling_log1+
     &        diff1*diff3/(-xdiff12*xdiff23)*cooling_log2 +
     &        diff1*diff2/(xdiff13*xdiff23)*cooling_log3

         IF (ipos.EQ.2 .AND. logxnH2.GT.3.0) THEN
            xpos4 = Goldsmith2001_table2(3,ipos+3)
            diff4 = logxnH2-xpos4
            xdiff24 = xpos2 - xpos4
            xdiff34 = xpos3 - xpos4
            cooling_log4 = Goldsmith2001_table2(1,ipos+3) +
     &           Goldsmith2001_table2(2,ipos+3)*logT
            cooling_line_rate = (1.0-diff2)*cooling_line_rate + diff2*(
     &           diff3*diff4/(xdiff23*xdiff24)*cooling_log2+
     &           diff2*diff4/(-xdiff23*xdiff34)*cooling_log3 +
     &           diff2*diff3/(xdiff24*xdiff34)*cooling_log4)

         ENDIF

         cooling_line_rate = EXP(cooling_line_rate)

         IF (logxnH2.LT.1.) cooling_line_rate = 0.

c         IF (logxnH2.GT.7.0 .AND. logxnH2.LE.7.3) cooling_line_rate =
c     &        cooling_line_rate * (7.3-logxnH2)/(0.3)
c
c--If gas is depleted, use Goldsmith (2001) depletion
c
         x_Cplus = carbon_chemistry(ipart, xnH2, gas_temp, depletion, 
     &        x_CO)
c
c--This table: index 1: log alpha or beta value
c              index 2: depletion factor: 1, 3, 10, 30, 100
c              index 3: nH2 density: 10^3, 10^4, 10^5, 10^6
c
c     If values lie outside these ranges, use undepleted values for
c     high densities and cooling rates that depend linearly on the depletion
c     for low densities
c
c         GOTO 777

         IF (logxnH2.GT.7.) THEN
c
c--Don't do anything since depletion doesn't change the cooling much 
c     according to Goldsmith (2001)
c
         ELSEIF (logxnH2.LT.2.) THEN
c
c--Assume that cooling rate scales linearly with depletion
c
            cooling_line_rate = cooling_line_rate*depletion
         ELSE
c
c--Use cooling from Table 4 of Goldmsith (2001)
c
            Goldsmith2001_depletion(1,1,1) = -24.0
            Goldsmith2001_depletion(1,2,1) = -24.2
            Goldsmith2001_depletion(1,3,1) = -24.5
            Goldsmith2001_depletion(1,4,1) = -24.8
            Goldsmith2001_depletion(1,5,1) = -25.1

            Goldsmith2001_depletion(1,1,2) = -23.3
            Goldsmith2001_depletion(1,2,2) = -23.4
            Goldsmith2001_depletion(1,3,2) = -23.5
            Goldsmith2001_depletion(1,4,2) = -23.7
            Goldsmith2001_depletion(1,5,2) = -23.8

            Goldsmith2001_depletion(1,1,3) = -22.6
            Goldsmith2001_depletion(1,2,3) = -22.8
            Goldsmith2001_depletion(1,3,3) = -22.9
            Goldsmith2001_depletion(1,4,3) = -23.0
            Goldsmith2001_depletion(1,5,3) = -23.2

            Goldsmith2001_depletion(1,1,4) = -22.3
            Goldsmith2001_depletion(1,2,4) = -22.4
            Goldsmith2001_depletion(1,3,4) = -22.4
            Goldsmith2001_depletion(1,4,4) = -22.5
            Goldsmith2001_depletion(1,5,4) = -22.6

            Goldsmith2001_depletion(2,1,1) = 2.4
            Goldsmith2001_depletion(2,2,1) = 2.1
            Goldsmith2001_depletion(2,3,1) = 1.9
            Goldsmith2001_depletion(2,4,1) = 1.8
            Goldsmith2001_depletion(2,5,1) = 1.8

            Goldsmith2001_depletion(2,1,2) = 2.7
            Goldsmith2001_depletion(2,2,2) = 2.7
            Goldsmith2001_depletion(2,3,2) = 2.6
            Goldsmith2001_depletion(2,4,2) = 2.5
            Goldsmith2001_depletion(2,5,2) = 2.2

            Goldsmith2001_depletion(2,1,3) = 3.0
            Goldsmith2001_depletion(2,2,3) = 2.9
            Goldsmith2001_depletion(2,3,3) = 2.8
            Goldsmith2001_depletion(2,4,3) = 2.8
            Goldsmith2001_depletion(2,5,3) = 2.8

            Goldsmith2001_depletion(2,1,4) = 3.4
            Goldsmith2001_depletion(2,2,4) = 3.3
            Goldsmith2001_depletion(2,3,4) = 3.2
            Goldsmith2001_depletion(2,4,4) = 3.2
            Goldsmith2001_depletion(2,5,4) = 3.1

            ipos = logxnH2 - 2
            IF (ipos.LT.1) THEN
               ipos = 1
            ELSEIF (ipos.GT.3) THEN
               ipos = 3
            ENDIF
            log_depletion = LOG10(depletion)
            ipos_depletion = 1 - log_depletion*2.0
            IF (ipos_depletion.LT.1) THEN
               ipos_depletion = 1
            ELSEIF (ipos_depletion.GT.4) THEN
               ipos_depletion = 4
            ENDIF
            xpos1 = ipos + 2.0
            xpos2 = xpos1 + 1.0
            diff1 = logxnH2 - xpos1
            diff2 = xpos2 - logxnH2
            xdiff12 = xpos2 - xpos1

            ypos1 = -0.5*(ipos_depletion - 1.0)
            ypos2 = -0.5*(ipos_depletion)
            ydiff1 = log_depletion - ypos1
            ydiff2 = ypos2 - log_depletion
            ydiff12 = ypos2 - ypos1
c
c--Values are given as log10, but interpolation is done using natural logs.
c
            tenlog = LOG(10.)

            cooling_log1 = tenlog*
     &           Goldsmith2001_depletion(1,ipos_depletion,ipos) + 
     &           Goldsmith2001_depletion(2,ipos_depletion,ipos)*logT
            cooling_log2 = tenlog*
     &           Goldsmith2001_depletion(1,ipos_depletion,ipos+1) + 
     &           Goldsmith2001_depletion(2,ipos_depletion,ipos+1)*logT
            cooling_log3 = tenlog*
     &           Goldsmith2001_depletion(1,ipos_depletion+1,ipos) + 
     &           Goldsmith2001_depletion(2,ipos_depletion+1,ipos)*logT
            cooling_log4 = tenlog*
     &           Goldsmith2001_depletion(1,ipos_depletion+1,ipos+1) + 
     &           Goldsmith2001_depletion(2,ipos_depletion+1,ipos+1)*logT

            cooling_line_rate_dep = (cooling_log1*diff2*ydiff2 + 
     &           cooling_log2*diff1*ydiff2 + 
     &           cooling_log3*diff2*ydiff1 +
     &           cooling_log4*diff1*ydiff1)/(xdiff12*ydiff12)
c
c--If outside of range of Table 4 of Goldsmith (2001), use linear interpolation
c     in log(nH2) to have smooth transition between depleted and undepleted
c     cooling rates
c
            IF (logxnH2.LT.3.0) THEN
               diff1 = logxnH2 - 2.0
               cooling_line_rate = EXP(
     &              LOG(cooling_line_rate*depletion)*(1.0-diff1) +
     &              cooling_line_rate_dep*diff1 )
            ELSEIF (logxnH2.GT.6.0) THEN
               diff1 = 7.0 - logxnH2
               cooling_line_rate = EXP(
     &              LOG(cooling_line_rate)*(1.0-diff1) +
     &              cooling_line_rate_dep*diff1 )
            ELSE
               cooling_line_rate = EXP(cooling_line_rate_dep)
            ENDIF
         ENDIF
c 777     CONTINUE
         cooling_line_rate = cooling_line_rate * metallicity
c     &        * 3.

c         cooling_line_rate = 0.

c
c--Add C+ cooling from Tielens (2005, ISM book)
c     NOTE: 0.14 is 1.4E-4/1.0E-3, which is the total carbon abundance
c     divided by the 1.0E-3 value in the Tielens formula
c
         xnH = 2.0 * xnH2

         cooling_line_rate = cooling_line_rate + 
     &        x_Cplus *
     &        3E-27*EXP(-92./gas_temp) * xnH**2 * metallicity
c
c--Add oxygen OI cooling from equation C3 of Wolfire et al. (2003)
c     Assume that OI abundance scales with x_CO abundance (not
c     including depletion).
c
         oxygen = 2.5E-27 * xnH**2 * (gas_temp/100.0)**0.4 * metallicity
         oxygen = oxygen*(1.0-x_CO)
c         oxygen = oxygen*(1.0-2.0*h2frac)
         IF (gas_temp.GT.4.0) THEN
            oxygen = oxygen * EXP(-228.0/gas_temp)
         ELSE
            oxygen = 0.
         ENDIF
c
c--Tapper oxygen cooling to zero from n_H2 = 1000 to 10000 
c
c         IF (logxnH2.GT.3.0 .AND. logxnH2.LE.4.0) oxygen =
c     &        oxygen * (4.0-logxnH2)
c         oxygen = 0.

         cooling_line_rate = cooling_line_rate + oxygen
c
c--Add recombination cooling (Glover & MacLow 2007; Wolfire et al 2003;
c     Bakes & Tielens 1994)
         xne = electron_fraction(xnH,phiPAH)
c
c--Needs the local UV ISR, but the recombination cooling should only
c     be significant when the extinction is very low anyway.
c
         Gr = G0*heatingISR(2,ipart)
c         Gr = 1.0
         brackets = Gr*SQRT(gas_temp)/xne/phiPAH

         recombination = 4.65E-30 * gas_temp**0.94 *
     &        brackets**(0.74/gas_temp**0.068)*xne*xnH/phiPAH *
     &        metallicity
c         recombination = 0.

         cooling_line_rate = cooling_line_rate + recombination

      ELSE  ! Cooling turned off completely
         cooling_line_rate = 0.
      ENDIF

 1000 RETURN
      END


c-----------------------------------------------------------
      FUNCTION h2_formation(ipart,gas_temp,dust_temp,xnH2)
c
c--Based on Glover et al. 2010 (rate 165 in appendix)
c     
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rbnd'

      INTEGER ipart
      REAL gas_temp,dust_temp,xnH2,grain_formation_rate
      REAL h2_formation,fA,fB

      REAL xnH
c
c--nH is twice nH2, where nH is actually the number density of protons from
c     hydrogen (i.e. it does not depend on the fractions of H and H_2 )
c
      xnH = 2.0*xnH2
c
c--Partial rate from Glover et al. 2010 (rate 165 in appendix)
c
      IF (dust_temp.GT.20.0) THEN
         fA = 1./(1.0 + 1.E+4*EXP(-600.0/dust_temp))
      ELSE
         fA = 1.0
      ENDIF
      fB = 1.0/(1.0 + 0.04*SQRT(gas_temp + dust_temp) + 0.002*gas_temp
     &     + 8.0E-06*gas_temp**2)
      grain_formation_rate = 3.0E-18*SQRT(gas_temp) * fA * fB *
     &     ((1.0-2.0*h2frac(ipart))*xnH)**2 * metallicity

      h2_formation = grain_formation_rate
c
c--Energy released is 4.48eV times the formation rate (eV in erg),
c     with H_2 formation *heating* the gas
c
      RETURN
      END

c-----------------------------------------------------------
      FUNCTION h2_destruction(ipart,xnH2,icosmic_ray)
c
c--Based on table B2 of Glover et al. (2010).
c     
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/physcon'

      INTEGER ipart
      REAL xnH2
      REAL h2_destruction
c
c--The passed logical "icosmic_ray" allows the code to turn off the
c     cosmic ray component of H_2 destruction when including heating
c     in the thermal energy equation, but turn it on for evolving the
c     abundance of H_2 itself.  This is because the H_2 destruction
c     heating by photodissociation includes an amount for the destruction
c     itself, but also a significant amount due to H_2 pumping from UV
c     photons that do not result in H_2 destruction but do heat the gas.
c     This extra pumping term is assumed to be proportional to the H_2
c     destruction rate (which should not include cosmic ray destruction).
c
c     In practice the cosmic ray destruction is usually very small, 
c     so it actually doesn't matter whether it is turned
c     on or off for the thermal behaviour.
c
      LOGICAL icosmic_ray

      REAL xnH
c
c--nH is twice nH2, where nH is actually the number density of protons from
c     hydrogen (i.e. it does not depend on the fractions of H and H_2 )
c
      xnH = 2.0*xnH2
c
c--Table B2 of Glover et al. (2010), and Section 2.2.
c     Includes cosmic ray dissociation.
c     Self-shielding and dust extinction is included in the 
c     photodissociation through the "heatingISR(4,i)" term
c     which is calculated via the TreeCol routine (needs H_2 
c     column density as well as dust extinction).
c
      IF (icosmic_ray) THEN
         h2_destruction = 1.2E-17
      ELSE
         h2_destruction = 0.
      ENDIF
      h2_destruction = h2_destruction + 5.6E-11*G0*heatingISR(4,ipart)
      h2_destruction = h2_destruction * xnH*h2frac(ipart)

      RETURN
      END

c-----------------------------------------------------------
      FUNCTION criticaln(ipart, gas_temp)
c
c--Based on Glover & MacLow (2007a)
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/interstellar'

      INTEGER ipart
      REAL gas_temp, criticaln
      REAL t4, crit_H, crit_H2

      t4 = gas_temp / 1.0E+4
      crit_H = 10**(3.0-0.461*t4 - 0.327*t4**2)
      crit_H2 = 10**(4.845 - 1.3*t4 + 1.62*t4**2)
      criticaln = 1.0/( (1.0-2.0*h2frac(ipart))/crit_H + 
     &     2*h2frac(ipart)/crit_H2 )

      RETURN
      END

c-----------------------------------------------------------
      FUNCTION carbon_chemistry(ipart, xnH2, gas_temp, depletion, x_CO)
c
c--Based on Keto & Caselli (2008)
c     
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rbnd'

      INTEGER ipart
      REAL carbon_chemistry,xnH2,gas_temp,depletion
      REAL x_CO_Cplus,x_C_Cplus,x_Cplus,x_C,x_CO
      REAL Vtherm,tau_on,tau_off

      IF (G0*heatingISR(3,ipart).LT.1.0E-10) THEN
         x_Cplus = 0.
         x_C = 23333./350000. * heatingISR(3,ipart)**0.6
         x_CO = 1.
      ELSE
         x_CO_Cplus = 1./23333. * xnH2 / (G0*heatingISR(3,ipart)**3.2)
         x_C_Cplus = 1./350000. * xnH2 / (G0*heatingISR(3,ipart)**2.6)

         x_Cplus = 1.0 / (1. + x_CO_Cplus + x_C_Cplus)
         x_C = x_C_Cplus * x_Cplus
         x_CO = x_CO_Cplus * x_Cplus
      ENDIF

      IF (iCHEM_depletion) THEN
c
c--Need to calculate depletion of CO from gas to grains
c
         Vtherm = SQRT(8.0*boltzmannk*gas_temp/(pi*30*mH))
c--Tau is in seconds
         tau_on = 1.0 / (1.0*4.0E-10*xnH2*3.4E-12*Vtherm)
c--The 1/3 is assuming a cosmic ray ionisation rate of 3x10^-17
         tau_off = 1.04E+14 * (1/3.)

         depletion = tau_on / (tau_on + tau_off)
      ELSE
         depletion = 1.0
      ENDIF
c
c--Set carbon fractions
c
      chemistry(1,ipart) = x_Cplus * metallicity
      chemistry(2,ipart) = x_C * metallicity
      chemistry(3,ipart) = x_CO * depletion * metallicity

      carbon_chemistry = x_Cplus
c      carbon_chemistry = 0.

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION electron_fraction(xnH,phiPAH)
c
c--Need to approximate electron density (approximation from Fig 10 of
c     Wolfire et al. 2003).  Note that this needs n(H) rather than
c     n(H2).
c
      REAL electron_fraction,xnH,phiPAH
      REAL xne_over_nH
c
c--Need to set phiPAH (Wolfire et al. use 0.5, but I use 0.55 to get
c     closer to 10^4 K at very low ISM densities (<1 cm^-3)
c
      phiPAH = 0.55

      xne_over_nH = 0.008/xnH
      IF (xne_over_nH.LT.1.0E-4) xne_over_nH = 1.0E-4
      IF (xne_over_nH.GT.1.) xne_over_nH = 1.
      electron_fraction = xne_over_nH*xnH

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION get_heatingISR(ext)
c
c--NOTE: This gives the heating rate assuming a standard ISRF (G0=1).
c     The scaling for G0!=1 is done elsewhere.
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/interstellar'

      REAL ext, xpos, heating, get_heatingISR
      INTEGER ipos

      IF (ext.LE.1.0E-3) THEN
         xpos = 1.
         ipos = 1
      ELSE
         xpos = LOG(ext/1.0E-3)/xISRlogfac + 1
         ipos = xpos + 1
      ENDIF
      IF (ipos.LT.nISRpoints) THEN
         heating = ISR_table(ipos) +
     &        (ISR_table(ipos+1)-ISR_table(ipos))*(xpos-ipos)
      ELSE
         heating = ISR_table(nISRpoints)
      ENDIF

      get_heatingISR = heating

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION get_Gphotoelectric(ext)
c
c--Returns the photoelectric heating for a given extinction value.
c     Uses interpolation from precalculated table of values.
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/interstellar'

      REAL ext, xpos, heating, get_Gphotoelectric
      INTEGER ipos

      IF (ext.LE.1.0E-3) THEN
         xpos = 1.
         ipos = 1
      ELSE
         xpos = LOG(ext/1.0E-3)/xISRlogfac + 1
         ipos = xpos + 1
      ENDIF
      IF (ipos.LT.nISRpoints) THEN
         heating = Ge_ISR_table(ipos) + 
     &        (Ge_ISR_table(ipos+1)-Ge_ISR_table(ipos))*(xpos-ipos)
      ELSE
         heating = Ge_ISR_table(nISRpoints)
      ENDIF

      get_Gphotoelectric = heating

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION gas_dust_collisional_term(xnH2,metallicity,gas_temp)
c
c--Returns the first bit of the gas-dust collisional heating/cooling 
c     term used by Keto & Field (2005) which is
c
c     Lambda_gd = 1.0E-33 * n_H2**2 * SQRT(T_gas) * (T_gas - T_dust)
c
c     or the coupling term used by Glover & Clark (2012) which is
c
c     Lambda_gd = 1.5E-32 * n_H2**2 * SQRT(T_gas) * (T_gas - T_dust) *
c                   (1.0 + 0.8*EXP(-75.0/T_gas)
c
c     in erg/cm^3/s.  This function also includes a metallicity term
c     which essentially assumes that the dust number density depends
c     linearly on the metallicity.
c
c     The value is calculated in cgs units and excludes the 
c     (T_gas-T_dust) term because of the need for derivatives w.r.t.
c     T_gas in solving the system of equations.
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/interstellar'

      REAL*8 gas_dust_collisional_term, xnH2, metallicity, gas_temp
      REAL*8 expfac

      IF (iDRT_gasdust.EQ.1) THEN
         gas_dust_collisional_term = 1.0E-33*(xnH2**2)*SQRT(gas_temp)*
     &              metallicity
      ELSEIF (iDRT_gasdust.EQ.2) THEN
         IF (gas_temp.GT.4.0) THEN
            expfac = 1.0-0.8*EXP(-75.0/gas_temp)
         ELSE
            expfac = 1.0
         ENDIF
         gas_dust_collisional_term = 1.5E-32*(xnH2**2)*SQRT(gas_temp)*
     &              expfac*metallicity
      ELSE
         gas_dust_collisional_term = 0.
      ENDIF

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION cosmic_ray_heating(xnH2)
c
c--Following Goldsmith et al. (2001), we set the cosmic ray heating
c     rate to be
c
c     Gamma_cr = 1.0E-27 * n_H2
c
c     in erg/cm^3/s.
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/interstellar'

      REAL*8 cosmic_ray_heating,xnH2

      IF (iDRT_cosmic_ray) THEN
         cosmic_ray_heating = 1.0E-27*xnH2
      ELSE
         cosmic_ray_heating = 0.
      ENDIF

      RETURN
      END

c-----------------------------------------------------------

      FUNCTION photoelectric_heating(ipart,gas_temp,xnH2,metallicity)
c
c--Following Young et al. (2004) we set the photoelectric heating rate
c     of the gas by electrons liberated from grains as
c     
c     Gamma_pe = 1E-24 * 0.05 * G(r) * n_H
c
c     and we take n_H = 2*n_H2 + 4*n_He = 2.66 nH2
c
c     The units are erg/cm^3/s and the function G(r) is the attenuation
c     factor for high energy photons (=1 for no extinction).
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/interstellar'

      REAL*8 photoelectric_heating,xnH2,gas_temp,metallicity
      REAL xne,phiPAH,xnH,Gr,ep,electron_fraction
      INTEGER ipart

      xnH = 2.0 * xnH2
      xne = electron_fraction(xnH,phiPAH)
      Gr = G0*heatingISR(2,ipart)
c      Gr = 1.

      ep = 0.049/(1.0+0.004*(Gr*SQRT(gas_temp)/xne/phiPAH)**0.73) +
     &     0.037*(gas_temp/1.E+4)**0.7/
     &     (1.0+2.0E-4*(Gr*SQRT(gas_temp)/xne/phiPAH))      
c      ep = 0.05

      IF (iDRT_photoelectric) THEN
         photoelectric_heating = 1.3E-24*ep*Gr
     &        *xnH  *metallicity
      ELSE
         photoelectric_heating = 0.
      ENDIF

      RETURN
      END

c-----------------------------------------------------------

