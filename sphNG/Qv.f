      FUNCTION Qv(v,Qabs)
c
c--This returns the total opacity (scattering + absorption) of the dust 
c     with frequency.  It also gives Qabs, which is the absorption 
c     efficiency of the dust with frequency.
c     If MOD(iopmodel,10).LE.1 the two are the same, but if 
c     MOD(iopmodel,10).EQ.2 they are not.
c
c--For MOD(iopmodel,10)=0 or 1 is based on Zucconi et al. (2001)
c     For MOD(iopmodel,10)=2, a table is read from file/created in 
c     preset.F of opacity vs frequency
c
c--NOTE: The values are in cgs units, not code units.
c     The units are cm^2/H_2 and it assumes a standard 100:1 gas to dust
c     ratio.  Therefore, to allow different metallicities, we multiply
c     the result by metallicity which is in solar units.
c
      IMPLICIT NONE

      REAL*8 v,Qv,Qabs
      REAL*8 value,xlog10_lambda
      INTEGER ipos

      INCLUDE 'idim'

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/optbl'

      IF (MOD(iopmodel,10).LE.1) THEN
         IF (c/v .LT. 1e-3) THEN
            value = 3.9E-22 * (1.0E-4/(c/v))**1.4
         ELSEIF (c/v .GT. 0.04) THEN
            value = 3.3E-26 * (1.06E-1/(c/v))**2.0
         ELSE
            value = 1.5E-24 * (1.4E-2/(c/v))**1.6
         ENDIF

         Qabs = value

      ELSEIF (MOD(iopmodel,10).EQ.2) THEN
         xlog10_lambda = LOG10(c/v*1.0E+04)
         ipos = MIN(nopacity_lambda-1.0,
     &        MAX(1.0,(xlog10_lambda+1.0)/0.02885581 + 1.))
         xlog10_lambda = xlog10_lambda - ((ipos-1)*0.02885581 - 1.0)

         value = xlog10_lambda*(opacity_vs_lambda(2,ipos+1) - 
     &        opacity_vs_lambda(2,ipos))/0.02885581 + 
     &        opacity_vs_lambda(2,ipos)
         Qabs = xlog10_lambda*(opacity_vs_lambda(3,ipos+1) - 
     &        opacity_vs_lambda(3,ipos))/0.02885581 + 
     &        opacity_vs_lambda(3,ipos)

         value = 10**value  *(gmw*mH)
         Qabs = 10**Qabs *(gmw*mH)

      ELSE
         WRITE (*,*) 'ERROR - opacity model not recognised ',iopmodel
         CALL quit(0)
      ENDIF

      Qabs = Qabs * metallicity
      Qv = value * metallicity
c
c--For fixed value
c
c      Qabs = 2.381*mH  * 2.0 * 3.0E+11
c      Qv = Qabs

      RETURN
      END

c-----------------------------------------------------------

