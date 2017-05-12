      FUNCTION planck(v,T)
c
c--Planck function (in cgs units)
c
      IMPLICIT NONE

      INCLUDE 'COMMONS/physcon'

      REAL*8 planck, v, T, hoverk, val2, val3

      hoverk = planckconst/boltzmannk

      val3 = hoverk*v/T
      IF (val3.LT.600.) THEN
         val2 = EXP(val3) - 1.0
         IF (val2.GT.0.) THEN
            planck = (2.0*((planckconst*v)*v)*v/c**2) / val2
         ELSE
            planck = 0.
         ENDIF
      ELSE
         planck = 0.
      ENDIF

      RETURN
      END

c-----------------------------------------------------------
