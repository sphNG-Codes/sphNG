      REAL FUNCTION thermal_velocity(ipart)
c
c--Returns thermal speed of gas particles for use with dust/gas mixtures
c     In code units.
c
      IMPLICIT NONE

      INCLUDE 'idim'      

      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/cgas'

      INTEGER ipart

      thermal_velocity = vthermal_coef*vsound(ipart)

      RETURN

      END FUNCTION thermal_velocity
