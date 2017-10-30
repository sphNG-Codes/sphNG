      FUNCTION stoppingtime(ipart)
c************************************************************
c                                                           *
c  One-fluid dust:  This function returns the stopping      *
c     time of particle ipart.                               *
c                                                           *
c************************************************************

      IMPLICIT NONE

      REAL stoppingtime
      INTEGER ipart

c      stoppingtime = 0.1
      stoppingtime = 0.0005

      RETURN
      END
