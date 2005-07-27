      SUBROUTINE scaling(gt, rs, drdt, dlnrdt)
c************************************************************
c                                                           *
c  Subroutine to compute the scaling factor and its various *
c     derivatives.                                          *
c                                                           *
c************************************************************

      INCLUDE 'COMMONS/expan'

      rs = 1.0 + vexpan*gt
      drdt = vexpan
      dlnrdt = drdt/rs

      RETURN
      END
