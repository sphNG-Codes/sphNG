      SUBROUTINE homexp(ipart, ti, vxyzu, fx, fy, fz)
c************************************************************
c                                                           *
c  This subroutine computes the correction to the momentum  *
c     equation due to homologous expansion or contraction.  *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/logun'

      DIMENSION vxyzu(4,idim)
c
c--Scaling factors
c
      CALL scaling(ti, rscale, drdt, dlnrdt)

      dlnrdt2 = 2.*dlnrdt
      rscale3 = 1./rscale**3
      fx = rscale3*fx - vxyzu(1,ipart)*dlnrdt2
      fy = rscale3*fy - vxyzu(2,ipart)*dlnrdt2
      fz = rscale3*fz - vxyzu(3,ipart)*dlnrdt2

      RETURN
      END
