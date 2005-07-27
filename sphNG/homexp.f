      SUBROUTINE homexp(ipart, ti, vxyzu, fxyzu)
c************************************************************
c                                                           *
c  This subroutine computes the correction to the momentum  *
c     equation due to homologous expansion or contraction.  *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/logun'

      DIMENSION vxyzu(4,idim), fxyzu(3,idim)
c
c--Scaling factors
c
      CALL scaling(ti, rscale, drdt, dlnrdt)

      dlnrdt2 = 2.*dlnrdt
      rscale3 = 1./rscale**3
      fxyzu(1,ipart) = rscale3*fxyzu(1,ipart) - vxyzu(1,ipart)*dlnrdt2
      fxyzu(2,ipart) = rscale3*fxyzu(2,ipart) - vxyzu(2,ipart)*dlnrdt2
      fxyzu(3,ipart) = rscale3*fxyzu(3,ipart) - vxyzu(3,ipart)*dlnrdt2

      RETURN
      END
