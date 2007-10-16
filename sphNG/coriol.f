      SUBROUTINE coriol(ipart, ntot, ti, xyzmh, vxyzu, fxyzu)
c************************************************************
c                                                           *
c  This routine adds centrifugal and coriolis forces        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim), vxyzu(4,idim), fxyzu(4,idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/typef'
c
c--Needed for MPI code
c
      IF (ipart.GT.ntot) THEN
         iparttree = ipart + ntot + 2
      ELSE
         iparttree = ipart
      ENDIF
c
c--Scaling factors
c
      CALL scaling(ti, rscale, drdt, dlnrdt)
      a3 = rscale**3
      a32 = 2.*omeg0*a3
      aadot = 2.*omeg0*drdt*rscale**2
      omeg2 = omeg0**2
      IF (ifcor.EQ.1) THEN
c
c--Rotation around z axis
c
c--Coriolis
c
         fcorcx = a32*vxyzu(2,ipart) + aadot*xyzmh(2,iparttree)
         fcorcy = -a32*vxyzu(1,ipart) - aadot*xyzmh(1,iparttree)
c
c--Centrifugal
c
         fcentx = a3*omeg2*xyzmh(1,iparttree)
         fcenty = a3*omeg2*xyzmh(2,iparttree)

         fxyzu(1,ipart) = fxyzu(1,ipart) + fcorcx + fcentx
         fxyzu(2,ipart) = fxyzu(2,ipart) + fcorcy + fcenty

         IF (idebug(1:6).EQ.'coriol') THEN
            WRITE (iprint, 99002) fxyzu(1,ipart)
            WRITE (iprint, 99002) fxyzu(2,ipart)
         ENDIF

      ELSEIF (ifcor.EQ.2) THEN
c    
c--Rotation end over end around x axis
c
c--Coriolis
c
         fcorcy = a32*vxyzu(3,ipart) + aadot*xyzmh(3,iparttree)
         fcorcz = -a32*vxyzu(2,ipart) - aadot*xyzmh(2,iparttree)
c
c--Centrifugal
c
         fcenty = a3*omeg2*xyzmh(2,iparttree)
         fcentz = a3*omeg2*xyzmh(3,iparttree)

         fxyzu(2,ipart) = fxyzu(2,ipart) + fcorcy + fcenty
         fxyzu(3,ipart) = fxyzu(3,ipart) + fcorcz + fcentz

         IF (idebug(1:6).EQ.'coriol') THEN
            WRITE (iprint, 99002) fxyzu(2,ipart)
            WRITE (iprint, 99002) fxyzu(3,ipart)
         ENDIF

      ENDIF

99002 FORMAT (1PE12.5)

      RETURN
      END

