      SUBROUTINE sinksurface(ipart,xyzmh,rsurface,fx,fy,fz)
c************************************************************
c                                                           *
c  This subroutine computes the boundary force              *
c  for a sink particle.                                     *
c                                                           *
c************************************************************
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/ptmass'

      INTEGER ipart
      REAL xyzmh(5,idim)
      REAL rsurface, fx, fy, fz
      REAL xi, yi, zi, d, d2
      REAL runix, runiy, runiz
      REAL fsurface, ftotal

      xi = xyzmh(1,ipart) - xyzmh(1,listpm(1))
      yi = xyzmh(2,ipart) - xyzmh(2,listpm(1))
      zi = xyzmh(3,ipart) - xyzmh(3,listpm(1))

      d2 = (xi*xi + yi*yi + zi*zi + tiny)
      d = SQRT(d2)
      
      runix = xi/d
      runiy = yi/d
      runiz = zi/d

      IF (d.LE.(2.*rsurface)) THEN
         fsurface = (((2.*rsurface)-d)/
     &        (rsurface))**4
      ELSE
         fsurface = 0.0
      ENDIF
      
      ftotal = -fsurface*xyzmh(4,listpm(1))/d2

      fx = fx - ftotal*runix
      fy = fy - ftotal*runiy
      fz = fz - ftotal*runiz

      END SUBROUTINE sinksurface
