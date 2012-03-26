      subroutine sphericalpos(xi,yi,zi,radius,phi,theta)

c************************************************************
c                                                           *
c  Routine to find a position in spherical coordinates      *
c  about a planet, keeping the angles relative to the line  *
c  connecting the star and planet.                          *
c                                                           *
c  Star-Planet line at +-pi on planet's surface.            *
c                                                           *
c  Ben Ayliffe, 17th November 2011.                         *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/pxpy'
      
      REAL xi, yi, zi
      REAL radius, phi, theta
      REAL dx, dy, dz
      REAL xr, yr, dxr, dyr
      REAL radius_p, phi_p, theta_p

      radius_p = sqrt(px*px + py*py + pz*pz)
      phi_p = ATAN2(py,px)
      theta_p = ACOS(pz/radius_p)

      dx = xi - px
      dy = yi - py
      dz = zi - pz

      radius = sqrt(dx*dx + dy*dy + dz*dz)

      xr =  xi*cos(phi_p) + yi*sin(phi_p)
      yr = -xi*sin(phi_p) + yi*cos(phi_p)
      dxr = xr - sqrt(px*px + py*py)
      dyr = yr

      phi = ATAN2(dyr,dxr)

      theta = ACOS(dz/radius)

      end subroutine sphericalpos
