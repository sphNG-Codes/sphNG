      SUBROUTINE externf(ipart, xyzmh, fxyzu, iexf)
c************************************************************
c                                                           *
c  This subroutine computes the effect of an external       *
c     force.                                                *
c                                                           *
c        (1) Vertical gravitational field                   *
c        (2) -------                                        *
c        (3) Accretion disk                                 *
c        (4) Rotating cylinder                              *
c        (5) Central point mass                             *
c        (6) Distant point mass                             *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim), fxyzu(4,idim)

      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/xforce'
c
c--Unit angular momentum
c
      uang = udist**2/utime
c
c--Gravitational field
c
      IF (iexf.EQ.1) THEN
         grav = 1.0E4*utime**2/udist

         fxyzu(3,ipart) = fxyzu(3,ipart) - grav
c
c--Accretion disk
c
      ELSEIF (iexf.EQ.3) THEN
         t1 = 1.
         q = 2.
         h0 = 0.8
         r0 = h0**2
         r02 = r0**2
         omega = h0/r02
         d2 = xyzmh(1,ipart)**2 + xyzmh(2,ipart)**2
         d = SQRT(d2)
         d3 = d2*d
         omeg = omega*(r0/d)**q
         h1 = omeg*d2
         h2 = h1**2
         r2 = d2 + xyzmh(3,ipart)**2
         r = SQRT(r2)

         runix = xyzmh(1,ipart)/r
         runiy = xyzmh(2,ipart)/r
         runiz = xyzmh(3,ipart)/r
         dunix = xyzmh(1,ipart)/d
         duniy = xyzmh(2,ipart)/d

         fxyzu(1,ipart) = fxyzu(1,ipart) - 0.9999*runix/r2 + 
     &        t1*h2*dunix/d3
         fxyzu(2,ipart) = fxyzu(2,ipart) - 0.9999*runiy/r2 + 
     &        t1*h2*duniy/d3
         fxyzu(3,ipart) = fxyzu(3,ipart) - 0.9999*runiz/r2
c
c--Rotating cylinders
c
      ELSEIF (iexf.EQ.4) THEN
         omega = 0.6
         d2 = xyzmh(1,ipart)**2 + xyzmh(2,ipart)**2
         d = SQRT(d2)
         d3 = d2*d
         omeg = omega
         IF (d.GT.1.) omeg = 0.
         h1 = d2*omeg
         h2 = h1**2
         dunix = xyzmh(1,ipart)/d
         duniy = xyzmh(2,ipart)/d

         fxyzu(1,ipart) = fxyzu(1,ipart) + h2*dunix/d3
         fxyzu(2,ipart) = fxyzu(2,ipart) + h2*duniy/d3
c
c--Central point mass
c
      ELSEIF (iexf.EQ.5) THEN
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         d2 = (xi*xi + yi*yi + zi*zi + tiny)
         d = SQRT(d2)
         runix = xi/d
         runiy = yi/d
         runiz = zi/d

         fxyzu(1,ipart) = fxyzu(1,ipart) - xmass*runix/d2
         fxyzu(2,ipart) = fxyzu(2,ipart) - xmass*runiy/d2
         fxyzu(3,ipart) = fxyzu(3,ipart) - xmass*runiz/d2
         poten(ipart) = poten(ipart) - xmass/d
c
c--Distant point mass
c
      ELSEIF (iexf.EQ.6) THEN
         zdist = 40.
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = zdist - xyzmh(3,ipart)
         d2 = (xi*xi + yi*yi + zi*zi)
         d = SQRT(d2)
         runix = xi/d
         runiy = yi/d
         runiz = zi/d

         fxyzu(1,ipart) = fxyzu(1,ipart) - xmass*runix/d2
         fxyzu(2,ipart) = fxyzu(2,ipart) - xmass*runiy/d2
         fxyzu(3,ipart) = fxyzu(3,ipart) - xmass*runiz/d2
         poten(ipart) = poten(ipart) - xmass/d
      ENDIF

      RETURN

      END
