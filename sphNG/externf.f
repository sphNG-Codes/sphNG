      SUBROUTINE externf(ipart, ti, xyzmh, fxyzu, rho, iexf)
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
c        (7) Central point mass and planet                  *
c        (8) Galactic spiral potential (C. Dobbs)           *
c        (9) Tokamak potential (D. Price)                   *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim), fxyzu(4,idim)
      REAL*4 rho(idim)

      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/potent'
      INCLUDE 'COMMONS/tokamak'
c
c--Unit angular momentum
c
      uang = udist**2/utime

c   z potential constant
       zq=0.7

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
c
c--Central point mass and planet
c
      ELSEIF (iexf.EQ.7) THEN
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
c--Assumes planet at location (1,0,0) in code units
c
         xi = xi - 1.0
c
c--Radius of earth is 8.2e-6 in 5.2au code units
c--Radius of 30 M_e solid core would be 2.5e-5 in 5.2 au code units
c
         rplanet = 2.5E-05
         rplanet = 1.0E-03
         d2 = (xi*xi + yi*yi + zi*zi + (rplanet)**2)
         d = SQRT(d2)
         runix = xi/d
         runiy = yi/d
         runiz = zi/d

         fxyzu(1,ipart) = fxyzu(1,ipart) - planetmass*runix/d2
         fxyzu(2,ipart) = fxyzu(2,ipart) - planetmass*runiy/d2
         fxyzu(3,ipart) = fxyzu(3,ipart) - planetmass*runiz/d2
         poten(ipart) = poten(ipart) - planetmass/d
      ELSEIF (iexf.EQ.8) THEN
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         hi = xyzmh(5,ipart)
         dhi = hi/1000.
         d2 = (xi*xi + yi*yi)
         r=SQRT(d2+zi*zi)

c contribution of logarithmic potential and halo

         fxyzu(1,ipart)=fxyzu(1,ipart)-2.*Co*xi/(Rc**2.+d2+(zi/zq)**2.)
     &          -p1*(rc2/r)**2.*(r-rc2*atan(r/rc2))*xi/r

         fxyzu(2,ipart)=fxyzu(2,ipart)-2.*Co*yi/(Rc**2.+d2+(zi/zq)**2.)
     &          -p1*(rc2/r)**2.*(r-rc2*atan(r/rc2))*yi/r

         fxyzu(3,ipart)=fxyzu(3,ipart)-2.*Co*zi/((Rc**2.+d2+(zi/zq)**2.)
     &             *zq**2.)

c spiral perturbation

         call potential(xi,yi,zi,ti,potent1)

         call potential(xi+dhi,yi,zi,ti,potent2)
         fxyzu(1,ipart) = fxyzu(1,ipart) - (potent2-potent1)/dhi

         call potential(xi,yi+dhi,zi,ti,potent2)
         fxyzu(2,ipart) = fxyzu(2,ipart) - (potent2-potent1)/dhi

         call potential(xi,yi,zi+dhi,ti,potent2)
         fxyzu(3,ipart) = fxyzu(3,ipart) - (potent2-potent1)/dhi
      ELSEIF (iexf.EQ.9) THEN
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
c get cylindrical r
         rcyl = sqrt(xi**2 + yi**2)
         IF (rcyl.GT.tiny) THEN
            drcyl = 1./rcyl
         ELSE
            drcyl = 0.
         ENDIF
c rintorus is radius from centre of torus
         rintorus2 = (rcyl - Rtorus)**2 + zi**2
         rintorus = SQRT(rintorus2)
         IF (rintorus.GT.tiny) THEN
            drintorus = 1./rintorus
         ELSE
            drintorus = 0.
         ENDIF
cc         IF (rintorus2.LT.atorus**2) THEN
c current in torus "z" direction
            term = 1. - rintorus2*da2
            currjz = currj0*term**nutorus
c Bfield in torus "theta" direction
            Btheta = currj0*atorus**2/(2.*(nutorus+1))*
     &           (1. - term**(nutorus+1))*drintorus
c force is in torus "r" direction  (J X B)
            frtorus = -Btheta*currjz
cc         ELSE
ccc outside torus, use a f proportional to distance
cc            frtorus = -(rintorus - atorus)
cc         ENDIF
c get force in cylindrical co-ordinates
         frcyl = frtorus*(rcyl - Rtorus)*drintorus
         fz = frtorus*zi*drintorus
c translate to cartesians
         fxyzu(1,ipart) = fxyzu(1,ipart) + frcyl*xi*drcyl/rho(ipart)
         fxyzu(2,ipart) = fxyzu(2,ipart) + frcyl*yi*drcyl/rho(ipart)
         fxyzu(3,ipart) = fxyzu(3,ipart) + fz/rho(ipart)
      ENDIF

      RETURN

      END
