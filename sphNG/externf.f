      SUBROUTINE externf(ipart,ntot,ti,xyzmh,vxyzu,fx,fy,fz,rho,iexf,
     &     ibound, irotpot)
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
c        (9) potential from external B field                *
c           (D. Price/C.Toniolo)                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim),vxyzu(4,idim)
      REAL rhoi,fextx,fexty,fextz, fbound, rxyplane, drxyplane
      REAL*4 rho(idim)

      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/potent'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/cgas'
c      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/cylinder'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/prdrag'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/pxpy'
c
c--Needed for MPI code
c
      IF (ipart.GT.ntot) RETURN
c
c--Unit angular momentum
c
      uang = udist**2/utime
c
c--z potential constant
c
      zq=0.7
c
c--Gravitational field
c
      IF (iexf.EQ.1) THEN
         grav = 1.0E4*utime**2/udist

         fz = fz - grav
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

         fx = fx - 0.9999*runix/r2 + 
     &        t1*h2*dunix/d3
         fy = fy - 0.9999*runiy/r2 + 
     &        t1*h2*duniy/d3
         fz = fz - 0.9999*runiz/r2
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

         fx = fx + h2*dunix/d3
         fy = fy + h2*duniy/d3
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

         fx = fx - xmass*runix/d2
         fy = fy - xmass*runiy/d2
         fz = fz - xmass*runiz/d2
         poten(ipart) = poten(ipart) - xmass/d
c
c--Forces related to radiation pressure and PR drag
c
         IF (prswitch) THEN
            IF (iphase(ipart).EQ.11) THEN
               clight = c*utime/udist
               vx = vxyzu(1,ipart)
               vy = vxyzu(2,ipart)
               vz = vxyzu(3,ipart)
               vr = vx*runix + vy*runiy + vz*runiz
               
               prcoeff = xmass*pbeta/d2
               fx = fx + ((1.-vr/clight)*runix - vx/clight)*prcoeff
               fy = fy + ((1.-vr/clight)*runiy - vy/clight)*prcoeff
               fz = fz + ((1.-vr/clight)*runiz - vz/clight)*prcoeff
            ENDIF
         ENDIF
c
c--Inner disc edge boundary which acts only on gas
c
         IF (ibound.EQ.102 .AND. iphase(ipart).EQ.0) THEN
            range = 0.005
            rlmin = rmind - (2.*range)
            
            IF (d-rlmin.LE.2.*range) THEN
               fsurface = (((2.*range) + rlmin - d)/(range))**4
            ELSE
               fsurface = 0.0
            ENDIF
            
            fx = fx + fsurface*xmass*runix/d2
            fy = fy + fsurface*xmass*runiy/d2
            fz = fz + fsurface*xmass*runiz/d2
         ENDIF

c
c--Wyatt style migration, now useable for sinks and potentials.
c  Odd use of runiy, runix due to force acting in phi direction.
c
         IF (imigrate.EQ.1 .AND. iphase(ipart).GE.1 .AND.
     &        iphase(ipart).LT.10 .AND. ti.LT.rorbitmax) THEN
            fwyatt = 0.5*pmrate*sqrt(xmass/d**3)
            fx = fx - fwyatt*runiy
            fy = fy + fwyatt*runix
         ENDIF
c
c--Orbiting potential case, where gravity and surface forces are
c  applied here.
c
         IF (irotpot.EQ.1) THEN
            CALL planetpotential (px, py, pz, imigrate, rorbitmax,
     &           pmrate, rorbit_orig, ti)

            xi = xi - px
            yi = yi - py

            d2 = (xi*xi + yi*yi + zi*zi + tiny)
            d = SQRT(d2)

            runix = xi/d
            runiy = yi/d
            runiz = zi/d

            IF (d.LE.(2.*rplanet*pradfac(1)) .AND.
     &           iphase(ipart).EQ.0) THEN
               fsurface = (((2.*rplanet*pradfac(1))-d)/
     &              (rplanet*pradfac(1)))**4
            ELSE
               fsurface = 0.0
            ENDIF
            
            ftotal = planetmass/d2*(1.0 - fsurface)

            fx = fx - ftotal*runix
            fy = fy - ftotal*runiy
            fz = fz - ftotal*runiz
            poten(ipart) = poten(ipart) - planetmass/d
         ENDIF

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

         fx = fx - xmass*runix/d2
         fy = fy - xmass*runiy/d2
         fz = fz - xmass*runiz/d2
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

         fx = fx - xmass*runix/d2
         fy = fy - xmass*runiy/d2
         fz = fz - xmass*runiz/d2
         poten(ipart) = poten(ipart) - xmass/d
c
c--Assumes planet at location (1,0,0) in code units
c
         xi = xi - 1.0
c
c--Radius of earth is 8.2e-6 in 5.2au code units
c--Radius of 30 M_e solid core would be 2.5e-5 in 5.2 au code units
c
c         rplanet = 2.5E-05
c         rplanet = 1.0E-03
c         d2 = (xi*xi + yi*yi + zi*zi + (rplanet)**2)
c         d = SQRT(d2)

         d2 = (xi*xi + yi*yi + zi*zi + tiny)
         d = SQRT(d2)
         rp = d - (rplanet*pradfac(1))

         runix = xi/d
         runiy = yi/d
         runiz = zi/d
c
c--Forces relating to planet's surface
c
         IF (d.LE.(2.*rplanet*pradfac(1))) THEN
            fsurface = (((2.*rplanet*pradfac(1))-d)/
     &           (rplanet*pradfac(1)))**4
         ELSE
            fsurface = 0.0
         ENDIF

         ftotal = planetmass/d2*(1.0 - fsurface)

         fx = fx - ftotal*runix
         fy = fy - ftotal*runiy
         fz = fz - ftotal*runiz
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

         fx=fx-2.*Co*xi/(Rc**2.+d2+(zi/zq)**2.)
     &          -p1*(rc2/r)**2.*(r-rc2*atan(r/rc2))*xi/r

         fy=fy-2.*Co*yi/(Rc**2.+d2+(zi/zq)**2.)
     &          -p1*(rc2/r)**2.*(r-rc2*atan(r/rc2))*yi/r

         fz=fz-2.*Co*zi/((Rc**2.+d2+(zi/zq)**2.)
     &             *zq**2.)

c spiral perturbation

         call potential(xi,yi,zi,ti,potent1)

         call potential(xi+dhi,yi,zi,ti,potent2)
         fx = fx - (potent2-potent1)/dhi

         call potential(xi,yi+dhi,zi,ti,potent2)
         fy = fy - (potent2-potent1)/dhi

         call potential(xi,yi,zi+dhi,ti,potent2)
         fz = fz - (potent2-potent1)/dhi
      ELSEIF (iexf.EQ.9) THEN
c
c--external force due to an assumed external B field
c  (with non-zero curl)
c         
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
         hi = xyzmh(5,ipart)
         rhoi = rho(ipart)
         CALL fexternalB(xi,yi,zi,hi,rhoi,fextx,fexty,fextz)
         fx = fx + fextx
         fy = fy + fexty
         fz = fz + fextz
c
c External force to re-inject particles that try to escape from the cylinder
c considered as a solid boundary
c
      ELSEIF (iexf.EQ.10) THEN
         xi = xyzmh(1,ipart)
         yi = xyzmh(2,ipart)
         zi = xyzmh(3,ipart)
c get cylindrical r
         rxyplane = sqrt(xi**2 + yi**2)
         v2 = gamma*2./3.*RK2*(rhozero)**(gamma-1.)
c         v2 = RK2
c Correct estimate of velocity with Alfven speed when B field on   
         IF (imhd.EQ.idim) THEN
            valfven2 = ampl*ampl/rhozero
            v2 = v2 + valfven2
         ENDIF   
c         if (rxyplane.ge.0.9) print *, v2, ampl*ampl/rhozero, rxyplane, 
c     &    radius,hzero
         frcyl = fbound(rxyplane,radius,hzero,10.*v2)
         IF (rxyplane.GT.tiny) THEN
            drxyplane = 1./rxyplane
         ELSE
            drxyplane = 0.
         ENDIF
c get force in cartesian co-ordinates
         fru = sqrt(fx*fx+
     &             fy*fy)
c         if (rxyplane.ge.0.9) print *, frcyl
c         if (abs(frcyl).gt.0.) print *,frcyl/fru,frcyl,radius-rxyplane
         fx = fx+frcyl*xi*drxyplane
         fy = fy+frcyl*yi*drxyplane  
         
      ENDIF

      RETURN

      END

      FUNCTION fbound(rr,rmax,hi,v2)
c************************************************************
c                                                           *
c  This subroutine computes the boundary force              *
c  as in Monaghan (2005)                                    *
c                                                           *
c************************************************************
      IMPLICIT NONE
      REAL fbound
      REAL rr,rmax,hi,v2
      REAL gamfac,yy,qfac
      REAL cnormk,tiny
      PARAMETER (cnormk = 1.) !/3.1415926536)
      PARAMETER (tiny = 1.e-14)

      IF (hi.LE.0.) THEN
         WRITE (*,*) 'Stop, hzero.LE.0 in fbound'
         CALL quit
      ENDIF
      yy = abs(rr-rmax)
      qfac = yy/hi
      IF (qfac.LT.2./3.) THEN
         gamfac = cnormk*2./3.
      ELSEIF (qfac.LT.1.) THEN
         gamfac = cnormk*(2.*qfac - 1.5*qfac**2)
      ELSEIF (qfac.LT.2.) THEN
         gamfac = cnormk*(0.5*(2.-qfac)**2)
      ELSE
         gamfac = 0.
      ENDIF
   
      fbound = -0.1*v2*gamfac/(yy+tiny)
      IF (rr.gt.rmax) THEN
         WRITE(*,*) 'ERROR! particle crossed boundary'
         WRITE(*,*) 'rr = ',rr,rmax,v2,fbound
         CALL quit
      ENDIF
      
      END FUNCTION fbound
