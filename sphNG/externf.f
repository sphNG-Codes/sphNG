      SUBROUTINE externf(ipart,ntot,ti,xyzmh,fxyzu,rho,iexf)
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

      DIMENSION xyzmh(5,idim), fxyzu(4,idim)
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
         fxyzu(1,ipart) = fxyzu(1,ipart) + fextx
         fxyzu(2,ipart) = fxyzu(2,ipart) + fexty
         fxyzu(3,ipart) = fxyzu(3,ipart) + fextz
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
         fru = sqrt(fxyzu(1,ipart)*fxyzu(1,ipart)+
     &             fxyzu(2,ipart)*fxyzu(2,ipart))
c         if (rxyplane.ge.0.9) print *, frcyl
c         if (abs(frcyl).gt.0.) print *,frcyl/fru,frcyl,radius-rxyplane
         fxyzu(1,ipart) = fxyzu(1,ipart)+frcyl*xi*drxyplane
         fxyzu(2,ipart) = fxyzu(2,ipart)+frcyl*yi*drxyplane  
         
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
