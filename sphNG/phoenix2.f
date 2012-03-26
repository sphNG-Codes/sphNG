      SUBROUTINE phoenix2(i, idtsyn, itime, iphasein)
c***********************************************************************
c                                                                      *
c  Assigns an accreted particle or particle that goes outside the      *
c     boundary a new position, velocity etc to allow gas flow into a   *
c     box surrounding a planet to stablise.                            *
c                                                                      *
c***********************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/planetesimal'

      PARAMETER (numinjectmax = 10001)
      COMMON /zeustab/ xinjecttable(6,numinjectmax), dtheta, numinject

      REAL*4 rhoreal4
c
c--Assumes planet is at radius=1.  Also assumes disc
c     rotating anticlockwise and calculation done in rotating reference 
c     frame so planet is fixed.  Also assumes disc is thin.
c
      IF (iphasein.EQ.11) THEN
         CALL orbital_elements (i, eccentricity_p, inclination_p,
     &        1.0-variation, 1.0+variation)
c
c--Section injection, rotate position to be at phibound
c
         rxy2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
         radius = sqrt(rxy2 + xyzmh(3,i)**2)
         thetaadd = -phibound*0.999
         IF (radius.LT.1.0) thetaadd = phibound*0.999
         theta  = ATAN2(xyzmh(2,i), xyzmh(1,i)) + thetaadd
         xi = xyzmh(1,i)*cos(theta) + xyzmh(2,i)*sin(theta)
         yi = -xyzmh(1,i)*sin(theta) + xyzmh(2,i)*cos(theta)
c
c--Section injection, rotate velocities to phibound
c
         rr = sqrt(rxy2)
         vphi = (-xyzmh(2,i)*vxyzu(1,i) + vxyzu(2,i)*xyzmh(1,i))/rr
         vr = (vxyzu(1,i)*xyzmh(1,i) + vxyzu(2,i)*xyzmh(2,i))/rr
         omega = vphi/rr
         omega = omega - sqrt(xmass)
         vphi = omega*rr
         phi = ATAN(xyzmh(2,i)/xyzmh(1,i))

         xyzmh(1,i) = xi
         xyzmh(2,i) = yi
         vxyzu(1,i) = vr*cos(phi) - vphi*sin(phi)
         vxyzu(2,i) = vr*sin(phi) + vphi*cos(phi)

c         write(99,99000) xyzmh(1,i), xyzmh(2,i), xyzmh(3,i),
c     &        vxyzu(1,i),vxyzu(2,i),theta
c99000    FORMAT (6(1PE12.5,2X))
         xyzmh(4,i) = planetesimalmass
         goto 222
      ENDIF

c
c--Pick r and z from random standard accretion disc
c     iflag = 0
c--Pick from ZEUS
c
      iflag = 1
      
      IF (iflag.EQ.1 .AND. xmass.NE.1.0) THEN
         print *, 'Welcome to Phoenix2.f'
         print *, 'You have a star with a mass ne 1, and are using'
         print *, 'tables based upon ZEUS calculations.'
         print *, 'This code now uses xmass for the star mass,'
         print *, 'but the tables were made assuming a solar mass,'
         print *, 'so I am going to stop you to make sure you'
         print *, 'have engaged your brain. Comment this out to'
         print *, 'continue.'
         stop
      ENDIF

      IF (iflag.EQ.0) THEN
         radinject = findradius()
 10      zrandom = hoverr*gasdev(1)
         IF (zrandom.GT.zmax .OR. zrandom.LT.zmin) GOTO 10
         zinject = zrandom
      ELSE
c
c--Pick r and z from ZEUS tables of r and theta
c
 77      rand = ran1(1)
         ipos = numinject/2
         imove = ipos/2
c         write (*,*) rand
         DO j = 1, INT(LOG(1.0*numinject)/LOG(2.0))+1
c            write (*,*) ipos,xinjecttable(1,ipos),imove
            IF (rand.GT.xinjecttable(1,ipos)) THEN
               ipos = ipos+imove
            ELSE
               ipos = ipos-imove
            ENDIF
            ipos = MAX(1,MIN(numinject,ipos))
            imove = imove/2
         END DO
 110     IF (rand.GT.xinjecttable(1,ipos) .AND. ipos.LT.numinject) THEN
            ipos = ipos+1
c            write (*,*) 'Corr ',ipos,xinjecttable(1,ipos)
            GOTO 110
         ENDIF
 120     IF (ipos.GT.1) THEN
            IF (rand.LT.xinjecttable(1,ipos-1)) THEN
               ipos = MAX(1,ipos-1)
c               write (*,*) 'Corr ',ipos,xinjecttable(1,ipos)
               GOTO 120
            ENDIF
         ENDIF
c         WRITE (*,*) 'rand ',rand,xinjecttable(1,MAX(1,ipos-1)),
c     &        xinjecttable(1,ipos)

         IF (ipos.NE.numinject) THEN
            iposother = ipos+1
            dradius = ABS(ABS(xinjecttable(2,iposother))-
     &           ABS(xinjecttable(2,ipos)))
            IF (dradius.GT.0.1) THEN
               iposother = MAX(1,ipos-1)
               dradius = ABS(ABS(xinjecttable(2,ipos))-
     &              ABS(xinjecttable(2,iposother)))
            ENDIF
         ELSE
            iposother = MAX(1,ipos-1)
            dradius = ABS(ABS(xinjecttable(2,ipos))-
     &           ABS(xinjecttable(2,iposother)))
         ENDIF

         radinject = ABS(xinjecttable(2,ipos)) +
     &        (0.5*dradius) + (ran1(1)-0.5)*dradius

         IF (radinject.GT.(1.0+variation).OR.radinject.LT.
     &        (1.0-variation)) GOTO 77

c         radinject = MIN(ABS(xinjecttable(2,ipos)) +
c     &        (0.5*dradius) + (ran1(1)-0.5)*dradius,
c     &        1.0+variation)

c         radinject = MAX(ABS(xinjecttable(2,ipos)) +
c     &        (0.5*dradius) + (ran1(1)-0.5)*dradius,
c     &        1.0-variation)

         zinject = xinjecttable(3,ipos)+(ran1(1)-0.5)*radinject*dtheta

         IF (radinject**2+zinject**2.LT.(1.0-variation)**2.OR.
     &        radinject**2+zinject**2.GT.(1.0+variation)**2) GOTO 77

c     
c--Tables only occupy half of disc - need to populate above and below midplane
c
         IF (ran1(1).LT.0.5) zinject = -zinject
      ENDIF

c      WRITE (*,*) 'rad, z ',radinject,zinject
c
c--Set Cartesian coordinates
c
      IF (iflag.EQ.0) THEN
         IF (radinject.LT.1.0) THEN
            phiinject = -0.999*phibound
         ELSE
            phiinject = 0.999*phibound
         ENDIF
      ELSE
         IF (xinjecttable(2,ipos).LT.0.0) THEN
            phiinject = -0.999*phibound
         ELSE
            phiinject = 0.999*phibound
         ENDIF
      ENDIF
      cosangle = COS(phiinject)
      sinangle = SIN(phiinject)
      xyzmh(1,i) = radinject*cosangle
      xyzmh(2,i) = radinject*sinangle
      xyzmh(3,i) = zinject
c
c--Velocity equals Keplerian velocity (but in frame rotating with planet)
c
      IF (iflag.EQ.0) THEN
         velmag = xmass/SQRT(radinject)
         vxyzu(1,i) = - velmag*sinangle + xyzmh(2,i)
         vxyzu(2,i) = velmag*cosangle - xyzmh(1,i)
         vxyzu(3,i) = 0.0
      ELSE
c
c--Or velocity set from ZEUS tables of r and theta
c
         lastbin = 1
         if(iposother.LT.ipos) lastbin = -1

         rtop = abs(abs(radinject) - abs(xinjecttable(2,ipos)))
         rbot = abs(xinjecttable(2,iposother) - xinjecttable(2,ipos))

         dvelx = xinjecttable(4,iposother) - xinjecttable(4,ipos)
         vxyzu(1,i) = xinjecttable(4,ipos) + (rtop/rbot)*dvelx*lastbin

         dvely = xinjecttable(5,iposother) - xinjecttable(5,ipos)
         vxyzu(2,i) = xinjecttable(5,ipos) + (rtop/rbot)*dvely*lastbin

         dvelz = xinjecttable(6,iposother) - xinjecttable(6,ipos)
         vxyzu(3,i) = xinjecttable(6,ipos) + (rtop/rbot)*dvelz*lastbin
      ENDIF
c
c--Other quantities
c
      vxyzu(4,i) = hoverr**2/(radinject*gamma*(gamma-1.0))
      IF (encal.EQ.'r') THEN
         uradconst = radconst/uergcc
c
c--Assumes 75 g/cm^2 surface density at planet's radius, multiplied by signorm
c
         profile = 0.5
         rhoreal4 = signorm*75.0/(SQRT(2.0*pi)*hoverr*radinject*udist)/
     &        (radinject**profile)*EXP(-xyzmh(3,i)**2/
     &        (2.0*(hoverr*radinject)**2))/umass*udist**3

         ekcle(3,i) = getcv(rhoreal4,vxyzu(4,i))
         ekcle(1,i) = uradconst*(vxyzu(4,i)/ekcle(3,i))**4/rhoreal4
         ekcle(2,i) = getkappa(vxyzu(4,i),ekcle(3,i),rhoreal4)
         dumekcle(1,i) = ekcle(1,i)
         dumekcle(2,i) = ekcle(2,i)
         dumekcle(3,i) = ekcle(3,i)
      ENDIF

      xyzmh(4,i) = partm/gapfac
 222  IF (iplanetesimals.EQ.2) THEN
         rndiv = REAL(npart/2.)
      ELSE
         rndiv = REAL(npart)
      ENDIF
      xyzmh(5,i) = (4.0*hoverr*phibound*variation/(rndiv))**
     &     (1.0/3.0)

      DO j = 1, 5
         dumxyzmh(j,i) = xyzmh(j,i)
      END DO
      DO j = 1, 4
         dumvxyzu(j,i) = vxyzu(j,i)
         f1vxyzu(j,i) = 0.0
      END DO

      f1ha(1,i) = 0.0
      f1ha(2,i) = 0.0

      dgrav(i) = 0.0
      hasghost(i) = .FALSE.

      itry = istepmax/2
ccc      itry = istepmin
ccc      itry = imaxstep
 100  IF (MOD(idtsyn, itry).NE.0) THEN
         itry = itry/2
         GOTO 100
      ELSE
         isteps(i) = itry
      ENDIF

      it0(i) = itime
      it1(i) = it0(i) + isteps(i)/2
      it2(i) = it0(i) + isteps(i)
      iphase(i) = iphasein

      nreassign = nreassign + 1

      RETURN
      END

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
 1       v1=2.*ran1(idum)-1.
         v2=2.*ran1(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
      END
c
c--Find radius
c
      FUNCTION findradius()

      PARAMETER (nradtable = 10001)
      COMMON /cumradtab/ radiustable(nradtable)

      random = ran1(1)
      ipos = MIN(INT(random*(nradtable-1))+1, nradtable-1)
      x1 = (ipos-1)/REAL(nradtable-1)
      xoff = (random-x1)*(nradtable-1)
      findradius = radiustable(ipos)+xoff*(radiustable(ipos+1)-
     &     radiustable(ipos))

      RETURN
      END
c
c--Make table for finding radius
c
      SUBROUTINE cumradtable

      INCLUDE 'idim'

      PARAMETER (nradtable = 10001)
      COMMON /cumradtab/ radiustable(nradtable)
      INCLUDE 'COMMONS/rbnd'

      DIMENSION radtable(nradtable)
c
c--Set surface density profile
c
      profile = 0.5
c
c--Make table
c 
      pr1 = profile+1
      pr2 = profile+2
      pr3 = profile+3
      pr4 = profile+4
      rad = 1.0-variation
      radtable0 = -0.5*(rad**pr3)/pr3 + 0.5*(rad**pr2)/pr2 +
     &        0.375*(rad**pr4)/pr4 - 0.75*(rad**pr3)/pr3 +
     &        0.375*(rad**pr2)/pr2
      rad = 1.0
      radtable1 = -0.5*(rad**pr3)/pr3 + 0.5*(rad**pr2)/pr2 +
     &        0.375*(rad**pr4)/pr4 - 0.75*(rad**pr3)/pr3 +
     &        0.375*(rad**pr2)/pr2
      DO i = 1, nradtable
         rad = 1.0 + ((i-nradtable/2-1)*variation)/(nradtable/2)
         radtable(i) = -0.5*(rad**pr3)/pr3 + 0.5*(rad**pr2)/pr2 + 
     &        0.375*(rad**pr4)/pr4 - 0.75*(rad**pr3)/pr3 + 
     &        0.375*(rad**pr2)/pr2
         IF (rad.LT.1.0) THEN
            radtable(i) = radtable(i) - radtable0
         ELSE
            radtable(i) = radtable(nradtable/2)+(radtable1-radtable(i))
         ENDIF
      END DO
      DO i = 1, nradtable-1
         radtable(i) = radtable(i)/radtable(nradtable)
      END DO
      radtable(nradtable) = 1.0
c
c--Make new table equi-spaced in 'x' (random number between 0 and 1)
c
      DO j = 1, nradtable
         x = (j-1)/REAL(nradtable)
         DO i = 1, nradtable
            IF (x.LT.radtable(i)) THEN
               ipos = i
               GOTO 100
            ENDIF
         END DO
 100     ipos = MIN(ipos, nradtable-1)
         xoff = (x-radtable(ipos))/(radtable(ipos+1)-radtable(ipos))
         rad1 = 1.0 + ((ipos-nradtable/2-1)*variation)/(nradtable/2)
         rad2 = 1.0 + ((ipos+1-nradtable/2-1)*variation)/(nradtable/2)
         radiustable(j) = rad1 + xoff*(rad2-rad1)
      END DO

      RETURN
      END
c
c--Make table for finding injection information from ZEUS output
c
      SUBROUTINE zeustable

      INCLUDE 'idim'

      PARAMETER (numinjectmax = 10001)
      COMMON /zeustab/ xinjecttable(6,numinjectmax),dtheta,numinject
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/zeus'

      REAL muint,interp_func
c
c--Find which zeus dump applies to planet's Hill mass
c
      izeus = 1
      DO i = 6, 1,-1
         IF (hmass.GE.zmasses(i)) THEN
            izeus = i
            GOTO 200
         ENDIF
      END DO
 200  CONTINUE

      IF (izeus.NE.6) THEN
         muint = (hmass - zmasses(izeus))/
     &        (zmasses(izeus+1) - zmasses(izeus))
      ENDIF
c
c--Modify so that radii and thetas are zone-centred
c
      DO i=1,nradius-1
         drad(i) = radii(i+1)-radii(i)
         radiic(i) = (radii(i+1)+radii(i))/2.0
      END DO

      dtheta = thetas(2)-thetas(1)
      DO j=1,ntheta
         IF (j.EQ.ntheta) THEN
            thetasc(j) = (pi/2.0+thetas(j))/2.0
         ELSE
            thetasc(j) = (thetas(j+1)+thetas(j))/2.0
         ENDIF
      END DO
c
c--Need to produce cummulative table
c
      cum = 0.0
      numinject = 0
      DO k=1,2 
         DO j=1,ntheta
            DO i=1,nradius-1
c
c--Interpolate to find values for planet's current Hill mass
c
               IF (izeus.GE.5) THEN
                  densityint = extrap(0.0003,0.001,density(5,i,j,k),
     &                 density(6,i,j,k),(hmass-0.0003))
                  vrint = extrap(0.0003,0.001,vr(5,i,j,k),
     &                 vr(6,i,j,k),(hmass-0.0003))
                  vtint = extrap(0.0003,0.001,vt(5,i,j,k),
     &                 vt(6,i,j,k),(hmass-0.0003))
                  vpint = extrap(0.0003,0.001,vp(5,i,j,k),
     &                 vp(6,i,j,k),(hmass-0.0003))
               ELSEIF(izeus.LE.1) THEN
                  densityint = extrap(0.000003,0.00001,density(1,i,j,k),
     &                 density(2,i,j,k),(hmass-0.000003))
                  vrint = extrap(0.000003,0.00001,vr(1,i,j,k),
     &                 vr(2,i,j,k),(hmass-0.000003))
                  vtint = extrap(0.000003,0.00001,vt(1,i,j,k),
     &                 vt(2,i,j,k),(hmass-0.000003))
                  vpint = extrap(0.000003,0.00001,vp(1,i,j,k),
     &                 vp(2,i,j,k),(hmass-0.000003))

               ELSE
                  densityint = interp_func(log(density(izeus-1,i,j,k)),
     &                 log(density(izeus,i,j,k)),
     &                 log(density(izeus+1,i,j,k)),
     &                 log(density(izeus+2,i,j,k)),muint)

                  densityint = exp(densityint)

                  vrint = interp_func(vr(izeus-1,i,j,k),
     &                 vr(izeus,i,j,k),vr(izeus+1,i,j,k),
     &                 vr(izeus+2,i,j,k),muint)

                  vtint = interp_func(vt(izeus-1,i,j,k),
     &                 vt(izeus,i,j,k),vt(izeus+1,i,j,k),
     &                 vt(izeus+2,i,j,k),muint)

                  vpint = interp_func(vp(izeus-1,i,j,k),
     &                 vp(izeus,i,j,k),vp(izeus+1,i,j,k),
     &                 vp(izeus+2,i,j,k),muint)
               ENDIF

               IF (vpint*(k-1.5).LT.0.0) THEN
                  numinject = numinject + 1
                  IF (numinject.GT.numinjectmax) THEN
                     WRITE (*,*) 'numinjectmax too small'
                     STOP
                  ENDIF

                  cum = cum + densityint*radiic(i)*drad(i)*
     &                 ABS(vpint)
                  IF(densityint.LT.0.0) print *, 'dens',densityint
                  IF(radiic(i).LT.0.0) print *, 'rad', radiic(i)
                  IF(drad(i).LT.0.0) print *, 'dr',drad(i),i
                  sintheta = sin(thetasc(j))
                  costheta = cos(thetasc(j))
                  xinjecttable(1,numinject) = cum
c
c--xinjecttable(2,*) is negative for k=1 (inner radius)
c
              xinjecttable(2,numinject)=2.0*(k-1.5)*radiic(i)*sintheta

                  xinjecttable(3,numinject) = radiic(i)*costheta
                  IF (k-1.5 .LT. 0.0) THEN
                     phi = -phibound
                  ELSE
                     phi = phibound
                  ENDIF
                  sinphi = sin(phi)
                  cosphi = cos(phi)

                  xinjecttable(4,numinject) = vrint*cosphi*sintheta
     &         - vpint*sinphi*sintheta + vtint*costheta*cosphi
                  xinjecttable(5,numinject) = vrint*sinphi*sintheta
     &         + vpint*cosphi*sintheta + vtint*costheta*sinphi
                  xinjecttable(6,numinject) = vrint*costheta - 
     &                 vtint*sintheta

               ENDIF
            END DO
         END DO
      END DO
c
c--Normalise cummulative probability
c
      DO i=1,numinject
         xinjecttable(1,i)=xinjecttable(1,i)/cum
         write (33,*) xinjecttable(1,i)
      END DO

      WRITE (*,*) 'Number inject points ',numinject

      RETURN
      END
c
c--Interpolation function
c
      FUNCTION interp_func(y0,y1,y2,y3,mu)

      REAL y0,y1,y2,y3,mu,a0,a1,a2,a3,interp_func,mu2

      interp_func = 0.0

      a0 = y3 - y2 - y0 + y1
      a1 = y0 - y1 - a0
      a2 = y2 - y0
      a3 = y1

      mu2 = mu*mu

      interp_func = a0*mu*mu2+a1*mu2+a2*mu+a3

      END FUNCTION interp_func
c
c--Extrapolation function
c
      FUNCTION extrap(x0,x1,y0,y1,mu)

      REAL x0,x1,y0,y1,mu

      extrap = 0.0
      extrap = ((y1-y0)/(x1-x0))*mu + y0

      END FUNCTION extrap
