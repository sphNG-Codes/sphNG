      SUBROUTINE phoenix2(i, idtsyn, itime)
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

      PARAMETER (numinjectmax = 10001)
      COMMON /zeustab/ xinjecttable(6,numinjectmax), dtheta, numinject
c
c--Assumes planet is at radius=1 and stellar mass is 1.  Also assumes disc
c     rotating anticlockwise and calculation done in rotating reference 
c     frame so planet is fixed.  Also assumes disc is thin.
c
c
c--Pick r and z from random standard accretion disc
c     iflag = 0
c--Pick from ZEUS
c     iflag = 1
c
      iflag = 1

      IF (iflag.EQ.0) THEN
         radinject = findradius()
 10      zrandom = hoverr*gasdev(1)
         IF (zrandom.GT.zmax .OR. zrandom.LT.zmin) GOTO 10
         zinject = zrandom
      ELSE
c
c--Pick r and z from ZEUS tables of r and theta
c
         rand = ran1(1)
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
         dradius = ABS(ABS(xinjecttable(2,ipos+1))-
     &        ABS(xinjecttable(2,ipos)))
         IF (dradius.GT.0.1) dradius = ABS(ABS(xinjecttable(2,ipos))-
     &        ABS(xinjecttable(2,MAX(1,ipos-1))))
         IF (dradius.GT.0.1) dradius = 0.0
         
         radinject = ABS(xinjecttable(2,ipos)) + (ran1(1)-0.5)*dradius
         zinject = xinjecttable(3,ipos)+(ran1(1)-0.5)*radinject*dtheta
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
            phiinject = -0.98*phibound
         ELSE
            phiinject = 0.98*phibound
         ENDIF
      ELSE
         IF (xinjecttable(2,ipos).LT.0.0) THEN
            phiinject = -0.98*phibound
         ELSE
            phiinject = 0.98*phibound
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
         velmag = 1.0/SQRT(radinject)
         vxyzu(1,i) = - velmag*sinangle + xyzmh(2,i)
         vxyzu(2,i) = velmag*cosangle - xyzmh(1,i)
         vxyzu(3,i) = 0.0
      ELSE
c
c--Or velocity set from ZEUS tables of r and theta
c
         vxyzu(1,i) = xinjecttable(4,ipos)
         vxyzu(2,i) = xinjecttable(5,ipos)
         vxyzu(3,i) = xinjecttable(6,ipos)
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
         rho = signorm*75.0/(SQRT(2.0*pi)*hoverr*radinject*udist)/
     &        (radinject**profile)*EXP(-xyzmh(3,i)**2/
     &        (2.0*(hoverr*radinject)**2))/umass*udist**3
         ekcle(3,i) = getcv(rho,vxyzu(4,i))
         ekcle(1,i) = uradconst*(vxyzu(4,i)/ekcle(3,i))**4/rho
         ekcle(2,i) = getkappa(vxyzu(4,i),ekcle(3,i),rho)
      ENDIF

      xyzmh(5,i) = (4.0*hoverr*phibound*variation/(FLOAT(npart)))**
     &     (1.0/3.0)
      xyzmh(4,i) = partm

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
      iphase(i) = 0

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

      PARAMETER (nradiusmax=100)
      PARAMETER (nthetamax=36)
      DIMENSION radii(nradiusmax),thetas(nthetamax),drad(nradiusmax)
      DIMENSION density(nradiusmax,nthetamax,2),vr(nradiusmax,
     &     nthetamax,2),vt(nradiusmax,nthetamax,2),
     &     vp(nradiusmax,nthetamax,2)
c
c--Read in from ZEUS/IDL file
c
      OPEN (44,FILE='outputgrid.dat')
      READ (44,*) nradius, ntheta
      IF (nradius.GT.nradiusmax .OR. ntheta.GT.nthetamax) THEN
         WRITE (*,*) 'Grid too big in outputgrid.dat'
         STOP
      ENDIF
      READ (44,*) (radii(i), i=1, nradius)
      READ (44,*) (thetas(i), i=1, ntheta)
      READ (44,*) ((density(i,j,1), j=1, ntheta), i=1, nradius)
      READ (44,*) ((vr(i,j,1), j=1, ntheta), i=1, nradius)
      READ (44,*) ((vt(i,j,1), j=1, ntheta), i=1, nradius)
      READ (44,*) ((vp(i,j,1), j=1, ntheta), i=1, nradius)
      READ (44,*) ((density(i,j,2), j=1, ntheta), i=1, nradius)
      READ (44,*) ((vr(i,j,2), j=1, ntheta), i=1, nradius)
      READ (44,*) ((vt(i,j,2), j=1, ntheta), i=1, nradius)
      READ (44,*) ((vp(i,j,2), j=1, ntheta), i=1, nradius)
      CLOSE (44)
c
c--Modify so that radii and thetas are zone-centred
c
      DO i=1,nradius-1
         drad(i) = radii(i+1)-radii(i)
         radii(i) = (radii(i+1)+radii(i))/2.0
c         write (44,*) radii(i),vp(i,36,1),vp(i,36,2)
      END DO
      nradius = nradius - 1 
      dtheta = thetas(2)-thetas(1)
      DO j=1,ntheta
         IF (j.EQ.ntheta) THEN
            thetas(j) = (pi/2.0+thetas(j))/2.0
         ELSE
            thetas(j) = (thetas(j+1)+thetas(j))/2.0
         ENDIF
      END DO
c
c--Need to produce cummulative table
c
      cum = 0.0
      numinject = 0
      DO k=1,2 
         DO j=1,ntheta
            DO i=1,nradius
               IF (vp(i,j,k)*(k-1.5).LT.0.0) THEN
                  numinject = numinject + 1
                  IF (numinject.GT.numinjectmax) THEN
                     WRITE (*,*) 'numinjectmax too small'
                     STOP
                  ENDIF

                  cum = cum + density(i,j,k)*radii(i)*drad(i)*
     &                 ABS(radii(i)*vp(i,j,k))
                  sintheta = sin(thetas(j))
                  costheta = cos(thetas(j))
                  xinjecttable(1,numinject) = cum
c
c--xinjecttable(2,*) is negative for k=1 (inner radius)
c
               xinjecttable(2,numinject)=2.0*(k-1.5)*radii(i)*sintheta
                  xinjecttable(3,numinject) = radii(i)*costheta
                  IF (k-1.5 .LT. 0.0) THEN
                     phi = -phibound
                  ELSE
                     phi = phibound
                  ENDIF
                  sinphi = sin(phi)
                  cosphi = cos(phi)
                  xinjecttable(4,numinject) = vr(i,j,k)*cosphi*sintheta
     &         - vp(i,j,k)*sinphi*sintheta + vt(i,j,k)*costheta*cosphi
                  xinjecttable(5,numinject) = vr(i,j,k)*sinphi*sintheta
     &         + vp(i,j,k)*cosphi*sintheta + vt(i,j,k)*costheta*sinphi
                  xinjecttable(6,numinject) = vr(i,j,k)*costheta - 
     &                 vt(i,j,k)*sintheta
c                  write (*,*) xinjecttable(6,numinject),costheta,
c     &                 radii(i),thetas(j),vr(i,j,k),vt(i,j,k)
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
