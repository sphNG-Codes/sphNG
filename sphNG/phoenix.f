      SUBROUTINE phoenix(i, idtsyn, itime)
c***********************************************************************
c                                                                      *
c  Assigns an accreted particle or particle that goes outside the      *
c     dead boundary a new position, velocity etc to allow accretion    *
c     flow to stablise.                                                *
c                                                                      *
c***********************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/timei'
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
      INCLUDE 'COMMONS/timeextra'

      radinject = 0.98*deadbound
      radinject2 = radinject**2

 10   rnd1 = ran1(1)
      rnd2 = ran1(1)
      rnd3 = ran1(1)

      ang1 = 2*pi*rnd1
      ang2 = pi*rnd2
      s2 = SIN(ang2)
c
c--Make UNIFORM distribution over the surface of a SPHERE.
c
      IF (rnd3.GT.s2) GOTO 10

      xyzmh(1,i) = radinject*COS(ang1)*s2
      xyzmh(2,i) = radinject*SIN(ang1)*s2
      xyzmh(3,i) = radinject*COS(ang2)

      rxy2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
      rxy = SQRT(rxy2)
      r2 = rxy2 + xyzmh(3,i)**2
c
c--If ibound=90, then inject all particles with same specific angular momentum
c     ibound=91, then inject particles with uniform angular velocity at
c                    injection radius, radinject
c
      IF (ibound.EQ.90) THEN
         fractanhere = fractan
      ELSEIF (ibound.EQ.91) THEN
         fractanhere = fractan*(rxy2/radinject2)
      ELSE
         fractanhere = fractan
         WRITE (iprint,*) 'ERROR in phoenix'
         CALL quit(1)
      ENDIF

      accelcent = rxy*(fractanhere*specang/rxy2)**2
      accelpt = ptmassin/r2*SIN(ATAN2(rxy,ABS(xyzmh(3,i))))
      IF (accelcent.GT.accelpt) GOTO 10

      r = SQRT(r2)
c
c--Velocity equals tangential + radial contributions.  The tangential
c     velocity is defined as a fraction (fractanhere) of the specific angular 
c     momentum (specang) of a particle in circular orbit about the total 
c     mass of the binary at the semi-major axis of the system (=1 in 
c     code units).
c     The radial velocity is defined as a fraction (fracradial) of the
c     orbital velocity at the semi-major axis.
c
      tanmag = fractanhere*specang/rxy
c
c--Maximum radial velocity gives kinetic energy of particle to be that of a 
c     particle which has fallen from infinity, accounting for tangential 
c     velocity.  Note that specang=SQRT(Mtot).
c
      radmag = fracradial*SQRT(2.0*specang*specang/r - tanmag*tanmag)
      vxyzu(1,i) = - radmag*xyzmh(1,i)/r - tanmag*xyzmh(2,i)/rxy
      vxyzu(2,i) = - radmag*xyzmh(2,i)/r + tanmag*xyzmh(1,i)/rxy
      vxyzu(3,i) = - radmag*xyzmh(3,i)/r
      vxyzu(4,i) = vxyzu(4,i)
      xyzmh(5,i) = deadbound/(FLOAT(npart)**(1.0/3.0))
      xyzmh(4,i) = xyzmh(4,3)

      DO j = 1, 5
         dumxyzmh(j,i) = xyzmh(j,i)
      END DO
      DO j = 1, 4
         dumvxyzu(j,i) = vxyzu(j,i)
      END DO

      f1vxyzu(1,i) = 0.0
      f1vxyzu(2,i) = 0.0
      f1vxyzu(3,i) = 0.0
      f1vxyzu(4,i) = 0.0
      f1ha(1,i) = 0.0

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
