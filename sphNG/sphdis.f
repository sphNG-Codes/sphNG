      SUBROUTINE sphdis(igeom, idist, np, h1, facx, facy, facz, 
     &                   delx, dely, nx, ny, nz)
c************************************************************
c                                                           *
c  This subroutine positions particles in a spherical       *
c     distribution                                          *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/setlocal'

      REAL*4 velmax

      CHARACTER*20 filename
      CHARACTER*1 prof
c
c--Allow for tracing flow
c
      third = 1./3.
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine sphdis ')
c
c--Set Centrally Condensed Spherical Distribution
c
      IF ((idist.EQ.1) .OR. (idist.EQ.2)) THEN
         CALL unifdis(igeom, idist, np, h1, facx, facy, facz, 
     &                   delx, dely, nx, ny, nz, 0)
         rcyl2 = rcyl * rcyl 
         rmind2 = rmind * rmind
         rmax2 = rmax * rmax
         WRITE (*, 88001)
88001    FORMAT (' enter profile along r ',/,
     &           '                  r^-1 : 1',/,
     &           '                  r^-2 : 2',/,
     &           '           exponential : e',/,
     &           '             velocity  : v',/,  
     &           ' cos(mtheta) perturbn. : c')
         READ (*, 99004) prof
         iprofr = 0
         IF (prof.NE.'e') THEN
            IF (prof.EQ.'1') iprofr = 1
            IF (prof.EQ.'2') iprofr = 2
         ENDIF
         profr = iprofr/2.
c
c--Deform particle grid 
c
         IF (iprofr.EQ.1) THEN
            DO i = 1, npart
               xi = xyzmh(1,i)
               yi = xyzmh(2,i)
               zi = xyzmh(3,i)
               rxy2 = xi*xi + yi*yi
               r2 = rxy2 + zi*zi
               theta1 = ATAN2(yi,xi)
               theta2 = ATAN2(zi,SQRT(rxy2))
               r1 = rmax*(r2/rmax2)**0.75
               xyzmh(1,i) = r1*COS(theta1)*COS(theta2)
               xyzmh(2,i) = r1*SIN(theta1)*COS(theta2)
               xyzmh(3,i) = r1*SIN(theta2)
            END DO
         ELSEIF (prof.EQ.'c') THEN
c
c--linear cos(mtheta) perturbation by deforming grid
c
            WRITE (*,*) 'Warning: perturbation correct only for linear'
            WRITE (*,99034)
99034       FORMAT(' Enter m, and density contrast ')
            READ(*,*) m, densc
            IF (m .LE. 0) STOP 'ERROR: m must be > 0'
            WRITE (*,*) 'Setting perturbation m = ',m,' ampl = ',densc
            DO i = nptmass+1, npart
               xi = xyzmh(1,i)
               yi = xyzmh(2,i)
               zi = xyzmh(3,i)
               rxy = sqrt(xi*xi + yi*yi)
               theta1 = ATAN2(yi,xi)
               dtheta = -densc*SIN(m*theta1)/m
               theta1 = theta1 + dtheta
               xyzmh(1,i) = rxy*COS(theta1)
               xyzmh(2,i) = rxy*SIN(theta1)
            END DO
         ELSEIF (iprofr.EQ.2 .OR. prof.EQ.'e' .OR. prof.EQ.'v') THEN
            WRITE (*,99002)
99002       FORMAT ('NOT IMPLEMENTED')
            CALL quit
         ENDIF
      ELSE
         npart = np + nptmass
         rcyl2 = rcyl * rcyl 
         rmind2 = rmind * rmind
         rmax2 = rmax * rmax
         WRITE (*, 88001)
         READ (*, 99004) prof
         iprofr = 0
         IF (prof.NE.'e' .AND. prof.NE.'v') THEN
            IF (prof.EQ.'1') iprofr = 1
            IF (prof.EQ.'2') iprofr = 2
         ENDIF
         profr = iprofr/2.

         IF (prof.EQ.'v') THEN
            WRITE (*, 88401)
88401       FORMAT ('Enter name of file')
            READ (*,99010) filename
99010       FORMAT(A20)

            WRITE (*,*) 'Enter size of velocity files (e.g.N=32^3)'
            READ (*,*) nspace

            OPEN (45,FILE=filename,FORM='unformatted')
            READ (45) (((velx(i,j,k), i=1,nspace),j=1,nspace), 
     &           k=1,nspace)
            CLOSE (45)

            velmax = 0.
            DO k = 1, nspace
               DO j = 1, nspace
                  DO i = 1, nspace
                     velmax = MAX(velmax, velx(i,j,k))
                  END DO
              END DO
            END DO
            WRITE (*,*) 'velmax ',velmax
         ENDIF

         rm20 = rmax2 * 0.58 * 0.58
         xmax5 = xmax * 0.5
         ymax5 = ymax * 0.5
         zmax5 = zmax * 0.5
         probavn = 0.

         DO i = nptmass + 1, npart
            IF (MOD(i,1000).EQ.0) write (*,*) i
 100        x1 = 2.*(xmax*ran1(1) - xmax5)
            y1 = 2.*(ymax*ran1(1) - ymax5)
            z1 = 2.*(zmax*ran1(1) - zmax5)

            IF (prof.EQ.'v') THEN
               weight = velweight(nspace,x1,y1,z1,velx)/velmax
               IF (ran1(1).GT.weight) GOTO 100
            ENDIF

            rc2 = x1*x1 + y1*y1
            r2 = rc2 + z1*z1
            IF ((igeom.EQ.2) .AND. (rc2.GT.rcyl2)) GOTO 100
            IF ((igeom.EQ.3) .AND. (r2.GT.rmax2)) GOTO 100
            IF ((igeom.EQ.4) .AND. ((rc2.GT.rcyl2) .OR. 
     &           (rc2.LT.rmind2))) GOTO 100
            IF ((igeom.EQ.7) .AND. ((r2.GT.rmax2) 
     &           .OR. (r2.LT.rmind2))) GOTO 100

            IF (prof.NE.'e') THEN
               probr = 0.01 * (rmax2/r2) ** profr
            ELSE
               probr = 0.05 * 20.0 * EXP(-1.0*r2/rm20)
            ENDIF

            IF (ran1(1).GT.probr) GOTO 100

            probavn = probavn + probr
            xyzmh(1,i) = x1
            xyzmh(2,i) = y1
            xyzmh(3,i) = z1
            xyzmh(5,i) = 1.0 !!(1.0 / zprobr) ** third
         END DO

         probavn = probavn / FLOAT(npart - nptmass)
         WRITE (*,*) probavn
         DO i = nptmass + 1, npart
            xyzmh(5,i) = xyzmh(5,i) * h1 * probavn ** third
         END DO
      ENDIF

      DO i = nptmass + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END



      FUNCTION velweight(nspace, x, y, z, vel)

      INCLUDE 'idim'

      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/setlocal'

      REAL*4 vel(nvelmax,nvelmax,nvelmax)
      
      deli = 1.0/REAL(nspace/2)

      iposx = INT(x/rmax*(nspace/2)+(nspace/2))
      iposy = INT(y/rmax*(nspace/2)+(nspace/2))
      iposz = INT(z/rmax*(nspace/2)+(nspace/2))

      delx = x - (iposx-(nspace/2))/REAL(nspace/2)*rmax
      dely = y - (iposy-(nspace/2))/REAL(nspace/2)*rmax
      delz = z - (iposz-(nspace/2))/REAL(nspace/2)*rmax
c     
c--Find interpolated velocity
c
      velx1 = vel(iposx,iposy,iposz) + delx/deli*
     &     (vel(iposx+1,iposy,iposz)-vel(iposx,iposy,iposz))
      velx2 = vel(iposx,iposy+1,iposz) + delx/deli*
     &     (vel(iposx+1,iposy+1,iposz)-vel(iposx,iposy+1,iposz))
      vely1 = velx1 + dely/deli*(velx2-velx1)
      
      velx1 = vel(iposx,iposy,iposz+1) + delx/deli*
     &     (vel(iposx+1,iposy,iposz+1)-vel(iposx,iposy,iposz+1))
      velx2 = vel(iposx,iposy+1,iposz+1) + delx/deli*
     &     (vel(iposx+1,iposy+1,iposz+1)-vel(iposx,iposy+1,iposz+1))
      vely2 = velx1 + dely/deli*(velx2-velx1)
      
      velweight = vely1 + delz/deli*(vely2-vely1)

      RETURN
      END
