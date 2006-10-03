      SUBROUTINE unifdis(igeom, idist, np, h1, facx, facy, facz,
     &                   delx, dely, nx, ny, nz, ibound)
c************************************************************
c                                                           *
c  This subroutine positions particles in a uniform         *
c     coordinate distribution                               *
c                                                           *
c  Changes by D.Price: added hfact so that the smoothing    *
c  length is hfact * particle spacing (previous was unity)  *
c  This determines the number of neighbours                 *
c  Also added iexcludepart function to neaten things up     *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/flag'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine unifdis ')

      hfact = 1.0
      delz = 0.
c
c--Set uniform particle distribution
c
      notdone = 0
      rmax2 = rmax * rmax
      rcyl2 = rcyl * rcyl
      rmind2 = rmind * rmind
      zmax2 = zmax*zmax
      npartold = max(npart,nptmass)

      IF (idist.EQ.3) THEN
         xmax5 = 0.5*xmax
         ymax5 = 0.5*ymax
         zmax5 = 0.5*zmax
         npart = npartold + np
         DO i = npartold + 1, npart
 100        xi = 2.0*(xmax*ran1(1) - xmax5)
            yi = 2.0*(ymax*ran1(1) - ymax5)
            zi = 2.0*(zmax*ran1(1) - zmax5)
            d2 = xi**2 + yi**2 + zi**2
            r2 = d2 - zi*zi
            
            IF (iexcludepart(igeom,ibound,xi,yi,zi,r2,d2,xmax,ymax,zmax,
     &          rcyl2,rmax2,zmax2,rmind2).NE.0) GOTO 100

            xyzmh(1,i) = xi
            xyzmh(2,i) = yi
            xyzmh(3,i) = zi
            xyzmh(5,i) = hfact*h1
         END DO
c
c--Rings - so that m-modes are all zero
c
      ELSEIF (idist.EQ.6) THEN
         npart = npartold

         ntheta1 = 11
         steptheta = 0.99999*2.0*pi/ntheta1

         stepr = h1
         nr = INT(rmax/stepr)
         nz = INT(2.0*rmax/stepr)
         stepz = 2.0*rmax/nz

         WRITE (*,*) 'nr = ',nr,' nz = ',nz

         DO ir = 1, nr
            WRITE (*,*) 'ir = ',ir
c            radius = (ir - 0.5)*stepr
            radius = ir*stepr
            ntheta = INT(2.0*pi*radius/(steptheta*stepr))

            WRITE (*,*) 'ntheta = ',ntheta

            DO iz = 1, nz
               zi = (iz - 0.5)*stepz - rmax

               DO itheta = 1, ntheta
c                  theta = (itheta - 1./3.*MOD(iz,3))*steptheta
                  theta = (itheta - 0.5)*(steptheta/ir)
                  xi = radius*COS(theta)
                  yi = radius*SIN(theta)
                  d2 = xi**2 + yi**2 + zi**2
                  r2 = d2 - zi*zi
                  
                  IF (iexcludepart(igeom,ibound,xi,yi,zi,r2,d2,xmax,
     &               ymax,zmax,rcyl2,rmax2,zmax2,rmind2).EQ.0) THEN
                     npart = npart + 1
                     xyzmh(1,npart) = xi
                     xyzmh(2,npart) = yi
                     xyzmh(3,npart) = zi
                     xyzmh(5,npart) = hfact*stepr
                  ENDIF
               END DO
            END DO
         END DO
c
c--Custom
c
      ELSEIF (idist.EQ.5) THEN
         stepx = h1*facx
         stepy = h1*facy
         stepz = h1*facz
         k = 0
         l = 1
         m = 1
         npart = npartold
         DO i = npartold+ 1, npartold+np/2
            k = k + 1
            IF (k.GT.nx/2) THEN
               k = 1
               l = l + 1
               IF (l.GT.ny) THEN
                  l = 1
                  m = m + 1
                  IF (m.GT.nz) THEN
                     k = 1
                     l = 1
                     nz = nz + 1
                  ENDIF
               ENDIF
            ENDIF

            xstart = xmin + 0.5*stepx
            ystart = ymin + 0.5*stepy
            zstart = zmin + 0.5*stepz

            xi = xstart + (k - 1)*stepx
            yi = ystart + (l - 1)*stepy
            zi = zstart + (m - 1)*stepz 
            r2 = xi*xi + yi*yi
            d2 = r2 + zi*zi
            
            IF (iexcludepart(igeom,ibound,xi,yi,zi,r2,d2,xmax,ymax,
     &          zmax,rcyl2,rmax2,zmax2,rmind2).EQ.0) THEN
               npart = npart + 1
               IF (npart.GT.idim) THEN
                  WRITE (*, 99002)
99002             FORMAT (' dimensions too small ! ', /,
     &            ' choices  : 1  =  change dimensions and recompile'
     &            , /, '            2  =  change input parameters')
                  READ (*, *) ichoice
                  IF (ichoice.EQ.1) CALL quit
                  notdone = 1
                  GOTO 300
               ENDIF
               IF (xi.GT.xmax) xi = xmax
               xyzmh(1,npart) = xi
               IF (yi.GT.ymax) yi = ymax
               xyzmh(2,npart) = yi
               IF (zi.GT.zmax) zi = zmax
               xyzmh(3,npart) = zi
               xyzmh(5,npart) = hfact*h1
            ENDIF
         END DO

         k = 0
         l = 1
         m = 1
         DO i = 1, np/2
            k = k + 1
            IF (k.GT.nx/2) THEN
               k = 1
               l = l + 1
               IF (l.GT.ny) THEN
                  l = 1
                  m = m + 1
                  IF (m.GT.nz) THEN
                     k = 1
                     l = 1
                     nz = nz + 1
                  ENDIF
               ENDIF
            ENDIF
            
            xstart = (xmax + xmin)/2.0 + 0.5*stepx
            ystart = ymin
            zstart = zmin

            xi = xstart + (k - 1)*stepx
            yi = ystart + (l - 1)*stepy
            zi = zstart + (m - 1)*stepz 
            r2 = xi*xi + yi*yi
            d2 = r2 + zi*zi
            IF (iexcludepart(igeom,ibound,xi,yi,zi,r2,d2,xmax,ymax,
     &          zmax,rcyl2,rmax2,zmax2,rmind2).EQ.0) THEN
               npart = npart + 1
               IF (npart.GT.idim) THEN
                  WRITE (*, 99002)
                  READ (*, *) ichoice
                  IF (ichoice.EQ.1) CALL quit
                  notdone = 1
                  GOTO 300
               ENDIF
               IF (xi.GT.xmax) xi = xmax
               xyzmh(1,npart) = xi
               IF (yi.GT.ymax) yi = ymax
               xyzmh(2,npart) = yi
               IF (zi.GT.zmax) zi = zmax
               xyzmh(3,npart) = zi
               xyzmh(5,npart) = hfact*h1
            ENDIF
         END DO

      ELSE
         stepx = h1*facx
         stepy = h1*facy
         stepz = h1*facz
         k = 0
         l = 1
         m = 1
         npart = npartold
         DO 200 i = npartold + 1, npartold+np
            k = k + 1
            IF (k.GT.nx) THEN
               k = 1
               l = l + 1
               IF (l.GT.ny) THEN
                  l = 1
                  m = m + 1
                  IF (m.GT.nz) THEN
c                     k = 1
c                     l = 1
                     nz = nz + 1
                     print*,'nz =',nz
                  ENDIF
               ENDIF
            ENDIF
            
            IF (idist.EQ.2) THEN
c	          xstart = (xmin+xmax)/2.0
c     &               -(INT(ABS(0.5*(xmax-xmin)/delx)+0.5))*delx
c     &               + 0.5*delx	
c                  ystart = (ymin+ymax)/2.0
c     &              -(INT(ABS(0.5*(ymax-ymin)/dely)+0.5))*dely
c     &               + 0.5*dely
c                  zstart = (zmin+zmax)/2.0
c     &              -(INT(ABS(0.5*(zmax-zmin)/stepz)+0.5))*stepz
c     &               + 0.5*stepz
c               xstart = FLOAT(INT(xmin/stepx))*stepx + 0.5*delx
c               ystart = FLOAT(INT(ymin/stepy))*stepy + 0.5*dely
c     &              + stepy/2.0 
c     &              - dely/2.0
c               zstart = FLOAT(INT(zmin/stepz))*stepz + 0.5*stepz
               xstart = xmin + 0.125*delx
	       ystart = ymin + 0.25*dely
	       zstart = zmin + 0.5*stepz
               jy = MOD(l, 2)
               jz = MOD(m, 3)
	       IF (jz.EQ.0) THEN	! 3rd layer
		  ystart = ystart + 2.*dely
	          IF (jy.EQ.0) xstart = xstart + delx		   
	       ELSEIF (jz.EQ.2) THEN	! 2nd layer	  
	          ystart = ystart + dely
		  IF (jy.EQ.1) xstart = xstart + delx
	       ELSEIF (jy.EQ.0) THEN    ! first layer	  
	          xstart = xstart + delx
	       ENDIF

            ELSEIF (idist.EQ.1) THEN
	          xstart = (xmin+xmax)/2.0
     &               -(INT(ABS(0.5*(xmax-xmin)/stepx)+0.5))*stepx
     &               + 0.5*stepx	
                  ystart = (ymin+ymax)/2.0
     &              -(INT(ABS(0.5*(ymax-ymin)/stepy)+0.5))*stepy
     &               + 0.5*stepy
                  zstart = (zmin+zmax)/2.0
     &              -(INT(ABS(0.5*(zmax-zmin)/stepz)+0.5))*stepz
     &               + 0.5*stepz

!               xstart = (xmin+xmax)/2.0-(INT(ABS(xmin/stepx)+0.5))*stepx
!     &              + 0.5*stepx
!               ystart = (ymin+ymax)/2.0-(INT(ABS(ymin/stepy)+0.5))*stepy
!     &              + 0.5*stepy
!               zstart = (zmin+zmax)/2.0-(INT(ABS(zmin/stepz)+0.5))*stepz
!     &              + 0.5*stepz
            ELSEIF (idist.EQ.4) THEN
               jx = MOD(k, 2)
               IF (jx.EQ.1) THEN
                  xstart = xmin + 0.5*stepx
                  ystart = ymin + dely
                  zstart = zmin + delz
               ELSE
                  xstart = xmin + 0.5*stepx
                  ystart = ymin + 0.5*stepy + dely
                  zstart = zmin + 0.5*stepz + delz
               ENDIF
            ELSE
               WRITE(*,*) 'weird choice of idist in unifdis'
               WRITE(*,*) 'means xstart is used uninitialised'
               xstart = 0. ! to stop compiler warnings
               ystart = 0.
               zstart = 0.
               CALL quit
            ENDIF
            xi = xstart + FLOAT(k - 1)*stepx
            yi = ystart + FLOAT(l - 1)*stepy
            zi = zstart + FLOAT(m - 1)*stepz 
            r2 = xi*xi + yi*yi
            d2 = r2 + zi*zi
            IF (iexcludepart(igeom,ibound,xi,yi,zi,r2,d2,xmax,ymax,
     &          zmax,rcyl2,rmax2,zmax2,rmind2).EQ.0) THEN
               npart = npart + 1
               IF (npart.GT.idim) THEN
                  WRITE (*, 99002)
                  READ (*, *) ichoice
                  IF (ichoice.EQ.1) CALL quit
                  notdone = 1
                  GOTO 300
               ENDIF
               IF (xi.GT.xmax) xi = xmax
               xyzmh(1,npart) = xi
               IF (yi.GT.ymax) yi = ymax
               xyzmh(2,npart) = yi
               IF (zi.GT.zmax) zi = zmax
               xyzmh(3,npart) = zi
               xyzmh(5,npart) = hfact*h1
            ENDIF
 200     CONTINUE
      ENDIF

 300  DO i = npartold + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END
c
c--this function is to neaten up unifdis
c  determines whether or not particle is within prescribed boundaries
c  for the given geometry
c  DJP 02.10.06
c     
      FUNCTION iexcludepart(igeom,ibound,xi,yi,zi,r2,d2,xmax,ymax,zmax,
     &                      rcyl2,rmax2,zmax2,rmind2)

      INCLUDE 'COMMONS/setbin'
      INTEGER iexcludepart

      DATA small/0.0001/
        
      iexcludepart = 0      
      
      IF ((igeom.EQ.1) .AND. ((xi.GT.xmax + small) .OR.
     &    (yi.GT.ymax + small) .OR. (zi.GT.zmax + small))) THEN
         iexcludepart = 1
      ELSEIF ((igeom.EQ.2) .AND. (r2.GT.rcyl2)) THEN
         iexcludepart = 1
      ELSEIF ((igeom.EQ.3) .AND. (d2.GT.rmax2)) THEN
         iexcludepart = 1
      ELSEIF ((igeom.EQ.4) .AND. ((r2.LT.rmind2).OR.(r2.GT.rcyl2))) THEN
         iexcludepart = 1
      ELSEIF ((igeom.EQ.5) .AND. ((r2/rcyl2 + zi*zi/zmax2.GT.1.0)
     &   .OR. (r2/rmind2 + zi*zi/(zmax2*(rmind2/rcyl2)).LT. 1.0))) THEN
         iexcludepart = 1
      ELSEIF (igeom.EQ.6) THEN
         IF (d2.LE.rmax2 .AND. d2.GE.rmind2) THEN
            r1 = SQRT(r2)
            IF (ibound.EQ.90) THEN
               angvelhere = angvel
            ELSEIF (ibound.EQ.91) THEN
               angvelhere = angvel*(r2/d2)
            ELSE
               angvelhere = 0.
            ENDIF
            accelcent = r1*(angvelhere/r2)**2
            accelpt = totptmass/d2*SIN(ATAN2(r1,ABS(zi)))
            IF (accelcent.GT.accelpt) iexcludepart = 1
         ELSE
            iexcludepart = 1
         ENDIF
      ELSEIF ((igeom.EQ.7) .AND. ((d2.GT.rmax2).OR.(d2.LT.rmind2))) THEN
         iexcludepart = 1
      ELSEIF ((igeom.EQ.8) .AND. ((d2.LE.rmax2) .OR. 
     &    ((xi.GT.xmax + small) .OR. (yi.GT.ymax + small) 
     &     .OR. (zi.GT.zmax + small)))) THEN
         iexcludepart = 1
      ENDIF
      
      RETURN
      END
