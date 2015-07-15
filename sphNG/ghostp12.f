      SUBROUTINE ghostp12(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c************************************************************
c                                                           *
c  This subroutine computes the list of ghost particles for *
c     treating the boundaries.                              *
c                                                           *
c************************************************************

      
      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      DIMENSION vxyzu(4,idim)
      DIMENSION ekcle(4,iradtrans)
      DIMENSION Bevolxyz(imhdevol,imhd)

      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/cylinder'
      INCLUDE 'COMMONS/eosq'

      CHARACTER*7 where

      DATA where/'ghostp12'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp12')

      nghost = 0
      uradconst = radconst/uergcc
c
c--Find z-limits from cylinder common block
c
      zmax = length/2.
      zmin = - zmax
      
c 
c--define r-limits from cylinder common block
c
      rcyl = radius
      rcyl2 = radius*radius

c     Force-free magnetic field configuration, in a cylinder
c     Muff = 2.5 set in Cylinder common

      nr = 20
      nz = int(nr*length/radius)    

      deltar = radius/nr
      deltaz = length/nz     
C      
C      do iz = 1, nz+6
C         do ir = nr+2, nr +6
C            rgr = deltar*(ir-1)
C            nc = 2*pi*rgr/deltar
C            angle = 2*pi/nc
C            do ic = 1, nc     
C               nghost = nghost+1
C               nptot = MIN0(npart + nghost, idim)
C               xyzmh(1,nptot) = rgr*COS(angle*ic)
C               xyzmh(2,nptot) = rgr*SIN(angle*ic)
C               xyzmh(3,nptot) = -length/2.+deltaz/2.+deltaz*(iz-4)
C               xyzmh(4,nptot) = xyzmh(4,1)
C               xyzmh(4,nptot) = 1.2*deltar
C               vxyzu(1,nptot) = 0.
C               vxyzu(2,nptot) = 0.
C               vxyzu(3,nptot) = 0. 
C               vxyzu(3,nptot) = 1.5 
C               rho(nptot) = 0.001
C               iphase(nptot) = 0                    
C               Bevolxyz(1,nptot) = 0.!-0.001*dbesj1(Muff*rr)
C     &                                     *xyzmh(2,ipart)/rgr
C               Bevolxyz(2,nptot) = 0.!0.001*dbesj1(Muff*rr)
C     &                                    *xyzmh(1,ipart)/rgr
C               Bevolxyz(3,nptot) = 0.!0.001*dbesj0(Muff*rr) 
C            enddo           
C         enddo
C      enddo         

     
c
c--Find ghost particles (for all particles within radkernel*h of boundary)
c
c
c  For periodic boundary conditions, safest is to use hmax
c
      hmax = 0.
      DO i = 1, npart
         hmax = max(1.1*xyzmh(5,i),hmax)
      ENDDO

      
      DO 200 i = 1, npart
         nghostold = nghost
         hasghost(i) = .FALSE.
         IF (iphase(i).NE.0) GOTO 200
         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         pmassi = xyzmh(4,i)
         hi = xyzmh(5,i)

         vxi = vxyzu(1,i)
         vyi = vxyzu(2,i)
         vzi = vxyzu(3,i)
         ui = vxyzu(4,i)
         rhoi = rho(i)
         vsoundi = vsound(i)
         presi = pr(i)

c
c--Z axis
c
         delta = 0.
         dzmin = (zi - zmin)/hmax
         IF (dzmin.GE.delta .AND. dzmin.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi + (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
            IF (imhd.EQ.idim) THEN
               DO k = 1, imhdevol
                  Bevolxyz(k,nptot) = Bevolxyz(k,i)
               END DO
            ENDIF
         ENDIF

         dzmax = (zmax - zi)/hmax
         IF (dzmax.GT.delta .AND. dzmax.LT.radkernel) THEN
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            xyzmh(1,nptot) = xi
            xyzmh(2,nptot) = yi
            xyzmh(3,nptot) = zi - (zmax-zmin)
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
            vxyzu(1,nptot) = vxi
            vxyzu(2,nptot) = vyi
            vxyzu(3,nptot) = vzi
            vxyzu(4,nptot) = ui
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
            IF (imhd.EQ.idim) THEN
               DO k = 1, imhdevol
                  Bevolxyz(k,nptot) = Bevolxyz(k,i)
               END DO
            ENDIF
         ENDIF

 200  CONTINUE
 
      ntot = npart + nghost     

c
c--If ibound = 12 then we have cylindrical boundaries
c--now find ghost for cylindrical boundaries
c
      DO 300 i = 1, ntot
         IF (iphase(i).NE.0) GOTO 300
         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         pmassi = xyzmh(4,i)
         hi = xyzmh(5,i)

         vxi = vxyzu(1,i)
         vyi = vxyzu(2,i)
         vzi = vxyzu(3,i)
         ui = vxyzu(4,i)
         rhoi = rho(i)      
         vsoundi = vsound(i)
         presi = pr(i)

         delta = 0.001*hi
         hi2 = hi*hi
         r2 = xi**2 + yi**2
         delta2 = delta*delta
         drmin = rcyl - radkernel*hi
         drmin2 = drmin*drmin
         drmax2 = rcyl2 + delta2 - 2*rcyl*delta
         IF (r2.GT.drmin2 .AND. r2.LT.drmax2.OR. drmin.LT.0.) THEN  
            hasghost(i) = .TRUE.
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            r = SQRT(r2)
            rgr = 2*rcyl - r
            xyzmh(1,nptot) = xi*rgr / r
            xyzmh(2,nptot) = yi*rgr / r
            xyzmh(3,nptot) = zi
            xyzmh(4,nptot) = pmassi
            xyzmh(5,nptot) = hi
c            v2 = vxi*vxi + vyi*vyi
c            vg = SQRT(v2)
c            ang1 = ATAN2(yi,xi)
c            ang2 = ATAN2(vyi,vxi)
c            angr = pi - ang2 + ang1
            rho(nptot) = rhoi
            vsound(nptot) = vsoundi
            pr(nptot) = presi
            iphase(nptot) = 0
            IF (imhd.EQ.idim) THEN
               Bxi = -ampl*dbesj1(Muff*rgr)*xyzmh(2,nptot)/rgr
               Byi = ampl*dbesj1(Muff*rgr)*xyzmh(1,nptot)/rgr
               Bzi = ampl*dbesj0(Muff*rgr)
               IF (varmhd(1:4).EQ.'Bvol') THEN              
                  Bevolxyz(1,nptot) = Bxi
                  Bevolxyz(2,nptot) = Byi
                  Bevolxyz(3,nptot) = Bzi
               ELSEIF (varmhd(1:4).EQ.'Brho') THEN
                  Bevolxyz(1,nptot) = Bxi/rhoi
                  Bevolxyz(2,nptot) = Byi/rhoi
                  Bevolxyz(3,nptot) = Bzi/rhoi
               ELSE
                  STOP 'ERROR: ghostp12+euler not implemented'
               ENDIF
c This gives a velocity field which follows magnetic field lines               
               vxyzu(1,nptot) = vampl*Bxi
               vxyzu(2,nptot) = vampl*Byi
               vxyzu(3,nptot) = vampl*Bzi  
            ELSE
               vxyzu(1,nptot) = 0.!vg*COS(angr)
               vxyzu(2,nptot) = 0.!vg*SIN(angr)
               vxyzu(3,nptot) = vzi            
            ENDIF            
            vxyzu(4,nptot) = ui            
         ENDIF
 300  CONTINUE         

      ntot = npart + nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      RETURN
      END
