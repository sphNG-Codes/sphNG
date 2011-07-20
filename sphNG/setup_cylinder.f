      SUBROUTINE setup_cylinder(xyzmh,vxyzu,Bevolxyz,rhozero,
     &           RK2,npart,ntot)
c************************************************************
c                                                           *
c     Cylindrical setup for Force-Free in cylinder          *     
c                                                           *
c************************************************************


      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/xforce' 
      INCLUDE 'COMMONS/cylinder'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/typef'
      
      REAL xyzmh(5,idim),vxyzu(4,idim), Bevolxyz(imhdevol,imhd)      
      REAL pi, gama1, rhozero, RK2, Slund, valfven
      PARAMETER (pi=3.141592653589)
      
      INTEGER npart,ntot
      INTEGER ipart,nr,nz,nc,ir,iz,ic,i, ipfree
      REAL massp,totvol,totmass
      REAL deltar,deltaz,rr,angle
      REAL offset, dbesj0, dbesj1, ran1, zz

c     Set Boundary To be Periodic in z
      ibound = 12
233   WRITE (*,*) 'Do you want pressure free cylinder? (0/1)'
      READ (*,*) ipfree
      
      IF (ipfree.EQ.0) THEN
         encal = 'p'
         gamma = 5./3.
         iexf = 10
         print *, 'encal =', encal, 'gamma =', gamma, 'iexf =',iexf 
      ELSEIF (ipfree.EQ.1) THEN
         encal = 'i'
         gamma = 1.
         iexf = 0
         print *, 'encal =', encal, 'gamma =', gamma, 'iexf =',iexf
      ELSE 
         GO TO 233
      ENDIF
      
      rcyl = radius     
      ipart = 0
      nr = 20
      nz = int(nr*length/radius)    

      deltar = radius/(nr+0.5)
      deltaz = length/nz     
      
c     Force-free magnetic field configuration, in a cylinder
c     Muff = 2.5 set in Cylinder common
      
      do iz = 1, nz
         zz = -length/2.+deltaz/2.+deltaz*(iz-1) 
         print *, 'zz =', zz, iz, deltaz
         do ir = 1, nr +1
            rr = deltar*(ir-1) 
            nc = 2*pi*rr/deltar
c            write (*,*) 'nc =', nc
c            if (ir.eq.nr+1) write (*,*) rr, radius-rr, 0.5*deltar
            if (nc.eq.0) then
              ipart = ipart +1
              IF (ipart.GT.idim) STOP 'setup_cyl: dims too small'
              xyzmh(1,ipart) = 0.d0
              xyzmh(2,ipart) = 0.d0
              xyzmh(3,ipart) = zz
              if (imhd.eq.idim) then
                 Bevolxyz(1,ipart) =  0.d0
                 Bevolxyz(2,ipart) =  0.d0
                 Bevolxyz(3,ipart) =  ampl*dbesj0(Muff*rr)
              endif
            else
              offset = ran1(1)
              do ic = 1, nc     
                 angle = 2*pi/nc
                 ipart = ipart+1
                 IF (ipart.GT.idim) STOP 'setup_cyl: dims too small'
                 xyzmh(1,ipart) = rr*COS(angle*(ic-0.5)+offset*angle)
                 xyzmh(2,ipart) = rr*SIN(angle*(ic-0.5)+offset*angle)
                 xyzmh(3,ipart) = -length/2.+deltaz/2.+deltaz*(iz-1)    
                 if (imhd.eq.idim) then 
                    Bevolxyz(1,ipart) = -ampl*xyzmh(2,ipart)/rr
     &                                     *dbesj1(Muff*rr)
                    Bevolxyz(2,ipart) = ampl*xyzmh(1,ipart)/rr
     &                                     *dbesj1(Muff*rr)                
                    Bevolxyz(3,ipart) = ampl*dbesj0(Muff*rr) 
                 endif   
              enddo
            endif                
         enddo
      enddo              
      
      npart = ipart
      ntot = ipart
!     choose an initial density
      rhozero = 1.e-3
      hzero = 1.2*deltar
      totvol = (pi*radius**2)*length
      totmass = rhozero*totvol
      massp = totmass/npart   
      RK2 = 1.5
      gama1 = gamma-1.
      DO i=1,npart
         vxyzu(1,i) = 0.
         vxyzu(2,i) = 0.
         vxyzu(3,i) = 0.
         IF (imhd.EQ.idim) THEN
            vxyzu(1,i) = vxyzu(1,i) + vampl*Bevolxyz(1,i)
            vxyzu(2,i) = vxyzu(2,i) + vampl*Bevolxyz(2,i)
            vxyzu(3,i) = vxyzu(3,i) + vampl*Bevolxyz(3,i)   
         ENDIF   
         IF (ipfree.EQ.0) THEN
            vxyzu(4,i) = RK2*(rhozero**gama1)
         ELSE
            vxyzu(4,i) = 0.
         ENDIF
         xyzmh(4,i) = massp
         xyzmh(5,i) = hzero
c-- Setup an initial uniform density         
         rho(i) = rhozero
      ENDDO

c-- Estimate Lundquist number
      valfven = ampl**2/rhozero
      Slund =radius*valfven/etamhd 
      
      WRITE(*,*) 'npart = ',npart,' particle mass = ',massp,
     &           ' denszero = ',rhozero
      WRITE(*,*) 'The Lundquist number = t_eta/t_alfven =', Slund
      
      RETURN      
      
      END



C      nr = 20
C      deltar = radius/nr
C      nz = int(2.0*radius/deltar)    
C      deltaz = 2.0*radius/nz     
C      nc1 = 5
C      steptheta = 0.99999*2.0*pi/nc1
C      
C      do ir=1, nr
C         write (*,*) 'ir=',ir
C         rr = ir*deltar
C         nc = int(2.0*pi*rr/(steptheta*deltar))
C         
C         write(*,*) 'ntheta',nc
C         
C         do iz =1, nz
C            zi = (iz-0.5)*deltaz - radius
C            
C            do ic = 1, nc
C               ipart = ipart+1
C               IF (ipart.GT.idim) STOP 'setup_cyl: dims too small'            
C               angle = (ic-0.5)*steptheta/ir
C               xyzmh(1,ipart) = rr*COS(angle)
C               xyzmh(2,ipart) = rr*SIN(angle)
C               xyzmh(3,ipart) = zi
C               Bevolxyz(1,ipart) =  -0.001*dbesj1(Muff*rr)
C     &                                     *xyzmh(2,ipart)/rr
C               Bevolxyz(2,ipart) =  0.001*dbesj1(Muff*rr)
C     &                                     *xyzmh(1,ipart)/rr
C               Bevolxyz(3,ipart) =  0.001*dbesj0(Muff*rr) 
C            enddo
C         enddo
C      enddo

      
      
C      do iz = 1,nz
C         do ix = 0, 2*nr-1
C            xcyl = - radius + deltar/2.+ix*deltar
C            do iy = 0,2*nr-1
C               ycyl = -radius + deltar/2.+iy*deltar
C               rincyl2 = xcyl*xcyl+ycyl*ycyl
C               IF (rincyl2.LE.radius**2) THEN
C                 ipart = ipart +1      
C                 IF (ipart.GT.idim) STOP 'setup_cyl: dims too small'
C                 xyzmh(1,ipart) = xcyl
C                 xyzmh(2,ipart) = ycyl
C                 xyzmh(3,ipart) = -length/2.+deltaz/2.+deltaz*(iz-1)     
C                 Bevolxyz(1,ipart) = 0. 
C                 Bevolxyz(2,ipart) = 0. 
C                 Bevolxyz(3,ipart) = 0.            
C               ENDIF
C            enddo
C         enddo   
C      enddo
          
