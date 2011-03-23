      SUBROUTINE setup_torus(xyzmh,vxyzu,rhozero,RK2, npart,ntot)
c************************************************************
c                                                           *
c     setup for Tokamak torus                               *
c                                                           *
c************************************************************

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/xforce'     
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/tokamak'

      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL pi
      PARAMETER (pi=3.141592653589)
      INTEGER npart,ntot
      REAL rhozero, RK2, gama1
      REAL deltaz,deltax,deltaphi,rcyl,rintorus2,phi
      REAL xtorus,ztorus,massp,totmass,totvol
      REAL phitorus, randphi, ran1, thetator, r1, rho_c
      REAL rintorus, xi, yi, zi, rfactor, kappa, ra2, pri
      INTEGER ipart,nx,nz,nphi,k,j,i

      nx = 15
      nz = 15
      deltaz = 2.*atorus/nz
      deltax = 2.*atorus/nx
      ipart = 0
      encal  = 'i'
      gamma = 5./3.
      
c Particle and eos properties/constants
      rhozero = 1.
      hzero = 0.6*deltax
      ra2 = da2*(atorus-hzero)**2      
      pri = currj0*currj0*atorus*atorus*(47.-12.*ra2**5+
     &      75.*ra2**4-200.*ra2**3+270.*ra2**2-180.*ra2)/720.      
      RK2 = pri/rhozero*gamma
      gama1 = gamma-1.
      
      DO k=1,nz
         ztorus = (k-0.5)*deltaz - atorus
         DO j=1,nx        
            xtorus = (j-0.5)*deltax - atorus
            
            rintorus2 = xtorus**2 + ztorus**2
            
            IF (rintorus2.LT.atorus**2) THEN
               rcyl = xtorus + Rtorus

               deltaphi = deltaz*Rtorus/rcyl
               nphi = int(2.*pi/deltaphi)
               deltaphi = 2.*pi/nphi
               randphi  = ran1(1)*deltaphi
            !--make ring of particles at r=rcyl
               DO i=1,nphi                   
                  ipart = ipart + 1
                  IF (ipart.GT.idim) STOP 'setup_torus: dims too small'
                  phi = (i-1)*deltaphi
                  xyzmh(1,ipart) = rcyl*COS(phi+randphi)
                  xyzmh(2,ipart) = rcyl*SIN(phi+randphi)
                  xyzmh(3,ipart) = ztorus
                  ra2 = rintorus2*da2
                  pri = currj0*currj0*atorus*atorus*(47.-12.*ra2**5+
     &                 75.*ra2**4-200.*ra2**3+270.*ra2**2-180.*ra2)/720.
c                 zero velocities
                  vxyzu(1,i) = 0.
                  vxyzu(2,i) = 0.
                  vxyzu(3,i) = 0.       
                  vxyzu(4,ipart) = pri/gama1/rhozero                 
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      
      npart = ipart
      ntot = ipart
c!     choose an initial density
      totvol = (pi*atorus**2)*2.*pi*Rtorus      
      totmass = totvol*rhozero !(pi*currj0*atorus**2)**2*Rtorus/36.
      massp = totmass/npart
c      rhozero = totmass/totvol
      
      DO i=1,npart
         xyzmh(4,i) = massp
         xyzmh(5,i) = hzero*2.
c        zero velocities
c         vxyzu(1,i) = 0.
c         vxyzu(2,i) = 0.
c         vxyzu(3,i) = 0.         
c         vxyzu(4,i) = RK2*(rhozero**gama1)
      ENDDO

c
c--Stretching the spatial distribution to have exponential density distribution
c      
c      kappa = 4.
c      rho_c = rhozero*kappa/(1.- EXP(-kappa))
c      DO i = 1, npart
c         xi = xyzmh(1,i)
c         yi = xyzmh(2,i)
c         zi = xyzmh(3,i)
c         rcyl = SQRT(xi*xi+yi*yi)
c         rintorus2 = (rcyl - Rtorus)**2 + zi*zi
c         rintorus = sqrt(rintorus2)
c         IF (rintorus.gt.tiny) THEN
c            phitorus = ATAN2(yi,xi)
c            thetator = ATAN2(zi,rcyl-Rtorus)
c            rfactor = 1.- kappa*rhozero/rho_c*(rintorus/atorus)**2
c            r1 = atorus*SQRT(-LOG(rfactor)/kappa)
c            xyzmh(1,i) = (Rtorus+r1*COS(thetator))*COS(phitorus)
c            xyzmh(2,i) = (Rtorus+r1*COS(thetator))*SIN(phitorus)
c            xyzmh(3,i) = r1*SIN(thetator)
c           ENDIF   
c        END DO        
        
      WRITE(*,*) 'npart = ',npart,' particle mass = ',massp,
     &           ' denszero = ',rhozero
      
      RETURN

      END
