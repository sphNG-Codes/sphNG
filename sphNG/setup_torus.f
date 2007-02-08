      SUBROUTINE setup_torus(xyzmh,vxyzu,rhozero,npart,ntot)
c************************************************************
c                                                           *
c     setup for Tokamak torus                               *
c                                                           *
c************************************************************

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/tokamak'

      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL pi
      PARAMETER (pi=3.141592653589)
      INTEGER npart,ntot
      REAL rhozero
      REAL deltaz,deltax,deltaphi,rcyl,rintorus2,phi
      REAL xtorus,ztorus,massp,totmass,totvol
      REAL phitorus, thetator, r1, rho_c
      REAL rintorus, xi, yi, zi, rfactor, kappa
      INTEGER ipart,nx,nz,nphi,k,j,i

      nx = 10
      nz = 10
      deltaz = 2.*atorus/nz
      deltax = 2.*atorus/nx
      ipart = 0

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
            !--make ring of particles at r=rcyl
               DO i=1,nphi                   
                  ipart = ipart + 1
                  IF (ipart.GT.idim) STOP 'setup_torus: dims too small'
                  phi = (i-1)*deltaphi
                  xyzmh(1,ipart) = rcyl*COS(phi)
                  xyzmh(2,ipart) = rcyl*SIN(phi)
                  xyzmh(3,ipart) = ztorus
c                 zero velocities
                  vxyzu(1,ipart) = 0.
                  vxyzu(2,ipart) = 0.
                  vxyzu(3,ipart) = 0.
                  vxyzu(4,ipart) = 1.5 ! so that p = rho in the equation of state routine
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      
      npart = ipart
      ntot = ipart
!     choose an initial density
      totvol = (pi*atorus**2)*2.*pi*Rtorus
      totmass = (pi*currj0*atorus**2)**2*Rtorus/36.
      massp = totmass/npart
      rhozero = totmass/totvol
      
      DO i=1,npart
         xyzmh(4,i) = massp
         xyzmh(5,i) = 1.2*deltax
      ENDDO

c
c--Stretching the spatial distribution to have exponential density distribution
c
      kappa = 4.
      rho_c = rhozero*kappa/(1.- EXP(-kappa))
      DO i = 1, npart
         xi = xyzmh(1,i)
         yi = xyzmh(2,i)
         zi = xyzmh(3,i)
         rcyl = SQRT(xi*xi+yi*yi)
         rintorus2 = (rcyl - Rtorus)**2 + zi*zi
         rintorus = sqrt(rintorus2)
         IF (rintorus.gt.tiny) THEN
            phitorus = ATAN2(yi,xi)
            thetator = ATAN2(zi,rcyl-Rtorus)
            rfactor = 1.- kappa*rhozero/rho_c*(rintorus/atorus)**2
            r1 = atorus*SQRT(-LOG(rfactor)/kappa)
            xyzmh(1,i) = (Rtorus+r1*COS(thetator))*COS(phitorus)
            xyzmh(2,i) = (Rtorus+r1*COS(thetator))*SIN(phitorus)
	        xyzmh(3,i) = r1*SIN(thetator)
	     ENDIF   
	  END DO        
	  
      WRITE(*,*) 'npart = ',npart,' particle mass = ',massp,
     &           ' denszero = ',rhozero
      
      RETURN

      END
