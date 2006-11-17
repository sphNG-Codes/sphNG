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
      INTEGER ipart,nx,nz,nphi,k,j,i

      nx = 10
      nz = 10
      deltaz = 2.*atorus/nz
      deltax = 2.*atorus/nx
      ipart = 0

      DO k=1,nz
         ztorus = (k-1)*deltaz - atorus
         DO j=1,nx        
            xtorus = (j-1)*deltax - atorus
            
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
      rhozero = 0.001
      totvol = (pi*atorus**2)*2.*pi*Rtorus
      totmass = rhozero*totvol
      massp = totmass/npart
      
      DO i=1,npart
         xyzmh(4,i) = massp
         xyzmh(5,i) = 1.2*deltax
      ENDDO
      
      WRITE(*,*) 'npart = ',npart,' particle mass = ',massp,
     &           ' denszero = ',rhozero
      
      RETURN

      END
