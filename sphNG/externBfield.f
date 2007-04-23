      SUBROUTINE externBfield(xi,yi,zi,hi,vxi,vyi,vzi,rhoi,
     &                        Bintx, Binty, Bintz,
     &                        currJintx, currJinty, currJintz,
     &                        Bextx, Bexty, Bextz, 
     &                        fextx, fexty, fextz, 
     &                        vdotgradBx, vdotgradBy, vdotgradBz,
     &                        string)
c************************************************************
c                                                           *
c    This subroutine and associated routines handles        *
c    everything to do with external B fields                *
c    (with non-zero spatial derivatives)                    *
c                                                           *
c************************************************************      

      IMPLICIT NONE
      REAL tiny
      PARAMETER (tiny=1.e-14)
      INCLUDE 'COMMONS/tokamak'
c IN      
      REAL xi,yi,zi,hi,vxi,vyi,vzi,rhoi
      REAL Bintx,Binty,Bintz
      REAL currJintx,currJinty,currJintz
      CHARACTER*(*) string
c OUT
      REAL Bextx,Bexty,Bextz,fextx,fexty,fextz
      REAL vdotgradBx,vdotgradBy,vdotgradBz
c LOCAL
      REAL rcyl,drcyl,rintorus2,rintorus,ra2,term
      REAL costheta,sintheta,cosphi,sinphi,drintorus
      REAL Bextr,Bexttheta,Bextphi,dBthetadr
      REAL currJextr,currJexttheta,currJextphi
      REAL currJextx,currJexty,currJextz
      REAL frtorus,fbound,drhoi,valfven2
      REAL currjintbextx,currjintbexty,currjintbextz
      REAL currjextbintx,currjextbinty,currjextbintz
c
c--get coordinate factors in torus coordinate system
c
      CALL get_factors(xi,yi,zi,costheta,sintheta,cosphi,sinphi,
     &                 rcyl,drcyl,rintorus,rintorus2,drintorus)
c
c--get 1/rho
c
      IF (rhoi.GT.tiny) THEN 
         drhoi = 1./rhoi
      ELSE
         drhoi = 0.
         WRITE(*,*) ' rho <= tiny in externalBfield !!!'
         CALL quit
      ENDIF

c============ choose external B field profile ===============
c
c     iBext is a parameter which allows one to have different
c     choices for the external field profile 
c     (set in COMMONS/tokamak)
c      
      IF (iBext.EQ.1) THEN
c
c     these are the profiles as described in Wesson, Tokamaks
c
         ra2 = rintorus2*da2
         term = 1. - ra2
c
c        external currents
c
         currJextr = 0.
         currJexttheta = 0.
         currJextphi = currj0*term**nutorus
c
c        Bfield in torus "theta" direction
c
         Bextr = 0.
         Bexttheta = 0.5*currj0*atorus**2/(nutorus+1)*
     &            (1. - term**(nutorus+1))*drintorus
         Bextphi = Bphi
c
c        derivative of Btheta with respect to torus 'r'
c
         dBthetadr = 0.5*currj0*atorus**2/(nutorus+1)*
     &     (((2*nutorus + 1)*ra2 + 1.)*term**nutorus - 1.)*drintorus**2
c
c        external force is in torus "r" direction  (J X B)
c
         frtorus = -Bexttheta*currJextphi*drhoi
c
c        estimate of alfven speed used to give right order
c        of magnitude to boundary force
c        
         valfven2 = (Bexttheta*Bexttheta + Bextphi*Bextphi)*drhoi

      ELSEIF (iBext.EQ.2) THEN
c
c     this is a Gaussian profile for J designed so that
c     J does not go to zero at the edge but rather tails off gently
c
         ra2 = rintorus2*da2
         term = exp(-4.*ra2)
         currJextr = 0.
         currJexttheta = 0.
         currJextphi = currj0*term
         
         Bextr = 0.
         Bexttheta = 0.125*currj0*atorus**2*(1.-term)*drintorus
         Bextphi = 0.
         
         dBthetadr = -0.125*currj0*atorus**2*(1.-term)*drintorus**2
     &               +currj0*term
         
         frtorus = -Bexttheta*currJextphi*drhoi
         
         valfven2 = (Bexttheta*Bexttheta + Bextphi*Bextphi)*drhoi
      ELSE
         STOP 'unknown/unimplemented iBext in COMMONS/tokamak'
      ENDIF
c============================================================

      IF (string(1:6).eq.'Bfield') THEN
c
c--in this case return only the external B field
c  is returned in cartesian co-ordinates
c      
         CALL vec_rthetaphi_to_xyz(Bextr,Bexttheta,Bextphi,
     &                             Bextx,Bexty,Bextz,
     &                             costheta,sintheta,cosphi,sinphi)
      
      ELSEIF (string(1:4).eq.'fext') THEN
c
c--add boundary force to force in torus 'r' direction
c
         frtorus = frtorus + fbound(rintorus,hi,valfven2)
c
c--get J_ext x B_ext as the external force
c        
         CALL vec_rthetaphi_to_xyz(frtorus,0.,0.,
     &                             fextx,fexty,fextz,
     &                             costheta,sintheta,cosphi,sinphi)
      
      ELSEIF (string(1:3).eq.'all') THEN
c
c--in this case return the external force
c  including mixed Jint x Bext terms and also
c  the term -v.grad Bext needed in the B evolution equation
c      
         CALL vec_rthetaphi_to_xyz(Bextr,Bexttheta,Bextphi,
     &                             Bextx,Bexty,Bextz,
     &                             costheta,sintheta,cosphi,sinphi)
c
c--external currents
c
         CALL vec_rthetaphi_to_xyz(currJextr,currJexttheta,currJextphi,
     &                             currJextx,currJexty,currJextz, 
     &                             costheta,sintheta,cosphi,sinphi)
c
c--Add  J_int x B_ext
c
         currjintbextx = currJinty*Bextz - currJintz*Bexty
         currjintbexty = currJintz*Bextx - currJintx*Bextz
         currjintbextz = currJintx*Bexty - currJinty*Bextx
c
c--Add  J_ext x B_int
c
         currjextbintx = currJexty*Bintz - currJextz*Binty
         currjextbinty = currJextz*Bintx - currJextx*Bintz
         currjextbintz = currJextx*Binty - currJexty*Bintx
c
c--add boundary force to force in torus 'r' direction
c
         frtorus = frtorus + fbound(rintorus,hi,valfven2)
c
c--get J_ext x B_ext
c
         CALL vec_rthetaphi_to_xyz(frtorus,0.,0.,
     &                             fextx,fexty,fextz,
     &                             costheta,sintheta,cosphi,sinphi)
c
c--construct total external force
c         
         fextx = fextx + (currjintbextx + currjextbintx)*drhoi
         fexty = fexty + (currjintbexty + currjextbinty)*drhoi
         fextz = fextz + (currjintbextz + currjextbintz)*drhoi

c
c--next calculate the advection terms
c      
         CALL get_advectionterm(Bextr,Bexttheta,Bextphi,dBthetadr,
     &            vxi,vyi,vzi,vdotgradBx,vdotgradBy,vdotgradBz,
     &            drintorus,drcyl,costheta,sintheta,cosphi,sinphi)

      ELSE
         STOP 'ERROR: unknown string in call to externBfield'
      ENDIF
      
      RETURN
      END SUBROUTINE externBfield
      
c************************************************************
c                                                           *
c  This subroutine computes the v.grad B term               *
c  assuming only Btheta has a gradient                      *
c                                                           *
c************************************************************      
      SUBROUTINE get_advectionterm(Br,Btheta,Bphi,dBthetadr,
     &           vx,vy,vz,vdotgradBx,vdotgradBy,vdotgradBz,
     &           drintorus,drcyl,costheta,sintheta,cosphi,sinphi)
      IMPLICIT NONE
      REAL Br,Btheta,Bphi,dBthetadr,vx,vy,vz
      REAL drintorus,drcyl,costheta,sintheta,cosphi,sinphi
      REAL vdotgradBx,vdotgradBy,vdotgradBz
      
      REAL dBxdr,dBxdtheta,dBxdphi
      REAL dBydr,dBydtheta,dBydphi
      REAL dBzdr,dBzdtheta,dBzdphi
      REAL drdx,dthetadx,dphidx
      REAL drdy,dthetady,dphidy
      REAL drdz,dthetadz !,dphidz
      REAL dBxdx,dBxdy,dBxdz
      REAL dBydx,dBydy,dBydz
      REAL dBzdx,dBzdy,dBzdz
            
c
c--get gradients of cartesian components of B
c  with respect to toroidal coordinates
c
      dBxdr = -dBthetadr*sintheta*cosphi
      dBxdtheta = -Br*sintheta*cosphi - Btheta*costheta*cosphi
      dBxdphi = -Br*costheta*sinphi +Btheta*sintheta*sinphi -Bphi*cosphi
      
      dBydr = -dBthetadr*sintheta*sinphi
      dBydtheta = -Br*sintheta*sinphi - Btheta*costheta*sinphi
      dBydphi = Br*costheta*cosphi - Btheta*sintheta*cosphi -Bphi*sinphi
      
      dBzdr = dBthetadr*costheta
      dBzdtheta = Br*costheta - Btheta*sintheta
      dBzdphi = 0.

c
c--set transformation factors
c
      drdx = costheta*cosphi
      drdy = costheta*sinphi
      drdz = sintheta
      
      dthetadx = -sintheta*cosphi*drintorus
      dthetady = -sintheta*sinphi*drintorus
      dthetadz = costheta*drintorus
      
      dphidx = -sinphi*drcyl
      dphidy = cosphi*drcyl
      !dphidz = 0.

c
c--translate to get grad B in cartesians
c
      dBxdx = dBxdr*drdx + dBxdtheta*dthetadx + dBxdphi*dphidx
      dBxdy = dBxdr*drdy + dBxdtheta*dthetady + dBxdphi*dphidy
      dBxdz = dBxdr*drdz + dBxdtheta*dthetadz !+ dBxdphi*dphidz
      
      dBydx = dBydr*drdx + dBydtheta*dthetadx + dBydphi*dphidx
      dBydy = dBydr*drdy + dBydtheta*dthetady + dBydphi*dphidy
      dBydz = dBydr*drdz + dBydtheta*dthetadz !+ dBydphi*dphidz

      dBzdx = dBzdr*drdx + dBzdtheta*dthetadx + dBzdphi*dphidx
      dBzdy = dBzdr*drdy + dBzdtheta*dthetady + dBzdphi*dphidy
      dBzdz = dBzdr*drdz + dBzdtheta*dthetadz !+ dBzdphi*dphidz

c
c--get v.grad B in cartesians
c
      vdotgradBx = vx*dBxdx + vy*dBxdy + vz*dBxdz
      vdotgradBy = vx*dBydx + vy*dBydy + vz*dBydz
      vdotgradBz = vx*dBzdx + vy*dBzdy + vz*dBzdz
     
      RETURN     
      END SUBROUTINE

      SUBROUTINE vec_rthetaphi_to_xyz(vr,vtheta,vphi,vx,vy,vz,
     &           costheta,sintheta,cosphi,sinphi)
c************************************************************
c                                                           *
c  This subroutine transforms a vector in torus coordinates *
c  back to cartesian co-ordinates                           *
c                                                           *
c************************************************************      
      IMPLICIT NONE
      REAL vr,vtheta,vphi,vx,vy,vz
      REAL costheta,sintheta,cosphi,sinphi
      
      vx = vr*costheta*cosphi - vtheta*sintheta*cosphi - vphi*sinphi
      vy = vr*costheta*sinphi - vtheta*sintheta*sinphi + vphi*cosphi
      vz = vr*sintheta + vtheta*costheta
      
      RETURN           
      END SUBROUTINE vec_rthetaphi_to_xyz

      SUBROUTINE vec_xyz_to_rthetaphi(vr,vtheta,vphi,vx,vy,vz,
     &           costheta,sintheta,cosphi,sinphi)
c************************************************************
c                                                           *
c  This subroutine transforms a vector in cartesian         *
c  co-ordinates into torus coordinates                      *
c                                                           *
c  NOTE: interface is THE SAME as for the previous function *
c  so that both functions are called in exactly the same    *
c  way.                                                     *
c                                                           *
c************************************************************      
      IMPLICIT NONE      
      REAL vr,vtheta,vphi,vx,vy,vz
      REAL costheta,sintheta,cosphi,sinphi
      
      vr     =  vx*costheta*cosphi + vy*costheta*sinphi + vz*sintheta
      vtheta = -vx*sintheta*cosphi - vy*sintheta*sinphi + vz*costheta
      vphi   = -vx*sinphi + vy*cosphi
      
      RETURN
      END SUBROUTINE vec_xyz_to_rthetaphi
      
      SUBROUTINE get_factors(xi,yi,zi,costheta,sintheta,cosphi,sinphi,
     &                       rcyl,drcyl,rintorus,rintorus2,drintorus)
c************************************************************
c                                                           *
c  This subroutine deals with the coordinate part of the    *
c  transformation (returns costheta,sintheta,cosphi,sinphi) *
c                                                           *
c************************************************************      
      IMPLICIT NONE
      INCLUDE 'COMMONS/tokamak'
      
      REAL tiny
      PARAMETER (tiny=1.e-14)      
      REAL xi,yi,zi
      REAL rcyl,drcyl,rintorus2,rintorus,drintorus
      REAL sintheta,costheta,cosphi,sinphi

      rcyl = sqrt(xi**2 + yi**2)
      IF (rcyl.GT.tiny) THEN
         drcyl = 1./rcyl
      ELSE
         drcyl = 0.
      ENDIF
c rintorus is radius from centre of torus
      rintorus2 = (rcyl - Rtorus)**2 + zi**2
      rintorus = SQRT(rintorus2)
      IF (rintorus.GT.tiny) THEN
         drintorus = 1./rintorus
      ELSE
         drintorus = 0.
      ENDIF
      sintheta = zi*drintorus
      costheta = (rcyl-Rtorus)*drintorus
      cosphi = xi*drcyl
      sinphi = yi*drcyl
      
      RETURN
      END SUBROUTINE get_factors


      FUNCTION Bexternal(xcoord,ycoord,zcoord,icomp)      
c************************************************************
c                                                           *
c  This function acts as an interface to externBfield       *
c  cutting out the dummy arguments needed on the first call *
c                                                           *
c************************************************************  
      IMPLICIT NONE
      INTEGER icomp
      REAL xcoord, ycoord, zcoord, Bexternal
      REAL dumx,dumy,dumz,Bextx,Bexty,Bextz
      REAL dumh,dumrhoi
          
      dumrhoi = 1.
      CALL externBfield(xcoord,ycoord,zcoord,dumh,dumx,dumy,dumz,
     &                  dumrhoi,dumx, dumy, dumz,
     &                  dumx, dumy, dumz,
     &                  Bextx, Bexty, Bextz, 
     &                  dumx, dumy, dumz, 
     &                  dumx, dumy, dumz,'Bfield')
          
      IF (icomp.EQ.1) THEN
         Bexternal= Bextx
      ELSEIF (icomp.EQ.2) THEN
         Bexternal= Bexty
      ELSEIF (icomp.EQ.3) THEN
         Bexternal= Bextz
      ELSE
         Bexternal = 0.
         STOP 'error in Bexternal call'
      ENDIF
      
      RETURN
      
      END FUNCTION Bexternal

      SUBROUTINE fexternalB(xi,yi,zi,hi,rhoi,fextx,fexty,fextz)
c************************************************************
c                                                           *
c  This subroutine acts as an interface to externBfield     *
c  which returns only the external force component          *
c  (ie. can be called from externf for doing                *
c   hydro relaxation runs in the external tokamak potential)*
c                                                           *
c************************************************************
      IMPLICIT NONE
      REAL xi,yi,zi,hi,rhoi
      REAL dumx,dumy,dumz,fextx,fexty,fextz
      
      CALL externBfield(xi,yi,zi,hi,dumx,dumy,dumz,rhoi,
     &                  dumx, dumy, dumz,
     &                  dumx, dumy, dumz,
     &                  dumx, dumy, dumz, 
     &                  fextx, fexty, fextz, 
     &                  dumx, dumy, dumz, 'fext')
     
      RETURN
      END SUBROUTINE fexternalB
      
      FUNCTION fbound(rintorus,hi,valfven2)
c************************************************************
c                                                           *
c  This subroutine computes the boundary force              *
c  as in Monaghan (2005)                                    *
c                                                           *
c************************************************************
      IMPLICIT NONE
      REAL fbound
      REAL rintorus,hi,valfven2
      INCLUDE 'COMMONS/tokamak'
      REAL gamfac,yy,qfac
      REAL cnormk,tiny
      PARAMETER (cnormk = 1./3.1415926536)
      PARAMETER (tiny = 1.e-14)

      yy = abs(rintorus-atorus)
      qfac = yy/hi
      IF (qfac.LT.2./3.) THEN
         gamfac = cnormk*2./3.
      ELSEIF (qfac.LT.1.) THEN
         gamfac = cnormk*(2.*qfac - 1.5*qfac**2)
      ELSEIF (qfac.LT.2.) THEN
         gamfac = cnormk*(0.5*(2.-qfac)**2)
      ELSE
         gamfac = 0.
      ENDIF
      
      fbound = -0.1*valfven2*gamfac/(yy + tiny)
      IF (rintorus.gt.atorus) THEN
         WRITE(*,*) 'ERROR! particle crossed boundary'
         WRITE(*,*) 'rintorus = ',rintorus,atorus,valfven2,fbound
         STOP 'particle crossed boundary'
      ENDIF
      
      END FUNCTION fbound
