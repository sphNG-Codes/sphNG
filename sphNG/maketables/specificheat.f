c-------------------------------------------------------------------------
c  GENERATEU is now a simplified interface to the more general routine
c  to maintain backwards compatibility with the maketables program
c-------------------------------------------------------------------------
      SUBROUTINE GENERATEU(rho,tm,specific,uoverT,mu) !,delta)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: rho,tm
      REAL*8, INTENT(OUT) :: specific,uoverT,mu
      REAL*8 :: X,Y,Z,ne
c
c--Set metallicity and H/He fractions.  These are consistent with opacity
c     tables used.
c
      X = 0.700
      Y = 0.28
      Z = 0.02
      
      CALL ionisation(rho,tm,X,Y,Z,specific,uoverT,mu,ne)
      
      END SUBROUTINE GENERATEU

c-------------------------------------------------------------------------
c  This is the most general routine, with modified
c  interface to return the electron number density
c  and take in an assumed composition
c  (DJP 23/5/11)
c-------------------------------------------------------------------------
      SUBROUTINE ionisation(rho,tm,X,Y,Z,specific,uoverT,mu,ne)
c
c--Originally written by Stuart Whitehouse following Black & Bodenheimer 1975.
c     However, Boley et al. 2007 pointed out that Black & Bodenheimer's 
c     assumption that u=c_v*T is incorrect when c_v is not a constant and
c     this lead to an incorrect u(T) function.  Also, BB75's partition function
c     for ortho/para H_2 mix (3:1) was incorrect as T->0.
c
c     So, on 31/10/07 Matthew Bate modified this code to follow Boley et al.'s
c     equations for H_2 equation of state.
c
c     DJP (23/5/11): General cleanup, more general interface so can return
c     electron number density and change input composition, implicit none and
c     parameters set as parameters. References to units removed,
c     as everything should be in cgs anyway. Additional comments added.
c     Also incorrect value of Pi fixed!
c
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: rho,tm,X,Y,Z
      REAL*8, INTENT(OUT) :: specific,uoverT,mu,ne
      REAL*8 aionx1,aionx2,EH2,ftm
      REAL*8 aiony,aionx,aionz1,aionz2
      REAL*8 rhsy,nh,nhe,nH2,zrot
      REAL*8 aiony1,aiony2,delta,aionz21,aionz22 !,y0,yp
      REAL*8 aionz11,aionz12,rhsx,rhsz1,rhsz2,thetaRot,thetaVib
      REAL*8, PARAMETER :: pi = 3.1415926536
      REAL*8, PARAMETER :: Rg = 8.3145d7
      REAL*8, PARAMETER :: h  = 6.626d-27
      REAL*8, PARAMETER :: me = 9.109d-28
      REAL*8, PARAMETER :: k  = 1.381d-16
      REAL*8, PARAMETER :: mH = 1.6733d-24
      !1eV = 1.602176462E-12 ergs
      REAL*8, PARAMETER :: eV = 1.602176462E-12
      !1eV = 1.602176462E-12 ergs
      REAL*8, PARAMETER :: IH   = 13.5984*1.602176462d-12
      !-13.6eV ionisation potential of hydrogen
      REAL*8, PARAMETER :: IH2  = 4.73*1.602176462d-12
      !First ionisation pot of Helium
      REAL*8, PARAMETER :: IHe1 = 24.58*1.602176462d-12
      !Second one
      REAL*8, PARAMETER :: IHe2 = 54.416*1.602176462d-12
c
c--Set delta for numerical integration of equations
c
!      delta=1.47572e-6
      delta=tm*1.0d-5
c
c--Set Number density of molecular hydrogen
c
      nH2=X*rho/(2.0*mH)
c
c--Right Hand Side of B&B1975 eqn 10
c
c--This switch caused a discontinuity at log(rho) = -2.0 for lowish
c  temperatures. Left here intact for reference.

c      IF(rho.LT.1.0d-2) THEN
c         rhsy=2.11/(rho*X)*exp(-MIN(200.,52490.0/tm))
c      ELSE
c         rhsy=1.0/(nH2*h**3)*(2*pi*mH*k*tm)**1.5*(exp(-MIN(200.,IH2/
c     $        (k*tm))))
c      ENDIF

c--Instead use the equation previously reserved for rho < 0.01,
c  throughout the density space.

      rhsy=2.11/(rho*X)*exp(-MIN(200.,52490.0/tm))
c
c--Which makes degree of dissocation of molecular hydrogen "y" equal to either
c
      IF((log10(rhsy)).GT.10.0) THEN
         aiony1=-1.0!SQRT(rhsy)
         aiony2=1.0!sqrt(rhsy)
      ELSE
         aiony1=-rhsy/2.0-sqrt((4.0*rhsy+rhsy**2))/2.0
         aiony2=sqrt((4.0*rhsy+rhsy**2))/2.0-rhsy/2.0
      ENDIF

!      PRINT *,"aiony",aiony1,aiony2
      IF(aiony1.LT.0.0.AND.aiony2.GE.0.0) THEN
         aiony=aiony2
      ELSE IF(aiony2.LT.0.0.AND.aiony1.GE.0.0) THEN
         aiony=aiony1
      ELSEIF(aiony1.EQ.0.0.AND.aiony2.EQ.0.0) THEN
         aiony=0.0
      ELSE
         PRINT *,"SPECIFIC error in dissociation y",aiony1,aiony2
         STOP
      END IF
c
c--Number density of hydrogen
c
      nh=X*rho/(mH)
c
c--Number density of helium
c      
      nHe=Y*rho/(4.002*mH)

!      PRINT *,"Number density of H,He"
!      PRINT *,nh,nHe

c      PRINT *,'Degree of dissociation of H_2:',aiony

c
c--Other ionisation - ionisation of hydrogen
c
      rhsx=1.0/(aiony*nh*h**3)*(2*pi*me*k*tm)**1.5*(exp(-MIN(200.,IH/
     $     (k*tm))))
      IF((log10(rhsx)).GT.8.0) THEN
         aionx1=-1.0
         aionx2=1.0
      ELSE
         aionx1=-rhsx/2.0-sqrt((4.0*rhsx+rhsx**2))/2.0
         aionx2=sqrt((4.0*rhsx+rhsx**2))/2.0-rhsx/2.0
      ENDIF
!      PRINT *,"aionx",aionx1,aionx2,rhsx
      IF(aionx1.LT.0.0.AND.aionx2.GE.0.0) THEN
         aionx=aionx2
      ELSE IF(aionx2.LT.0.0.AND.aionx1.GE.0.0) THEN
         aionx=aionx1
      ELSEIF(aionx1.EQ.0.0.AND.aionx2.EQ.0.0) THEN
         aionx=0.0
      ELSE
         PRINT *,"SPECIFIC error in ionisation x",aionx1,aionx2
         STOP
      END IF

c      PRINT *,'Ionisation of H:',aionx
c
c--Degree of single ionisation of He
c
      rhsz1=1.0/(nhe*h**3)*(2*pi*me*k*tm)**1.5*
     $     exp(-MIN(200.,IHe1/(k*tm)))
      IF((log10(rhsz1)).GT.8.0) THEN
         aionz11=-1.0
         aionz12=1.0
      ELSE
         aionz11=-rhsz1/2.0-sqrt((4.0*rhsz1+rhsz1**2))/2.0
         aionz12=sqrt((4.0*rhsz1+rhsz1**2))/2.0-rhsz1/2.0
      ENDIF
      IF(aionz11.LT.0.0.AND.aionz12.GE.0.0) THEN
         aionz1=aionz12
      ELSE IF(aionz12.LT.0.0.AND.aionz11.GE.0.0) THEN
         aionz1=aionz11
      ELSEIF(aionz12.EQ.0.0.AND.aionz11.EQ.0.0) THEN
         aionz1=0.0
      ELSE
         PRINT *,"SPECIFIC error in ionisation z1",aionz11,aionz12
         STOP
      END IF
c
c--Degree of double ionisation of He
c
      IF(aionz1.EQ.0.0) THEN
         aionz2=0.0
      ELSE
         rhsz2=1.0/(nhe*aionz1*h**3)*(2*pi*me*k*tm)**1.5*
     $        exp(-MIN(200.,IHe2/(k*tm)))
         IF((log10(rhsz2)).GT.8.0) THEN
            aionz21=-1.0
            aionz22=1.0
         ELSE
            aionz21=-rhsz2/2.0-sqrt((4.0*rhsz2+rhsz2**2))/2.0
            aionz22=sqrt((4.0*rhsz2+rhsz2**2))/2.0-rhsz2/2.0
         ENDIF
!      PRINT *,"aionz2:",aionz21,aionz22
         IF(aionz21.LT.0.0.AND.aionz22.GE.0.0) THEN
            aionz2=aionz22
         ELSE IF(aionz22.LT.0.0.AND.aionz21.GE.0.0) THEN
            aionz2=aionz21
         ELSEIF(aionz21.EQ.0.0.AND.aionz22.EQ.0.0) THEN
            aionz2=0.0
         ELSE
            PRINT *,"SPECIFIC error in ionisation z2",aionz21,aionz22
            PRINT *,rho,tm,aionz1
            STOP
         END IF
      ENDIF

c      PRINT *,'Single & Double He ionization: ',aionz1,aionz2

c
c--Temperatures for rotational and vibrational excitation of molecular H
c     Rot: from Black & Bodenheimer 1975
c     Vib: updated to Draine et al. 1983 (from BB75's value of 6100 K)
c
      thetaVib=5987.0
      thetaRot=85.4
c
c--Assuming 3:1 ortho/para ratio of molecular hydrogen
c
!      PRINT *,tm,tm+delta,tm-delta
!      PRINT *,y0(tm+delta),y0(tm),y0(tm-delta)
!      PRINT *,yp(tm+delta),yp(tm),yp(tm-delta)
c
c--Old equations based on BB75 (that require an annoying second derivative)
c
c      y0tmpd=(tm+delta)*log(y0(tm+delta))
c      y0tm=(tm)*log(y0(tm))
c      y0tmmd=(tm-delta)*log(y0(tm-delta))
c
c      yptmpd=(tm+delta)*log(yp(tm+delta))
c      yptm=(tm)*log(yp(tm))
c      yptmmd=(tm-delta)*log(yp(tm-delta))
c  
c      dz0=((tm+delta)*log(y0(tm+delta))-tm*log(y0(tm)))/delta
c      dzp=((tm+delta)*log(yp(tm+delta))-tm*log(yp(tm)))/delta
c
c      dz0m=((tm)*log(y0(tm))-(tm-delta)*log(y0(tm-delta)))/delta
c      dzpm=((tm)*log(yp(tm))-(tm-delta)*log(yp(tm-delta)))/delta
c
c      d2z0=(y0tmpd-2.0*y0tm+y0tmmd)/delta**2
c      d2zp=(yptmpd-2.0*yptm+yptmmd)/delta**2
c
c      d2z0=(dz0-dz0m)/delta
c      d2zp=(dzp-dzpm)/delta
c
c
c!      PRINT *,"z0,zp: ",y0(tm),yp(tm)
c!      PRINT *,"First derivatives, z0, zp:",dz0,dzp
c!      PRINT *,"Second derivatives, z0, zp:",d2z0,d2zp
c
c      ftm1=yp(tm)/(yp(tm)+3.0*y0(tm))
c      ftm2=3.0*y0(tm)/(yp(tm)+3.0*y0(tm))
c
c!      PRINT *,"ftm12",ftm1,ftm2
c
c      ftm=tm*(ftm1*d2zp+ftm2*d2z0)
c
c--New equations for 2nd term in Boley et al. 2007
c     which is T^2 d(ln z)/dT
c
      ftm = tm**2*(LOG(zrot(tm+delta,thetaRot))
     $            -LOG(zrot(tm-delta,thetaRot)))/(2.0*delta)

!      PRINT *,"f(tm):",ftm
c
c--Definition of specific internal energy of H_2 (following Boley et al., eqn. 5)
c     Note that the derivative of this with respect to T gives the equation
c     in BB75 because BB75 has the right equation for c_v=du/dT but then
c     an incorrect equation u=c_v*T.
c     The equation (11) in BB75 is incorrect in that the second term has
c     an extra term: (thetaRot/tm)**2 which needs to be removed.
c
      EH2=Rg/2.0*(1.5*tm + ftm + 
     $     thetaVib/(exp(MIN(500.,thetaVib/tm))-1.0))
c
c--INCORRECT OLD EQUATIONS (for incorrect version of cv for H_2)
c
c      EH2=Rg/2.0*(1.5+ftm+(thetaVib/tm)**2*
c     $     exp(MIN(500.,thetaVib/tm))/(exp(MIN(500.,thetaVib/tm))
c     $     -1.0)**2)

!      PRINT *,"EH2: ",EH2
!      PRINT *,"plot quant: ",EH2*2.0/Rg
!      WRITE(66,*) tm,EH2*2.0/Rg,ftm
c
c--Mu
c
      mu=(1.0/((2.*X*(1.+aiony+aionx*aiony*2.0)+
     $     Y*(1.+aionz1+aionz1*aionz2))/4.0))
c
c--Final value of specific internal energy
c
      specific = X*(1.0-aiony)*EH2 + (1.5*X*(1 + aionx)*aiony+0.375*Y*
     $     (1.0 + aionz1 + aionz1*aionz2))*Rg*Tm
     $     + X*(1.304e13*aionx + 2.143e12)*aiony
     $     + Y*(5.888e12*(1.0-aionz2) + 1.892e13*aionz2)*aionz1
c
c--INCORRECT OLD EQUATIONS
c
c      cv=((X*(1.0-aiony)*EH2/Rg+1.5*X*(1+aionx)*aiony+0.375*Y*
c    $     (1.0+aionz1+aionz1*aionz2))*Rg*Tm+X*(1.304e13*aionx+2.143e12)
c     $     *aiony+Y*(5.888e12*(1.0-aionz2)+1.892e13*aionz2)*aionz1)/tm
c
c--Set ratio of specific internal energy to temperature (note that this is
c     NOT equal to c_v because c_v is not a constant.  c_v = du/dT.
c
      uoverT = specific/tm
c
c--Return number density of electrons
c  from the mean molecular weight, we have n = rho/m_h * 1/mu
c  where 1/mu = X/2*(1 + aiony + 2*aiony*aionx)
c             + Y/4*(1 + aionz1 + aionz1*aionz2) 
c  so, noting that n = ni + ne, and hence 1/mu = 1/mu_i + 1/mu_e
c  ne are just the electron parts making up 1/mu_e
c
      ne = nh*(aionx*aiony*2.)
     $   + nHe*(aionz1 + aionz1*aionz2)

      IF(uoverT.LE.1) THEN
         PRINT *,(X*(1.0-aiony)*EH2)
         PRINT *,1.5*X*(1+aionx)*aiony+0.375*Y*(1.0+aionz1+aionz1*
     $        aionz2)*Rg
         PRINT *,X*(1.304e13*aionx+2.143e12)/tm*aiony
         PRINT *,Y*(5.888e12*(1.0-aionz2)+1.892e13*aionz2)*aionz1/tm
 109     FORMAT(9((1pd12.5),1x))
         WRITE(80,109) tm,uoverT,aionx,aiony,aionz1,aionz2,rhsy,
     $        uoverT/Rg,EH2/Rg*mu
         PRINT *,uoverT
         STOP 
      ENDIF

      RETURN
      END SUBROUTINE

c-------------------------------------------------------------------------
c Overall partition function assuming a fixed 3:1 ortho/para ratio (b:a)
c from Boley et al. 2007:
c     zrot = zp^(a/(a+b))*zodash^(b/(a+b))
c     where zodash = zortho*exp(2*thetarot/T)
c     Here a and b are hardwired to 1 and 3 as above.
c-------------------------------------------------------------------------
      FUNCTION zrot(tm,theta_rot)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: tm,theta_rot
      REAL*8 zrot,ypara,yortho

      zrot = ypara(tm,theta_rot)**0.25*(yortho(tm,theta_rot)
     $        *EXP(2.0*theta_rot/tm))**0.75

      RETURN
      END

c-------------------------------------------------------------------------
c Partition function for para-Hydrogen
c (see Boley et al. 2007)
c-------------------------------------------------------------------------
      FUNCTION ypara(tm,thetaRot)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: tm,thetaRot
      REAL*8 ypara,zpara !,zpold
      INTEGER I,J

      zpara=0.0
      !zpold=-9999.9
      DO I=2,200,2
         J=I-2
         !zpold=zp
         zpara=zpara+(2*J+1)*exp(-J*(J+1)*thetaRot/tm)
!         PRINT *,zpara
   !      IF(ABS((zpold-zp)/zp).GE.0.99999.AND.ABS((zpold-zp)/zp).LE.
   !  $        1.000001) THEN
   !         GOTO 101
   !      END IF
      END DO
!      STOP
! 100  PRINT *,"zp",tm,zp,I
! 101  ypara=tm*log(zpara)
! 101  CONTINUE
      ypara=(zpara)

      RETURN
      END

c-------------------------------------------------------------------------
c Partition function for ortho-Hydrogen
c (see Boley et al. 2007)
c-------------------------------------------------------------------------
      FUNCTION yortho(tm,thetaRot)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: tm,thetaRot
      REAL*8 yortho,zortho !,z0old
      INTEGER I,J
      !z0old=-9999.9
      zortho=0.0
      DO I=1,199,2
         J=I
         !z0old=zortho
         zortho=zortho+(2*J+1)*exp(-J*(J+1)*thetaRot/tm)
!         PRINT *,zortho
  !       IF(ABS((z0old-z0)/z0).GE.0.99999.AND.ABS((z0old-z0)/z0).LE.
  !   $        1.000001) THEN
  !          GOTO 101
  !       END IF
      END DO
!      STOP
! 100  PRINT *,"z0",tm,z0,I
! 101  y0=tm*log(z0)
! 101  CONTINUE
      yortho=(zortho)

      RETURN
      END
