      SUBROUTINE GENERATEU(rho,tm,specific,uoverT,mu) !,delta)
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
      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 rho,tm,aiony,aionx,aionz1,aionz2
      REAL*8 rhsy,specific,nh,h,IH,IHe1,k,me,nhe,IHe2,mu,IH2,nH2
      REAL*8 aiony1,aiony2,delta,y0,yp,aionz21,aionz22
      REAL*8 aionz11,aionz12,mH

c
c--Set metallicity and H/He fractions.  These are consistent with opacity
c     tables used.
c
      Z=0.02
      Y=0.28
      X=0.700
c      Z=0.0
c      Y=0.0
c      X=1.0
c
c--Set delta for numerical integration of equations
c
!      delta=1.47572e-6
      delta=tm*1.0d-5
c
c--Calculations done in cgs units
c
      utime=1.0
      umass=1.0
      uerg=1.0
c
c--Set other constants
c
      h=6.626d-27!/uerg/utime
      me=9.109d-28!/umass
      k=1.381d-16!/uerg
      pi=3.14259265
      Rg=8.3145d7
      mH=1.6733d-24

      !1eV = 1.602176462E-12 ergs
      IH=13.5984*1.602176462d-12!/uerg
      !-13.6eV ionisation potential of hydrogen in code units
		IH2=4.73*1.602176462d-12!/uerg
      
      IHe1=24.58*1.602176462d-12!/uerg  !First ionisation pot of Helium
      IHe2=54.416*1.602176462d-12!/uerg !Second one
c
c--Set Number density of molecular hydrogen
c
      nH2=X*rho/(2.0*mH)!*(umass)
c
c--Right Hand Side of B&B1975 eqn 10
c
      IF(rho.LT.1.0d-2) THEN
         rhsy=2.11/(rho*X)*exp(-MIN(200.,52490.0/tm))
      ELSE
         rhsy=1.0/(nH2*h**3)*(2*pi*mH*k*tm)**1.5*(exp(-MIN(200.,IH2/
     $        (k*tm))))
      ENDIF
!      PRINT *,rhsy
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
      nh=X*rho/(1.6733d-24)!*(umass)
c
c--Number density of helium
c      
      nHe=Y*rho/(4.002*1.6733d-24)

!      PRINT *,"Number density of H,He"
!      PRINT *,nh,nHe

!      PRINT *,"Degree of disociation of H_2:",aiony

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


!      PRINT *,"Ionisation of H:",aionx

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
c--Degree of double ionisatoin of He
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

!      PRINT *,"Single & Double He ionization: ",aionz1,aionz2

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
      ftm = tm**2*(LOG(zrot(tm+delta)) -LOG(zrot(tm-delta)))/(2.0*delta)

!      PRINT *,"f(tm):",ftm
c
c--Definition of specific internal energy of H_2 (following Boley et al.)
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
      mu=(1.0/((2*X*(1+aiony+aionx*aiony*2.0)+
     $     Y*(1+aionz1+aionz1*aionz2))/4.0))
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
      END

c-------------------------------------------------------------------------

      FUNCTION zrot(tm)

      IMPLICIT REAL*8 (A-H,O-Z)

      theta_rot = 85.4

      zrot = yp(tm)**0.25*(y0(tm)*EXP(2.0*theta_rot/tm))**0.75

      RETURN
      END

c-------------------------------------------------------------------------

      FUNCTION yp(tm)

      REAL*8 yp,tm,zp,thetaRot,zpold
      INTEGER I

      thetaRot=85.4
      zp=0.0
      zpold=-9999.9
      DO I=2,200,2
         J=I-2
         zpold=zp
         zp=zp+(2*J+1)*exp(-J*(J+1)*thetaRot/tm)
!         PRINT *,zp
   !      IF(ABS((zpold-zp)/zp).GE.0.99999.AND.ABS((zpold-zp)/zp).LE.
   !  $        1.000001) THEN
   !         GOTO 101
   !      END IF
      END DO
!      STOP
! 100  PRINT *,"zp",tm,zp,I
! 101  yp=tm*log(zp)
 101  yp=(zp)

      RETURN
      END

c-------------------------------------------------------------------------

      FUNCTION y0(tm)

      REAL*8 y0,tm,z0,thetaRot,z0old
      INTEGER I
      z0old=-9999.9
      thetaRot=85.4
      z0=0.0
      DO I=1,199,2
         J=I
         z0old=z0
         z0=z0+(2*J+1)*exp(-J*(J+1)*thetaRot/tm)
!         PRINT *,z0
  !       IF(ABS((z0old-z0)/z0).GE.0.99999.AND.ABS((z0old-z0)/z0).LE.
  !   $        1.000001) THEN
  !          GOTO 101
  !       END IF
      END DO
!      STOP
! 100  PRINT *,"z0",tm,z0,I
! 101  y0=tm*log(z0)
 101  y0=(z0)

      RETURN
      END

c-------------------------------------------------------------------------

      PROGRAM MAKETABLE

      REAL*8 ltm,lrho,tm,rho,specific,uoverT,cvv,muu,mu,Rg
      DIMENSION cvv(4602,1802),muu(4602,1802)
      INTEGER I,J

!		PRINT *,"Rho, Tg"
!		READ *,rho,tm
!		CALL GENERATECV(rho,tm,cv,mu)
!		PRINT *,"CV: ",cv,cv/8.3144d7
!		PRINT *,"Mu: ",mu
!		STOP

      Rg=8.3145d7

!     OPEN(UNIT=55,FILE='specheattbl')

      OPEN(UNIT=8,FILE='specheattbl',FORM='unformatted')
!      OPEN(UNIT=10,FILE='molmasstbl',FORM='unformatted')

      DO nrho=-20000,3000,5!-20,0,0.005

         IF(MOD(nrho,500).EQ.0) PRINT *,nrho
         lrho=nrho/1000.0
         i=nrho/5+4001
         rho=10.0**lrho

         DO ntm=0,9000,5!0.05,4,0.005
            ltm=ntm/1000.0
            j=ntm/5+1
            tm=10.0**ltm
            CALL GENERATEU(rho,tm,specific,uoverT,mu)
            write (95,99001) tm,specific/Rg,uoverT/Rg,mu
99001       FORMAT(4(1PE12.5,1X))
            cvv(I,J)=log10(uoverT)
            muu(I,J)=log10(mu)
         END DO

c         STOP

         IF(MOD(nrho,500).EQ.0) PRINT *,cvv(I,900)

      ENDDO
      PRINT *,"Doing making logcv table. Writing to disk..."
!J is temperature J = 1 to 1002 rows
!I is density I = 1 to 4602 columns

      DO j=1,1802
         WRITE(8) (cvv(i,j), i=1, 4601)
!         WRITE(10) (muu(i,j), i=1, 4601)
         PRINT *,cvv(2096,j),muu(2096,j)
         WRITE(66,*) cvv(4600,j)
      END DO
!      PRINT *,cvv(20,50)
      PRINT *,"Complete"
      CLOSE(55)

      END 
            
c=========================================================================
