      
      SUBROUTINE GENERATECV(rho,tm,cv,mu)!,delta)
      
      REAL*8 rho,tm,aiony,aionx,aionz1,aionz2
      REAL*8 rhsy,cv,nh,h,IH,IHe1,k,me,nhe,IHe2,mu,IH2,nH2
      REAL*8 aiony1,aiony2,delta,mu,y0,yp,aionz21,aionz22
      REAL*8 aionz11,aionz12,mH

      tm=DBLE(tm)
!      INCLUDE 'COMMONS/massfrac'
!      INCLUDE 'COMMONS/units'
      Z=0.02
      Y=0.28
      X=0.700
!

!      delta=1.47572e-6
      delta=tm*1.0d-5

      utime=1.0
      umass=1.0
      uerg=1.0

      h=6.626d-27!/uerg/utime
      me=9.109d-28!/umass
      k=1.381d-16!/uerg
      pi=3.142
      Rg=8.3145d7
		mH=1.6733d-24

      !1eV = 1.602176462E-12 ergs
      IH=13.5984*1.602176462d-12!/uerg
      !-13.6eV ionisation potential of hydrogen in code units
		IH2=4.73*1.602176462d-12!/uerg
      
      IHe1=24.58*1.602176462d-12!/uerg  !First ionisation pot of Helium
      IHe2=54.416*1.602176462d-12!/uerg !Second one

      !Number density of molecular hydrogen
      nH2=X*rho/(2.0*mH)!*(umass)
    

      !Right Hand Side of B&B1975 eqn 10
		IF(rho.LT.1.0d-2) THEN
      	rhsy=2.11/(rho*X)*exp(-MIN(200.,52490.0/tm))
		ELSE
      	rhsy=1.0/(nH2*h**3)*(2*pi*mH*k*tm)**1.5*(exp(-MIN(200.,IH2/
     $     (k*tm))))
		ENDIF
!      PRINT *,rhsy
      !Which makes degree of dissocation of molecular hydrogen 
      !"y" equal to either
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
         PRINT *,"CV error in dissociation y",aiony1,aiony2
         STOP
      END IF
      

      !Number density of hydrogen
      nh=X*rho/(1.6733d-24)!*(umass)
      
      !And similarly of Helium
      nHe=Y*rho/(4.002*1.6733d-24)

!      PRINT *,"Number density of H,He"
!      PRINT *,nh,nHe

!      PRINT *,"Degree of disociation of H_2:",aiony


      !Other ionisation - ionisation of hydrogen
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
         PRINT *,"CV error in ionisation x",aionx1,aionx2
         STOP
      END IF


!      PRINT *,"Ionisation of H:",aionx

      
      !Degree of single ionisation of He
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
         PRINT *,"CV error in ionisation z1",aionz11,aionz12
         STOP
      END IF


      !Degree of double ionisatoin of He
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
        PRINT *,"CV error in ionisation z2",aionz21,aionz22
          PRINT *,rho,tm,aionz1
         STOP
      END IF
      ENDIF


!      PRINT *,"Single & Double He ionization: ",aionz1,aionz2

  

      !Temperatures
      thetaVib=6100.0
      thetaRot=85.4

      !assuming eqm mix of ortho- and para-hydrogen
!      PRINT *,tm,tm+delta,tm-delta
!      PRINT *,y0(tm+delta),y0(tm),y0(tm-delta)
!      PRINT *,yp(tm+delta),yp(tm),yp(tm-delta)
      !that annoying second derivative

      y0tmpd=(tm+delta)*log(y0(tm+delta))
      y0tm=(tm)*log(y0(tm))
      y0tmmd=(tm-delta)*log(y0(tm-delta))

      yptmpd=(tm+delta)*log(yp(tm+delta))
      yptm=(tm)*log(yp(tm))
      yptmmd=(tm-delta)*log(yp(tm-delta))

  
      dz0=((tm+delta)*log(y0(tm+delta))-tm*log(y0(tm)))/delta
      dzp=((tm+delta)*log(yp(tm+delta))-tm*log(yp(tm)))/delta

      dz0m=((tm)*log(y0(tm))-(tm-delta)*log(y0(tm-delta)))/delta
      dzpm=((tm)*log(yp(tm))-(tm-delta)*log(yp(tm-delta)))/delta

      d2z0=(y0tmpd-2.0*y0tm+y0tmmd)/delta**2
      d2zp=(yptmpd-2.0*yptm+yptmmd)/delta**2

      d2z0=(dz0-dz0m)/delta
      d2zp=(dzp-dzpm)/delta


!      d2z0=(y0(tm+delta)-2.0*y0(tm)+y0(tm-delta))/delta**2
!      d2zp=(yp(tm+delta)-2.0*yp(tm)+yp(tm-delta))/delta**2


!      PRINT *,"z0,zp: ",y0(tm),yp(tm)
!      PRINT *,"First derivatives, z0, zp:",dz0,dzp
!      PRINT *,"Second derivatives, z0, zp:",d2z0,d2zp

      ftm1=yp(tm)/(yp(tm)+3.0*y0(tm))
      ftm2=3.0*y0(tm)/(yp(tm)+3.0*y0(tm))

!      PRINT *,"ftm12",ftm1,ftm2

      ftm=tm*(ftm1*d2zp+ftm2*d2z0)

!      PRINT *,"f(tm):",ftm
!      PRINT *,"Term 2,3:",ftm,(thetaVib/tm)**2*
!     $     exp(thetaVib/tm)/(exp(thetaVib/tm)-1.0)**2

!(thetaRot/tm)**2*
      EH2=Rg/2.0*(1.5+ftm+(thetaVib/tm)**2*
     $     exp(MIN(500.,thetaVib/tm))/(exp(MIN(500.,thetaVib/tm))
     $     -1.0)**2)

!      PRINT *,"EH2: ",EH2
!      PRINT *,"plot quant: ",EH2*2.0/Rg
!      WRITE(66,*) tm,EH2*2.0/Rg,ftm

      !Mu 
      mu=(1.0/((2*X*(1+aiony+aionx*aiony*2.0)+
     $     Y*(1+aionz1+aionz1*aionz2))/4.0))

      !Final value

      cv=((X*(1.0-aiony)*EH2/Rg+1.5*X*(1+aionx)*aiony+0.375*Y*
     $     (1.0+aionz1+aionz1*aionz2))*Rg*Tm+X*(1.304e13*aionx+2.143e12)
     $     *aiony+Y*(5.888e12*(1.0-aionz2)+1.892e13*aionz2)*aionz1)/tm
	
		IF(cv.LE.1) THEN

			
		PRINT *,(X*(1.0-aiony)*EH2)
     	 PRINT *,1.5*X*(1+aionx)*aiony+0.375*Y*(1.0+aionz1+aionz1*
     $aionz2)*Rg
		PRINT *,X*(1.304e13*aionx+2.143e12)/tm*aiony
      PRINT *,Y*(5.888e12*(1.0-aionz2)+1.892e13*aionz2)*aionz1/tm
 109  FORMAT(9((1pd12.5),1x))
      WRITE(80,109) tm,cv,aionx,aiony,aionz1,aionz2,rhsy,cv/Rg,EH2/Rg*mu
      PRINT *,cv
		STOP 
		ENDIF
!      cv=((X*(1-aiony)*EH2/Rg+1.5*X*(1+aionx)*aiony+0.375*Y*(1+aionz1+aionz1*aionz2))*Rg+X*(1.304e13*aionx+2.143e12)*aiony+Y*(5.888e12*(1-aionz2)+1.892e13*aionz2)*aionz1/tm) !Back-up copy without newlines

      END

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
      END

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
      END




      PROGRAM MAKETABLE

      REAL*8 ltm,lrho,tm,rho,cv,cvv,muu,mu
      DIMENSION cvv(4602,1802),muu(4602,1802)
      INTEGER I,J

!		PRINT *,"Rho, Tg"
!		READ *,rho,tm
!		CALL GENERATECV(rho,tm,cv,mu)
!		PRINT *,"CV: ",cv,cv/8.3144d7
!		PRINT *,"Mu: ",mu
!		STOP

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
				CALL GENERATECV(rho,tm,cv,mu)
            cvv(I,J)=log10(cv)
				muu(I,J)=log10(mu)
         ENDDO
			
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
            
            
         

      
