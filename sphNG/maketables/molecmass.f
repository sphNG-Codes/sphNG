      
      FUNCTION MU(rho,tm)!,delta)

		IMPLICIT NONE      
      REAL*8 rho,tm,aiony,aionx,aionz1,aionz2
      REAL*8 rhsy,cv,nh,h,IH,IHe1,me,nhe,IHe2,IH2,nH2
      REAL*8 aiony1,aiony2,delta,mu,y0,yp,aionz21,aionz22
      REAL*8 aionz11,aionz12,K,aionx1,aionx2,rhsz1,rhsz2,mH
		REAL*8 rhsx,pi,Y,X,Rg,Z

      tm=DBLE(tm)
!      INCLUDE 'COMMONS/massfrac'
!      INCLUDE 'COMMONS/units'
      Z=0.02
      Y=0.28
      X=0.700
!

!      delta=1.47572e-6
      delta=tm*1.0d-5

!      utime=1.0
!      umass=1.0
!      uerg=1.0

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

!    	PRINT *,"Density: ",rho, " Temperature: ",tm
      !Number density of molecular hydrogen
      nH2=X*rho/(2.0*mH)!*(umass)


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
         PRINT *,"MU: error in dissociation y",aiony1,aiony2
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
 !     PRINT *,"aionx",aionx1,aionx2,rhsx
      IF(aionx1.LT.0.0.AND.aionx2.GE.0.0) THEN
         aionx=aionx2
      ELSE IF(aionx2.LT.0.0.AND.aionx1.GE.0.0) THEN
         aionx=aionx1
      ELSEIF(aionx1.EQ.0.0.AND.aionx2.EQ.0.0) THEN
         aionx=0.0
      ELSE
         PRINT *,"mu: error in ionisation x",aionx1,aionx2,rhsx,tm,
     &        aiony,IH/(k*tm),IH
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

  

      !Mu 
      mu=(((2*X*(1+aiony+aionx*aiony*2.0)+
     $     Y*(1+aionz1+aionz1*aionz2))/4.0))

 !     print *,mu
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

      REAL*8 ltm,lrho,rho8,tg,mu,muu,startu,endu,cv,u,K,y1,y2
      DIMENSION muu(4602,1800)
      INTEGER I,J,R

      INCLUDE '../COMMONS/tgtbl'
      INCLUDE '../COMMONS/mutbl'
      INCLUDE '../COMMONS/units'

      REAL*4 rho

      umass = 1.0
      udist = 1.0
      udens = 1.0
      utime = 1.0
      uergg = 1.0
      uergcc = 1.0


      OPEN(UNIT=8,FILE='gasttbl',FORM='unformatted')
      DO i=1, tgmxu
         READ(8) (tgtable(i,j), j=1, tgmxrh)
      ENDDO
      CLOSE(8)

!		PRINT *,mu(0.9,10.0**14.71/getcv(0.9,10.0**14.71))
!!		PRINT *,10.0**14.71/getcv(0.9,10.0**14.71)
!		STOP

      OPEN(UNIT=10,FILE='molmasstbl',FORM='unformatted')

      DO nrho=-20000,2999,5!-20,0,0.005

         IF(MOD(nrho,500).EQ.0) PRINT *,nrho
         lrho=nrho/1000.0
         i=nrho/5+4001
         rho=10.0**lrho

	 		startu=7.725
         DO R=1,1798
            startu=startu+0.005
            u=10.0**REAL(startu)
            cv = getcv(rho,u)
            tg = u/cv
            j=R
!            print *,"Making ",I,J," as ",rho,tg,u,cv
            rho8 = rho
            muu(I,J)=log10(mu(rho8,tg))
!            print *,"Done"
         ENDDO

         IF(MOD(nrho,500).EQ.0) PRINT *,muu(I,900),10.0**muu(I,900)

      ENDDO
      PRINT *,"Doing making logmu table. Writing to disk..."
!J is specific gas energy J = 1 to 1002 rows
!I is density I = 1 to 4602 columns

      DO j=1,1799
         WRITE(10) (muu(i,j), i=1, 4601)
         PRINT *,muu(2096,j)
      END DO
!      PRINT *,cvv(20,50)
      PRINT *,"Complete"
      CLOSE(10)

      OPEN(UNIT=10,FILE='molmasstbl',FORM='unformatted')

	DO I=1,mumxu
           READ(10) (mutable(i,j), j=1,4601)
        END DO

	PRINT *,"Read table"
      K=7.725
      DO I=1,mumxu-1
         K=K+0.005
         u=10.0**(REAL(K))
         rho=1e-18			
!     print *,1e-10,u
         PRINT *,mutable(I,4600),muu(4600,I)			
         y1=get1overmu(rho,u)
         y2=getcv(rho,u)
c         WRITE(57,*)  u/y2,y2,1.0/y1
         rho = 0.9
         y1=get1overmu(rho,u)
         y2=getcv(rho,u)
c         WRITE(58,*)  u/y2,y2,1.0/y1
      ENDDO

	

      CLOSE(10)

      END 
            
            
         

      
