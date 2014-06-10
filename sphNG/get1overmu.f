      FUNCTION GET1OVERMU(rho2,u2)

!      IMPLICIT NONE

      REAL*4 rho2
      REAL rho,u2,u,lrho,lu,rhoval1,rhoval2,uval1,uval2
      REAL y1,y2,y3,y4,w,v,rmu,get1overmu
      INTEGER nkrho1,nkrho2,nku1,nku2

      INCLUDE 'COMMONS/eostbl'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      
c      get1overmu = 1.0/gmw     
c      RETURN
      
      rho=rho2*umass/udist**3!udens
      u=u2*uergg


      lrho = log10(rho)
c
c--If goes to very high density, use last values in table
c
      IF (lrho.GT.3.0) lrho = 3.0

      lu=log10(u)

      nkrho1=INT(lrho/0.005)+eostbl_rho1
      IF(lrho.LT.0.0) nkrho1=nkrho1-1
      IF(nkrho1.GE.mumxrh) nkrho1 = mumxrh - 1
      nku1=INT(lu/0.005)-1545
      IF(nkrho1.LT.1.0) nkrho1=1
      IF(nku1.LT.1.0) nku1=1
      
      nkrho2=nkrho1+1
      nku2=nku1+1

      rhoval1=(nkrho1-eostbl_rho1)*0.005
      rhoval2=(nkrho2-eostbl_rho1)*0.005
      uval1=(nku1+1545)*0.005
      uval2=(nku2+1545)*0.005



c      IF(nkrho1.LT.1.0.OR.nkrho1.GE.mumxrh.OR.nku1.LT.1)
      IF(nkrho1.LT.1.0.OR.nku1.LT.1)
     $     CALL FAILED3(0,lu,lrho,nkrho1,nkrho2,nku1,nku2)
      
      IF(nku2.GE.mumxu) THEN
         get1overmu = 1.61
         RETURN
      ENDIF

!     Bilinear interpolation from Numerical Recipes
!     Interpolation is in log T_g
       
 
      !Interpolation in log10 space
         y1=(mutable(nku1,nkrho1))
         y2=(mutable(nku2,nkrho1))
         y3=(mutable(nku2,nkrho2))
         y4=(mutable(nku1,nkrho2))
       
       
         w=(lrho-rhoval1)/(rhoval2-rhoval1)
         v=(lu-uval1)/(uval2-uval1)

         !     Final value of log10 temperature
         rmu=(1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4

         get1overmu=(10.0**rmu)


         IF(get1overmu.EQ.0.00) CALL FAILED3(5,lu,lrho,nkrho1,nkrho2,
     $        nku1,nku2)

         END

      SUBROUTINE FAILED3(i,ltg,lrho,nkrho1,nkrho2,nktg1,nktg2)
 
      INTEGER i
      REAL ltg,lrho
  
 
      IF(i.EQ.0) THEN
         PRINT *,"GET1OVERMU: Error. Estimated array locations"
         PRINT *,"outside bounds. Ending."
         PRINT *,"Log gas specific energy is ",ltg,10.0**ltg
         PRINT *,"Log density is ",lrho,10.0**lrho
         PRINT *,nkrho1,nkrho2,nktg1,nktg2
         PRINT *,i
      ELSE
         PRINT *,"GET1OVERMU: Molecular mass returned is zero!"
         PRINT *,"Log gas temperature is ",ltg
         PRINT *,"Log density is ",lrho
      END IF
       

      STOP
      END
      
