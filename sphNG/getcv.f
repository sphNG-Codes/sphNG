      REAL FUNCTION GETCV(rho2,u2)

      IMPLICIT NONE

      REAL*4 rho2
      REAL rho,u2,u,lrho,lu,rhoval1,rhoval2,uval1,uval2
      REAL y1,y2,y3,y4,w,v,rtg
      INTEGER nkrho1,nkrho2,nku1,nku2

      INCLUDE 'COMMONS/eostbl'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'

c      getcv = 1.5*Rg/uergg     
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
      IF(nkrho1.GE.tgmxrh) nkrho1 = tgmxrh - 1
      nku1=INT(lu/0.005)-1545
      IF(nkrho1.LT.1) nkrho1=1
      IF(nku1.LT.1) nku1=1
      
      nkrho2=nkrho1+1
      nku2=nku1+1

      rhoval1=(nkrho1-eostbl_rho1)*0.005
      rhoval2=(nkrho2-eostbl_rho1)*0.005
      uval1=(nku1+1545)*0.005
      uval2=(nku2+1545)*0.005

c      IF(nkrho1.LT.1.0.OR.nkrho1.GE.tgmxrh.OR.nku1.LT.1)
      IF(nkrho1.LT.1.0.OR.nku1.LT.1)
     $     CALL FAILED2(0,lu,lrho,nkrho1,nkrho2,nku1,nku2)
      
      IF(nku2.GE.tgmxu) THEN
!It's reached the stage where c_v is constant for all rho at high T
         getcv = 198561558.8045447/uergg
         RETURN
      ENDIF

!     Bilinear interpolation from Numerical Recipes
!     Interpolation is in log T_g

      !Interpolation in log10 space
         y1=(tgtable(nku1,nkrho1))
         y2=(tgtable(nku2,nkrho1))
         y3=(tgtable(nku2,nkrho2))
         y4=(tgtable(nku1,nkrho2))
       
       
         w=(lrho-rhoval1)/(rhoval2-rhoval1)
         v=(lu-uval1)/(uval2-uval1)

         !     Final value of log10 temperature
         rtg=(1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4
!                  PRINT *,u,rtg
!                  PRINT *,y1,y2,y3,y4
         getcv=u/(10.0**rtg)/uergg
      

         IF(getcv.EQ.0.00) CALL FAILED2(5,lu,lrho,nkrho1,nkrho2,
     $        nku1,nku2)

         END

      SUBROUTINE FAILED2(i,ltg,lrho,nkrho1,nkrho2,nktg1,nktg2)
 
      INTEGER i
      REAL ltg,lrho
  
 
      IF(i.EQ.0) THEN
         PRINT *,"GETCV: Error. Estimated array locations"
         PRINT *,"outside bounds. Ending."
         PRINT *,"Log gas specific energy is ",ltg,10.0**ltg
         PRINT *,"Log density is ",lrho,10.0**lrho
         PRINT *,nkrho1,nkrho2,nktg1,nktg2
      ELSE
         PRINT *,"GETCV: Specific heat capacity returned is zero!"
         PRINT *,"Log gas temperature is ",ltg
         PRINT *,"Log density is ",lrho
      END IF
       

      CALL quit(1) 
      END
      
