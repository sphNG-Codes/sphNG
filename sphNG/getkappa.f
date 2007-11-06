      FUNCTION GETKAPPA(u2,cv2,rho2)
  
!      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/optbl'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/rbnd'

      COMMON /getkap/ iflag

      REAL*4 rho2
      REAL lTg,lrho,getkappa,u2,u,cv,cv2,rho
      REAL rkappa,y1,y2,y3,y4,w,v
      INTEGER nkrho1,nkrho2,nktg1,nktg2,iflag
      LOGICAL y1alex,y2alex,y3alex,y4alex,tryalex

c     SHOCK!!!!!!!!!
!      getkappa=4.0d2*umass/udist**2
!      RETURN

      nkrho1=2
      nkrho2=3
      nktg1=2
      nktg2=3

c      write (*,*) 'Enter kap'

      u=u2*uergg
      rho=rho2*umass/udist**3
      cv=cv2*uergg

      IF(u.EQ.0.0.OR.rho.EQ.0.0.OR.cv.EQ.0.0) CALL FAILED(7,ltg,lrho)

!     Call for each particle
!     rkappa (intent=out) is opacity of particle 
!     u (intent=in) is specific gas energy
!     cv (intent=in) is specfic heat capacity
!     rho (intent=in) is density

!     optable is the table of opacities in the common block optbl
!     Organised by log10 gas temperature (first array index) and
!     log10 density (second array index).
!     Actual values start at 2, and dimensions are optable(56:9) 
!     so table runs from optable(2:2) to optable(56:9)


!     optable(1:x) is log10 temperature
!     optable(x:1) is log10 density
      tryalex=.FALSE.

      ltg = log10(u/cv)
      lrho = log10(rho)

      !First get array positions for density
      nkrho1=MIN(INT(lrho)+15,15)

!     This bit extrapolates the table to lower densities as being the same
! as 1d-14

      IF(nkrho1.LE.1) THEN
         nkrho1=2
         IF (ltg.LT.3.35) lrho = -14
      ENDIF

      nkrho2=nkrho1+1
 
      IF(u/cv.GE.10.0.AND.u/cv.LT.500.0) THEN
         nktg1=(INT(0.1*u/cv))+2
      ELSEIF(u/cv.Ge.500.0.AND.u/cv.LT.1500.0) THEN
         nktg1=INT(0.02*u/cv)+42
      ELSEIF(u/cv.GE.1500.AND.u/cv.LT.2500) THEN 
         nktg1=INT(0.01*u/cv)+57
      ELSEIF(u/cv.GE.1.0.AND.u/cv.LT.10.0) THEN
         nktg1=2
      ELSEIF(u/cv.GE.2500.0.AND.u/cv.Le.10000.0) THEN
         nktg1=INT(20*ltg)+15
      ELSEIF(u/cv.GT.10000.0 .AND. lrho.GE.-14.0) THEN
         val1l = optable(94,nkrho1)
         val1h = optable(94,nkrho2)
         val2l = optable(95,nkrho1)
         val2h = optable(95,nkrho2)
         
         w=(lrho-optable(1,nkrho1))/
     $        (optable(1,nkrho2)-optable(1,nkrho1))

         val1 = (1.0 - w) * val1l + w * val1h
         val2 = (1.0 - w) * val2l + w * val2h
c         IF(val2.GT.(1.2512E+22* rho * (10000.0)**(-3.5))) THEN
c            getkappa = 0.4*umass/udist**2
c            RETURN
c!     PRINT *,"GETKAPPA: Warning, low density-high temperature"
c!     PRINT *,"How to join tables together?"
c!     STOP
c         ENDIF

         grad = (val2 - val1)/0.05
         valup = (ltg - 4.0) * grad + val2
!Kramer's opacity + electron scattering (0.4 cm^2/g, according to 
!page two of astro-ph/0410343
         valdown = (1.2512E+22* rho *(10.0**ltg)**(-3.5)) + 0.4
!     PRINT *,valup,valdown, (10.0**ltg)**(-3.5)
         getkappa = MIN(valup,valdown)*umass/udist**2
c      write (*,*) 'Exit kap'
         RETURN
      ENDIF
      
      nktg2=nktg1+1

      IF(nkrho1.LE.1.OR.nkrho1.GE.opmxrh.OR.nktg1.LE.1.OR.
     $     nktg1.GE.opmxtg)
     $     CALL FAILED(0,ltg,lrho)

!     Bilinear interpolation from Numerical Recipes
!     Pollack table gives kappa, Alexander gives log10 kappa
!     Interpolation is in log kappa
      
      IF(optable(nktg1,nkrho1).EQ.0.00) CALL FAILED(1,ltg,lrho)
      IF(optable(nktg2,nkrho1).EQ.0.00) CALL FAILED(2,ltg,lrho)
      IF(optable(nktg2,nkrho2).EQ.0.00) CALL FAILED(3,ltg,lrho)
      IF(optable(nktg1,nkrho2).EQ.0.00) CALL FAILED(4,ltg,lrho)

!TEMP CODE
!      IF(optable(nktg1,nkrho1).EQ.0.00.OR.optable(nktg2,nkrho1).EQ.0.00.
!     $     OR.optable(nktg2,nkrho2).EQ.0.00.OR.
!     $     optable(nktg1,nkrho2).EQ.0.00) THEN
!         getkappa=0.0
!         RETURN
!      ENDIF

      !Interpolation in log10 space
      y1=log10(optable(nktg1,nkrho1))
      y2=log10(optable(nktg2,nkrho1))
      y3=log10(optable(nktg2,nkrho2))
      y4=log10(optable(nktg1,nkrho2))

      w=(lrho-optable(1,nkrho1))/
     $     (optable(1,nkrho2)-optable(1,nkrho1))
      v=(ltg-optable(nktg1,1))/
     $     (optable(nktg2,1)-optable(nktg1,1))

!     Final value of log10 opacity
      rkappa=(1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4

c
c--For very low densities and temperatures > 2000K (use Bell & Lin 1994)
c
      IF (lrho.LT.-14 .AND. ltg.GT.3.35) THEN
         rkappa = MAX(-8.0 + lrho*2.0/3.0 + ltg*3.0,
     &        -36.0 + lrho/3.0 + ltg*10.0)
         rkappa = MIN(rkappa,LOG10(0.4))
      ENDIF
      
!     And the actual opacity returned
      rkappa=(10.0**rkappa)*umass/udist**2
!     PRINT *,getkappa,umass,udist,rkappa,v,w


      getkappa=rkappa!/10.0
!      getkappa=MAX(rkappa,1.0/(10.0*rmax*rho2))
!      getkappa=MAX(rkappa,1.0/(rmax*4.0*0.01))
!      getkappa=15.0*umass/udist**2


      IF(getkappa.EQ.0.00) CALL FAILED(5,ltg,lrho)

c      WRITE (*,*) 'Exit kap'
      RETURN

      END !SUBROUTINE GETKAPPA
      



      SUBROUTINE FAILED(i,ltg,lrho)

      INTEGER i
      REAL ltg,lrho
 

      IF(i.GE.1.AND.i.LE.4 .OR. i.EQ.6) THEN
         PRINT *,"GETKAPPA: Warning, out of range values of"
         PRINT *,"temperature/density in Pollack and Alexander tables."
         PRINT *,"Log gas temperature is ",ltg
         PRINT *,"Log density is ",lrho
         PRINT *,"Trying to interpolate in point ",i
      ELSEIF(i.EQ.0) THEN
         PRINT *,"GETKAPPA: Error. Estimated array locations"
         PRINT *,"outside bounds. Ending."
         PRINT *,"Log gas temperature is ",ltg
         PRINT *,"Log density is ",lrho
      ELSEIF(i.EQ.5) THEN
         PRINT *,"GETKAPPA: Attention. Opacity returned would have "
         PRINT *,"been zero. This will cause problems in the future."
         PRINT *,"(c.f. divide by zero). Aborting."
         PRINT *,"Abnormal program termination done."
         PRINT *,"Log gas temperature is ",ltg
         PRINT *,"Log density is ",lrho
      ELSE
         PRINT *,"GETKAPPA: Warning!"
         PRINT *,"One or more of the parameters sent to GETKAPPA"
         PRINT *,"are zero!"
      END IF
      
!      RETURN
      STOP
      END



