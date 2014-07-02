      FUNCTION GETKAPPA(u2,cv2,rho2)
  
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/optbl'
      INCLUDE 'COMMONS/dusttbl'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/rbnd'

      COMMON /getkap/ iflag

      REAL*4 rho2
      REAL lTg,lrho,getkappa,u2,u,cv,cv2,rho
      REAL rkappa,y1,y2,y3,y4,w,v
      REAL minop, mu, mintemp, maxtemp, maxop
      REAL grad,val1,val1h,val1l,val2,val2h,val2l,valdown,valup
      REAL xlogR,rkappa_ferguson,xlogR_bound,ltg_bound,dustfac,uval
      INTEGER ndrho,nkrho1,nkrho2,nktg1,nktg2,iflag,nkR1
      LOGICAL tryalex

      REAL trilinear
c
c--Constant opacity for shock testing
c
c      getkappa=4.0d2*umass/udist**2
c      RETURN
c
c--Dust opacity can be reduced through "opdenom>1" or from "metallicity<1"
c
      dustfac = opdenom/metallicity
c
c--Convert to dimensionfull quantities
c
      u=u2*uergg
      rho=rho2*umass/udist**3
      cv=cv2*uergg

      IF(u.EQ.0.0.OR.rho.EQ.0.0.OR.cv.EQ.0.0) CALL FAILED(7,ltg,lrho)
c
c--Log10 of gas temperature and density
c
      ltg = log10(u/cv)
      lrho = log10(rho)
c
c--Original opacities (Whitehouse & Bate 2006 to Bate 2012)
c     These use Pollack et al. (1985) for low temperatures and
c     Alexander (1975) for high temperatures.
c
      nkrho1=2
      nkrho2=3
      nktg1=2
      nktg2=3

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
c!     CALL quit(1)
c         ENDIF

         grad = (val2 - val1)/0.05
         valup = (ltg - 4.0) * grad + val2
!Kramer's opacity + electron scattering (0.4 cm^2/g, according to 
!page two of astro-ph/0410343
         valdown = (1.2512E+22* rho *(10.0**ltg)**(-3.5)) + 0.4
!     PRINT *,valup,valdown, (10.0**ltg)**(-3.5)
         rkappa = MIN(valup,valdown)*umass/udist**2
c      write (*,*) 'Exit kap'
         GOTO 666
c         RETURN
      ENDIF
      
      nktg2=nktg1+1

      IF(nkrho1.LE.1.OR.nkrho1.GE.opmxrh.OR.nktg1.LE.1.OR.
     $     nktg1.GE.opmxtg)
     $     CALL FAILED(0,ltg,lrho)

!     Bilinear interpolation from Numerical Recipes
!     Pollack table gives kappa, Alexander gives log10 kappa
!     Interpolation is in log kappa
      
      IF (iopmodel/10.EQ.0) THEN
         IF(optable(nktg1,nkrho1).EQ.0.00) CALL FAILED(1,ltg,lrho)
         IF(optable(nktg2,nkrho1).EQ.0.00) CALL FAILED(2,ltg,lrho)
         IF(optable(nktg2,nkrho2).EQ.0.00) CALL FAILED(3,ltg,lrho)
         IF(optable(nktg1,nkrho2).EQ.0.00) CALL FAILED(4,ltg,lrho)
      ENDIF
!TEMP CODE
!      IF(optable(nktg1,nkrho1).EQ.0.00.OR.optable(nktg2,nkrho1).EQ.0.00.
!     $     OR.optable(nktg2,nkrho2).EQ.0.00.OR.
!     $     optable(nktg1,nkrho2).EQ.0.00) THEN
!         getkappa=0.0
!         RETURN
!      ENDIF

      !Interpolation in log10 space
      y1=log10(MAX(optable(nktg1,nkrho1),1.0E-20))
      y2=log10(MAX(optable(nktg2,nkrho1),1.0E-20))
      y3=log10(MAX(optable(nktg2,nkrho2),1.0E-20))
      y4=log10(MAX(optable(nktg1,nkrho2),1.0E-20))

      w=(lrho-optable(1,nkrho1))/
     $     (optable(1,nkrho2)-optable(1,nkrho1))
      v=(ltg-optable(nktg1,1))/
     $     (optable(nktg2,1)-optable(nktg1,1))

!     Final value of log10 opacity
      rkappa=(1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4
c
c--For mixed tables, make sure interpolation does not go past end of Pollack
c     table
c
      IF (iopmodel/10.EQ.1) THEN
         IF (y1.LT.-19.0 .OR. y2.LT.-19.0 .OR. y3.LT.-19.0 .OR. 
     &        y4.LT.-19.0) rkappa = -20.0
      ENDIF
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

      IF (dustfac.NE.1.0) THEN
         IF (iopmodel/10.EQ.0) THEN
c--New index, ndrho, required for new dusttbl now with different
c  dimensions to optable.
            ndrho = INT(((lrho + 14)/0.014)+2)
            IF (ndrho.LE.1) ndrho = 2
            IF (ndrho.GT.duslen) ndrho = 1001

            mu = (lrho-dusttable(ndrho,1))/
     &           (dusttable(ndrho-1,1) - dusttable(ndrho,1))

            minop = log10(dusttable(ndrho,3)) + mu*(log10(
     &           dusttable(ndrho-1,3)) -
     &           log10(dusttable(ndrho,3)))

            mintemp = dusttable(ndrho,2) + mu*(dusttable(ndrho-1,2) -
     &           dusttable(ndrho,2))

            maxop = log10(dusttable(ndrho,5)) + mu*(log10(
     &           dusttable(ndrho-1,5)) -
     &           log10(dusttable(ndrho,5)))

            maxtemp = dusttable(ndrho,4) + mu*(dusttable(ndrho-1,4) -
     &           dusttable(ndrho,4))

            minop = 10**(minop)
            minop = minop*umass/udist**2
            maxop = 10**(maxop)
            maxop = maxop*umass/udist**2

            IF (ltg.le.mintemp .AND. ltg.ge.maxtemp) THEN
               rkappa = MAX(rkappa/dustfac, minop)
            ELSEIF (ltg.le.mintemp) THEN
               rkappa = rkappa/dustfac
            ELSE
               rkappa = rkappa
            ENDIF
         ELSE
            rkappa = rkappa/dustfac
         ENDIF
      ENDIF

 666  CONTINUE
c
c--Below here are the additions for using the Ferguson et al. (2005) opacities
c     at higher temperatures (beyond dust sublimation).
c
      IF (iopmodel/10.EQ.1) THEN
c
c--Use Ferguson et al. (2005) for high temperatures
c
         IF (ltg.GT.3.0) THEN
c
c--Value from Ferguson
c
            xlogR = lrho - 3.0*(ltg - 6.0)

            xlogR_bound = MAX(Rmin_ferguson,MIN(Rmax_ferguson,xlogR))
            ltg_bound = MAX(3.0,MIN(Tmax_ferguson,ltg))
c
c--Ensure that log10(R) is within Ferguson table bounds
c
c     Otherwise, use normal opacities
c
c--Ferguson table indices
c
            IF (ltg.GT.3.5) THEN
               nktg1 = (ltg_bound - 3.5)/0.05 + 66
            ELSEIF (ltg.LT.3.5) THEN
               nktg1 = (ltg_bound - 3.0)/0.01 + 16
            ENDIF
c
c--Don't use temperatures below 1000 K from Ferguson
c
            nktg1 = MAX(16,MIN(opmxtg_ferguson-1,nktg1))

            nkR1 = (xlogR_bound - Rmin_ferguson)/Rdelta_ferguson + 2
            nkR1 = MAX(2,MIN(opmxR_ferguson-1,nkR1))
c
c--Use trilinear interpolation in log space because table values are logs
c
            uval=(xlogZ-metaltable(metaltable1))/
     &           (metaltable(metaltable1+1)-metaltable(metaltable1))

            rkappa_ferguson = trilinear(0,ltg_bound,xlogR_bound,uval,
     &           nktg1,nkR1,metaltable1,
     &           optable_ferguson,opmxtg_ferguson,opmxR_ferguson,
     &           opmxMet_ferguson)

            rkappa_ferguson = 10**rkappa_ferguson*umass/udist**2
c
c--Now choose between Pollack et al. values at low temperatures and Ferguson
c     at high temperatures for intermediate values
c
            IF (rkappa.LE.1.0E-9 .OR. ltg.GT.3.3979 .OR. 
     &           lrho.LT.-14 .AND. ltg.GT.3.35) 
     &           rkappa = rkappa_ferguson
         ENDIF
      ENDIF

      getkappa = rkappa
!      getkappa=rkappa!/1000.0
!      getkappa=MAX(rkappa,1.0/(10.0*rmax*rho2))
!      getkappa=MAX(rkappa,1.0/(rmax*4.0*0.01))
!      getkappa=15.0*umass/udist**2
      IF(getkappa.EQ.0.00) THEN
         print *,'Info: ',u/cv,lrho,val1l,val1h,val2l,val2h
         print *,w,val1,val2,grad,valup,ltg,valdown
         CALL FAILED(5,ltg,lrho)
      ENDIF

      RETURN

      END !SUBROUTINE GETKAPPA
      

      FUNCTION trilinear(itype, xval, yval, zval, ix, iy, iz, table, 
     &     idim1, idim2, idim3)
c
c--Performs trilinear interpolation on a 3-D table of values
c
      INTEGER itype, idim1, idim2, idim3
      REAL trilinear,xval,yval,zval,table(idim1, idim2, idim3)
      REAL value,y1,y2,y3,y4,y5,y6,y7,y8,bilinear1,bilinear2
c
c--Performs trilinear interpolation in the table
c     itype=0 is linear, itype=1 is in log space
c
      y1 = table(ix,iy,iz)
      y2 = table(ix+1,iy,iz)
      y3 = table(ix+1,iy+1,iz)
      y4 = table(ix,iy+1,iz)

      y5 = table(ix,iy,iz+1)
      y6 = table(ix+1,iy,iz+1)
      y7 = table(ix+1,iy+1,iz+1)
      y8 = table(ix,iy+1,iz+1)
      IF (itype.EQ.1) THEN
         y1=log10(y1)
         y2=log10(y2)
         y3=log10(y3)
         y4=log10(y4)
         y5=log10(y5)
         y6=log10(y6)
         y7=log10(y7)
         y8=log10(y8)
      ELSEIF (itype.NE.0) THEN
         WRITE (*,*) 'ERROR - trilinear interpolation ',itype
         CALL quit(1)
      ENDIF
      
c      print *,xval,table(ix,1), yval,table(1,iy)

      w=(yval-table(1,iy,iz))/
     &     (table(1,iy+1,iz)-table(1,iy,iz))
      v=(xval-table(ix,1,iz))/
     &     (table(ix+1,1,iz)-table(ix,1,iz))

      bilinear1 = (1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4
      bilinear2 = (1.0-v)*(1.0-w)*y5+v*(1.0-w)*y6+w*v*y7+(1.0-v)*w*y8

      value = (1.0-zval)*bilinear1 + zval*bilinear2

      IF (itype.EQ.1) value = 10**value

c      print *,'trilinear ',value,y1,y2,y3,y4,y5,y6,y7,y8

      trilinear = value

      RETURN
      END


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
      
      CALL quit(1)
      END



