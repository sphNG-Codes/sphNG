      subroutine uset
      
      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/phase'
      
      INTEGER i
      REAL radius, boundtempl
      REAL getu,getcv

      DO i=1, npart
         IF (iphase(i).EQ.0) THEN
           radius = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
            boundtempl = gmw*hoverr**2*radius**
     &           (tprof+1)/((Rg/uergg)*gamma*radius)
            
            vxyzu(4,i) = getu(rho(i), boundtempl)
            ekcle(3,i) = getcv(rho(i),vxyzu(4,i))
         ENDIF
      ENDDO

      RETURN
      END

      FUNCTION getu (rho, boundtempl)

c-- Pass in rho in code units.

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/eostbl'

      REAL ltemp, lrho, getu
      REAL rhoval1, rhoval2, tval1, tval2, rug
      REAL*4 rho
      REAL tinc
      INTEGER nkrho1, nkrho2, nkt1, nkt2

      tinc = 0.005

      lrho = log10(rho*umass/udist**3)
      ltemp = log10(boundtempl)

      nkt1 = 1 + INT(ltemp/tinc)
      IF(nkt1.GE.umxt) nkt1 = umxt - 1
      nkt2 = nkt1 + 1

      nkrho1=INT(lrho/0.005)+eostbl_rho1
      IF(lrho.LT.0.0) nkrho1=nkrho1-1
      IF(nkrho1.GE.tgmxrh) nkrho1 = tgmxrh - 1
      IF(nkrho1.LT.1) nkrho1 = 1
      nkrho2=nkrho1+1

      rhoval1=(nkrho1-eostbl_rho1)*0.005
      rhoval2=(nkrho2-eostbl_rho1)*0.005
      tval1 = (nkt1 - 1)*tinc
      tval2 = (nkt2 - 1)*tinc

c      write(*,*) rho, umass, udist, lrho
c      write(*,*) nkrho1, nkrho2, nkt1, nkt2
c      write(*,468) rhoval1, rhoval2, tval1, tval2, ltemp
c 468  FORMAT (5(1PE12.5,1X))
c      STOP
c      IF(nkrho1.LT.1.0.OR.nkt1.LT.1)
c     $     CALL FAILED2(0,lu,lrho,nkrho1,nkrho2,nkt1,nkt2)
      
!     Bilinear interpolation from Numerical Recipes
      
                                !Interpolation in log10 space
      y1=(utable(nkt1,nkrho1))
      y2=(utable(nkt2,nkrho1))
      y3=(utable(nkt2,nkrho2))
      y4=(utable(nkt1,nkrho2))
      
      w=(lrho-rhoval1)/(rhoval2-rhoval1)
      v=(ltemp-tval1)/(tval2-tval1)
                                !     Final value of log10 energies
      rug=(1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4
      getu = 10.0**rug/uergg
      
      RETURN

      END FUNCTION getu
