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

      FUNCTION getu(rho, temp_in)

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
      ltemp = log10(temp_in)

      nkt1 = 1 + INT(ltemp/tinc)
      IF (nkt1.GE.umxt) THEN
         nkt1 = umxt - 1
      ELSEIF (nkt1.LT.1) THEN
         nkt1 = 1
      ENDIF
      nkt2 = nkt1 + 1

      nkrho1=INT(lrho/0.005)+eostbl_rho1
      IF (lrho.LT.0.0) nkrho1=nkrho1-1
      IF (nkrho1.GE.tgmxrh) THEN
         nkrho1 = tgmxrh - 1
      ELSEIF (nkrho1.LT.1) THEN
         nkrho1 = 1
      ENDIF
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
      u_got = 10.0**rug/uergg

c
c--Then uses Newton-Raphson iteration to get the value consistent with
c     getcv().  Without this, the values can differ by ~1%
c
      iteration = 0
 100  u_last = u_got
      t_check = u_got/getcv(rho, u_got)
      IF (u_got.LT.100.) THEN
         u_plus = u_got*1.01
      ELSE
         u_plus = u_got + 1.0
      ENDIF
      u_delta = u_plus-u_got
      t_plus = u_plus/getcv(rho, u_plus)

      func = t_check - temp_in
      derivative = (t_plus-t_check)/u_delta
      u_got = u_got - func/derivative

      IF (ABS(u_got-u_last).GT.1.E-5) THEN
         iteration = iteration + 1
         IF (iteration.GT.100) THEN
            WRITE (*,*) 'ERROR - getu failed ',u_got,temp_in
            CALL quit(1)
         ENDIF
         GOTO 100
      ENDIF

      getu = u_got

      RETURN

      END FUNCTION getu
