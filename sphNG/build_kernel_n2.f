      SUBROUTINE build_kernel_n2
c
c--Calculates the integral along a line of sight through a smoothing
c     kernel of weight^2 (for density^2 or number density^2),
c     and the integral of weight^2 times y^2
c
      PARAMETER (ntable=1000)
      COMMON /columnkerneln2/ coltable_n2(ntable), coltable_n2y2(ntable)

      DATA pi/3.14159265/
c
c--Table of v^2=x^2/h^2
c
      DO i=1,ntable
         r2=(i-1)/(ntable/4.0)
         dist=SQRT(4.0-r2)
         step=dist/(4.0*ntable)
         ypos=0.

         xintegral1=0.0
         xintegral2=0.0
         DO j=1,4*ntable
            v=SQRT(r2+ypos*ypos)
            IF (v.LT.1.0) THEN
               v2=v*v
               val=1.0-1.5*v2+0.75*v2*v
            ELSE
               v2m=2.0-v
               val=0.25*v2m*v2m*v2m
            ENDIF
            xintegral1=xintegral1+val**2*step
            xintegral2=xintegral2+val**2*ypos**2*step
            ypos=ypos+step
         END DO
         coltable_n2(i)=2.0*xintegral1/pi*1260./491.
         coltable_n2y2(i)=2.0*xintegral2/pi*1260/491.
      END DO

      check=0.0
      r_old = 0.
      DO i=1,ntable
         r=SQRT((i-1)/(ntable/4.0))
         dr = r-r_old
         check=check+2*pi*r*coltable_n2(i)*dr
         r_old = r
      END DO
      WRITE(*,*) 'Integral over kernel^2 column ',check
c      STOP

      RETURN
      END
