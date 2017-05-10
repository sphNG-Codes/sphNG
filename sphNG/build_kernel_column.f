      SUBROUTINE build_kernel_column
c
c--Builds a table of the column density through the cubic SPH kernel
c
      INCLUDE 'idim'

      INCLUDE 'COMMONS/stellarradiation'

      PARAMETER (ntable=1000)
      COMMON /columnkernel/ coltable(ntable)

      DATA pi/3.14159265/
c
c--Table of v^2=x^2/h^2
c
      DO i=1,ntable
         r2=(i-1)/(ntable/4.0)
         dist=SQRT(4.0-r2)
         step=dist/(4.0*ntable)
         ypos=0.
         
         coldens=0.0
         DO j=1,4*ntable
            v=SQRT(r2+ypos*ypos)
            IF (v.LT.1.0) THEN
               v2=v*v
               val=1.0-1.5*v2+0.75*v2*v
               coldens=coldens+val*step
            ELSE
               v2m=2.0-v
               val=0.25*v2m*v2m*v2m
               coldens=coldens+val*step        
            ENDIF
            ypos=ypos+step
         END DO
         coltable(i)=2.0*coldens/pi
      END DO
c
c--Check values by integrating over table
c
      check=0.0
      r_old = 0.
      DO i=1,ntable
         r=SQRT((i-1)/(ntable/4.0))
         dr = r-r_old
         check=check+2*pi*r*coltable(i)*dr
         r_old = r
      END DO
      WRITE(*,*) 'Integral over kernel column density ',check
c      STOP

      RETURN
      END
