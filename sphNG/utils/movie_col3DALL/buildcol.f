      SUBROUTINE buildcol

      COMMON /columntab/ coltable(1000)

      DATA pi/3.14159265/
c
c--Table of v=x/h
c
c      DO i=1,1000
c         r=(i-1)/500.
c         dist=SQRT(4.0-r*r)
c         step=dist/4000.0
c         ypos=0.
c         
c         coldens=0.0
c         DO j=1,4000
c            v=SQRT(r*r+ypos*ypos)
c            IF (v.LT.1.0) THEN
c               v2=v*v
c               val=1.0-1.5*v2+0.75*v2*v
c               coldens=coldens+val*step
c            ELSE
c               v2m=2.0-v
c               val=0.25*v2m*v2m*v2m
c               coldens=coldens+val*step        
c            ENDIF
c            ypos=ypos+step
c         END DO
c         coltable(i)=2.0*coldens/pi
c      END DO
c
c--Table of v^2=x^2/h^2
c
      DO i=1,1000
         r2=(i-1)/250.
         dist=SQRT(4.0-r2)
         step=dist/4000.0
         ypos=0.
         
         coldens=0.0
         DO j=1,4000
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

c      check=0.0
c      DO i=1,1000
c         r=(i-1)/500.
c         WRITE(*,*) r, table(i)
c         check=check+2*pi*r*table(i)*(1.0/500.)
c      END DO
c      WRITE(*,*)check

      RETURN
      END
