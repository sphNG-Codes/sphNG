      SUBROUTINE value (xi,yi,zi,hi,gx,gy,val,pmassi)
c***********************************************************
c                                                          *
c  this subroutine computes the value of the quantities at *
c  one given point.                                        *
c                                                          *
c***********************************************************
c
      INCLUDE 'idim'

      COMMON /columntab/ coltable(1000)
c
      pi=3.14159265
c
      val = 0.0
c
      vr2 = ((xi-gx)**2 + (yi-gy)**2)/hi**2

      IF (vr2.GE.4.0) RETURN

      vr = SQRT(vr2)
c
c--find nearest table points
c
      pos=INT(vr*500.0)+1
      IF (pos.EQ.1000) THEN
         pos1=999
         pos2=1000
      ELSE
         pos1=pos
         pos2=pos1+1
      ENDIF
      r1=(pos1-1)/500.
      r2=(pos2-1)/500.
      valuemean = (vr-r1)/(r2-r1)*(coltable(pos2)-coltable(pos1))
      valuemean = valuemean + coltable(pos1)
c
      val = valuemean*pmassi/hi/hi
c
      RETURN
      END

