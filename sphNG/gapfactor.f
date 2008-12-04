      subroutine gapfactor(variation, hmass, gapfac)
      
      INCLUDE 'COMMONS/gaptbl'

c      print *, '------------'
c      print *, gaplen, gapwid
c         DO i=1, gaplen
c            write(*,*) (gaptable(i,j), j=1, gapwid)
c         ENDDO

      DO i = gapwid, 1, -1
         IF (variation.GE.gaptable(1,i)) THEN
            ldr = i
            mdr = i+1
            exit
         ENDIF
      END DO
      
      DO i = gaplen, 1, -1
         IF (hmass.GE.gaptable(i,1)) THEN
            lpmass = i
            mpmass = i+1
            exit
         ENDIF
      END DO

      IF (lpmass.EQ.gaplen) THEN
         lpmass = gaplen -1
         mpmass = gaplen
      ENDIF

      IF (ldr.EQ.gapwid) THEN
         ldr = gapwid -1
         mdr = gapwid
      ENDIF
      
      y1=log10(gaptable(lpmass,ldr))
      y2=log10(gaptable(mpmass,ldr))
      y3=log10(gaptable(mpmass,mdr))
      y4=log10(gaptable(lpmass,mdr))
      
      w=(variation-gaptable(1,ldr))/
     $     (gaptable(1,mdr)-gaptable(1,ldr))
      v=(hmass-gaptable(lpmass,1))/
     $     (gaptable(mpmass,1)-gaptable(lpmass,1))
      gapfac = (1.0-v)*(1.0-w)*y1+v*(1.0-w)*y2+w*v*y3+(1.0-v)*w*y4
      
      gapfac = 10**gapfac
      
c      print *, ldr,mdr, lpmass,mpmass
c      print *, hmass, variation
c      print *, gaptable(1,ldr),gaptable(lpmass,1)
c      print *, gaptable(lpmass,ldr),gaptable(mpmass,mdr)
c      print *, 'Gap Factor =', gapfac

c      STOP
      RETURN
      END
