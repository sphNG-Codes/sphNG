      subroutine semimajoraxis(i, e, a)

      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/xforce'

      REAL e, a
      REAL mass1, mass2, radius
      REAL amx, amy, amz, angmom, angmom2
      REAL vel2, energy
      INTEGER i

      mass1 = xmass
      mass2 = xyzmh(4,i)
      
      radius = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2)

      amx = xyzmh(2,i)*vxyzu(3,i) - xyzmh(3,i)*vxyzu(2,i)
      amy = xyzmh(3,i)*vxyzu(1,i) - xyzmh(1,i)*vxyzu(3,i)
      amz = xyzmh(1,i)*vxyzu(2,i) - xyzmh(2,i)*vxyzu(1,i)
      
      angmom = amx + amy + amz
      angmom2 = amx**2 + amy**2 + amz**2
      
      vel2 = vxyzu(1,i)**2 + vxyzu(2,i)**2 + vxyzu(3,i)**2
      energy = vel2/2. - (mass1+mass2)/radius
      
      e = SQRT(1.0 + (2.*energy*angmom2/
     &     (mass1+mass2)**2))
      a = angmom2/((mass1+mass2)*(1.-e**2))

      RETURN
      end subroutine semimajoraxis
