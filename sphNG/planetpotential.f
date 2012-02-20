      subroutine planetpotential (px, py, pz, imigrate, rorbitmax,
     &     pmrate, rorbit_orig, ti)

      IMPLICIT NONE
      INCLUDE 'idim'
      INCLUDE 'COMMONS/xforce'

      INTEGER imigrate
      REAL rorbit, ang, angend
      REAL px, py, pz
      REAL rorbitmax, rorbit_orig, pmrate
      REAL ti

c--Allow for migration of potential at prescribed rate.
      IF (imigrate.EQ.1 .AND. ti.LT.rorbitmax) THEN
         rorbit = rorbit_orig + (pmrate*ti)
         ang = (2.*rorbit_orig*sqrt(xmass/rorbit_orig**3) - 
     &        2.*(pmrate*ti+rorbit_orig)*sqrt(xmass/(pmrate*
     &        ti + rorbit_orig)**3))/pmrate
      ELSEIF (imigrate.EQ.1 .AND. ti.GE.rorbitmax) THEN
         rorbit = rorbit_orig + (pmrate*rorbitmax)
c--This only needs to be calculated once really, consider moving.
         angend = (2.*rorbit_orig*sqrt(xmass/rorbit_orig**3) - 
     &        2.*(pmrate*rorbitmax+rorbit_orig)*sqrt(xmass/(pmrate*
     &        rorbitmax + rorbit_orig)**3))/pmrate
         
         ang = angend + (ti-rorbitmax)/sqrt(rorbit**3/xmass)
      ELSE
         rorbit = rorbit_orig
         ang = ti*sqrt(xmass/rorbit**3)
      ENDIF
      
      px = rorbit*COS(ang)
      py = rorbit*SIN(ang)
      pz = 0.0

      RETURN

      end subroutine planetpotential
