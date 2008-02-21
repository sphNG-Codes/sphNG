      SUBROUTINE zeusread

      INCLUDE 'COMMONS/zeus'
c
c--Read in from ZEUS/IDL file
c
      OPEN (44,FILE='outputgrid.dat')
      READ (44,*) nradius, ntheta
      IF (nradius.GT.nradiusmax .OR. ntheta.GT.nthetamax) THEN
         WRITE (*,*) 'Grid too big in outputgrid.dat'
         STOP
      ENDIF
      DO l = 1, 6
         READ (44,*) (radii(i), i=1, nradius)
         READ (44,*) (thetas(i), i=1, ntheta)
         READ (44,*) ((density(l,i,j,1), j=1, ntheta), i=1, nradius)
         READ (44,*) ((vr(l,i,j,1), j=1, ntheta), i=1, nradius)
         READ (44,*) ((vt(l,i,j,1), j=1, ntheta), i=1, nradius)
         READ (44,*) ((vp(l,i,j,1), j=1, ntheta), i=1, nradius)
         READ (44,*) ((density(l,i,j,2), j=1, ntheta), i=1, nradius)
         READ (44,*) ((vr(l,i,j,2), j=1, ntheta), i=1, nradius)
         READ (44,*) ((vt(l,i,j,2), j=1, ntheta), i=1, nradius)
         READ (44,*) ((vp(l,i,j,2), j=1, ntheta), i=1, nradius)
      ENDDO

      zmasses(1) = 0.000003
      zmasses(2) = 0.00001
      zmasses(3) = 0.00003
      zmasses(4) = 0.0001
      zmasses(5) = 0.0003
      zmasses(6) = 0.001

      CLOSE (44)

      RETURN

      END
