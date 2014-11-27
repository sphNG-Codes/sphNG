      PROGRAM healpixtable
c
c--Writes out an ASCII table of points in 3-D for a Healpix sphere
c     Needs to be compiled with the Healpix library
c
c     Compile command (using ifort) is:
c     ifort healpixtable.f -openmp -I/sw/Healpix/Healpix_3.11/include -L/sw/Healpix/Healpix_3.11/lib -L/sw/fits/lib -lhealpix -lhpxgif -lcfitsio
c
      use pix_tools

      INTEGER nside, npoints
      REAL*8 vector(3)

      WRITE (*,*) 'Enter nside'
      WRITE (*,*) '    1: 12 points'
      WRITE (*,*) '    2: 48 points'
      WRITE (*,*) '    4: 192 points'
      WRITE (*,*) '    8: 768 points'
      READ (*,*) nside

      IF (nside.EQ.1) THEN
         npoints = 12
      ELSEIF (nside.EQ.2) THEN
         npoints = 48
      ELSEIF (nside.EQ.4) THEN
         npoints = 192
      ELSEIF (nside.EQ.8) THEN
         npoints = 768
      ENDIF

      OPEN (11,FILE='healpixtbl')
      WRITE (11,*) npoints
      DO i = 0, npoints-1
         CALL pix2vec_ring(nside,i,vector)
         WRITE (11,99001) vector
      END DO
99001 FORMAT(3(1PE12.5,1X))
      CLOSE(11)

      CALL pix2vec_ring(nside,npoints,vector)

      END
