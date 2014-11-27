      SUBROUTINE buildsightlines
c
c--Reads in a table of sightlines on a sphere (stored in Cartesian unit
c     vectors).
c     These should have been generated using the HEALPix library 
c     functions using the program in maketables/healpixtable.f
c
      INCLUDE 'idim'

      INCLUDE 'COMMONS/sightlines'
      INCLUDE 'COMMONS/physcon'

      CHARACTER*100 homedir, filetable

      CALL GET_ENVIRONMENT_VARIABLE('SPH_HOME',homedir)
c
c--Choose how many HEALPix vectors to use to cover the sphere
c     (12, 48, 192, 768).  Typical choice would be 48.
c
c      filetable = TRIM(homedir)//'tables/healpixtbl_12'
      filetable = TRIM(homedir)//'tables/healpixtbl_48'
c      filetable = TRIM(homedir)//'tables/healpixtbl_192'

      OPEN(UNIT=8,FILE=filetable)
      READ (8,*) nsightlines
      IF (nsightlines.GT.nsightlinesmax) THEN
         WRITE (*,*) 'ERROR - nsightlines.GT.nsightlinesmax'
         CALL quit(0)
      ENDIF
      DO i = 1, nsightlines
         READ (8,*,END=100) sightline(1,i),sightline(2,i),sightline(3,i)
      END DO
      CLOSE(8)
      GOTO 200
 100  WRITE (*,*) 'ERROR - buildsightlines end of file'
      CALL quit(0)
 200  CONTINUE
c
c--Approximate values in the absence of Healpix (but these are not uniformly
c     distributed on a sphere)
c
c     nsightlines = 100
c      DO i=1,10
c         theta = i*pi/10 - pi/20
c         DO j = 1,10
c            phi = j*2.0*pi/10 - pi/10
c            k = (i-1)*10 + j
c            sightline(1,k) = sin(theta)*cos(phi)
c            sightline(2,k) = sin(theta)*sin(phi)
c            sightline(3,k) = cos(theta)
c         END DO
c      END DO

      RETURN
      END
