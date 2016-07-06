c************************************************************
c                                                           *
c  These subroutine contains all machine dependant code     *
c                                                           *
c************************************************************

c
c--Get date
c
      SUBROUTINE getdat(id, im, iy)
      IMPLICIT NONE

      INCLUDE 'COMMONS/ibmcom'

      INTEGER id,im,iy,ihour,imins,isec,imilli
      INTEGER ivalues

      CHARACTER*10 time
      CHARACTER*8 date
      CHARACTER*5 zone

      LOGICAL ifirst

      DIMENSION ivalues(8)

      DATA ifirst/.TRUE./

      CALL DATE_AND_TIME(date,time,zone,ivalues)
      id = ivalues(3)
      im = ivalues(2)
      iy = ivalues(1)
      IF (ifirst) THEN
         ifirst = .FALSE.
         ihour = ivalues(5)
         imins = ivalues(6)
         isec = ivalues(7)
         imilli = ivalues(8)
         istarttime(1) = iy
         istarttime(2) = im
         istarttime(3) = id
         istarttime(4) = ihour
         istarttime(5) = imins
         istarttime(6) = isec
         starttime = id*86400. + ihour*3600. + imins*60. + isec +
     &        imilli*0.001
      ENDIF
      RETURN
      END SUBROUTINE getdat
c
c--Get time
c
      SUBROUTINE getime(ih, imin, is, fhour)
      IMPLICIT NONE 

      INTEGER ih,imin,is,ivalues
      REAL fhour

      CHARACTER*10 time
      CHARACTER*8 date
      CHARACTER*5 zone

      DIMENSION ivalues(8)

      CALL DATE_AND_TIME(date,time,zone,ivalues)
      ih = ivalues(5)
      imin = ivalues(6)
      is = ivalues(7)
      fhour = ih + imin/60. + is/3600.
      RETURN
      ENDSUBROUTINE getime
c
c--Get time used since begining
c
      SUBROUTINE getused(tused)
      IMPLICIT NONE

      INCLUDE 'COMMONS/ibmcom'

      INTEGER i,idu,ihu,iminu,isu,imilliu,ivalues,imu
      REAL tused

      CHARACTER*10 time
      CHARACTER*8 date
      CHARACTER*5 zone

      DIMENSION ivalues(8)

      CALL DATE_AND_TIME(date,time,zone,ivalues)
      idu = ivalues(3)
      ihu = ivalues(5)
      iminu = ivalues(6)
      isu = ivalues(7)
      imilliu = ivalues(8)

      IF (ivalues(2).LT.istarttime(2)) THEN
         ivalues(2) = ivalues(2) + 12
      ENDIF
      DO i = istarttime(2), ivalues(2) - 1
         imu = MOD(i,12)
         IF (imu.EQ.4 .OR. imu.EQ.6 .OR. imu.EQ.9. .OR. imu.EQ.11) THEN
            idu = idu + 30
         ELSEIF (imu.EQ.2) THEN
            idu = idu + 28
         ELSE
            idu = idu + 31
         ENDIF
      END DO
      tused = idu*86400. + ihu*3600. + iminu*60. + isu + imilliu*0.001
     &     - starttime
      RETURN
      ENDSUBROUTINE getused
c
c--Check for file status
c
      SUBROUTINE statfile(file, ifsize)
      IMPLICIT NONE
      CHARACTER*7 file
      INTEGER ifsize
c      ifile = lstat(file, istatb)
c      ifsize = istatb(8)
      ifsize = 0
      RETURN
      ENDSUBROUTINE statfile
c
c--Get argument on command line
c
      SUBROUTINE getcom(job, inname)
      IMPLICIT NONE
      CHARACTER*11 job, inname
      CALL getarg(1, job)
      CALL getarg(2, inname)
      WRITE (*,*) 'Job: ', job, inname
      RETURN
      ENDSUBROUTINE getcom
