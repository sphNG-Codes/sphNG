      SUBROUTINE zzibm
c************************************************************
c                                                           *
c  This routine contains all machine dependant routines     *
c                                                           *
c************************************************************

c      use f90_unix
      
      INCLUDE 'COMMONS/ibmcom'

      DIMENSION ivalues(8), istatb(19)

      CHARACTER*11 job, inname
      CHARACTER*10 time
      CHARACTER*8 date
      CHARACTER*5 zone
      CHARACTER*7 file
      LOGICAL ifirst

      DATA ifirst/.TRUE./
c
c--Get date
c
      ENTRY getdat(id, im, iy)
      CALL DATE_AND_TIME(date,time,zone,ivalues)
      id = ivalues(3)
      im = ivalues(2)
      iy = ivalues(1)
      IF (ifirst) THEN
         ifirst = .FALSE.
         ihour = ivalues(5)
         imins = ivalues(6)
         isec = ivalues(7)
         istarttime(1) = iy
         istarttime(2) = im
         istarttime(3) = id
         istarttime(4) = ihour
         istarttime(5) = imins
         istarttime(6) = isec
         starttime = id*86400. + ihour*3600. + imins*60. + isec
      ENDIF
      RETURN
c
c--Get time
c
      ENTRY getime(ih, imin, is, fhour)
      CALL DATE_AND_TIME(date,time,zone,ivalues)
      ih = ivalues(5)
      imin = ivalues(6)
      is = ivalues(7)
      fhour = ih + imin/60. + is/3600.
      RETURN
c
c--Get time used since begining
c
      ENTRY getused(tused)
      CALL DATE_AND_TIME(date,time,zone,ivalues)
      idu = ivalues(3)
      ihu = ivalues(5)
      iminu = ivalues(6)
      isu = ivalues(7)

      IF (ivalues(2).LT.istarttime(2)) THEN
         ivalues(2) = ivalues(2) + 12
      ENDIF
      DO i = istarttime(2), ivalues(2) - 1
         imu = MOD(i,12)
         IF (imu.EQ.4 .OR. imu.EQ.6 .OR. imu.EQ.9. .OR. imu.EQ.11) THEN
            idu = idu + 30
         ELSEIF (im.EQ.2) THEN
            idu = idu + 28
         ELSE
            idu = idu + 31
         ENDIF
      END DO
      tused = idu*86400. + ihu*3600. + iminu*60. + isu - starttime
      RETURN
c
c--Check for file status
c
      ENTRY statfile(file, ifsize)
c      ifile = lstat(file, istatb)
c      ifsize = istatb(8)
      ifsize = 0
      RETURN
c
c--Get argument on command line
c
      ENTRY getcom(job, inname)
      CALL getarg(1, job)
      CALL getarg(2, inname)
      WRITE (*,*) 'Job: ', job, inname
      RETURN
c
c--Get argument on command line
c
      ENTRY flush(iunit)
      CALL FLUSH_(iunit)
      RETURN

      END
