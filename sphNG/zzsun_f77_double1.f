      SUBROUTINE zzsun

c      USE f90_unix

c************************************************************
c                                                           *
c  This routine contains all machine dependant routines     *
c                                                           *
c************************************************************

      INTEGER*4 iarray(3), istatb(13)
      INTEGER*4 ival1, ival2
      INTEGER*4 Dtime(8)

c      REAL*4 tarray(2)

      CHARACTER*11 job, inname
      CHARACTER*7 file
c
c--Get date
c
      ENTRY getdat(id, im, iy)
c      CALL idate(iarray)
c      id = iarray(1)
c      im = iarray(2)
c      iy = iarray(3)
      CALL DATE_AND_TIME(VALUES=Dtime)
      id = Dtime(3)
      im = Dtime(2)
      iy = Dtime(1)
      RETURN
c
c--Get time
c
      ENTRY getime(ih, im, is, fhour)
c      CALL itime(iarray)
c      ih = iarray(1)
c      im = iarray(2)
c      is = iarray(3)
      CALL DATE_AND_TIME(VALUES=Dtime)
      ih = Dtime(5)
      im = Dtime(6)
      is = Dtime(7)
      fhour = ih + im/60. + is/3600.
      RETURN
c
c--Get time used since begining
c
c      ENTRY getused(tused)
c      tused = etime(tarray)
c      RETURN
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
      ival1 = 1
      ival2 = 2
      CALL getarg(ival1, job)
      CALL getarg(ival2, inname)
      WRITE (*,*) 'Job: ', job, inname
      RETURN

      END
