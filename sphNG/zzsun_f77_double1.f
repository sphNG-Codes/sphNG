      SUBROUTINE zzsun

c      USE f90_unix

c************************************************************
c                                                           *
c  This routine contains all machine dependant routines     *
c                                                           *
c************************************************************

      INTEGER*4 iarray(3), istatb(13)
      INTEGER*4 ival1, ival2

c      REAL*4 tarray(2)

      CHARACTER*11 job, inname
      CHARACTER*7 file
c
c--Get date
c
      ENTRY getdat(id, im, iy)
      CALL idate(iarray)
      id = iarray(1)
      im = iarray(2)
      iy = iarray(3)
      RETURN
c
c--Get time
c
      ENTRY getime(ih, im, is, fhour)
      CALL itime(iarray)
      ih = iarray(1)
      im = iarray(2)
      is = iarray(3)
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
      ifsize = istatb(8)
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
