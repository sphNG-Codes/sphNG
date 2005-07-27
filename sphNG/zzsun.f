      SUBROUTINE zzsun

      USE f90_unix

c************************************************************
c                                                           *
c  This routine contains all machine dependant routines     *
c                                                           *
c************************************************************

      DIMENSION iarray(3), tarray(2), istatb(13)

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
      ENTRY getused(tused)
      tused = etime(tarray)
      RETURN
c
c--Check for file status
c
      ENTRY statfile(file, ifsize)
      ifile = lstat(file, istatb)
      ifsize = istatb(8)
      RETURN
c
c--Get argument on command line
c
      ENTRY getcom(job, inname)
      CALL getarg(1, job)
      CALL getarg(2, inname)
      WRITE (*,*) 'Job: ', job, inname
      RETURN

      END
