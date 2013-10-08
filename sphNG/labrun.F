      SUBROUTINE labrun
c************************************************************
c                                                           *
c  This routine writes the title page                       *
c                                                           *
c************************************************************

      INCLUDE 'igrape'
      INCLUDE 'svn_version.h'

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/infor'

c--Get time and date
c
      CALL getdat(iday, imon, iyear)
      CALL getime(ih, im, is, heure)
c
c--Write title page
c
      WRITE (iprint, 99001)
99001 FORMAT ('1', /, 1X,
     &        '**********************************************', /, 1X,
     &        '*                                            *', /, 1X,
     &        '*      S P H - 3 D - TIMESTEPS - GRAPE       *', /, 1X,
     &        '*              WITH POINT MASSES             *', /, 1X,
     &        '*                                            *', /, 1X,
     &        '**********************************************')
      WRITE (iprint, 99002) svnversion
99002 FORMAT (/, 1X, 'Current version : ', A20)

      IF (igrape.EQ.0) THEN
         WRITE (iprint,*) '                  RUNNING WITH BINARY TREE'
      ELSEIF (igrape.EQ.1) THEN
         WRITE (iprint,*) '                  RUNNING WITH GRAPE'
      ENDIF

#ifdef MPICOPY
      WRITE (iprint, *) 'Compiled using MPICOPY'
#endif
#ifdef MPI
      WRITE (iprint, *) 'Compiled using MPI'
#endif

      WRITE (iprint, 99003)namerun, job, iday, imon, iyear, ih, im, is
      WRITE (*, 99003)namerun, job, iday, imon, iyear, ih, im, is
99003 FORMAT (//, ' SPH run ', A20, ' is running option : ', A11, /,
     &        ' Started on : ', I2, '/', I2, '/', I4, ' at ', I2,
     &        ' h. ', I2, ' min. ', I2, ' sec. ')

      RETURN
      END
