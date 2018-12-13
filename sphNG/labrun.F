      SUBROUTINE labrun
c************************************************************
c                                                           *
c  This routine writes the title page                       *
c                                                           *
c************************************************************
#ifdef _OPENMP
      USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS
#endif

      INCLUDE 'idim'
      INCLUDE 'igrape'
      INCLUDE 'git_version.h'

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/infor'
#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif

#ifdef _OPENMP
      INTEGER numthread
#endif

c--Get time and date
c
      CALL getdat(iday, imon, iyear)
      CALL getime(ih, im, is, heure)
c
c--Write title page
c
      WRITE (iprint, 99001)
99001 FORMAT (1X,
     &        '**********************************************', /, 1X,
     &        '*                                            *', /, 1X,
     &        '*            sphNG - 3 D - RT - MHD          *', /, 1X,
     &        '*              WITH POINT MASSES             *', /, 1X,
     &        '*                                            *', /, 1X,
     &        '**********************************************')
      WRITE (iprint, 99002) gitversion
99002 FORMAT (/, 1X, 'Current version : ', A20)

      WRITE (iprint, 99010) __DATE__, __TIME__
99010 FORMAT (1X, 'Built           : ', 2A12)

      WRITE (iprint,*) '                  RUNNING WITH BINARY TREE'

#ifdef MPICOPY
      WRITE (iprint, *) 'Compiled using MPICOPY'
#endif
#ifdef MPI
      WRITE (iprint, *) 'Compiled using MPI'
#endif
#ifdef MPIALL
      WRITE (iprint, 99100) numproc
99100 FORMAT (' Running with ',I5,' MPI processes')
#endif
#ifdef _OPENMP
C$OMP PARALLEL
C$OMP&DEFAULT (none)
C$OMP&SHARED(iprint, numthread)
C$OMP MASTER
      numthread = OMP_GET_NUM_THREADS()
      WRITE (iprint, 99200)
      WRITE (iprint, 99201) numthread
C$OMP END MASTER
99200 FORMAT (' Compiled using openMP')
99201 FORMAT (' Running with ',I5,' openMP threads')
C$OMP END PARALLEL
#endif
      IF (iprint.NE.istdout) THEN
        WRITE (iprint, 99003)namerun, job, iday, imon, iyear, ih, im, is
      ENDIF
      WRITE (*, 99003)namerun, job, iday, imon, iyear, ih, im, is
99003 FORMAT (//, ' SPH run ', A20, ' is running option : ', A11, /,
     &        ' Started on : ', I2, '/', I2, '/', I4, ' at ', I2,
     &        ' h. ', I2, ' min. ', I2, ' sec. ')

      RETURN
      END
