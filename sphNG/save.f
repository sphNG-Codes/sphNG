      SUBROUTINE save
c************************************************************
c                                                           *
c  This routine determines whether or not the time has come *
c     to write a dump on disk and print out a detailed      *
c     state of the run                                      *
c                                                           *
c************************************************************

      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/stop'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/binary'

      CHARACTER*7 where

      DATA where/'save'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine save')
c
c--Increment time step counter
c
      ncount = ncount + 1
      nbuild = nbuild + 1
c
c--Get remaining time for the job
c
      CALL getused(tused)
      tleft = 60.*tmax - tused
c
c--Evaluate time needed for next timestep
c
      nleft = nstep - ncount
      tneed = 1.0*tstep + 1.0
c
c--If enough time and no dump required go on with integration
c
      IF (tleft.LT.tneed) istop = 1
      IF (gt.GT.tstop) istop = 1
      IF (ibound.EQ.8 .AND. naccrete.GT.nstop) istop = 1
      IF (nleft.LE.0) THEN
         ncount = 0
         IF (ibound.EQ.8 .AND. naccrete.GT.nstop-nfastd) nstep = 1
      ENDIF
c
c--Transform into original frame of reference
c
      CALL chanref(1)
c
c--Compute and print out the state of the system
c
      CALL inform(where)
c
c--Transform into expanding frame of reference
c
      CALL chanref(2)

      IF (idebug(1:4).EQ.'save') THEN
         WRITE (iprint, 99002) nstep, ncount, nleft
99002    FORMAT (1X, 3(I4,1X))
         WRITE (iprint, 99003) tmax, tleft, tstep, tneed
99003    FORMAT (1X, 5(1PE12.5,1X))
      ENDIF
c
c--Check for run termination
c
      IF (istop.EQ.1) CALL endrun

      RETURN
      END
