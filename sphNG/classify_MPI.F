      FUNCTION classify(mpitype,xyzmh)

      IMPLICIT NONE

      INCLUDE 'idim'

      CHARACTER*1 mpitype
      REAL xyzmh(5,1)

      INTEGER classify_binary,classify_cartesian,classify_disc,
     &     classify_sphere

      INTEGER classify,ibelong

      IF (mpitype.EQ.'b') THEN
         ibelong = classify_binary(xyzmh)
      ELSEIF (mpitype.EQ.'c') THEN
         ibelong = classify_cartesian(xyzmh)
      ELSEIF (mpitype.EQ.'d') THEN
         ibelong = classify_disc(xyzmh)
      ELSEIF (mpitype.EQ.'p') THEN
         ibelong = classify_disc(xyzmh)
      ELSEIF (mpitype.EQ.'r') THEN
         ibelong = classify_disc(xyzmh)
      ELSEIF (mpitype.EQ.'s') THEN
         ibelong = classify_sphere(xyzmh)
      ELSE
         WRITE (*,*) 'ERROR - mpitype ',mpitype,' not defined'
         CALL quit(1)
      ENDIF

      classify = ibelong

      RETURN
      END


c---------------------------------------------------------------------------

      FUNCTION classify_disc(xyzmh)
c************************************************************
c                                                           *
c  Subroutine to determine which MPI process a particle     *
c     should be on.                                         *
c     This version is applicable to disc geometries.        *
c     The split is done on z (assumed rotation axis) first. *
c     The split is located at xlocationmpi, ylocationmpi,   *
c     zlocationmpi.
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/rbnd'

      DIMENSION xyzmh(5,1)

      INTEGER classify_disc

      ival = 0

#ifdef MPI
      xcen = xyzmh(1,1) - xlocationmpi
      ycen = xyzmh(2,1) - ylocationmpi
      zcen = xyzmh(3,1) - zlocationmpi

      IF (mpitype.EQ.'p') THEN
         IF (xcen**2+ycen**2+zcen**2.GT.(0.5*phibound)**2) GOTO 100
      ENDIF

      IF (zcen.LT.0.0) THEN
         ival = 1
      ELSE
         ival = 0
      ENDIF

      phiincrement = 2.0*pi/(numproc/2)

      IF (xcen.EQ.0. .AND. ycen.EQ.0.) THEN
         ival = 0
      ELSE
         ival = ival*numproc/2 +
     &     INT((pi + ATAN2(ycen,xcen))/phiincrement)
      ENDIF

 100  IF (ival.GE.numproc) ival = numproc - 1
#endif

      classify_disc = ival

      RETURN
      END


c---------------------------------------------------------------------------

      FUNCTION classify_sphere(xyzmh)
c************************************************************
c                                                           *
c  Subroutine to determine which MPI process a particle     *
c     should be on.                                         *
c     This version is applicable to disc geometries.        *
c     The split is done on z (assumed rotation axis) first. *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'

      DIMENSION xyzmh(5,1)

      INTEGER classify_sphere

      ival = 0

      classify_sphere = ival

      RETURN
      END


c---------------------------------------------------------------------------

      FUNCTION classify_cartesian(xyzmh)
c************************************************************
c                                                           *
c  Subroutine to determine which MPI process a particle     *
c     should be on.                                         *
c     This version is applicable to cartesian geometries.   *
c     The split is done on x then y then z, then back to x. *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'

      DIMENSION xyzmh(5,1)

      INTEGER classify_cartesian

      ival = 0

#ifdef MPI
      DO idimension = 1, 3
         DO number = 1, nsplit(idimension) - 1
            IF (xyzmh(idimension,1).LT.
     &           splitcoords(number,idimension)) THEN
               ival = ival + number - 1
               GOTO 100
            ENDIF
         END DO
         ival = ival + nsplit(idimension) - 1
 100     IF (idimension.LT.3) ival = ival*nsplit(idimension+1)
      END DO

      IF (ival.GE.numproc) THEN
         WRITE (*,*) 'ERROR - classify_cartesian ',numproc,ival
         CALL quit(1)
      ENDIF
#endif

      classify_cartesian = ival

      RETURN
      END

c---------------------------------------------------------------------------

      FUNCTION classify_binary(xyzmh)
c************************************************************
c                                                           *
c  Subroutine to determine which MPI process a particle     *
c     should be on.                                         *
c     This version is applicable to cartesian geometries.   *
c     It does binary splits.                                *
c     The splits are done on x then y, z, then back to x.   *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'

      DIMENSION xyzmh(5,1)

      INTEGER classify_binary

      ival = 0

#ifdef MPI
      nsplits = INT(LOG(REAL(numproc)+0.1)/LOG(2.0))
      DO i = 1, nsplits
         idimension = MOD(i-1,3)+1

         IF (i.EQ.1) THEN
            ipos = 1
         ELSE
            ipos = MOD(ival,2**(i-1)) + 2**(i-1)
         ENDIF

         IF (xyzmh(idimension,1).GE.binarysplits(ipos))
     &        ival = ival + 2**(i-1)
      END DO

      IF (ival.GE.numproc) THEN
         WRITE (*,*) 'ERROR - classify_binary ',numproc,ival
         CALL quit(1)
      ENDIF
#endif

      classify_binary = ival

      RETURN
      END
