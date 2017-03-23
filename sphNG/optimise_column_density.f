      SUBROUTINE optimise_column_density(m,itime,ntot)
c
c********************************************************************* 
c                                                                    *
c     For use with column_density.F subroutine.                      *
c     Optimisation so that if column density doesn't change much     *
c     from step to step, don't recalculate it every particle time.   *
c                                                                    *
c********************************************************************* 
c
      IMPLICIT NONE

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INTEGER m,itime,ntot

      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/sort'

      REAL fracchange

      IF (idustRT.NE.1 .OR. ioptimise_column.NE.1) THEN
         WRITE (*,*) 'ERROR - idustRT.NE.1 .OR. ioptimise_column.NE.1'
         CALL quit(0)
      ENDIF

      IF (m.GT.ntot) THEN
         WRITE (*,*) 'ERROR - m.GT.ntot'
         CALL quit(0)
      ENDIF

      IF (itime.EQ.0) THEN
         icolumnsteps(m) = isteps(m)
      ELSE
         IF (heatingISRold(1,m).LE.0. .OR. heatingISRold(2,m).LT.0.)
     &        THEN
            WRITE (*,*) 'ERROR - heatingISRold(1 or 2,m)<=0 ',itime,m,
     &           heatingISRold(1,m),heatingISRold(2,m),iorig(m),
     &           iunique(iorig(m))
            CALL quit(1)
         ENDIF
         fracchange = MAX(
     &        ABS(heatingISR(1,m)/heatingISRold(1,m) - 1.0),
     &        ABS(heatingISR(4,m)/heatingISRold(2,m) - 1.0))
         IF (fracchange.LT.0.02) THEN
c
c--Make sure that if increasing the time between recalculations, that an
c     whole number of recalculations can be done between now and the
c     synchronisation timestep.
c
            IF (MOD(imaxstep-itime,icolumnsteps(m)*2).EQ.0) THEN
               icolumnsteps(m) = icolumnsteps(m)*2
            ENDIF
         ELSEIF (fracchange.GT.0.05) THEN
            icolumnsteps(m) = icolumnsteps(m)/2
         ENDIF
      ENDIF
      icolumnsteps(m) = MIN(imaxstep,MAX(icolumnsteps(m),isteps(m)))
      icolumnnext(m) = itime + icolumnsteps(m)
      IF (icolumnnext(m).GT.imaxstep) icolumnnext(m) = 
     &     MOD(icolumnnext(m),imaxstep)
      IF (icolumnnext(m).EQ.0) icolumnnext(m) = imaxstep

c      IF (m.LT.5) THEN
c         print *,itime,m,icolumnnext(m),icolumnsteps(m),isteps(m),
c     &        fracchange, heatingISR(1,m), heatingISRold(1,m),
c     &        heatingISR(4,m), heatingISRold(2,m)
c      ENDIF

      heatingISRold(1,m) = heatingISR(1,m)
      heatingISRold(2,m) = heatingISR(4,m)

      RETURN
      END
c
c------------------------------------------------------------------------
c
