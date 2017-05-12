      SUBROUTINE stellar_radiation(nlst_in,nlst_end,nlst_tot,npart,ntot,
     &     itime,xyzmh)
c************************************************************
c                                                           *
c  Subroutine to compute the stellar radiation from sinks.  *
c     For MPI calculations, this only does local particles. *
c     Contributions from other MPI processes are done in    *
c     densityiterate.                                       *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/stellarradiation'

      DIMENSION xyzmh(5,mmax2)

c      print *,iproc,': entered stellar_radiation ',nlst_tot

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine stellar_radiation')
c
c--First calculate column densities between active particles and sinks
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,nlst_tot,itime,ntot)
C$OMP& shared(ivar,xyzmh,iphase)
C$OMP& private(n,ipart)
      DO n = nlst_in, nlst_end
         ipart = ivar(3,n)
         IF (iphase(ipart).EQ.0) CALL column_los(ipart,itime,ntot,xyzmh)
      END DO
C$OMP END PARALLEL DO

      IF (itrace.EQ.'all') WRITE (iprint,300)
 300  FORMAT ('exit subroutine stellar_radiation')

c      print *,'Exited stellar_radiation ',itime

      RETURN
      END
