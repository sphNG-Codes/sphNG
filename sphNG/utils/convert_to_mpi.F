      PROGRAM convert_to_mpi
c
c--Converts full sphNG scalar dump to an MPI dump file (i.e. multiple blocks)
c
c--NOTE: To make a small memory footprint you need to:
c     1) Set idim_MPI=1 in idim
c     2) Comment out COMMONS/treecom_P in wdump.F and the associated 
c           EQUIVALENCE statement
c
      INCLUDE '../idim'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      PARAMETER (maxnumber = 10000)

      INCLUDE '../COMMONS/tming'
      INCLUDE '../COMMONS/units'
      INCLUDE '../COMMONS/sort'
      INCLUDE '../COMMONS/stepopt'
      INCLUDE '../COMMONS/phase'
      INCLUDE '../COMMONS/part'
      INCLUDE '../COMMONS/ptmass'
      INCLUDE '../COMMONS/varmhd'
      INCLUDE '../COMMONS/files'
      INCLUDE '../COMMONS/physcon'
      INCLUDE '../COMMONS/logun'
      INCLUDE '../COMMONS/actio'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      CHARACTER*1 ianswer
      CHARACTER*7 fullname(maxnumber)

#ifdef MPIALL
      iloop = 0
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)
#endif

      inname = 'ipb07'
      CALL lunit

      CALL options

      nfullstep = 1

      mpitype = 'd'
      xlocationmpi = 0.
      ylocationmpi = 0.
      zlocationmpi = 0.
      pi = 3.141592654

      WRITE (*,*) 'Do you wish to overwrite the input dump file?'
      WRITE (*,*) '    (If you answer anything but y, the output dump'
      WRITE (*,*) '    will have the original name appended with .mpi)'
c      READ (*,88001) ianswer
      ianswer = 'n'
88001 FORMAT(A1)

      IF (imhd.EQ.idim) THEN
50    WRITE(*, 88002)
88002 FORMAT (/,' Which MHD variable did you evolve?', /,
     &   ' B     (flux density per unit volume)         : Bvol b)', /,
     &   ' B/rho (flux density per unit mass)           : Brho r)', /,
     &   ' Euler potentials (B = grad alpha x grad beta): eulr e)', /,
     &   ' Vector potential (B = curl A)                : vecp v)')
c      READ (*, 88003) varmhd
      varmhd = 'r'
88003 FORMAT (A4)
      IF (varmhd.EQ.'b' .or. varmhd.EQ.'B') varmhd = 'Bvol'
      IF (varmhd.EQ.'r') varmhd = 'Brho'
      IF (varmhd.EQ.'e') varmhd = 'eulr'
      IF (varmhd.EQ.'v') varmhd = 'vecp'
      IF (varmhd.NE.'Bvol' .and. varmhd.NE.'Brho'
     &    .and. varmhd.NE.'eulr' .and. varmhd.NE.'vecp') GOTO 50
      END IF

      WRITE (*,*) 'How many full dumps?'
c      READ (*,*) numfulldumps
      numfulldumps = 1
      IF (numfulldumps.GT.maxnumber) THEN
         WRITE (*,*) 'Filename array too small'
         STOP
      ENDIF

      WRITE (*,*) 'Enter names of dump files'
      DO i = 1, numfulldumps
c         READ (*,*) fullname(i)
         fullname(i) = 'MF50000'
      END DO

      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO

      DO iloop = 1, numfulldumps

c      WRITE (*,*) 'Enter name of (full) file to convert'
c      READ (*,*) fullname
         IF (ianswer.EQ.'y') THEN
            file1 = fullname(iloop)
         ELSE
            file1 = 'MF50002'
         ENDIF
          
         OPEN (11,FILE=fullname(iloop),FORM='unformatted')

         CALL rdump(11,ichkl,1)

         CLOSE (11)
         IF (ichkl.EQ.1) THEN
            WRITE (*,*) 'ERROR reading input file: ichkl ',ichkl
            STOP
         ENDIF

         umass = umassi
         udist = udisti
         utime = utimei
         umagfd = umagfdi
         dtmax = dtmaxdp

         IF (iproc.EQ.0) THEN
            CALL file
            OPEN (11,FILE=file1,FORM='unformatted')
         ENDIF

         print *,iproc,': Calling wdump ',npart

         CALL wdump(11)

         print *,iproc,': Done wdump'

      END DO

      CALL quit

      END


      SUBROUTINE quit
      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'

      CALL MPI_FINALIZE(ierr)
#endif
      STOP
      END


      SUBROUTINE endrun
      CALL quit
      END

      SUBROUTINE ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)

      STOP
      END

      SUBROUTINE insulate(ival, ntot, npart, dumxyzmh, f1vxyzu)

      STOP
      END


      SUBROUTINE labrun
      RETURN
      END
