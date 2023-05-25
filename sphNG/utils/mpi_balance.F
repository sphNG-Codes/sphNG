      PROGRAM mpi_balance_d_to_b
c
c--Converts an MPI dump file (i.e. multiple blocks) that was created
c     with type "d" (disc) domain decomposition to type "b" (binary)
c     domain decomposition.
c
c     NOTE: idim must be set to be at LEAST twice as large as the
c     maximum number of particles per domain in the input dump file
c     to allow particles to be moved between MPI domains.
c     
c     NOTE: idim CANNOT be larger than about 20,000,000 otherwise the
c     MPI data transfers will exceed 2^31 (INTEGER*4) in size and the
c     subroutine balance.F will not work properly.
c
c     This can be done simply by running sphNG to the first dump time
c     requesting "b" in the i-file.  However, in order to do this,
c     the code has to be compiled with "idim" twice as large as that
c     required to hold all the particles to allow particles to be
c     moved between the MPI processes.  However, this means that the
c     memory required can be large.
c
c     This programme gives a start alone code that can be used to
c     convert an existing dump file with "d" decomposition to "b"
c     decomposition (the number of MPI processes and domains remains 
c     the same).
c
c     The input files are called "DUM0000" and "ifile".  The output
c     file is "DUM0002".
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
      INCLUDE '../COMMONS/active'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      CHARACTER*7 fullname

#ifdef MPIALL
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)

      lblocklengths = 1  ! array operation
#endif

      inname = 'ifile'
      iprint = istdout 
      iread  = istdin

      CALL options

      nfullstep = 1
c
c--Set type of output domain decomposition
c
      mpitype = 'b'
      xlocationmpi = 0.
      ylocationmpi = 0.
      zlocationmpi = 0.
c
c--Handle type of MHD
c
      IF (imhd.EQ.idim) THEN
 50      WRITE(*, 88002)
88002    FORMAT (/,' Which MHD variable did you evolve?', /,
     &   ' B     (flux density per unit volume)         : Bvol b)', /,
     &   ' B/rho (flux density per unit mass)           : Brho r)', /,
     &   ' Euler potentials (B = grad alpha x grad beta): eulr e)', /,
     &   ' Vector potential (B = curl A)                : vecp v)')
c      READ (*, 88003) varmhd
         varmhd = 'r'
88003    FORMAT (A4)
         IF (varmhd.EQ.'b' .or. varmhd.EQ.'B') varmhd = 'Bvol'
         IF (varmhd.EQ.'r') varmhd = 'Brho'
         IF (varmhd.EQ.'e') varmhd = 'eulr'
         IF (varmhd.EQ.'v') varmhd = 'vecp'
         IF (varmhd.NE.'Bvol' .and. varmhd.NE.'Brho'
     &        .and. varmhd.NE.'eulr' .and. varmhd.NE.'vecp') GOTO 50
      END IF
c
c--Input filename
c
      fullname = 'DUM0000'

      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO
c
c--Set output file (number is one less than it will actually be)
c
      file1 = 'DUM0001'
c
c--Read in the existing dump file
c
      OPEN (11,FILE=fullname,FORM='unformatted')

      CALL rdump(11,ichkl,1)

      CLOSE (11)
      IF (ichkl.EQ.1) THEN
         WRITE (*,*) 'ERROR reading input file: ichkl ',ichkl
         STOP
      ENDIF
c
c--Set the units to the units from the input dump file.
c
      umass = umassi
      udist = udisti
      utime = utimei
      umagfd = umagfdi
      dtmax = dtmaxdp
c
c--Open output file
c
      IF (iproc.EQ.0) THEN
         CALL file
         OPEN (11,FILE=file1,FORM='unformatted')
      ENDIF
c
c--Change the domain decomposition.
c--Need to define nactive for balance 'b' to work properly.
c     
      nactive = npart

      CALL balance(0)
c
c--Write out new dump file
c
      print *,iproc,': Calling wdump ',npart

      CALL wdump(11)

      print *,iproc,': Done wdump'

      CALL quit(0)

      END
c
c--Define the "quit" subroutine that is used by sphNG subroutines.
c
      SUBROUTINE quit(i)
      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'

      CALL MPI_FINALIZE(ierr)
#endif
      STOP
      END
c
c--Dummy routines that allow sphNG subroutines to work properly.
c
      SUBROUTINE endrun
      CALL quit(0)
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
