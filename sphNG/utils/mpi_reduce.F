      PROGRAM mpi_reduce_b
c
c--Converts an MPI dump file (i.e. multiple blocks) that has type "b"
c     (binary) domain decomposition over N domains, to one that has
c     half the number of domains (i.e. N/2, for example from 128->64).
c
c--NOTE: The input dump file MUST have type "b" domain decomposition,
c     otherwise the new domains will have very weird shapes!  There is
c     currently NO WAY to check that the incoming file actually has "b".
c
c--NOTE: idim can be set to the smallest number that will hold all the
c     particles in one of the NEW domains.  In this respect, it requires
c     a lot less memory than balance.F requires to change from one
c     domain type to another.
c
c     "idim" could be larger than 20,000,000 because this code does not
c     use MPI to transfer particles between nodes (it just adds together
c     the particles from the appropriate blocks in the input file).
c     However, the sphNG code itself cannot handle more than about
c     idim=20,000,000 for each MPI process otherwise the MPI data 
c     transfers will exceed 2^31 (INTEGER*4) in size and the subroutine
c     balance.F (and probably others too) will not work properly.
c
c     The input files are called "BIN0000" and "ifile".  The output
c     file is "BIN0002".
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
      fullname = 'BIN0000'

      isort(1:idim) = i
      iorig(1:idim) = i
c
c--Set output file (number is one less than it will actually be)
c
      file1 = 'BIN0001'
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
