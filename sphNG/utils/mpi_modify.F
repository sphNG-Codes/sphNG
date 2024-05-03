      PROGRAM mpi_modify
c
c--SUMMARY:
c     This utility converts an existing dump file to one or more
c     new dump files, each with a specified number of MPI blocks.
c     All output dumps will have the binary decomposition.
c     Written by TB May 2024.
c
c--COMPILE:
c     make mpi_modify mpimod=yes
c     if there is spare memory you can add openmp=yes (see note 4)
c     Do not compile with mpi=yes.
c
c--INPUTS:
c     (1) a serial or MPI dump (any number of blocks), with any mpitype.
c     (2) a file named mod_param.txt, of the following format,
c            -the input dump file name \newline
c            -the output dump file name \newline
c            -the number of output dumps required \newline
c            -a space separated list of the desired number MPI ranks for
c             each required output dump (binary numbers only).
c
c--OUTPUTS:
c     (1) file_names.txt, which contains a list of output dump files and
c         the numper of MPI ranks they are split into.
c     (2) The dump files in list mentioned above. They will all have
c         the binary decomposition type, mpitype='b'.
c
c--NOTES:
c     (1) idim must be at least equal to the total number of particles
c         in the input dump.
c     (2) While this method requires a lot of memory for large particle
c         sets, it is never compiled with mpi=yes. This means that the
c         idim limit is 2^31 (~2 billion), rather than typical limit of
c         ~20 million. This typical limit is imposed by the maximum
c         size of one MPI message.
c     (3) For idim > 85 million you need to change the type of imaxrec
c         in COMMONS/logun to INTEGER*8. This requires a compiler newer
c         than 2018. Tested with intel 2020.
c     (4) With openMP and large idim, the OMP_STACKSIZE needs to be
c         high. 1 billion particles works with, export OMP_STACKSIZE=4G
c         Note that this means that openMP will use 4 gigabytes of
c         additional memory for each openMP thread so do not set
c         OMP_NUM_THREADS too high! Testing on a 128 core node suggests
c         there is no point using more than 32 threads. Running in
c         serial may be the best option, lots of the time is spent in
c         serial I/O anyway.
c
#ifdef _OPENMP
      INCLUDE 'omp_lib.h'
#endif
      INCLUDE '../idim'
      INCLUDE '../COMMONS/mpimod'
      INCLUDE '../COMMONS/units'
      INCLUDE '../COMMONS/sort'
      INCLUDE '../COMMONS/part'
      INCLUDE '../COMMONS/phase'
      INCLUDE '../COMMONS/varmhd'
      INCLUDE '../COMMONS/files'
      INCLUDE '../COMMONS/logun'
      INCLUDE '../COMMONS/actio'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      INTEGER, ALLOCATABLE :: imap(:,:),numprocall(:)
      INTEGER, ALLOCATABLE :: isplitposkey(:),npartkey(:)
      INTEGER*2 numprocprint,nfile

      CALL getused(tstart)
      WRITE(*,*) "This is the mpi_modify utility"
#ifndef MPIMOD
      PRINT*, "cannot run mpi_modify without mpimod=yes"
      CALL quit(0)
#endif
#ifdef MPIALL
      PRINT*, "Do not compile this utility with mpi=yes"
      CALL quit(0)
#endif
      IF (istellarfeedback.GT.0 .AND. nptmass.GT.0) THEN
         WRITE (*,*) 'ERROR - Cannot start MPI calculation ',
     &        'from scalar calculation if it'
         WRITE (*,*) 'has stellar feedback and exisiting ',
     &        'sink particles'
         CALL quit(0)
      ENDIF

      iprint = istdout
      iread  = istdin

#ifdef _OPENMP
C$OMP PARALLEL
      IF (OMP_GET_THREAD_NUM() .eq. 0) THEN
         WRITE(*,*) "Using", OMP_GET_NUM_THREADS(), "openMP threads."
      END IF
C$OMP END PARALLEL
#endif
c
c--Set type of output domain decomposition
c
      mpitype = 'b'
c
c--Handle type of MHD
c     Hardcoded to varmhd='r', for other options see utils/mpi_reduce.F
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'r'
         PRINT*,"varmhd set to r, for other values change mpi_modify.F"
      END IF
c
c--Read parameters file
c
99003 FORMAT (A7)
      iparam = 20
      OPEN (iparam,FILE='mod_param.txt',STATUS='old',ACTION='read')
      READ (iparam,99003) file2
      READ (iparam,99003) file1
      READ (iparam,*) noutput
      ALLOCATE (numprocall(noutput))
      READ (iparam,*) (numprocall(i),i=1,noutput)
      CLOSE (iparam)
      PRINT*,"Making files with the following numprocs",numprocall
c
c--Check provided values of numproc are binary numbers and
c     find maximum number of splits.
c
      nsplits_max = 0
      DO i=1,noutput
         numproc = numprocall(i)
         xtemp = LOG(REAL(numproc))/LOG(2.0)
         IF (xtemp - INT(xtemp+1.0E-4).GT.1.0E-4) THEN
            WRITE (*,*) 'ERROR - numproc is not a power of two ',
     &           i,numproc,xtemp
            CALL quit(1)
         ENDIF
         nsplits = INT(xtemp+1.0E-4)
         nsplits_max = MAX(nsplits,nsplits_max)
      ENDDO
c
c--Read in the existing dump file
c
      CALL getused(t1)
      OPEN (idisk2,FILE=file2,FORM='unformatted')
      CALL rdump(idisk2,ichkl,1)
      IF (ichkl.EQ.1) THEN
         WRITE (*,*) 'ERROR reading input file: ichkl ',ichkl
         STOP
      ENDIF
      CLOSE (idisk2)
      nparttotorig = npart
      CALL getused(t2)
      print*,"rdump time",t2-t1,npart
c
c--Set the units to the units from the input dump file.
c
      umass = umassi
      udist = udisti
      utime = utimei
      umagfd = umagfdi
      dtmax = dtmaxdp
c
c--Alocate memory for binary decomposition storage arrays.
c
      nparttot = 0
      DO i = 1,nparttotorig
         IF (iphase(i) .GE. 0) THEN
            nparttot = nparttot + 1
         ENDIF
      END DO
      ALLOCATE(imap(nsplits_max,nparttot))
      nkey = 2**(nsplits_max+1)-1
      ALLOCATE(isplitposkey(nkey))
      ALLOCATE(npartkey(nkey))
c
c--Split particles onto their MPI ranks and save in imap array.
c
      CALL getused(t1)
      CALL binary_split(imap,isplitposkey,npartkey)
      CALL getused(t2)
      PRINT*,"binary_split, total time",t2-t1
c
c--Loop over all output files.
c
      OPEN (idisk3,FILE='file_names.txt')
      PRINT*,""
      PRINT*,"Now writing to files"
      DO nfile=1,noutput
         numproc = numprocall(nfile)
         numprocprint = numproc
c
c--Find the number of binary splits.
c
         xtemp = LOG(REAL(numproc))/LOG(2.0)
         nsplits = INT(xtemp+1.0E-4)
         PRINT*,""
         PRINT*,"Writing ",file1," numproc ",numprocprint
         PRINT*,': binary number of splits is ',nsplits
c
c--Open new output file and write numproc blocks to it.
c
         CALL file
         CALL getused(t1)
         CALL write_mpi_dump(nsplits,imap,isplitposkey,npartkey)
         CALL getused(t2)
         PRINT*,"File, ",file1," complete, total write time ",t2-t1

         WRITE (idisk3,*) "numproc ",numprocprint," in ",file1
      ENDDO
      CLOSE (idisk3)
      CALL getused(tend)
      PRINT*,"mpi_modify, total time",tend-tstart

      END PROGRAM mpi_modify

c
c--Define the "quit" subroutine that is used by sphNG subroutines.
c
      SUBROUTINE quit(i)
      STOP
      END

c
c--Dummy routines that allow sphNG subroutines to work properly.
c
      SUBROUTINE endrun
      CALL quit(0)
      END

c
c--This routine builds a map of all binary domain decompositions.
c
      SUBROUTINE binary_split(imap,isplitposkey,npartkey)

      INCLUDE 'idim'
      INCLUDE 'igrape'
      INCLUDE 'COMMONS/mpimod'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'

      INTEGER, INTENT(INOUT) :: imap(nsplits_max,nparttot)
      INTEGER, INTENT(INOUT) :: isplitposkey(nkey),npartkey(nkey)
      INTEGER, ALLOCATABLE :: iorder(:),list(:)

      ALLOCATE(iorder(nparttot))
      ALLOCATE(list(nparttot))

      PRINT*, "Entering binary_split ",nsplits_max,numproc,nkey
      imap = 0

      DO nb=1,nsplits_max
         CALL getused(t1)
         iposstart = 2**(nb-1)
         iposend = 2**(nb)-1
         ncdim = MOD(nb-1,3)+1
C$OMP PARALLEL DO SCHEDULE(static,1) default(none)
C$OMP& SHARED(npartkey,isplitposkey,imap)
C$OMP& SHARED(iphase,xyzmh)
C$OMP& SHARED(nb,ncdim,iposstart,iposend,nkey)
C$OMP& SHARED(nparttot,nparttotorig,npartactive)
C$OMP& PRIVATE(list,iorder)
C$OMP& PRIVATE(ipos,i,isplitpos,npproc)
         DO ipos=iposstart,iposend
            IF (nb .EQ. 1) THEN
c--First time around need to build list of active particles.
               isplitposkey(1) = 0
               isplitpos = 0
               npartkey(1) = nparttot
               npproc = nparttot
               npartactive = 0
               DO i = 1,nparttotorig
                  IF (iphase(i) .GE. 0) THEN
                     npartactive = npartactive + 1
                     list(npartactive) = i
                  ENDIF
               END DO
               IF (nparttot .NE. npartactive) THEN
                  PRINT*,"Active particle mismatch",npartactive,nparttot
                  STOP
               ENDIF
            ELSE
c
c--Other times copy part of the previous row of imap to get the list.
c
               isplitpos = isplitposkey(ipos)
               npproc = npartkey(ipos)
               DO i=1,npproc
                  list(i) = imap(nb-1,i+isplitpos)
               ENDDO
            ENDIF
c
c--Sort xyzmh along ncdim axis (i.e. x,y,z,x,y,z,...)
c
            CALL indexx3(npproc,list,xyzmh,ncdim,iorder)
            DO ipart=1,npproc
               imap(nb,ipart+isplitpos) = list(iorder(ipart))
            ENDDO
            isplitposkey(2*ipos) = isplitpos
            isplitposkey(2*ipos+1) = isplitpos + (npproc / 2)
            npartkey(2*ipos) = npproc / 2
            npartkey(2*ipos+1) = npproc - (npproc / 2)
            IF (ipos .GT. nkey) THEN
               PRINT*,"ipos > nkey",ipos,nkey
               STOP
            ENDIF
         END DO
C$OMP END PARALLEL DO
         CALL getused(t2)
         PRINT*,nb,"split time",t2-t1
      ENDDO
      DEALLOCATE(iorder)
      DEALLOCATE(list)
c
c--Check for bad values in imap array.
c
C$OMP PARALLEL DO SCHEDULE(static,1) default(none)
C$OMP& SHARED(imap)
C$OMP& SHARED(nsplits_max,nparttot)
C$OMP& PRIVATE(nb,i)
      DO nb=1,nsplits_max
         DO i=1,nparttot
            IF (imap(nb,i) .LE. 0) THEN
               PRINT*,"imap zero",nb,i,imap(nb,i)
               STOP
            ENDIF
         ENDDO
      ENDDO
C$OMP END PARALLEL DO

      END SUBROUTINE binary_split

c
c--This routine writes a binary MPI decomposition to a dump given imap.
c
      SUBROUTINE write_mpi_dump(nsplits,imap,isplitposkey,npartkey)

      INCLUDE 'idim'
      INCLUDE 'COMMONS/mpimod'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/rbnd'

      INTEGER, INTENT(IN) :: nsplits
      INTEGER, INTENT(INOUT) :: imap(nsplits_max,nparttot)
      INTEGER, INTENT(INOUT) :: isplitposkey(nkey),npartkey(nkey)
c
c--With a bit of thinking this could be done without the inner loop.
c     However, the loops are not expensive and this approach
c     is the same as in the split() routine in rdump.F.
c
      DO iproc=0,numproc-1
         ipos = 1
         DO nb=1,nsplits
            ipris = iproc
            ndiv = numproc
            DO i=1,nsplits-nb
               ndiv = ndiv/2
               IF (ipris .GE. ndiv) ipris = ipris-ndiv
            ENDDO
c--SPLIT x,y,z,x,y,z etc...
            ncdim = MOD(nb-1,3)+1
c--Decide which half of particles needed for this iproc.
            IF (ipris .LT. ndiv/2) THEN
               ipos = 2*ipos
               PRINT*,"Top half",nb,ipos
            ELSE
               ipos = 2*ipos+1
               PRINT*,"Bottom half",nb,ipos
            ENDIF
         ENDDO
         isplitpos = isplitposkey(ipos)
         npart = npartkey(ipos)
c--Check no block has too many particles.
         IF (npart .GT. 20000000 .AND. numproc .GT. 1) THEN
            PRINT*,"npart too large to evolve with MPI",npart
            STOP
         ENDIF
         DO ipart=1,npart
            isort(ipart) = imap(nsplits,ipart+isplitpos)
         ENDDO
c
c--Write this block to file.
c
         nlistinactive = 0
         CALL getused(t1)
         PRINT *,iproc,': Calling wdump ',npart,file1
         CALL wdump(idisk1,iproc,numproc)
         PRINT *,iproc,': Done wdump'
         CALL getused(t2)
         PRINT*,"block ",iproc+1," wdump time ",t2-t1
      ENDDO

      END SUBROUTINE write_mpi_dump
