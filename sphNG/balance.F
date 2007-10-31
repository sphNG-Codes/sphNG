      SUBROUTINE balance
c************************************************************
c                                                           *
c  Subroutine to load balance job across MPI processes      *
c                                                           *
c     Code written by MRB (23/10/2007).                     *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/f2'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/bodys'

      DIMENSION listtrans(idim), numlist(nummaxproc)
c
c--Make derived types
c
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(isizealphaMM, MPI_REAL4, 
     &     ialphaREAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(isizealphaMM+1, MPI_REAL4, 
     &     ialphaREAL41, ierr)

      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)
      CALL MPI_TYPE_COMMIT(ialphaREAL4,ierr)
      CALL MPI_TYPE_COMMIT(ialphaREAL41,ierr)
c
c--Find particles that need to be transferred to each MPI process.
c     Returns an INTEGER*2 array with the process number that each
c     particle belongs to and the number of particles that need to be moved
c     to each process.
c
      CALL maketransferlist(npart,numlist,listtrans,xyzmh,iphase)

      DO i = 1, npart
         IF (listtrans(i).EQ.iproc) write (40+iproc,*) xyzmh(1,i),
     &              xyzmh(2,i),xyzmh(3,i)
         IF (listtrans(i).NE.iproc) write (50+iproc,*) xyzmh(1,i),
     &              xyzmh(2,i),xyzmh(3,i)
      END DO

#ifdef MPIDEBUGT
      DO i = 1, numproc
         print *,iproc,': BELONG ',i-1,numlist(i)
      END DO
#endif

c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
      npartnew = npart



c      GOTO 100



      DO ii = 1, numproc - 1
         iahead = MOD(iproc+ii,numproc)
         ibehind = MOD(numproc+iproc-ii,numproc)
c
c--Make list of particles going to process "iahead"
c
         IF (numlist(iahead+1).GT.0) THEN
            ipos = 0
            DO j = 1, npart
               IF (iahead.EQ.listtrans(j)) THEN
                  ipos = ipos + 1
                  lsendlist(ipos) = j - 1
                  IF (ipos.GE.numlist(iahead+1)) GOTO 10
               ENDIF
            END DO
 10         CONTINUE


            DO i = 1, numlist(iahead+1)
               write (30+iproc,*) xyzmh(1,lsendlist(i)+1),
     &              xyzmh(2,lsendlist(i)+1),xyzmh(3,lsendlist(i)+1)
            END DO

         ENDIF
c
c--Make derived types required for data transfer:
c     Need to send: xyzmh, vxyzu, alphaMM(1,i), gradhs, isteps, iphase
c        it0, it1, it2, f1vxyzu, f2vxyzu, rho, dumrho
c     For sinks: all sink arrays: listpm, spinx,y,z, angadd, spinad
c     For RT: ekcle
c     For MHD: Bevolxyz, Bxyz, alphaMM(isizealphaMM,i),
c        f1Bxyz, f2Bxyz
c
c     Other points to note: ??dumxyzmh, ??dumvxyzu, it1bin, it2bin,
c        nlstbins, listbins, iscurrent
c
c     Also need to define iorig() and isort()
c
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead), 1, 
c     &       lsendlist, i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &        lsendlist, i5REAL8, indexMPI5, ierr)
         CALL MPI_TYPE_COMMIT(indexMPI5,ierr)
c
c--Send to node ahead, receive from node behind
c
#ifdef MPIDEBUGT
         print *,iproc,': Transferring ',numlist(iahead+1),
     &        ' particles to ',iahead,' Receiving from ',
     &        ibehind
#endif

         CALL MPI_SENDRECV(xyzmh,1,indexMPI5,iahead,100,
     &        xyzmh(1,npartnew+1), idim, i5REAL8, ibehind,
     &        100, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
         igotnumber = inumber

#ifdef MPIDEBUGT
         print *,iproc,': Received xyzmh transfer of ',inumber,
     &        ' and wants to send ',numlist(iahead+1)
#endif

ccc         IF (igotnumber.NE.0 .OR. numlist(iahead+1).NE.0) THEN

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, i4REAL8, indexMPI4, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, i4REAL8, indexMPI4, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI4,ierr)
            
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, i3REAL8, indexMPI3, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, i3REAL8, indexMPI3, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI3,ierr)
            
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, MPI_INTEGER, indexMPI_INT, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, MPI_INTEGER, indexMPI_INT, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI_INT,ierr)
            
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, MPI_INTEGER1, indexMPI_INT1, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, MPI_INTEGER1, indexMPI_INT1, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)
            
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, i2REAL4, indexMPI2, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, i2REAL4, indexMPI2, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI2,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, MPI_REAL4, indexMPI1, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, MPI_REAL4, indexMPI1, ierr)
            CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, ialphaREAL4, indexMPIalpha, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, ialphaREAL4, indexMPIalpha, ierr)
            CALL MPI_TYPE_COMMIT(indexMPIalpha,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(iahead+1), 1, 
c     &       lsendlist, ialphaREAL41, indexMPIalpha1, ierr)
            CALL MPI_TYPE_INDEXED(numlist(iahead+1), lblocklengths, 
     &           lsendlist, ialphaREAL41, indexMPIalpha1, ierr)
            CALL MPI_TYPE_COMMIT(indexMPIalpha1,ierr)

            
            CALL MPI_SENDRECV(vxyzu,1,indexMPI4,iahead,101,
     &           vxyzu(1,npartnew+1), idim, i4REAL8, ibehind,
     &           101, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received vxyzu transfer of ',inumber
#endif

            CALL MPI_SENDRECV(f1vxyzu,1,indexMPI4,iahead,102,
     &           f1vxyzu(1,npartnew+1), idim, i4REAL8, ibehind,
     &           102, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received f1vxyzu transfer of ',inumber
#endif

            CALL MPI_SENDRECV(f2vxyzu,1,indexMPI4,iahead,103,
     &           f2vxyzu(1,npartnew+1), idim, i4REAL8, ibehind,
     &           103, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received f2vxyzu transfer of ',inumber
#endif

            CALL MPI_SENDRECV(isteps,1,indexMPI_INT,iahead,104,
     &           isteps(npartnew+1), idim, MPI_INTEGER, ibehind,
     &           104, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received isteps transfer of ',inumber,
     &           isteps(3033),isteps(3034),isteps(3035)
#endif

            CALL MPI_SENDRECV(it0,1,indexMPI_INT,iahead,105,
     &           it0(npartnew+1), idim, MPI_INTEGER, ibehind,
     &           105, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received it0 transfer of ',inumber
#endif

            CALL MPI_SENDRECV(it1,1,indexMPI_INT,iahead,106,
     &           it1(npartnew+1), idim, MPI_INTEGER, ibehind,
     &           106, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received it1 transfer of ',inumber
#endif

            CALL MPI_SENDRECV(it2,1,indexMPI_INT,iahead,107,
     &           it2(npartnew+1), idim, MPI_INTEGER, ibehind,
     &           107, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received it2 transfer of ',inumber
#endif

            CALL MPI_SENDRECV(iphase,1,indexMPI_INT1,iahead,108,
     &           iphase(npartnew+1), idim, MPI_INTEGER1, ibehind,
     &           108, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER1, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received iphase transfer of ',inumber
#endif

            CALL MPI_SENDRECV(rho,1,indexMPI1,iahead,109,
     &           rho(npartnew+1), idim, MPI_REAL4, ibehind,
     &           109, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received rho transfer of ',inumber
#endif
            DO ipart = npartnew + 1, npartnew + igotnumber
               dumrho(ipart) = rho(ipart)
            END DO

            CALL MPI_SENDRECV(poten,1,indexMPI1,iahead,110,
     &           poten(npartnew+1), idim, MPI_REAL4, ibehind,
     &           110, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received poten transfer of ',inumber
#endif

            IF (nlmax.EQ.1) THEN
               CALL MPI_SENDRECV(gradhs,1,indexMPI2,iahead,111,
     &              gradhs(1,npartnew+1), idim, i2REAL4, ibehind,
     &              111, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i2REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit
               ENDIF

#ifdef MPIDEBUGT
               print *,iproc,': Received gradhs transfer of ',inumber
#endif
            ENDIF

            CALL MPI_SENDRECV(alphaMM,1,indexMPIalpha,iahead,112,
     &           alphaMM(1,npartnew+1), idim, ialphaREAL4, ibehind,
     &           112, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, ialphaREAL4, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received alphaMM transfer of ',inumber
#endif

            CALL MPI_SENDRECV(f1ha,1,indexMPIalpha1,iahead,113,
     &           f1ha(1,npartnew+1), idim, ialphaREAL41, ibehind,
     &           113, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, ialphaREAL41, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received f1ha transfer of ',inumber
#endif

            CALL MPI_SENDRECV(f2ha,1,indexMPIalpha1,iahead,114,
     &           f2ha(1,npartnew+1), idim, ialphaREAL41, ibehind,
     &           114, MPI_COMM_WORLD, istatus, ierr)
            CALL MPI_GET_COUNT(istatus, ialphaREAL41, inumber, ierr)
            IF (inumber.NE.igotnumber) THEN
               WRITE (*,*) 'ERROR - balance igotnumber'
               CALL quit
            ENDIF

#ifdef MPIDEBUGT
            print *,iproc,': Received f2ha transfer of ',inumber
#endif

            IF (imhd.EQ.idim) THEN
               CALL MPI_SENDRECV(Bevolxyz,1,indexMPI3,iahead,115,
     &              Bevolxyz(1,npartnew+1), idim, i3REAL8, ibehind,
     &              115, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i3REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit
               ENDIF

#ifdef MPIDEBUGT
               print *,iproc,': Received Bevol transfer of ',inumber
#endif

               CALL MPI_SENDRECV(Bxyz,1,indexMPI3,iahead,116,
     &              Bxyz(1,npartnew+1), idim, i3REAL8, ibehind,
     &              116, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i3REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit
               ENDIF
                  
#ifdef MPIDEBUGT
               print *,iproc,': Received Bxyz transfer of ',inumber
#endif

               CALL MPI_SENDRECV(f1Bxyz,1,indexMPI3,iahead,117,
     &              f1Bxyz(1,npartnew+1), idim, i3REAL8, ibehind,
     &              117, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i3REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit
               ENDIF
                  
#ifdef MPIDEBUGT
               print *,iproc,': Received f1Bxyz transfer of ',inumber
#endif

               CALL MPI_SENDRECV(f2Bxyz,1,indexMPI3,iahead,118,
     &              f2Bxyz(1,npartnew+1), idim, i3REAL8, ibehind,
     &              118, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i3REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit
               ENDIF

#ifdef MPIDEBUGT
               print *,iproc,': Received f2Bxyz transfer of ',inumber
#endif
            ENDIF
            CALL MPI_TYPE_FREE(indexMPIalpha1,ierr)
            CALL MPI_TYPE_FREE(indexMPIalpha,ierr)
            CALL MPI_TYPE_FREE(indexMPI1,ierr)
            CALL MPI_TYPE_FREE(indexMPI2,ierr)
            CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
            CALL MPI_TYPE_FREE(indexMPI_INT,ierr)
            CALL MPI_TYPE_FREE(indexMPI3,ierr)
            CALL MPI_TYPE_FREE(indexMPI4,ierr)
ccc         ENDIF
         CALL MPI_TYPE_FREE(indexMPI5,ierr)

         npartnew = npartnew + igotnumber
         nactive = nactive + igotnumber

#ifdef MPIDEBUGT
         print *,iproc,': Finished transfer to ',iahead
         print *,iproc,': New value of npart = ',npartnew, npart
#endif
      END DO

 100  CALL MPI_TYPE_FREE(ialphaREAL41,ierr)
      CALL MPI_TYPE_FREE(ialphaREAL4,ierr)
      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)
c
c--End circular send/receive
c
c
c--Need to set particle indexxing
c
      DO i = npart+1, npartnew
         iscurrent(i) = .FALSE.
         isort(i) = i
         iorig(i) = i
      END DO
c
c--Now remove particles that have been transferred to another MPI process
c
      numtransfer = 0
      DO i = 1, numproc
         numtransfer = numtransfer + numlist(i)
      END DO
      numtransfer = numtransfer - numlist(iproc+1)
#ifdef MPIDEBUGT
      print *,iproc,': Total number transferred ',numtransfer
#endif
      IF (numtransfer.GT.0) THEN
         DO ipart = 1, npart
            IF (listtrans(ipart).GE.0 .AND. 
     &           listtrans(ipart).NE.iproc) THEN
               iphase(ipart) = -1
               nactive = nactive - 1
            ENDIF
         END DO
      ENDIF
      npart = npartnew
      n1 = npart

#ifdef MPIDEBUGT
      print *,iproc,': New npart and nactive ',npart,nactive
#endif
c
c--Need to reset lists of particles in each timestep bin
c
      DO i = 1, nbinmax
         nlstbins(i) = 0
         it2bin(i) = 2**i
      END DO

      xlog2 = 0.30103

c      goto 200

      DO i = 1, npart
         IF (iphase(i).NE.-1) THEN
            it0(i) = 0
            it2(i) = isteps(i)

            ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
            IF (ibin.GT.nbinmax) THEN
               WRITE (*,*) 'ERROR - ibin.GT.nbinmax B'
               WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax B'
               CALL quit
            ENDIF
            nlstbins(ibin) = nlstbins(ibin) + 1
            listbins(nlstbins(ibin),ibin) = i
            IF (it2bin(ibin).NE.it2(i)) THEN
               WRITE (*,*) 'ERROR - it2bin B'
               CALL quit
            ENDIF
         ENDIF
      END DO
c
c--End of MPI-only section
c
#endif


c 200  CALL quit



      RETURN
      END

c---------------------------------------------------------------------------

      SUBROUTINE maketransferlist(npart,numlist,listtrans,xyzmh,iphase)
c************************************************************
c                                                           *
c  Subroutine to load balance job across MPI processes      *
c                                                           *
c     Code written by MRB (23/10/2007).                     *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'

      DIMENSION listtrans(idim), numlist(nummaxproc), xyzmh(5,idim)
      INTEGER*1 iphase(idim)
      INTEGER classify_cartesian, classify_disc, classify_sphere
c
c--Identify which MPI process each particle should belong to
c
      DO i = 1, numproc
         numlist(i) = 0
      END DO
      DO ipart = 1, npart
         IF (iphase(ipart).GE.0) THEN
            IF (mpitype.EQ.'c') THEN
               ibelong = classify_cartesian(ipart,xyzmh)
               ibelong1 = ibelong + 1
               numlist(ibelong1) = numlist(ibelong1) + 1
               listtrans(ipart) = ibelong
            ELSEIF (mpitype.EQ.'d') THEN
               ibelong = classify_disc(ipart,xyzmh)
               ibelong1 = ibelong + 1
               numlist(ibelong1) = numlist(ibelong1) + 1            
               listtrans(ipart) = ibelong
            ELSEIF (mpitype.EQ.'s') THEN
               ibelong = classify_sphere(ipart,xyzmh)
               ibelong1 = ibelong + 1
               numlist(ibelong1) = numlist(ibelong1) + 1
               listtrans(ipart) = ibelong
            ENDIF
         ELSE
            listtrans(ipart) = -1
         ENDIF
      END DO
c
c--End of MPI-only section
c
#endif
      RETURN
      END

c---------------------------------------------------------------------------

      FUNCTION classify_disc(ipart,xyzmh)
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
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'

      DIMENSION xyzmh(5,idim)

      INTEGER classify_disc

      ival = 0

#ifdef MPI
      IF (xyzmh(3,ipart).LT.0.0) THEN
         ival = 1
      ELSE
         ival = 0
      ENDIF

      phiincrement = 2.0*pi/(numproc/2)

      IF (xyzmh(1,ipart).EQ.0. .AND. xyzmh(2,ipart).EQ.0.) THEN
         ival = 0
      ELSE
         ival = ival*numproc/2 +
     &     INT((pi + ATAN2(xyzmh(2,ipart),xyzmh(1,ipart)))/phiincrement)
      ENDIF
#endif

      classify_disc = ival

      RETURN
      END


c---------------------------------------------------------------------------

      FUNCTION classify_sphere(ipart,xyzmh)
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
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'

      DIMENSION xyzmh(5,idim)

      INTEGER classify_sphere

      ival = 0

      classify_sphere = ival

      RETURN
      END


c---------------------------------------------------------------------------

      FUNCTION classify_cartesian(ipart,xyzmh)
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
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpidebug'
#endif
      INCLUDE 'COMMONS/physcon'

      DIMENSION xyzmh(5,idim)

      INTEGER classify_cartesian

      ival = 0


      classify_cartesian = ival

      RETURN
      END

c==========================================================================
