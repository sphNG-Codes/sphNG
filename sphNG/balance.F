      SUBROUTINE balance
c************************************************************
c                                                           *
c  Subroutine to load balance job across MPI processes      *
c                                                           *
c     Code written by MRB (23/10/2007).                     *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
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
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/accnum'
      INCLUDE 'COMMONS/accurpt'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/radsink'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/dustimplicit'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/recor'

      DIMENSION itempsort(idim)
c      EQUIVALENCE (itempsort, next1)

      DIMENSION listtrans(idim), numlist(nummaxproc),numrecv(nummaxproc)
      DIMENSION listptmass(iptdim), listpmtrans(iptdim)

      CHARACTER*13 filename_split
      CHARACTER*3 idumpnumber
      CHARACTER*3 iprocnumber

      IF (iproc.EQ.0) WRITE (iprint,*) 'BALANCING'
c
c--Make derived types
c
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL8, i5REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(imhdevol, MPI_REAL8, imhdevolREAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL8, i3REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL8, i2REAL8, ierr)
      CALL MPI_TYPE_CONTIGUOUS(3, MPI_REAL4, i3REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_REAL4, i2REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(nheatingISR, MPI_REAL4, iHREAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(nchemistry, MPI_REAL4, iCREAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(5, MPI_REAL4, i5REAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(isizealphaMM, MPI_REAL4, 
     &     ialphaREAL4, ierr)
      CALL MPI_TYPE_CONTIGUOUS(isizealphaMM+1, MPI_REAL4, 
     &     ialphaREAL41, ierr)

      CALL MPI_TYPE_COMMIT(i5REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8,ierr)
      CALL MPI_TYPE_COMMIT(imhdevolREAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL8,ierr)
      CALL MPI_TYPE_COMMIT(i3REAL4,ierr)
      CALL MPI_TYPE_COMMIT(i2REAL4,ierr)
      CALL MPI_TYPE_COMMIT(iHREAL4,ierr)
      CALL MPI_TYPE_COMMIT(iCREAL4,ierr)
      CALL MPI_TYPE_COMMIT(i5REAL4,ierr)
      CALL MPI_TYPE_COMMIT(ialphaREAL4,ierr)
      CALL MPI_TYPE_COMMIT(ialphaREAL41,ierr)

      IF (encal.EQ.'r') THEN
         CALL MPI_TYPE_CONTIGUOUS(nmaxexposed, MPI_INTEGER, iNINT, ierr)
         CALL MPI_TYPE_COMMIT(iNINT,ierr)
      ENDIF

      IF (istellar_ionisation) THEN
         CALL MPI_TYPE_CONTIGUOUS(4*iptdim,MPI_REAL8,istellarREAL8,ierr)
         CALL MPI_TYPE_COMMIT(istellarREAL8,ierr)
      ENDIF
c
c--Find particles that need to be transferred to each MPI process.
c     Returns an INTEGER*2 array with the process number that each
c     particle belongs to and the number of particles that need to be moved
c     to each process.
c
      CALL maketransferlist(npart,numlist,listtrans,xyzmh,iphase)

      numlistmax = 0
      DO ii = 1, numproc
         IF (ii.NE.iproc) numlistmax = MAX(numlistmax,numlist(ii))
         numrecv(ii) = 0
      END DO

c      DO i = 1, npart
c         IF (listtrans(i).EQ.iproc) write (40+iproc,*) xyzmh(1,i),
c     &              xyzmh(2,i),xyzmh(3,i)
c         IF (listtrans(i).NE.iproc) write (50+iproc,*) xyzmh(1,i),
c     &              xyzmh(2,i),xyzmh(3,i)
c      END DO

#ifdef MPIDEBUGT
      DO i = 1, numproc
         print *,iproc,': BELONG ',i-1,numlist(i)
      END DO
#endif

c
c--IMPLEMENTATION USING CIRCULAR SEND_RECV
c
      npartnew = npart
      nptmassnew = nptmass
c
c--Need to send number of particles to be sent (so that other process knows
c     whether to expect information or not).
c
      inumtoget = 0
      numrecvtotal = 0
      DO ii = 1, numproc - 1
         iahead = MOD(iproc+ii,numproc)
         ibehind = MOD(numproc+iproc-ii,numproc)

         CALL MPI_SENDRECV(numlist(iahead+1),1,MPI_INTEGER,iahead,120,
     &        numrecv(ibehind+1), 1, MPI_INTEGER, ibehind,
     &        120, MPI_COMM_WORLD, istatus, ierr)
         CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
         IF (inumber.NE.1) THEN
            WRITE (*,*) 'ERROR - balance inumber.NE.1 ',inumber,iahead,
     &           ibehind
            CALL quit(1)
         ENDIF
         numrecvtotal = numrecvtotal + numrecv(ibehind+1)
         IF (numrecv(ibehind+1).GT.0) inumtoget = inumtoget + 1
      END DO
#ifdef MPIDEBUGT
      print *,iproc,': inumtoget ',inumtoget,(numrecv(i),i=1,numproc)

      print *,iproc,': test values ',npartnew,numrecvtotal,idim
#endif

      IF (npartnew + numrecvtotal.GT.idim) THEN
         WRITE (*,*) 'ERROR - npartnew + numrecvtotal.GT.idim ',
     &        npartnew, numrecvtotal, idim
         CALL quit(1)
      ENDIF

#ifdef MPIDEBUGT
      print *,iproc,': test done'
#endif
c
c--Now transfer particles
c
      inumgot = 0
      DO ii = 0, numproc - 1
         IF (iproc.EQ.ii) THEN
            IF (numlistmax.GT.0) THEN
c
c--Otherwise this process does not need to send any particles!
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
                     IF (numlist(j+1).GT.0) THEN
c
c--Make list of particles going to process "j"
c
                        ipos = 0
                        npttrans = 0
                        DO ipart = 1, npart
                           IF (j.EQ.listtrans(ipart)) THEN
                              ipos = ipos + 1
                              lsendlist(ipos) = ipart - 1
                              nlistinactive = nlistinactive + 1
                              listinactive(nlistinactive) = ipart
c                                 IF (ipos.LE.5) print *,'To ',
c     &                              ipart,iunique(iorig(ipart))
                              IF (iphase(ipart).GT.0 .AND. 
     &                             iphase(ipart).LT.10) THEN
                                 npttrans = npttrans + 1
                                 DO k = 1, nptmass
                                    IF (listpm(k).EQ.ipart) THEN
                                       listptmass(npttrans) = k - 1
                                       listpmtrans(npttrans) = ipos
                                       GOTO 5
                                    ENDIF
                                 END DO
                              ENDIF
 5                            CONTINUE
                              IF (ipos.GE.numlist(j+1)) GOTO 10
                           ENDIF
                        END DO
 10                     CONTINUE

c                        DO i = 1, numlist(j+1)
c                          write (30+iproc,*) xyzmh(1,lsendlist(i)+1),
c     &                          xyzmh(2,lsendlist(i)+1),
c     &                          xyzmh(3,lsendlist(i)+1)
c                        END DO
c
c--Make derived types required for data transfer:
c     Need to send: xyzmh, vxyzu, alphaMM(1,i), gradhs, isteps, iphase
c        it0, it1, it2, f1vxyzu, f2vxyzu, rho, dumrho
c     For sinks: all sink arrays: listpm, spinx,y,z, angadd, spinad
c     For RT: ekcle
c     For MHD: Bevolxyz, Bxyz, alphaMM(isizealphaMM,i),
c        f1Bxyz, f2Bxyz, divcurlB
c
c     Other points to note: ??dumxyzmh, ??dumvxyzu, it1bin, it2bin,
c        nlstbins, listbins, iscurrent
c
c     Also need to define iorig() and isort()
c
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j), 1, 
c     &       lsendlist, i5REAL8, indexMPI5, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i5REAL8, 
     &                       indexMPI5, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI5,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, i4REAL8, indexMPI4, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i4REAL8, 
     &                       indexMPI4, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI4,ierr)

                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, imhdevolREAL8, 
     &                       indexMPIimhdevol, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIimhdevol,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, i3REAL8, indexMPI3, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i3REAL8, 
     &                       indexMPI3, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI3,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, i2REAL8, indexMPI2, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i2REAL8, 
     &                       indexMPI2, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI2,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, MPI_REAL8, indexMPIR, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, MPI_REAL8, 
     &                       indexMPIR, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIR,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, MPI_INTEGER, indexMPI_INT, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, MPI_INTEGER, 
     &                       indexMPI_INT, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI_INT,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, MPI_INTEGER1, indexMPI_INT1, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, MPI_INTEGER1, 
     &                       indexMPI_INT1, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI_INT1,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, MPI_LOGICAL, indexMPI_L, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, MPI_LOGICAL, 
     &                       indexMPI_L, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI_L,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, i3REAL4, indexMPI34, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i3REAL4, 
     &                       indexMPI34, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI34,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, i2REAL4, indexMPI24, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i2REAL4, 
     &                       indexMPI24, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI24,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, iHREAL4, indexMPIH4, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, iHREAL4, 
     &                       indexMPIH4, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIH4,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, iCREAL4, indexMPIC4, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, iCREAL4, 
     &                       indexMPIC4, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIC4,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, i5REAL4, indexMPI54, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, i5REAL4, 
     &                       indexMPI54, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI54,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, MPI_REAL4, indexMPI1, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, MPI_REAL4, 
     &                       indexMPI1, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI1,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, ialphaREAL4, indexMPIalpha, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, ialphaREAL4, 
     &                       indexMPIalpha, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIalpha,ierr)

c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, ialphaREAL41, indexMPIalpha1, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, ialphaREAL41, 
     &                       indexMPIalpha1, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIalpha1,ierr)

                        IF (istellar_ionisation) THEN
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, istellarREAL8, indexMPIstellar, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, lsendlist, istellarREAL8, 
     &                       indexMPIstellar, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPIstellar,ierr)
                        ENDIF
c
c--NOTE: itempsort is used to send iunique indexes of particles which 
c     are not sorted like all the other arrays
c
                        IF (numlist(j+1).GT.idim/2) THEN
                           WRITE (*,*) 'ERROR - itempsort not large',
     &                          ' enough ',numlist(j+1),j
                           CALL quit(1)
                        ENDIF
                        DO i = 1, numlist(j+1)
                           itempsort(i) = iorig(lsendlist(i) + 1) - 1
                        END DO
c         CALL MPI_TYPE_CREATE_INDEXED_BLOCK(numlist(j+1), 1, 
c     &       lsendlist, MPI_INTEGER8, indexMPI_INT8, ierr)
                        CALL MPI_TYPE_INDEXED(numlist(j+1), 
     &                       lblocklengths, itempsort, MPI_INTEGER8, 
     &                       indexMPI_INT8, ierr)
                        CALL MPI_TYPE_COMMIT(indexMPI_INT8,ierr)
c
c--Now do sends
c
#ifdef MPIDEBUGT
                        print *,iproc,': Transferring ',numlist(j+1),
     &                        ' particles to ',j
#endif

                        CALL MPI_SEND(xyzmh,1,indexMPI5,j,100,
     &                       MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUGT
                        print *,iproc,': Sent xyzmh to ',j
#endif
                        CALL MPI_SEND(vxyzu,1,indexMPI4,j,101,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(f1vxyzu,1,indexMPI4,j,102,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(f2vxyzu,1,indexMPI4,j,103,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(isteps,1,indexMPI_INT,j,104,
     &                       MPI_COMM_WORLD, ierr)

                        IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) 
     &                       THEN
                           CALL MPI_SEND(icolumnsteps,1,indexMPI_INT,j,
     &                          168,MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(icolumnnext,1,indexMPI_INT,j,
     &                          169,MPI_COMM_WORLD, ierr)
                        ENDIF

                        CALL MPI_SEND(it0,1,indexMPI_INT,j,105,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(it1,1,indexMPI_INT,j,106,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(it2,1,indexMPI_INT,j,107,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(iphase,1,indexMPI_INT1,j,108,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(iscurrent,1,indexMPI_L,j,151,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(notacc,1,indexMPI_L,j,152,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(iunique,1,indexMPI_INT8,j,109,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(rho,1,indexMPI1,j,110,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(divv,1,indexMPI1,j,153,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(curlv,1,indexMPI1,j,154,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(ddv,1,indexMPI1,j,155,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(pr,1,indexMPI1,j,156,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(vsound,1,indexMPI1,j,157,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(poten,1,indexMPI1,j,111,
     &                       MPI_COMM_WORLD, ierr)

                        IF (idim_h2.EQ.idim) THEN
c                           CALL MPI_SEND(h2mol,1,indexMPIR,j,158,
c     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(h2ratio,1,indexMPIR,j,159,
     &                          MPI_COMM_WORLD, ierr)
                        ENDIF

                        IF (idustIMPL.GT.0) THEN
                           CALL MPI_SEND(dustnorm,1,indexMPIR,j,171,
     &                          MPI_COMM_WORLD, ierr)
                        ENDIF

                        IF (nlmax.EQ.1) THEN
                           CALL MPI_SEND(gradhs,1,indexMPI24,j,112,
     &                          MPI_COMM_WORLD, ierr)
                        ENDIF

                        CALL MPI_SEND(alphaMM,1,indexMPIalpha,j,113,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(f1ha,1,indexMPIalpha1,j,114,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_SEND(f2ha,1,indexMPIalpha1,j,115,
     &                       MPI_COMM_WORLD, ierr)

                        IF (imhd.EQ.idim) THEN
                           CALL MPI_SEND(Bevolxyz,1,indexMPIimhdevol,
     &                          j,116,MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(dumBevolxyz,1,indexMPIimhdevol,
     &                          j,162,MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(Bxyz,1,indexMPI3,j,117,
     &                          MPI_COMM_WORLD, ierr)
                  
                           CALL MPI_SEND(f1Bxyz,1,indexMPIimhdevol,
     &                          j,118,MPI_COMM_WORLD, ierr)
                  
                           CALL MPI_SEND(f2Bxyz,1,indexMPIimhdevol,
     &                          j,119,MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(divcurlB,1,indexMPI54,j,160,
     &                          MPI_COMM_WORLD, ierr)
#ifdef NONIDEAL
                           CALL MPI_SEND(jcurrent,1,indexMPI3,j,172,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(n_R,1,indexMPI4,j,173,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(n_electronT,1,indexMPIR,j,174,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(eta_nimhd,1,indexMPI4,j,179,
     &                          MPI_COMM_WORLD, ierr)
#endif
                        ENDIF

                        IF (encal.EQ.'r') THEN
                           CALL MPI_SEND(ekcle,1,indexMPI5,j,142,
     &                          MPI_COMM_WORLD, ierr)
                           IF (idustRT.GT.0) THEN
                              CALL MPI_SEND(dust_tk,1,indexMPI2,j,163,
     &                             MPI_COMM_WORLD, ierr)
                              CALL MPI_SEND(heatingISR,1,indexMPIH4,j,
     &                             164,MPI_COMM_WORLD, ierr)
                              IF (ioptimise_column.EQ.1) THEN
                               CALL MPI_SEND(heatingISRold,1,indexMPI24,
     &                                j,170,MPI_COMM_WORLD, ierr)
                              ENDIF
                              CALL MPI_SEND(chemistry,1,indexMPIC4,j,
     &                             165,MPI_COMM_WORLD, ierr)
                              CALL MPI_SEND(h2frac,1,indexMPI1,j,
     &                             166,MPI_COMM_WORLD, ierr)
                              CALL MPI_SEND(dh2dt,1,indexMPI24,j,
     &                             167,MPI_COMM_WORLD, ierr)
                           ENDIF
                        ENDIF

                        IF (istellar_ionisation) THEN
                           CALL MPI_SEND(HIIion,1,indexMPIR,
     &                          j,175,MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(dHIIdt,1,indexMPI2,
     &                          j,176,MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(stellarrad,1,indexMPIstellar,
     &                          j,177,MPI_COMM_WORLD, ierr)
                        ENDIF
c
c--Deallocate types
c
                        IF (istellar_ionisation) THEN
                           CALL MPI_TYPE_FREE(indexMPIstellar,ierr)
                        ENDIF
                        CALL MPI_TYPE_FREE(indexMPIalpha1,ierr)
                        CALL MPI_TYPE_FREE(indexMPIalpha,ierr)
                        CALL MPI_TYPE_FREE(indexMPI1,ierr)
                        CALL MPI_TYPE_FREE(indexMPI54,ierr)
                        CALL MPI_TYPE_FREE(indexMPIH4,ierr)
                        CALL MPI_TYPE_FREE(indexMPIC4,ierr)
                        CALL MPI_TYPE_FREE(indexMPI24,ierr)
                        CALL MPI_TYPE_FREE(indexMPI34,ierr)
                        CALL MPI_TYPE_FREE(indexMPI_INT8,ierr)
                        CALL MPI_TYPE_FREE(indexMPI_INT1,ierr)
                        CALL MPI_TYPE_FREE(indexMPI_L,ierr)
                        CALL MPI_TYPE_FREE(indexMPI_INT,ierr)
                        CALL MPI_TYPE_FREE(indexMPIR,ierr)
                        CALL MPI_TYPE_FREE(indexMPI2,ierr)
                        CALL MPI_TYPE_FREE(indexMPI3,ierr)
                        CALL MPI_TYPE_FREE(indexMPIimhdevol,ierr)
                        CALL MPI_TYPE_FREE(indexMPI4,ierr)
                        CALL MPI_TYPE_FREE(indexMPI5,ierr)
#ifdef MPIDEBUGT
                print *,iproc,': Sent transferred particles info to ',j
#endif
c
c--Transfer sink particle information
c
                        CALL MPI_SEND(npttrans,1,MPI_INTEGER,j,121, 
     &                       MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUGT
                print *,iproc,': Send TRANS ',npttrans,' SINK to ',j
#endif

                        IF (npttrans.GT.0) THEN
                           CALL MPI_TYPE_INDEXED(npttrans, 
     &                          lblocklengths, listptmass, 
     &                          MPI_REAL8, indexMPI1, ierr)
                           CALL MPI_TYPE_COMMIT(indexMPI1,ierr)
                              
                           CALL MPI_TYPE_INDEXED(npttrans, 
     &                          lblocklengths, listptmass, 
     &                          MPI_INTEGER, indexMPI_INT, ierr)
                           CALL MPI_TYPE_COMMIT(indexMPI_INT,ierr)
                 
                           CALL MPI_TYPE_INDEXED(npttrans, 
     &                          lblocklengths, listptmass, 
     &                          MPI_LOGICAL, indexMPI_L, ierr)
                           CALL MPI_TYPE_COMMIT(indexMPI_L,ierr)

                           IF (encal.EQ.'r') THEN
                              CALL MPI_TYPE_INDEXED(npttrans, 
     &                             lblocklengths, listptmass, 
     &                             iNINT, indexMPI_LIST, ierr)
                              CALL MPI_TYPE_COMMIT(indexMPI_LIST,ierr)

                              IF (istellarfeedback.GE.5) THEN
                                 CALL MPI_TYPE_INDEXED(npttrans, 
     &                                lblocklengths, listptmass, 
     &                                MPI_INTEGER2, indexMPI_INT2, ierr)
                                CALL MPI_TYPE_COMMIT(indexMPI_INT2,ierr)
                              ENDIF
                           ENDIF

                           IF (istellar_ionisation) THEN
                              CALL MPI_TYPE_INDEXED(npttrans, 
     &                             lblocklengths, listptmass, 
     &                             i3REAL8, indexMPI_S3R8, ierr)
                              CALL MPI_TYPE_COMMIT(indexMPI_S3R8,ierr)
                           ENDIF

                           CALL MPI_SEND(spinx,1,indexMPI1,j,122,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(spiny,1,indexMPI1,j,123,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(spinz,1,indexMPI1,j,124,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(angaddx,1,indexMPI1,j,125,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(angaddy,1,indexMPI1,j,126,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(angaddz,1,indexMPI1,j,127,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(spinadx,1,indexMPI1,j,128,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(spinady,1,indexMPI1,j,129,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(spinadz,1,indexMPI1,j,130,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(ptmsyn,1,indexMPI1,j,131,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(ptmadd,1,indexMPI1,j,132,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(xmomsyn,1,indexMPI1,j,133,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(ymomsyn,1,indexMPI1,j,134,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(zmomsyn,1,indexMPI1,j,135,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(xmomadd,1,indexMPI1,j,136,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(ymomadd,1,indexMPI1,j,137,
     &                          MPI_COMM_WORLD, ierr)
                           CALL MPI_SEND(zmomadd,1,indexMPI1,j,138,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(ptmassinner,1,indexMPI1,j,139,
     &                          MPI_COMM_WORLD, ierr)

                           CALL MPI_SEND(nactotal,1,indexMPI_INT,j,140,
     &                          MPI_COMM_WORLD, ierr)

                           IF (encal.EQ.'r') THEN
                              CALL MPI_SEND(nexposed,1,indexMPI_INT,
     &                             j, 143, MPI_COMM_WORLD, ierr)

                              CALL MPI_SEND(nexposedold,1,indexMPI_INT,
     &                             j, 144, MPI_COMM_WORLD, ierr)

                              CALL MPI_SEND(ptmassluminosity,1,
     &                             indexMPI1,j,145,MPI_COMM_WORLD,ierr)

                              CALL MPI_SEND(cummulativeenergy,1,
     &                             indexMPI1,j,146,MPI_COMM_WORLD,ierr)

                              CALL MPI_SEND(actualcumenergy,1,
     &                             indexMPI1,j,147,MPI_COMM_WORLD,ierr)

                              CALL MPI_SEND(iexposedold,1,
     &                             indexMPI_L,j,148,MPI_COMM_WORLD,ierr)

                              CALL MPI_SEND(listexposed,1,indexMPI_LIST,
     &                             j,149,MPI_COMM_WORLD,ierr)

                           CALL MPI_SEND(listexposedold,1,indexMPI_LIST,
     &                             j,150,MPI_COMM_WORLD,ierr)

                              IF (istellarfeedback.GE.5) THEN
                                 CALL MPI_SEND(ptmassform,1,indexMPI1,
     &                                j,180,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(diskinit,1,indexMPI1,
     &                                j,181,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(stellartime,1,indexMPI1,
     &                                j,182,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(stellarmass,1,indexMPI1,
     &                                j,183,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(stellarradius,1,
     &                            indexMPI1,j,184,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(stellarluminosity,1,
     &                            indexMPI1,j,185,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(lyontime,1,indexMPI1,
     &                                j,186,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(lyonmass,1,indexMPI1,
     &                                j,187,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(lyonradius,1,indexMPI1,
     &                                j,188,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(lyonluminosity,1,
     &                            indexMPI1,j,189,MPI_COMM_WORLD, ierr)

                                 CALL MPI_SEND(listlyon,1,indexMPI_INT2,
     &                                j,190,MPI_COMM_WORLD, ierr)

                              ENDIF
                           ENDIF

                           IF (istellar_ionisation) THEN
                              CALL MPI_SEND(stellar_radtempion,1,
     &                             indexMPI_S3R8,
     &                             j,178,MPI_COMM_WORLD,ierr)
                           ENDIF
c
c--Send information about which particle matches each sink so that listpm
c     can be re-constructed on the new MPI process
c
                           CALL MPI_SEND(listpmtrans,npttrans,
     &                          MPI_INTEGER, j, 141, 
     &                          MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUGT
                ipart = lsendlist(listpmtrans(1)) + 1
                print *,iproc,': TRANS SINK INFO ',listpmtrans(1),
     &               ipart,listpm(listptmass(1)+1),iphase(ipart),
     &               xyzmh(1,ipart),xyzmh(2,ipart),xyzmh(3,ipart),
     &               xyzmh(4,ipart),xyzmh(5,ipart),
     &               ptmsyn(listptmass(1)+1),xmomsyn(listptmass(1)+1),
     &               xmomadd(listptmass(1)+1)
#endif
                           IF (istellar_ionisation) THEN
                              CALL MPI_TYPE_FREE(indexMPI_S3R8,ierr)
                           ENDIF
                           IF (encal.EQ.'r') THEN
                              CALL MPI_TYPE_FREE(indexMPI_LIST,ierr)
                              IF (istellarfeedback.GE.5) THEN
                                 CALL MPI_TYPE_FREE(indexMPI_INT2,ierr)
                              ENDIF
                           ENDIF
                           CALL MPI_TYPE_FREE(indexMPI_L,ierr)
                           CALL MPI_TYPE_FREE(indexMPI_INT,ierr)
                           CALL MPI_TYPE_FREE(indexMPI1,ierr)
                        ENDIF
#ifdef MPIDEBUGT
                        print *,iproc,': Finished transfer to ',j
#endif
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
#ifdef MPIDEBUGT
        print *,iproc,': TEST ',inumgot,inumtoget,numrecv(ii+1)
#endif
            IF (inumgot.LT.inumtoget .AND. numrecv(ii+1).GT.0) THEN
               inumgot = inumgot + 1

#ifdef MPIDEBUGT
        print *,iproc,': Receiving transfer of ',numrecv(ii+1),
     &              ' from ',ii
#endif

               CALL MPI_RECV(xyzmh(1,npartnew+1), idim, i5REAL8, 
     &              ii, 100, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ii.NE.iprocrec) THEN
                  WRITE (*,*) 'ERROR - ii.NE.iprocrec ',ii,iprocrec
                  CALL quit(1)
               ENDIF
               igotnumber = inumber
               IF (inumber.NE.numrecv(ii+1)) THEN
                  WRITE (*,*) 'ERROR - inumber.NE.numrecv(ii+1) ',
     &                 inumber,numrecv(ii+1),ii,iprocrec
                  CALL quit(1)
               ENDIF

#ifdef MPIDEBUGT
               print *,iproc,': Received transfer of ',igotnumber,
     &              ' from ',iprocrec
#endif
               CALL MPI_RECV(vxyzu(1,npartnew+1), idim, i4REAL8,
     &              iprocrec, 101, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(f1vxyzu(1,npartnew+1), idim, i4REAL8,
     &              iprocrec, 102, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(f2vxyzu(1,npartnew+1), idim, i4REAL8,
     &              iprocrec, 103, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(isteps(npartnew+1), idim, MPI_INTEGER,
     &              iprocrec, 104, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) THEN
                  CALL MPI_RECV(icolumnsteps(npartnew+1), idim, 
     &                 MPI_INTEGER,
     &                 iprocrec, 168, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(icolumnnext(npartnew+1), idim, 
     &                 MPI_INTEGER,
     &                 iprocrec, 169, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF
               ENDIF

               CALL MPI_RECV(it0(npartnew+1), idim, MPI_INTEGER,
     &              iprocrec, 105, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(it1(npartnew+1), idim, MPI_INTEGER,
     &              iprocrec, 106, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(it2(npartnew+1), idim, MPI_INTEGER,
     &              iprocrec, 107, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(iphase(npartnew+1), idim, MPI_INTEGER1,
     &              iprocrec, 108, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_INTEGER1,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF
#ifdef MPIDEBUGT
c               print *,iproc,': iPHASE ',(iphase(ipart),
c     &              ipart= npartnew + 1, npartnew + igotnumber),
c     &              ' from ',iprocrec
#endif               

               CALL MPI_RECV(iscurrent(npartnew+1), idim, MPI_LOGICAL,
     &              iprocrec, 151, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_LOGICAL,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber iscurrent'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(notacc(npartnew+1), idim, MPI_LOGICAL,
     &              iprocrec, 152, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,MPI_LOGICAL,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber notacc'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(iunique(npartnew+1), idim, MPI_INTEGER8,
     &              iprocrec, 109, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER8,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

c                  DO iii = npartnew+1,npartnew+MIN(5,igotnumber)
c                     print *,'REC ',iproc,iii,iunique(iii)
c                  END DO

               CALL MPI_RECV(rho(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 110, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               DO ipart = npartnew + 1, npartnew + igotnumber
                  dumrho(ipart) = rho(ipart)
               END DO

               CALL MPI_RECV(divv(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 153, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber divv'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(curlv(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 154, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber curlv'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(ddv(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 155, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber ddv'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(pr(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 156, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber pr'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(vsound(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 157, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber vsound'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(poten(npartnew+1), idim, MPI_REAL4,
     &              iprocrec, 111, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL4, inumber, ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               IF (idim_h2.EQ.idim) THEN
c                  CALL MPI_RECV(h2mol(npartnew+1), idim, MPI_REAL8,
c     &                 iprocrec, 158, MPI_COMM_WORLD, istatus, ierr)
c                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
c                  IF (inumber.NE.igotnumber) THEN
c                     WRITE (*,*) 'ERROR - balance igotnumber h2mol'
c                     CALL quit(1)
c                  ENDIF

                  CALL MPI_RECV(h2ratio(npartnew+1), idim, MPI_REAL8,
     &                 iprocrec, 159, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber h2ratio'
                     CALL quit(1)
                  ENDIF
               ENDIF

               IF (idustIMPL.GT.0) THEN
                  CALL MPI_RECV(dustnorm(npartnew+1), idim, MPI_REAL8,
     &                 iprocrec, 171, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber dustnorm'
                     CALL quit(1)
                  ENDIF
               ENDIF

               IF (nlmax.EQ.1) THEN
                  CALL MPI_RECV(gradhs(1,npartnew+1), idim, i2REAL4,
     &                 iprocrec, 112, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i2REAL4, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF
               ENDIF

               CALL MPI_RECV(alphaMM(1,npartnew+1),idim,ialphaREAL4,
     &              iprocrec, 113, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, ialphaREAL4,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(f1ha(1,npartnew+1), idim, ialphaREAL41,
     &              iprocrec, 114, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus,ialphaREAL41,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               CALL MPI_RECV(f2ha(1,npartnew+1), idim, ialphaREAL41,
     &              iprocrec, 115, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, ialphaREAL41,inumber,ierr)
               IF (inumber.NE.igotnumber) THEN
                  WRITE (*,*) 'ERROR - balance igotnumber'
                  CALL quit(1)
               ENDIF

               IF (imhd.EQ.idim) THEN
                  CALL MPI_RECV(Bevolxyz(1,npartnew+1), idim, 
     &                 imhdevolREAL8,
     &                 iprocrec, 116, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, imhdevolREAL8, 
     &                 inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(dumBevolxyz(1,npartnew+1), idim, 
     &                 imhdevolREAL8,
     &                 iprocrec, 162, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, imhdevolREAL8, 
     &                 inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(Bxyz(1,npartnew+1), idim, i3REAL8,
     &                 iprocrec, 117, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i3REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(f1Bxyz(1,npartnew+1), idim, 
     &                 imhdevolREAL8,
     &                 iprocrec, 118, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, imhdevolREAL8, 
     &                 inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF
                  
                  CALL MPI_RECV(f2Bxyz(1,npartnew+1), idim, 
     &                 imhdevolREAL8,
     &                 iprocrec, 119, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, imhdevolREAL8, 
     &                 inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(divcurlB(1,npartnew+1), idim, i5REAL4,
     &                 iprocrec, 160, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i5REAL4, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber'
                     CALL quit(1)
                  ENDIF
#ifdef NONIDEAL
                  CALL MPI_RECV(jcurrent(1,npartnew+1), idim, i3REAL8,
     &                 iprocrec, 172, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i3REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber: jcurrent'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(n_R(1,npartnew+1), idim, i4REAL8,
     &                 iprocrec, 173, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber: n_R'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(n_electronT(npartnew+1), idim,MPI_REAL8,
     &                 iprocrec, 174, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber: n_elecT'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(eta_nimhd(1,npartnew+1), idim, i4REAL8,
     &                 iprocrec, 179, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i4REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber: eta_nimhd'
                     CALL quit(1)
                  ENDIF
#endif
               ENDIF

               IF (encal.EQ.'r') THEN
                  CALL MPI_RECV(ekcle(1,npartnew+1), idim, i5REAL8,
     &                 iprocrec, 142, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, i5REAL8, inumber, ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber ekcle'
                     CALL quit(1)
                  ENDIF

                  IF (idustRT.GT.0) THEN
                     CALL MPI_RECV(dust_tk(1,npartnew+1), idim, i2REAL8,
     &                    iprocrec, 163, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, i2REAL8, inumber, ierr)
                     IF (inumber.NE.igotnumber) THEN
                        WRITE (*,*) 'ERROR - balance igotnumber dust'
                        CALL quit
                     ENDIF

                     CALL MPI_RECV(heatingISR(1,npartnew+1), idim,
     &                    iHREAL4, iprocrec, 164, MPI_COMM_WORLD, 
     &                    istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, iHREAL4, inumber, ierr)
                     IF (inumber.NE.igotnumber) THEN
                        WRITE (*,*) 'ERROR - balance igotnumber heating'
                        CALL quit
                     ENDIF

                     IF (ioptimise_column.EQ.1) THEN
                        CALL MPI_RECV(heatingISRold(1,npartnew+1), idim,
     &                       i2REAL4, iprocrec, 170, MPI_COMM_WORLD, 
     &                       istatus, ierr)
                        CALL MPI_GET_COUNT(istatus,i2REAL4,inumber,ierr)
                        IF (inumber.NE.igotnumber) THEN
                        WRITE (*,*) 'ERROR - balance igotnumber heating'
                           CALL quit
                        ENDIF
                     ENDIF

                     CALL MPI_RECV(chemistry(1,npartnew+1), idim,
     &                    iCREAL4, iprocrec, 165, MPI_COMM_WORLD, 
     &                    istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, iCREAL4, inumber, ierr)
                     IF (inumber.NE.igotnumber) THEN
                        WRITE (*,*) 'ERROR - balance igotnumber chem'
                        CALL quit
                     ENDIF

                     CALL MPI_RECV(h2frac(npartnew+1), idim, MPI_REAL4,
     &                    iprocrec, 166, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, MPI_REAL4, 
     &                    inumber, ierr)
                     IF (inumber.NE.igotnumber) THEN
                        WRITE (*,*) 'ERROR - balance igotnumber h2frac'
                        CALL quit
                     ENDIF

                     CALL MPI_RECV(dh2dt(1,npartnew+1), idim, i2REAL4,
     &                    iprocrec, 167, MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus, i2REAL4, 
     &                    inumber, ierr)
                     IF (inumber.NE.igotnumber) THEN
                        WRITE (*,*) 'ERROR - balance igotnumber dh2dt'
                        CALL quit
                     ENDIF

                  ENDIF

               ENDIF

               IF (istellar_ionisation) THEN
                  CALL MPI_RECV(HIIion(npartnew+1), idim, MPI_REAL8,
     &                 iprocrec, 175, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber HIIion'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(dHIIdt(1,npartnew+1), idim, i2REAL8,
     &                 iprocrec, 176, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,i2REAL8,inumber,ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber dHIIdt'
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(stellarrad(1,1,npartnew+1), idim,
     &                 istellarREAL8,
     &                 iprocrec, 177, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,istellarREAL8,inumber,ierr)
                  IF (inumber.NE.igotnumber) THEN
                     WRITE (*,*) 'ERROR - balance igotnumber stellarrad'
                     CALL quit(1)
                  ENDIF
               ENDIF

#ifdef MPIDEBUGT
               print *,iproc,': Received transfer of ',inumber,' from ',
     &              iprocrec
#endif
c
c--Transfer sink particle information
c
               CALL MPI_RECV(nptrecv, 1, MPI_INTEGER,
     &              iprocrec, 121, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber, ierr)
               IF (inumber.NE.1) THEN
                  WRITE (*,*) 'ERROR - received no sink info ',iproc,
     &                 iprocrec
                  CALL quit(1)
               ENDIF
               IF (nptmassnew + nptrecv.GT.iptdim) THEN
                  WRITE (*,*) 'ERROR - nptmassnew+nptrecv.GT.iptdim',
     &                 nptmass, nptmassnew, nptrecv, iprocrec
               ENDIF
                  
#ifdef MPIDEBUGT
                print *,iproc,': Recv TRANS ',nptrecv,' SINK from ',
     &              iprocrec
#endif

               IF (nptrecv.GT.0) THEN
                  CALL MPI_RECV(spinx(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 122, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv1 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(spiny(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 123, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv2 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(spinz(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 124, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv3 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF
                  
                  CALL MPI_RECV(angaddx(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 125, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv4 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(angaddy(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 126, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv5 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(angaddz(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 127, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv6 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(spinadx(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 128, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv7 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(spinady(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 129, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv8 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(spinadz(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 130, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecv9 ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(ptmsyn(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 131, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvA ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(ptmadd(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 132, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvB ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(xmomsyn(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 133, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvC ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(ymomsyn(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 134, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvD ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(zmomsyn(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 135, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvE ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(xmomadd(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 136, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvF ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(ymomadd(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 137, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvG ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(zmomadd(nptmassnew+1),nptrecv,MPI_REAL8,
     &                 iprocrec, 138, MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvH ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF


                  CALL MPI_RECV(ptmassinner(nptmassnew+1),nptrecv,
     &                 MPI_REAL8, iprocrec, 139, MPI_COMM_WORLD, 
     &                 istatus, ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvI ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  CALL MPI_RECV(nactotal(nptmassnew+1),nptrecv,
     &                 MPI_INTEGER, iprocrec, 140, MPI_COMM_WORLD,
     &                 istatus,ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvJ ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF

                  IF (encal.EQ.'r') THEN
                     CALL MPI_RECV(nexposed(nptmassnew+1),nptrecv,
     &                    MPI_INTEGER, iprocrec, 143, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,MPI_INTEGER,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvK ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                     CALL MPI_RECV(nexposedold(nptmassnew+1),nptrecv,
     &                    MPI_INTEGER, iprocrec, 144, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,MPI_INTEGER,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvL ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                   CALL MPI_RECV(ptmassluminosity(nptmassnew+1),nptrecv,
     &                    MPI_REAL8, iprocrec, 145, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvM ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                  CALL MPI_RECV(cummulativeenergy(nptmassnew+1),nptrecv,
     &                    MPI_REAL8, iprocrec, 146, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvN ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                    CALL MPI_RECV(actualcumenergy(nptmassnew+1),nptrecv,
     &                    MPI_REAL8, iprocrec, 147, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvO ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                     CALL MPI_RECV(iexposedold(nptmassnew+1),nptrecv,
     &                    MPI_LOGICAL, iprocrec, 148, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,MPI_LOGICAL,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvP ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                     CALL MPI_RECV(listexposed(1,nptmassnew+1),nptrecv,
     &                    iNINT, iprocrec, 149, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,iNINT,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvQ ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                   CALL MPI_RECV(listexposedold(1,nptmassnew+1),nptrecv,
     &                    iNINT, iprocrec, 150, MPI_COMM_WORLD,
     &                    istatus,ierr)
                    CALL MPI_GET_COUNT(istatus,iNINT,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvR ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF

                     IF (istellarfeedback.GE.5) THEN

                        CALL MPI_RECV(ptmassform(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 180, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvS ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(diskinit(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 181, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvT ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(stellartime(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 182, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvU ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(stellarmass(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 183, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvV ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(stellarradius(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 184, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvW ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(stellarluminosity(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 185, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvX ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(lyontime(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 186, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvY ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(lyonmass(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 187, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvY ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(lyonradius(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 188, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvZ ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(lyonluminosity(nptmassnew+1),
     &                       nptrecv, MPI_REAL8, iprocrec, 189, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_REAL8,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvAA ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                        CALL MPI_RECV(listlyon(nptmassnew+1),
     &                       nptrecv, MPI_INTEGER2, iprocrec, 190, 
     &                       MPI_COMM_WORLD, istatus, ierr)
                  CALL MPI_GET_COUNT(istatus,MPI_INTEGER2,inumber,ierr)
                        IF (inumber.NE.nptrecv) THEN
                           WRITE (*,*) 'ERROR - inumber.NE.nptrecvAB ',
     &                       inumber,nptrecv
                           CALL quit(1)
                        ENDIF

                     ENDIF ! istellarfeedback
                  ENDIF ! encal

                  IF (istellar_ionisation) THEN
                     CALL MPI_RECV(stellar_radtempion(1,nptmassnew+1),
     &                    nptrecv,
     &                    i3REAL8, iprocrec, 178, MPI_COMM_WORLD,
     &                    istatus,ierr)
                     CALL MPI_GET_COUNT(istatus,i3REAL8,inumber,ierr)
                     IF (inumber.NE.nptrecv) THEN
                        WRITE (*,*) 'ERROR - inumber.NE.nptrecvSR ',
     &                       inumber,nptrecv
                        CALL quit(1)
                     ENDIF                     
                  ENDIF

                  CALL MPI_RECV(listpm(nptmassnew+1),nptrecv,
     &                 MPI_INTEGER, iprocrec, 141, MPI_COMM_WORLD,
     &                 istatus,ierr)
                  CALL MPI_GET_COUNT(istatus, MPI_INTEGER, inumber,ierr)
                  IF (inumber.NE.nptrecv) THEN
                     WRITE (*,*) 'ERROR - inumber.NE.nptrecvK ',inumber,
     &                    nptrecv
                     CALL quit(1)
                  ENDIF
c
c--Set indexing for sink particles
c
                  DO k = nptmassnew + 1, nptmassnew + nptrecv
                     listpm(k) = listpm(k) + npartnew
                     listrealpm(listpm(k)) = k
                  END DO

#ifdef MPIDEBUGT
                  ipart = listpm(nptmassnew + 1)
                  print *,iproc,': RECV SINK INFO ',ipart,
     &                 xyzmh(1,ipart),xyzmh(2,ipart),xyzmh(3,ipart),
     &                 xyzmh(4,ipart),xyzmh(5,ipart),
     &                 ptmsyn(nptmassnew + 1),xmomsyn(nptmassnew + 1),
     &                 xmomadd(nptmassnew + 1)
#endif

               ENDIF
               npartnew = npartnew + igotnumber
               nactive = nactive + igotnumber
               nptmassnew = nptmassnew + nptrecv
            ENDIF
         ENDIF
      END DO

      print *,iproc,': New value of npart, nptmass = ',npartnew, npart, 
     &     nptmassnew, nptmass

      IF (istellar_ionisation) THEN
         CALL MPI_TYPE_FREE(istellarREAL8,ierr)
      ENDIF

      IF (encal.EQ.'r') THEN
         CALL MPI_TYPE_FREE(iNINT,ierr)
      ENDIF

      CALL MPI_TYPE_FREE(ialphaREAL41,ierr)
      CALL MPI_TYPE_FREE(ialphaREAL4,ierr)
      CALL MPI_TYPE_FREE(i5REAL4,ierr)
      CALL MPI_TYPE_FREE(iCREAL4,ierr)
      CALL MPI_TYPE_FREE(iHREAL4,ierr)
      CALL MPI_TYPE_FREE(i2REAL4,ierr)
      CALL MPI_TYPE_FREE(i3REAL4,ierr)
      CALL MPI_TYPE_FREE(i2REAL8,ierr)
      CALL MPI_TYPE_FREE(i3REAL8,ierr)
      CALL MPI_TYPE_FREE(imhdevolREAL8,ierr)
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
      CALL MPI_TYPE_FREE(i5REAL8,ierr)
c
c--Need to set particle indexxing
c
      DO i = npart+1, npartnew
         iscurrent(i) = .FALSE.
         isort(i) = i
         iorig(i) = i
         DO k = 1, 5
            dumxyzmh(k,i) = xyzmh(k,i)
         END DO
         DO k = 1, 4
            dumvxyzu(k,i) = vxyzu(k,i)
         END DO
         DO k = 1, isizealphaMM
            dumalpha(k,i) = alphaMM(k,i)
         END DO
         IF (encal.EQ.'r') THEN
            DO k = 1, 5
               dumekcle(k,i) = ekcle(k,i)
            END DO
         ENDIF
         IF (imhd.EQ.idim) THEN
            DO k = 1, 3
               dumBevolxyz(k,i) = Bevolxyz(k,i)
            END DO
         ENDIF
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
      npttrans = 0
      IF (numtransfer.GT.0) THEN
         DO ipart = 1, npart
            IF (listtrans(ipart).GE.0 .AND. 
     &           listtrans(ipart).NE.iproc) THEN
               IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
                  npttrans = npttrans + 1
                  listpm(listrealpm(ipart)) = 0
                  listrealpm(ipart) = 0
               ENDIF
               iphase(ipart) = -1
               nactive = nactive - 1
            ENDIF
         END DO
      ENDIF
      IF (npttrans.GT.0) THEN
         nptmass = 0
         DO ipt = 1, nptmassnew
            IF (listpm(ipt).NE.0) THEN
               nptmass = nptmass + 1
               listpm(nptmass) = listpm(ipt)
               listrealpm(listpm(nptmass)) = nptmass

               spinx(nptmass) = spinx(ipt)
               spiny(nptmass) = spiny(ipt)
               spinz(nptmass) = spinz(ipt)

               angaddx(nptmass) = angaddx(ipt)
               angaddy(nptmass) = angaddy(ipt)
               angaddz(nptmass) = angaddz(ipt)

               spinadx(nptmass) = spinadx(ipt)
               spinady(nptmass) = spinady(ipt)
               spinadz(nptmass) = spinadz(ipt)

               ptmassinner(nptmass) = ptmassinner(ipt)
               nactotal(nptmass) = nactotal(ipt)

               ptmsyn(nptmass) = ptmsyn(ipt)
               ptmadd(nptmass) = ptmadd(ipt)

               xmomsyn(nptmass) = xmomsyn(ipt)
               ymomsyn(nptmass) = ymomsyn(ipt)
               zmomsyn(nptmass) = zmomsyn(ipt)

               xmomadd(nptmass) = xmomadd(ipt)
               ymomadd(nptmass) = ymomadd(ipt)
               zmomadd(nptmass) = zmomadd(ipt)

               IF (istellarfeedback.GE.5) THEN
                  ptmassform(nptmass) = ptmassform(ipt)
                  diskinit(nptmass) = diskinit(ipt)
                  stellartime(nptmass) = stellartime(ipt)
                  stellarmass(nptmass) = stellarmass(ipt)
                  stellarradius(nptmass) = stellarradius(ipt)
                  stellarluminosity(nptmas) = stellarluminosity(ipt)
                  lyontime(nptmass) = lyontime(ipt)
                  lyonmass(nptmass) = lyonmass(ipt)
                  lyonradius(nptmass) = lyonradius(ipt)
                  lyonluminosity(nptmass) = lyonluminosity(ipt)
                  listlyon(nptmass) = listlyon(ipt)
               ENDIF
            ENDIF
         END DO
      ELSE
         nptmass = nptmassnew
      ENDIF

      npartold = npart
      npart = npartnew
c      n1 = npart
      nptmasslast = nptmass

      CALL MPI_ALLREDUCE(nptmass, nptmasstot, 1, MPI_INTEGER, MPI_SUM,
     &     MPI_COMM_WORLD,ierr)

#ifdef MPIDEBUGT
      print *,iproc,': New npart and nactive ',npart,nactive,nptmass,
     &     nptmasstot
      print *,iproc,': PTMASS ID ',(listpm(i),i=1,nptmass),
     &     (iphase(listpm(i)),i=1,nptmass)
#endif
c
c--Now move new particles around within MPI process to fill in blanks left
c     by killed particles or by particles moved to other MPI processes.
c
      print *,iproc,': MOVE ',npart-npartold,nlistinactive,npart
      IF (iproc.EQ.0) WRITE (iprint,*) 'MOVE ',npart-npartold,
     &     nlistinactive,npart
      nlistinactiveold = nlistinactive

      DO i = npart, npartold + 1, -1
         IF (nlistinactive.GT.0) THEN
            CALL move_particle(i,listinactive(nlistinactive))
            nlistinactive = nlistinactive - 1
         ENDIF
      END DO
      npart = npart - (nlistinactiveold - nlistinactive)
c      n1 = n1 - (nlistinactiveold - nlistinactive)
      IF (iproc.EQ.0) WRITE (iprint,*) 'MOVE-POST ',npart,
     &     nlistinactive,nlistinactiveold
c
c--Test MPI decomposition: writes out particles assigned to each process
c
      IF (.FALSE.) THEN
         IF (irec.LT.10) THEN
            WRITE (idumpnumber,99200) irec
99200       FORMAT('00',I1)
         ELSEIF (irec.LT.100) THEN
            WRITE (idumpnumber,99201) irec
99201       FORMAT('0',I2)
         ELSEIF (irec.LT.1000) THEN
            WRITE (idumpnumber,99202) irec
99202       FORMAT(I3)
         ENDIF

         IF (iproc.LT.10) THEN
            WRITE (iprocnumber,99300) iproc
99300       FORMAT('0',I1)
         ELSEIF (iproc.LT.100) THEN
            WRITE (iprocnumber,99301) iproc
99301       FORMAT(I2)
         ENDIF

         filename_split = 'TESTMPI'//idumpnumber//iprocnumber
         OPEN (60+iproc,FILE=filename_split)
         DO i = 1, npart
            IF (iphase(i).NE.-1) THEN
               WRITE (60+iproc,777) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)
 777           FORMAT(3(1PE12.5,1X))
            ENDIF
         
         END DO

         CLOSE (60+iproc)
      ENDIF
c
c--End of test of MPI decomposition
c
c
c--Need to reset lists of particles in each timestep bin
c
      DO i = 1, nbinmax
         nlstbins(i) = 0
         it1bin(i) = 2**(i-1)
         it2bin(i) = 2**i
      END DO

      xlog2 = 0.30103

      DO i = 1, npart
         IF (iphase(i).NE.-1) THEN
            it0(i) = 0
            it1(i) = isteps(i)/2
            it2(i) = isteps(i)

            IF (idustRT.GT.0  .AND. ioptimise_column.EQ.1 
     &           .AND. iphase(i).EQ.0) THEN
               icolumnsteps(i) = MAX(icolumnsteps(i),isteps(i))
               icolumnnext(i) = icolumnsteps(i)
            ENDIF

            ibin = INT(LOG10(REAL(isteps(i)))/xlog2+0.5)
            IF (ibin.GT.nbinmax .OR. ibin.GT.30 .OR. ibin.LT.0) THEN
               WRITE (*,*) 'ERROR - ibin.GT.nbinmax B ',ibin
               WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax B ',ibin
               CALL quit(1)
            ENDIF
            nlstbins(ibin) = nlstbins(ibin) + 1
            listbins(nlstbins(ibin),ibin) = i
            IF (it1bin(ibin).NE.it1(i)) THEN
               WRITE (*,*) 'ERROR - it1bin B'
               CALL quit(1)
            ENDIF
            IF (it2bin(ibin).NE.it2(i)) THEN
               WRITE (*,*) 'ERROR - it2bin B'
               CALL quit(1)
            ENDIF
         ENDIF
      END DO
c
c--Need to re-make ghosts and re-make tree
c
      IF (ibound.GT.0) THEN
         CALL ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)

         DO i = npart+1, ntot
            DO k = 1, 5
               dumxyzmh(k,i) = xyzmh(k,i)
            END DO
            DO k = 1, 4
               dumvxyzu(k,i) = vxyzu(k,i)
            END DO
            DO k = 1, isizealphaMM
               dumalpha(k,i) = alphaMM(k,i)
            END DO
            IF (encal.EQ.'r') THEN
               DO k = 1, 5
                  dumekcle(k,i) = ekcle(k,i)
               END DO
            ENDIF
            IF (imhd.EQ.idim) THEN
               DO k = 1, 3
                  dumBevolxyz(k,i) = Bevolxyz(k,i)
               END DO
            ENDIF
         END DO
      ELSE
         nghost = 0
         ntot = npart + nghost
      ENDIF
#ifdef MPIDEBUGT
      print *,iproc,': Balance is making tree ',ntot, npart
#endif
      CALL insulate(1, 0, ntot, npart, dumxyzmh, f1vxyzu)
c
c--End of MPI-only section
c
#endif

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

      LOGICAL ifirst
      DATA ifirst/.TRUE./

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
      INCLUDE 'COMMONS/mpisink'
      INCLUDE 'COMMONS/mpidebug'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/astrcon'

      DIMENSION listtrans(idim), numlist(nummaxproc), xyzmh(5,idim)
      INTEGER classify
      INTEGER*1 iphase(idim2)
      REAL*4 rhomaxlocal,rhomaxtot

c
c--Identify which MPI process each particle should belong to
c
      DO i = 1, numproc
         numlist(i) = 0
      END DO
c
c--For balance type 'r' need to find maximum density (rho) particle across all
c     MPI processes because this is where the domain decomposition is centred.
c
      IF (mpitype.EQ.'r') THEN
         rhomaxlocal = 0.
         irhomax = 0
         DO ipart = 1, npart
            IF (iphase(ipart).EQ.0 .AND. rho(ipart).GT.rhomaxlocal) THEN
               irhomax = ipart
               rhomaxlocal = rho(ipart)
            ENDIF
         END DO
         rhomaxlocal = rho(irhomax)
         print *,iproc,': BAL found ',irhomax,rho(irhomax)
         print *,iproc,': BAL loc ',xyzmh(1,irhomax),xyzmh(2,irhomax),
     &        xyzmh(3,irhomax)
         CALL MPI_ALLREDUCE(rhomaxlocal,rhomaxtot,1,MPI_REAL4,MPI_MAX,
     &        MPI_COMM_WORLD, ierr)

         IF (rhomaxtot.GT.densmpi) THEN
            IF (rhomaxtot.EQ.rhomaxlocal) THEN
               inum = iproc + 1
            ELSE
               inum = 0
            ENDIF
            print *,iproc,': BAL inum ',inum,rhomaxtot
            CALL MPI_ALLREDUCE(inum,inumkeep,1,MPI_INTEGER,MPI_MAX,
     &           MPI_COMM_WORLD, ierr)
            inumkeep = inumkeep - 1
            print *,iproc,': BAL inumkeep ',inumkeep
            IF (inumkeep.LT.0 .OR. inumkeep.GE.numproc) THEN
               WRITE (*,*) 'ERROR - Balance: Invalid inumkeep ',inumkeep
               CALL quit(1)
            ENDIF

            IF (inumkeep.EQ.iproc) THEN
               xlocationmpi = xyzmh(1,irhomax)
               ylocationmpi = xyzmh(2,irhomax)
               zlocationmpi = xyzmh(3,irhomax)
            ENDIF
            CALL MPI_BCAST(xlocationmpi,1,MPI_REAL8,inumkeep,
     &           MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(ylocationmpi,1,MPI_REAL8,inumkeep,
     &           MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(zlocationmpi,1,MPI_REAL8,inumkeep,
     &           MPI_COMM_WORLD, ierr)
            print *,iproc,': BALANCE location ',xlocationmpi,
     &           ylocationmpi,zlocationmpi
         ELSE
            print *,iproc,': BALANCE: location unchanged'
         ENDIF
      ENDIF
c
c--For balance type '1', if there is only one sink particle, then centre
c     the MPI domain decomposition on it.  If there is more than one,
c     then leave the centre unchanged from last time.
c
      IF (mpitype.EQ.'1') THEN
         IF (nptmasstot.EQ.1) THEN
            IF (nptmass.EQ.1) THEN
               iprocsinklocal = iproc
            ELSE
               iprocsinklocal = 0
            ENDIF
            CALL MPI_ALLREDUCE(iprocsinklocal,iprocsink,1,MPI_INTEGER,
     &           MPI_MAX,MPI_COMM_WORLD, ierr)

            IF (iprocsink.EQ.iproc) THEN
               xlocationmpi = xyzmh(1,listpm(1))
               ylocationmpi = xyzmh(2,listpm(1))
               zlocationmpi = xyzmh(3,listpm(1))
            ENDIF
            CALL MPI_BCAST(xlocationmpi,1,MPI_REAL8,iprocsink,
     &           MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(ylocationmpi,1,MPI_REAL8,iprocsink,
     &           MPI_COMM_WORLD, ierr)
            CALL MPI_BCAST(zlocationmpi,1,MPI_REAL8,iprocsink,
     &           MPI_COMM_WORLD, ierr)
            print *,iproc,': BALANCE location for mpitype.EQ.1 ',
     &           xlocationmpi,ylocationmpi,zlocationmpi
         ELSE
            print *,iproc,': BALANCE: location unchanged'
         ENDIF
      ENDIF
c
c--For Cartesian geometry, process 0 finds out how many particles are on each
c     node so that it can figure out where the new boundaries should be.  It 
c     then broadcasts the new node boundaries to the other MPI processes.
c
      IF (mpitype.EQ.'c') THEN
         print *,iproc,' BAL entered mpitype.EQ.c'

         xtemp = LOG(REAL(numproc))/LOG(2.0)
         IF (xtemp - INT(xtemp+1.0E-4).GT.1.0E-4) THEN
            WRITE (*,*) 'ERROR - numproc is not a power of two ',
     &           numproc,xtemp
            CALL quit(1)
         ENDIF
c
c--Find total number of active particles on all MPI processes
c
         CALL MPI_ALLREDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,ierr)

         print *,iproc,' BAL nactivetot ',nactivetot
c
c--Determine number of divisions in x, y, z
c
         nsplit(1) = (INT(((REAL(numproc)-0.1)**(1.0/3.0))/2.0)+1)*2
         nsplit(2) = (INT(((REAL(numproc/nsplit(1))-0.1)**(1.0/2.0))/
     &        2.0)+1)*2
         nsplit(3) = numproc/nsplit(1)/nsplit(2)

         print *,iproc,' BAL nsplit ',(nsplit(i),i=1,3)

         DO idimension = 1, 3
            nsplittemp = nsplit(idimension)
            CALL MPI_BCAST(nsplittemp,1,MPI_INTEGER,0,
     &           MPI_COMM_WORLD,ierr)
            nsplit(idimension) = nsplittemp
         END DO
c
c--Need to find extent of particle distribution if not done before or if
c     numproc is less than 64 (since then x,y, or z will have only 2 domains
c     and the code below won't know how to change the split coords).
c
         DO idimension = 1, 3
            xyzmax(idimension) = -1.0E-30
            xyzmin(idimension) = 1.0E-30
            DO i = 1, npart
               xyzmax(idimension) = MAX(xyzmax(idimension),
     &              xyzmh(idimension,i))
               xyzmin(idimension) = MIN(xyzmin(idimension),
     &              xyzmh(idimension,i))
            END DO
         END DO
c
c--Find global bounds of active particles
c
         DO idimension = 1, 3
            xsend = xyzmax(idimension)
            CALL MPI_ALLREDUCE(xsend,xsendtot,1,MPI_REAL8,MPI_MAX,
     &           MPI_COMM_WORLD, ierr)
            xyzmax(idimension) = xsendtot
            splitcoords(nsplit(idimension),idimension) = 
     &           xyzmax(idimension)

            xsend = xyzmin(idimension)
            CALL MPI_ALLREDUCE(xsend,xsendtot,1,MPI_REAL8,MPI_MIN,
     &           MPI_COMM_WORLD, ierr)
            xyzmin(idimension) = xsendtot
         END DO
c
c--If not done before, need to guess splitcoords()
c
         IF (ifirst) THEN
            DO idimension = 1, 3
               DO number = 1, nsplit(idimension) - 1
                  splitcoords(number,idimension) = xyzmin(idimension) +
     &                 ((xyzmax(idimension)-xyzmin(idimension))/
     &                 nsplit(idimension))*number
               END DO
            END DO
            print *,iproc,' BAL: Guess coords ',
     &           ((splitcoords(number,idimension),
     &           number=1,nsplit(idimension)), idimension=1,3)
         ENDIF
c
c--Determine new coordinates of MPI domains
c
         DO idimension = 1, 3
            DO number = 1, nsplit(idimension) - 1
               idone = 0
               nfoundlast = 0
c
c--Target value for all processes and the number of x,y,z splits
c
               IF (iproc.EQ.0) THEN
                  ntarget = nactivetot/nsplit(idimension)

                  print *,iproc,' BAL ntarget ',ntarget
c
c--Use bisection
c
                  IF (number.EQ.1) THEN
                     dx = (splitcoords(number,idimension) -
     &                    xyzmin(idimension))/2.0
                  ELSE
                     dx = (splitcoords(number,idimension) - 
     &                    splitcoords(number-1,idimension))/2.
                  ENDIF
               ENDIF

 500           splitpos = splitcoords(number,idimension)
               CALL MPI_BCAST(splitpos,1,MPI_REAL8,0,
     &              MPI_COMM_WORLD,ierr)
               splitcoords(number,idimension) = splitpos

               IF (idone.EQ.1) GOTO 600
c
c--All MPI processes need to find how many active particles are in the 
c     specified zone
c
               nfound = 0
               IF (number.EQ.1) THEN
                  DO i = 1, npart
                     IF (iphase(i).GE.0 .AND. 
     & xyzmh(idimension,i).LT.splitcoords(number,idimension)) THEN
                        nfound = nfound + 1
                     ENDIF
                  END DO
               ELSEIF (number.EQ.nsplit(idimension)-1) THEN
                  DO i = 1, npart
                     IF (iphase(i).GE.0 .AND. 
     & xyzmh(idimension,i).GE.splitcoords(number,idimension)) THEN
                        nfound = nfound + 1
                     ENDIF
                  END DO
               ELSE
                  DO i = 1, npart
                     IF (iphase(i).GE.0 .AND. 
     & xyzmh(idimension,i).GE.splitcoords(number-1,idimension) .AND.
     & xyzmh(idimension,i).LT.splitcoords(number,idimension)) THEN
                        nfound = nfound + 1
                     ENDIF
                  END DO
               ENDIF

               print *,iproc,' BAL found ',nfound,' for ',number,
     &              idimension
c
c--Find total actual number of particles in this zone
c
               IF (iproc.EQ.0) THEN
                  DO i = 1, numproc-1
                     CALL MPI_RECV(nfoundremote, 1, MPI_INTEGER,
     &                    MPI_ANY_SOURCE, 161, MPI_COMM_WORLD, 
     &                    istatus, ierr)
                     nfound = nfound + nfoundremote
                  END DO

                  print *,iproc,' BAL nfound tot ',nfound,nfoundlast,dx
               ELSE
                  CALL MPI_SEND(nfound,1,MPI_INTEGER,0,161,
     &                 MPI_COMM_WORLD,ierr)
               ENDIF
c
c--Modify region to try and get closer to target
c
               IF (iproc.EQ.0) THEN
                  IF (nfound.GT.ntarget) THEN
                     splitcoords(number,idimension) = 
     &                    splitcoords(number,idimension) - dx
                  ELSE
                     splitcoords(number,idimension) =
     &                    splitcoords(number,idimension) + dx
                  ENDIF
                  dx = dx/2.0

                  IF ((nfound-nfoundlast).EQ.0) idone = 1
                  nfoundlast = nfound

                  print *,iproc,' BAL new ',splitcoords(number,
     &                 idimension),number,idimension
               ENDIF

               CALL MPI_BCAST(idone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

               GOTO 500

 600           CONTINUE
            END DO
         END DO

         print *,iproc,' BAL: New coords ',
     &        ((splitcoords(number,idimension),
     &        number=1,nsplit(idimension)), idimension=1,3)
      ENDIF

c
c--For successive binary splitting, process 0 controls the process of defining
c     where the new boundaries should be.
c     It broadcasts the new node boundaries to the other MPI processes.
c
      IF (mpitype.EQ.'b') THEN
         print *,iproc,' BAL entered mpitype.EQ.b'

         xtemp = LOG(REAL(numproc))/LOG(2.0)
         IF (xtemp - INT(xtemp+1.0E-4).GT.1.0E-4) THEN
            WRITE (*,*) 'ERROR - numproc is not a power of two ',
     &           numproc,xtemp
            CALL quit(1)
         ENDIF
         nsplits = INT(xtemp+1.0E-4)

         print *,iproc,': BAL number of splits is ',nsplits
c
c--Find total number of active particles on all MPI processes
c
         CALL MPI_ALLREDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,ierr)

         print *,iproc,' BAL nactivetot ',nactivetot
c
c--Need to find extent of particle distribution.
c
         DO idimension = 1, 3
            xyzmax(idimension) = -1.0E-30
            xyzmin(idimension) = 1.0E-30
            DO i = 1, npart
               xyzmax(idimension) = MAX(xyzmax(idimension),
     &              xyzmh(idimension,i))
               xyzmin(idimension) = MIN(xyzmin(idimension),
     &              xyzmh(idimension,i))
            END DO
         END DO
c
c--Find global bounds of active particles
c
         DO idimension = 1, 3
            xsend = xyzmax(idimension)
            CALL MPI_ALLREDUCE(xsend,xsendtot,1,MPI_REAL8,MPI_MAX,
     &           MPI_COMM_WORLD, ierr)
            xyzmax(idimension) = xsendtot

            xsend = xyzmin(idimension)
            CALL MPI_ALLREDUCE(xsend,xsendtot,1,MPI_REAL8,MPI_MIN,
     &           MPI_COMM_WORLD, ierr)
            xyzmin(idimension) = xsendtot
         END DO
c
c--Determine new coordinates of MPI domains
c
         nsplitbeingdone = 0
         DO i = 1, nsplits
            DO number = 1, 2**(i-1)
               nsplitbeingdone = nsplitbeingdone + 1
               idone = 0
               nfoundlast = 0
c
c--Find direction (x,y,z) of split
c
               idimension_now = MOD(i-1,3)+1
c
c--Set search space
c
               DO idimension = 1,3
                  searchxyzmax(idimension) = xyzmax(idimension)
                  searchxyzmin(idimension) = xyzmin(idimension)
               END DO

               DO i2 = 1, i-1
                  idimension = MOD(i2-1,3)+1

                  ipos = MOD( (nsplitbeingdone - 2**i2), 2**(i2-1)) + 
     &                 2**(i2-1)
                  
                  zeroorone = MOD((nsplitbeingdone-2**i2)/(2**(i2-1)),2)

                  IF (zeroorone.EQ.0) THEN
                     searchxyzmax(idimension) = binarysplits(ipos)
                  ELSE
                     searchxyzmin(idimension) = binarysplits(ipos)
                  ENDIF
               END DO
c
c--Set starting limits for bisection and dx value
c
               IF (iproc.EQ.0) THEN
                  boundupper = searchxyzmax(idimension_now)
                  boundlower = searchxyzmin(idimension_now)
                  dx = (boundupper - boundlower)/2.0
                  searchxyzmin(idimension_now) = boundlower + dx
                  dx = dx/2.0
c
c--Target value for all processes and the number of x,y,z splits
c
                  ntarget = nactivetot/(2**i)

                  print *,iproc,' BAL ntarget ',ntarget, 
     &                 idimension_now, dx, boundupper, boundlower
               ENDIF

 700           binarysplits(nsplitbeingdone) = 
     &              searchxyzmin(idimension_now)
               splitpos = binarysplits(nsplitbeingdone)
               CALL MPI_BCAST(splitpos,1,MPI_REAL8,0,
     &              MPI_COMM_WORLD,ierr)
               binarysplits(nsplitbeingdone) = splitpos
               searchxyzmin(idimension_now) = splitpos

               IF (iproc.EQ.0) THEN
                  print *,iproc,' BAL splitpos ',splitpos

                  print *,iproc,' BAL search ',searchxyzmin(1),
     &                 searchxyzmax(1),searchxyzmin(2),searchxyzmax(2),
     &                 searchxyzmin(3),searchxyzmax(3)
               ENDIF

               IF (idone.EQ.1) GOTO 800
c
c--All MPI processes need to find how many active particles are in the 
c     specified zone
c
               nfound = 0
               DO ipart = 1, npart
                  IF (iphase(ipart).GE.0) THEN
                     DO idimension = 1, 3
                        IF (xyzmh(idimension,ipart).LT.
     &                       searchxyzmin(idimension) .OR. 
     &                       xyzmh(idimension,ipart).GE.
     &                       searchxyzmax(idimension)) GOTO 750
                     END DO
                     nfound = nfound + 1
                  ENDIF
 750              CONTINUE
               END DO

               print *,iproc,' BAL found ',nfound,' for ',number,
     &              idimension
c
c--Find total actual number of particles in this zone
c
               IF (iproc.EQ.0) THEN
                  DO iloop = 1, numproc-1
                     CALL MPI_RECV(nfoundremote, 1, MPI_INTEGER,
     &                    MPI_ANY_SOURCE, 161, MPI_COMM_WORLD, 
     &                    istatus, ierr)
                     nfound = nfound + nfoundremote
                  END DO

                  print *,iproc,' BAL nfound tot ',nfound,nfoundlast,dx
               ELSE
                  CALL MPI_SEND(nfound,1,MPI_INTEGER,0,161,
     &                 MPI_COMM_WORLD,ierr)
               ENDIF
c
c--Modify region to try and get closer to target
c
               IF (iproc.EQ.0) THEN
                  IF (nfound.GT.ntarget) THEN
                     searchxyzmin(idimension_now) = 
     &                    searchxyzmin(idimension_now) + dx
                  ELSE
                     searchxyzmin(idimension_now) = 
     &                    searchxyzmin(idimension_now) - dx
                  ENDIF
                  dx = dx/2.0

                  IF (nfound.EQ.ntarget .OR. dx.LT.1.0E-10) idone = 1
                  nfoundlast = nfound

                  print *,iproc,' BAL new',
     &                 searchxyzmin(idimension_now),idimension_now
               ENDIF

               CALL MPI_BCAST(idone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

               GOTO 700

 800           CONTINUE
            END DO
         END DO

         print *,iproc,' BAL: New coords ',
     &        (binarysplits(i), i=1,nsplitbeingdone)

      ENDIF
c
c--With sink particle radiative feedback and MPI, need to make sure
c     that sink particles that are close to each other are on the same
c     MPI process.  This needs to be done iteratively until convergence.
c
      IF (istellarfeedback.GT.0) THEN
         DO i = 1, nptmasstot
            ibelong_sink(i) = classify(mpitype,stellar_xyzmhrti(1,i),0)
            PRINT 99100,iproc,i,ibelong_sink(i),stellar_xyzmhrti(1:3,i)
99100       FORMAT(I3,1X,' BAL: Initial sink classification ',2(I3,1X),
     &           3(1PE12.5,1x))
         END DO
 900     nchanged = 0
         DO i = 1, nptmasstot
            DO j = i+1, nptmasstot
               IF (ibelong_sink(i).NE.ibelong_sink(j)) THEN
                  distance2 = 
     &                 (stellar_xyzmhrti(1,i)-stellar_xyzmhrti(1,j))**2+
     &                 (stellar_xyzmhrti(2,i)-stellar_xyzmhrti(2,j))**2+
     &                 (stellar_xyzmhrti(3,i)-stellar_xyzmhrti(3,j))**2
                 PRINT *,' BAL: sink dist ',i,j,SQRT(distance2)*udist/au
                 IF (distance2.LT.(1000.0*au/udist)**2) THEN
                    ibelong_new = MIN(ibelong_sink(i),ibelong_sink(j))
                    ibelong_sink(i) = ibelong_new
                    ibelong_sink(j) = ibelong_new
                    nchanged = nchanged + 1
                 ENDIF
              ENDIF
            END DO
         END DO
         IF (nchanged.NE.0) GOTO 900

         DO i = 1, nptmasstot
            print *,iproc,' BAL: Final sink classification ',
     &           i,ibelong_sink(i)
         END DO
      ENDIF
c
c--Now classify every particle as to which MPI process it should belong to
c
c     If istellarfeedback.GT.0, need to make sure that sink particles
c     near an MPI boundary have sufficient particles around them on
c     the same MPI process as them to give their stellar feedback too.
c
      DO ipart = 1, npart
         IF (istellarfeedback.GT.0) THEN
            IF (iphase(ipart).GE.1 .AND. iphase(ipart).LE.9) THEN
               ibelong = ibelong_sink(liststellarID(listrealpm(ipart)))
               PRINT *,iproc,' BAL: ibelong of sink ',ipart,
     &              listrealpm(ipart),liststellarID(listrealpm(ipart)),
     &              ibelong

               ibelong1 = ibelong + 1
               numlist(ibelong1) = numlist(ibelong1) + 1
               listtrans(ipart) = ibelong
            ELSEIF (iphase(ipart).GE.0) THEN
               ibelong = classify(mpitype,xyzmh(1,ipart),1)
               
               ibelong1 = ibelong + 1
               numlist(ibelong1) = numlist(ibelong1) + 1
               listtrans(ipart) = ibelong
            ELSE
               listtrans(ipart) = -1
            ENDIF
         ELSE
            IF (iphase(ipart).GE.0) THEN
               ibelong = classify(mpitype,xyzmh(1,ipart),0)
               
               ibelong1 = ibelong + 1
               numlist(ibelong1) = numlist(ibelong1) + 1
               listtrans(ipart) = ibelong
            ELSE
               listtrans(ipart) = -1
            ENDIF
         ENDIF
      END DO
c
c--End of MPI-only section
c
#endif

      IF (ifirst) THEN
         ifirst = .FALSE.
      ENDIF

      RETURN
      END

c---------------------------------------------------------------------------

      SUBROUTINE move_particle(ifrom,ito)
c************************************************************
c                                                           *
c     Moves a particle from location ifrom to ito (within   *
c     an MPI process, not between MPI processes)            *
c                                                           *
c     NOTE: Does not alter listbins for timesteps or ireal  *
c     or hasghost for ghost particles.                      *
c     Does not alter tree.                                  *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

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
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/accnum'
      INCLUDE 'COMMONS/accurpt'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/radsink'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/dustimplicit'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif
      INCLUDE 'COMMONS/stellarradiation'

      DO i = 1, 5
         xyzmh(i,ito) = xyzmh(i,ifrom)
         dumxyzmh(i,ito) = dumxyzmh(i,ifrom)
      END DO
      DO i = 1, 4
         vxyzu(i,ito) = vxyzu(i,ifrom)
         dumvxyzu(i,ito) = dumvxyzu(i,ifrom)
         f1vxyzu(i,ito) = f1vxyzu(i,ifrom)
         f2vxyzu(i,ito) = f2vxyzu(i,ifrom)         
      END DO
      isteps(ito) = isteps(ifrom)
      IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) THEN
         icolumnsteps(ito) = icolumnsteps(ifrom)
         icolumnnext(ito) = icolumnnext(ifrom)
      ENDIF
      it0(ito) = it0(ifrom)
      it1(ito) = it1(ifrom)
      it2(ito) = it2(ifrom)
      iphase(ito) = iphase(ifrom)
      IF (iphase(ifrom).GE.1 .AND. iphase(ifrom).LT.10) THEN
         ipt = listrealpm(ifrom)
         IF (ipt.LE.0 .OR. ipt.GT.nptmass) THEN
            WRITE (*,*) 'ERROR - ipt.LE.0 .OR. ipt.GT.nptmass'
            CALL quit(1)
         ENDIF
         listrealpm(ito) = ipt
         listpm(ipt) = ito
         listrealpm(ifrom) = 0
      ENDIF
      iphase(ifrom) = -1
c      iorig(ito) = ito
c      isort(ito) = ito
c      iorig(ifrom) = 0
c      isort(ifrom) = 0
      iunique(iorig(ito)) = iunique(iorig(ifrom))
      iunique(iorig(ifrom)) = 0
      rho(ito) = rho(ifrom)
      dumrho(ito) = dumrho(ifrom)
      divv(ito) = divv(ifrom)
      curlv(ito) = curlv(ifrom)
      ddv(ito) = ddv(ifrom)
      pr(ito) = pr(ifrom)
      vsound(ito) = vsound(ifrom)
      poten(ito) = poten(ifrom)
      IF (idim_h2.EQ.idim) THEN
c         h2mol(ito) = h2mol(ifrom)
         h2ratio(ito) = h2ratio(ifrom)
      ENDIF
      IF (idustIMPL.GT.0) dustnorm(ito) = dustnorm(ifrom)
      iscurrent(ito) = iscurrent(ifrom)
      notacc(ito) = notacc(ifrom)
      DO i = 1, 2
         gradhs(i,ito) = gradhs(i,ifrom)
      END DO
      DO i = 1, isizealphaMM
         alphaMM(i,ito) = alphaMM(i,ifrom)
         dumalpha(i,ito) = dumalpha(i,ifrom)
      END DO
      DO i = 1, 1+isizealphaMM
         f1ha(i,ito) = f1ha(i,ifrom)
         f2ha(i,ito) = f2ha(i,ifrom)
      END DO
      IF (imhd.EQ.idim) THEN
         DO i = 1, imhdevol
            Bevolxyz(i,ito) = Bevolxyz(i,ifrom)
            dumBevolxyz(i,ito) = dumBevolxyz(i,ifrom)
            f1Bxyz(i,ito) = f1Bxyz(i,ifrom)
            f2Bxyz(i,ito) = f2Bxyz(i,ifrom)
         END DO
         DO i = 1, 3
            Bxyz(i,ito) = Bxyz(i,ifrom)
         END DO
         DO i = 1, 5
            divcurlB(i,ito) = divcurlB(i,ifrom)
         END DO
#ifdef NONIDEAL
         DO i = 1, 3
            jcurrent(i,ito) = jcurrent(i,ifrom)
         END DO
         DO i = 1,4
            n_R(i,ito) = n_R(i,ifrom)
         END DO
         n_electronT(ito) = n_electronT(ifrom)
         DO i = 1,4
            eta_nimhd(i,ito) = eta_nimhd(i,ifrom)
         END DO
#endif
      ENDIF
      IF (iradtrans.EQ.idim) THEN
         DO i = 1, 5
            ekcle(i,ito) = ekcle(i,ifrom)
            dumekcle(i,ito) = dumekcle(i,ifrom)
         END DO
         IF (idustRT.GT.0) THEN
            DO i = 1, 2
               dust_tk(i,ito) = dust_tk(i,ifrom)
            END DO
            DO i = 1, nheatingISR
               heatingISR(i,ito) = heatingISR(i,ifrom)
            END DO
            IF (ioptimise_column.EQ.1) THEN
               DO i = 1, 2
                  heatingISRold(i,ito) = heatingISRold(i,ifrom)
               END DO
            ENDIF
            DO i = 1, nchemistry
               chemistry(i,ito) = chemistry(i,ifrom)
            END DO
            h2frac(ito) = h2frac(ifrom)
            DO i = 1, 2
               dh2dt(i,ito) = dh2dt(i,ifrom)
            END DO
         ENDIF
      ENDIF
      IF (istellar_ionisation) THEN
         HIIion(ito) = HIIion(ifrom)
         dHIIdt(1:2,ito) = dHIIdt(1:2,ifrom)
         stellarrad(1:4,1:iptdim,ito) =stellarrad(1:4,1:iptdim,ifrom)
      ENDIF

      RETURN
      END

c===========================================================================
