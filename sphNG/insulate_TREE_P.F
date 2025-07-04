      SUBROUTINE insulate(itype, itime, ntot, npart, xyzmh, fxyzu)
c************************************************************
c                                                           *
c  The purpose of this subroutine is to insulate the main   *
c     SPH code from being able to tell whether it uses a    *
c     binary tree or the GRAPE to get gravity forces and    *
c     neighbours.                                           *
c                                                           *
c     Written by M R Bate 2/2/96                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax2), fxyzu(4,idim3)

      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/active'

      CHARACTER*7 where
      LOGICAL*1 ireset
      DATA where/'insul'/

      IF (itiming) CALL getused(tins1)

c------------------------------------------------------------
c--Call mtree to make the binary tree
c------------------------------------------------------------

      IF (itype.EQ.1) THEN

         IF (itiming) CALL getused(tmtree1)

         CALL mtree(ntot, npart, xyzmh)

         itreeupdate = .TRUE.

         IF (itiming) THEN
            CALL getused(tmtree2)
            tmtree = tmtree + (tmtree2 - tmtree1)
         ENDIF

c------------------------------------------------------------
c--Call revtree to revise the binary tree
c------------------------------------------------------------

      ELSEIF (itype.EQ.2) THEN

         IF (itiming) CALL getused(trevt1)

         CALL revtree(ntot, npart, xyzmh)

         itreeupdate = .TRUE.

         IF (itiming) THEN
            CALL getused(trevt2)
            trevt = trevt + (trevt2 - trevt1)
         ENDIF

c------------------------------------------------------------
c--Get gravity forces and neighbours using the binary tree
c------------------------------------------------------------

      ELSEIF (itype.EQ.3) THEN
         IF (itiming) CALL getused(ttreef1)
c
c--Compute the neighbour indexes & gravitational forces of the distant 
c     particles for all the particles in the list
c
c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , I1'
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , I1'
c            WRITE (*,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            WRITE (iprint,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            CALL quit(1)
c         ENDIF
c      END DO

         IF (nlst.EQ.nactive) THEN
            ireset = .TRUE.
            numoverflowmax = MAX(numoverflowmax,numoverflow)
            numoverflow = 0
         ELSE
            ireset = .FALSE.
         ENDIF

C$OMP PARALLEL default(none), shared(nlst,npart,ntot,acc,igphi)
C$OMP& shared(llist,xyzmh,fxyzu,poten,iphase)
C$OMP& shared(ireset,neighb,itime)
C$OMP& private(i,ipart,fsx,fsy,fsz,epot)
C$OMP DO SCHEDULE(runtime)
         DO 100 i = 1, nlst
            ipart = llist(i)
            IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10 .AND. 
     &           iptintree.EQ.0) GOTO 100
            IF (ireset) neighb(nlmax,ipart) = 0
c
c--Walk through the tree to get the neighbours, the acceleration
c     and the potential due to outside 2h particles
c
            CALL treef(ipart,itime,npart,ntot,xyzmh,acc,igphi,
     &           fsx,fsy,fsz,epot)

            IF (igphi.GT.0 .OR. (iphase(ipart).GE.1 .AND. 
     &               iphase(ipart).LT.10 .AND. iptintree.EQ.1)) THEN
               fxyzu(1,ipart) = fsx
               fxyzu(2,ipart) = fsy
               fxyzu(3,ipart) = fsz
               poten(ipart) = epot
            ELSE
               fxyzu(1,ipart) = 0.
               fxyzu(2,ipart) = 0.
               fxyzu(3,ipart) = 0.
               poten(ipart) = 0.
            ENDIF
 100     CONTINUE
C$OMP END DO
C$OMP END PARALLEL

c      DO i = 1, nlstbins(29)-1
c         IF (listbins(i,29).NE.i+94296) THEN
c            WRITE (*,*) 'listbins(i,29).NE.i+94296 , I2'
c            WRITE (iprint,*) 'listbins(i,29).NE.i+94296 , I2'
c            WRITE (*,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            WRITE (iprint,*) i,nlstbins(29),listbins(i-1,29),
c     &           listbins(i,29),listbins(i+1,29),listbins(i+2,29)
c            CALL quit(1)
c         ENDIF
c      END DO

         IF (itiming) THEN
            CALL getused(ttreef2)
            ttreef = ttreef + (ttreef2 - ttreef1)
         ENDIF

c------------------------------------------------------------
c--Get gravity forces and neighbours using the binary tree
c------------------------------------------------------------

      ELSEIF (itype.EQ.5) THEN
         IF (itiming) CALL getused(ttreef1)
c
c--Compute the neighbour indexes & gravitational forces of the distant
c     particles for all the particles in the list
c

C$OMP PARALLEL default(none), shared(nlst,npart,ntot,acc,igphi)
C$OMP& shared(llist,xyzmh,iphase,itime)
C$OMP& private(i,ipart,fsx,fsy,fsz,epot)
C$OMP DO SCHEDULE(runtime)
         DO 200 i = 1, nlst
            ipart = llist(i)
            IF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10 .AND. 
     &           iptintree.EQ.0) GOTO 200
c
c--Walk through the tree to get the neighbours, the acceleration
c     and the potential due to outside 2h particles
c
            CALL treef(ipart,itime,npart,ntot,xyzmh,acc,igphi,
     &           fsx,fsy,fsz,epot)

 200     CONTINUE
C$OMP END DO
C$OMP END PARALLEL

         IF (itiming) THEN
            CALL getused(ttreef2)
            ttreef = ttreef + (ttreef2 - ttreef1)
         ENDIF

      ELSE
         CALL error(where,2)
      ENDIF

      IF (itiming) THEN
         CALL getused(tins2)
         tins = tins + (tins2 - tins1)
      ENDIF

      RETURN
      END


c------------------------------------------------------------

      SUBROUTINE quit(iflag)
c************************************************************
c                                                           *
c  This subroutine should be used instead of 'STOP'         *
c                                                           *
c************************************************************

#ifdef __INTEL_COMPILER
c Needed so that ABORT works correctly
      USE IFLPORT, ONLY: ABORT
#endif
#ifdef MPIALL
#define IMPLICITNONE
#include "mpi_sup.h"
#endif
      INCLUDE 'idim'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/actio'

      INTEGER iflag

#ifdef MPIALL
      INTEGER iret
      INCLUDE 'COMMONS/mpiall'
#endif

#ifndef MPIALL
      LOGICAL ionce
      DATA ionce/.TRUE./
#endif

c First, if iflag < 0 this is a graceful exit
      IF (iflag.LT.0) THEN
#ifdef MPIALL
c Shut down the MPI system
        CALL MPI_FINALIZE(ierr)
#endif
        GOTO  168
#ifndef MPIALL
      ELSEIF (iflag.GT.0.AND.ionce.AND.job(1:9).EQ.'evolution') THEN
c Prevent infinite loops
         ionce = .FALSE.
         WRITE(*,99999)
         WRITE(iprint,99999)
99999    FORMAT(' Quitting. Trying to dump state-of-play to FAILDUMP')

c Make sure the unit number is free
         CLOSE(UNIT=idisk2)
         OPEN(UNIT=idisk2, ERR=166, FILE="FAILDUMP",
     &        FORM='unformatted', RECL=imaxrec,
     &        STATUS='replace', ACTION='write')
         CALL wdump(idisk2)
         CLOSE(UNIT=idisk2)

         GOTO 167

 166     WRITE(*,99998)
         WRITE(iprint,99998)
99998    FORMAT(' Failed to write dump. Sorry.')
#endif
#ifdef MPIALL
      ELSE
c iret should probably be set to something 'meaningful'
         iret = 65
         CALL MPI_ABORT(MPI_COMM_WORLD, iret, ierr)
#endif
      ENDIF
c ABORT allows for tracebacks and core files (sometimes)
 167  CALL ABORT
c If by some error we make it this far, really really STOP
 168  STOP 
      END
