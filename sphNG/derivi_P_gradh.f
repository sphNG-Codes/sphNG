      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &     dvxyzu,dha,npart,ntot,ireal,alphaMM,ekcle,Bxyz,dBxyz)
c************************************************************
c                                                           *
c  This subroutine drives the computation of the forces on  *
c     every particle on the list.                           *
c                                                           *
c  MRB 14/12/2005:                                          *
c                                                           *
c  The non-gradh version calculates neighbours and gravity  *
c     with neighbours stored in a list, then computes a     *
c     list of particles whose values are required to be     *
c     interpolated (density, pressure, sound speed) and     *
c     does the interpolation, then updates ghosts, the      *
c     does implicit radiative transfer, and finally calls   *
c     forcei to calculate the rest of the forces.           *
c                                                           *
c  In this grad-h version of the code, the structure is     *
c     very different because of the interations required to *
c     set density and the fact that enormous numbers of     *
c     neighbours can be obtained to make density and h      *
c     consistent (overflowing any global neighbour store).  *
c     Neighbours are no longer stored. Rather, they are got *
c     each time for density and forcei. Note that this will *
c     not work for radiative transfer because of the        *
c     iterative solve - would need to calculate neighbours  *
c     each iteration !!!!  In the grad-h version of the     *
c     code, density and neighbours (h's) are calculated     *
c     simulataneously in densityiterate_gradh.f.  This code *
c     also calculates divv() and other quantities as in the *
c     normal density.f code.  It also interpolates the      *
c     density, pressure and sound speed of neighbouring     *
c     particles that are required for forces at the same    *
c     time, rather than constructing a list and doing it    *
c     later.  Thus, it combines density.f and the middle    *
c     part of derivi.  Forcei_gradh then calculates all the *
c     forces, including all gravity forces when it re-finds *
c     the neighbours of each particle.  Thus, it combines   *
c     the first part of derivi.f and forcei.f .             *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax),vxyzu(4,idim),dvxyzu(4,idim)
      REAL*4 dha(2,idim),alphaMM(idim)
      DIMENSION ireal(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bxyz(3,imhd),dBxyz(3,imhd)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/hagain'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/avail'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/gtime'

      CHARACTER*7 where
      DIMENSION dedxyz(3,iradtrans)

      DATA where/'derivi'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine derivi_gradh')
      IF (igrape.NE.0) THEN
         WRITE (iprint,*) 'ERROR: derivi_P_gradh must have igrape.EQ.0'
         CALL quit
      ENDIF
      IF (nlmax.NE.1) THEN
         WRITE (iprint,*) 'ERROR: derivi_P_gradh must have nlmax.EQ.1'
         CALL quit
      ENDIF
c     
c--Set constants first time around
c
      uradconst = radconst/uergcc
      nlst_in = 1
      nlst_end = nlst
      IF (itrace.EQ.'all') WRITE (iprint, 99002) nlst_in, nlst_end
99002 FORMAT(' derivi_gradh ',I8, I8)
c
c--Find self-consistent density and smoothing length for all the particles
c     in the list.  Only calculate density etc for non-sinks.
c
      IF (nlst_end.GT.nptmass) THEN
         IF (itiming) CALL getused(tdens1)

         CALL densityiterate_gradh(dt,npart,ntot,xyzmh,vxyzu,
     &        nlst_in,nlst_end,llist,itime,ekcle)

         IF (itiming) THEN
            CALL getused(tdens2)
            tdens = tdens + (tdens2 - tdens1)
         ENDIF
c     
c--Compute implicit radiative transfer
c     
         IF (itiming) CALL getused(tass1)

         IF(encal.EQ.'r') THEN
            WRITE (*,*) 'Calling ass at realtime ',dt*itime/imaxstep+gt
            WRITE(*,*) 'radiation/gradh not yet implemented'
            STOP

c            WRITE (*,*) 'nlstall is ',nlstall

c            CALL ASS(nlst_in,nlst_end,nlstall,llist,dt,itime,npart,
c     &           xyzmh,vxyzu,ekcle,dumrho,vsound,dedxyz)

cC$OMP PARALLEL DO SCHEDULE(runtime) default(none)
cC$OMP& shared(nlstall,vxyzu,dumrho,pr,vsound,llist,iphase,ekcle)
cC$OMP& private(i,ipart)
c            DO i = 1, nlstall
c               ipart = llist(i)
c               IF (iphase(ipart).EQ.0) 
c     &              CALL eospg(ipart,vxyzu,dumrho,pr,vsound,ekcle)
c            END DO
cC$OMP END PARALLEL DO
         END IF

         IF (itiming) THEN
            CALL getused(tass2)
            tass = tass + (tass2 - tass1)
         ENDIF
      ENDIF
c
c--Compute forces on EACH particle
c
      IF (itiming) CALL getused(tforce1)

      CALL forcei(nlst_in,nlst_end,llist,dt,itime,npart,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,
     &     dedxyz,Bxyz,dBxyz)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

      RETURN
      END
