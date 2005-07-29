      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &              dvxyzu,dha,npart,ntot,ireal,alphaMM)
c************************************************************
c                                                           *
c  This subroutine drives the computation of the forces on  *
c     every particle on the list.                           *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,idim),vxyzu(4,idim),dvxyzu(4,idim)
      REAL*4 dha(2,idim)
      DIMENSION ireal(idim)
      DIMENSION alphaMM(idim)

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

c      INCLUDE 'COMMONS/timeextra'

      CHARACTER*7 where

      DATA where/'derivi'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine derivi')
c
c--Set constants first time around
c
      nlst_in = 1
      nlst_end = nlst
      IF (itrace.EQ.'all') WRITE (iprint, 99002) nlst_in, nlst_end
99002 FORMAT(' derivi ',I8, I8)
c
c--Compute the neighbor indexes & gravitational forces of the distant 
c     particles for all the particles in the list
c
      IF (igrape.EQ.0) THEN
         IF (nlst_end.GT.nptmass) THEN
            CALL insulate(3, ntot, npart, xyzmh, dvxyzu)
         ELSE
c--Don't bother to update gravity from gas particles acting on sinks
c    Also means potential energy and neighbours of sinks are not updated
            DO i = nlst_in, nlst_end
               ipart = llist(i)
               dvxyzu(1,ipart) = gravxyzstore(1,ipart)
               dvxyzu(2,ipart) = gravxyzstore(2,ipart)
               dvxyzu(3,ipart) = gravxyzstore(3,ipart)
               poten(ipart) = potenstore(ipart)
            END DO
         ENDIF
      ELSEIF (igrape.EQ.1) THEN
         CALL insulate(4, ntot, npart, xyzmh, dvxyzu)
      ELSE
         CALL error(where,1)
      ENDIF
c
c--Find neighbours and calculate gravity on and from point masses
c     The point masses are no longer in the TREE/GRAPE - done separately for
c     higher accuracy
c
      IF (nptmass.GT.0) THEN
         IF (iptintree.EQ.0) THEN
            CALL gptall(xyzmh,npart)
         ELSEIF (iptintree.EQ.1) THEN
            CALL gforspt(xyzmh,dvxyzu)
         ENDIF
      ENDIF

c
c--Only want to do density and updates for non-sinks
c
      IF (nlst_end.GT.nptmass) THEN
c
c--Compute the pressure, divv, etc. on list particles but
c     do not search twice the list particles
c
         IF (itiming) CALL getused(tdens1)

         jneigh = 0

c         CALL densityi(jneigh,npart,xyzmh,vxyzu,
c     &        nlst_in,nlst_end,llist,itime)
         CALL iterate_density(jneigh,npart,xyzmh,vxyzu,
     &        nlst_in,nlst_end,llist,itime)

         IF (itiming) THEN
            CALL getused(tdens2)
            tdens = tdens + (tdens2 - tdens1)
         ENDIF
c
c--Predict the pressure, divv, etc. on list-particle neighbors
c
         IF (itiming) CALL getused(td11)
         nlistavail = 0
         DO i = 1, nlst
            ipart = llist(i)
            DO j = 1, nneigh(ipart)
               IF (j.GE.nlmax) THEN
                  jpart = neighover(j-nlmax+1,ABS(neighb(nlmax,ipart))) 
               ELSE
                  jpart = neighb(j,ipart)
               ENDIF
               IF (jpart.GT.idim) THEN
                  WRITE (*,*) 'jpart.GT.idim'
                  WRITE (iprint,*) 'jpart.GT.idim'
                  CALL quit
               ENDIF
               IF (jpart.LT.1) THEN
                  WRITE (*,*) 'jpart.LT.1'
                  WRITE (iprint,*) 'jpart.LT.1'
                  CALL quit
               ENDIF
               IF (iavail(jpart).EQ.0 .AND.(.NOT.iscurrent(jpart))) THEN
                  iavail(jpart) = 1
                  nlistavail = nlistavail + 1
                  IF (nlistavail+nlst.GT.idim) THEN
                     WRITE (*,*) 'nlistavail+nlst.GT.idim'
                     WRITE (iprint,*) 'nlistavail+nlst.GT.idim'
                     CALL quit
                  ENDIF
                  llist(nlst+nlistavail) = jpart
               ENDIF
            END DO
         END DO
         IF (itiming) THEN
            CALL getused(td12)
            td1 = td1 + (td12 - td11)
         ENDIF

         IF (itiming) CALL getused(td21)
C$OMP PARALLEL default(none)
C$OMP& shared(npart,it1,imax,itime,it0,isteps,dt,imaxstep)
C$OMP& shared(divv,rho,dumrho,vxyzu,pr,vsound,ntot,ireal)
C$OMP& shared(iphase,nlst,llist,nlistavail,iavail)
C$OMP& private(i,j,deltat,deltarho,ipart)
C$OMP DO SCHEDULE (runtime)
         DO i = 1, nlst+nlistavail
            ipart = llist(i)
            IF (i.LE.nlst) THEN
               dumrho(ipart) = rho(ipart)
            ELSE
               iavail(ipart) = 0

               IF (ipart.LE.npart) THEN
                  IF (it1(ipart).EQ.imax) THEN
            deltat = dt*(itime - it0(ipart) - isteps(ipart)/2)/imaxstep
                  ELSE
                     deltat = dt*(itime - it0(ipart))/imaxstep
                  ENDIF
c
c--Update the density value at neighbor's locations
c--Avoid, though, abrupt changes in density
c
                  deltarho = -deltat*divv(ipart)
                  IF (ABS(deltarho).GT.rho(ipart)/2.) THEN
                     deltarho = SIGN(1.0,deltarho)*rho(ipart)/2.0
                  ENDIF
                  dumrho(ipart) = rho(ipart) + deltarho

                  CALL eospg(ipart,vxyzu,dumrho,pr,vsound)
               ENDIF
            ENDIF
         END DO
C$OMP END DO

C$OMP DO SCHEDULE (runtime)
         DO i = npart + 1, ntot
            j = ireal(i)
            dumrho(i) = dumrho(j)
            pr(i) = pr(j)
            vsound(i) = vsound(j)
            divv(i) = 0.
            vxyzu(4,i) = vxyzu(4,j)
         END DO
C$OMP END DO
C$OMP END PARALLEL

         IF (itiming) THEN
            CALL getused(td22)
            td2 = td2 + (td22 - td21)
         ENDIF
c
c--End if for nlst>nptmass
c
      ENDIF
c
c--Compute forces on EACH particle
c
      jneigh = 0

      IF (itiming) CALL getused(tforce1)

      CALL forcei(jneigh,nlst_in,nlst_end,llist,dt,itime,npart,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

      RETURN
      END
