      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &     dvxyzu,dha,npart,ntot,ireal,alphaMM,ekcle,Bevolxyz,dBxyz)
c************************************************************
c                                                           *
c  This subroutine drives the computation of the forces on  *
c     every particle on the list.                           *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax),vxyzu(4,idim),dvxyzu(4,idim)
      REAL*4 dha(2,idim),alphaMM(idim)
      DIMENSION ireal(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevolxyz(3,imhd)
      DIMENSION dBxyz(3,imhd),Bxyz(3,imhd)

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

c      INCLUDE 'COMMONS/timeextra'

      CHARACTER*7 where
      DIMENSION dedxyz(3,iradtrans)

      DATA where/'derivi'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine derivi')
c
c--Set constants first time around
c
      uradconst = radconst/uergcc

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

c         CALL densityi(npart,xyzmh,vxyzu,
c     &        nlst_in,nlst_end,llist,itime)
         CALL iterate_density(npart,xyzmh,vxyzu,dha,
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
         nlstall = nlst + nlistavail
         IF (itiming) THEN
            CALL getused(td12)
            td1 = td1 + (td12 - td11)
         ENDIF

         IF (itiming) CALL getused(td21)
C$OMP PARALLEL default(none)
C$OMP& shared(npart,it1,imax,itime,it0,isteps,dt,imaxstep,ekcle)
C$OMP& shared(divv,rho,dumrho,vxyzu,pr,vsound,ntot,ireal,ibound,icall)
C$OMP& shared(iphase,nlst,llist,nlistavail,iavail,encal,uradconst)
C$OMP& private(i,j,deltat,deltarho,ipart)
C$OMP DO SCHEDULE (runtime)
         DO i = 1, nlst+nlistavail
            ipart = llist(i)
            IF (i.LE.nlst) THEN
               dumrho(ipart) = rho(ipart)
c
c--Set e(i) for the first time around because density only set here
c     This sets radiation and matter to have same initial temperature
c
               IF (icall.EQ.1 .AND. encal.EQ.'r') THEN
                  IF (ekcle(1,ipart).EQ.0.0) THEN
                     ekcle(3,ipart) = getcv(rho(ipart),vxyzu(4,ipart))
                     ekcle(1,ipart) = uradconst*(vxyzu(4,ipart)/
     &                    ekcle(3,ipart))**4/rho(ipart)
                  ENDIF
               ENDIF
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
            IF (encal.EQ.'r' .AND. ibound.EQ.100)
     &           ekcle(1,i) = ekcle(1,j)
         END DO
C$OMP END DO
c
c--make sure that B is sent in to force (instead of the evolved variable B/rho)
c
      IF (imhd.EQ.idim) THEN
C$OMP    DO SCHEDULE (runtime)
         DO i=1,ntot
            DO j=1,3
               Bxyz(j,i) = Bevolxyz(j,i)*dumrho(i)
            ENDDO
         ENDDO
C$OMP    END DO
cc      ELSE
ccC$OMP    DO SCHEDULE (runtime)
cc         DO i=1,ntot
cc            DO j=1,3
cc               Bxyz(j,i) = Bevolxyz(j,i)
cc            ENDDO
cc         ENDDO
ccC$OMP    END DO
      ENDIF
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
c--Compute implicit radiative transfer
c
      IF (itiming) CALL getused(tass1)

      IF(encal.EQ.'r') THEN
         WRITE (*,*) 'Calling ass at realtime ',dt*itime/imaxstep+gt
c         CALL ASS(nlst_in,nlst_end,nlstall,llist,dt,itime,npart,
c     &        xyzmh,vxyzu,ekcle,dumrho,vsound,dedxyz)

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlstall,vxyzu,dumrho,pr,vsound,llist)
C$OMP& private(i,ipart)
         DO i = 1, nlstall
            ipart = llist(i)
            CALL eospg(ipart,vxyzu,dumrho,pr,vsound)
         END DO
C$OMP END PARALLEL DO
      END IF

      IF (itiming) THEN
        CALL getused(tass2)
        tass = tass + (tass2 - tass1)
      ENDIF
c
c--Compute forces on EACH particle
c
      IF (itiming) CALL getused(tforce1)

      CALL forcei(nlst_in,nlst_end,llist,dt,itime,npart,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,
     &     ekcle,dedxyz,Bxyz,dBxyz)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

      RETURN
      END
