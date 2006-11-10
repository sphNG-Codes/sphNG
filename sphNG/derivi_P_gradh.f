      SUBROUTINE derivi (dt,itime,xyzmh,vxyzu,
     &     dvxyzu,dha,npart,ntot,ireal,alphaMM,ekcle,Bevolxyz,dBevolxyz)
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
c  MHD note: (DJP 6.1.06)                                   *
c   The MHD quantities passed to this routine               *
c     are the *evolved* MHD variables. These could be       *
c     B, B/rho (usual option) or the Euler potentials.      *
c     However, we send just the magnetic field into the     *
c     force routines. The dBevolxyz returned by derivi      *
c     is the derivative required for evolving the magnetic  *
c     field.                                                *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      DIMENSION xyzmh(5,mmax),vxyzu(4,idim),dvxyzu(4,idim)
      REAL*4 dha(1+isizealphaMM,idim),alphaMM(isizealphaMM,idim)
      DIMENSION ireal(idim)
      DIMENSION ekcle(5,iradtrans)
      DIMENSION Bevolxyz(3,imhd),dBevolxyz(3,imhd)

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
c     Bxyz is stored here for calculation of energy & writing to dump file
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/vsmooth'
      INCLUDE 'COMMONS/divcurlB'

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
c     for Euler potentials, this also means get B from grad alpha x grad beta
c
      IF (nlst_end.GT.nptmass) THEN
         IF (itiming) CALL getused(tdens1)

         CALL densityiterate_gradh(dt,npart,ntot,xyzmh,vxyzu,
     &        nlst_in,nlst_end,llist,itime,ekcle,Bevolxyz,Bxyz)

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
c--Calculate B from the evolved magnetic field variable
c
      IF (imhd.EQ.idim) THEN
         IF (varmhd.eq.'Bvol') THEN
            DO i=nlst_in,nlst_end
               ipart = llist(i)
               Bxyz(1,ipart) = Bevolxyz(1,ipart)
               Bxyz(2,ipart) = Bevolxyz(2,ipart)
               Bxyz(3,ipart) = Bevolxyz(3,ipart)
c
c--smooth magnetic field variable
c               
c               Bevolxyz(1,ipart) = Bsmooth(1,ipart)
c               Bevolxyz(2,ipart) = Bsmooth(2,ipart)
c               Bevolxyz(3,ipart) = Bsmooth(3,ipart)
            ENDDO
         ELSEIF (varmhd.EQ.'Brho') THEN
            DO i=nlst_in,nlst_end
               ipart = llist(i)
c
c--smooth magnetic field variable
c               
c               Bevolxyz(1,ipart) = Bsmooth(1,ipart)
c               Bevolxyz(2,ipart) = Bsmooth(2,ipart)
c               Bevolxyz(3,ipart) = Bsmooth(3,ipart)
c
c--used smoothed B in dBevol/dt *and* in forces
c
               Bxyz(1,ipart) = Bevolxyz(1,ipart)*dumrho(ipart)
               Bxyz(2,ipart) = Bevolxyz(2,ipart)*dumrho(ipart)
               Bxyz(3,ipart) = Bevolxyz(3,ipart)*dumrho(ipart)

c               Bsmooth(1,ipart) = Bsmooth(1,ipart)*dumrho(ipart)
c               Bsmooth(2,ipart) = Bsmooth(2,ipart)*dumrho(ipart)
c               Bsmooth(3,ipart) = Bsmooth(3,ipart)*dumrho(ipart)
            ENDDO
         ELSEIF (varmhd.NE.'eulr') THEN
            STOP 'unknown MHD variable in derivi'
         ENDIF
      ENDIF
c
c--Implicit hyperdiffusion of div B
c
      IF (itiming) CALL getused(tass1)

c      IF(imhd.EQ.idim) THEN
      IF(.FALSE.) THEN
         WRITE (*,*) 'Calling Hyper at realtime ',dt*itime/imaxstep+gt
         CALL divBdiffuse(dt,nlst_in,nlst_end,npart,llist,
     &        xyzmh,dumrho,Bxyz)

         IF (varmhd.EQ.'Bvol') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,Bevolxyz,Bxyz,llist)
C$OMP& private(i,ipart)
            DO i = nlst_in,nlst_end
               ipart = llist(i)
               Bevolxyz(1,ipart) = Bxyz(1,ipart)
               Bevolxyz(2,ipart) = Bxyz(2,ipart)
               Bevolxyz(3,ipart) = Bxyz(3,ipart)
            END DO
C$OMP END PARALLEL DO
         ELSEIF (varmhd.EQ.'Brho') THEN
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,Bevolxyz,Bxyz,llist,dumrho)
C$OMP& private(i,ipart)
            DO i = nlst_in,nlst_end
               ipart = llist(i)
               Bevolxyz(1,ipart) = Bxyz(1,ipart)/dumrho(ipart)
               Bevolxyz(2,ipart) = Bxyz(2,ipart)/dumrho(ipart)
               Bevolxyz(3,ipart) = Bxyz(3,ipart)/dumrho(ipart)
            END DO
C$OMP END PARALLEL DO
         ENDIF
      END IF
c
c--div B projection
c
c      IF (imhd.EQ.idim) THEN
c         IF (varmhd.EQ.'Brho' .OR. varmhd.EQ.'Bvol') THEN
c            divBmax = 0.
c            DO i=1,npart
c               divBmax = max(divcurlB(1,i),divBmax)
c            ENDDO
c            WRITE(iprint,*) 'div B max = ',divBmax
c            CALL divBclean(nlst_in,nlst_end,npart,ntot,llist,
c     &                     xyzmh,rho,Bevolxyz)
c         ENDIF
c      ENDIF
      
      IF (itiming) THEN
        CALL getused(tass2)
        tass = tass + (tass2 - tass1)
      ENDIF

c
c--Compute forces on EACH particle
c
      IF (itiming) CALL getused(tforce1)

      CALL forcei(nlst_in,nlst_end,llist,dt,itime,npart,
     &     xyzmh,vxyzu,dvxyzu,dha,dumrho,pr,vsound,alphaMM,ekcle,
     &     dedxyz,Bxyz,dBevolxyz,Bevolxyz)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

      RETURN
      END
