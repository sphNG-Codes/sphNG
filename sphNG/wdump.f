      SUBROUTINE wdump(idisk1)
c************************************************************
c                                                           *
c  This routine writes a dump on disk                       *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
c      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/treecom_P'

      DIMENSION itempsort(idim), tempsort(idim)
      EQUIVALENCE (itempsort, next1), (tempsort, key)

      CHARACTER*7 where

      DATA where/'wdump'/
c
c--Write
c
      irec = irec + 1
      iresort = iresort + 1

      WRITE (idisk1, ERR=100) udist, umass, utime,
     &     npart, n1, n2, gt, gamma, rhozero, RK2,
     &     escap, tkin, tgrav, tterm
      DO j = 1, 5
         WRITE (idisk1, ERR=100) (xyzmh(j,isort(i)), i=1, npart)
      END DO
      DO j = 1, 4
         WRITE (idisk1, ERR=100) (vxyzu(j,isort(i)), i=1, npart)
      END DO
      WRITE (idisk1, ERR=100) (rho(isort(i)), i=1, npart)
      WRITE (idisk1, ERR=100) (dgrav(isort(i)), i=1, npart)
      WRITE (idisk1, ERR=100) dtmax, (isteps(isort(i)), i=1, npart)
      WRITE (idisk1, ERR=100) (iphase(isort(i)), i=1, npart)
      WRITE (idisk1, ERR=100) nptmass, (iorig(listpm(i)), i=1,nptmass),
     &     (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass),
     &     (spinz(i),i=1,nptmass)
     &     ,(angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass),
     &     (angaddz(i),i=1,nptmass),
     &     anglostx, anglosty, anglostz,
     &     nreassign, naccrete, nkill, specang, ptmassin,
     &     (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass),
     &     (spinadz(i),i=1,nptmass)
c     &     ,(alphaMM(isort(i)), i=1, npart)

c
c--Sort particles to ensure most efficient running.  Note that this 
c     should not be visible to the outside observer.  In other words,
c     an array must be kept of original list of particles and this
c     must be used to index *ANY* value from an array which is written 
c     to the outside.  This requires modification to almost every output
c     line in the code.  Done 21 Nov 2000.
c
      IF (gt.EQ.0. .OR. iresort.LT.5) GOTO 777

      iresort = 0

      xminimum = 1.0E+30
      xmaximum = -1.0E+30
      DO i = 1, npart
         xmaximum = MAX(xmaximum, xyzmh(1,i))
         xminimum = MIN(xminimum, xyzmh(1,i))
      END DO
      xrange = xmaximum-xminimum
      istepmin = imax
      istepmax = 0
      DO i = 1, npart
         llist(i) = i
         IF (iphase(i).EQ.-1) THEN
            tempsort(i) = LOG(REAL(imax))/LOG(2.0)
            istepmax = imax
         ELSE
            tempsort(i) = LOG(REAL(isteps(i)))/LOG(2.0) + 
     &           (xyzmh(1,i)-xminimum)/xrange
            istepmin = MIN(istepmin, isteps(i))
            istepmax = MAX(istepmax, isteps(i))
         ENDIF
      END DO
c
c--Sort particles based on their individual timesteps and x
c
      CALL indexx(npart, llist, tempsort, itempsort)
c
c--Move into llist.  This now stores 'iorignew'
c
      DO i = 1, npart
         llist(i) = itempsort(i)
      END DO
c
c--Now re-order arrays
c
      DO i = 1, npart
         tempsort(i) = xyzmh(5,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(5,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(1,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(2,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(3,llist(i))
       END DO
      DO i = 1, npart
         xyzmh(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
        tempsort(i) = vxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         vxyzu(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(4,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = rho(llist(i))
      END DO
      DO i = 1, npart
         rho(i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dgrav(llist(i))
      END DO
      DO i = 1, npart
         dgrav(i) = tempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = isteps(llist(i))
      END DO
      DO i = 1, npart
         isteps(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = iphase(llist(i))
      END DO
      DO i = 1, npart
         iphase(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = it0(llist(i))
      END DO
      DO i = 1, npart
         it0(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = it1(llist(i))
      END DO
      DO i = 1, npart
         it1(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = it2(llist(i))
      END DO
      DO i = 1, npart
         it2(i) = itempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1vxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         f1vxyzu(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f1ha(1,llist(i))
      END DO
      DO i = 1, npart
         f1ha(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumrho(llist(i))
      END DO
      DO i = 1, npart
         dumrho(i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = alphaMM(llist(i))
      END DO
      DO i = 1, npart
         alphaMM(i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = ddv(llist(i))
      END DO
      DO i = 1, npart
         ddv(i) = tempsort(i)
      END DO
c
c--iorig
c
      DO i = 1, npart
         itempsort(i) = iorig(llist(i))
      END DO
      DO i = 1, npart
         iorig(i) = itempsort(i)
      END DO
c
c--isortnew
c
      DO i = 1, npart
         itempsort(llist(i)) = i
      END DO
      DO i = 1, nptmass
         listpm(i) = itempsort(listpm(i))
      END DO
      DO i = 1, nlstacc
         listacc(i) = itempsort(listacc(i))
      END DO
      DO i = 1, nbinmax
         DO j = 1, nlstbins(i)
            listbins(j,i) = itempsort(listbins(j,i))
         END DO
      END DO
      DO i = 1, npart
         isort(i) = itempsort(isort(i))
      END DO
c
c--Must rebuild the ghosts and make the tree again after sorting
c
      nghost = 0
      IF (ibound.EQ.1) CALL ghostp1(npart,xyzmh,vxyzu)
      IF (ibound.EQ.2) CALL ghostp2(npart,xyzmh,vxyzu)
      IF (ibound.EQ.3 .OR. ibound.EQ.8 .OR. ibound.GE.90)
     &     CALL ghostp3(npart,xyzmh,vxyzu)
      IF (ibound.EQ.11) CALL ghostp11(npart,xyzmh,vxyzu)

      ntot = npart + nghost

      DO i = 1, ntot
         DO j = 1, 5
            dumxyzmh(j,i) = xyzmh(j,i)
         END DO
         DO j = 1, 4
            dumvxyzu(j,i) = vxyzu(j,i)
         END DO
      END DO

      CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)
 777  CONTINUE
c
c--Zero torques
c
c      DO i = 1, npart
c         torqt(i) = 0.
c         torqg(i) = 0.
c         torqp(i) = 0.
c         torqv(i) = 0.
c         torqc(i) = 0.
c      END DO

      ENDFILE idisk1
      BACKSPACE idisk1

      ipos = irec

      RETURN
c
c--An error as occured while writing
c
 100  CALL error(where, 1)

      RETURN
      END
