      SUBROUTINE wdump_wrapper(idisk1)
c************************************************************
c                                                           *
c  This routine writes a dump on disk                       *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener1'
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
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/savernd'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/f2'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/perform'

      DIMENSION itempsort(idim)
      EQUIVALENCE (itempsort, next1)

      CHARACTER*7 where
      CHARACTER*100 fileident
      INTEGER*4 int1, int2
      INTEGER*8 number8
      DIMENSION nums(8)
c
c--Standard numbers
c
      PARAMETER (int1=690706,int2=780806)

      DATA where/'wdump_w'/
c
c--Do the actual write
c
      CALL wdump(idisk1)
c
c--Sort particles to ensure most efficient running.  Note that this 
c     should not be visible to the outside observer.  In other words,
c     an array must be kept of original list of particles and this
c     must be used to index *ANY* value from an array which is written 
c     to the outside.  This requires modification to almost every output
c     line in the code.  Done 21 Nov 2000.
c
      IF (itiming) CALL getused(tsort1)

      iresort = iresort + 1

c      IF (gt.EQ.0. .OR. iresort.LT.2) GOTO 777
      IF (gt.EQ.0.) GOTO 777
      GOTO 777

      iresort = 0

      IF (nptmass.GT.0) THEN
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
     &              (xyzmh(1,i)-xminimum)/xrange
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
      ELSE
c
c--Sort particles based on a walk through the tree
c
         CALL treelist(npart, llist, dumxyzmh)
      ENDIF
c
c--Now re-order arrays
c
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
         tempsort(i) = xyzmh(4,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(5,llist(i))
      END DO
      DO i = 1, npart
         xyzmh(5,i) = tempsort(i)
      END DO
c
c--vxyzu
c
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
c
c--dumxyzmh
c
      DO i = 1, npart
         tempsort(i) = dumxyzmh(1,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(2,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(3,llist(i))
       END DO
      DO i = 1, npart
         dumxyzmh(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(4,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumxyzmh(5,llist(i))
      END DO
      DO i = 1, npart
         dumxyzmh(5,i) = tempsort(i)
      END DO
c
c--dumvxyzu
c
      DO i = 1, npart
        tempsort(i) = dumvxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumvxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumvxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dumvxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         dumvxyzu(4,i) = tempsort(i)
      END DO
c
c--rho
c
      DO i = 1, npart
         tempsort(i) = rho(llist(i))
      END DO
      DO i = 1, npart
         rho(i) = tempsort(i)
      END DO
c
c--divv
c
      DO i = 1, npart
         tempsort(i) = divv(llist(i))
      END DO
      DO i = 1, npart
         divv(i) = tempsort(i)
      END DO
c
c--curlv
c
      DO i = 1, npart
         tempsort(i) = curlv(llist(i))
      END DO
      DO i = 1, npart
         curlv(i) = tempsort(i)
      END DO
c
c--pr
c
      DO i = 1, npart
         tempsort(i) = pr(llist(i))
      END DO
      DO i = 1, npart
         pr(i) = tempsort(i)
      END DO
c
c--vsound
c
      DO i = 1, npart
         tempsort(i) = vsound(llist(i))
      END DO
      DO i = 1, npart
         vsound(i) = tempsort(i)
      END DO
c
c--poten
c
      DO i = 1, npart
         tempsort(i) = poten(llist(i))
      END DO
      DO i = 1, npart
         poten(i) = tempsort(i)
      END DO
c
c--dgrav
c
      DO i = 1, npart
         tempsort(i) = dgrav(llist(i))
      END DO
      DO i = 1, npart
         dgrav(i) = tempsort(i)
      END DO
c
c--gradhs
c
      DO k = 1, 2
         DO i = 1, npart
            tempsort(i) = gradhs(k,llist(i))
         END DO
         DO i = 1, npart
            gradhs(k,i) = tempsort(i)
         END DO
      END DO
c
c--isteps
c
      DO i = 1, npart
         itempsort(i) = isteps(llist(i))
      END DO
      DO i = 1, npart
         isteps(i) = itempsort(i)
      END DO
c
c--iphase
c
      DO i = 1, npart
         itempsort(i) = iphase(llist(i))
      END DO
      DO i = 1, npart
         iphase(i) = itempsort(i)
      END DO
c
c--it0
c
      DO i = 1, npart
         itempsort(i) = it0(llist(i))
      END DO
      DO i = 1, npart
         it0(i) = itempsort(i)
      END DO
c
c--it1
c
      DO i = 1, npart
         itempsort(i) = it1(llist(i))
      END DO
      DO i = 1, npart
         it1(i) = itempsort(i)
      END DO
c
c--it2
c
      DO i = 1, npart
         itempsort(i) = it2(llist(i))
      END DO
      DO i = 1, npart
         it2(i) = itempsort(i)
      END DO
c
c--f1vxyzu
c
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
c
c--f2vxyzu
c
      DO i = 1, npart
         tempsort(i) = f2vxyzu(1,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f2vxyzu(2,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f2vxyzu(3,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = f2vxyzu(4,llist(i))
      END DO
      DO i = 1, npart
         f2vxyzu(4,i) = tempsort(i)
      END DO
c
c--f1ha
c
      DO k = 1, 1+isizealphaMM
         DO i = 1, npart
            tempsort(i) = f1ha(k,llist(i))
         END DO
         DO i = 1, npart
            f1ha(k,i) = tempsort(i)
         END DO
c
c--f2ha
c
         DO i = 1, npart
            tempsort(i) = f2ha(k,llist(i))
         END DO
         DO i = 1, npart
            f2ha(k,i) = tempsort(i)
         END DO
      END DO
c
c--dumrho
c
      DO i = 1, npart
         tempsort(i) = dumrho(llist(i))
      END DO
      DO i = 1, npart
         dumrho(i) = tempsort(i)
      END DO
c
c--alphaMM
c
      DO k = 1, isizealphaMM
         DO i = 1, npart
            tempsort(i) = alphaMM(k,llist(i))
         END DO
         DO i = 1, npart
            alphaMM(k,i) = tempsort(i)
         END DO
c
c--dumalpha
c
         DO i = 1, npart
            tempsort(i) = dumalpha(k,llist(i))
         END DO
         DO i = 1, npart
            dumalpha(k,i) = tempsort(i)
         END DO
      END DO
c
c--ddv
c
      DO i = 1, npart
         tempsort(i) = ddv(llist(i))
      END DO
      DO i = 1, npart
         ddv(i) = tempsort(i)
      END DO
c
c--chemistry
c
      IF (idim_h2.EQ.idim) THEN
         DO i = 1, npart
            tempsort(i) = h2ratio(llist(i))
         END DO

         DO i = 1, npart
            h2ratio(i) = tempsort(i)
         END DO

         DO i = 1, npart
            tempsort(i) = abHIq(llist(i))
         END DO

         DO i = 1, npart
            abHIq(i) = tempsort(i)
         END DO

         DO i = 1, npart
            tempsort(i) = abhpq(llist(i))
         END DO

         DO i = 1, npart
            abhpq(i) = tempsort(i)
         END DO

         DO i = 1, npart
            tempsort(i) = abeq(llist(i))
         END DO

         DO i = 1, npart
            abeq(i) = tempsort(i)
         END DO
 
         DO i = 1, npart
            tempsort(i) = abco(llist(i))
         END DO

         DO i = 1, npart
            abco(i) = tempsort(i)
         END DO
      ENDIF
c
c--ekcle, dumekcle
c
      IF (encal.EQ.'r') THEN
         DO k = 1, 5
            DO i = 1, npart
               tempsort(i) = ekcle(k,llist(i))
            END DO
            DO i = 1, npart
               ekcle(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = dumekcle(k,llist(i))
            END DO
            DO i = 1, npart
               dumekcle(k,i) = tempsort(i)
            END DO
         END DO
      ENDIF
c
c--Bevolxyz, dumBevolxyz, Bxyz, f1Bxyz, f2Bxyz
c
      IF (imhd.EQ.idim) THEN
         DO k = 1, 3
            DO i = 1, npart
               tempsort(i) = Bevolxyz(k,llist(i))
            END DO
            DO i = 1, npart
               Bevolxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = dumBevolxyz(k,llist(i))
            END DO
            DO i = 1, npart
               dumBevolxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = Bxyz(k,llist(i))
            END DO
            DO i = 1, npart
               Bxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = f1Bxyz(k,llist(i))
            END DO
            DO i = 1, npart
               f1Bxyz(k,i) = tempsort(i)
            END DO

            DO i = 1, npart
               tempsort(i) = f2Bxyz(k,llist(i))
            END DO
            DO i = 1, npart
               f2Bxyz(k,i) = tempsort(i)
            END DO
         END DO
      ENDIF
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
         listrealpm(listpm(i)) = i
      END DO
      DO i = 1, nlstacc
         listacc(i) = itempsort(listacc(i))
      END DO
      DO i = 1, nlistinactive
         listinactive(i) = itempsort(listinactive(i))
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
      CALL ghostp(ntot, npart, xyzmh, vxyzu, ekcle, Bevolxyz)

      DO i = 1, ntot
         DO j = 1, 5
            dumxyzmh(j,i) = xyzmh(j,i)
         END DO
         DO j = 1, 4
            dumvxyzu(j,i) = vxyzu(j,i)
         END DO
      END DO

      IF (itiming) CALL getused(tsortmtree1)

      CALL insulate(1,ntot,npart,dumxyzmh,f1vxyzu)

      IF (itiming) THEN
         CALL getused(tsortmtree2)
         tsortmtree = tsortmtree + (tsortmtree2 - tsortmtree1)
      ENDIF
 777  CONTINUE

      IF (itiming) THEN
         CALL getused(tsort2)
         tsort = tsort + (tsort2 - tsort1)
      ENDIF
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

      ipos = irec

      RETURN
c
c--An error as occured while writing
c
 100  CALL error(where, 1)

      RETURN
      END
