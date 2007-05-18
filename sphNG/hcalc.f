      SUBROUTINE hcalc
c************************************************************
c                                                           *
c  This subroutine computes the derivative of the smoothing *
c     length.                                               *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine hcalc')

      third = 1./3.
      neimin = 30
      neimax = 70
c      neimin = 80
c      neimax = 120
      acc = 0.5
      idumg = 0
      idummy = 0

      neimean = (neimax + neimin) / 2
      neirange = (neimax - neimin) / 4
      neisup = neimean + neirange
      neiinf = neimean - neirange

      ikount = 0

      WRITE(*,*) neimax, neimean, neimin, neisup, neiinf, neirange

   5  CONTINUE

      nlst = npart
      DO i = 1, npart
         iscurrent(i) = .TRUE.
         llist(i) = i
      END DO

      CALL ghostp(ntot, npart, xyzmh, vxyzu, ekcle, Bevolxyz)

      IF (igrape.EQ.0) THEN
         DO i = 1, ntot
            DO j = 1, 5
               dumxyzmh(j,i) = xyzmh(j,i)
            END DO
         END DO
         WRITE(*,*) ' Making tree'
         CALL insulate(1, ntot, npart, dumxyzmh, f1vxyzu)
      ENDIF

      icount = 0

   10 CONTINUE

      icount = icount + 1
      ikount = ikount + 1
      ioutinf1 = 0
      ioutsup1 = 0
      ioutinf2 = 0
      ioutsup2 = 0
      iagain = 0
      WRITE(*,*) ' Calculating neighbour changes'
      icall = 1
c
c--Get neighbours
c
      IF (igrape.EQ.0) THEN
         CALL insulate(5, ntot, npart, dumxyzmh, f1vxyzu)
      ELSEIF (igrape.EQ.1) THEN
         CALL insulate(4, ntot, npart, dumxyzmh, f1vxyzu)
      ENDIF

      DO ipart = 1, npart

         IF(iphase(ipart).GE.1) GOTO 12

         numneigh = nneigh(ipart)
c         print *, ipart, numneigh

         IF (FLOAT(numneigh)/FLOAT(neimean).LT.0.1)
     &                           numneigh = neimean/10
         IF (numneigh.LT.neiinf) THEN
            IF (numneigh.LT.neimin) THEN
               iagain = iagain + 1
               ioutinf1 = ioutinf1 + 1
               xyzmh(5,ipart) = xyzmh(5,ipart) *
     &              (FLOAT(neimean) / FLOAT(numneigh + 1))**third
            ELSE 
               ioutinf2 = ioutinf2 + 1
               xyzmh(5,ipart) = xyzmh(5,ipart) *
     &              (FLOAT(neimean) / FLOAT(numneigh + 1))**third
            ENDIF
         ELSEIF (numneigh.GT.neisup) THEN
            IF (numneigh.GT.neimax) THEN
               iagain = iagain + 1
               ioutsup1 = ioutsup1 + 1
               xyzmh(5,ipart) = xyzmh(5,ipart) *
     &              (FLOAT(neimean) / FLOAT(numneigh + 1))**third
            ELSE
               ioutsup2 = ioutsup2 + 1
               xyzmh(5,ipart) = xyzmh(5,ipart) *
     &              (FLOAT(neimean) / FLOAT(numneigh + 1))**third
            ENDIF
         ENDIF
 12      CONTINUE
      END DO

      IF (ioutinf1.NE.0)
     &     WRITE (*,*) 'h too small ', ioutinf1, ' times'
      IF (ioutsup1.NE.0)
     &     WRITE (*,*) 'h too big ', ioutsup1, ' times'
      IF (ioutinf2.NE.0)
     &     WRITE (*,*) 'h near lower limit ', ioutinf2, ' times'
      IF (ioutsup2.NE.0)
     &     WRITE (*,*) 'h near upper limit ', ioutsup2, ' times'

      IF (ikount.GT.10) GOTO 15
ccc      IF (iagain.GT.0 .AND. icount.GT.2) GOTO 5
ccc      IF (iagain.GT.0) GOTO 10
      IF (iagain.GT.0) GOTO 5
c
c--get density
c
      CALL densityi(npart,xyzmh,vxyzu,ekcle,
     &        1,npart,llist,itime)


 15   DO i = 1, npart
         iscurrent(i) = .FALSE.
      END DO

      RETURN
      END
