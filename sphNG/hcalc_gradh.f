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
      INCLUDE 'COMMONS/polyk2'
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
      acc = 0.9

      neimean = (neimax + neimin) / 2
      neirange = (neimax - neimin) / 4
      neisup = neimean + neirange
      neiinf = neimean - neirange

      ikount = 0

      WRITE(*,*) neimax, neimean, neimin, neisup, neiinf, neirange

      nlst = npart
      DO i = 1, npart
         iscurrent(i) = .TRUE.
         llist(i) = i
      END DO

      IF (ibound.EQ.1)
     &      CALL ghostp1(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
      IF (ibound.EQ.2)
     &      CALL ghostp2(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
      IF (ibound.EQ.3 .OR. ibound.EQ.8 .OR. ibound/10.EQ.9)
     &      CALL ghostp3(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
      IF (ibound.EQ.100)
     &     CALL ghostp100(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
      IF (ibound.EQ.11)
     &      CALL ghostp11(npart, xyzmh, vxyzu, ekcle, Bevolxyz)
      IF (ibound.EQ.0) nghost = 0

      ntot = npart + nghost

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

      icount = icount + 1
      ikount = ikount + 1
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
      WRITE(*,*) ' Got neighbours from tree'
c
c--calculate density
c
      WRITE(*,*) ' Calculating density and h... '
      nlst = 0
      DO i = 1, npart
         IF (iphase(i).NE.-1) THEN
            nlst = nlst + 1
            llist(nlst) = i
            iscurrent(i) = .TRUE.
         ENDIF
      END DO    
      itime = 0
      nlst_in = 1
      nlst_end = nlst
      hfact = 1.2
      
      DO ipart=1,npart
         dumxyzmh(5,ipart) = hfact*(dumxyzmh(4,ipart)/rhozero)**third
      ENDDO      

c      CALL iterate_density(npart,dumxyzmh,vxyzu,f1ha,
c     &                     nlst_in,nlst_end,llist,itime)
      CALL derivi(dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,npart,ntot,
     &            ireal,alphaMM,ekcle,Bevolxyz,f1Bxyz)

      DO ipart=1,npart
         xyzmh(5,ipart) = dumxyzmh(5,ipart)
      ENDDO
c
c--get h max and min and print to screen
c
      hhmin = 1.e12
      hhmax = -1.e12
      hhav = 0.
      DO ipart=1,npart
         IF (xyzmh(5,ipart).LT.hhmin) hhmin = dumxyzmh(5,ipart)
         IF (xyzmh(5,ipart).GT.hhmax) hhmax = dumxyzmh(5,ipart)
         hhav = hhav + xyzmh(5,ipart)
      ENDDO
      hhav = hhav/real(npart)
      WRITE(*,*) ' min h = ',hhmin,' max h = ',hhmax,' av h =',hhav

      DO i = 1, npart
         iscurrent(i) = .FALSE.
      END DO

      RETURN
      END
