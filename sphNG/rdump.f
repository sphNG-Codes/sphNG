      SUBROUTINE rdump(idisk1, ichkl)
c************************************************************
c                                                           *
c  This routine reads a dump into memory                    *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      REAL*8 umassi, udisti, utimei

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/gtdble'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
c      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/treecom_P'

      DIMENSION itempsort(idim), tempsort(idim)
      EQUIVALENCE (itempsort, next1), (tempsort, key)

      CHARACTER*7 where
      CHARACTER*100 fileident
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      DIMENSION nums(8)

      DATA icall/2/
      DATA where/'rdump'/
c
c--Read
c
      IF (itrace.EQ.'all') WRITE (*, 99001)
99001 FORMAT (' entry subroutine rdump')

      ichkl = 0
c
c--Dump file
c-------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Write ouput file
c
      READ (idisk1, END=100) int1i,r1i,int2i,i1i,int3i
      IF (int1i.NE.int1) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         CALL quit
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default integer size wrong'
         CALL quit
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default real size wrong'
         CALL quit
      ENDIF
      READ (idisk1, END=100) fileident
c
c--Single values
c
c--Default int
      READ (idisk1, END=100) number
      IF (number.NE.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump'
         CALL quit
      ENDIF
      READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (idisk1, END=100) number
      END DO
c--Default real
      READ (idisk1, END=100) number
      IF (number.NE.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump'
         CALL quit
      ENDIF
      READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &     escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &     specang, ptmassin
c--real*4
      READ (idisk1, END=100) number
c--real*8
      READ (idisk1, END=100) number
      IF (number.NE.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump'
         CALL quit
      ENDIF
      READ (idisk1, END=100) udisti, umassi, utimei
c
c--Arrays
c
c--Number of array lengths
c
      READ (idisk1, END=100) number
      IF (number.NE.2) THEN
         WRITE (*,*) 'ERROR 7 in rdump'
         CALL quit
      ENDIF
c
c--Read array type 1 header
c
      READ (idisk1, END=100) number8, (nums(i), i=1,8)
      IF (number8.NE.npart) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong'
         CALL quit
      ENDIF
      npart = number8
c
c--Read array type 2 header
c
      READ (idisk1, END=100) number8, (nums(i), i=1,8)
      nptmass = number8
c
c--Read array type 1 arrays
c
c--Default int
      READ (idisk1, END=100) (isteps(i), i=1, npart)
c--int*1
      READ (idisk1, END=100) (iphase(i), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      DO j = 1, 5
         READ (idisk1, END=100) (xyzmh(j,i), i=1, npart)
      END DO
      DO j = 1, 4
         READ (idisk1, END=100) (vxyzu(j,i), i=1, npart)
      END DO      
c--real*4
      READ (idisk1, END=100) (rho(i), i=1, npart)
      READ (idisk1, END=100) (dgrav(i), i=1, npart)
c     READ (idisk1, END=100) (alphaMM(i), i=1, npart)
c--real*8

c
c--Read array type 2 arrays
c
c--Default int
      READ (idisk1, END=100) (listpm(i), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      READ (idisk1, END=100) (spinx(i),i=1,nptmass)
      READ (idisk1, END=100) (spiny(i),i=1,nptmass)
      READ (idisk1, END=100) (spinz(i),i=1,nptmass)
      READ (idisk1, END=100) (angaddx(i),i=1,nptmass)
      READ (idisk1, END=100) (angaddy(i),i=1,nptmass)
      READ (idisk1, END=100) (angaddz(i),i=1,nptmass)
      READ (idisk1, END=100) (spinadx(i),i=1,nptmass)
      READ (idisk1, END=100) (spinady(i),i=1,nptmass)
      READ (idisk1, END=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

c
c--End reading of dump file
c--------------------------
c
      gtdouble = DBLE(gt)
c
c--Sort particles to ensure most efficient running.  Note that this 
c     should not be visible to the outside observer.  In other words,
c     an array must be kept of original list of particles and this
c     must be used to index *ANY* value from an array which is written 
c     to the outside.  This requires modification to almost every output
c     line in the code.  Done 21 Nov 2000.
c
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
            IF (gt.EQ.0.0 .OR. isteps(i).EQ.0) THEN
               tempsort(i) = (xyzmh(1,i)-xminimum)/xrange
            ELSE
               tempsort(i) = LOG(REAL(isteps(i)))/LOG(2.0) +
     &           (xyzmh(1,i)-xminimum)/xrange
            ENDIF
            istepmin = MIN(istepmin, isteps(i))
            istepmax = MAX(istepmax, isteps(i))
         ENDIF
      END DO
c
c--Sort particles based on their individual timesteps and x
c

      CALL indexx(npart, llist, tempsort, iorig)

      DO i = 1, npart
         isort(iorig(i)) = i
         tempsort(i) = xyzmh(5,iorig(i))
      END DO
      DO i = 1, npart
         xyzmh(5,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(1,iorig(i))
      END DO
      DO i = 1, npart
         xyzmh(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(2,iorig(i))
      END DO
      DO i = 1, npart
         xyzmh(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(3,iorig(i))
      END DO
      DO i = 1, npart
         xyzmh(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(1,iorig(i))
      END DO
      DO i = 1, npart
         vxyzu(1,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(2,iorig(i))
      END DO
      DO i = 1, npart
         vxyzu(2,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(3,iorig(i))
      END DO
      DO i = 1, npart
         vxyzu(3,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = vxyzu(4,iorig(i))
      END DO
      DO i = 1, npart
         vxyzu(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = xyzmh(4,iorig(i))
      END DO
      DO i = 1, npart
         xyzmh(4,i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = rho(iorig(i))
      END DO
      DO i = 1, npart
         rho(i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dgrav(iorig(i))
      END DO
      DO i = 1, npart
         dgrav(i) = tempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = isteps(iorig(i))
      END DO
      DO i = 1, npart
         isteps(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = iphase(iorig(i))
      END DO
      DO i = 1, npart
         iphase(i) = itempsort(i)
      END DO
      DO i = 1, nptmass
         listpm(i) = isort(listpm(i))
      END DO
c
c--Zero torques
c
c      DO i = 1, idim
c         torqt(i) = 0.
c         torqg(i) = 0.
c         torqp(i) = 0.
c         torqv(i) = 0.
c         torqc(i) = 0.
c      END DO
c
c--Check units in file the same as in the code!
c
      IF (udisti.LT.0.99999*udist .OR. udisti.GT.1.00001*udist) THEN
         CALL error(where,1)
      ELSEIF (umassi.LT.0.99999*umass .OR.umassi.GT.1.00001*umass) THEN
         CALL error(where,2)
      ENDIF
      IF (npart.GT.idim) THEN
         CALL error(where,3)
      ENDIF
c
c--Check that dtmax times are the same.  If not, modify isteps(i) as in mesop.f
c
ccc      GOTO 50

      IF (gt.NE.0.0 .AND. 
     &     (dtmaxdp.LT.0.9999*dtmax .OR. dtmaxdp.GT.1.0001*dtmax)) THEN
         ipower = INT(LOG10(dtmax/dtmaxdp)/LOG10(2.0))

         ifactor = 2**ABS(ipower)
         imaxnew = imaxstep/ifactor
         iminnew = 2*ifactor

         IF (ipower.LT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  isteps(j) = MIN(isteps(j), imaxnew)
                  isteps(j) = isteps(j)*ifactor
               ENDIF
            END DO
         ELSEIF (ipower.GT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  IF (isteps(j)/ifactor .LT. 2) CALL error(where, 4)
                  isteps(j) = isteps(j)/ifactor
               ENDIF
            END DO
         ENDIF
      ENDIF
c
c--Change reference frame
c
 50   IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2.AND.gt.NE.0.0)) THEN
c      IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2)) THEN
         CALL chanref(icall)
      ELSEIF (ifcor.GT.2) THEN
         ifcor = ifcor - 2
      ENDIF

      IF (itrace.EQ.'all') WRITE (*, 99002)
99002 FORMAT (' exit subroutine rdump')
      RETURN

 100  ichkl = 1

      RETURN
      END
