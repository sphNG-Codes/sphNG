      PROGRAM minimise
c************************************************************
c                                                           *
c  This routine reads a dump in old format and outputs NG   *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/gtdble'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener1'
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

      CHARACTER*1 ians1, ians2
      CHARACTER*20 filename
      CHARACTER*100 fileident
      INTEGER*4 int1, int2
      INTEGER*8 number8
      DIMENSION nums(8)

      DATA icall/2/
c
c--Read
c
      WRITE (*,*) 'Enter filename'
      READ (*,99001) filename
99001 FORMAT(A20)

      WRITE (*,*) 'Original dump file or semi-new (o/s)?'
      READ (*,99002) ians1
99002 FORMAT(A1)

      WRITE (*,*) 'Full or small NG dump (f/s)?'
      READ (*,99002) ians2

      OPEN (11,FILE=filename,FORM='unformatted')

      IF (ians1.EQ.'o') THEN
         READ (11, END=100) udist, umass, utime,
     &        npart, n1, n2, gt, gamma, rhozero, RK2,
     &        (xyzmh(5,i), i=1, npart), escap, tkin, tgrav, tterm,
     &        (xyzmh(1,i), i=1, npart), (xyzmh(2,i), i=1, npart),
     &        (xyzmh(3,i), i=1, npart), (vxyzu(1,i), i=1, npart),
     &        (vxyzu(2,i), i=1, npart), (vxyzu(3,i), i=1, npart),
     &        (vxyzu(4,i), i=1, npart), (xyzmh(4,i), i=1, npart),
     &        (rho(i), i=1, npart), (dgrav(i), i=1, npart),
     &        dtmax, (isteps(i), i=1, npart)
     &        ,(iphase(i), i=1, npart),
     &        nptmass, (listpm(i), i=1, nptmass),
     &        (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass),
     &        (spinz(i),i=1,nptmass)
     &        ,(angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass),
     &        (angaddz(i),i=1,nptmass),
     &        anglostx, anglosty, anglostz,
     &        nreassign, naccrete, nkill, specang, ptmassin,
     &        (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass),
     &        (spinadz(i),i=1,nptmass)
c     &     ,(alphaMM(i), i=1, npart)

      ELSE

         READ (11, END=100) udist, umass, utime,
     &        npart, n1, n2, gt, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm
         DO j = 1, 5
            READ (11, END=100) (xyzmh(j,i), i=1, npart)
         END DO
         DO j = 1, 4
            READ (11, END=100) (vxyzu(j,i), i=1, npart)
         END DO
         READ (11, END=100) (rho(i), i=1, npart)
         READ (11, END=100) (dgrav(i), i=1, npart)
         READ (11, END=100) dtmax, (isteps(i), i=1, npart)
         READ (11, END=100) (iphase(i), i=1, npart)
         READ (11, END=100) nptmass, (listpm(i), i=1, nptmass),
     &        (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass),
     &        (spinz(i),i=1,nptmass)
     &        ,(angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass),
     &        (angaddz(i),i=1,nptmass),
     &        anglostx, anglosty, anglostz,
     &        nreassign, naccrete, nkill, specang, ptmassin,
     &        (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass),
     &        (spinadz(i),i=1,nptmass)
c     &     ,(alphaMM(i), i=1, npart)
      ENDIF
      CLOSE(11)

      OPEN (12,FILE=filename,FORM='unformatted')

      IF (ians2.EQ.'f') THEN
c
c--Write full dump file
c----------------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
      i1 = int1
      r1 = real(int2)
c
c--Write output file
c
      WRITE (12, ERR=100) int1,r1,int2,i1,int1
      fileident = 'FHydro1'
      WRITE (12, ERR=100) fileident
c
c--Single values
c
c--Default int
      number = 6
      WRITE (12, ERR=100) number
      WRITE (12, ERR=100) npart,n1,n2,nreassign,naccrete,nkill
c--int*1, int*2, int*4, int*8
      number = 0
      DO i = 1, 4
         WRITE (12, ERR=100) number
      END DO
c--Default real
      number = 14
      WRITE (12, ERR=100) number
      WRITE (12, ERR=100) gt, dtmax, gamma, rhozero, RK2,
     &     escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &     specang, ptmassin
c--real*4
      number = 0
      WRITE (12, ERR=100) number
c--real*8
      number = 3
      WRITE (12, ERR=100) number
      WRITE (12, ERR=100) udist, umass, utime
c
c--Arrays
c
c--Number of array lengths
c
      number = 2
      WRITE (12, ERR=100) number
c
c--Array length 1
c
      number8 = npart
      nums(1) = 1
      nums(2) = 1
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 9
      nums(7) = 2
      nums(8) = 0
      WRITE (12, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 2
c
      number8 = nptmass
      nums(1) = 1
      nums(2) = 0
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 9
      nums(7) = 0
      nums(8) = 0
      WRITE (12, ERR=100) number8, (nums(i), i=1,8)
c
c--Arrays of length 1
c
c--Default int
      WRITE (12, ERR=100) (isteps(i), i=1, npart)
c--int*1
      WRITE (12, ERR=100) (iphase(i), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      DO j = 1, 5
         WRITE (12, ERR=100) (xyzmh(j,i), i=1, npart)
      END DO
      DO j = 1, 4
         WRITE (12, ERR=100) (vxyzu(j,i), i=1, npart)
      END DO      
c--real*4
      WRITE (12, ERR=100) (rho(i), i=1, npart)
      WRITE (12, ERR=100) (dgrav(i), i=1, npart)      
c     WRITE (12, ERR=100) (alphaMM(i), i=1, npart)
c--real*8

c
c--Arrays of length 2
c
c--Default int
      WRITE (12, ERR=100) (listpm(i), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      WRITE (12, ERR=100) (spinx(i),i=1,nptmass)
      WRITE (12, ERR=100) (spiny(i),i=1,nptmass)
      WRITE (12, ERR=100) (spinz(i),i=1,nptmass)
      WRITE (12, ERR=100) (angaddx(i),i=1,nptmass)
      WRITE (12, ERR=100) (angaddy(i),i=1,nptmass)
      WRITE (12, ERR=100) (angaddz(i),i=1,nptmass)
      WRITE (12, ERR=100) (spinadx(i),i=1,nptmass)
      WRITE (12, ERR=100) (spinady(i),i=1,nptmass)
      WRITE (12, ERR=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

c
c--End writing of full dump file
c-------------------------------
c
      ELSE
c
c--Write small dump file
c----------------------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
      i1 = int1
      r1 = i2
c
c--Write ouput file
c
      WRITE (12, ERR=100) int1,i1,int2,r1,int1
      fileident = 'SHydro1'
      WRITE (12, ERR=100) fileident
c
c--Single values
c
c--Default int
      number = 6
      WRITE (12, ERR=100) number
      WRITE (12, ERR=100) npart,n1,n2,nreassign,naccrete,nkill
c--int*1, int*2, int*4, int*8
      number = 0
      DO i = 1, 4
         WRITE (12, ERR=100) number
      END DO
c--Default real
      number = 15
      DO i = 1, npart
         IF (iphase(i).EQ.0) THEN
            pmassinitial = xyzmh(4,i)
            GOTO 40
         ENDIF
      END DO
 40      WRITE (12, ERR=100) number
      WRITE (12, ERR=100) gt, dtmax, gamma, rhozero, RK2,
     &     escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &     specang, ptmassin, pmassinitial
c--real*4
      number = 0
      WRITE (12, ERR=100) number
c--real*8
      number = 3
      WRITE (12, ERR=100) number
      WRITE (12, ERR=100) udist, umass, utime
c
c--Arrays
c
c--Number of array lengths
c
      number = 2
      WRITE (12, ERR=100) number
c
c--Array length 1
c
      number8 = npart
      nums(1) = 0
      nums(2) = 1
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 5
      nums(7) = 2
      nums(8) = 0
      WRITE (12, ERR=100) number8, (nums(i), i=1,8)
c
c--Array length 2
c
      number8 = nptmass
      nums(1) = 1
      nums(2) = 0
      nums(3) = 0
      nums(4) = 0
      nums(5) = 0
      nums(6) = 1
      nums(7) = 0
      nums(8) = 0
      WRITE (12, ERR=100) number8, (nums(i), i=1,8)
c
c--Arrays of length 1
c
c--Default int
c      WRITE (12, ERR=100) (isteps(i), i=1, npart)
c--int*1
      WRITE (12, ERR=100) (iphase(i), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      DO j = 1, 3
         WRITE (12, ERR=100) (xyzmh(j,i), i=1, npart)
      END DO
c      DO j = 1, 4
c         WRITE (12, ERR=100) (vxyzu(j,i), i=1, npart)
c      END DO      
c--real*4
      WRITE (12, ERR=100) (rho(i), i=1, npart)
      DO i = 1, npart
         dq(i) = xyzmh(5,i)
      END DO
      WRITE (12, ERR=100) (dq(i), i=1, npart)      
c      WRITE (12, ERR=100) (dgrav(i), i=1, npart)      
c     WRITE (12, ERR=100) (alphaMM(i), i=1, npart)
c--real*8

c
c--Arrays of length 2
c
c--Default int
      WRITE (12, ERR=100) (listpm(i), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
      WRITE (12, ERR=100) (xyzmh(4,listpm(i)), i=1,nptmass)
c      WRITE (12, ERR=100) (spinx(i),i=1,nptmass)
c      WRITE (12, ERR=100) (spiny(i),i=1,nptmass)
c      WRITE (12, ERR=100) (spinz(i),i=1,nptmass)
c      WRITE (12, ERR=100) (angaddx(i),i=1,nptmass)
c      WRITE (12, ERR=100) (angaddy(i),i=1,nptmass)
c      WRITE (12, ERR=100) (angaddz(i),i=1,nptmass)
c      WRITE (12, ERR=100) (spinadx(i),i=1,nptmass)
c      WRITE (12, ERR=100) (spinady(i),i=1,nptmass)
c      WRITE (12, ERR=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

c
c--End writing of small dump file
c--------------------------------
c
      ENDIF

      CLOSE(12)
      RETURN

 100  WRITE (*,*) 'ERROR'
      
      RETURN
      END
