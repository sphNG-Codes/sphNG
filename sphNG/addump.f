      SUBROUTINE addump
c************************************************************
c                                                           *
c  This routine reads two existing dumps and create a third *
c  one out of it.                                           *
c                                                           *
c************************************************************

      INCLUDE 'idim' 

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/new'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'

      CHARACTER*7 where
      CHARACTER*1 iok
      CHARACTER*100 fileident
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      DIMENSION nums(8)

      DATA where/'addump'/
c
c--allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine addump')

      fmas1 = 0.
      fmas2 = 0.
      WRITE (*, 99002)
99002 FORMAT (' PARTICLE SET UP', //)
      WRITE (*, 99003)
99003 FORMAT (' self-gravity included? (y/n)')
      READ (*, 99004) iok
99004 FORMAT (A1)
      igphi = 0
      IF ( iok.EQ.'y' ) igphi = 1
      WRITE (*, 99005)
99005 FORMAT (' give the position and the velocity of the cm of body 1')
      READ (*, *) xx1, yy1, zz1m, vvx1, vvy1, vvz1
      WRITE (*, 99006)
99006 FORMAT (' give the position and the velocity of the cm of body 2')
      READ (*, *) xx2, yy2, zz2, vvx2, vvy2, vvz2
      WRITE (*, 99007)
99007 FORMAT (' do you want to rotate body 2? (y/n)')
      READ (*, 99004) iok
      rotang = 0.
      IF ( iok.EQ.'y' ) THEN
         WRITE (*, 99008)
99008    FORMAT (' give rotation angle (rad.)')
         READ (*, *) rotang
      ENDIF
c
c--read first object
c
      WRITE (*, 99009)
99009 FORMAT (' reading first file and storing data')
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Write ouput file
c
      READ (idisk1, END=300) int1i,r1i,int2i,i1i,int3i
      IF (int1i.NE.int1) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         CALL quit(0)
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default integer size wrong'
         CALL quit(0)
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default real size wrong'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) fileident
c
c--Single values
c
c--Default int
      READ (idisk1, END=300) number
      IF (number.NE.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) np1,n11,n12
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (idisk1, END=300) number
      END DO
c--Default real
      READ (idisk1, END=300) number
      IF (number.NE.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) t1, dtmax1, gamma1, rhozero1, RK21,
     &     escap1, tkin1, tgrav1, tterm1
c--real*4
      READ (idisk1, END=300) number
c--real*8
      READ (idisk1, END=300) number
      IF (number.NE.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) udist, umass, utime
c
c--Arrays
c
c--Number of array lengths
c
      READ (idisk1, END=300) number
      IF (number.LT.2 .OR. number.GT.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump'
         CALL quit(0)
      ENDIF
c
c--Read array type 1 header
c
      READ (idisk1, END=300) number8, (nums(i), i=1,8)
      IF (number8.NE.np1) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong'
         CALL quit(0)
      ENDIF
      np1 = number8
c
c--Read array type 2 header
c
      READ (idisk1, END=300) number8, (nums(i), i=1,8)
      nptmass1 = number8
c
c--Read array type 3 header
c
      nradtrans = 0
      IF (number.GE.3) THEN
         READ (idisk1, END=300) number8, (nums(i), i=1,8)
         IF (number8.GT.iradtrans .OR. number8.NE.1 .AND. 
     &        number8.NE.np1) THEN
            WRITE (*,*) 'ERROR 9 in rdump: iradtrans wrong ',number8,
     &           iradtrans,np1
            CALL quit(0)
         ENDIF
         nradtrans = number8
      ENDIF
c
c--Read array type 4 header
c
      nmhd = 0
      IF (number.GE.4) THEN
         READ (idisk1, END=300) number8, (nums(i), i=1,8)
         IF (number8.GT.imhd .OR. number8.NE.1 .AND. 
     &        number8.NE.np1) THEN
            WRITE (*,*) 'ERROR 10 in rdump: imhd wrong ',number8,
     &           imhd,np1
            CALL quit(0)
         ENDIF
         nmhd = number8
      ENDIF
c
c--Read array type 1 arrays
c
c--Default int
      READ (idisk1, END=300) (isteps(i), i=1, np1)
c--int*1
      READ (idisk1, END=300) (iphase(i), i=1, np1)
c--int*2,int*4,int*8
c--Default real
      DO j = 1, 5
         READ (idisk1, END=300) (xyzmh(j,i), i=1, np1)
      END DO
      DO j = 1, 4
         READ (idisk1, END=300) (vxyzu(j,i), i=1, np1)
      END DO      
c--real*4
      READ (idisk1, END=300) (rho(i), i=1, np1)
      READ (idisk1, END=300)
c     READ (idisk1, END=300) (alphaMM(1,i), i=1, np1)
c--real*8

c
c--Read array type 2 arrays
c
      IF (nptmass1.GT.0) THEN
c--Default int
         READ (idisk1, END=300) 
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 9
            READ (idisk1, END=300)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (number.GE.3 .AND. nradtrans.GT.1) THEN
c
c--Array length 3 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 5
            READ (idisk1, END=300) (ekcle(j,i), i=1, np1)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (number.GE.4 .AND. nmhd.GT.1) THEN
c
c--Array length 4 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 3
            READ (idisk1, END=300) (Bevolxyz(j,i), i=1, np1)
         END DO
c--real*4

c--real*8

      ENDIF
c
c--End reading of dump file
c--------------------------

c
c--put first object into place
c
      DO 100 j = 1, np1
         fmas1 = fmas1 + xyzmh(4,j)
         xyzmh(1,j) = xyzmh(1,j) + xx1
         xyzmh(2,j) = xyzmh(2,j) + yy1
         xyzmh(3,j) = xyzmh(3,j) + zz1
         vxyzu(1,j) = vvx1
         vxyzu(2,j) = vvy1
         vxyzu(3,j) = vvz1
         dgrav(j) = 0.
 100  CONTINUE
c
c--read second dump
c
      WRITE (*, 99010)
99010 FORMAT (' reading second file and storing data')
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Write ouput file
c
      READ (idisk1, END=300) int1i,r1i,int2i,i1i,int3i
      IF (int1i.NE.int1) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         CALL quit(0)
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default integer size wrong'
         CALL quit(0)
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default real size wrong'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) fileident
c
c--Single values
c
c--Default int
      READ (idisk1, END=300) number
      IF (number.NE.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) np2,n21,n22
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (idisk1, END=300) number
      END DO
c--Default real
      READ (idisk1, END=300) number
      IF (number.NE.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) t2, dtmax2, gamma2, rhozero2, RK22,
     &     escap2, tkin2, tgrav2, tterm2
c--real*4
      READ (idisk1, END=300) number
c--real*8
      READ (idisk1, END=300) number
      IF (number.NE.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump'
         CALL quit(0)
      ENDIF
      READ (idisk1, END=300) udist, umass, utime
c
c--Arrays
c
c--Number of array lengths
c
      READ (idisk1, END=300) number
      IF (number.LT.2 .OR. number.GT.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump'
         CALL quit(0)
      ENDIF
c
c--Read array type 1 header
c
      READ (idisk1, END=300) number8, (nums(i), i=1,8)
      IF (number8.NE.np2) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong'
         CALL quit(0)
      ENDIF
      np2 = number8
c
c--Read array type 2 header
c
      READ (idisk1, END=300) number8, (nums(i), i=1,8)
      nptmass2 = number8
c
c--Read array type 3 header
c
      IF (number.GE.3) THEN
         READ (idisk1, END=300) number8, (nums(i), i=1,8)
         IF (number8.GT.iradtrans .OR. number8.NE.1 .AND. 
     &        number8.NE.np2) THEN
            WRITE (*,*) 'ERROR 9 in rdump: iradtrans wrong ',number8,
     &           iradtrans,np2
            CALL quit(0)
         ENDIF
         nradtrans = number8
      ENDIF
c
c--Read array type 4 header
c
      IF (number.GE.4) THEN
         READ (idisk1, END=300) number8, (nums(i), i=1,8)
         IF (number8.GT.imhd .OR. number8.NE.1 .AND. 
     &        number8.NE.np2) THEN
            WRITE (*,*) 'ERROR 10 in rdump: imhd wrong ',number8,
     &           imhd,np2
            CALL quit(0)
         ENDIF
         nmhd = number8
      ENDIF
c
c--Read array type 1 arrays
c
c--Default int
      READ (idisk1, END=300) (isteps(i+np1), i=1, np2)
c--int*1
      READ (idisk1, END=300) (iphase(i+np1), i=1, np2)
c--int*2,int*4,int*8
c--Default real
      DO j = 1, 5
         READ (idisk1, END=300) (xyzmh(j,i+np1), i=1, np2)
      END DO
      DO j = 1, 4
         READ (idisk1, END=300) (vxyzu(j,i+np1), i=1, np2)
      END DO      
c--real*4
      READ (idisk1, END=300) (rho(i+np1), i=1, np2)
      READ (idisk1, END=300)
c     READ (idisk1, END=300) (alphaMM(1,i+np1), i=1, np1)
c--real*8

c
c--Read array type 2 arrays
c
      IF (nptmass1.GT.0) THEN
c--Default int
         READ (idisk1, END=300) 
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 9
            READ (idisk1, END=300)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (number.GE.3 .AND. nradtrans.GT.1) THEN
c
c--Array length 3 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 5
            READ (idisk1, END=300) (ekcle(j,i+np1), i=1, np1)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (number.GE.4 .AND. nmhd.GT.1) THEN
c
c--Array length 4 arrays
c      
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 3
            READ (idisk1, END=300) (Bevolxyz(j,i+np1), i=1, np1)
         END DO
c--real*4

c--real*8

      ENDIF
c
c--End reading of dump file
c--------------------------
c
c--rotate second object if needed
c
      crotang = cos(rotang)
      srotang = sin(rotang)
      DO 200 j = np1 + 1, np1 + np2
         fmas2 = fmas2 + xyzmh(4,j)
         xyzmh(1,j) = crotang*xyzmh(1,j) + srotang*xyzmh(2,j) + xx2
         xyzmh(2,j) = -srotang*xyzmh(1,j) + crotang*xyzmh(2,j) + yy2
         xyzmh(3,j) = xyzmh(3,j) + zz2
         vxyzu(1,j) = vvx2
         vxyzu(2,j) = vvy2
         vxyzu(3,j) = vvz2
         dgrav(j) = 0.
 200  CONTINUE

      WRITE (*, 99011)
99011 FORMAT (' what is the equation of state variable:', /,
     &        ' specific internal energy :  intener', /,
     &        ' specific entropy         :  entropy')
      READ (*, 99012) varsta
99012 FORMAT (A7)
      WRITE (*, 99013)
99013 FORMAT (' is the equation of state a gamma-law? (y/n)')
      READ (*, 99004) iok
      IF ( iok.EQ.'y' ) THEN
         WRITE (*, 99014) gamma1, gamma2
99014    FORMAT (' what is the value of the adiabatic index gamma?', /,
     &           ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) gamma
         WRITE (*, 98014) RK21, RK22
98014    FORMAT (' what is the value of the adiabatic constant RK2?', /,
     &           ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) RK2
         WRITE (*, 98015) rhozero1, rhozero2
98015    FORMAT (' what is the value of the  constant rhozero?', /,
     &           ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) rhozero
      ENDIF
      npart = np1 + np2
      n1 = np1
      n2 = np2
      WRITE (*, 99015) n1, n2, npart
99015 FORMAT (' set up completed', /,
     &        ' particles in object 1          :', I6, /,
     &        ' particles in object 2          :', I6, /,
     &        ' total number of particles used :', I6)
      WRITE (*, 99016)
99016 FORMAT (//, ' END SET UP')

      IF ( idebug(1:6).EQ.'addump' ) THEN
         WRITE (iprint, 99017) (xyzmh(1,i), i=1, npart)
         WRITE (iprint, 99017) (xyzmh(2,i), i=1, npart)
         WRITE (iprint, 99017) (xyzmh(3,i), i=1, npart)
         WRITE (iprint, 99017) (vxyzu(1,i), i=1, npart)
         WRITE (iprint, 99017) (vxyzu(2,i), i=1, npart)
         WRITE (iprint, 99017) (vxyzu(3,i), i=1, npart)
         WRITE (iprint, 99017) (vxyzu(4,i), i=1, npart)
         WRITE (iprint, 99017) (xyzmh(4,i), i=1, npart)
         WRITE (iprint, 99017) (rho(i), i=1, npart)
99017    FORMAT (1X, 5(1PE12.5,1X))
      ENDIF
      GOTO 400

 300  CALL error(where, 1)

 400  RETURN
      END
