      SUBROUTINE readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)
c************************************************************
c                                                           *
c  This is rdump but as a standalone subroutine             *
c  (just essential hydro quantities+sinks at the moment)    *
c  DJP 27.09.06                                             *
c                                                           *
c************************************************************

      IMPLICIT NONE ! those sweet, sweet words
      INTEGER i,j,i1i
      INTEGER idim,iptdim,idisk1,ierr
      INTEGER npart,n1,n2,nreassign,naccrete,nkill,nptmass
      INTEGER nhydro,number
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      INTEGER nums1(8),nums2(8),nums3(8),nums4(8)
      INTEGER*1 iphase(idim)
      INTEGER isteps(idim),listpm(iptdim)
      REAL r1i,hacc
      REAL gt, dtmaxdp, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz
      REAL specang, ptmassin
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim)
      REAL*4 rho(idim)
      REAL*8 umassi, udisti, utimei, umagfdi
      CHARACTER*100 fileident
      CHARACTER(*) filename
c
c--unit to open file on
c     
      ierr = 0
      idisk1 = 8
      OPEN(UNIT=idisk1,FILE=filename,STATUS='old',FORM='unformatted')
c
c--Dump file
c-------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Read ouput file
c
      READ (idisk1, END=100, iostat=ierr) int1i,r1i,int2i,i1i,int3i
      IF (int1i.NE.int1 .OR. ierr.NE.0) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         ierr = 1
         RETURN
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default real size wrong'
         PRINT*,int1i,r1i,int2i
         ierr = 2
         RETURN
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default integer size wrong'
         ierr = 3
         RETURN
      ENDIF
      READ (idisk1, END=100) fileident
      WRITE(*,*) fileident
c
c--Single values
c
c--Default int
      READ (idisk1, END=100) number
      IF (number.LT.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump: not enough default ints'
         ierr = 4
         RETURN
      ENDIF
      READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill
      IF (npart.GT.idim) THEN
         WRITE (*,*) 'ERROR in rdump: npart>idim'
         STOP
      ELSE
         WRITE (*,*) 'npart = ',npart
      ENDIF
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (idisk1, END=100) number
      END DO
c--Default real
      READ (idisk1, END=100) number
      IF (number.LT.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump: not enough default reals'
         ierr = 5
         RETURN
      ENDIF
         READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin
      WRITE(*,*) 'time = ',gt
c--real*4
      READ (idisk1, END=100) number
c--real*8
      READ (idisk1, END=100) number
      IF (number.LT.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump: nreal8 too small header'
         ierr = 6
         RETURN
      ENDIF
      READ (idisk1, END=100) udisti, umassi, utimei

c
c--Arrays
c
c--Number of array lengths
c
      READ (idisk1, END=100) number
      IF (number.LT.2 .OR. number.GT.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump'
         ierr = 7
         RETURN
      ENDIF
c
c--Read array type 1 header
c
      READ (idisk1, END=100) number8, (nums1(i), i=1,8)
      IF (number8.LT.npart) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong',number8,
     &           npart
         RETURN
      ENDIF
      nhydro = number8
c
c--Read array type 2 header
c
      READ (idisk1, END=100) number8, (nums2(i), i=1,8)
      nptmass = number8
      WRITE(*,*) 'nptmass = ',nptmass
c
c--Read array type 3 header
c
      IF (number.GE.3) THEN
         READ (idisk1, END=100) number8
      ENDIF
c
c--Read array type 4 header
c
      IF (number.GE.4) THEN
         READ (idisk1, END=100)
      ENDIF
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
c--skip unnecessary reals
      IF (nums1(6).GT.9) THEN
         DO j=1,nums1(6)-9
            READ (idisk1, END=100)
         ENDDO
      ENDIF    
c--real*4
      READ (idisk1, END=100) (rho(i), i=1, npart)     
c--skip unnecessary real*4's
      IF (nums1(7).GT.1) THEN
         DO j=1,nums1(7)-1
            READ (idisk1, END=100)
         ENDDO
      ENDIF
c     READ (idisk1, END=100) (alphaMM(1,i), i=1, npart)
c--real*8

c
c--Read array type 2 arrays
c
c--Default int
      READ (idisk1, END=100) (listpm(i), i=1,nptmass)

      DO i = 1, nptmass
         WRITE(*,*) 'Setting ',i,listpm(i),xyzmh(5,listpm(i)),hacc
         xyzmh(5,listpm(i)) = hacc
      END DO


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

      WRITE(*,*) 'finished reading (hydro) file'
      CLOSE(idisk1)
      RETURN

 100  CLOSE (idisk1)
      WRITE(*,*) ' end of file reached in data read'
      RETURN
      
      END SUBROUTINE
