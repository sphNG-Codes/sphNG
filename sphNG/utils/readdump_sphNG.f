      SUBROUTINE readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,rho,
     &   iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,udisti,umassi,utimei,ierr)
c************************************************************
c                                                           *
c  This is rdump but as a standalone subroutine             *
c  (just essential hydro quantities+sinks at the moment)    *
c  DJP 27.09.06                                             *
c  DJP 02.11.06 support for small dumps added               *
c                                                           *
c************************************************************

      IMPLICIT NONE ! those sweet, sweet words
      INTEGER i,j
      INTEGER idim,iptdim,idisk1,ierr
      INTEGER npart,n1,n2,nreassign,naccrete,nkill,nptmass
      INTEGER nhydro,jlen
      INTEGER nums1(8),nums2(8),nums3(8),nums4(8)
      INTEGER*1 iphase(idim)
      INTEGER isteps(idim),listpm(iptdim)
      REAL hacc, gt, dtmaxdp, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz
      REAL specang, ptmassin, tiny, pmassinitial
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim)
      REAL*4 rho(idim),dummyr4(idim)
      REAL*8 umassi, udisti, utimei
      PARAMETER (tiny=1.e-12)
      LOGICAL smalldump
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
      CALL read_header(idisk1,ierr,smalldump,idim,
     &   n1,n2,nreassign,naccrete,nkill,gt,dtmaxdp,gamma,rhozero,RK2,
     &   escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &   specang, ptmassin, pmassinitial, udisti, umassi, utimei,
     &   npart,nhydro,nptmass,nums1,nums2,nums3,nums4)
c
c--Read array type 1 arrays
c
c--Default int
      IF (nums1(1).GE.1) THEN
         READ (idisk1, END=100) (isteps(i), i=1, npart)
         DO j=1,nums1(1)-1
            READ (idisk1, END=100)
         ENDDO
      ENDIF
c--int*1
      READ (idisk1, END=100) (iphase(i), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      IF (smalldump) THEN
c--set masses for equal mass particles in small dumps
         IF (abs(pmassinitial).GT.tiny) THEN
            jlen = 3
            WRITE(*,*) 'setting masses = ',pmassinitial
            DO i=1,npart
               xyzmh(4,i) = pmassinitial
            ENDDO
            IF (nums1(6).GT.3) 
     &         WRITE(*,*) '??? but masses are dumped ???',nums1(6)
         ELSE
            jlen = 4
         ENDIF
         IF (nums1(6).LT.jlen) THEN
            WRITE(*,*) 'ERROR in rdump: not enough reals (small dump)'
            ierr = 9
            RETURN
         ENDIF
         DO j = 1, jlen
            READ (idisk1, END=100) (xyzmh(j,i), i=1, npart)
         END DO
         IF (nums1(6).GE.jlen+4) THEN
c--read velocity and u if present
            DO j = 1, 4
               READ (idisk1, END=100) (vxyzu(j,i), i=1, npart)
            END DO         
         ELSE
c--skip unnecessary reals
            DO j=1,nums1(6)-jlen
               READ (idisk1, END=100)
            ENDDO
         ENDIF
c--real*4
         IF (nums1(7).LT.2) THEN
            WRITE(*,*) 'ERROR in rdump: not enough real*4s (small dump)'
            ierr = 10
            RETURN
         ELSE
            READ (idisk1, END=100) (rho(i), i=1, npart)
c--smoothing length (convert back from real*4)
            READ (idisk1, END=100) (dummyr4(i), i=1, npart)
            DO i=1,npart
               xyzmh(5,i) = dummyr4(i)
            ENDDO
         ENDIF  
c--skip unnecessary real*4's
         IF (nums1(7).GT.2) THEN
            DO j=1,nums1(7)-1
               READ (idisk1, END=100)
            ENDDO
         ENDIF

      ELSE
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
      ENDIF
c--real*8 : skip all
      DO j=1,nums1(8)
         READ(idisk1, END=100)
      ENDDO
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
      IF (smalldump) THEN
         READ (idisk1, END=100) (xyzmh(4,listpm(i)),i=1,nptmass)
      ELSE
         READ (idisk1, END=100) (spinx(i),i=1,nptmass)
         READ (idisk1, END=100) (spiny(i),i=1,nptmass)
         READ (idisk1, END=100) (spinz(i),i=1,nptmass)
         READ (idisk1, END=100) (angaddx(i),i=1,nptmass)
         READ (idisk1, END=100) (angaddy(i),i=1,nptmass)
         READ (idisk1, END=100) (angaddz(i),i=1,nptmass)
         READ (idisk1, END=100) (spinadx(i),i=1,nptmass)
         READ (idisk1, END=100) (spinady(i),i=1,nptmass)
         READ (idisk1, END=100) (spinadz(i),i=1,nptmass)
      ENDIF
c--real*4

c--real*8

      WRITE(*,*) 'finished reading (hydro) file'
      CLOSE(idisk1)
      RETURN

 100  CLOSE (idisk1)
      WRITE(*,*) ' end of file reached in data read'
      ierr = 666
      RETURN
      
      END SUBROUTINE

c************************************************************
c                                                           *
c  Separate routine for reading the header only so that     *
c  this can be done as a separate operation                 *
c  DJP 30/5/07                                              *
c                                                           *
c************************************************************

      SUBROUTINE read_header(idisk1,ierr,smalldump,idim,
     &   n1,n2,nreassign,naccrete,nkill,gt,dtmaxdp,gamma,rhozero,RK2,
     &   escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &   specang, ptmassin, pmassinitial, udisti, umassi, utimei,
     &   npart,nhydro,nptmass,nums1,nums2,nums3,nums4)
      IMPLICIT NONE
      INTEGER idisk1,ierr,i1i,idim
      INTEGER npart,n1,n2,nreassign,naccrete,nkill,nptmass
      INTEGER nhydro
      INTEGER nums1(8),nums2(8),nums3(8),nums4(8)
      LOGICAL smalldump
      REAL gt, dtmaxdp, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz
      REAL specang, ptmassin, tiny, pmassinitial
      REAL*8 umassi, udisti, utimei

      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER i,number
      INTEGER*8 number8
      REAL r1i
      CHARACTER*100 fileident
      PARAMETER (tiny=1.e-12)
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
c--if file is a small dump, return an error code but still read what
c  can be read from a small dump
c
      IF (fileident(1:1).EQ.'S') THEN
         smalldump = .true.
      ELSE
         smalldump = .false.
      ENDIF
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
      IF (smalldump .and. number.GE.15) THEN
         READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin, pmassinitial
      ELSE
         READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &        escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &        specang, ptmassin
         pmassinitial = 0.
      ENDIF
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
         ierr = 8
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

      RETURN

100   CONTINUE
    
      CLOSE (idisk1)
      WRITE(*,*) ' end of file reached in header read'
      ierr = 666
      RETURN      
      
      END SUBROUTINE read_header
