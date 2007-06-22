      PROGRAM masswithtime
c************************************************************
c                                                           *
c  This program extracts sink particle properties           *
c   from the 'P' files                                      *
c                                                           *
c  DJP May 2007                                             *
c************************************************************
      IMPLICIT NONE
      INTEGER, PARAMETER :: idim=36000000
      INTEGER, PARAMETER :: nout=100
      INTEGER, PARAMETER :: iptmass=100
      INTEGER i,j,nfiles,iargc,ierr,ioutput
      INTEGER nptmass,ipart,iout,listpm(iptmass)
      REAL pmasstemp(idim)
      REAL pmass(iptmass),timeff,time,dumi
      CHARACTER*100 filename,fileout(iptmass)
      LOGICAL pfiles
!--stuff for reading from dump files
      INTEGER npart,n1,n2,nreassign,naccrete,nkill
      INTEGER nhydro,idisk1
      INTEGER nums1(8),nums2(8),nums3(8),nums4(8)
      LOGICAL smalldump
      REAL dtmaxdp, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz
      REAL specang, ptmassin, pmassinitial, totmass
      REAL*8 umassi, udisti, utimei
      REAL, PARAMETER :: pi=3.1415926536

      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: masswithtime Pfile(s)'
      ENDIF

      pfiles = .false.
      CALL getarg(1,filename)
      IF (filename(1:1).EQ.'P') THEN
         pfiles = .true.
      ENDIF
      idisk1 = 15
c
c--open file for j vs time plots
c
c      PRINT*,' 0: sink particle masses vs time'
c      PRINT*,' 1: not implemented'
c      WRITE(*,*) 'Please select output:'
c      READ*,ioutput
c      IF (ioutput.GT.1 .OR. ioutput.LT.0) STOP 'unknown output choice'
      ioutput = 0
c
c--open file
c
      IF (ioutput.EQ.0) THEN
         iout = 55
         DO i=1,nout
            IF (pfiles) THEN
               WRITE(fileout(i),"(a,i3.3,a)") 'sink',i,'.out'
            ELSE
               WRITE(fileout(i),"(a,i3.3,a)") 'sinkd',i,'.out'
            ENDIF
            PRINT*,' opening ',fileout(i)
            OPEN(UNIT=iout+i,FILE=fileout(i),STATUS='replace',
     &           FORM='formatted')
         ENDDO
         OPEN(UNIT=iout+nout+1,FILE='sinkstot.out',STATUS='replace',
     &        FORM='formatted')
      ENDIF
      
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         OPEN (idisk1, FILE=filename, STATUS='unknown', 
     &         FORM='unformatted')
         
         ierr = 0
         DO WHILE (ierr.EQ.0)
            IF (pfiles) THEN
c
c--get sink properties from sink P file
c
               READ(idisk1,iostat=ierr) timeff,time,nptmass
               PRINT*,' t_ff= ',timeff,' t= ',time,' nptmass= ',nptmass
               IF (ierr /= 0) THEN
                  PRINT*,'ERROR READING TIME FROM FILE'
               ELSE
                  pmass = 0.
                  DO i=1,nptmass
                     READ(idisk1,iostat=ierr) ipart,dumi,dumi,dumi,
     &                                       dumi,dumi,dumi,pmass(i)
                     !PRINT*,ipart,pmass(i)
                  ENDDO
                  totmass = 0.
                  DO i=1,nptmass
                     totmass = totmass + pmass(i)
                     WRITE(iout+i,*) timeff,time,pmass(i)
                  ENDDO
                  IF (nptmass.GE.1) THEN
                     WRITE(iout+nout+1) timeff,time,totmass
                  ENDIF

                  IF (ierr /= 0) PRINT*,'ERROR READING SINKS FROM FILE'
               ENDIF
            ELSE
c
c--get sink properties from dump file
c
               CALL read_header(idisk1,ierr,smalldump,idim,
     &              n1,n2,nreassign,naccrete,nkill,time,dtmaxdp,gamma,
     &              rhozero,RK2,escap, tkin, tgrav, tterm, anglostx,
     &              anglosty, anglostz,specang, ptmassin, pmassinitial,
     &              udisti, umassi, utimei,
     &              npart,nhydro,nptmass,nums1,nums2,nums3,nums4)

               timeff = time/(SQRT((3. * pi) / (32. * rhozero)))
               PRINT*,' t_ff= ',timeff,' t= ',time,' nptmass= ',nptmass
c
c--skip normal particle block
c
               IF (nptmass.GT.0) THEN
                  IF (smalldump) THEN
                     DO i=1,SUM(nums1(1:8))
                        READ(idisk1,iostat=ierr)
                     ENDDO
                  ELSE
                     DO i=1,SUM(nums1(1:5))+3
                        READ(idisk1,iostat=ierr)
                     ENDDO
                     IF (ierr /= 0) PRINT*,'ERROR SKIPPING BLOCKS'
c
c--read masses of all particles from full dumps
c
                     READ(idisk1,iostat=ierr) (pmasstemp(i),i=1,npart)
                     DO i=1,nums1(6)-4 + SUM(nums1(7:8))
                        READ(idisk1,iostat=ierr)
                     ENDDO
                  ENDIF
                  IF (ierr /= 0) PRINT*,'ERROR READING FROM FILE'
c
c--read sink masses from block 2 of the dump file
c
c
c--Read array type 2 arrays
c
c--Default int
                  READ (idisk1,iostat=ierr) (listpm(i), i=1,nptmass)
                  IF (smalldump) THEN
                     READ (idisk1,iostat=ierr) (pmass(i),i=1,nptmass)
                  ELSE
                     DO i = 1,nptmass
                        pmass(i) = pmasstemp(listpm(i))
                     ENDDO
                  ENDIF
                  totmass = 0.
                  DO i=1,nptmass
                     totmass = totmass + pmass(i)
                     WRITE(iout+i,*) timeff,time,pmass(i)
                  ENDDO
                  IF (nptmass.GE.1) THEN
                     WRITE(iout+nout+1) timeff,time,totmass
                  ENDIF
               ENDIF
               ierr = -1 ! move on to next dump file
           ENDIF
         ENDDO
         
         CLOSE(idisk1)

      ENDDO
      
      DO i=1,nout
         IF (i.GT.nptmass) THEN
      !--delete empty files
            CLOSE(iout+i,STATUS='DELETE')
         ELSE
            CLOSE(iout+i)
         ENDIF
         CLOSE(iout+nout+1)
      ENDDO
            
      END PROGRAM masswithtime
