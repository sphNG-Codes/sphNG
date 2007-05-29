      PROGRAM masswithtime
c************************************************************
c                                                           *
c  This program extracts sink particle properties           *
c   from the 'P' files                                      *
c                                                           *
c  DJP May 2007                                             *
c************************************************************
      IMPLICIT NONE
      INTEGER, PARAMETER :: idim=100
      INTEGER, PARAMETER :: nout=100
      INTEGER i,j,nfiles,iargc,ierr,ioutput
      INTEGER nptmass,ipart,iout
      REAL pmass(idim),timeff,time,dumi
      CHARACTER*100 filename,fileout(idim)

      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: masswithtime Pfile(s)'
      ENDIF
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
            WRITE(fileout(i),"(a,i3.3,a)") 'sink',i,'.out'
            PRINT*,' opening ',fileout(i)
            OPEN(UNIT=iout+i,FILE=fileout(i),STATUS='replace',
     &           FORM='formatted')
         ENDDO
      ENDIF
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         OPEN (15, FILE=filename, STATUS='unknown', FORM='unformatted')
         
         ierr = 0
         
         DO WHILE (ierr.EQ.0)
            READ(15,iostat=ierr) timeff,time,nptmass
            PRINT*,' t_ff = ',timeff,' t = ',time,' nptmass = ',nptmass
            IF (ierr /= 0) THEN
                PRINT*,'ERROR READING TIME FROM FILE'
            ELSE
               pmass = 0.
               DO i=1,nptmass
                  READ(15,iostat=ierr) ipart,dumi,dumi,dumi,
     &                                       dumi,dumi,dumi,pmass(i)
                  !PRINT*,ipart,pmass(i)
               ENDDO
               DO i=1,nptmass
                  WRITE(iout+i,*) timeff,time,pmass(i)
               ENDDO
         
               IF (ierr /= 0) PRINT*,'ERROR READING SINKS FROM FILE'
            ENDIF
         ENDDO
         
         CLOSE(15)

      ENDDO
      
      DO i=1,nout
         IF (i.GT.nptmass) THEN
      !--delete empty files
            CLOSE(iout+i,STATUS='DELETE')
         ELSE
            CLOSE(iout+i)
         ENDIF
      ENDDO
            
      END PROGRAM masswithtime
