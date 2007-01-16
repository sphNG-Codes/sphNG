      SUBROUTINE writedump_ascii(fileout,idim,iptdim,
     &   npart,gt,gamma,xyzmh,vxyzu,rho,iphase,nptmass,ierr)
      IMPLICIT NONE
      CHARACTER(*) fileout
      INTEGER idim,iptdim,npart,ierr,nptmass,i
      REAL gt,gamma
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL*4 rho(idim)
      INTEGER*1 iphase(idim)

      OPEN(UNIT=23,FILE=fileout,STATUS='REPLACE',FORM='FORMATTED')
c
c--write header line
c
      WRITE(23,*) gt,npart,nptmass,gamma
c
c--write ascii body
c
      DO i=1,npart
         WRITE(23,*) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),xyzmh(4,i),
     &    xyzmh(5,i),rho(i),
     &    vxyzu(1,i),vxyzu(2,i),vxyzu(3,i),vxyzu(4,i),iphase(i)
      ENDDO
      CLOSE(UNIT=23)
      WRITE(*,*) fileout//' written successfully'
      
      RETURN
      
      END SUBROUTINE writedump_ascii
