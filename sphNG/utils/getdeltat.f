      PROGRAM getdeltat 
c************************************************************
c                                                           *
c  gets delta t between dump files DJP 02.11.06             *
c                                                           *
c************************************************************
      IMPLICIT NONE
      INCLUDE '../idim'
      INTEGER npart,n1,n2,nptmass
      INTEGER*1 iphase(idim)
      INTEGER isteps(idim),listpm(iptdim)
      INTEGER i,j,nfiles,iargc,ierr,i1,i2
      REAL gt, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL*4 rho(idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim)
      CHARACTER*100 filename
      REAL tprev
c      LOGICAL smalldump
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: getdeltat dumpfile(s)'
      ENDIF
      
      tprev = 0.
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         ierr = 666
         CALL readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,rho,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)
         !IF (ierr.GT.0) STOP 'aborting...'
c         IF (ierr.EQ.-1) THEN
c            PRINT*,'WARNING: small dump file: AM not present'
c            smalldump = .true.
c         ELSE
c            smalldump = .false.
c         ENDIF
         print*,'time = ',gt,' delta t = ',gt-tprev 
         tprev = gt
      ENDDO
      
      CLOSE(UNIT=55)
      STOP

      END PROGRAM
