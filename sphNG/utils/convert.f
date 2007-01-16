      PROGRAM convert
c************************************************************
c                                                           *
c  This program converts from sphNG format to a variety     *
c  (one) of other formats                                   *
c                                                           *
c  DJP & BAA 16.01.07                                       *
c                                                           *
c************************************************************
      IMPLICIT NONE
      INCLUDE '../idim'
      INTEGER npart,n1,n2,nptmass
      INTEGER*1 iphase(idim)
      INTEGER isteps(idim),listpm(iptdim)
      INTEGER nfiles,iargc,ierr,i1,i2,iformat
      REAL gt, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL*4 rho(idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim)
      CHARACTER*100 filename,fileout
c      LOGICAL smalldump
      nfiles = iargc()
      IF (nfiles.LT.2) THEN
         WRITE(*,*) 'Converts sphNG dump format to other formats'
         WRITE(*,*) 'Usage: convert inputdump outputdump'
         STOP
      ELSE
         CALL getarg(1,filename)
         CALL getarg(2,fileout)
      ENDIF
      
      i1 = 0
      i2 = 0
      
      CALL readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,rho,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)      
c
c--write new full dump
c
      write(*,*) '0: same format (e.g. to change endian using compiler)'
      write(*,*) '1: ascii'
      write(*,*)
      write(*,*) 'which format do you want to write?'
      read(*,*) iformat
      write(*,*) 'writing to dump ',fileout

      if (iformat.eq.0) then
         STOP 'NOT IMPLEMENTED'    
      elseif (iformat.eq.1) then
         CALL writedump_ascii(fileout,idim,iptdim,
     &   npart,gt,gamma,xyzmh,vxyzu,rho,iphase,nptmass,ierr)
      else
         stop 'unknown output format: doing nothing'
      endif

      WRITE(*,*) '...program terminated normally'
 
      END PROGRAM convert
