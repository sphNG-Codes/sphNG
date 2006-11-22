      PROGRAM discmass
c************************************************************
c                                                           *
c  This program extracts properties of a disc               *
c  (mass, radius, angular momentum)                         *
c  from a series of dump files                              *
c                                                           *
c  DJP Nov 2006                                             *
c************************************************************
      IMPLICIT NONE
      INCLUDE '../idim'
      INTEGER npart,n1,n2,nptmass
      INTEGER*1 iphase(idim)
      INTEGER isteps(idim),listpm(iptdim)
      INTEGER i,j,nfiles,iargc,ierr
      REAL gt, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim) 
      REAL xcm1,ycm1,zcm1,dx,dy,dz
      REAL pmassi
      REAL discmas,rdisc,dr,totmass
      REAL radius(idim),fracmass,frac
      REAL*4 rho(idim)
      REAL discr50pc,discr90pc,discr95pc,discr99pc
      CHARACTER*100 filename,fileout
      INTEGER imax,ioutput,ndisc,ipart,list(idim),listr(idim)
      LOGICAL smalldump
      INTEGER nbinsize,ibindim,nlastbin,ibin,nbins
      PARAMETER(nbinsize=100)
      PARAMETER(ibindim=idim/nbinsize)
      REAL rmean(ibindim),rhomean(ibindim),rhomax(ibindim)
      REAL drhodr
            
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: massvsangmom dumpfile(s)'
      ENDIF
c
c--open file for j vs time plots
c
      PRINT*,' 0: disc mass, radius etc vs time'
      PRINT*,' 1: disc particles selected for a given dump'
      PRINT*,' 2: disc rho vs r'
      WRITE(*,*) 'Please select output:'
      READ*,ioutput
      IF (ioutput.GT.2 .OR. ioutput.LT.0) STOP 'unknown output choice'
c
c--open file
c
      IF (ioutput.EQ.0) THEN
         fileout = 'disc.out'
         PRINT*,' opening ',fileout
         OPEN(UNIT=55,FILE=fileout,STATUS='replace',FORM='formatted')
      ENDIF
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         ierr = -1 ! means do not read small dump files
         CALL readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,rho,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)
         IF (ierr.GT.0) STOP 'aborting...'
         IF (ierr.EQ.-1) THEN
c            PRINT*,'WARNING: small dump file: AM not present'
            smalldump = .true.
c            PRINT*,'skipping small dump file...'
c            GOTO 100
         ELSE
            smalldump = .false.
         ENDIF
         
         IF (nptmass.GT.1) THEN
            STOP 'not implemented for more than one point mass'
         ELSEIF (nptmass.EQ.1) THEN
            PRINT*,'Point mass found: calculating disc properties'
            imax = listpm(1)
         ELSE
            PRINT*,'No point masses: using density maximum to get disc'
            PRINT*, '!!! not yet implemented'
            GOTO 100
         ENDIF
         xcm1 = xyzmh(1,imax)
         ycm1 = xyzmh(2,imax)
         zcm1 = xyzmh(3,imax)
         discmas = 0. !xyzmh(4,imax)
         totmass = xyzmh(4,imax)
         rdisc = 0.
         ndisc = 0
         discr50pc = 0.
         discr90pc = 0.
         discr95pc = 0.
         discr99pc = 0.
         
         IF (ioutput.EQ.1) OPEN(UNIT=10,FILE=trim(filename)//'.disc',
     &        STATUS='replace',FORM='formatted')
         IF (ioutput.EQ.2) OPEN(UNIT=66,FILE=trim(filename)//'.discr',
     &        STATUS='replace',FORM='formatted')
                  
         npart = 0
         DO i=1,n1
            IF (iphase(i).EQ.0) THEN
               dx = xyzmh(1,i) - xcm1
               dy = xyzmh(2,i) - ycm1
               dz = xyzmh(3,i) - zcm1
               pmassi = xyzmh(4,i)
               dr = sqrt(dx**2 + dy**2 + dz**2)
               npart = npart + 1
               list(npart) = i
               radius(npart) = dr
               totmass = totmass + pmassi
            ENDIF
         ENDDO
c
c--sort particles by radius from point mass
c
         CALL indexx2(npart,radius,listr)
c
c--now take bins of 100 particles and calculate radius and mean density for bin
c         
         ibin = 1
         rmean(1) = 0.
         rhomean(1) = 0.
         rhomax(1) = 0.
         nlastbin = 0
         DO ipart=1,npart
            IF (ipart.GT.1 .and. mod(ipart-1,nbinsize).EQ.0) THEN
!               print*,rmean(ibin),rhomean(ibin)
               rhomean(ibin) = rhomean(ibin)/nbinsize
               rmean(ibin) = rmean(ibin)/nbinsize
!               print*,'end of bin ',ibin,nlastbin,nbinsize,
!     &                'r = ',rmean(ibin),' rhomean = ',rhomean(ibin)
               ibin = ibin + 1
               rmean(ibin) = 0.
               rhomean(ibin) = 0.
               rhomax(ibin) = 0.
               nlastbin = 0
            ENDIF
            i = list(listr(ipart))
            nlastbin = nlastbin + 1
            rmean(ibin) = rmean(ibin) + radius(listr(ipart))
            rhomean(ibin) = rhomean(ibin) + rho(i)
            rhomax(ibin) = max(rhomax(ibin),rho(i))
         ENDDO
         rhomean(ibin) = rhomean(ibin)/nlastbin
         rmean(ibin) = rmean(ibin)/nlastbin
         nbins = ibin
!         print*,'end of bin ',ibin,nlastbin,nbinsize,
!     &          'r = ',rmean(ibin),' rhomean = ',rhomean(ibin)
c
c--extract disc size based on where the gradient of log(rho) vs log(r)
c  has a steep drop-off
c        
         drhodr = 0.
         nlastbin = 0
         DO ibin=1,nbins
            IF (ibin.GT.1) THEN
               drhodr =(log(rhomean(ibin)) - log(rhomean(ibin-1)))
     &                  /(log(rmean(ibin)) - log(rmean(ibin-1)))
            ENDIF
            IF (nlastbin.EQ.0 .AND. (rhomax(ibin).LT.5 
     &          .OR. drhodr .LT. -5.)) nlastbin = ibin
            IF (ioutput.EQ.2) THEN
               WRITE(66,*) rmean(ibin),rhomean(ibin),drhodr
            ENDIF
         ENDDO
         IF (nlastbin.EQ.0) THEN
            WRITE(*,*) 'COULD NOT EXTRACT DISC!!!!'
            GOTO 100
         ENDIF
         
         ndisc = nbinsize*nlastbin
c
c--calculate total mass in disc
c
         discmas = 0.
         rdisc = 0.
         DO ipart=1,ndisc
            i = list(listr(ipart))
            discmas = discmas + xyzmh(4,i)
            rdisc = max(rdisc,radius(listr(ipart)))
            IF (ioutput.EQ.1) THEN
               WRITE(10,*) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),
     &                     xyzmh(4,i),xyzmh(5,i),rho(i)
            ENDIF
         ENDDO
c
c--calculate disc radii for half/90%/95%/99% of disc mass
c
         fracmass = 0. !!xyzmh(4,imax)
         DO ipart=1,ndisc
            i = listr(ipart)
            fracmass = fracmass + xyzmh(4,list(i))
            frac = fracmass/discmas
            IF (frac.LT.0.5) THEN
               discr50pc = radius(i)
            ENDIF
            IF (frac.LT.0.9) THEN
               discr90pc = radius(i)
            ENDIF
            IF (frac.LT.0.95) THEN
               discr95pc = radius(i)
            ENDIF
            IF (frac.LT.0.99) THEN
               discr99pc = radius(i)        
            ENDIF           
         ENDDO
         IF (ioutput.EQ.1) CLOSE(UNIT=10)
         IF (ioutput.EQ.2) CLOSE(UNIT=66)
         IF (ioutput.EQ.0) THEN
            WRITE(55,*) gt,discmas,discmas/totmass,discr50pc,discr90pc,
     &                 discr95pc,discr99pc,rdisc
         ENDIF
         WRITE(*,*) ' npart in disc = ',ndisc
         WRITE(*,*) ' point mass / (disc + pt mass) = ',
     &               xyzmh(4,imax)/(discmas + xyzmh(4,imax))
         WRITE(*,*) 't = ',gt,'mass in disc = ',discmas,discmas/totmass,
     &        'mtot=',totmass,'disc radius=',discr50pc,discr90pc,
     &        discr95pc,discr99pc,rdisc

100   END DO
      
      IF (ioutput.EQ.0) CLOSE(UNIT=55)
               
      END PROGRAM discmass
