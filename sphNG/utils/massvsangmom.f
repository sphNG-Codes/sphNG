      PROGRAM massvsangmom
c************************************************************
c                                                           *
c  This subroutine gets 
c  cumulative mass vs specific angular momentum             *
c  from a series of dump files
c************************************************************
      IMPLICIT NONE
      INCLUDE '../idim'
      INTEGER npart,n1,n2,nptmass
      INTEGER*1 iphase(idim)
      INTEGER isteps(idim),listpm(iptdim),list(idim)
      INTEGER i,j,k,nfiles,iargc,ierr,ipart
      REAL gt, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim) 
      REAL angx,angy,angz,angxi,angyi,angzi,angto
      REAL angx1,angy1,angz1,angto1
      REAL angm(idim),radius(idim),spinxi,spinyi,spinzi
      REAL totmass,mass,totmass1
      CHARACTER*100 filename,fileout
      INTEGER ioutput,njvals,nrvals
      PARAMETER (njvals=6,nrvals=9)
      REAL massbelowJ(njvals),angmval(njvals)
      REAL masswithinR(nrvals),rval(nrvals)
      LOGICAL smalldump
      
c      DATA angmval /1.e-3,1.e-2,0.1,1.0,5.0/
      DATA angmval /1.e-3,1.e-2,0.1,0.2,0.5,1.0/
      DATA rval /0.02,0.05,0.1,0.2,0.5,1.0,2.0,3.0,4.0/
      
      
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: massvsangmom dumpfile(s)'
      ENDIF
c
c--open file for j vs time plots
c
      PRINT*,' 0: mass frac vs ang mom'
      PRINT*,' 1: mass fraction below J value vs time'
      PRINT*,' 2: mass fraction vs radius'
      PRINT*,' 3: mass fraction within a given radius vs time'
      WRITE(*,*) 'Please select output:'
      READ*,ioutput
      IF (ioutput.GT.3 .OR. ioutput.LT.0) STOP 'unknown output choice'
      
      IF (ioutput.EQ.1) THEN
         fileout = 'mbelowj.out'
         PRINT*,' opening ',fileout
         OPEN(UNIT=55,FILE=fileout,STATUS='replace',FORM='formatted')
      ELSEIF (ioutput.EQ.3) THEN
         fileout = 'mwithinr.out'
         PRINT*,' opening ',fileout
         OPEN(UNIT=55,FILE=fileout,STATUS='replace',FORM='formatted')
      ENDIF
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         CALL readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)
         IF (ierr.GT.0) STOP 'aborting...'
         IF (ierr.EQ.-1) THEN
            PRINT*,'WARNING: small dump file: AM not present'
            smalldump = .true.
            IF (ioutput.LE.1) THEN
               PRINT*,'skipping small dump file...'
               GOTO 100
            ENDIF
         ELSE
            smalldump = .false.
         ENDIF
c
c--calculate specific AM for each particle
c  also total angular momentum
c
         angx = 0.
         angy = 0.
         angz = 0.
         angx1 = 0.
         angy1 = 0.
         angz1 = 0.
         totmass = 0.
         totmass1 = 0.
         print*,'npart=',npart,'n1=',n1,'n2=',n2
         DO i = 1, npart
            IF (iphase(i).GE.0) THEN
               radius(i) = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 
     &                                        + xyzmh(3,i)**2)
               angxi = (xyzmh(2,i)*vxyzu(3,i) - vxyzu(2,i)*xyzmh(3,i))
               angyi = (vxyzu(1,i)*xyzmh(3,i) - xyzmh(1,i)*vxyzu(3,i))
               angzi = (xyzmh(1,i)*vxyzu(2,i) - vxyzu(1,i)*xyzmh(2,i))
               angm(i) = angxi**2 + angyi**2 + angzi**2
               angx = angx + xyzmh(4,i)*angxi
               angy = angy + xyzmh(4,i)*angyi
               angz = angz + xyzmh(4,i)*angzi
               totmass = totmass + xyzmh(4,i)
               IF (i.EQ.n1) THEN
                  angx1 = angx
                  angy1 = angy
                  angz1 = angz
                  totmass1 = totmass
               ENDIF
            ELSE
               angm(i) = 0.
            ENDIF
         END DO
c
c--Add spin angular momentum of point masses
c
         IF (.not.smalldump) THEN
            DO i = 1, nptmass
               spinxi = spinx(i)
               spinyi = spiny(i)
               spinzi = spinz(i)
               angx = angx + spinxi
               angy = angy + spinyi
               angz = angz + spinzi
               ipart = listpm(i)
               angm(ipart) = angm(ipart) + 
     &             (spinxi**2 + spinyi**2 + spinzi**2)/xyzmh(4,ipart)**2
            END DO
            angto = SQRT(angx**2 + angy**2 + angz**2)
            angto1 = SQRT(angx1**2 + angy1**2 + angz1**2)
            PRINT*,' total angular momentum = ',angto,'n1=',angto1,
     &           angto1/angto
         ENDIF
         PRINT*,' total mass = ',totmass

         IF (ioutput.EQ.0) THEN
c
c--now sort particles by specific AM
c
            CALL indexx2(n1,angm,list)
c
c--print out in list order, cumulative mass vs specific AM
c
            fileout = trim(filename)//'.am'
            PRINT*,' writing to file ',fileout
            OPEN(UNIT=66,FILE=fileout,STATUS='replace',FORM='formatted')
            mass = 0.
            DO ipart=1, n1
               i = list(ipart)
               IF (iphase(i).GE.0) THEN
                  mass = mass + xyzmh(4,i)
                  WRITE(66,*) mass/totmass1,sqrt(angm(i))
               ENDIF
            ENDDO
            CLOSE(UNIT=66)
         ELSEIF (ioutput.EQ.1) THEN
c
c--print out time slices : mass fraction below a certain J
c
            DO k=1,njvals
               massbelowJ(k) = 0.
            ENDDO
            DO i=1,n1
               IF (iphase(i).GE.0) THEN
                  DO k=1,njvals
                     IF (sqrt(angm(i)).LE.angmval(k)) THEN
                        massbelowJ(k) = massbelowJ(k) + xyzmh(4,i)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            DO k=1,njvals
               massbelowJ(k) = massbelowJ(k)/totmass1
            ENDDO
            
            WRITE(*,*) 't=',gt,'AM=',angto1,angto1/angto,
     &                 (massbelowJ(k),k=1,njvals)
            WRITE(55,*) gt,angto1,angto1/angto,
     &                 (massbelowJ(k),k=1,njvals)
         
         ELSEIF (ioutput.EQ.2) THEN
c
c--sort particles by radius
c  print out cumulative mass vs spherical radius
c  (ie mass contained within a radius r, for various radius values
c
            CALL indexx2(n1,radius,list)
c
c--print out in list order, cumulative mass vs radius
c
            fileout = trim(filename)//'.mvsr'
            PRINT*,' writing to file ',fileout
            OPEN(UNIT=66,FILE=fileout,STATUS='replace',FORM='formatted')
            mass = 0.
            DO ipart=1, n1
               i = list(ipart)
               IF (iphase(i).GE.0) THEN
                  mass = mass + xyzmh(4,i)
                  WRITE(66,*) radius(i),mass/totmass1,sqrt(angm(i))
               ENDIF
            ENDDO
            CLOSE(UNIT=66)
         
         ELSEIF (ioutput.EQ.3) THEN
c
c--print out time slices : mass fraction contained within a certain R
c
            DO k=1,nrvals
               masswithinR(k) = 0.
            ENDDO
            DO i=1,n1
               IF (iphase(i).GE.0) THEN
                  DO k=1,nrvals
                     IF (radius(i).LE.rval(k)) THEN
                        masswithinR(k) = masswithinR(k) + xyzmh(4,i)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            DO k=1,nrvals
               masswithinR(k) = masswithinR(k)/totmass1
            ENDDO
            
            WRITE(*,*) 't=',gt,'M(<R)=',totmass1,
     &                 (masswithinR(k),k=1,nrvals)
            WRITE(55,*) gt,totmass1,
     &                 (masswithinR(k),k=1,nrvals)         
         ENDIF
100   END DO
      
      IF (ioutput.EQ.1 .OR. ioutput.EQ.3) CLOSE(UNIT=55)
               
      END PROGRAM massvsangmom
