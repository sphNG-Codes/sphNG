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
      INTEGER i,j,nfiles,iargc,ierr,ipart
      REAL gt, gamma, rhozero, RK2
      REAL escap, tkin, tgrav, tterm
      REAL xyzmh(5,idim),vxyzu(4,idim)
      REAL spinx(iptdim),spiny(iptdim),spinz(iptdim)
      REAL angaddx(iptdim),angaddy(iptdim),angaddz(iptdim)
      REAL spinadx(iptdim),spinady(iptdim),spinadz(iptdim) 
      REAL angx,angy,angz,angxi,angyi,angzi,angto
      REAL angm(idim),spinxi,spinyi,spinzi
      REAL totmass,mass,angmom,totang
      CHARACTER*100 filename,fileout
      
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: massvsangmom dumpfile(s)'
      ENDIF
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         CALL readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)
         IF (ierr /= 0) STOP 'aborting...'
c
c--calculate specific AM for each particle
c  also total angular momentum
c
         angx = 0.
         angy = 0.
         angz = 0.
         totmass = 0.
         totang = 0.
         DO i = 1, npart
            IF (iphase(i).GE.0) THEN
               angxi = (xyzmh(2,i)*vxyzu(3,i) - vxyzu(2,i)*xyzmh(3,i))
               angyi = (vxyzu(1,i)*xyzmh(3,i) - xyzmh(1,i)*vxyzu(3,i))
               angzi = (xyzmh(1,i)*vxyzu(2,i) - vxyzu(1,i)*xyzmh(2,i))
               angm(i) = SQRT(angxi**2 + angyi**2 + angzi**2)
               totang = totang + angm(i)
               angx = angx + xyzmh(4,i)*angxi
               angy = angy + xyzmh(4,i)*angyi
               angz = angz + xyzmh(4,i)*angzi
               totmass = totmass + xyzmh(4,i)
            ELSE
               angm(i) = 0.
            ENDIF
         END DO
c
c--Add spin angular momentum of point masses
c
         DO i = 1, nptmass
            spinxi = spinx(i)
            spinyi = spiny(i)
            spinzi = spinz(i)
            angx = angx + spinxi
            angy = angy + spinyi
            angz = angz + spinzi
            ipart = listpm(i)
            angm(ipart) = SQRT(spinxi**2 + spinyi**2 + spinzi**2)
     &                     /xyzmh(4,ipart)
            totang = totang + angm(ipart)
         END DO
         angto = SQRT(angx**2 + angy**2 + angz**2)
         PRINT*,' total angular momentum = ',angto,angto/totmass
         PRINT*,' total mass = ',totmass
c
c--now sort particles by specific AM
c
         CALL indexx2(npart,angm,list)
c
c--print out in list order, cumulative mass vs specific AM
c
         fileout = trim(filename)//'.am'
         PRINT*,' writing to file ',fileout
         OPEN(UNIT=66,FILE=fileout,STATUS='replace',FORM='formatted')
         mass = 0.
         angmom = 0.
         DO ipart=1, npart
            i = list(ipart)
            IF (iphase(i).GE.0) THEN
               mass = mass + xyzmh(4,i)
               angmom = angmom + angm(i)
               WRITE(66,*) mass/totmass,angm(i),angmom/totang  !!!/(angto/totmass)
            ENDIF
         ENDDO
         CLOSE(UNIT=66)
      END DO
      
      END PROGRAM massvsangmom
