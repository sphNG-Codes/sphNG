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
      REAL xcm1,ycm1,zcm1,dx,dy,dz,vcmx,vcmy,vcmz
      REAL pmassi,dv2,vpotdif,angxtot,angytot,angztot
      REAL dvx,dvy,dvz,dangx,dangy,dangz,dang2,denergy,decc
      REAL discmas,rdisc,amdisc,dpmasstot,da,dr,totmass
      REAL radius(idim),fracmass,frac
      REAL*4 rho(idim)
      REAL discr50pc,discr90pc,discr95pc,discr99pc
      CHARACTER*100 filename,fileout
      INTEGER imax,ioutput,ndisc,ipart,list(idim),listr(idim)
      LOGICAL smalldump
      INTEGER nbinsize,ibindim,nlastbin,ibin,nbins
      PARAMETER(nbinsize=100)
      PARAMETER(ibindim=idim/nbinsize)
      REAL rmean(ibindim),rhomean(ibindim)
            
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: massvsangmom dumpfile(s)'
      ENDIF
c
c--open file for j vs time plots
c
      PRINT*,' 0: disc mass, radius etc vs time'
      PRINT*,' 1: disc particles selected for a given dump'
      WRITE(*,*) 'Please select output:'
      READ*,ioutput
      IF (ioutput.GT.1 .OR. ioutput.LT.0) STOP 'unknown output choice'
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
            PRINT*,'WARNING: small dump file: AM not present'
            smalldump = .true.
            PRINT*,'skipping small dump file...'
            GOTO 100
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
         vcmx = vxyzu(1,imax)
         vcmy = vxyzu(2,imax)
         vcmz = vxyzu(3,imax)
         discmas = 0. !xyzmh(4,imax)
         totmass = xyzmh(4,imax)
         rdisc = 0.
         angxtot = spinx(1)
         angytot = spiny(1)
         angztot = spinz(1)
         ndisc = 0
         discr50pc = 0.
         discr90pc = 0.
         discr95pc = 0.
         discr99pc = 0.
         
         IF (ioutput.EQ.1) OPEN(UNIT=10,FILE=trim(filename)//'.disc',
     &        STATUS='replace',FORM='formatted')
                  
         npart = 0
         DO i=1,n1
            IF (iphase(i).EQ.0) THEN
               dx = xyzmh(1,i) - xcm1
               dy = xyzmh(2,i) - ycm1
               dz = xyzmh(3,i) - zcm1
               pmassi = xyzmh(4,i)
               dr = sqrt(dx**2 + dy**2 + dz**2)
               radius(i) = dr
               npart = npart + 1
               list(npart) = i
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
         nlastbin = 0
         DO ipart=1,npart
            IF (mod(ipart,nbinsize).EQ.0) THEN
               print*,rmean(ibin),rhomean(ibin)
               rhomean(ibin) = rhomean(ibin)/nbinsize
               rmean(ibin) = rmean(ibin)/nbinsize
               print*,'end of bin ',ibin,nlastbin,nbinsize,
     &                'r = ',rmean(ibin),' rhomean = ',rhomean(ibin)
               ibin = ibin + 1
               rmean(ibin) = 0.
               rhomean(ibin) = 0.
               nlastbin = 0
            ENDIF
            i = list(listr(ipart))
            nlastbin = nlastbin + 1
            if (ibin.eq.1) then
               print*,i,radius(i)
            endif
            rmean(ibin) = rmean(ibin) + radius(i)
            rhomean(ibin) = rhomean(ibin) + rho(i)
         ENDDO
         rhomean(ibin) = rhomean(ibin)/nlastbin
         rmean(ibin) = rmean(ibin)/nlastbin
         nbins = ibin
         DO ibin=1,nbins
            WRITE(55,*) rmean(ibin),rhomean(ibin)
         ENDDO
         print*,'end of bin ',ibin,nlastbin,nbinsize,
     &          'r = ',rmean(ibin),' rhomean = ',rhomean(ibin)
         
         STOP
         DO i=1,n1
            IF (iphase(i).EQ.0) THEN
               dx = xyzmh(1,i) - xcm1
               dy = xyzmh(2,i) - ycm1
               dz = xyzmh(3,i) - zcm1
               pmassi = xyzmh(4,i)
               dr = sqrt(dx**2 + dy**2 + dz**2)
               dvx = vxyzu(1,i) - vcmx
               dvy = vxyzu(2,i) - vcmy
               dvz = vxyzu(3,i) - vcmz
               dv2 = dvx*dvx + dvy*dvy + dvz*dvz
               totmass = totmass + pmassi
c
c--work out which particles are in the disc
c
c--check if particle is bound to point mass
c  (ie. v^2/r < GM/r^2)
               vpotdif = -xyzmh(4,imax)/dr + dv2
c
c--check if specific ang mom < that required to form circular orbit
c
c               vrad = (dvx*dx + dvy*dy + dvz*dz)/dr
c               vtan2 = dv2 - vrad*vrad
c               specangmom2 = vtan2*dr*dr
c               specangmomr2 = xyzmh(4,imax)*dr
c
c--check if eccentricity is small
c               
               dpmasstot = pmassi + xyzmh(4,imax)

               denergy = - dpmasstot/dr
     &                 + (dvx*dvx + dvy*dvy + dvz*dvz)/2.0
               dangx = dvy*dz - dvz*dy
               dangy = dvz*dx - dvx*dz
               dangz = dvx*dy - dvy*dx
               dang2 = (dangx*dangx + dangy*dangy + dangz*dangz)

               da = - dpmasstot/2.0/denergy
               decc = 1.0D+0 - dang2/dpmasstot/da
               IF (decc.LT.-0.000001) THEN
                  WRITE (*,*) 'ERROR: ecc^2 = ',decc
               ELSEIF (decc.LT.0.0) THEN
                  decc = 0.0
               ENDIF
               decc = SQRT(decc)
c
c--disc particles are those on nearly circular, bound orbits around centre of mass
c  (defined as eccentricity < 0.2)
c
               IF (vpotdif.LE.0. .AND. decc.LT.0.5) THEN
c                  PRINT*,i,' is in disc, decc = ',decc,'vpotdif = ',
c     &             vpotdif,'m = ',xyzmh(4,i),' r = ',dr,rdisc
                  ndisc = ndisc + 1
                  radius(ndisc) = dr
                  list(ndisc) = i
                  discmas = discmas + xyzmh(4,i)
                  rdisc = MAX(rdisc,dr)
                  angxtot = angxtot + dangx
                  angytot = angytot + dangy
                  angztot = angztot + dangz
                  IF (ioutput.EQ.1) THEN
                     WRITE(10,*) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),
     &                        xyzmh(4,i),xyzmh(5,i),rho(i)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
c
c--calculate disc radii (radii within which 50%, 90% and 99% of mass is contained)
c
c  sort disc particles by radius
         CALL indexx2(ndisc,radius,listr)
         fracmass = 0. !!xyzmh(4,imax)
         DO ipart=1,ndisc
            i = listr(ipart)
            fracmass = fracmass + xyzmh(4,list(i))
            frac = fracmass/discmas
c            print*,list(i),radius(i),xyzmh(4,list(i)),
c     &  iphase(list(i)),frac
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
c
c--spit out disc mass and radius for this time
c
         amdisc = SQRT(angxtot**2 + angytot**2 + angztot**2)
         IF (ioutput.EQ.0) THEN
            WRITE(55,*) gt,discmas,discmas/totmass,discr50pc,discr90pc,
     &                 discr95pc,discr99pc,rdisc,amdisc
         ENDIF
         WRITE(*,*) ' npart in disc = ',ndisc
         WRITE(*,*) ' point mass / (disc + pt mass) = ',
     &               xyzmh(4,imax)/(discmas + xyzmh(4,imax))
         WRITE(*,*) 't = ',gt,'mass in disc = ',discmas,discmas/totmass,
     &        'mtot=',totmass,'disc radius=',discr50pc,discr90pc,
     &        discr95pc,discr99pc,rdisc,' disc AM = ',amdisc
100   END DO
      
      IF (ioutput.EQ.0) CLOSE(UNIT=55)
               
      END PROGRAM discmass
