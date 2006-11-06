      PROGRAM binaryproperties
c************************************************************
c                                                           *
c  This program extracts properties of a binary star system *
c  such as binary separation and position angle             *
c  from a series of dump files                              *
c                                                           *
c  If 2 sink particles present, uses sink properties        *
c                                                           *
c  If < 2 sinks looks for two density maxima in opposite    *
c  hemispheres.                                             *
c                                                           *
c  DJP 02.11.06                                             *
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
      REAL dx,dy,dz,dr,r1,r2,x1,y1,z1,x2,y2,z2,rhomax1,rhomax2
      REAL r2dotr1,angle,pi,rtol,ri
      PARAMETER (pi=3.1415926536)
      PARAMETER (rtol=0.05)
      CHARACTER*100 filename,fileout
c      LOGICAL smalldump
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: massvsangmom dumpfile(s)'
      ENDIF
      
      i1 = 0
      i2 = 0
      x1 = 0.
      y1 = 0.
      z1 = 0.
      r1 = -1.0
      rhomax1 = 0.
      fileout = 'binary.out'
      PRINT*,' opening ',fileout
      OPEN(UNIT=55,FILE=fileout,STATUS='replace',FORM='formatted')
      
      DO j=1,nfiles
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
         CALL readdump_sphNG(filename,idim,iptdim,
     &   npart,n1,n2,gt,gamma,rhozero,RK2,
     &   escap,tkin,tgrav,tterm,xyzmh,vxyzu,rho,iphase,isteps,nptmass,
     &   listpm,spinx,spiny,spinz,angaddx,angaddy,angaddz,
     &   spinadx,spinady,spinadz,ierr)
         IF (ierr.GT.0) STOP 'aborting...'
c         IF (ierr.EQ.-1) THEN
c            PRINT*,'WARNING: small dump file: AM not present'
c            smalldump = .true.
c         ELSE
c            smalldump = .false.
c         ENDIF
         
         IF (nptmass.GE.2) THEN
            PRINT*,'WARNING: more than two sink particles found '
            IF (nptmass.GT.2 .and. 
     &         (listpm(1).ne.i1 .or. listpm(2).ne.i2)) THEN
               PRINT*,'BEWARE: BINARY IDENTITY MAY HAVE CHANGED!'
            ELSE
               PRINT*,'Binary found -- using sink particle properties'
            ENDIF
            i1 = listpm(1)
            i2 = listpm(2)
            dx = xyzmh(1,i1) - xyzmh(1,i2)
            dy = xyzmh(2,i1) - xyzmh(2,i2)
            dz = xyzmh(3,i1) - xyzmh(3,i2)
            dr = SQRT(dx*dx + dy*dy + dz*dz)
            r1 = SQRT(xyzmh(1,i1)**2 + xyzmh(2,i1)**2 + xyzmh(3,i1)**2)
            r2 = SQRT(xyzmh(1,i2)**2 + xyzmh(2,i2)**2 + xyzmh(3,i2)**2)
            
c            dvx = vxyzu(1,i1) - vxyzu(1,i2)
c            dvy = vxyzu(2,i1) - vxyzu(2,i2)
c            dvz = vxyzu(3,i1) - vxyzu(3,i2)
                        
            angle = ACOS(dx/dr)
            angle = MAX(angle,pi-angle)
            WRITE(*,*) 't =',gt,' r1 = ',r1,' r2 = ',r2,' sep = ',dr
            WRITE(55,*) gt,r1,r2,dr,angle
         ELSE
            PRINT*,'Binary not found -- must look at density maxima'
            rhomax1 = 0.5*rhomax1

            DO i=1,npart
               IF (iphase(i).GE.0) THEN
                  IF (rho(i).GT.rhomax1) THEN
c
c--only accept a density maximum that is within a certain radius of previous maximum
c  (relies on using all dumps at once in the correct order)
c
                     dr = SQRT((xyzmh(1,i)-x1)**2 + (xyzmh(2,i)-y1)**2
     &                        +(xyzmh(3,i)-z1)**2)
                     IF (dr/r1 .LT. 0.01 .or. j.eq.1) THEN
                        rhomax1 = rho(i)
                        i1 = i
c                     ELSE
c                        PRINT*,'rejecting ',i,
c     &                         'rho=',rho(i),(ri-r1)/r1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            IF (i1.EQ.0) STOP 'density maximum not found'
            x1 = xyzmh(1,i1)
            y1 = xyzmh(2,i1)
            z1 = xyzmh(3,i1)
            r1 = SQRT(x1**2 + y1**2 + z1**2)
            PRINT*,'Maximum density on particle ',i1,' rho = ',rho(i1),
     &      'r1 = ',r1
c
c           look for density maximum in the opposite hemisphere
c           roughly equidistant from origin (ie r2 = r1 to 10%)
c
            rhomax2 = 0.
            DO i=1,npart
               IF (iphase(i).GE.0) THEN
                  IF (rho(i).GT.rhomax2 .and. i.NE.i1) THEN
                     !--check hemisphere
                     r2dotr1 = xyzmh(1,i)*x1 + xyzmh(2,i)*y1 
     &                       + xyzmh(3,i)*z1
                     IF (r2dotr1 .le.0.) THEN
                        x2 = xyzmh(1,i)
                        y2 = xyzmh(2,i)
                        z2 = xyzmh(3,i)
                        r2 = SQRT(x2**2 + y2**2 + z2**2)
c                        IF (abs(r2 - r1)/r1.LT.0.2) THEN
                           rhomax2 = rho(i)
                           i2 = i
c                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            x2 = xyzmh(1,i2)
            y2 = xyzmh(2,i2)
            z2 = xyzmh(3,i2)
            r2 = SQRT(x2**2 + y2**2 + z2**2)
            dr = SQRT((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            
            PRINT*,'Maximum density on particle ',i2,' rho = ',rho(i2),
     &      'r2 = ',r2

            angle = ACOS((x2-x1)/dr)
            angle = MAX(angle,pi-angle)

            WRITE(*,*) 't =',gt,' r1 = ',r1,' r2 = ',r2,' sep = ',dr,
     &                 'pos angle = ',angle
            WRITE(55,*) gt,r1,r2,dr,angle

         ENDIF        
      ENDDO
c200   CONTINUE
      
      CLOSE(UNIT=55)
      STOP

      END PROGRAM binaryproperties
