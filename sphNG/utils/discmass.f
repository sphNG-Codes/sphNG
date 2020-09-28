      PROGRAM disc1
c************************************************************
c                                                           *
c  This program extracts properties of a disc               *
c  (mass, radius, angular momentum)                         *
c  from a series of dump files                              *
c                                                           *
c  DJP Nov 2006; MRB Sept 2020                              *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'

      INTEGER i,j,k,nfiles,iargc,ierr,ichkl
      INTEGER imax,ioutput,ndisc,ipart,list(idim),listr(idim)

      REAL xcm1,ycm1,zcm1,dx,dy,dz
      REAL vxcm1,vycm1,vzcm1
      REAL pmassi
      REAL discmass,rdisc,dr,totmass
      REAL etot,semimajor,eccentricity
      REAL radius(idim),fracmass,frac
      REAL discr50pc,discr90pc,discr95pc,discr99pc
      
      REAL xyzm(4,2),vxyz(3,2)
      
      CHARACTER*100 filename,fileout
      LOGICAL smalldump
c
c--See how many dump files have been listed (error if none)
c            
      nfiles = iargc()
      IF (nfiles.LE.0) THEN
         STOP 'Usage: massvsangmom dumpfile(s)'
      ENDIF
c
c--Open file for j vs time plots
c
      PRINT*,'0: disc mass, radius etc vs time'
      PRINT*,'1: write disc particles selected from given dump to file'
      WRITE(*,*) 'Please select output:'
      READ*,ioutput
      IF (ioutput.GT.1 .OR. ioutput.LT.0) STOP 'unknown output choice'
c
c--Open file to write sequence of disc parameters from each dump file
c
      IF (ioutput.EQ.0) THEN
         fileout = 'disc.out'
         PRINT*,' opening ',fileout
         OPEN(UNIT=55,FILE=fileout,STATUS='replace',FORM='formatted')
      ENDIF

      DO j = 1, nfiles
c
c--Get next input file name from argument list of command
c
         CALL getarg(j,filename)

         PRINT*,'opening ',filename
c
c--Read dump file
c
         OPEN (UNIT = 11, FILE = filename, FORM = 'unformatted')
         CALL rdump(11, ichkl, 0)
         CLOSE (11)
         IF (ichkl.NE.0) STOP 'Dump file error - aborting...'
c
c--Set the properties of the effective central mass
c     If no sink particles, give error message.
c     If only one sink particle, use its xyz location and mass,
c       and its vx, vy, vz velocity.
c     If more than one sink particle, take centre of mass position,
c       total sink particle mass, and centre of mass velocity.
c
         IF (nptmass.EQ.1) THEN
            PRINT*,'Point mass found: calculating disc properties'
            i = listpm(1)
            xyzm(1:4,1) = xyzmh(1:4,i)
            vxyz(1:3,1) = vxyzu(1:3,i)
         ELSEIF (nptmass.EQ.0) THEN
            PRINT*,'No point masses - undefined central mass'
            STOP
         ELSE
            PRINT*,'nptmass>1: using centre of mass of point masses'
            xcm1 = 0.
            ycm1 = 0.
            zcm1 = 0.
            vxcm1 = 0.
            vycm1 = 0.
            vzcm1 = 0.
            totmass = 0.
            DO k = 1, nptmass
               i = listpm(k)
               pmassi = xyzmh(4,i)
               xcm1 = xcm1 + pmassi*xyzmh(1,i)
               ycm1 = ycm1 + pmassi*xyzmh(2,i)
               zcm1 = zcm1 + pmassi*xyzmh(3,i)
               totmass = totmass + pmassi
               vxcm1 = vxcm1 + pmassi*vxyzu(1,i)
               vycm1 = vycm1 + pmassi*vxyzu(2,i)
               vzcm1 = vzcm1 + pmassi*vxyzu(3,i)
            END DO
            xyzm(1,1) = xcm1/totmass
            xyzm(2,1) = ycm1/totmass
            xyzm(3,1) = zcm1/totmass
            xyzm(4,1) = totmass
            vxyz(1,1) = vxcm1/totmass
            vxyz(2,1) = vycm1/totmass
            vxyz(3,1) = vzcm1/totmass
         ENDIF
c
c--Zero quantities that need to be determined for each dump file
c
         rdisc = 0.
         ndisc = 0
         discr50pc = 0.
         discr90pc = 0.
         discr95pc = 0.
         discr99pc = 0.
c
c--If it has been requested to write out the particles that are selected
c     then open the file that will contain the properties of these
c     particles (e.g. positions, masses, etc).
c
         IF (ioutput.EQ.1) OPEN(UNIT=10,FILE=trim(filename)//'.disc',
     &        STATUS='replace',FORM='formatted')
c
c--Zero the number of particles selected and the associated disc mass.
c     Then loop over all particles to select disc particles.
c
         ndisc = 0
         discmass = 0.
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN  ! if particle is gas
c
c--These lines load the x,y,z positions and mass of the gas particle
c     into xyzm(*,2) and its x,y,z velocities into vxyz(*,2) 
c
               xyzm(1:4,2) = xyzmh(1:4,i)
               vxyz(1:3,2) = vxyzu(1:3,i)
c
c--Find the properties of the orbit that the gas particle has around the
c     central object.
c
               CALL find_orbit(xyzm,vxyz,etot,semimajor,eccentricity)
c
c--Take a particle to be a disc particle if its orbit has an eccentricity
c     e<0.3.
c
               IF (eccentricity .LT. 0.3) THEN
c
c--Compute its radius from the central mass, and the maximum radius
c
                  dx = xyzmh(1,i) - xyzm(1,1)
                  dy = xyzmh(2,i) - xyzm(2,1)
                  dz = xyzmh(3,i) - xyzm(3,1)
                  dr = sqrt(dx**2 + dy**2 + dz**2)
                  rdisc = MAX(rdisc,dr)
c
c--Add to a list of particles so that particles can be sorted in radius
c
                  ndisc = ndisc + 1
                  list(ndisc) = i
                  radius(ndisc) = dr
c
c--Keep track of total mass of particles that have been selected
c
                  discmass = discmass + xyzmh(4,i)
c
c--Write out the particle properties if requested (ioutput=1)
c
                  IF (ioutput.EQ.1) THEN
                     WRITE(10,99000) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i),
     &                     xyzmh(4,i),xyzmh(5,i),rho(i)*udens
99000                FORMAT(6(1PE12.5,1X))
                  ENDIF
               ENDIF
            ENDIF
         END DO
c
c--Close output file listing selected particles
c
         IF (ioutput.EQ.1) CLOSE(10)
c
c--Sort disc particles by radius from central mass
c
         CALL indexx2(ndisc,radius,listr)
c
c--Calculate disc radii containing half/90%/95%/99% of total disc mass
c
         fracmass = 0.
         DO ipart = 1, ndisc
            i = listr(ipart)
            fracmass = fracmass + xyzmh(4,list(i))
            frac = fracmass/discmass
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
         END DO
c
c--Write disc parameters to a file (unit 55), and also to the screen
c
         IF (ioutput.EQ.0) THEN
            WRITE(55,99001) gt,discmass,discmass/totmass,discr50pc,
     &                 discr90pc,discr95pc,discr99pc,rdisc
99001       FORMAT(8(1PE12.5,1X))
         ENDIF
         WRITE(*,*) ' npart in disc = ',ndisc
         WRITE(*,*) ' point mass / (disc + pt mass) = ',
     &               xyzmh(4,imax)/(discmass + xyzmh(4,imax))
         WRITE(*,*) 't = ',gt
         WRITE(*,*)'mass in disc = ',discmass,
     &        discmass/totmass
         WRITE(*,*)
     &        'mtot=',totmass,'disc radius=',discr50pc,discr90pc,
     &        discr95pc,discr99pc,rdisc
c
c--Loop again for next dump file
c
      END DO
c
c--Close summary output file
c
      IF (ioutput.EQ.0) CLOSE(55)
               
      END PROGRAM disc1
c
c========================================================================
c
      SUBROUTINE find_orbit(xyzm,vxyz,etot,semimajor,eccentricity)
c
c--Assuming object 2 is in orbit around object 1, determine orbital
c     parameters.  Object 1 positions and mass are stored in xyzm(*,1)
c     and Object 2 in xyzm(*,2).  Similarly for vxyz().
c
c     NOTE: Only determines values for e<=1 (not hyperbolic orbits).
c
      IMPLICIT NONE

      REAL xyzm(4,2),vxyz(3,2)
      REAL etot,semimajor,eccentricity

      REAL totalmass,rx,ry,rz,vxdiff,vydiff,vzdiff,r2
      REAL ekin,epot,dangx,dangy,dangz,dang2,value
c
c--Set total mass to consider whether two objects are bound or not
c
      totalmass = xyzm(4,1) + xyzm(4,2)
      IF (totalmass.LE.0.) THEN
         PRINT *,'ERROR - find_orbit'
         STOP
      ENDIF
c
c--Find relative positions and velocities
c
      rx = xyzm(1,1) - xyzm(1,2)
      ry = xyzm(2,1) - xyzm(2,2)
      rz = xyzm(3,1) - xyzm(3,2)
      vxdiff = vxyz(1,1) - vxyz(1,2)
      vydiff = vxyz(2,1) - vxyz(2,2)
      vzdiff = vxyz(3,1) - vxyz(3,2)
      r2 = (rx**2 + ry**2 + rz**2)
c
c--Energies and angular momenta
c
      ekin = 0.5*(vxdiff**2 + vydiff**2 + vzdiff**2)
      epot = -totalmass/SQRT(r2)
      etot = ekin + epot
      dangx = vydiff*rz - vzdiff*ry
      dangy = vzdiff*rx - vxdiff*rz
      dangz = vxdiff*ry - vydiff*rx
      dang2 = (dangx**2 + dangy**2 + dangz**2)
c
c--Orbital parameters
c
      IF (etot.LE.0.) THEN
         semimajor = - totalmass/(2.0*etot)
         value = 1.0 - dang2/totalmass/semimajor
         IF (value.GE.0.0) THEN
            eccentricity = SQRT(value)
         ELSE
            eccentricity = 1.0
         ENDIF
      ELSE
         semimajor = 1.0E+30
         eccentricity = 1.0
      ENDIF

      RETURN
      END
c
c========================================================================
c
      SUBROUTINE quit
      STOP
      END
c
c========================================================================
c
      SUBROUTINE endrun
      CALL quit
      END

