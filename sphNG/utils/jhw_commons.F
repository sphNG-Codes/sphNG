c This contains generic functions to be used across many jhw* analysis routines
c  I/O: Input/outputs
c  CoM: Shift to the centre of mass
c  Miscellaneous required to run
c----------------------------------------------------------------------!
c     Gets common in/out from the interactive sessions
      subroutine get_in_give_out(prefix,outfix,need_out,dumprange,
     &                           maxfiles,istart,istop,iarray)
      character(len=3) :: prefix,outfix
      logical :: need_out,dumprange
      integer :: maxfiles,istart,istop,iarray(maxfiles)
      integer :: iread,numfulldumps
      
      iread = 5
      WRITE (*,*) 'Enter the input file prefix'
      READ (iread,*) prefix

      IF (need_out) THEN
         print*, 'if moddump_CoM:'
         print*, ' all/snk: provides options to shift to sinks'
         print*, ' fil:     will shift and rotate primary filament'
         print*, ' m??:     will not untilt model'
         print*, ' d??:     shift to densest gas particle; untilted'
         print*, ' else:    shift to CoM sinks or densest gas particle'
         WRITE (*,*) 'Enter the output file prefix now'
         READ (iread,*) outfix
      ENDIF

      WRITE (*,*) 'How many full dumps (< 1 to select range)?'
      READ (iread,*) numfulldumps
      IF (numfulldumps < 1) THEN
         dumprange = .true.
      ELSE
         dumprange = .false.
      ENDIF
      IF (numfulldumps.GT.maxfiles) THEN
         WRITE (*,*) 'Filename array too small'
         STOP
      ENDIF

      IF (dumprange) THEN
         WRITE (*,*) 'Enter first dump to convert'
         READ (iread,*) istart
         WRITE (*,*) 'Enter final dump to convert'
         READ (iread,*) istop
      ELSE
         WRITE (*,*) 'Enter dump numbers'
         DO i = 1, numfulldumps
            READ (iread,*) iarray(i)
         END DO
         istart = 1
         istop  = numfulldumps
      ENDIF

      end subroutine get_in_give_out
c----------------------------------------------------------------------!
c     Will create an incremental filename
      subroutine get_filename(i,filename,prefix,suffix)
      integer           :: i
      character(len=*)  :: filename,prefix,suffix
      logical           :: iexists

      iexists = .true.
      i = -1
      do while (iexists .and. i < 98)
         i = i + 1
         write(filename,'(a,I2.2,a)') trim(prefix),i,trim(suffix)
         inquire(file=trim(filename),exist=iexists)
      enddo
      if (iexists) then
         print*,'Cannot determine unique output filename.  Aborting.'
         stop
      endif

      end subroutine get_filename
c----------------------------------------------------------------------!
c     Determine if dumps must be from splash.filenames; if so, get the list
      subroutine get_splashdumps(ks,ke,ns,maxn,prefix,splashnames)
      integer       :: ks,ke,ns,maxn,iio
      character*3   :: prefix
      character*256 :: splashnames(maxn),splashname,splashin,folder
      logical       :: fexist

      ns = 0
      splashnames = ""
      if (ks < 0 .or. ke < 0) then
         ks  = abs(ks)
         ke  = abs(ke)
         iio = 0
         inquire(file='Folder.dat',exist=fexist)
         if (fexist) then
            open(unit=16,file="Folder.dat")
            read(16,'(a)') folder
            close(16)
         else
            folder = prefix
         endif
         splashin = 'splashList.filenames'
         inquire(file=trim(splashin),exist=fexist)
         if (.not. fexist) then
            splashin = 'splash.filenames'
            inquire(file=trim(splashin),exist=fexist)
         endif
         if (fexist) then
            print*, 'Reading from ',trim(splashin)
            open(unit=16,file=trim(splashin))
            do while (iio.EQ.0)
               read(16,'(a)',iostat=iio) splashname
               if (iio.eq.0 .and. 
     &         index( trim(splashname),trim(folder) ) > 0) then
                  ns = ns + 1
                  splashnames(ns) = trim(splashname)
               endif
            enddo
            close(16)
         endif
      endif

      end subroutine get_splashdumps
c----------------------------------------------------------------------!
c     Determine if the file is on the splash list
      function is_splash_dump(ns,maxn,splashnames,dump)
      integer      :: ns,maxn,j
      character*7  :: dump
      character*256:: splashnames(maxn)
      logical      :: kexists,is_splash_dump

      kexists = .false.
      j = 1
      do while(j.le.ns .and. .not.kexists)
        ! if dump is in splashnames(j)
        if (index( trim(splashnames(j)),trim(dump) )> 0) kexists=.true.
        j = j + 1
      enddo
      is_splash_dump = kexists

      end function is_splash_dump
c----------------------------------------------------------------------!
c     Moves the particles to the centre of gas mass or to the sink particle
      subroutine centre_of_mass(idim,idim2,iptdim,npart,nptmass,xyzmh,
     &              vxyzu,rho,iphase,listpm,umassi,udisti,au,isink_in,
     &              xcen,ycen,zcen)
      integer   :: idim,iptdim,npart,nptmass,listpm(iptdim)
      integer*1 :: iphase(idim2)
      real      :: mcen,xyzmh(5,idim),vxyzu(4,idim)
      real*4    :: rho(idim2)
      real*8    :: au,umassi,udisti,udistau,udensity
      integer   :: i,ipt,irhomax,ncen,isink,isink_in
      real      :: xcen,ycen,zcen,rcen,rhothresh
      real      :: vxcen,vycen,vzcen
c
      print*, 'Modifying CoM',isink_in
      ncen  = 0
      xcen  = 0.
      ycen  = 0.
      zcen  = 0.
      rcen  = 0.
      vxcen = 0.
      vycen = 0.
      vzcen = 0.
      mcen  = 0.
      isink = isink_in
      udensity = umassi/udisti**3
      udistau  = udisti/au
      IF (isink==0 .and. nptmass==1) isink = 1

      IF (isink > 0) THEN
         ipt    = listpm(isink)
         print*, "Shifting to sink particle ",isink
         ncen   = 0
         xcen   = xyzmh(1,ipt)
         ycen   = xyzmh(2,ipt)
         zcen   = xyzmh(3,ipt)
         vxcen  = vxyzu(1,ipt)
         vycen  = vxyzu(2,ipt)
         vzcen  = vxyzu(3,ipt)
         irhomax = ipt
      ELSE IF (nptmass > 1) THEN
         print*, "Shifting to CoM of massive sink particles"
         DO i = 1,nptmass
            ipt  = listpm(i)
            mcen = max(mcen,xyzmh(4,ipt))
         ENDDO
         DO i = 1,nptmass
            ipt = listpm(i)
            IF (xyzmh(4,ipt).GT.0.66*mcen) THEN
               ncen  = ncen + 1
               xcen  = xcen + xyzmh(1,ipt)*xyzmh(4,ipt)
               ycen  = ycen + xyzmh(2,ipt)*xyzmh(4,ipt)
               zcen  = zcen + xyzmh(3,ipt)*xyzmh(4,ipt)
               rcen  = rcen + xyzmh(4,ipt)
            ENDIF
         ENDDO
         call normalise_cen(xcen,ycen,zcen,rcen)
      ELSE IF (nptmass== -1) THEN
         print*, "Shifting to CoM of dense gas"
         rhothresh = 1.0d-14/udensity
         DO i = 1, npart
            IF (iphase(i).EQ.0 .and. rho(i).GE.rhothresh) THEN
               ncen = ncen + 1
               xcen = xcen + xyzmh(1,i)*rho(i)
               ycen = ycen + xyzmh(2,i)*rho(i)
               zcen = zcen + xyzmh(3,i)*rho(i)
               rcen = rcen + rho(i)
            ENDIF
         END DO
         call normalise_cen(xcen,ycen,zcen,rcen)
      ELSE
         print*, "Shifting to CoM of gas particle"
         irhomax = 1
         DO i = 1, npart
            IF (iphase(i).EQ.0 .and. rho(i).GT.rho(irhomax)) irhomax = i
         END DO
         WRITE (*,'(a,I8)')      'irhomax:', irhomax
         WRITE (*,'(a,2Es16.4)') '  density (code,cgs):',
     &               rho(irhomax),rho(irhomax)*udensity
         WRITE (*,'(a,3Es16.4)') '  location (code): ',
     &               xyzmh(1:3,irhomax)
         WRITE (*,'(a,3Es16.4)') '  location (au):   ',
     &               xyzmh(1:3,irhomax)*udistau
         rhothresh = 0.98*rho(irhomax)
         r2thresh  = (16.0*xyzmh(5,irhomax))**2
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               r2 = (xyzmh(1,i)-xyzmh(1,irhomax))**2
     &            + (xyzmh(2,i)-xyzmh(2,irhomax))**2
     &            + (xyzmh(3,i)-xyzmh(3,irhomax))**2
               IF (r2.LT.r2thresh .and. rho(i).GE.rhothresh) THEN
                  ncen = ncen + 1
                  xcen = xcen + xyzmh(1,i)*rho(i)
                  ycen = ycen + xyzmh(2,i)*rho(i)
                  zcen = zcen + xyzmh(3,i)*rho(i)
                  rcen = rcen + rho(i)
               ENDIF
            ENDIF
         END DO
         IF (rho(irhomax)*udensity.LT.1.d-14) THEN
            rcen = 0.0
            write(*,'(a,Es16.4,a)') 
     &      "WARNING! rho_max = ", rho(irhomax)*udensity,
     &      " < 1d-14g/cm^3, therefore not moving to CoM"
         ENDIF
         call normalise_cen(xcen,ycen,zcen,rcen)
      ENDIF

      WRITE (*,'(a,3Es16.4,2I8)') '  new centre:      ',
     &    xcen,ycen,zcen,ncen,isink
      WRITE (*,'(a,3Es16.4,2I8)') '  new centre (au): ',
     &    xcen*udistau,ycen*udistau,zcen*udistau,ncen,isink

      rhothresh = 1.0d-13/udensity
      DO i = 1, npart
         IF (iphase(i).EQ.0 .and. rho(i).GE.rhothresh) THEN
            ncen  = ncen  + 1
            vxcen = vxcen + vxyzu(1,i)
            vycen = vycen + vxyzu(2,i)
            vzcen = vzcen + vxyzu(3,i)
         ENDIF
      ENDDO
c
ccc      ncen = 0
      IF (ncen > 0) THEN
         print*, 'ncen = ',ncen
         print*, 'vcen = ',vxcen/ncen,vycen/ncen,vzcen/ncen
      ENDIF
      DO i = 1, npart
         xyzmh(1,i) = xyzmh(1,i) - xcen
         xyzmh(2,i) = xyzmh(2,i) - ycen
         xyzmh(3,i) = xyzmh(3,i) - zcen
         IF (ncen > 0) THEN
            vxyzu(1,i) = vxyzu(1,i) - vxcen/ncen
            vxyzu(2,i) = vxyzu(2,i) - vycen/ncen
            vxyzu(3,i) = vxyzu(3,i) - vzcen/ncen
         ENDIF
      END DO

      end subroutine
c----------------------------------------------------------------------!
      subroutine normalise_cen(xcen,ycen,zcen,rcen)
      real :: xcen,ycen,zcen,rcen

      if (rcen > 0.0) then
         xcen = xcen/rcen
         ycen = ycen/rcen
         zcen = zcen/rcen
      else
         xcen = 0.
         ycen = 0.
         zcen = 0.
      endif

      end subroutine
c----------------------------------------------------------------------!
c     Calculate the rotation angles of the disc
      SUBROUTINE get_disc_angles(npart,idim,idim2,imhd2,xyzmh,vxyzu,
     &                           Bxyz,jcur,rho,iphase,
     &                           rhothresh_cgs,rmax,slice,
     &                           umassi,udisti,pi)
      integer   :: npart,idim,idim2,imhd2,i,j,p
      real      :: xyzmh(5,idim),vxyzu(4,idim),Bxyz(3,imhd2),
     &             rad(0:idim)
      real      :: rhothresh_cgs,rmax,pi
      real      :: jcur(3,imhd2)
      real*4    :: rho(idim2)
      real*8    :: umassi,udisti,udensity
      integer*1 :: iphase(idim2)
      logical   :: slice,include_part
      DIMENSION indx(idim),list(idim)

c     Tilt the disc (for testing purposes only)
      if (.false.) then
         aalpha = 73.0
         agamma = 32.0
         cosa   = cos(aalpha*pi/180.0)
         sina   = sin(aalpha*pi/180.0)
         cosg   = cos(agamma*pi/180.0)
         sing   = sin(agamma*pi/180.0)
         print*,aalpha,cosa,sina
         print*,agamma,cosg,sing
         do i = 1,npart
            call rotate_vec(xyzmh(1:3,i),cosa,sina,cosg,sing)
            call rotate_vec(vxyzu(1:3,i),cosa,sina,cosg,sing)
            call rotate_vec(Bxyz (1:3,i),cosa,sina,cosg,sing)
            call rotate_vec(jcur (1:3,i),cosa,sina,cosg,sing)
         enddo
      endif

c     Align the disc
      udensity  = umassi/udisti**3
      rhothresh = rhothresh_cgs/udensity
      print*, 'Beging de-rotation'
c     sort the particles by radius from the centre.  Make list of only the dense particles
      ik  = 0
      rad = 0.0
      do i = 1,npart
         if (iphase(i).EQ.0 .and. rho(i).GT.rhothresh) then
            ik      = ik + 1
            rad(ik) = sqrt(xyzmh(1,i)**2+xyzmh(2,i)**2+xyzmh(3,i)**2)
            list(ik)= i
            if (rad(ik) > rmax) ik = ik - 1
         endif
      enddo
      CALL indexx(ik, list, rad, indx)
      print*, 'There are ',ik,' particles that meet the criteria'

      ang1 = 1.0
      ang2 = 1.0
      ictr = 0
      angthresh = 0.1
      do while ((ik.ge.1 ).and. 
     &          (abs(ang1) > angthresh .or. abs(ang2) > angthresh) )
         ictr = ictr + 1
         do p = 1,2
            angave = 0.
            nctr   = 0
            ! These particles are guaranteed to be acticve gas particles that are dense enough
            im1 = list(indx(1))
            do j = 2,ik
               i = list(indx(j))
               if (p.eq.1) then
                  axis1 = xyzmh(1,i)
                  axis2 = xyzmh(2,i)
               elseif (p.eq.2) then
                  axis1 = xyzmh(2,i)
                  axis2 = xyzmh(1,i)
               endif
               if (slice) then
                  if (2.0*xyzmh(5,i).GE.abs(axis2) .AND. 
     &            2.0*(xyzmh(5,i)+xyzmh(5,im1)).GE.rad(i)-rad(im1).AND.
     &            abs(axis1).GT.0.0) THEN
                     include_part = .true.
                  else
                     include_part = .false.
                  endif
               else
                  include_part = .true.
               endif
               if (include_part .and. abs(axis1).GT.0.0) then
                  angave = angave + atan(xyzmh(3,i)/axis1)
                  nctr   = nctr + 1
               endif
               im1 = i
            enddo
            if (nctr > 0) then
               aalpha = angave/nctr
            else
               aalpha = 0.0
            endif
            if (p.eq.1) then
               ang1 = aalpha*180.0/pi
               cosa = cos(aalpha)
               sina = sin(aalpha)
               do i = 1,npart
                  call update_vec(xyzmh(1,i),xyzmh(3,i),cosa,sina)
                  call update_vec(vxyzu(1,i),vxyzu(3,i),cosa,sina)
                  call update_vec( Bxyz(1,i), Bxyz(3,i),cosa,sina)
                  call update_vec( jcur(1,i), jcur(3,i),cosa,sina)
               enddo
            elseif (p.eq.2) then
               ang2 = aalpha*180.0/pi
               cosa = cos(aalpha)
               sina = sin(aalpha)
               do i = 1,npart
                  call update_vec(xyzmh(2,i),xyzmh(3,i),cosa,sina)
                  call update_vec(vxyzu(2,i),vxyzu(3,i),cosa,sina)
                  call update_vec( Bxyz(2,i), Bxyz(3,i),cosa,sina)
                  call update_vec( jcur(2,i), jcur(3,i),cosa,sina)
               enddo
            endif
         enddo
         write(*,'(a,I4,2F10.4,I8)')
     &       "iteration, rotation angles (degrees), particles: ",
     &       ictr,ang1,ang2,nctr
      enddo
      print*, 'Required ',ictr,' iterations.'

      end subroutine
c----------------------------------------------------------------------!
      subroutine rotate_vec(xyz,cosa,sina,cosg,sing)
      real :: xyz(3),xyztmp(3),cosa,sina,cosg,sing

      xyztmp(1) =  xyz(1)*cosg      + xyz(2)*sing
      xyztmp(2) = -xyz(1)*sing*cosa + xyz(2)*cosg*cosa + xyz(3)*sina
      xyztmp(3) =  xyz(1)*sing*sina - xyz(2)*cosg*sina + xyz(3)*cosa
      xyz       =  xyztmp

      end subroutine
c----------------------------------------------------------------------!
      subroutine update_vec(r1,r2,cosa,sina)
      real :: r1,r2,r1tmp,r2tmp,cosa,sina

      r1tmp =  r1*cosa + r2*sina
      r2tmp = -r1*sina + r2*cosa
      r1 = r1tmp
      r2 = r2tmp

      end subroutine
c----------------------------------------------------------------------!
c     Read in the SinkList:
      subroutine sinklist_read(iptdim,iptmax,idsinkarray)
      integer :: i,io,iptmax,iptdim,itmp
      integer :: idsinkarray(iptdim)
      logical :: kexists
      character*124 :: cdummy

      idsinkarray = 0
      iptmax      = 0
      inquire(file='SinkList.log',exist=kexists)
      if (kexists) then
         open(unit=11,file='SinkList.log')
         read(11,*) cdummy
         io = 0
         iptmax = 0
         do while (io==0)
            read(11,*,iostat=io) i, itmp
            if (io.eq.0) then
               if (i > iptdim) then
                  print*, 'iptdim=',iptdim,' is too small! abort!'
                  call quit(0)
               endif
               iptmax = i
               idsinkarray(iptmax) = itmp
            endif
         enddo
         close(11)
      endif

      print*, "Number of sink id's read in: ",iptmax
      if (iptmax.GT.0) then
         print*, "Sink id's:"
         do i = 1,iptmax
            print*, i,idsinkarray(i)
         enddo
      endif

      end subroutine sinklist_read
c----------------------------------------------------------------------!
c     Write the SinkList:
      subroutine sinklist_write(iptdim,iptmax,idsinkarray)
      integer :: iptmax,iptdim,idsinkarray(iptdim)

      if (iptmax > 0) then
         open(unit=11,file='SinkList.log')
         write(11,'(a)') 'Sink id in output names, Sink particle id'
         do i = 1,iptmax
            write(11, '(2I10)') i, idsinkarray(i)
         enddo
         close(11)
      endif

      end subroutine sinklist_write
c----------------------------------------------------------------------!
c     Get matching sink id
      subroutine sinklist_getid(isink,iptdim,iptmax,idsink,idsinkarray)
      integer :: j,iptmax,iptdim,idsink,idsinkarray(iptdim)

      isink   = -1
      do j = 1,iptmax
         if (idsinkarray(j).EQ.idsink) isink = j
      enddo
      if (isink==-1) then
         iptmax = iptmax + 1
         isink  = iptmax
         idsinkarray(iptmax) = idsink
         print*, "Appending new sink; ids: ",iptmax,isink,idsink
      endif

      end subroutine sinklist_getid
c----------------------------------------------------------------------!
c     From the Properties file, determine if this is a cluster or no
      function is_cluster()
      logical            :: is_cluster,iexists,foundval
      character(len=128) :: cline

      is_cluster = .false.
      foundval   = .false.
      inquire(file="Properties",exist=iexists)
      if (iexists) then
         open(unit=1,file="Properties")
         iio = 0
         do while(iio==0)
            read(1,*,iostat=iio) cline
            iindx = index(cline,'is_cluster')
            if (iindx > 0) then
               foundval = .true.
               iindx = index(cline,'is_cluster=yes')
               if (iindx > 0) is_cluster = .true.
            endif
         enddo
         close(1)
      else
         print*, 'Properties file not found.  Aborting.'
         call quit(0)
      endif
      if (.not. foundval) then
         print*, 'Unable to determine if this is a cluster or not.'
         print*, 'aborting'
         call quit(0)
      endif

      end function is_cluster
c----------------------------------------------------------------------!
c     Miscellaneous subroutines originally at the end of profile.F
      SUBROUTINE quit(i)
      INCLUDE 'idim'
#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      CALL MPI_FINALIZE(ierr)
#endif
      STOP
      END

      SUBROUTINE endrun
      CALL quit(0)
      END

      SUBROUTINE ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)
      STOP
      END

      SUBROUTINE labrun
      RETURN
      END
