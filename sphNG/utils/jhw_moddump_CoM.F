      PROGRAM jhw_moddump_CoM
c***********************************************************
c  Recentres a dumpfile.  Options include moving it to the *
c  CoM of the particles near the densest particle, the CoM *
c  of particles above a density threshold, or to a         *
c  selected sink particle.                                 *
c  in sphNG/, compile as make jhw_moddump_CoM              *
c     author: James Wurster                                *
c***********************************************************

#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/logun'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
      INCLUDE 'COMMONS/vargam'
#endif
#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      PARAMETER (maxfiles = 200)
      PARAMETER (itests   =  11)
      PARAMETER (nmaxx    = 300)

      INTEGER      :: i,j,k,p,istart,istop,iptmax0
      INTEGER      :: iarray(maxfiles),ptest(itests)
      INTEGER      :: indx(idim),list(idim),idsinkarray(iptdim)
      INTEGER      :: ibin(nmaxx,nmaxx,nmaxx)
      REAL         :: kyr,rhothresh_cgs,rmax,xyz_au(3)
      REAL         :: qrat(nmaxx,nmaxx),den(3000,3000)
      REAL         :: q3(nmaxx,nmaxx,nmaxx),std(nmaxx,nmaxx,nmaxx)
      REAL         :: hproj,mi,xcen,ycen,zcen
      REAL*8       :: udistau,dummy(3,imhd2)
      CHARACTER*3  :: prefix,outfix
      CHARACTER*7  :: filein,fileout
      CHARACTER*11 :: fileouta
      LOGICAL      :: dumprange,printascii,printdump,lineone
      LOGICAL      :: untilt,com_all_sinks,slice,is_cluster
      LOGICAL      :: update_sink_list,make_null,make_null_vr
      LOGICAL      :: make_fpfk_neigh,recentre

#ifdef MPIALL
      k = 0
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)
#endif

      iprint     = istdout 
      iread      = 5
      printdump  = .true.
      printascii = .false.
      untilt     = .true.
      isink      = 0
      is0        = 0
      rhothresh_cgs = 1.0d-13
      slice         = .true.
      maxrec        = 100*idim
      rmax_au       = 400. ! max distance to search when untilting disc
      recentre      = .false.  ! manually set the origin
      make_null     = .false.
      make_null_vr  = .false. !Set the velocity to have small vr, else just negative Cartesian values
      make_fpfk_neigh = .false. ! use pre-defined grids and remove outiers

#ifndef MPI
      iproc = 0
#endif

      if (make_null) then
         print*, 'Resetting values to make null file for Splashing'
      endif
      call get_in_give_out(prefix,outfix,.true.,dumprange,
     &                           maxfiles,istart,istop,iarray)
      WRITE(inname,'(2a)')'ifile_',prefix
      CALL options
      if (outfix=="mal") then
         untilt = .false.
         outfix = "all"
      else if (outfix=="msn") then
         untilt = .false.
         outfix = "snk"
      endif
      if (outfix=="all" .or. outfix=="snk") then
         call sinklist_read(iptdim,iptmax,idsinkarray)
         if (outfix=="snk") then
            WRITE (*,*) 'Select sink number (0 for all)'
            READ (iread,*) is0
            if (is0 > iptmax) then
               print*, "That sink does not exist.  Aborting."
               call quit(0)
            endif
         endif
      elseif (outfix=="fil") then
         nptmass_in    = -1
         rhothresh_cgs = 1.0d-14
         rmax_au       = 1.0d30
         slice         = .false.
      elseif (outfix(1:1)=="m") then
         untilt        = .false.
      elseif (outfix(1:1)=="d") then
         nptmass_in    = 0
         untilt        = .false.
      elseif (outfix=="rec") then
         recentre      = .true.
      endif
      if (is_cluster()) then
         ! these values are required due to lower resolution
         print*, "This is cluster data"
         rhothresh_cgs = 1.0d-14
         rmax_au       = 60. ! max distance to search when untilting disc
         slice         = .false.
      endif
      if (recentre) then
         write (*,*) 'Will re-centre on the given input coordinates'
         write (*,*) 'Enter the x-coordinate (au)'
         read (iread,*) xcen_au
         write (*,*) 'Enter the y-coordinate (au)'
         read (iread,*) ycen_au
         write (*,*) 'Enter the z-coordinate (au)'
         read (iread,*) zcen_au
      endif

#ifdef MPI
      mpitype = 'd'
#endif
      xlocationmpi = 0.
      ylocationmpi = 0.
      zlocationmpi = 0.
      varmhd       = 'Brho'
      kyr          = 3.15576d10
      ptest( 1)    = 1000
      ptest( 2)    = 10000
      ptest( 3)    = 100000
      ptest( 4)    = 1000000
      ptest( 5)    = 3000
      ptest( 6)    = 30000
      ptest( 7)    = 300000
      ptest( 8)    = 3000000
      ptest( 9)    = 6000
      ptest(10)    = 60000
      ptest(11)    = 600000
c     Not sure why, but this loop is necessary
      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO

      DO k = istart, istop
         IF (dumprange) THEN
            iloop = k
         ELSE
            iloop = iarray(k)
         ENDIF

c        Read in file
         WRITE(filein,  '(a,I4.4)')prefix,iloop
         print*, 'Reading file ',trim(filein)
c         OPEN (11,FILE=trim(filein),FORM='unformatted', recl=maxrec)
         OPEN (11,FILE=trim(filein),FORM='unformatted')
c         CALL rdump(11,ichkl,1)
         CALL rdump(11,ichkl,0)
          CLOSE (11)
         IF (ichkl.EQ.1) THEN
            WRITE (*,*) 'ERROR reading input file: ichkl ',ichkl
            STOP
         ENDIF
         IF (outfix=="fil") then
            nptmass_in = -1
         ELSEIF (outfix(1:1)=="d") then
            nptmass_in =  0
         ELSE
            nptmass_in = nptmass
         ENDIF
         umass  = umassi
         udist  = udisti
         utime  = utimei
         umagfd = umagfdi
         dtmax  = dtmaxdp
         udistau= udisti/au
         rmax   = rmax_au/udistau
         xcen   = xcen_au/udistau
         ycen   = ycen_au/udistau
         zcen   = zcen_au/udistau
         print*, udist,udisti
         print*, utime,utimei

c        Ensure all the sinks have been read in
         update_sink_list = .false.
         do j = 1,nptmass
            iptmax0 = iptmax
            call sinklist_getid(idummy,iptdim,iptmax,
     &                          iunique(listpm(j)),idsinkarray)
            if (iptmax0.ne.iptmax) update_sink_list = .true.
         enddo
         if (update_sink_list) then
c           Print updated sink list
            call sinklist_write(iptdim,iptmax,idsinkarray)
         endif

c        Determine loops
         isink_in = 0
         if (is0 > 0) then
            do j = 1,iptmax
              if (idsinkarray(is0).EQ.iunique(listpm(j))) isink_in = j
            enddo
         endif
         print*, "active sink ids: local & unique"
         do j = 1,nptmass
              print*, listpm(j),iunique(listpm(j))
         enddo
         if (outfix=="all") then
            js = 0
            je = nptmass
         else if (outfix=="snk") then
            if (isink_in==0) then
               js = 1
               je = nptmass
            else
               js = isink_in
               je = isink_in
            endif
         else
            js = 0
            je = 0
         endif

c        !---Print new dump-file
         do j = js,je
c           Get output filename
            if (j==0) then
               WRITE(fileout, '(a,I4.4)')outfix,iloop
            else
               call sinklist_getid(isink,iptdim,iptmax,
     &              iunique(listpm(j)),idsinkarray)
               WRITE(fileout, '(a,I2.2,I4.4)')'s',isink,iloop
            endif

            if (recentre) then
c              Shift to the input coordinates
               do i = 1,npart
                  xyzmh(1,i) = xyzmh(1,i) - xcen
                  xyzmh(2,i) = xyzmh(2,i) - ycen
                  xyzmh(3,i) = xyzmh(3,i) - zcen
               enddo
            else
c              Shift to the centre of mass
               call centre_of_mass(idim,idim2,iptdim,npart,nptmass_in,
     &                    xyzmh,vxyzu,rho,iphase,listpm,umassi,udisti,
     &                    au,j,xcen,ycen,zcen)
            endif

c           Un-tilt the disc if requested
            if (untilt) then
               dummy = 0.0
               call get_disc_angles(npart,idim,idim2,imhd2,xyzmh,vxyzu,
     &                             Bxyz,dummy,rho,iphase,rhothresh_cgs,
     &                             rmax,slice,umassi,udisti,pi)
               do i = 1,npart
                  Bevolxyz(1:3,i) = Bxyz(1:3,i)/rho(i)
               enddo
            endif
            if (make_null) then
               do i = 1,npart
                  rho(i) = 1.d-30*udisti**3/umassi
                  Bxyz(:,i) = 1.d-20
                  vxyzu(4,i) = epsilon(vxyzu(4,i))
                  Bevolxyz(1:3,i) = Bxyz(1:3,i)/rho(i)
                  if (make_null_vr) then
                     if (xyzmh(1,i) > 0.) then
                        vxyzu(1,i) = -100.d5*utimei/udisti
                     else
                        vxyzu(1,i) =  100.d5*utimei/udisti
                     endif
                     if (xyzmh(2,i) > 0.) then
                        vxyzu(2,i) = -100.d5*utimei/udisti
       	       	     else
                        vxyzu(2,i) =  100.d5*utimei/udisti
       	       	     endif
                     if (xyzmh(3,i) > 0.) then
                        vxyzu(3,i) = -100.d5*utimei/udisti
       	       	     else
                        vxyzu(3,i) =  100.d5*utimei/udisti
       	       	     endif
                  else
                     vxyzu(1:3,i) = -100.d5*utimei/udisti
                  endif
               enddo
            endif

            if (make_fpfk_neigh) then
               if (.false.) then
                  dxz = 5.0/udistau  ! FHC
                  x0 = -130./udistau
                  y0 = -130./udistau
                  z0 = -130./udistau
               else
                  dxz = 0.33/udistau ! SHC
                  x0 = -6./udistau
                  y0 = -6./udistau
                  z0 = -6./udistau
               endif
               nmax = int(-2*x0/dxz)
               if (nmax > nmaxx) then
                  print*, 'recompile with nmaxx > ',nmaxx
                  stop
               endif
               q3   = 0.
               ibin = 0.
               std  = 0.
               do i = 1,npart 
                  ! calcualte Q & bin particles
                  if (iphase(i)==0) then
                     xi = xyzmh(1,i)
                     yi = xyzmh(2,i)
                     zi = xyzmh(3,i)
                     mi = xyzmh(4,i)
                     hi = xyzmh(5,i)
                     Bx = Bxyz(1,i)
                     By = Bxyz(2,i)
                     Bz = Bxyz(3,i)
                     vx = vxyzu(1,i)
                     vy = vxyzu(2,i)
                     vz = vxyzu(3,i)
                     v2 = vx*vx+vy*vy+vz*vz
                   vxyzu(4,i) = 2.*(Bx*Bx+By*By - Bz*(Bx*vx+By*vy)/vz)
     &                       /(rho(i)*v2)
                     if (vxyzu(4,i) < 0.) vxyzu(4,i) = -vxyzu(4,i)
                     if (vxyzu(4,i) > 1.d3 .or. vxyzu(4,i)< 0.1) then
                        iphase(i) = -1
                     else
                        ix = int( (xi - x0) / dxz) + 1
                        iy = int( (yi - y0) / dxz) + 1
                        iz = int( (zi - z0) / dxz) + 1
                        if (0 < ix .and. ix < nmax .and.
     &                      0 < iy .and. iy < nmax .and.
     &                      0 < iz .and. iz < nmax) then
                           q3(ix,iy,iz) = q3(ix,iy,iz)+vxyzu(4,i)
                           ibin(ix,iy,iz) = ibin(ix,iy,iz)+1
                        else
                           iphase(i) = -1
                        endif
                     endif
                   endif
                enddo
                ! calculate average
                do iz = 1,nmax
                  do iy = 1,nmax
                     do ix = 1,nmax
                        if (ibin(ix,iy,iz)> 0) then
                           q3(ix,iy,iz) = q3(ix,iy,iz)/ibin(ix,iy,iz)
                         endif
                      enddo
                   enddo
                enddo
                ! loop thopugh again to calculate standard deviation
                do i = 1,npart 
                  if (iphase(i)==0) then
                     xi = xyzmh(1,i)
                     yi = xyzmh(2,i)
                     zi = xyzmh(3,i)
                     ix = int( (xi - x0) / dxz) + 1
                     iy = int( (yi - y0) / dxz) + 1
                     iz = int( (zi - z0) / dxz) + 1
                     std(ix,iy,iz) = std(ix,iy,iz)
     &                             + (vxyzu(4,i)-q3(ix,iy,iz))**2
                  endif
                enddo
                do iz = 1,nmax
                  do iy = 1,nmax
                     do ix = 1,nmax
                        if (ibin(ix,iy,iz)> 0) then
                           std(ix,iy,iz) = sqrt(std(ix,iy,iz)
     &                                   / ibin(ix,iy,iz))
                         endif
                      enddo
                   enddo
                enddo
                ! final loop to remove all particle outide of 5 standard deviations
                ictr = 0
                do i = 1,npart 
                  if (iphase(i)==0) then
                     xi = xyzmh(1,i)
                     yi = xyzmh(2,i)
                     zi = xyzmh(3,i)
                     ix = int( (xi - x0) / dxz) + 1
                     iy = int( (yi - y0) / dxz) + 1
                     iz = int( (zi - z0) / dxz) + 1
                     if(vxyzu(4,i) > q3(ix,iy,iz)+3.*std(ix,iy,iz).or.
     &                  vxyzu(4,i) < q3(ix,iy,iz)-3.*std(ix,iy,iz))then
                       ictr = ictr + 1
                       iphase(i) = -1
                       print*, 'removing ',vxyzu(4,i),' from ',ix,iy,iz
                     endif
                  endif
               enddo
               print*, 'removed ',ictr
            endif

            if (printdump) then
               print*, 'Writing to output file < ',trim(fileout),' >'
               file1 = trim(fileout)
               IF (iproc.EQ.0) THEN
c                  CALL file  ! this just makes a new filename, which is not what we want
                  OPEN (11,FILE=file1,FORM='unformatted')
               ENDIF
               print *,iproc,': Calling wdump ',npart
               CALL wdump(11)
c              CLOSE(11) ! included in wdump
               print *,iproc,': Done wdump'
            endif
         enddo

c        !---This will print out the particle positions and density to an ascii file
         if (printascii) then
            WRITE(fileouta,'(a,I4.4,a)')outfix,iloop,'.dat'
            open(23,file=trim(fileouta))
            write(23,"('#',6(1x,'[',i2.2,1x,a11,']',2x))")
     &      1,'id',
     &      2,'x (au)',
     &      3,'y (au)',
     &      4,'z (au)',
     &      5,'rho(g/cm^3)',
     &      6,'time (kyr)'

c           Identify and sort the particles
            DO i = 1,npart
               list(i) = i
            ENDDO
c           NOTE: Inside of sphng, use iunique(iorig(i))
c                 in analysis files, use iunique(i)
            CALL indexx(npart,list,real(iunique(1:npart)),indx)

c           Print them out in order 
            lineone = .true.
            ictr = 0
            DO i = 1,npart
               j = list(indx(i))
               p = iunique(j)
               if (iphase(j).EQ.0 .and. p.GT.0 .AND. p.LE.n1) then
                  ictr    = ictr + 1
                  xyz_au  = xyzmh(1:3,abs(j))*udistau
                  rho_cgs = rho(j)*umassi/udisti**3
                  IF (lineone) THEN
                     lineone = .false.
                      write(23,'(i18,1x,5(es18.4,1x))')
     &                iunique(j),xyz_au,rho_cgs,gt*utime/kyr
                  ELSE
                     write(23,'(i18,1x,4(es18.4,1x))')
     &                iunique(j),xyz_au,rho_cgs
                  ENDIF

c                 print values of test particles
                  do q = 1,itests
                    if (p==ptest(q)) call pnt_val(941+q,xyz_au,rho_cgs)
                  enddo

               ENDIF
            END DO
            close(23)
            print *,iproc,': Done ascii.  Npart = ',npart
            write(953,*) k,': Done ascii.  ictr = ',ictr
         ENDIF
      END DO

      END
c----------------------------------------------------------------------!
c     subroutine to print values of test particles
      subroutine pnt_val(iunit,xyz_au,rho_cgs)
      integer :: iunit
      real    :: xyz_au(3),rho_cgs

      write(iunit,'(4(es18.4,1x))') xyz_au, rho_cgs

      end subroutine pnt_val

