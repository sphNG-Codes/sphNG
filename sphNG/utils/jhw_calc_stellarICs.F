      PROGRAM jhw_calc_stellarICs
c***********************************************************
c  This is the first of two programmes required to         *
c  generate initial conditions for                         *
c  Wurster & Rowan (2023a,b); details given in WR2023a.    *
c  This programme craetes the 'associated list'.           *
c  This programme will read a series of dump files         *
c  (e.g., from Wurster, Bate & Price 2019) moving          *
c  backwards in time to create a list of particles         *
c  associated with each sink.  There are several options   *
c  included below the variable declarations; by default,   *
c  we export 4 files per sink, where each file has a       *
c  different extraction density.                           *
c  in sphNG/, compile as make jhw_calc_stellarICs          *
c     author: James Wurster                                *
c***********************************************************
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/sort'
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      INTEGER, PARAMETER :: maxfiles = 512
      INTEGER, PARAMETER :: nden     =   3  ! number of additional density outputs
      INTEGER      :: i,k,j,p,iloop,nneigh,nneigh_all,pstart,pend,ipt
      INTEGER      :: iarray(maxfiles),iouti(2,iptdim)
      INTEGER      :: isink_status(iptdim),idsinkarray(iptdim)
      INTEGER      :: icore(idim,iptdim),sinkformk(iptdim)
      INTEGER      :: icore0(idim,iptdim),itag(idim)
      INTEGER      :: icoresink(iptdim,iptdim),list(idim),indx(idim)
      integer*8    :: iuj,iui
      LOGICAL      :: kexists,dumprange
      real         :: xi,yi,zi,mi,xj,yj,zj,mj,mtot,mtot_all,Ecoef
      real         :: rthresh_au,rthresh,rthresh2,dt,gt_prev,sinkage
      real         :: sinkformt(iptdim),sinkfinal(3,iptdim)
      real         :: epot_array(iptdim)
      real         :: msink(2,iptdim),mcore(iptdim),mtmp(iptdim)
      real         :: uj,Tj,xyzmhj(3),vxyzj(3),Bxyzj(3)
      real         :: tcoresink(iptdim,iptdim),ttag(idim),dr(idim)
      real         :: tcore(idim,iptdim)
      real         :: rhothresh4(nden,iptdim)
      real*8       :: udensity,udist,umass,utime,umagfd,uvel,
     &                uergg,udistau
      real*4       :: rhoj
      character* 3 :: prefix,outfix
      character* 7 :: dump
      character*22 :: fileout
      character*25 :: fileoutAll
      logical      :: firstcall,debug,keep_searching,strict,most_bound
      logical      :: is_associ,is_assoc,is_bound,is_final
      logical      :: check_bound,j0,incl_bound_sinks,test_Emag
      logical      :: print_now,print_evol,print_evolAF,print_all,dave
      logical      :: modV,d4,mass_criteria,in_sink
c
c--read options
c
      call get_in_give_out(prefix,outfix,.false.,dumprange,
     &                           maxfiles,kstart,kstop,iarray)
c
c--initialise parameters (feel free to modify!)
c
      rthresh_au    =  1000.0 !250.   ! Outer radius where particles will be check for density and tighly bound
      rsink_au      =   0.5   ! Sink radius
      runcond_au    =  50.0   ! Inner radius where particles will only be checked for boundness & not density
      rhothresh_cgs = 1.0d-17 ! as soon as the max or average density drops below this value, we will call it a core!
      rhodense_cgs  = 1.0d-18 ! 1.0d-17 ! 1.0d-13 ! will only consider particles greater than this denstiy
      eccentricity0 = 0.5     !0.3     ! will only consider particles with less than this eccentricity
      check_bound   = .true.  ! only record particles that are bound
      print_evol    = .true.  ! print all the particles at all times
      print_evolAF  = .true.  ! print all the particles within the radius  (for *final* time if print_evol)
      print_all     = .false.  ! print all the particles within the radius (for *all*   time if print_evol)
      debug         = .false.  ! print additional debugging statements
      strict        = .false.   ! will only consider stronly bound particles
      most_bound    = .false.   ! at each step, will only associate a particle with the sink to which it is most bound
      incl_bound_sinks = .false. ! will determine what sinks are bound to one another
      dave          = .false.   ! if true, print both when average & max rho drop below rhothresh; else just dmax
      d4            = .true.    ! if true, will test 4 densities
      modV          = .false.    ! modify particle velocity to account for magnetosonic waves
      test_Emag     = .false.  ! check boundness using both Ekin & Emag individually
      mass_criteria = .true.   ! enforce a mass criteria

      print_evol    = .false.  ! print all the particles at all times
      print_evolAF  = .false.  ! print all the particles within the radius  (for *final* time if print_evol)
      print_all     = .false.  ! print all the particles within the radius (for *all* time if print_evol; at end time is not print_evol)
      !ifm          = 2         ! if mass_criteria, then will ensure the clump has at least the amount of mass equal to the final sink mass
      ifm           = 3         ! if mass_criteria, then will ensure the clump has at least the amount of mass equal to the final sink mass+final gas mass
c
c--initialise parameters (do not modify!)
c
      varmhd    = 'Brho'
      encal     = 'r'
      maxrec    = 100*idim
      firstcall = .true.
      icoresink = 0
      icore     = 0
      tcoresink = 0.
      tcore     = 0.
      iout      = 0
      iouti     = 0
      iptmax    = 0
      sinkfinal = -1.
      isink_status = 0     ! 0 = no sink
                           ! 1 = sink is present
                           ! 2 = we are at a time before the sink's birth
                           ! 3 = we are at a time before the core's birth
      if (strict .or. most_bound) check_bound = .true.
      if (most_bound) incl_bound_sinks = .false.
      if (strict) then ! bound if Ecoef*Ekin+Epot < 0
         Ecoef = 2.
      else
         Ecoef = 1.
      endif
      pstart = 223         ! id for the tagged particle output

      if (incl_bound_sinks) then
         print*, '*************************************************'
         print*, 'If A & B are bound B & C are bound,'
         print*, 'it does not guarantee that A & C are bound.'
         print*, 'Therefore perfect symmetry will never be reached.'
         print*, '*************************************************'
      endif
      print*, '*************************************************'
      print*, 'LOGICAL FLAGS:'
      print*, '  debug = ',debug
      print*, '  strict = ',strict
      print*, '  most_bound = ',most_bound
      print*, '  check_bound = ',check_bound
      print*, '  incl_bound_sinks = ',incl_bound_sinks
      print*, '  test_Emag = ',test_Emag
      print*, '  print_evol = ',print_evol
      print*, '  print_evolAF = ',print_evolAF
      print*, '  print_all = ',print_all
      print*, '  dave = ',dave
      print*, '  modV = ',modV
      print*, '  d4 = ',d4
      print*, '*************************************************'

c
c--Initialise dt
      gt = 0.
      if (.true.) then
         if (dumprange) then
            iloop = kstop-1
         else
            iloop = iarray(kstop-1)
         endif
         write(dump,'(a,I4.4)') prefix,iloop
         inquire(file=dump,exist=kexists)
         if (kexists) then
            open (unit=11, file=dump, form='unformatted', recl=maxrec)
            call rdump(11, ichkl, 0)
         endif
      endif
c
c--Process data files
c
      do k=kstop,kstart,-1
         IF (dumprange) THEN
            iloop = k
         ELSE
            iloop = iarray(k)
         ENDIF
         write(dump,'(a,I4.4)') prefix,iloop
         inquire(file=dump,exist=kexists)
         if (kexists) then
            gt_prev = gt
            print*, '**********************'
            print*, 'analysing file ',dump
            print*, '**********************'
            open (unit=11, file=dump, form='unformatted', recl=maxrec)
            call rdump(11, ichkl, 0)
            if (firstcall) then
               udensity  = umassi/udisti**3
               udist     = udisti
               umass     = umassi
               utime     = utimei
               umagfd    = umagfdi
               uvel      = udist/utime
               uergg     = uvel**2
               udistau   = udisti/au
               rthresh   = rthresh_au/udistau
               rthresh2  = rthresh**2
               rsink     = rsink_au/udistau
               rsink2    = rsink**2
               runcond   = runcond_au/udistau
               runcond2  = runcond**2
               rhothresh = rhothresh_cgs/udensity
               rhothresh4(1,:) = rhothresh*10.
               rhothresh4(2,:) = rhothresh*100.
               rhothresh4(3,:) = rhothresh*1000.
               rhodense  = rhodense_cgs/udensity
               iptmax    = 0
               msink     = 0.
               firstcall = .false.
            endif
            mcore = 0.
            dt = abs(gt_prev - gt) ! abs required for the first step since we are predicting dt
            print*, 'The timestep is ',dt

            !--ensure we have all the sinks; this is necessary every
            !  loop since sinks do merge, so a new one will appear part
            !  way through the sim
            isink_status = -isink_status
            iptmaxold = iptmax
            do j = 1,nptmass
               call sinklist_getid(isink,iptdim,iptmax,
     &              iunique(listpm(j)),idsinkarray)
               isink_status(isink) = 1    ! sink exists in this dump
               if (iptmaxold < iptmax) then
                  ! we have a new sink!
                  iptmaxold = iptmax
                  sinkfinal(1,isink) = gt
                  sinkfinal(2,isink) = xyzmh(4,listpm(j))
                  write(fileout,'(a,I2.2,a)') 'core_Masses',isink,'.dat'
                  open(unit=3000+isink,file=trim(fileout))
                  write(3000+isink,"('#',5(1x,'[',i2.2,1x,a11,']',2x))")
     &               1,'idump',
     &               2,'time',
     &               3,'CoreMassNow',  ! i.e. the particles tagged on just this loop
     &               4,'CoreMassAll',  ! i.e. all the particles ever tagged
     &               5,'Sink Mass'
               endif
               msink(2,isink) = msink(1,isink)     ! The previous sink mass
               msink(1,isink) = xyzmh(4,listpm(j)) ! The current sink mass
               print*, j, iunique(listpm(j)),xyzmh(1:3,listpm(j))
            enddo
            !--ensure all sinks from previous dumps are properly tagged
            do j = 1,iptmax
               if (    isink_status(j) == -3 ) then
                  isink_status(j)   = 3   ! the core of this sink has been found
               elseif (isink_status(j) == -2 ) then
                  isink_status(j)   = 2   ! still looking for the core
!               elseif (isink_status(j) == -1 ) then
!                  isink_status(j)   = 2   ! this is the first dump without the core 
!                  sinkformt(iptdim) = gt
!                  sinkformk(iptdim) = iloop
               endif
            enddo
            !
            !
            !--loop through all the particles to see who is close enough to a star
            !
            !
            if (mass_criteria .and. .false.) then
               ! algorithm to enforce our new mass criteria
               do j = 1,nptmass
                  ikount = 0
                  dr   = huge(dr(1))
                  xj   = xyzmh(1,listpm(j))
                  yj   = xyzmh(2,listpm(j))
                  zj   = xyzmh(3,listpm(j))
                  do i = 1,npart
                     if (iphase(i)==0) then
                        ikount       = ikount + 1
                        list(ikount) = i
                        xi   = xyzmh(1,i)
                        yi   = xyzmh(2,i)
                        zi   = xyzmh(3,i)
                        dr(i) = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                     endif
                  enddo
                  CALL indexx(ikount, list, dr, indx)

                  mtot   = msink(1,isink)
                  ii     = 1
                  idsink = iunique(listpm(j))
                  call sinklist_getid(isink,iptdim,iptmax,
     &                        idsink,idsinkarray)
                  keep_searching = .true.
                  do while (keep_searching .and. ii < ikount)
                     i      = list(indx(ii))
                     ii     = ii + 1
                     mtot   = mtot + xyzmh(4,i)
                     iui    = iunique(i)
                     icore(iui,isink) = 2
                     if (mtot > sinkfinal(2,isink)) then
                        keep_searching = .false.
                     endif
                  enddo
               enddo
            endif

            icore0 = 0
            mtmp   = 0.
            print*, 'Finding gas particles near sinks'
            ! find gas particles near sinks
cc$omp  parallel do default(none)
cc$omp& shared(npart,nptmass,iphase,xyzmh,rho,listpm,iunique,iptmax)
cc$omp& shared(vxyzu,Bxyz,ekcle,uergg,udistau,test_Emag,debug)
cc$omp& shared(idsinkarray,rthresh2,check_bound,udistsau,rhodense,modV)
cc$omp& shared(Ecoef,runcond2,most_bound,icore,tcore,eccentricity0,dt)
cc$omp& private(i,j,rad2,xi,xj,yi,yj,zi,zj,mi,mj,dvx,dvy,dvz,Epot,Ekin)
cc$omp& private(vxi,vyi,vzi,iui,vmhdxi,vmhdyi,vmhdzi,Emag)
cc$omp& private(idsink,isink,jsink,iu,is_bound,epot_array)
cc$omp& private(dv2x,dv2y,dv2z)
cc$omp& reduction(+:icore0,mtmp)
            do i = 1,npart
               !print*, 'i=',i
               if (iphase(i)==0) then
                  iui = iunique(i)
                  if (iui > idim) print*, 'PANIC'
                  xi  = xyzmh(1,i)
                  yi  = xyzmh(2,i)
                  zi  = xyzmh(3,i)
                  mi  = xyzmh(4,i)
                  vxi = vxyzu(1,i)
                  vyi = vxyzu(2,i)
                  vzi = vxyzu(3,i)
                  call v_magnetosonic(vmhdxi,vmhdyi,vmhdzi,Bxyz(1:3,i),
     &                                rho(i),Rg,vxyzu(4,i),ekcle(3,i),
     &                                gmw,uergg,modV)
                  Emag=0.5*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))/rho(i)
                  epot_array = 1.0
                  do j = 1,nptmass
                     !print*, 'start',i,j
                     idsink = iunique(listpm(j))
                     call sinklist_getid(isink,iptdim,iptmax,
     &                     idsink,idsinkarray)

                     !--only look at active sinks
                     xj   = xyzmh(1,listpm(j))
                     yj   = xyzmh(2,listpm(j))
                     zj   = xyzmh(3,listpm(j))
                     mj   = xyzmh(4,listpm(j))
                     rad2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                     if (rad2 < rthresh2) then
                        is_bound = .true.
                        if (check_bound .and. rad2 > 0.) then
                           ! relative velocity
                           dvx  = calc_dv(vxi,vxyzu(1,listpm(j)),vmhdxi)  !  = vxi - vxyzu(1,listpm(j)) if vmhdxi=0
                           dvy  = calc_dv(vyi,vxyzu(2,listpm(j)),vmhdyi)
                           dvz  = calc_dv(vzi,vxyzu(3,listpm(j)),vmhdzi)
                           Epot = -(mi+mj)/sqrt(rad2)
                           Ekin = 0.5*(dv2x + dv2y + dv2z)
                           if (test_Emag) Ekin = min(Ekin,Emag)
                           epot_array(j) = Ecoef*Ekin + Epot
                           !--apply density criteria for distant particles
                           if (rad2 > runcond2) then
                              if (rho(i) < rhodense) then
                                 epot_array(j) = abs(epot_array(j))
                              endif
                           endif
                           if (epot_array(j) > 0.) is_bound = .false.
                           if (is_bound) then
                              ! state unbound if high eccentricity
                              is_bound = low_ecc(mi+mj,Ekin,Epot,
     &                                       xi-xj,yi-yj,zi-zj,
     &                                       dvx,dvy,dvz,eccentricity0)
                           endif
                        endif
                        if (icore(iui,isink)==2) is_bound = .true.

                        if (is_bound .and. .not. most_bound) then
                           icore(iui,isink) = 1 ! identify this particle with this sink
                           tcore(iui,isink) = tcore(iui,isink) + dt
                           !mtmp(isink)      = mtmp(isink) + xyzmh(4,i)
                           if (debug) write(76,'(3Es18.6,3I8)')
     &                     xyzmh(1:3,i)*udistau,j,idsink,isink
                        endif
                     endif
                     ! determine the mass of all the active & tagged gas
                     if (icore(iui,isink)==1 .and. .true.) then
                        mtmp(isink) = mtmp(isink) + xyzmh(4,i)
                     endif
                     !print*, 'end  ',i,j
                  enddo
                  if (most_bound) then
                     !--find the most-bound particle
                     jsink = 1
                     do j = 2,nptmass
                        if (epot_array(j) < epot_array(jsink)) jsink=j
                     enddo
                     !--add the most bound particle to the array
                     if (epot_array(jsink) < 0.) then ! since it's possible that i is not associated with any sink
                        idsink = iunique(listpm(jsink))
                        call sinklist_getid(isink,iptdim,iptmax,
     &                                      idsink,idsinkarray)
                        iui = iunique(i)
                        if (iui > idim) print*, 'PANIC'
                        icore0(iui,isink) = icore0(iui,isink) + 1 ! identify this particle with this sink
                        tcore(iui,isink)  = tcore(iui,isink)  + dt
!                        mtmp(isink)       = mtmp(isink) + xyzmh(4,i)
                        if (debug) write(76,'(3Es18.6,3I8)')
     &                  xyzmh(1:3,i)*udistau,j,idsink,isink
                     endif
                  endif
               endif
            enddo
cc$omp end parallel do
            print*, 'Done finding gas particles near sinks'
            icore = icore + icore0
            mcore = mcore + mtmp

            if (mass_criteria .and. .true.) then
               ! it is likely that there is not enough tagged mass.  
               ! We will determine the total amount of active & tagged mass,
               ! and augment it with more if necessary.
               do j = 1,nptmass
                  idsink = iunique(listpm(j))
                  call sinklist_getid(isink,iptdim,iptmax,
     &                  idsink,idsinkarray)

                  if (sinkfinal(3,isink) < 0.) then
                     ! initialise as core +sink mass of final dump
                     sinkfinal(3,isink) = sinkfinal(2,isink)
     &                                  + mcore(isink)
                  endif

                  mtot = mcore(isink) + msink(1,isink)
                  if (mtot < sinkfinal(ifm,isink)) then
                     ikount = 0
                     dr   = huge(dr(1))
                     xj   = xyzmh(1,listpm(j))
                     yj   = xyzmh(2,listpm(j))
                     zj   = xyzmh(3,listpm(j))
                     do i = 1,npart
                        if (iphase(i)==0) then
                           iui = iunique(i)
                           if (icore(iui,isink)==0) then ! untagged
                              ikount       = ikount + 1
                              list(ikount) = i
                              xi    = xyzmh(1,i)
                              yi    = xyzmh(2,i)
                              zi    = xyzmh(3,i)
                              dr(i) = (xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                           endif
                        endif
                     enddo
                     CALL indexx(ikount, list, dr, indx)

                     ii = 1
                     keep_searching = .true.
                     do while (keep_searching .and. ii < ikount)
                        i      = list(indx(ii))
                        ii     = ii + 1
                        mtot   = mtot + xyzmh(4,i)
                        iui    = iunique(i)
                        icore(iui,isink) = 1
                        if (mtot > sinkfinal(ifm,isink)) then
                           keep_searching = .false.
                        endif
                     enddo
                     mcore(isink) = mtot - msink(1,isink)
                     print*, 'tagging ',ii,' particles to increase ',
     &                       'mass around isink=',isink
                  endif
               enddo
            endif
            !
            if (incl_bound_sinks) then
               print*, 'Finding sinks near sinks'
               ! find sinks near sinks
c$omp  parallel do default(none)
c$omp& shared(nptmass,xyzmh,vxyzu,listpm,rthresh2,Ecoef,icoresink,gt)
c$omp& shared(Bxyz,iunique,iptmax,idsinkarray,icore,tcoresink)
c$omp& shared(udistau,iloop,eccentricity0,dt,check_bound,debug)
c$omp& private(i,j,rad2,xi,yi,zi,mi,vxi,vyi,vzi,xj,yj,zj,mj,Ekin,Epot)
c$omp& private(is_bound,dvx,dvy,dvz,idsink,isink,jsink)
               do i = 1,nptmass-1
                  xi  = xyzmh(1,listpm(i))
                  yi  = xyzmh(2,listpm(i))
                  zi  = xyzmh(3,listpm(i))
                  mi  = xyzmh(4,listpm(i))
                  vxi = vxyzu(1,listpm(i))
                  vyi = vxyzu(2,listpm(i))
                  vzi = vxyzu(3,listpm(i))
                  do j = i+1,nptmass
                     xj   = xyzmh(1,listpm(j))
                     yj   = xyzmh(2,listpm(j))
                     zj   = xyzmh(3,listpm(j))
                     mj   = xyzmh(4,listpm(j))
                     rad2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                     if (rad2 < rthresh2) then
                        is_bound = .true.
                        if (check_bound .and. rad2 > 0.) then
                           dvx  = vxi - vxyzu(1,listpm(j))
                           dvy  = vyi - vxyzu(2,listpm(j))
                           dvz  = vzi - vxyzu(3,listpm(j))
                           Epot = -(mi+mj)/sqrt(rad2)
                           Ekin = 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)
                           if (Ecoef*Ekin+Epot > 0.) is_bound = .false.
c                          Ensure low eccentricity
                           if (is_bound) then
                              is_bound = low_ecc(mi+mj,Ekin,Epot,
     &                                       xi-xj,yi-yj,zi-zj,
     &                                       dvx,dvy,dvz,eccentricity0)
                           endif
                        endif
                        if (is_bound) then
                           idsink = iunique(listpm(i))
                           call sinklist_getid(isink,iptdim,iptmax,
     &                                      idsink,idsinkarray)
                           idsink = iunique(listpm(j))
                           call sinklist_getid(jsink,iptdim,iptmax,
     &                                      idsink,idsinkarray)
                           icoresink(isink,jsink) = 1
                           icoresink(jsink,isink) = 1
                           tcoresink(isink,jsink) =
     &                        tcoresink(isink,jsink) + dt
                           tcoresink(jsink,isink) =
     &                        tcoresink(jsink,isink) + dt
                           if (debug) then
                              write(777,'(a,Es18.4,5I18,4Es18.4)') 
     &                        'sink-sink bound: ',gt,iloop,isink,jsink,
     &                        iunique(listpm(i)),iunique(listpm(j)),
     &                        Ecoef*Ekin,Ekin,Epot,sqrt(rad2)*udistau
                              !print*, 'bound',i,j,xyzmh(1:3,listpm(i))
                              !print*, 'bound',i,j,xyzmh(1:3,listpm(j))
                              !print*, 'bound',i,j,Epot,Ekin
                              !print*, 'bound',i,j,isink,jsink
                           endif
                        endif
                     endif
                  enddo
               enddo
c$omp end parallel do
               print*, 'Done finding sinks near sinks'
            endif
            !
            !
            !--On the first dump without the sink, use it's ID to find all neighbours one last time
            !
            !
            print*, iptmax
            do j = 1,iptmax
               if (isink_status(j) == -1 ) then
                  print*,'first no-sink dump for ',j
                  ! set / reset values
                  isink_status(j) = 2
                  sinkformt(j) = gt
                  sinkformk(j) = iloop
                  ! find particle number
                  keep_searching = .true.
                  idsink = idsinkarray(j)
                  i = 1
                  do while (keep_searching)
                     if (iphase(i)==0 .and. iunique(i)==idsink) then
                        keep_searching = .false.
                     else
                        i = i + 1
                     endif
                     if (i > npart) then
                        print*, 'CANNOT FIND SINK!',k,j,idsink
                        keep_searching = .false.
                     endif
                  enddo
                  xj = xyzmh(1,i)
                  yj = xyzmh(2,i)
                  zj = xyzmh(3,i)

                  ! if checking boundedness, then find the properties the gas particles that make up the sink
                  if (check_bound) then
                     ! sort particles by radius from the target particle
                     ikount = 0
                     dr     = huge(dr(1))
                     do i = 1,npart
                        if (iphase(i)==0) then
                           ikount       = ikount + 1
                           list(ikount) = i
                           xi   = xyzmh(1,i)
                           yi   = xyzmh(2,i)
                           zi   = xyzmh(3,i)
                           dr(i) = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                        endif
                     enddo
                     CALL indexx(ikount, list, dr, indx)

                     nneigh = 0
                     mtot   = 0.
                     vx     = 0.
                     vy     = 0.
                     vz     = 0.
                     ii     = 1
                     keep_searching = .true.
                     do while (keep_searching .and. ii < ikount)
                        i      = list(indx(ii))
                        ii     = ii + 1
                        nneigh = nneigh + 1
                        mtot   = mtot + xyzmh(4,i)
                        vx     = vx   + vxyzu(1,i)
                        vy     = vy   + vxyzu(2,i)
                        vz     = vz   + vxyzu(3,i)
                        iui    = iunique(i)
                        ! tag all particles within the given mass or radius as part of the parent clump
                        if (mass_criteria) then
                           if (mtot > msink(1,j)) then
                              keep_searching = .false.
                           endif
                        else
                           if (dr(i)**2 > rsink2) then
                              keep_searching = .false.
                           endif
                        endif
                     enddo

                     if (nneigh > 0) then
                        vx = vx/nneigh
                        vy = vy/nneigh
                        vz = vz/nneigh
                     endif
                  endif
                  write(888,'(a,3I10,5Es18.6)') 'core created: ',
     &            k,j,nneigh,mtot,vx,vy,vz,gt
                  nctrA = 0
                  nctrB = 0
                  mtmp  = 0.
                  ! tag neighbouring particles
c$omp  parallel do default(none)
c$omp& shared(npart,iphase,xyzmh,rho,rhodense,iunique,icore)
c$omp& shared(idsinkarray,rthresh2,j,xj,yj,zj,Ecoef,runcond2)
c$omp& shared(most_bound,vx,vy,vz,mtot,nptmass,listpm,vxyzu)
c$omp& shared(list,indx,nneigh)
c$omp& shared(Bxyz,ekcle,uergg,modV,test_Emag,check_bound,tcore,dt)
c$omp& private(i,jj,jsink,rad2,xi,yi,zi,mi,mj,Ekin,Epot,Epotj,is_bound)
c$omp& private(dvx,dvy,dvz,epot_array,tcoretmp,vmhdxi,vmhdyi,vmhdzi)
c$omp& private(Emag,eccentricity0,iui,ii,in_sink)
c$omp& reduction(+:nctrA,nctrB,mtmp)
                  do i = 1,npart
                     if (iphase(i)==0) then
                        iui  = iunique(i)
                        xi   = xyzmh(1,i)
                        yi   = xyzmh(2,i)
                        zi   = xyzmh(3,i)
                        mi   = xyzmh(4,i)
                        rad2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                        Emag = 0.5*dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
     &                       / rho(i)
                        if (rad2 < rthresh2) then
                           is_bound = .true.
                           if (check_bound) then
                              call v_magnetosonic(vmhdxi,vmhdyi,vmhdzi,
     &                                Bxyz(1:3,i),rho(i),Rg,vxyzu(4,i),
     &                                ekcle(3,i),gmw,uergg,modV)
                              dvx  = calc_dv(vxyzu(1,i),vx,vmhdxi)  !  = vxyzu(1,i)-vx if vmhdxi=0
                              dvy  = calc_dv(vxyzu(2,i),vy,vmhdyi)
                              dvz  = calc_dv(vxyzu(3,i),vz,vmhdzi)
                              Epot = -(mi+mtot)/sqrt(rad2)
                              Ekin = 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)
                              if (test_Emag) Ekin = min(Ekin,Emag)
                              Epotj = Ecoef*Ekin + Epot
                              !--apply density criteria for distant particles
                              if (rad2 > runcond2) then
                                 if (rho(i) < rhodense) then
                                    Epotj = abs(Epotj)
                                 endif
                              endif
                              if (Epotj > 0.) is_bound = .false.
c                             Ensure low eccentricity
                              if (is_bound) then
                                 is_bound = low_ecc(mi+mtot,Ekin,Epot,
     &                                          xi-xj,yi-yj,zi-zj,
     &                                       dvx,dvy,dvz,eccentricity0)
                              endif
                           endif

                           in_sink = .false.
                           ! see if this particle should be bound for mass/radius criteria
                           do ii = 1,nneigh+1
                              if (i==list(indx(ii))) then
                                 is_bound = .true.
                                 in_sink  = .true.
                              endif
                           enddo

                           if (is_bound .and. .not. most_bound) then
                              tcoretmp = abs(tcore(iui,j))+dt
                              !mtmp(j) = mtmp(j) + xyzmh(4,i)
                              if (in_sink) then
                                 icore(iui,j) = -1 ! particle will likely merge with sink
                                 tcore(iui,j) = -tcoretmp
                                 nctrA = nctrA + 1
                              else
                                 icore(iui,j) =  1
                                 if (tcore(iui,j) < 0.) then
                                    tcore(iui,j) = -tcoretmp
                                 else
                                    tcore(iui,j) =  tcoretmp
                                 endif
                                 nctrB = nctrB + 1
                              endif
                           endif
                        else
                           is_bound = .false.
                        endif
                        if (abs(icore(iui,j))==1) then
                           mtmp(j) = mtmp(j) + xyzmh(4,i)
                        endif

                        ! ensure that i is most tightly bound to this soon to be sink
                        if (most_bound .and. is_bound) then ! (note that check_bound==.true.)
                           epot_array = 1.
                           do jj = 1,nptmass
                              xj  = xyzmh(1,listpm(jj))
                              yj  = xyzmh(2,listpm(jj))
                              zj  = xyzmh(3,listpm(jj))
                              mj  = xyzmh(4,listpm(jj))
                              rad2= (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                              if (rad2 < rthresh2) then
                                 dvx  = vxyzu(1,i)-vxyzu(1,listpm(jj))
                                 dvy  = vxyzu(2,i)-vxyzu(2,listpm(jj))
                                 dvz  = vxyzu(3,i)-vxyzu(3,listpm(jj))
                                 Epot = -(mi+mj)/sqrt(rad2)
                                 Ekin = 0.5*(dvx*dvx+dvy*dvy+dvz*dvz)
                                 epot_array(jj) = Ecoef*Ekin + Epot
                              endif
                           enddo
                           !--ensure that we are most bound to the soon to be formed sink
                           do jj = 1,nptmass
                              if (epot_array(jj)<Epotj) is_bound=.false.
                           enddo
                           if (is_bound) then
                              icore(iunique(i),j) = -1 ! unconditionally identify this particle with this sink
                              mtmp(j) = mtmp(j) + xyzmh(4,i)
                           endif
                        endif
                     endif
                  enddo
c$omp end parallel do
                  mcore = mcore + mtmp
                  write(888,'(a,4I10,1Es18.6)') 'core augmented: ',
     &            iloop,j,nctrA,nctrB,gt

                  if (mass_criteria .and. .true.) then
                     ii = nneigh+2 ! since we have previously forced ii < nneigh+1 to be tagged
                     do while (mcore(j)<sinkfinal(ifm,j).and.ii<ikount)
                        ! We are still missing mass, so start tagging particle outward from the candidate
                        i   = list(indx(ii))
                        ii  = ii + 1
                        iui = iunique(i)
                        if (icore(iui,j)==0) then ! add this particle to the count to ensure we have enough mass
                           mcore(j)     = mcore(j) + xyzmh(4,i)
                           icore(iui,j) = -1
                           tcore(iui,j) = -tcoretmp
                           nctrA        = nctrA + 1
                        endif
                     enddo
                  endif
                  write(888,'(a,4I10,1Es18.6)') 'core re-augmented: ',
     &            iloop,j,nctrA,nctrB,gt
 
                  ! check to see if there are any sinks nearby the group about to from
                  if (incl_bound_sinks) then ! since this will not be symmetric due to when sinks form
c$omp  parallel do default(none)
c$omp& shared(j,nptmass,xyzmh,vxyzu,listpm,rthresh2,Ecoef,icoresink)
c$omp& shared(xj,yj,zj,mtot,vx,vy,vz,iloop,check_bound,iunique,gt)
c$omp& shared(idsinkarray,iptmax)
c$omp& private(i,rad2,xi,yi,zi,mi,Ekin,Epot,is_bound,dvx,dvy,dvz)
c$omp& private(idsink,isink)
                     do i = 1,nptmass
                        xi   = xyzmh(1,listpm(i))
                        yi   = xyzmh(2,listpm(i))
                        zi   = xyzmh(3,listpm(i))
                        mi   = xyzmh(4,listpm(i))
                        rad2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                        if (rad2 < rthresh2) then
                           is_bound = .true.
                           if (check_bound) then
                              dvx  = vx - vxyzu(1,listpm(i))
                              dvy  = vy - vxyzu(2,listpm(i))
                              dvz  = vz - vxyzu(3,listpm(i))
                              Epot = -(mtot+mi)/sqrt(rad2)
                              Ekin = 0.5*(dvx*dvx + dvy*dvy + dvz*dvz)
                              if (Ecoef*Ekin+Epot > 0.) then
                                 is_bound = .false.
                              endif
                           endif
                           if (is_bound) then
                              idsink = iunique(listpm(i))
                              call sinklist_getid(isink,iptdim,iptmax,
     &                                         idsink,idsinkarray)
                              icoresink(isink,j) = 1
                              icoresink(j,isink) = 1
                              write(777,'(a,Es18.4,3I18)') 
     &                        'sink-core bound: ',gt,iloop,isink,j
                           endif
                        endif
                     enddo
c$omp end parallel do
                  endif
               endif
            enddo
            !
            !--Record sink & core properties; note that the core cannot gain mass in the next loop
            !
            print*,'printing bookkeeping for masses'
            do j = 1,iptmax
               mtot = 0.
c$omp  parallel do default(none)
c$omp& shared(j,npart,icore,xyzmh)
c$omp& private(i)
c$omp& reduction(+:mtot)
               do i = 1,npart
                  if (abs(icore(i,j))==1) then
                     mtot = mtot + xyzmh(4,i)
                  endif
               enddo
c$omp end parallel do
               write(3000+j,'((I18,1x),4(es18.4,1x))') iloop,gt,
     &              mcore(j),mtot,msink(1,j)
            enddo
            print*,'done printing bookkeeping for masses'
            !
            !
            !--If we are now before the star has formed, loop through the particles to determine if we have the initial core 
            !
            !
            do j = 1,iptmax
               if (debug) print*,'looping ',j,' of ',iptmax
               if (isink_status(j)==2 .or. print_evol) then
                  if (debug) print*,'searching for core for ',j
                  rhomax = 0.
                  rhoave = 0.
                  nneigh = 0
                  sinkage = sinkfinal(1,j)-sinkformt(j)
c$omp  parallel do default(none)
c$omp& shared(j,npart,iphase,iunique,icore,rho,xyzmh,most_bound)
c$omp& private(i,is_associ)
c$omp& reduction(max:rhomax)
c$omp& reduction(+:rhoave,nneigh)
                  do i = 1,npart
                     is_associ = is_assoc(j,icore(iunique(i),:),
     &                                    iptdim,most_bound)
                     if (iphase(i)==0 .and. is_associ) then
                        rhomax = max(rhomax,rho(i))
                        rhoave = rhoave + rho(i)
                        nneigh = nneigh + 1
                     endif
                  enddo
c$omp end parallel do
                  if (nneigh > 0) rhoave = rhoave/nneigh
                  if (debug) then
                     print*,'core ',j,' values: ',rhomax*udensity
     &               ,nneigh,mtot
                  endif

                  ! Determing if we are printing now
                  print_now = .false.
                  is_final  = .false.
                  if (isink_status(j)==2) then
                     if (d4) then
                        if (rhomax < rhothresh4(3,j)) then
                           is_final   = .true.
                           print_now  = .true.
                           write(fileout,'(a,I2.2,a)')
     &                     'core_vals',j,'_dmax14.dat'
                           write(fileoutAll,'(a,I2.2,a)')
     &                     'core_Allvals',j,'_dmax14.dat'
                           rhothresh4(3,j) = -rhothresh4(3,j)
                        elseif (rhomax < rhothresh4(2,j)) then
                           is_final   = .true.
                           print_now  = .true.
                           write(fileout,'(a,I2.2,a)')
     &                     'core_vals',j,'_dmax15.dat'
                           write(fileoutAll,'(a,I2.2,a)')
     &                     'core_Allvals',j,'_dmax15.dat'
                           rhothresh4(2,j) = -rhothresh4(2,j)
                        elseif (rhomax < rhothresh4(1,j)) then
                           is_final   = .true.
                           print_now  = .true.
                           write(fileout,'(a,I2.2,a)')
     &                     'core_vals',j,'_dmax16.dat'
                           write(fileoutAll,'(a,I2.2,a)')
     &                     'core_Allvals',j,'_dmax16.dat'
                           rhothresh4(1,j) = -rhothresh4(1,j)
                        endif
                     endif
                     if (iouti(1,j)==0 .and. rhomax<rhothresh) then
                        is_final   = .true.
                        print_now  = .true.
                        iout       = iout + 1
                        iouti(1,j) = 1
                        if (d4) then
                           write(fileout,'(a,I2.2,a)')
     &                     'core_vals',j,'_dmax17.dat'
                           write(fileoutAll,'(a,I2.2,a)')
     &                     'core_Allvals',j,'_dmax17.dat'
                        else
                           write(fileout,'(a,I2.2,a)')
     &                     'core_vals',j,'_dmax.dat'
                           write(fileoutAll,'(a,I2.2,a)')
     &                     'core_Allvals',j,'_dmax.dat'
                        endif
                     elseif (iouti(2,j)==0 .and. rhoave<rhothresh) then
                        is_final   = .true.
                        print_now  = .true.
                        iout       = iout + 1
                        iouti(2,j) = 1
                        write(fileout,'(a,I2.2,a)')
     &                  'core_vals',j,'_dave.dat'
                        write(fileoutAll,'(a,I2.2,a)')
     &                  'core_Allvals',j,'_dave.dat'
                        if (.not.dave) print_now = .false.
                     endif
                  endif
                  if (print_evol) then
                     print_now = .true.
                     ! overwrite file names if ploting evolution
                     write(fileout,'(a,I2.2,a,I4.4,a)')
     &               'core_vals',j,'_',iloop,'.dat'
                     write(fileoutAll,'(a,I2.2,a,I4.4,a)')
     &               'core_Allvals',j,'_',iloop,'.dat'
                  endif
                  if (print_evol .and. isink_status(j)==3) then
                     ! stop printing after the core has formed
                     print_now = .false.
                  endif
                  if (iouti(1,j)==1 .and. iouti(2,j)==1) then
                     isink_status(j) = 3
                  endif

                  ! Print the particle data
                  if (print_now) then
                     ! loop through all the particles to determine the density weighted centre of mass
                     nneigh = 0
                     xcom   = 0.
                     ycom   = 0.
                     zcom   = 0.
                     mtot   = 0.
                     rtot   = 0.
c$omp  parallel do default(none)
c$omp& shared(j,npart,iphase,xyzmh,rho,iunique,icore,most_bound)
c$omp& private(i,is_associ)
c$omp& reduction(+:nneigh,mtot,xcom,ycom,zcom,rtot)
                     do i = 1,npart
                        is_associ = is_assoc(j,icore(iunique(i),:),
     &                                       iptdim,most_bound)
                        if(iphase(i)==0 .and. is_associ) then
                           nneigh = nneigh + 1
                           xcom   = xcom + rho(i)*xyzmh(1,i)
                           ycom   = ycom + rho(i)*xyzmh(2,i)
                           zcom   = zcom + rho(i)*xyzmh(3,i)
                           mtot   = mtot + xyzmh(4,i)
                           rtot   = rtot + rho(i)
                        endif
                     enddo
c$omp end parallel do
                     if (nneigh > 0) then
                        xcom = xcom/rtot
                        ycom = ycom/rtot
                        zcom = zcom/rtot
                     endif
                     ! loop through all the particles to determine core radius
                     radius2 = 0.
                     itag    = 0
                     ttag    = 0.
c$omp  parallel do default(none)
c$omp& shared(j,npart,iphase,xyzmh,xcom,ycom,zcom,itag,ttag)
c$omp& shared(sinkage,icore,iunique,most_bound,tcore)
c$omp& private(i,rtmp2,is_associ)
c$omp& reduction(max:radius2)
                     do i = 1,npart
                        is_associ = is_assoc(j,icore(iunique(i),:),
     &                                       iptdim,most_bound)
                        if(iphase(i)==0 .and. is_associ) then
                           rtmp2 = (xyzmh(1,i) - xcom)**2
     &                           + (xyzmh(2,i) - ycom)**2
     &                           + (xyzmh(3,i) - zcom)**2
     &                           +  xyzmh(5,i)**2  ! to endure that we can tag all the neighbours
                           radius2 = max(radius2,rtmp2)
                           itag(i) = icore(iunique(i),j)
                           ttag(i) = tcore(iunique(i),j)/sinkage
                        endif
                     enddo
c$omp end parallel do
                     ! loop through all the particles to determine total values
                     mtot_all   = 0.
                     nneigh_all = 0
c$omp  parallel do default(none)
c$omp& shared(npart,iphase,xyzmh,xcom,ycom,zcom,radius2)
c$omp& private(i,rtmp2)
c$omp& reduction(+:mtot_all,nneigh_all)
                     do i = 1,npart
                        if (iphase(i)==0) then
                           rtmp2 = (xyzmh(1,i) - xcom)**2
     &                           + (xyzmh(2,i) - ycom)**2
     &                           + (xyzmh(3,i) - zcom)**2
                           if (rtmp2 <= radius2) then
                              mtot_all   = mtot_all   + xyzmh(4,i)
                              nneigh_all = nneigh_all + 1
                           endif
                        endif
                     enddo
c$omp end parallel do
                     ! loop through all the particles to add particles associated with bound sinks
                     if (incl_bound_sinks) then
                        radius2a = 0.
c$omp  parallel do default(none)
c$omp& shared(j,npart,iphase,xyzmh,xcom,ycom,zcom,itag,ttag,tcoresink)
c$omp& shared(icoresink,sinkage,iptmax,icore,iunique,most_bound)
c$omp& private(i,jj,rtmp2,is_associ)
c$omp& reduction(max:radius2a)
                        do i = 1,npart
                           if (iphase(i)==0 .and. itag(i)==0) then  ! active particle that is not associated with this sink
                              do jj = 1,iptmax
                                 if (icoresink(j,jj) == 1) then
                                    is_associ = is_assoc(jj,
     &                                          icore(iunique(i),:),
     &                                          iptdim,most_bound)
                                    if (is_associ) then ! is associated with sink jj but not j; add to list
                                       rtmp2 = (xyzmh(1,i) - xcom)**2
     &                                       + (xyzmh(2,i) - ycom)**2
     &                                       + (xyzmh(3,i) - zcom)**2
     &                                       +  xyzmh(5,i)**2  ! to endure that we can tag all the neighbours
                                       radius2a = max(radius2a,rtmp2)
                                       itag(i) = jj+100 ! this might be a race condition if associated with two other sinks, but the important thing is it becomes non-zero
                                       ttag(i) =tcoresink(j,jj)/sinkage
                                    endif
                                 endif
                              enddo
                           endif
                        enddo
c$omp end parallel do
                        radius2 = max(radius2a,radius2)
                     endif

                     ! open files and write header
                     open(unit=223,file=trim(fileout))
                     if ( print_all .or. 
     &                   (print_evolAF .and. is_final)) then
                        open(unit=224,file=trim(fileoutAll))
                        pend = 224
                     else
                        pend = pstart
                     endif
                     do p = pstart,pend
                        write(p,'(a,g18.6)')
     &                  'the simulation ends at t = ',sinkfinal(1,j)
                        write(p,'(a,g18.6)')
     &                  'the sinks final age is dt = ',sinkage
                        write(p,'(a,g18.6)') 
     &                  'the sinks final mass is m = ',sinkfinal(2,j)
                        write(p,'(a,g18.6,a,I5)')
     &                  'the sink forms shortly after t = ',
     &                  sinkformt(j),' and dump ',sinkformk(j)
                        write(p,'(a,g18.6,a,I5)')
     &                  'the core forms at            t = ',gt,
     &                  ' and dump ',iloop
                        write(p,'(a,es18.10)')
     &                  'The radius of the core is ',
     &                  sqrt(radius2)*udistau
                        write(p,'(a,I8,a,es18.10)')
     &                  'Using only tagged particles, there are ',
     &                  nneigh,
     &                  ' particles in the core, for a total mass of ',
     &                  mtot
                        write(p,'(a,I8,a,es18.10)')
     &                  'Using all particles, there are ',
     &                  nneigh_all,
     &                  ' particles in the core, for a total mass of ',
     &                  mtot_all
                        write(p,'(a,3E18.10)')
     &                 'the centre of mass is ',
     &                  xcom*udistau,ycom*udistau,zcom*udistau
                        write(p,'(a)') ' '
                        write(p,"('#',15(1x,'[',i2.2,1x,a11,']',2x))")
     &                     1,'x [au]',
     &                     2,'y [au]',
     &                     3,'z [au]',
     &                     4,'rho [g/cm^3]',
     &                     5,'vx [cm/s]',
     &                     6,'vy [cm/s]',
     &                     7,'vz [cm/s]',
     &                     8,'u [erg]',
     &                     9,'Bx [G]',
     &                    10,'By [G]',
     &                    11,'Bz [G]',
     &                    12,'T [K]',
     &                    13,'t(bnd)/age',
     &                    14,'tag',
     &                    15,'iunique'
                     enddo
                     ! print the location of all the associated sinks and their precursors
                     jj = j
                     j0 = .true.
                     do while (jj <= iptmax)
                        if(j0 .or. (jj/=j .and. icoresink(j,jj)==1))then
                           call get_sink(jj,idsinkarray,xyzmhj,rhoj,
     &                                   vxyzj,uj,Bxyzj,Tj,idj,iuj)
                           do p = pstart,pend
                              write(p,'(13(es18.4,1x),2(I18,1x))')
     &                            xyzmhj*udistau,rhoj*udensity,
     &                            vxyzj*uvel,uj*uergg,Bxyzj*umagfd,
     &                            Tj,1.0,idj,iuj
                           enddo
                        endif
                        if (j0) then
                           j0 = .false.
                           jj = 1
                        else
                           jj = jj + 1
                        endif
                     enddo

                     ! loop through all particle again and print to file
                     do i = 1,npart
                        if (iphase(i)==0) then
                           rtmp2 = (xyzmh(1,i) - xcom)**2
     &                           + (xyzmh(2,i) - ycom)**2
     &                           + (xyzmh(3,i) - zcom)**2
                           if (rtmp2 <= radius2) then
                              is_associ =is_assoc(j,icore(iunique(i),:),
     &                                            iptdim,most_bound)
                              if(incl_bound_sinks .and. itag(i)/=0)then
                                 is_associ = .true.
                              endif
                              do p = pstart,pend
                                 if (p==224 .or.
     &                              (p==223 .and. is_associ)) then
                                    write(p,'(13(es18.4,1x),2(I18,1x))')
     &                               xyzmh(1:3,i)*udistau,
     &                               rho(i)*udensity,
     &                               vxyzu(1:3,i)*uvel,
     &                               vxyzu(4,i)*uergg,
     &                               Bxyz(1:3,i)*umagfd,
     &                               vxyzu(4,i)/ekcle(3,i),
     &                               ttag(i),itag(i),iunique(i)
                                 endif
                              enddo
                           endif
                        endif
                     enddo
                     do p = pstart,pend
                        close(p)
                     enddo
                  endif
               endif
            enddo
            if (iout==2*iptmax) then
               print*, 'We have obtained outputs for all sinks'
               print*, 'Task completed'
               do i = 1,iptmax
                  close(3000+i)
               enddo
               stop
            endif
            do i = 1,iptmax
               do j = 1,iptmax
                  if (icoresink(i,j)/=icoresink(j,i)) then
                     print*, "BUG IN DUMP k"
                     print*, 'error with sinks ',i,j
                     do jj = 1,iptmax
                        print*,icoresink(jj,1:iptmas)
                     enddo
                     stop
                  endif
               enddo
            enddo

         endif
      enddo

      print*, 'Task completed, but likely not with all sinks'
      STOP
      END
c
c-----------------------------------------------------------------------
c     get dv, accounting for magnetosonic speed
      function calc_dv(vi,vj,vmhdi)
      real vi,vj,vmhdi,calc_dv
      real dvp,dvn

      dvn = (vi - vmhdi) - vj
      dvp = (vi + vmhdi) - vj
      if (abs(dvn) < abs(dvp)) then
         calc_dv = dvn
      else
         calc_dv = dvp
      endif

      end function calc_dv

c-----------------------------------------------------------------------
c     calculate the fast magnetosonic wave
      subroutine v_magnetosonic(vmhdxi,vmhdyi,vmhdzi,Bxyzi,rhoi,Rg,ui,
     &                          ekcle3,gmw,uergg,modV)
      real vmhdxi,vmhdyi,vmhdzi,Bxyzi(3),rhoi,Rg,ui,ekcle3,gmw,uergg
      logical modV

      if (modV) then
         va2x = Bxyzi(1)*Bxyzi(1)/rhoi
         va2y = Bxyzi(2)*Bxyzi(2)/rhoi
         va2z = Bxyzi(3)*Bxyzi(3)/rhoi
         va2  = va2x + va2y + va2z
         cs2  = 5.*Rg*ui/(3.*ekcle3*gmw*uergg)
         vacs = va2 + cs2
         vmhdxi = sqrt(0.5*(vacs + sqrt(vacs**2 - 4.*va2x*cs2)))
         vmhdyi = sqrt(0.5*(vacs + sqrt(vacs**2 - 4.*va2y*cs2)))
         vmhdzi = sqrt(0.5*(vacs + sqrt(vacs**2 - 4.*va2z*cs2)))
      else
         vmhdxi = 0.0
         vmhdyi = 0.0
         vmhdzi = 0.0
      endif

      end subroutine v_magnetosonic
c-----------------------------------------------------------------------
c     adjust boundness based upon eccentricity
      function low_ecc(mtot,Ekin,Epot,dx,dy,dz,dvx,dvy,dvz,eccen0)
      real mtot,Ekin,Epot,dx,dy,dz,dvx,dvy,dvz,eccen0
      real dangx,dangy,dangz,dang2,semimajor,val,eccentricity
      logical low_ecc

      dangx = dvy*dz - dvz*dy
      dangy = dvz*dx - dvx*dz
      dangz = dvx*dy - dvy*dx
      dang2 = (dangx**2 + dangy**2 + dangz**2)
      semimajor = -mtot/(2.0*(Ekin+Epot))
      val = 1.0 - dang2/mtot/semimajor
      if (val >= 0.0) then
         eccentricity = SQRT(val)
      else
         eccentricity = 1.0 + eccen0
      endif
      if (eccentricity < eccen0) then
         low_ecc = .true.
      else
         low_ecc = .false.
      endif

      end function low_ecc

c-----------------------------------------------------------------------
c     determine if the particle is associated with the star
      function is_assoc(j,icore,iptdim,most_bound)
      integer j,jj,icore(iptdim)
      logical is_assoc,most_bound

      if (most_bound) then
         if (icore(j) == 0) then
            is_assoc = .false. ! not associated
         elseif (icore(j) == -1) then
            is_assoc = .true. ! unconditional association
         else
            ! is associated only if associated more with this sink than any others
            is_assoc = .true.
            do jj = 1,iptdim
               if (icore(jj) > icore(j) .or. icore(jj)==-1) then
                  is_assoc = .false.
               endif
            enddo
         endif
      else
         is_assoc = (abs(icore(j))==1)
      endif

      end function is_assoc
c
c-----------------------------------------------------------------------
c     Get the sink properties
      subroutine get_sink(j,idsinkarray,xyzmhj,rhoj,
     &                    vxyzj,uj,Bxyzj,Tj,idj,iuj)
      INCLUDE 'idim'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/sort'
      integer   :: iii,idj,idsinkarray(iptdim)
      integer*8 :: iuj
      real      :: xyzmhj(3),vxyzj(3),uj,Bxyzj(3),Tj
      real*4    :: rhoj
      logical   :: keep_searching

      ! set defaults
      xyzmhj = 0.
      rhoj   = 0.
      vxyzj  = 0.
      uj     = 0.
      Bxyzj  = 0.
      Tj     = 0.
      idj    = 0
      iuj    = 0

      ! see if a sink exists
      keep_searching = .true.
      iii = 1
      do while (keep_searching .and. iii <= nptmass)
         if (iunique(listpm(iii))==idsinkarray(j))then
            keep_searching = .false.
         else
            iii = iii + 1
         endif
      enddo

      if (keep_searching) then
         ! sink does not exist; search for gas particle
         iii = 1
         do while (keep_searching .and. iii <=npart)
            if (iphase(iii)==0 .and. iunique(iii)==idsinkarray(j)) then
               keep_searching = .false.
            else
               iii = iii + 1
            endif
         enddo
         if (keep_searching) then
            print*,'This is an error that should not occur'
            idj = 0
         else
            ! gas particle precursor exists
            xyzmhj = xyzmh(1:3,iii)
            rhoj   = rho(iii)
            vxyzj  = vxyzu(1:3,iii)
            uj     = vxyzu(  4,iii)
            Bxyzj  = Bxyz( 1:3,iii)
            Tj     = vxyzu(4,iii)/ekcle(3,iii)
            idj    = j
            iuj    = -iunique(iii)
         endif
      else
         ! sink exists
         iii    = listpm(iii)
         xyzmhj = xyzmh(1:3,iii)
         vxyzj  = vxyzu(1:3,iii)
         idj    = -j
         iuj    = -iunique(iii)
      endif

      end subroutine get_sink
