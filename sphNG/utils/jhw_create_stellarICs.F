      PROGRAM jhw_create_stellarICs
c***********************************************************
c  This is the second of two programmes required to        *
c  generate initial conditions for                         *
c  Wurster & Rowan (2023a,b); details given in WR2023a.    *
c  This programme create the 'agumented clump' with the    *
c  background particles.                                   *
c  This programme will read the input from                 *
c  jhw_calc_stellarICs to create a new dump & ifile.       *
c  Developmental options are included below the variable   *
c  declarations, but I would suggest retaining the default *
c  options.                                                *
c  in sphNG/, compile as make jhw_create_stellarICs        *
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
      INCLUDE 'COMMONS/ener1'
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
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/interstellar'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
      INCLUDE 'COMMONS/vargam'
#endif

      REAL*8 umassi,udisti,utimei,umagfdi,uerggi
      COMMON /unitsin/ umassi,udisti,utimei,umagfdi,uerggi
      COMMON /dtmaxin/ dtmaxdp

      INTEGER, parameter :: ncore_max = 2000000
      INTEGER      :: i,p,io,idummy(2),itag(2,idim),ntag,ntag0,pmax
      INTEGER      :: ctr0,ctr1,ncore,lctr,np,iloop,ploop
      INTEGER      :: nctr,nctrtot
      INTEGER      :: ii,jj,kk,idx,ibdy(idim),list(idim),indx(idim)
      REAL         :: hbkg,pmass,rhomin,rhobkg,ffac,mach,have
      REAL         :: sinkage,sinkmass,coreage,coreform,sinkform
      REAL         :: hdbx,hdby,hdbz,timeend,wkern,rho18,rad,rtmp2
      REAL         :: r2xy,r2xz,r2yz,rkern
      REAL         :: xyzcom(3),vxyzcom(3),Bxyzcom(3),Bi(ncore_max)
      REAL         :: nimhdi(4),dr(3),xyz(3,idim),rmin(idim)
      REAL         :: rdummy(13),rhorange(7),Brange(4,7),Ba(3)
      REAL         :: veldisp(4),vave(4)
      REAL         :: Temp(ncore_max),Trange(7),Ta
      REAL         :: Bztheta(ncore_max),Bztrange(7),Bzta
      REAL         :: Lcore(3),Bcore(3),xcore(3),vcore(3),hcore2,BLta
      REAL         :: Lxyz(3),Lcloud(3),BLtheta(ncore_max),BLtrange(7)
      REAL         :: costheta,xyztmp(5,ncore_max),ell(3,2),phi
      REAL         :: cloudxyz(3),rad2_tmp(ncore_max),Erot(3)
      REAL*8       :: uvel,uerg,udensity,uang,udistpc
      CHARACTER*3  :: prefix,outfix
      CHARACTER*7  :: filein,fileout,dat
      CHARACTER*10 :: wfmt
      CHARACTER*11 :: debugfile
      CHARACTER*256:: filetag,filestats,cdummy,fileiunique
      LOGICAL      :: verbose,nooverlap,addparticle
      LOGICAL      :: get_stats,Bbkgave,ubkgave,mod_h
      LOGICAL      :: fexists,set_Beq0,Bbkgbdy
      LOGICAL      :: d4,find_neigh,keepsearching
      LOGICAL      :: casually_connected,old_rad,kp_orig_bkg
      LOGICAL      :: readline,write_bkg,updateB,randomborder
      LOGICAL      :: indiv_out,write_header,find_axes,recentre
      LOGICAL      :: print_iunique

c-----------------------------------------------------------------------
c--   Initialise variables
c-----------------------------------------------------------------------
      iprint = istdout 
      iread  = 5
      maxrec = 100*idim
      varmhd = 'Brho'
      frho   =  30.  ! The ratio between the cloud and the background densities & energies
      do i = 1, idim ! Not sure why, but this loop is necessary
         isort(i) = i
         iorig(i) = i
      enddo
      wfmt         = '(a,3F10.4)'
      verbose      = .false.
      nooverlap    = .true. ! if true, remove particels within h of an active particle
      get_stats    = .true.  ! if true, calculate the properties of the tagged particles
      Bbkgave      = .true. ! if true, use an average B-background; else copy properties from neighbours (accounting for perodic)
      ubkgave      = .true. ! if true, use an average u-background; else copy properties from neighbours (accounting for perodic)
      set_Beq0     = .false. ! will reset B = 0 to turn an MHD dump into an HD dump
      Bbkgbdy      = .true.  ! if true, will set B to average of its neighbours for boundary particles; average of all boundaries for remaining bkg particles
      nextra       =  4
      find_neigh   = .true.  ! will tag particles with nextra*h of a tagged particle
      d4           = .true.  ! set to match stellarICs.F90
      mod_h        = .true.  ! reduce h of the background particles close to tagged particles to get a better initial h guess
      casually_connected = .true. ! only keep particles that are casually connected to particles within 250au of the future sink !with rho > 1e-18g/cm^3
      old_rad      = .true.  ! if true, use flux-limited diffusion; else Bate & Keto RT
      write_bkg    = .false.  ! write the background particles to file for gnuplotting
      randomborder = .false.  ! particles between nextra/2*h & nextra*h will be kept with a probability of 2/frho
      recentre     = .true.  ! recentre particles & zero time
      print_iunique= .true.  ! print a list of iunique numbers to file

      if (set_Beq0 .or. Bbkgbdy) Bbkgave = .true.

      Bbkgbdy = .true.
      Bbkgave = .true.
      indiv_out = .false. ! if false, will append the core stats to an existing file

c-----------------------------------------------------------------------
c--   prompt for input
c-----------------------------------------------------------------------
      write (*,*) 'Input prefix'
      read (iread,*) prefix
      if (prefix=='hyd') Bbkgave = .true.

      write (*,*) 'Input outfix'
      read (iread,*) outfix

      write (*,*) 'Input dump number'
      read (iread,*) idump

      write (*,*) 'Input sink number'
      read (iread,*) isink

      if (d4) then
         write (*,*) 'Input desnsity id'
         read (iread,*) iden
         write(dat,'(I2.2,a)') iden,'.dat'
      else
         write(dat,'(a)')'.dat'
      endif

c     get filenames
      write(fileout,'(a,I4.4)') outfix,0
      write(filetag,'(a,I2.2,2a)') 'core_vals',isink,'_dmax',trim(dat)
      inquire(file=filetag,exist=fexists)
      if (.not.fexists) then
         write(filetag,'(a,I2.2,a,I4.4,a)')
     &       'core_valsA',isink,'_',idump,'.dat'
         inquire(file=filetag,exist=fexists)
      endif
      if (.not.fexists) then
         write(filetag,'(a,I2.2,a,I4.4,a)')
     &       'core_vals',isink,'_',idump,'.dat'
         inquire(file=filetag,exist=fexists)
      endif
      if (.not.fexists) then
         write(filetag,'(a,I2.2,2a)')
     &      'Core__not_most_bound/core_vals',isink,'_dmax',trim(dat)
         inquire(file=filetag,exist=fexists)
      endif
      if (.not.fexists) then
         write(filetag,'(a,I2.2,a,I4.4,a)')
     &       'Core__not_most_bound/core_vals',isink,'_',idump,'.dat'
         inquire(file=filetag,exist=fexists)
      endif
      if (fexists) then
         print*, 'Reading file ',trim(filetag)
      else
         print*, 'cannot find file.  aborting.'
         stop
      endif

c-----------------------------------------------------------------------
c--   Read files
c-----------------------------------------------------------------------
c     from the tagged file, get the iunique numbers
      print*, 'Reading tagged particles'
      open(12,file=filetag)
      do i = 1,11
         read(12,'(a)') cdummy       ! read and ignore header
         if (i==1) read(cdummy(30:  ),*) timeend
         if (i==2) read(cdummy(30:  ),*) sinkage
         if (i==3) read(cdummy(30:  ),*) sinkmass
         if (i==4) read(cdummy(35:50),*) sinkform
         if (i==5) read(cdummy(35:50),*) coreform
         if (i==5) read(cdummy(62:  ),*) idump
         coreage = timeend - coreform
      enddo
      print*, 'times', sinkage,sinkmass,coreage

      write(filein, '(a,I4.4)') prefix,idump
      print*, 'simulation end time is ',timeend
      print*, 'dumpfile to read: ',trim(filein)
      io    = 0
      ntag  = 0
      nsink = 0
      do while(io==0)
         read(12,*,iostat=io) rdummy,idummy ! read the tagged particles
         if (io==0 .and. rdummy(4) > 0.) then ! the latter conditional is to ignore distant particles
            ntag = ntag + 1
            itag(1,ntag) = abs(idummy(2))
            xyz(:,ntag)  = rdummy(1:3)
            if (idummy(2) < 0) nsink = nsink + 1
         endif
      enddo
      close(12)
      ntag0 = ntag
      print*, 'there are ',ntag,' tagged particles, including '
     &      ,nsink,' sinks'

c     Open the original dump file
      write(inname,'(2a)')'ifile_',prefix
      call options
      print*, 'Reading file ',trim(filein)
c     open (11,file=trim(filein),form='unformatted', recl=maxrec)
      open (11,file=trim(filein),form='unformatted')
c     call rdump(11,ichkl,1)
      call rdump(11,ichkl,0)
      close(11)
      if (ichkl.EQ.1) then
         write (*,*) 'ERROR reading input file: ichkl ',ichkl
         stop
      endif
      dtmax    = dtmaxdp
      umass    = umassi
      udist    = udisti
      udistpc  = udist/pc
      utime    = utimei
      umagfd   = umagfdi
      uvel     = udist/utime
      uang     = umass*udist**2/utime
      uergg    = (udist/utime)**2
      uerg     = umass*(udist/utime)**2
      uergcc   = umass/(udist*utime**2)
      udensity = umass/udist**3
      rho18    = 1.0d-18/udensity
      rad250   = 250.*au/udist
      xyz      = xyz*au/udist

c-----------------------------------------------------------------------
c--   Find all relevant particles
c-----------------------------------------------------------------------
c     Loop through all the particles, and kill any that are not in the tagged list
      print*, 'Condensing list'
      npartnew = 0
      nbdy     = 0
      ctr0     = 0
      ctr1     = 0
      have     = 0.
c$omp  parallel do default(none)
c$omp& shared(npart,ntag,xyzmh,vxyzu,iphase,iunique,itag,n1)
c$omp& shared(verbose,timeend)
c$omp& private(i,j,xi,yi,zi)
c$omp& reduction(+:npartnew,nbdy,have)
      do i = 1,npart
         if (iphase(i)==0) then
            j = minloc(abs(itag(1,1:ntag) - iunique(i)), 1)
            if (itag(1,j) - iunique(i) == 0) then
               if (iunique(i) <= n1) then
                  itag(2,j) = i
                  npartnew = npartnew + 1
                  xi = xyzmh(1,i)
                  yi = xyzmh(2,i)
                  zi = xyzmh(3,i)
                  have = have + xyzmh(4,i)
                  if (verbose) print*, i,j, itag(1,j),iunique(i)
               else
                  nbdy = nbdy + 1
                  iphase(i) = 66
               endif
            else
               if (iunique(i) <= n1) then
                  iphase(i) = 99
               else
                  iphase(i) = 66
               endif
            endif
         endif
      enddo
c$omp end parallel do
      have = have/npartnew
      print*, 'h_ave:  ',have
      print*, 'condensing partially complete ',ntag,npartnew,nbdy,nsink
      ntag = ntag - nsink - nbdy
      if (npartnew .ne. ntag) then
         print*, 'ntag .ne. npartnew: ',ntag,npartnew
         stop
      endif   

c     Untagging particles below a given density threshhold that are not casually connected
      if (casually_connected) then
c        Tag the particles within 250au ! rho > 1e-18
         nctr = 0
         nctrnd = 0
c$omp  parallel do default(none)
c$omp& shared(npart,iphase,rho,rho18,rad250,xyz,xyzmh)
c$omp& private(i,rad2)
c$omp& reduction(+:nctr,nctrnd)
         do i = 1,npart
            if (iphase(i)==0) then
               if (.false.) then
                  if (rho(i) < rho18) then
                     nctrnd = nctrnd + 1
                     iphase(i) = 50
                  else
                     nctr = nctr + 1
                  endif
               else
                  rad2 = (xyz(1,1)-xyzmh(1,i))**2
     &                 + (xyz(2,1)-xyzmh(2,i))**2
     &                 + (xyz(3,1)-xyzmh(3,i))**2
                  if (rad2 > rad250*rad250) then
                     nctrnd = nctrnd + 1
                     iphase(i) = 50
                  else
                     nctr = nctr + 1
                  endif
                endif
            endif
         enddo
c$omp end parallel do
         print*, 'There are ',nctr,
     &           ' dense particles associated with this core'
         print*, 'There are ',nctrnd,
     &           ' non-dense particles associated with this core'
c        keep looping until there are no more particles to associate
         if (nctr < 58) then
            print*, 'WARNING! there are not enough dense particles!'
            print*, 'Augmenting to 58 by distance!'
            ictr = 0
            do i = 1,npart
               if (iphase(i)==50) then
                  rad2 = (xyz(1,1)-xyzmh(1,i))**2
     &                 + (xyz(2,1)-xyzmh(2,i))**2
     &                 + (xyz(3,1)-xyzmh(3,i))**2
                  ictr = ictr + 1
                  list(ictr) = i
                  rad2_tmp(ictr) = rad2
               endif
            enddo

            CALL indexx(ictr, list, rad2_tmp, indx)
            ictr = 1
            do while (nctr < 58)
               iphase(list(indx(ictr))) = 0
               nctr = nctr + 1
               ictr = ictr + 1
               nctrnd = nctrnd - 1
            enddo
            print*, 'There are ',nctr,
     &              ' dense particles associated with this core'
            print*, 'There are ',nctrnd,
     &              ' non-dense particles associated with this core'
         endif

         p = 0
         nctrtot = nctr
         rkern   = 2.
         do while (nctr > 0)
c$omp parallel do default(none)
c$omp&shared(npart,ntag0,xyzmh,iphase,itag,rkern)
c$omp&private(i,j,ii,jj,xi,yi,zi,twoh2i,dx,dy,dz,twoh2j,rad2)
            do ii = 1,ntag0
               i = itag(2,ii)
               if (i > 0) then
                  if (iphase(i)==50) then  ! originally tagged, but now unjoined particle
                     xi = xyzmh(1,i)
                     yi = xyzmh(2,i)
                     zi = xyzmh(3,i)
                     twoh2i = (rkern*xyzmh(5,i))**2
                     jj = 1
                     do while (jj <= ntag0)
                        j = itag(2,jj)
                        if (j > 0) then
                           if (iphase(j)==0) then ! joined particle
                              dx = xi - xyzmh(1,j)
                              dy = yi - xyzmh(2,j)
                              dz = zi - xyzmh(3,j)
                              twoh2i = (rkern*xyzmh(5,j))**2
                              rad2 = dx*dx + dy*dy + dz*dz
                              if (rad2<twoh2i .or. rad2<twoh2j) then
                                 iphase(i) = 40
                                 j = npart
                              endif
                           endif
                        endif
                        jj = jj + 1
                     enddo
                  endif
               endif
            enddo
c$omp  end parallel do
            nctr = 0
c           reset the tags
c$omp parallel do default(none)
c$omp&shared(npart,iphase)
c$omp&private(i)
c$omp&reduction(+:nctr)
            do i = 1,npart
               if (iphase(i)==40) then
                  nctr = nctr + 1
                  iphase(i) = 0
               endif
            enddo
c$omp  end parallel do
            nctrtot = nctrtot + nctr
            print*, 'Loop',p,' added ',nctr,
     &              ' particles for a total of ',nctrtot
            p = p + 1

c           We need another criterion for diffuse associated cores.
c           We will require that at least 90% of the associate gas is in the final core
c            if (nctr.eq.0 .and. nctrtot < 0.954*ntag) then
            if (nctr.eq.0 .and. nctrtot < 0.900*ntag) then
c            if (nctr.eq.0 .and. nctrtot < 0.800*ntag) then
c            if (nctr.eq.0 .and. nctrtot < 0.682*ntag) then
               nctr  = 1           ! to kick back into the loop
               rkern = rkern + 2.  ! using additive values of two since this is the cubic kernel
               print*, 'increasing rkern to ',rkern
            else
               rkern = 2.
            endif
         enddo
c        remove tags of unused particles
         nctr = 0
c$omp parallel do default(none)
c$omp&shared(npart,iphase)
c$omp&private(i)
c$omp&reduction(+:nctr)
         do i = 1,npart
            if (iphase(i)==50) then
               nctr = nctr + 1
               !iphase(i) = 99
            endif
         enddo
c$omp  end parallel do
         print*, 'We removed ',nctr,' originally tagged particles'
         ntag = ntag - nctr
      endif

c     Add particles to the list that are within 2h of a tagged particle
      if (find_neigh) then
         ploop = 1
         iloop = 1
         do while (ploop > 0)
c$omp  parallel do default(none)
c$omp& shared(npart,ntag0,itag,xyzmh,iphase,iunique,n1,nextra)
c$omp& shared(randomborder,iloop,ploop,find_neigh50,find_neigh99)
c$omp& private(i,j,jj,xi,yi,zi,hi2,hj2,rtmp2)
         do i = 1,npart
            if (iphase(i)==99 .or. iphase(i)==50) then
               xi = xyzmh(1,i)
               yi = xyzmh(2,i)
               zi = xyzmh(3,i)
               hi2 = xyzmh(5,i)**2
               if (ploop==2 .and. iphase(i)==50) then
                  hi2 = (iloop*xyzmh(5,i))**2
               endif
               jj = 1
               do while (jj <= ntag0)
                  j = itag(2,jj)
                  if (j > 0) then
                     if (iphase(j)==0) then
                        rtmp2 = (xi-xyzmh(1,j))**2
     &                        + (yi-xyzmh(2,j))**2
     &                        + (zi-xyzmh(3,j))**2
                        hj2 = xyzmh(5,j)**2
                        if (randomborder) then
                           if (rtmp2 < max(hi2,hj2)) then
                              iphase(i) = 98
                              j = npart
                           elseif (rtmp2 <  4.*max(hi2,hj2)) then
                              iphase(i) = 97
                           elseif (rtmp2 < 16.*max(hi2,hj2)) then
                              iphase(i) = 96
                           endif
                        else
                           if (rtmp2 < nextra**2*max(hi2,hj2)) then
                              iphase(i) = 98
                           endif
                           j = npart
                        endif
                     endif
                  endif
                  jj = jj + 1
               enddo
            endif
         enddo
c$omp end parallel do
         nadd  = 0
         nadd1 = 0
         nadd2 = 0
         noadd = 0
         iseed = 1701
cc$omp  parallel do default(none)
cc$omp& shared(npart,iphase,iseed)
cc$omp& private(i,rannum)
cc$omp& reduction(+:nadd,nadd1,nadd2,noadd)
         do i = 1,npart
            if (iphase(i)==98) then
               iphase(i) = 0
               nadd = nadd + 1
            elseif (iphase(i)==96 .or. iphase(i)==97) then
               rannum = ran1(iseed)
               if (iphase(i)==97 .and. rannum < 0.5) then !0.25) then !0.06667) then
                  iphase(i) = 0
                  nadd1 = nadd1 + 1
               elseif (iphase(i)==96 .and. rannum < 0.25) then
                  iphase(i) = 0
                  nadd2 = nadd2 + 1
               else
                  iphase(i) = 20
                  noadd = noadd + 1
               endif
            endif
         enddo
cc$omp end parallel do
         print*, 'find_neigh: ',nadd,nadd1,nadd2,noadd,ntag
         ntag = ntag + nadd + nadd1 + nadd2
         ploop = ploop - 1
         if (ntag < ntag0) then
            print*, 'need more particles! ',ploop,iloop
            if (ploop==0) then
               ploop = 2
               iloop = iloop + 1
            endif
         endif
         enddo
      endif

c-----------------------------------------------------------------------
c--   Determine values (rho,B,u) the border particles
c-----------------------------------------------------------------------

      if (.true.) then
         print*, 'Determining values for background medium'
         ! Detemine the number of neighbours so that we can asses the particles near the boundary of the extracted cloud
         ibdy = 0
         nbdy = 0
         Bxave = 0.
         Byave = 0.
         Bzave = 0.
         nBave = 0
         uave  = 0.
         rhoave  = 0.
         rhoavel = 0.
c$omp  parallel do default(none)
c$omp& shared(npart,xyzmh,vxyzu,Bxyz,rho,iphase,ibdy)
c$omp& shared(dxbound,dybound,dzbound,hdbx,hdby,hdbz)
c$omp& private(i,j,xi,yi,zi,hi,twohi2,dx,dy,dz,nctr,rad2)
c$omp& reduction(+:nbdy,Bxave,Byave,Bzave,uave,rhoave,rhoavel)
         do i = 1,npart
            if (iphase(i)==0) then
               nctr = 0
               xi = xyzmh(1,i)
               yi = xyzmh(2,i)
               zi = xyzmh(3,i)
               hi = xyzmh(5,i)
               twohi2 = 4.0*hi*hi
               do j = 1,npart
                  if (iphase(j)==0 .and. i.ne.j) then
                     dx = xi - xyzmh(1,j)
                     dy = yi - xyzmh(2,j)
                     dz = zi - xyzmh(3,j)
                     if (abs(dx)>hdbx) dx = dx-dxbound*SIGN(1.0,dx)
                     if (abs(dy)>hdby) dy = dy-dybound*SIGN(1.0,dy)
                     if (abs(dz)>hdbz) dz = dz-dzbound*SIGN(1.0,dz)
                     rad2 = dx*dx + dy*dy + dz*dz
                     if (rad2 < twohi2) nctr = nctr + 1
                  endif
               enddo
               ibdy(i) = nctr
               if (nctr < 30) then
                  Bxave  = Bxave  + Bxyz( 1,i)
                  Byave  = Byave  + Bxyz( 2,i)
                  Bzave  = Bzave  + Bxyz( 3,i)
                  uave   = uave   + vxyzu(4,i)
                  rhoave = rhoave + rho(i)
                  rhoavel = rhoavel + log10(rho(i))
                  nbdy   = nbdy + 1
               endif
            endif
         enddo
c$omp end parallel do
         print*, 'There are ',nbdy,' boundary particles identified'

         if (nbdy > 0) then
            Bxave  = Bxave/nbdy
            Byave  = Byave/nbdy
            Bzave  = Bzave/nbdy
            rhoave = rhoave/nbdy
            rhoavel = 10**(rhoavel/nbdy)
            uave   = uave/nbdy
            print*, 'The boundary B-field is ',Bxave,Byave,Bzave
            print*, 'The boundary density is ',rhoave,rhoavel
            print*, 'The boundary thermal energy is ',uave
         else
            print*, 'FATAL! there are no particles associated with ',
     &              'the boundary!'
            do i = 1,npart
               write(666,*) i,ibdy(i),xyzmh(1:3,i)
            enddo
            stop
         endif

c        the background values
         rhobkg     = rhoavel/frho
         ubkg       = uave*frho
         uzero_n2   = ubkg
         pmass      = xyzmh(4,1)
         dxyz       = (pmass/rhobkg)**(1./3.)
         hbkg       = 1.2*dxyz
         twoh2bkg   = 4.0*hbkg*hbkg
         Bxyzcom(1) = Bxave
         Bxyzcom(2) = Byave
         Bxyzcom(3) = Bzave

c        Tag all particle with B > 5Bave within hbkg
         ibdy  = 0 ! retask to avoid duplicate calculations
         B2bkg = dot_product(Bxyzcom,Bxyzcom)
         nadd  = 1 ! to enter the loop

         nadd = 0 ! to turn off loop!

         do while (nadd > 0)
c$omp  parallel do default(none)
c$omp& shared(npart,xyzmh,iphase,Bxyz,B2bkg,ibdy)
c$omp& shared(twoh2bkg,hdbx,hdby,hdbz,dxbound,dybound,dzbound)
c$omp& private(i,j,xi,yi,zi,dx,dy,dz,rad2,B2i)
            do i = 1,npart
               if (iphase(i)==0 .and. ibdy(i) == 0) then
                  ibdy(i) = 1
                  xi = xyzmh(1,i)
                  yi = xyzmh(2,i)
                  zi = xyzmh(3,i)
                  B2i = dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
                  if (B2i > 25.*B2bkg) then
                     do j = 1,npart
                        if (iphase(j) == 99) then
                           dx = xi - xyzmh(1,j)
                           dy = yi - xyzmh(2,j)
                           dz = zi - xyzmh(3,j)
                           if (abs(dx)>hdbx) dx=dx-dxbound*SIGN(1.0,dx)
                           if (abs(dy)>hdby) dy=dy-dybound*SIGN(1.0,dy)
                           if (abs(dz)>hdbz) dz=dz-dzbound*SIGN(1.0,dz)
                           rad2 = dx*dx + dy*dy + dz*dz
                           if (rad2 < twoh2bkg) iphase(j) = 98
                        endif
                     enddo
                  endif
               endif
            enddo
c$omp end parallel do
            nadd = 0
c$omp  parallel do default(none)
c$omp& shared(npart,iphase)
c$omp& private(i)
c$omp& reduction(+:nadd)
            do i = 1,npart
               if (iphase(i)==98) then
                  iphase(i) = 0
                  nadd = nadd + 1
               endif
            enddo
c$omp end parallel do
            ntag = ntag + nadd
            print*, 'find_StrongB: ',nadd,ntag
         enddo


#ifdef NONIDEAL
c        Find background non-ideal coefficients for an initial guess
         keepsearching = .true.
         i = n1
         do while (keepsearching .and. i <= npart)
            nimhdi = eta_nimhd(:,i)
            if (iphase(i)==66) keepsearching = .false.
            i = i + 1
         enddo
#endif
      endif

c-----------------------------------------------------------------------
c--   Condense particle list to only inlclude iphase = 0 particles
c-----------------------------------------------------------------------
c     Print a list of particles, if required
      if (print_iunique) then
         write(fileiunique,'(a,I2.2,2a)') 
     &   'core_augmented',isink,'_dmax',trim(dat)
         open(unit=12664,file=trim(fileiunique))
         do i = 1,npart
            if (iphase(i)==0) write(12664,*) iunique(i)
         enddo
         close(12664)
      endif
    

c     We will manually copy the relevant arrays; the remainder should be automatically updated upon restart
      i = 1
      j = 1
      rhotot  = 0.
      xyzcom  = 0.
      vxyzcom = 0.

      do while (i <= npart .and. j < npart)
         if (iphase(i)==0) then
            i = i + 1
         else
c           an irrelevant particle we want to replace
            j = max(j,i + 1)
            do while (iphase(j).ne.0 .and. j<npart)
               j = j + 1
            enddo
            !print*, i,j,npart, iphase(j)
            if (iphase(j)==0) then ! likely necessary to deal with npart
c              Replace!
               if (verbose) print*, 'replacing', i, j
               xyzmh(:,i)    = xyzmh(:,j)
               vxyzu(:,i)    = vxyzu(:,j)
               ekcle(:,i)    = ekcle(:,j)
               Bxyz(:,i)     = Bxyz(:,j)
               Bevolxyz(:,i) = Bevolxyz(:,j)
               rho(i)        = rho(j)
               dgrav(i)      = dgrav(j)
               gradhs(:,i)   = gradhs(:,j)
               alphaMM(:,i)  = alphaMM(:,j)
               heatingISR(:,i) = heatingISR(:,j)
               h2frac(i)     = h2frac(j)
               dh2dt(:,i)    = dh2dt(:,j)
               dq(i)         = dq(j)
               poten(i)      = poten(j)
#ifdef NONIDEAL
               eta_nimhd(:,i)= eta_nimhd(:,j)
#endif
               isteps(i)     = 0 ! isteps(j)
               iunique(i)    = i            ! might as well reset
               isort(i)      = i
               iorig(i)      = i
               iphase(i)     = 0
               iphase(j)     = 1            ! since this particle has been relocated
c              find average properties
               rhotot        = rhotot  + rho(i)
               xyzcom        = xyzcom  + xyzmh(1:3,i)*rho(i)
               vxyzcom       = vxyzcom + vxyzu(1:3,i)*rho(i)

c              update index
               i = i + 1
            endif
         endif
      enddo
      print*, 'condensing complete',i-1,ntag
      xyzcom  =  xyzcom/rhotot
      vxyzcom = vxyzcom/rhotot

c     Update particle numbers
      npart   = ntag
      nactive = ntag
      n1      = ntag
      n2      = npart - n1
      iuniquemax = npart

c     Recentre particles & determine new borders
      if (.not. recentre) xyzcom = 0.
      xmin0 = xmax
      xmax0 = xmin
      ymin0 = ymax
      ymax0 = ymin
      zmin0 = zmax
      zmax0 = zmin
      write(*,wfmt) 'before: xmin0: ',xmin0,ymin0,zmin0
      write(*,wfmt) 'before: xmax0: ',xmax0,ymax0,zmax0
      write(*,wfmt) 'before: xmin:  ',xmin,ymin,zmin
      write(*,wfmt) 'before: xmax:  ',xmax,ymax,zmax
c$omp  parallel do default(none)
c$omp& shared(npart,xyzmh,vxyzu,vxyzcom,xyzcom)
c$omp& private(i)
c$omp& reduction(min:xmin0,ymin0,zmin0)
c$omp& reduction(max:xmax0,ymax0,zmax0)
      do i = 1,npart
         vxyzu(1:3,i) = vxyzu(1:3,i) - vxyzcom
         xyzmh(1:3,i) = xyzmh(1:3,i) -  xyzcom
         xmin0 = min(xmin0,xyzmh(1,i))
         xmax0 = max(xmax0,xyzmh(1,i))
         ymin0 = min(ymin0,xyzmh(2,i))
         ymax0 = max(ymax0,xyzmh(2,i))
         zmin0 = min(zmin0,xyzmh(3,i))
         zmax0 = max(zmax0,xyzmh(3,i))
      enddo
c$omp end parallel do
      write(*,wfmt) 'during: xmin0: ',xmin0,ymin0,zmin0
      write(*,wfmt) 'during: xmax0: ',xmax0,ymax0,zmax0
      write(*,wfmt) 'during: xmin:  ',xmin,ymin,zmin
      write(*,wfmt) 'during: xmax:  ',xmax,ymax,zmax
      ffac = 0.5 !1.0 !0.5
      xmin = xmin0 - ffac*(xmax0-xmin0)
      xmax = xmax0 + ffac*(xmax0-xmin0)
      ymin = ymin0 - ffac*(ymax0-ymin0)
      ymax = ymax0 + ffac*(ymax0-ymin0)
      zmin = zmin0 - ffac*(zmax0-zmin0)
      zmax = zmax0 + ffac*(zmax0-zmin0)
      write(*,wfmt) 'after: xmin0: ',xmin0,ymin0,zmin0
      write(*,wfmt) 'after: xmax0: ',xmax0,ymax0,zmax0
      write(*,wfmt) 'after: xmin:  ',xmin,ymin,zmin
      write(*,wfmt) 'after: xmax:  ',xmax,ymax,zmax

c-----------------------------------------------------------------------
c     Calculate properties of the cloud (diagnostic only)
c-----------------------------------------------------------------------
      if (get_stats) then
c---     Bulk properties
         Ekin    = 0.
         Epot    = 0.
         Emag    = 0.
         Etherm  = 0
         rhoavel = 0.
         Bavel   = 0.
         mach    = 0.
         Ta      = 0.
         Ba      = 0.
         Bzta    = 0.
         Bztheta = 0.
         list    = 0
         vave    = 0.
         np      = npart ! a shorter variable
         if (npart > ncore_max) then
            print*, 'npart = ',npart
            print*, 'ncore_max = ',ncore_max
            print*, 'recompile with larger ncore_max'
            stop
         endif
c$omp  parallel do default(none)
c$omp& shared(npart,xyzmh,vxyzu,Bxyz,ekcle,rho,Bi,uergg,list,Bztheta)
c$omp& shared(Temp)
c$omp& private(i,j,v2,B2,r2,dr,cs2)
c$omp& reduction(+:Ekin,Etherm,Epot,Emag,rhoavel,Bavel,mach,Ba,Bzta,Ta)
c$omp& reduction(+:vave)
         do i = 1,npart
            list(i)    = i
            v2         = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
            B2         = dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
            Bi(i)      = sqrt(B2)
            Temp(i)    = vxyzu(4,i)/ekcle(3,i)
            if (Bi(i).gt.0.0) then
               Bztheta(i) = Bxyz(3,i)/Bi(i) ! the cosine of the angle of the magnetic field with respect to z = (0,0,1)
            endif
            cs2        = 5.0*Rg*vxyzu(4,i)/(3.0*ekcle(3,i)*gmw*uergg)
            Ekin       = Ekin + 0.5*xyzmh(4,i)*v2
            Emag       = Emag + 0.5*xyzmh(4,i)*B2/rho(i)
            Etherm     = Etherm + xyzmh(4,i)*vxyzu(4,i)
            rhoavel    = rhoavel + log10(rho(i))
            Bavel      = Bavel   + log10(Bi(i))
            vave(1)    = vave(1)   + sqrt(v2)
            vave(2:4)  = vave(2:4) + vxyzu(1:3,i)
            mach       = mach    + v2/cs2
            Ba         = Ba      + Bxyz(1:3,i)
            Bzta       = Bzta    + Bztheta(i)
            Ta         = Ta      + Temp(i)
            do j = i+1, npart
               dr   = xyzmh(1:3,i)-xyzmh(1:3,j)
               r2   = dot_product(dr,dr)
               Epot = Epot + xyzmh(4,i)*xyzmh(4,j)/sqrt(r2)
            enddo
         enddo
c$omp end parallel do
         mach = sqrt(mach/npart)
         vave = vave/npart
         print*, 'Kinetic energy: ',Ekin*uerg
         print*, 'Thermal energy: ',Etherm*uerg
         if (.not.set_Beq0) print*, 'Magnetic energy: ',Emag*uerg
         print*, 'Potential energy: ',Epot*uerg
         print*, 'log-averaged density: ',10**(rhoavel/npart)*udensity
         print*, 'log-averaged B: ',10**(Bavel/npart)*umagfdi
         print*, 'mach number: ',mach

         CALL get_spread(np,list(1:np),real(rho(1:np)),
     &                   rhorange,rhoavel,0)
         rhorange = rhorange*udensity
         CALL get_spread(np,list(1:np),Temp(1:np),Trange,Ta,1)
         CALL get_spread(np,list(1:np),Bi(1:np),Brange(1,:),Bavel,0)
         CALL get_spread(np,list(1:np),Bxyz(1,1:np),Brange(2,:),
     &                   Ba(1),1)
         CALL get_spread(np,list(1:np),Bxyz(2,1:np),Brange(3,:),
     &                   Ba(2),1)
         CALL get_spread(np,list(1:np),Bxyz(3,1:np),Brange(4,:),
     &                   Ba(3),1)
         Brange = Brange*umagfdi
         CALL get_spread(np,list(1:np),Bztheta(1:np),Bztrange,Bzta,1)

c---     Find the core values
c        find the densest particle
         irho    = 1
         do i = 2,npart
            if (rho(i) > rho(irho)) irho = i
         enddo
         pmass    = xyzmh(4,irho)
         hcore2   = (2.*1.2*(pmass/rho(irho))**(1./3.))**2
         if (.false.) then
c           take the position & velocity of the densest particle
            xcore(1) = xyzmh(1,irho)
            xcore(2) = xyzmh(2,irho)
            xcore(3) = xyzmh(3,irho)
            vcore(1) = vxyzu(1,irho)
            vcore(2) = vxyzu(2,irho)
            vcore(3) = vxyzu(3,irho)
         else
c           take the position & velocity using the density-weighted centre of values of the particles within 2h of the densest particle
            xcore = 0.
            vcore = 0.
            rcore = 0.
            ncore = 0
c$omp  parallel do default(none)
c$omp& shared(npart,xyzmh,vxyzu,hcore2,rho,irho)
c$omp& private(i,dx,dy,dz,rtmp2)
c$omp& reduction(+:xcore,vcore,rcore,ncore)
            do i = 1,npart
               dx  = xyzmh(1,irho)-xyzmh(1,i)
               dy  = xyzmh(2,irho)-xyzmh(2,i)
               dz  = xyzmh(3,irho)-xyzmh(3,i)
               rtmp2 = dx**2 + dy**2 + dz**2
               if (rtmp2 < hcore2) then
                  xcore = xcore + xyzmh(1:3,i)*rho(i)
                  vcore = vcore + vxyzu(1:3,i)*rho(i)
                  rcore = rcore +              rho(i)
                  ncore = ncore + 1
               endif
            enddo
c$omp end parallel do
            xcore = xcore / rcore
            vcore = vcore / rcore
            print*, 'COM',rcore, ncore,hcore2,pmass
         endif

c---     Angular momentum
c        Assume the densest particle is the 'centre of momentum'
c        calculate the angular momentum & magnetic field of the core (i.e. all particles within 2h of the densest particle) & compare angles
c        calculate the total angular momentum of the cloud
c        on a particle-by-particle basis, determine the angle between B & L
c        Also calculate the velocity dispersion
         Lcore   = 0.
         Bcore   = 0.
         ncore   = 0
         Lcloud  = 0.
         BLtheta = 0.
         BLta    = 0.
         Erot    = 0.
         veldisp  = 0.
c$omp  parallel do default(none)
c$omp& shared(npart,xyzmh,vxyzu,Bxyz,BLtheta,pmass,hcore2,xcore,vcore)
c$omp& shared(vave)
c$omp& private(i,dx,dy,dz,dvx,dvy,dvz,Lxyz,rtmp2,v2,r2xy,r2xz,r2yz)
c$omp& reduction(+:Lcloud,Lcore,Bcore,ncore,BLta,veldisp,Erot)
         do i = 1,npart
            dx  = xcore(1)-xyzmh(1,i)
            dy  = xcore(2)-xyzmh(2,i)
            dz  = xcore(3)-xyzmh(3,i)
            dvx = vcore(1)-vxyzu(1,i)
            dvy = vcore(2)-vxyzu(2,i)
            dvz = vcore(3)-vxyzu(3,i)
            v2  = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))

            Lxyz(1) = dy*dvz - dz*dvy
            Lxyz(2) = dz*dvx - dx*dvz
            Lxyz(3) = dx*dvy - dy*dvx
            Lcloud  = Lcloud + Lxyz*pmass

            r2xy = dx*dx + dy*dy
            r2xz = dx*dx + dz*dz
            r2yz = dy*dy + dz*dz
            if (r2yz > 0.) Erot(1) = Erot(1) + pmass*Lxyz(1)**2/r2yz
            if (r2xz > 0.) Erot(2) = Erot(2) + pmass*Lxyz(2)**2/r2xz
            if (r2xy > 0.) Erot(3) = Erot(3) + pmass*Lxyz(3)**2/r2xy

            BLtheta(i) = costheta(Lxyz,Bxyz(1:3,i))
            BLta       = BLta + BLtheta(i)

            veldisp(1)   = veldisp(1)   + (sqrt(v2)     - vave(1))  **2
            veldisp(2:4) = veldisp(2:4) + (vxyzu(1:3,i) - vave(2:4))**2
            rtmp2 = dx**2 + dy**2 + dz**2
            if (rtmp2 < hcore2) then
               Lcore = Lcore + Lxyz
               Bcore = Bcore + Bxyz(1:3,i)
               ncore = ncore + 1
            endif
         enddo
c$omp end parallel do
         veldisp = sqrt(veldisp/npart)
         CALL get_spread(np,list(1:np),BLtheta(1:np),BLtrange,BLta,1)
         BLcore = costheta(Lcore,Bcore)

c        Get the size of the best-fit ellipse that describes the cloud
         lctr = 0
         find_axes = .true.
         xyztmp(1:5,1:npart) = xyzmh(1:5,1:npart)
         do while (lctr < 1000 .and. find_axes)
           lctr = lctr + 1
           find_axes = .false.
           call fit_ellipse(np,1,2,xyztmp(:,1:np),udistpc,ell(1,:),phi)
           if (phi > 1.e-3) find_axes = .true.
           call fit_ellipse(np,1,3,xyztmp(:,1:np),udistpc,ell(2,:),phi)
           if (phi > 1.e-3) find_axes = .true.
           call fit_ellipse(np,2,3,xyztmp(:,1:np),udistpc,ell(3,:),phi)
           if (phi > 1.e-3) find_axes = .true.
         enddo
         imin = minloc(ell(:,1),1)
         imax = maxloc(ell(:,1),1)
         cloudxyz(1) = ell(imax,1)*0.5  ! The semi-major axis
         cloudxyz(2) = ell(imin,1)*0.5  ! The semi-middle axis
         imin = minloc(ell(:,2),1)
         imax = maxloc(ell(:,2),1)
         ell(imin,2) = ell(imax,2)
         imin = minloc(ell(:,2),1)
         cloudxyz(3) = ell(imin,2)*0.5  ! The semi-minor axis

c        Print to file
         write_header = .true.
         if (indiv_out) then
            write(filestats,'(a,I2.2,a)') 'core_Stats',isink,'.dat'
            open(unit=67,file=trim(filestats))
         else
            write(filestats,'(a)') 'core_Stats.dat'
            inquire(file=filestats,exist=fexists)
            if (fexists) then
               write_header = .false.
               open(unit=67,file=trim(filestats),position='append')
            else
               open(unit=67,file=trim(filestats))
            endif
         endif
         if (write_header) then
            write(67,"(2a)") '# time in code units; M in Msun; ',
     &      'cloud semi axes in pc; rest in cgs'
            write(67,"('#',86(1x,'[',i2.2,1x,a11,']',2x))")
     &                     1,'isink',
     &                     2,'n core',
     &                     3,'cloud form',
     &                     4,'sink form',
     &                     5,'end time',
     &                     6,'cloud age',
     &                     7,'sink age',
     &                     8,'sink mass_f',
     &                     9,'M_tagged',
     &                    10,'M_total',
     &                    11,'Ekin',
     &                    12,'Etherm',
     &                    13,'Emag',
     &                    14,'Epot',
     &                    15,'Mach no',
     &                    16,'rho 0.2%',
     &                    17,'rho 2.3%',
     &                    18,'rho 15.9% ',
     &                    19,'rho logave',
     &                    20,'rho 84.1%',
     &                    21,'rho 97.7%',
     &                    22,'rho 99.8%',
     &                    23,'B 0.1%',
     &                    24,'B 2.3%',
     &                    25,'B 15.9%',
     &                    26,'B logave',
     &                    27,'B 84.1%',
     &                    28,'B 97.7%',
     &                    29,'B 99.8%',
     &                    30,'Bx 0.1%',
     &                    31,'Bx 2.3%',
     &                    32,'Bx 15.9%',
     &                    33,'Bx ave',
     &                    34,'Bx 84.1%',
     &                    35,'Bx 97.7%',
     &                    36,'Bx 99.8%',
     &                    37,'By 0.1%',
     &                    38,'By 2.3%',
     &                    39,'By 15.9%',
     &                    40,'By ave',
     &                    41,'By 84.1%',
     &                    42,'By 97.7%',
     &                    43,'By 99.8%',
     &                    44,'Bz 0.1%',
     &                    45,'Bz 2.3%',
     &                    46,'Bz 15.9%',
     &                    47,'Bz ave',
     &                    48,'Bz 84.1%',
     &                    49,'Bz 97.7%',
     &                    50,'Bz 99.8%',
     &                    51,'BZcos 0.1%',
     &                    52,'BZcos 2.3%',
     &                    53,'BZcos 15.9%',
     &                    54,'BZcos ave',
     &                    55,'BZcos 84.1%',
     &                    56,'BZcos 97.7%',
     &                    57,'BZcos 99.8%',
     &                    58,'Lx cloud',
     &                    59,'Ly cloud',
     &                    60,'Lz cloud',
     &                    61,'BLcos 0.1%',
     &                    62,'BLcos 2.3%',
     &                    63,'BLcos 15.9%',
     &                    64,'BLcos ave',
     &                    65,'BLcos 84.1%',
     &                    66,'BLcos 97.7%',
     &                    67,'BLcos 99.8%',
     &                    68,'Lx core',
     &                    69,'Ly core',
     &                    70,'Lz core',
     &                    71,'BL core',
     &                    72,'Cloud x',
     &                    73,'Cloud y',
     &                    74,'Cloud z',
     &                    75,'vel disp',
     &                    76,'vel dispX',
     &                    77,'vel dispY',
     &                    78,'vel dispZ',
     &                    79,'Erot',
     &                    80,'Temp 0.1%',
     &                    81,'Temp 2.3%',
     &                    82,'Temp 15.9%',
     &                    83,'Temp ave',
     &                    84,'Temp 84.1%',
     &                    85,'Temp 97.7%',
     &                    86,'Temp 99.8%'

         endif

         write(67,'(2(I18,1x),84(Es18.6,1x))') isink,ncore,
     &      coreform,sinkform,timeend,
     &      coreage,sinkage,sinkmass,ntag0*pmass,npart*pmass,
     &      Ekin*uerg,Etherm*uerg,Emag*uerg,Epot*uerg,mach,
     &      rhorange,Brange(1,:),Brange(2,:),Brange(3,:),Brange(4,:),
     &      Bztrange,Lcloud*uang,BLtrange,Lcore*uang,BLcore,
     &      cloudxyz*udistpc,veldisp*uvel,
     &      0.5*sqrt(dot_product(Erot,Erot))*uerg,Trange
         close(67)

c        Create histograms of density, magnetic field strength & angle
         print*, 'calling histograms'
         call get_histo( 32,npart,pmass,real(rho(1:np))*udensity
     &                     ,Bi(1:np)*umagfdi,Bztheta(1:np),Temp(1:np))
         call get_histo( 64,npart,pmass,real(rho(1:np))*udensity
     &                     ,Bi(1:np)*umagfdi,Bztheta(1:np),Temp(1:np))
         call get_histo(128,npart,pmass,real(rho(1:np))*udensity
     &                     ,Bi(1:np)*umagfdi,Bztheta(1:np),Temp(1:np))
         call get_histo(256,npart,pmass,real(rho(1:np))*udensity
     &                     ,Bi(1:np)*umagfdi,Bztheta(1:np),Temp(1:np))
         call get_histo(512,npart,pmass,real(rho(1:np))*udensity
     &                     ,Bi(1:np)*umagfdi,Bztheta(1:np),Temp(1:np))
      endif

c-----------------------------------------------------------------------
c--  Update additional values & arrays
c-----------------------------------------------------------------------
      nlistinactive = 0
      timeend = timeend - gt
      if (recentre) gt = 0. ! reset time

c-----------------------------------------------------------------------
c--   Add background medium
c-----------------------------------------------------------------------
      ekcle3 = 1.5*Rg/uergg  ! getcv(rhobkg,ubkg)
      !print*, 'prekappa',ubkg,ekcle3,rhobkg,frho,Rg,uergg 
      ekcle2 = getkappa(ubkg,ekcle3,real(rhobkg,kind=4))
      ekcle1 = radconst/uergcc*(ubkg/ekcle3)**4/rhobkg
      if (write_bkg) then
         write(debugfile,'(a,I2.2,a)') 'fort',isink,'.1233'
         open(unit=1233,file=trim(debugfile))
      endif
      pmax = 1
      n2 = 0
      do p = 1,pmax
         zi = zmin
         kk = 0
         do while (zi < zmax)
            yi = ymin
            jj = 0
            do while (yi < ymax)
               xi = xmin
               ii = 0
               do while (xi < xmax)
                  hbkg0 = hbkg
                  if (nooverlap .or. mod_h) then
                     addparticle = .true.
                     hbkg0 = -1.
c$omp  parallel do default(none)
c$omp& shared(npart,n1,xi,yi,zi,xyzmh,addparticle)
c$omp& private(i,rad2)
c$omp& reduction(max:hbkg0)
                     do i = 1,n1
                        rad2 = (xi-xyzmh(1,i))**2
     &                       + (yi-xyzmh(2,i))**2
     &                       + (zi-xyzmh(3,i))**2
                        if (rad2 < xyzmh(5,i)**2) then   ! ignore if within h of i
c                        if (rad2 < 0.25*xyzmh(5,i)**2) then   ! ignore if within h/2 of i
c                        if (rad2 < 4.0*xyzmh(5,i)**2) then   ! ignore if within 2h of i
                           addparticle = .false. ! acceptable race condition
                        endif
                        if (rad2 < 4.*xyzmh(5,i)**2) then
                           hbkg0 = max(hbkg0,xyzmh(5,i))
                        endif
                     enddo
c$omp end parallel do
                     if (.not. nooverlap) addparticle = .true.
                     if (mod_h) then
                        if (hbkg0 > 0.) then
                           hbkg0 = min(1.5*hbkg0,hbkg)
                        else
                           hbkg0 = hbkg
                        endif
                     else
                        hbkg0 = hbkg
                     endif
                  else
                     addparticle = .true.
                  endif
                  if (addparticle) then
                     npart = npart + 1
                     xyzmh(1,npart) = xi
                     xyzmh(2,npart) = yi
                     xyzmh(3,npart) = zi
                     xyzmh(4,npart) = pmass
                     xyzmh(5,npart) = hbkg0
                     vxyzu(1:3,npart) = 0.0
                     vxyzu(4,npart)   = ubkg
                     ekcle(:,npart)   = 0. ! dummy value
                     ekcle(3,npart)   = 1.462E-01
                     ekcle(2,npart)   = 3.186E-02
                     ekcle(1,npart)   = 4.292E-03
                     Bxyz(:,npart)    = Bxyzcom
                     Bevolxyz(1:3,npart) = Bxyzcom/rhobkg
                     rho(npart)          = rhobkg
                     dgrav(npart)        = 0.
                     gradhs(1,npart)     = 1.019E+00
                     gradhs(2,npart)     = -5.657E-04
                     alphaMM(1,npart)    = 0.1
                     alphaMM(2,npart)    = 0
                     heatingISR(:,npart) = 0. ! from icl0000, these are small, non-zero numbers with no pattern
                     h2frac(npart)       = 0.5
                     dh2dt(1,npart)      = -8.925E-05
                     dh2dt(2,npart)      = 0.
                     dq(npart)           = 0.
                     poten(npart)        = 0.
#ifdef NONIDEAL
                     eta_nimhd(:,npart)= nimhdi
#endif
                     iunique(npart)    = npart
                     isort(npart)      = npart
                     iorig(npart)      = npart
                     iphase(npart)     = 0
                     isteps(npart)     = 0
                     if (write_bkg) write(1233,*) xi,yi,zi
                  endif
                  xi = xi + dxyz
                  ii = ii + 1
               enddo
               yi = yi + dxyz
               jj = jj + 1
            enddo
            zi = zi + dxyz
            kk = kk + 1
            if (write_bkg) write(1233,*) ' '
         enddo
         if (p==1) then
            print*,'number of active particles',npart
         else
            print*,'number of total, inactive particles',npart,npart-n1
         endif
      enddo
      if (write_bkg) close(1233)
      nactive = npart
      n2      = npart - n1
      iuniquemax = npart
      write(*,wfmt) 'before final correction: xmin: ',xmin,ymin,zmin
      write(*,wfmt) 'before final correction: xmax: ',xmax,ymax,zmax
      xmin = xmin - 0.5*dxyz
      ymin = ymin - 0.5*dxyz
      zmin = zmin - 0.5*dxyz
      xmax = xmin + ii*dxyz
      ymax = ymin + jj*dxyz
      zmax = zmin + kk*dxyz
      write(*,wfmt) 'after  final correction: xmin: ',xmin,ymin,zmin
      write(*,wfmt) 'after  final correction: xmax: ',xmax,ymax,zmax
      dxbound = xmax - xmin
      dybound = ymax - ymin
      dzbound = zmax - zmin
      hdbx = 0.5*dxbound
      hdby = 0.5*dybound
      hdbz = 0.5*dzbound
c     verification loop
      do i = 1,npart
         if (xyzmh(1,i) < xmin) print*, 'xi < xmin for i = ',i
         if (xyzmh(1,i) > xmax) print*, 'xi > xmax for i = ',i
         if (xyzmh(2,i) < ymin) print*, 'yi < ymin for i = ',i
         if (xyzmh(2,i) > ymax) print*, 'yi > ymax for i = ',i
         if (xyzmh(3,i) < zmin) print*, 'zi < zmin for i = ',i
         if (xyzmh(3,i) > zmax) print*, 'zi > zmax for i = ',i
      enddo
c     remove magnetic fields, if requested
      if (set_Beq0) then
         print*, 'Removing magnetic fields from this dump'
         Bxyz     = 0.
         Bevolxyz = 0.
      endif

c-----------------------------------------------------------------------
c     Overwrite magentic background, using the local values of B for a smooth transition
c-----------------------------------------------------------------------     
      if (.false.) then
         print*, 'overwriting magnetic field or pressure'
         ! initialise tags of particle we need to update
         do i = n1+1,npart
            iunique(i) = -iunique(i)
            iphase(i)  = -1
         enddo
         twoh2bkg = 4.0*hbkg*hbkg
         twohbkg1 = 0.5/hbkg
         updateB  = .true.
         ffac = 0.5
         nt0 = npart - n1
         do while (updateB)
            nupdate = 0
            nmin = npart
            nmax = 0
            nave = 0
c$omp  parallel do default(none)
c$omp& shared(npart,n1,xyzmh,iphase,iunique,Bxyz,Bevolxyz,vxyzu,frho)
c$omp& shared(rhobkg,hdbx,hdby,hdbz,ffac,Bbkgave,ubkgave,twoh2bkg)
c$omp& shared(dxbound,dybound,dzbound,twohbkg1)
c$omp& private(i,j,xi,yi,zi,dx,dy,dz,hj,nctr,nneigh,q,wk,rad2,wctr)
c$omp& reduction(min:nmin)
c$omp& reduction(max:nmax)
c$omp& reduction(+:nave,nupdate,Bxave,Byave,Bzave,nBave)
            do i = n1+1,npart
               if (iphase(i)==-1) then ! particle needs updating
                  xi = xyzmh(1,i)
                  yi = xyzmh(2,i)
                  zi = xyzmh(3,i)
                  Bxyz(:,i)  = 0.
                  wctr   = 0.
                  nctr   = 0
                  nneigh = 0
                  ! find all potential sources
                  do j = 1,npart
                     dx = xi - xyzmh(1,j)
                     dy = yi - xyzmh(2,j)
                     dz = zi - xyzmh(3,j)
                     if (abs(dx)>hdbx) dx = dx-dxbound*SIGN(1.0,dx)
                     if (abs(dy)>hdby) dy = dy-dybound*SIGN(1.0,dy)
                     if (abs(dz)>hdbz) dz = dz-dzbound*SIGN(1.0,dz)
                     rad2 = dx*dx + dy*dy + dz*dz
                     if (rad2 < twoh2bkg) then
                        nneigh = nneigh + 1
                        if (iphase(j)==0) then ! particle can be source
                           q  = sqrt(rad2)*twohbkg1
                           wk = wkern(q)
                           wctr = wctr + wk
                           nctr = nctr + 1
                           Bxyz(1:3,i) = Bxyz(1:3,i) + Bxyz(1:3,j)*wk
                        endif
                     endif
                  enddo
                  if (nctr > int(ffac*nneigh)) then ! update B if there are sources
                     !print*, i,nctr,nneigh
                     nmin = min(nmin,nctr)
                     nmax = max(nmax,nctr)
                     nave = nave + nctr
                     nupdate = nupdate + 1
                     Bxyz(1:3,i)     = Bxyz(1:3,i)/wctr
                     Bevolxyz(1:3,i) = Bxyz(1:3,i)/rhobkg
                     iunique(i) = abs(iunique(i))    ! tag so can be identified as source candidate
                  endif
               endif
            enddo
c$omp end parallel do
            if (nupdate > 0) then
               nt0 = nt0 - nupdate
               write(*,'(a,5I8)')
     &         'resetting tags and starrting new loop:'
     &         ,nt0,nupdate,nmin,nave/nupdate,nmax
            endif
            ! reset tags
            updateB = .false.
            do i = n1+1,npart
               if (iunique(i) > 0) then
                  iphase(i) = 0
               else
                  updateB = .true. ! still have particles to correct
               endif
            enddo
            if (updateB .and. nupdate==0) then
               ffac = ffac * 0.9
               print*, 'Decreasing ffac to ',ffac
               if (ffac < 0.01) then
                  print*, 'ffac too small.  aborting'
                  stop
               endif
            endif
         enddo
         print*, 'Done updating B-field'
      endif

c-----------------------------------------------------------------------
c     Sanity checks & verifications
c-----------------------------------------------------------------------
c     Another sanity check
      nctr = 0
c$omp parallel do default(none)
c$omp&shared(npart,iphase)
c$omp&private(i)
c$omp&reduction(+:nctr)
      do i = 1,npart
         if (iphase(i) .ne. 0) nctr = nctr + 1
      enddo
c$omp  end parallel do
      if (nctr==0) then
         print*, 'Sanity check passed; iphase==0'
      else
         print*, 'Sanity check failed.  There are ',nctr,
     &           ' particles with iphase /= 0'
      endif


       if (.false.) then
         ! determine the minimum distance from an active particle to a background particle
         rmin = huge(rmin(1))
c$omp  parallel do default(none)
c$omp& shared(npart,n1,xyzmh,iphase,rmin)
c$omp& shared(hdbx,hdby,hdbz,dxbound,dybound,dzbound)
c$omp& private(i,j,xi,yi,zi,dx,dy,dz,rad2)
         do i = 1,n1
            if (iphase(i)==0) then
               xi = xyzmh(1,i)
               yi = xyzmh(2,i)
               zi = xyzmh(3,i)
               do j = n1+1,npart
                  if (iphase(j)==0) then
                     dx = xi - xyzmh(1,j)
                     dy = yi - xyzmh(2,j)
                     dz = zi - xyzmh(3,j)
                     if (abs(dx)>hdbx) dx = dx-dxbound*SIGN(1.0,dx)
                     if (abs(dy)>hdby) dy = dy-dybound*SIGN(1.0,dy)
                     if (abs(dz)>hdbz) dz = dz-dzbound*SIGN(1.0,dz)
                     rad2 = dx*dx + dy*dy + dz*dz
                     if (rad2 < rmin(i)) rmin(i) = rad2
                  endif
               enddo
            endif
         enddo
c$omp end parallel do
         rmin = sqrt(rmin)
         ! print results to file
         write(debugfile,'(a,I2.2,a)') 'fort',isink,'.4233'
         open(unit=4233,file=trim(debugfile))
         do i = 1,n1
            Bi = sqrt(dot_product(Bxyz(1:3,i),Bxyz(1:3,i)))
            Bave = sqrt(Bxave**2 + Byave**2 + Bzave**2)
            write(4233,'(12Es18.6)') rmin(i),Bxyz(1:3,i),
     &                               Bi,Bxave,Byave,Bzave,Bave,
     &                               rho(i),rhobkg,hbkg
         enddo
         close(4233)
      endif


c-----------------------------------------------------------------------
c--   Write dump file & new ifile
c-----------------------------------------------------------------------
c     Write new dump
      nptmass = 0
      print*, 'Writing to output file < ',trim(fileout),' >'
      file1 = trim(fileout)
      if (iproc.eq.0) open(11,file=file1,form='unformatted')
      print *,iproc,': Calling wdump ',npart
      call wdump(11)
c     close(11) ! included in wdump
      print *,iproc,': Done wdump'

c     Write new ifile
      open(unit=12,file='ifile_'//prefix)
      open(unit=13,file='ifile_'//outfix)
      read(12,*) cdummy
      read(12,*) cdummy
      write(13,'(2a)') outfix,'00'
      write(13,'(2a)') outfix,'0000'
      io = 0
      readline = .true.
      do while (io==0)
         read(12,'(a)',iostat=io) cdummy
         if (io==0 .and. readline) then
            if (old_rad .and. index(cdummy,'boundtemp') > 0) then
               write(13,'(4a)') 
     &         '   1.00000E-05  5.00000E+00  5.00000E-04  1.00000E+00',    ! Note that we changed boundarytemp from 13.8K to 5K
     &         ' 10  1.00000E+00 Implicit radiative transfer ',
     &         'tolerance, boundtemp, bounddens, opacity divisor, ',
     &         'opacity model, metallicity'
               readline = .false.
            elseif (index(cdummy,'Max Dynamic Time') > 0) then
               write(13,'(Es11.4,a)') timeend, '  Max Dynamic Time'
            elseif (index(cdummy,' Max Timestep') > 0 .and.
     &              index(cdummy,' Dump') == 0) then
c               write(13,'(a)') '   2.5205625-9  Max Timestep'
               write(13,'(a)') '   0.00025205625  Max Timestep'
c               write(13,'(a)') '   4.0329E-03  Max Timestep'
            elseif (index(cdummy,' Initial Timestep') > 0) then
               write(13,'(a)') '   4.03294E-10  Initial Timestep'
            elseif (index(cdummy,' Binary Dump Counter') > 0) then
               write(13,'(a)') '   1 Binary Dump Counter'
            elseif (index(cdummy,'Critical Density') > 0) then
               write(13,'(a)') trim(cdummy)
               write(13,'(8Es13.5)') -1.,0.,
     &                               xmin,xmax,ymin,ymax,zmin,zmax
               io = 1
            else
               write(13,'(a)') trim(cdummy)
            endif
         else
            readline = .true.
         endif
      enddo
      close(12)
      close(13)
      END

c ----------------------------------------------------------------------
      function wkern(q)
      real    :: wkern,q
      real, parameter   :: fourpi = 12.5663706144d0

      if (q.lt.0.0) then
         write(*,*) 'well, shit'
      else if (q.lt.1.0) then
         wkern = 4.0 - 6.0*q**2 + 3*q**3
      else if (q.lt.2.0) then
         wkern = (2.0 - q)**3
      else
         wkern = 0.0
      endif
      wkern = wkern / fourpi

      end function wkern

c ----------------------------------------------------------------------
      function costheta(vec1,vec2)
      real :: costheta,vec1(3),vec2(3)
      real :: vec12,vec11r,vec22r

      vec12  =      dot_product(vec1,vec2)
      vec11r = sqrt(dot_product(vec1,vec1))
      vec22r = sqrt(dot_product(vec2,vec2))
      if (vec11r > 0. .and. vec22r > 0.) then
         costheta = vec12/(vec11r*vec22r)
      else
         costheta = 0.
      endif

      end function costheta

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The next set of subroutine are required for ellipse fitting.
! First, we define a grid, then project the particles onto the grid
! Next, we use the grid to determine the ellipses, and rotate until the
! ellipses lies along the Cartesian axes.  This is necessary, since using
! the SPH particles directly is effectively a density-weighted fit and
! represents the core rather than the full cloud
!----------------------------------------------------------------------
      subroutine fit_ellipse(n,i1,i2,xyz,udistpc,axes,phi)
      integer, parameter   :: ngrid = 1000
      integer              :: n,i1,i2
      integer              :: i,j,p,nx1,nx2,imin,imax,jmin,jmax,hondx
      integer*1            :: igrid(ngrid,ngrid)
      real                 :: dxgrid0_pc,dxgrid0,dx0
      real                 :: h2,r2,x10,x20,phi,xc,yc
      real                 :: m00,m01,m10,m11,m20,m02,a,b,c,x1,x2
      real                 :: xyz(5,n),xmin(3),xmax(3),axes(2)
      real                 :: x1grid(ngrid),x2grid(ngrid)
      real*8               :: udistpc
      !
      !--Set the grid
      !
      ! set the resolution for the grid
      dxgrid0_pc = 0.001
      dxgrid0 = dxgrid0_pc/udistpc

      ! find extreme values
      xmin =  huge(xmin(1))
      xmax = -huge(xmax(1))
c$omp parallel do default(none)
c$omp&shared(n,xyz,i1,i2)
c$omp&private(i)
c$omp&reduction(min:xmin)
c$omp&reduction(max:xmax)
      do i = 1,n
         xmin(1) = min(xmin(1),xyz(i1,i))
         xmax(1) = max(xmax(1),xyz(i1,i))
         xmin(2) = min(xmin(2),xyz(i2,i))
         xmax(2) = max(xmax(2),xyz(i2,i))
         xmin(3) = min(xmin(3),xyz(5, i))
         xmax(3) = max(xmax(3),xyz(5, i))
      enddo
c$omp end parallel do

      ! Verify the grid size is large enough to encompass the cloud; otherwise reset dx
      nx1 = int((xmax(1) - xmin(1) + 4.*xmax(3))/dxgrid0) + 2
      nx2 = int((xmax(2) - xmin(2) + 4.*xmax(3))/dxgrid0) + 2
      if (nx1 > ngrid .or. nx2 > ngrid) then
         if (nx1 > nx2) then
            dx0 = (xmax(1) - xmin(1) + 4.*xmax(3))/(ngrid*0.9)
         else
            dx0 = (xmax(2) - xmin(2) + 4.*xmax(3))/(ngrid*0.9)
         endif
         print*, 'Cloud too big for grid.  Resetting dx = ',
     &           dx0, nx1,nx2,ngrid
      else
         dx0 = dxgrid0
      endif

      ! define the grid
      x1grid(ngrid/2) = 0.5*(xmax(1) + xmin(1))
      x2grid(ngrid/2) = 0.5*(xmax(2) + xmin(2))
      do i = ngrid/2+1,ngrid
         x1grid(i) = x1grid(i-1) + dx0
         x2grid(i) = x2grid(i-1) + dx0
      enddo
      do i = ngrid/2-1,1,-1
         x1grid(i) = x1grid(i+1) - dx0
         x2grid(i) = x2grid(i+1) - dx0
      enddo
      !
      !--Create a binary grid where the entry is one if a particle overlaps a cell
      !
      igrid = 0
c$omp parallel do default(none)
c$omp&shared(n,xyz,i1,i2,x1grid,x2grid,igrid,dx0)
c$omp&private(p,i,j,hondx,imin,imax,jmin,jmax,x10,x20,h2,r2)
      do p = 1,n
         i = 1
         do while (xyz(i1,p) > x1grid(i) .and. i < ngrid)
            i = i + 1
         enddo
         j = 1
         do while (xyz(i2,p) > x2grid(j) .and. j < ngrid)
            j = j + 1
         enddo
         if (i==1 .or. i==ngrid) print*, 
     &   'WARNING! SPH PARICLE CENTRE IS OUTSIDE OF THE GRID', 
     &   i,x1grid(1),xyz(i1,p),x1grid(ngrid)
         if (j==1 .or. j==ngrid) print*, 
     &   'WARNING! SPH PARICLE CENTRE IS OUTSIDE OF THE GRID', 
     &   j,x2grid(1),xyz(i2,p),x2grid(ngrid)

         igrid(i,j) = 1
         ! tag the surrounding cells to account for smoothing length
         hondx = int(xyz(5,p)/dx0)
         imin = max(1,    i-hondx - 1)
         imax = min(ngrid,i+hondx + 1)
         jmin = max(1,    j-hondx - 1)
         jmax = min(ngrid,j+hondx + 1)
         x10  = x1grid(i)
         x20  = x2grid(j)
         h2   = xyz(5,p)*xyz(5,p)
         do i = imin,imax
            do j = jmin,jmax
               r2 = (x1grid(i)-x10)**2 + (x2grid(j)-x20)**2
               if (h2 > r2) igrid(i,j) = 1
            enddo
         enddo
      enddo
c$omp end parallel do
      !
      !-- 2D ellipse fitting routine of Rocha,Velho & Carvalho (2002)
      !
      m00 = 0.
      m01 = 0.
      m10 = 0.
      m11 = 0.
      m02 = 0.
      m20 = 0.
c$omp parallel do default(none)
c$omp&shared(x1grid,x2grid,igrid)
c$omp&private(i,j)
c$omp&reduction(+:m00,m10,m01,m11,m20,m02)
      do i = 1,ngrid
         do j = 1,ngrid
            m00 = m00 +                     igrid(i,j)
            m10 = m10 + x1grid(i)          *igrid(i,j)
            m01 = m01 + x2grid(j)          *igrid(i,j)
            m11 = m11 + x1grid(i)*x2grid(j)*igrid(i,j)
            m20 = m20 + x1grid(i)**2       *igrid(i,j)
            m02 = m02 + x2grid(j)**2       *igrid(i,j)
         enddo
      enddo
c$omp end parallel do
 
      ! the centre of the ellipse
      xc = m10/m00
      yc = m01/m00
      ! intermediate calculations
      a = m20/m00 - xc**2
      b = 2.0*(m11/m00 - xc*yc)
      c = m02/m00 - yc**2
      ! the results: rotation angle, minor axis & major axis
      phi   = 0.5*atan2(b,a-c)
      axes(1) = sqrt(6.*(a+c + sqrt(b*b+(a-c)**2)))
      axes(2) = sqrt(6.*(a+c - sqrt(b*b+(a-c)**2)))
      !
      !-- Rotate the data to the Cartesian axes
      !
c$omp parallel do default(none)
c$omp&shared(n,xyz,i1,i2,xc,yc,phi)
c$omp&private(i,x1,x2)
c$omp&reduction(min:xmin)
c$omp&reduction(max:xmax)
      do i = 1,n
         x1 = xyz(i1,i) - xc
         x2 = xyz(i2,i) - yc
         xyz(i1,i) =  x1*cos(phi) + x2*sin(phi) + xc
         xyz(i2,i) = -x1*sin(phi) + x2*cos(phi) + yc
      enddo
c$omp end parallel do

      end subroutine fit_ellipse
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c     The next set of routines are to characterise data.  The first
c     is a generalised subroutine to calculate the spred of data,
c     while the latter will calculate histograms to plot the
c     distribution
c ----------------------------------------------------------------------
c     Get the range encompassing selected percentages:
      subroutine get_spread(npart,list,val,valrange,aval,itype)
      integer :: npart,list(npart),indx(npart),itype
      real    :: val(npart),valrange(7),aval
      real    :: std(3)

      std(1) = 0.341  ! all particles within one standard deviation    (skewed)
      std(2) = 0.477  ! all particles within two standard deviations   (skewed)
      std(3) = 0.498  ! all particles within three standard deviations (skewed)
      
      CALL indexx(npart, list, val, indx)
      valrange(1) = val(list(indx(int(npart*(0.5-std(3))))))
      valrange(2) = val(list(indx(int(npart*(0.5-std(2))))))
      valrange(3) = val(list(indx(int(npart*(0.5-std(1))))))
      valrange(4) = aval/npart
      if (itype.eq.0) valrange(4) = 10**valrange(4)
      valrange(5) = val(list(indx(int(npart*(0.5+std(1))))))
      valrange(6) = val(list(indx(int(npart*(0.5+std(2))))))
      valrange(7) = val(list(indx(int(npart*(0.5+std(3))))))
c     print*, 'value bracketing 99.6%:',valrange(1),valrange(9)
c     print*, 'value bracketing 95.4%:',valrange(2),valrange(8)
c     print*, 'value bracketing 68.2%:',valrange(4),valrange(6)

      end subroutine get_spread

c ----------------------------------------------------------------------
c     Create histograms of density, magnetic field strength & angle
      subroutine get_histo(nbins,npart,pmass,rho,Bi,theta,Temp)
      integer, parameter :: nproperty = 4
      integer      :: nbins,npart,i,j,ibins(nproperty,nbins)
      real         :: rho(npart),Bi(npart),theta(npart),Temp(npart)
      real         :: Dbins(nbins+1),Bbins(nbins+1)
      real         :: Tbins(nbins+1),Hbins(nbins+1)
      real         :: rbins(nproperty,nbins),area(nproperty)
      real         :: pmass
      real         :: Dmin,Dmax,Bmin,Bmax,Tmin,Tmax,Hmin,Hmax
      real         :: dlogD,dlogB,dT,dH
      logical      :: fexists
      character*64 :: filehisto

      print*, 'Making Histograms'
      Dmin = 1.0d-20
      Dmax = 1.0d-16
      Bmin = 1.0d-6
      Bmax = 1.0d-2
      Hmin = 0.
      Hmax = 150.
      Tmin = -1.
      Tmax =  1.
c     make bins
      dlogD = (log10(Dmax) - log10(Dmin))/nbins
      dlogB = (log10(Bmax) - log10(Bmin))/nbins
      dT    = (Tmax - Tmin)/nbins
      dH    = (Hmax - Hmin)/nbins
      write(*,'(a,I,4f18.6)') 'Bins:', nbins,dlogD,dlogB,dT,dH
      do i = 1,nbins+1
         Dbins(i) = 10**(log10(Dmin) +(i-1)*dlogD)
         Bbins(i) = 10**(log10(Bmin) +(i-1)*dlogB)
         Tbins(i) = Tmin + (i-1)*dT
         Hbins(i) = Hmin + (i-1)*dH
c        print*, i,Dbins(i),Bbins(i),Tbins(i)
      enddo

c     fill the bins
      ibins = 0
      do j = 1,npart 
         do i = 1,nbins
            if (Dbins(i) < rho(j) .and. rho(j) < Dbins(i+1) ) then
               ibins(1,i) = ibins(1,i) + 1
            endif
            if (Bbins(i) < Bi(j) .and. Bi(j) < Bbins(i+1) ) then
               ibins(2,i) = ibins(2,i) + 1
            endif
            if (Tbins(i) < theta(j) .and. theta(j) < Tbins(i+1) ) then
               ibins(3,i) = ibins(3,i) + 1
            endif
            if (Hbins(i) < Temp(j) .and. Temp(j) < Hbins(i+1) ) then
               ibins(4,i) = ibins(4,i) + 1
            endif
         enddo
      enddo
c     normalise the bins
      area  = 0.
      rbins = 0.
      do i = 1,nbins
         if (ibins(1,i) > 0) area(1) = area(1) + ibins(1,i)*dlogD
         if (ibins(2,i) > 0) area(2) = area(2) + ibins(2,i)*dlogB
         if (ibins(3,i) > 0) area(3) = area(3) + ibins(3,i)*dT
         if (ibins(4,i) > 0) area(4) = area(4) + ibins(4,i)*dH
      enddo
      do i = 1,nbins
         do j = 1,nproperty
            rbins(j,i) = ibins(j,i)/area(j)
         enddo
      enddo

c     Write to file
      write(filehisto,'(a,I3.3,a)') 'core_Histo',nbins,'.dat'
      inquire(file=filehisto,exist=fexists)
      if (fexists) then
         open(unit=67,file=trim(filehisto),position='append')
      else
         open(unit=67,file=trim(filehisto))
            write(67,"(2a)") '# units are cgs'
            write(67,"('#',13(1x,'[',i2.2,1x,a11,']',2x))")
     &                     1,'rho',
     &                     2,'n rho',
     &                     3,'n/A rho',
     &                     4,'B',
     &                     5,'n B',
     &                     6,'n/A B',
     &                     7,'theta',
     &                     8,'n theta',
     &                     9,'n/A theta',
     &                    10,'Temp',
     &                    11,'n Temp',
     &                    12,'n/A Temp',
     &                    13,'mass'
      endif

      do i = 1,nbins
         write(67,'(4(Es18.6,1x,I18,1x,Es18.6,1x),Es18.6,1x)') 
     &      0.5*(Dbins(i)+Dbins(i+1)),ibins(1,i),rbins(1,i),
     &      0.5*(Bbins(i)+Bbins(i+1)),ibins(2,i),rbins(2,i),
     &      0.5*(Tbins(i)+Tbins(i+1)),ibins(3,i),rbins(3,i),
     &      0.5*(Hbins(i)+Hbins(i+1)),ibins(4,i),rbins(4,i),
     &      npart*pmass
      enddo
      write(67, '(a)') ' '
      close(67)

      end subroutine get_histo
