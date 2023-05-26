      PROGRAM jhw_print_sinkorigin
c***********************************************************
c                                                          *
c  Sink order is not preserved, thus we cannot simply      *
c  select a sink in splash for the origin.  This programme *
c  will go through a series of dumps and compile           *
c  1. a list of positions of the first (i.e. oldest) sink  *
c  2. the xyz cooridinage of the first (i.e. oldest) sink  *
c  Splash has been update to read these files (use only    *
c  one) so that there will be a consistent                 *
c  origin/rest frame when making a movie                   *
c  The input is a splash.filenames file                    *
c     author: James Wurster                                *
c  in sphNG/, compile as make jhw_moddump_CoM              *
c                                                          *
c***********************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/varmhd'

      integer, parameter :: maxmodels = 18
      integer      :: i,io,irw,nmodels
      integer      :: idsink(maxmodels)
      logical      :: kexists,skipzero
      real         :: xyz(3)
      character* 3 :: prefix,outfix
      character*64 :: dump
      character*32 :: file_splashdumps,file_splashids,file_splashxyz
      character*32 :: file_sinklog
c
c--read options
      iread = 5
      write (*,*) 'Enter splash.filename: '
      read (iread, *) file_splashdumps
      write (*,*) 'The number of models per image: '
      read (iread, *) nmodels
      write (*,*) 'Skip dump 0: '
      read (iread, *) skipzero
      write (*,*) 'How would you like to interact with SinkOrderID.dat'
      write (*,*) '1: read only; 2:read+write; 3:write only'
      read (iread, *) irw
      if (nmodels > maxmodels) then
         print*, 'recompile with larger maxmodels.  aborting'
         stop
      endif
c
c--initialise parameters
c
      maxrec = 100*idim
      varmhd = 'Brho'
      write(file_splashids,'(2a)') 
     &      file_splashdumps(:index(file_splashdumps,'.')),'sinkpos'
      write(file_splashxyz,'(2a)') 
     &      file_splashdumps(:index(file_splashdumps,'.')),'origin'
      write(file_sinklog,'(a)') 'SinkOrderID.dat'
      idsink = 0
      if (irw==1 .or. irw==2) then
         ! sink file exists; read data
         open(unit=17,file=file_sinklog)
         do i = 1,nmodels
            read(17,*) idsink(i)
         enddo
         close(17)
      endif
c
c--Process data files
c
      inquire(file=file_splashdumps,exist=kexists)
      if (.not.kexists) then
         print*, 'Splash file does not exist.  aborting'
         stop
      endif

      io = 0
      k  = 0
      open(unit=17,file=file_splashdumps)
      open(unit=18,file=file_splashids)
      open(unit=19,file=file_splashxyz)
      do while (io==0)
         read(17,'(a)',iostat=io) dump
         isinkpos = 0
         xyz      = 0. 
         k        = k + 1
         if (k > nmodels) k = 1
         if (io==0) then
            inquire(file=trim(dump),exist=kexists)
            if (kexists .and. skipzero) then
               if (index(dump,'0000') > 0) kexists = .false.
            endif
         else
            kexists = .false.
         endif
         if (kexists) then
            print*, 'analysing file ',dump
            print*, 'searching for sink ',k,idsink(k)
            open (unit=11, file=dump, form='unformatted', recl=maxrec)
            call rdump(11, ichkl, 0)
            print*, 'nptmass = ', nptmass
            print*, 'listpm = ',listpm(1:nptmass)
            if (nptmass > 1 .or. irw==1) then
               ! find the location of the sink we care about
               isink = 0
               do i = 1,nptmass
                  if (idsink(k) == iunique(listpm(i))) isink = i
               enddo
               if (isink==0) then
                  print*, 'Aww crap: no match',idsink(k)
                  isinkpos = 0
                  xyz      = 0.
               else
                  print*, 'The sink we want is in position ',isink
                  ! determine the number of sinks with lower particle numbers
                  isinkpos = 1
                  do i = 1,nptmass
                     if (listpm(i) < listpm(isink)) isinkpos=isinkpos+1
                  enddo
                  xyz = xyzmh(1:3,listpm(isink))
               endif
             elseif (nptmass == 1) then
               idsink(k) = iunique(listpm(1))
               isinkpos  = 1
               xyz       = xyzmh(1:3,listpm(1))
            endif
         endif
         if (io==0) then
            write(18,*) isinkpos
            write(19,*) xyz
         endif
      enddo
      close(17)
      close(18)
      close(19)

c     Write sink ids to file
      open(unit=17,file=file_sinklog)
      do i = 1,nmodels
         write(17,*) idsink(i)
      enddo
      close(17)

      print*, 'analysis complete'

      stop
      end

c----------------------------------------------------------------------!
c     Miscellaneous subroutines originally at the end of profile.F
      SUBROUTINE quit
      INCLUDE 'idim'

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

