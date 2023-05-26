      PROGRAM jhw_calc_Mach
c***********************************************************
c                                                          *
c  This program reads sphNG dump files and calculates the  *
c  Mach number of the particles in object 1                *
c     author: James Wurster                                *
c  in sphNG/, compile as make jhw_calc_etaNden             *
c                                                          *
c***********************************************************
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/bodys'
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      INTEGER, PARAMETER :: maxfiles = 512
      INTEGER, PARAMETER :: nbin     = 256
      INTEGER     :: k,i,j,p,pm1,kstart,kstop,maxrec,iMach
      INTEGER     :: iarray(maxfiles)
      LOGICAL     :: kexists,dumprange,firstcall
      real        :: Mach(2),vrms,Tmin,Tave,Tmax
      real*8      :: uerggs
      character* 3 :: prefix,outfix
      character* 7 :: dump
      character*14 :: fileout
c
c--read options
c
      call get_in_give_out(prefix,outfix,.false.,dumprange,
     &                           maxfiles,kstart,kstop,iarray)
c
c--initialise parameters
c
      varmhd    = 'Brho'
      encal     = 'r'
      maxrec    = 100*idim
      firstcall = .true.
c
c--Process data files
c
      do k=kstart,kstop
         IF (dumprange) THEN
            iloop = k
         ELSE
            iloop = iarray(k)
         ENDIF
         write(dump,'(a,I4.4)') prefix,iloop
         inquire(file=dump,exist=kexists)
         if (kexists) then
            print*, '**********************'
            print*, 'analysing file ',dump
            print*, '**********************'
            open (unit=11, file=dump, form='unformatted', recl=maxrec)
            call rdump(11, ichkl, 0)
            close(11)

            if (firstcall) then
               uergg = udisti**2/utimei**2
               print*, uergg,Rg,gmw,encal
            endif
            Mach  = 0.
            iMach = 0
            vrms  = 0.
            Tmin  = huge(Tmin)
            Tave  = 0.
            Tmax  = 0.
            print*, 'Finding Mach number '
c$omp  parallel do default(none)
c$omp& shared(npart,vxyzu,ekcle,iphase,iunique,n1,uergg,encal)
c$omp& private(i,v2i,cs2i,temp)
c$omp& reduction(+:Mach,iMach,vrms,Tave)
c$omp& reduction(min:Tmin)
c$omp& reduction(max:Tmax)
            do i = 1,npart
               if (iphase(i)==0 .and. iunique(i).lt.n1) then
                  v2i = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
                  if (.false.) then
                     if (encal.EQ.'r') then
                        temp = vxyzu(4,i)/ekcle(3,i)
                        cs2i = (5.0*Rg*temp/(3.0*gmw*uergg))
                     else
                        cs2i = (2./3.*vxyzu(4,i))
                     endif
                     Mach  = Mach  + v2i/cs2i
                  else
                     temp    = vxyzu(4,i)/ekcle(3,i)
                     cs2i    = (5.0*Rg*temp/(3.0*gmw*uergg))
                     Mach(1) = Mach(1) + v2i/cs2i
                     cs2i    = (2./3.*vxyzu(4,i))
                     Mach(2) = Mach(2) + v2i/cs2i
                  endif
                  iMach = iMach + 1
                  vrms  = vrms  + v2i
                  Tmin  = min(Tmin,temp)
                  Tmax  = max(Tmax,temp)
                  Tave  =     Tave+temp
               endif
            enddo
c$omp end parallel do
            Mach = sqrt(Mach/iMach)
            vrms = sqrt(vrms/iMach)
            Tave = Tave/iMach
            print*, 'encalr: The Mach number of ',dump,' is ',Mach(1)
            print*, 'The Mach number of ',dump,' is ',Mach(2)
            if (firstcall) then
               firstcall = .false.
               write(fileout,'(2a)')prefix,'_MachNo.dat'
               open(unit=17,file=fileout)
               write(17,"('#',9(1x,'[',i2.2,1x,a11,']',2x))")
     &           1,'idump',
     &           2,'time',
     &           3,'Mach:encalr',
     &           4,'Mach:ideal',
     &           5,'vrms',
     &           6,'Tmin',
     &           7,'Tave',
     &           8,'Tmax',
     &           9,'npart'
            endif
            write(17,'(I18,1x,7(es18.10,1x),I18)')
     &      iloop,gt,Mach,vrms,Tmin,Tave,Tmax,iMach
         endif
      enddo
      close(17)

      stop
      end

