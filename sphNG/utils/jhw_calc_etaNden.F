      PROGRAM jhw_calc_etaNden
c***********************************************************
c  Calculates the distributions of non-ideal MHD           *
c  coefficients, eta, & ion densities                      *
c  This excludes the background gas from all calculations, *
c  thus is designed for sphere-in-box simulations.         *
c  This necessarily requires being compiled with           *
c  nonideal=yes                                            *
c  in sphNG/, compile as make jhw_calc_etaNden             *
c     author: James Wurster                                *
c***********************************************************

#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
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
      INCLUDE 'COMMONS/nonideal'
      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09rho'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi,udisti,utimei,umagfdi

      PARAMETER (maxfiles = 100)
      INTEGER     :: i,j,p,istart,istop,ic(3)
      INTEGER     :: iarray(maxfiles)
      INTEGER     :: ierrlist(n_warn),list(4,imhd2)
      REAL        :: nden_nimhdi(n_nden),data_out(n_data_out)
      REAL        :: nrat(4,imhd2)
      REAL        :: nave(8),ranges(8,7)
      REAL        :: ekclei,Bi,temperature
      REAL        :: n_total_gas,n_ion_gas,n_electron
      REAL        :: n_grain_pos,n_grain_neg,n_grain_neu,n_grain_tot
      REAL*8      :: unit_density,unit_eta
      CHARACTER*7 :: dump
      CHARACTER*3 :: prefix,outfix
      LOGICAL     :: dumprange,kexists,write_header

c     Get inputs
      iprint = istdout 
      iread  = 5
      call get_in_give_out(prefix,outfix,.false.,dumprange,
     &                           maxfiles,istart,istop,iarray)

      WRITE(inname,'(2a)')'ifile_',prefix
      CALL options

c     Other parameters
      varmhd   = 'Brho'
      encal    = 'r'
      maxrec   = 100*idim
      write_header = .true.

      do i = 1, idim
         isort(i) = i
         iorig(i) = i
      enddo

c     Process data files
      do kloop=istart,istop
         if (dumprange) then
            iloop = kloop
         else
            iloop = iarray(kloop)
         endif
         write(dump,'(a,I4.4)') prefix,iloop
         inquire(file=dump,exist=kexists)
         if (kexists) then
            print*, '!************************'
            print*, '!  analysing: ',dump
            print*, '!************************'
            open (unit=11, file=dump, form='unformatted', recl=maxrec)
            call rdump(11,ichkl,1)
            close (11)

            if (kloop.eq.istart) then
               umass  = umassi
               udist  = udisti
               utime  = utimei
               umagfd = umagfdi
               uergg  = uerggi
               unit_density = umass/udist**3
               unit_eta = udist**2/utime
               nicil_passdust = .FALSE.
               if (idimHY09.EQ.idim) then
               call nicil_initialise_dust(utime,udist,umass,umagfd,
     &                              HY09binsizes,HY09_ndust_bins,
     &                              nicil_passdust,ierr)

               else
                  call nicil_initialise(utime,udist,umass,umagfd,ierr)
               endif
               if (ierr/=0) call quit(0)
               call nimhd_init_gastemp(Rg,gamma,gam,gamdh,gamah,uergg
     &                                ,rhocrit,rhocrit2,rhocrit3)
            endif

c           Calculate ranges of non-ideal coefficients & number densities
            print*, "Calculating range of number densities"
            nrat = 1.e12
            nave = 0.
            ic   = 0
            list = 0
            itry = 0
            rho_max = 0.
c$omp  parallel default(none)
c$omp& shared(npart,encal,n1,n2,nicil_passdust,iden,gamma)
c$omp& shared(iunique,iphase,ekcle,vxyzu,Bxyz,rho)
c$omp& shared(eta_nimhd,HY09bin_rho,nrat,list)
c$omp& private(i,j,p,id2,ekclei,Bi,temperature,n_total_gas,itry)
c$omp& private(n_electron,n_grain_pos,n_grain_neu,n_grain_neg)
c$omp& private(n_grain_tot,n_ion_gas,ierrlist,nden_nimhdi,data_out)
c$omp& reduction(max:rho_max)
c$omp& reduction(+:nave,ic)
c$omp  do
            do i = 1,npart
               list(1,i) = i
               if (iphase(i)==0 .and. iunique(i).lt.n1) then
                  rho_max = max(rho_max,rho(i))
                  if (encal.EQ.'r') then
                     ekclei = ekcle(3,i)
                  else
                     ekclei = 0.0
                  endif
                  Bi = dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
                  Bi = sqrt(Bi)
                  temperature = nimhd_gastemp(encal,vxyzu(4,i),
     &                                ekclei,dble(rho(i)),gamma,
     &                                n1,n2,iunique(i))
                  ierrlist = 0
                  nden_nimhdi = 0.0
                  if (nicil_passdust) then
                     call nicil_update_nimhd(0,eta_nimhd(1,i),
     &                                 eta_nimhd(2,i),eta_nimhd(3,i),
     &                                 Bi,real(rho(i)),temperature,
     &                                 nden_nimhdi,ierrlist,data_out,
     &                                 itry,real(HY09bin_rho(:,i)))
                  else
                     call nicil_update_nimhd(0,eta_nimhd(1,i),
     &                                 eta_nimhd(2,i),eta_nimhd(3,i),
     &                                 Bi,real(rho(i)),temperature,
     &                                 nden_nimhdi,ierrlist,data_out)
                  endif
                  n_total_gas = data_out(5)
                  n_electron  = data_out(8)
                  n_ion_gas   = 0.
                  do j = 9,21 
                     n_ion_gas = n_ion_gas + data_out(j)
                  enddo
                  n_grain_pos = 0. 
                  n_grain_neu = 0. 
                  n_grain_neg = 0. 
                  do j = 1,na 
                     id2 = 22+3*(j-1)
                     n_grain_pos = n_grain_pos + data_out(id2  )
                     n_grain_neu = n_grain_neu + data_out(id2+1)
                     n_grain_neg = n_grain_neg + data_out(id2+2)
                  enddo
                  n_grain_tot = n_grain_pos +n_grain_neu +n_grain_neg

                  if (n_total_gas+n_grain_tot > 0. .and.
     &                n_electron > 0) then
                     nrat(1,i) = n_electron/(n_total_gas+n_grain_tot) ! electrons / total
                     nave(1)   = nave(1) + log10(nrat(1,i))
                  endif
                  if (n_total_gas > 0. .and. n_ion_gas > 0.) then
                     nrat(2,i) = n_ion_gas/n_total_gas                ! gas ion / gas total
                     nave(2)   = nave(2) + log10(nrat(2,i))
                  endif
                  if (n_grain_tot > 0. .and. n_grain_pos > 0) then
                     nrat(3,i) = n_grain_pos/n_grain_tot              ! positive dust / total dust
                     nave(3)   = nave(3) + log10(nrat(3,i))
                  endif
                  if (n_grain_tot > 0. .and. n_grain_neg > 0) then
                     nrat(4,i) = n_grain_neg/n_grain_tot              ! negative dust / total dust
                     nave(4)   = nave(4) + log10(nrat(4,i))
                  endif

               endif
            enddo
c$omp  enddo
c$omp  end parallel

            do j = 1,4
               call get_spread(npart,list(1,1:npart),nrat(j,1:npart),
     &                 ranges(j,1:7),nave(j),.true.,)
            enddo

            print*, "Calculating range of eta"
            ic = 0
            do i = 1,npart
               if (iphase(i)==0 .and. iunique(i).lt.n1) then
                  ic(1) = ic(1)+1
                  list(2,ic(1)) = i
                  nave(5) = nave(5) + log10(eta_nimhd(1,i))
                  nave(6) = nave(6) + log10(eta_nimhd(3,i))
                  if (eta_nimhd(2,i) > 0.) then
                     ic(2) = ic(2)+1
                     list(3,ic(2)) = i
                     nave(7) = nave(7) + log10(eta_nimhd(2,i))
                  elseif (eta_nimhd(2,i) < 0.) then
                     ic(3) = ic(3)+1
                     list(4,ic(3)) = i
                     nave(8) = nave(8) + log10(-eta_nimhd(2,i))
                  endif
               endif
            enddo

           call get_spread(ic(1),list(2,1:ic(1)),eta_nimhd(1,1:npart),
     &                 ranges(5,1:7),nave(5),.false.)
           call get_spread(ic(1),list(2,1:ic(1)),eta_nimhd(3,1:npart),
     &                 ranges(6,1:7),nave(6),.false.)
           call get_spread(ic(2),list(3,1:ic(2)),eta_nimhd(2,1:npart),
     &                 ranges(7,1:7),nave(7),.false.)
           call get_spread(ic(3),list(4,1:ic(3)),eta_nimhd(2,1:npart),
     &                 ranges(8,1:7),nave(8),.false.)
            ranges(8,4) = -ranges(8,4)

            if (write_header) then
               write_header = .false.
               open(unit=67,file='CoreIonisation.dat')
               write(67,"('#',59(1x,'[',i2.2,1x,a11,']',2x))")
     &                     1,'idump',
     &                     2,'time',
     &                     3,'rho_max',
     &                     4,'ne/nT 0.2%',
     &                     5,'ne/nT 2.3%',
     &                     6,'ne/nT 15.9% ',
     &                     7,'ne/nT logave',
     &                     8,'ne/nT 84.1%',
     &                     9,'ne/nT 97.7%',
     &                    10,'ne/nT 99.8%',
     &                    11,'ni/n 0.2%',
     &                    12,'ni/n 2.3%',
     &                    13,'ni/n 15.9% ',
     &                    14,'ni/n logave',
     &                    15,'ni/n 84.1%',
     &                    16,'ni/n 97.7%',
     &                    17,'ni/n 99.8%',
     &                    18,'ngp/nG 0.2%',
     &                    19,'ngp/nG 2.3%',
     &                    20,'ngp/nG 15.9% ',
     &                    21,'ngp/nG logave',
     &                    22,'ngp/nG 84.1%',
     &                    23,'ngp/nG 97.7%',
     &                    24,'ngp/nG 99.8%',
     &                    25,'ngn/nG 0.2%',
     &                    26,'ngn/nG 2.3%',
     &                    27,'ngn/nG 15.9% ',
     &                    28,'ngn/nG logave',
     &                    29,'ngn/nG 84.1%',
     &                    30,'ngn/nG 97.7%',
     &                    31,'ngn/nG 99.8%',
     &                    32,'eta_o 0.2%',
     &                    33,'eta_o 2.3%',
     &                    34,'eta_o 15.9% ',
     &                    35,'eta_o logave',
     &                    36,'eta_o 84.1%',
     &                    37,'eta_o 97.7%',
     &                    38,'eta_o 99.8%',
     &                    39,'eta_a 0.2%',
     &                    40,'eta_a 2.3%',
     &                    41,'eta_a 15.9% ',
     &                    42,'eta_a logave',
     &                    43,'eta_a 84.1%',
     &                    44,'eta_a 97.7%',
     &                    45,'eta_a 99.8%',
     &                    46,'eta_hp 0.2%',
     &                    47,'eta_hp 2.3%',
     &                    48,'eta_hp 15.9% ',
     &                    48,'eta_hp logave',
     &                    50,'eta_hp 84.1%',
     &                    51,'eta_hp 97.7%',
     &                    52,'eta_hp 99.8%',
     &                    53,'eta_hn 0.2%',
     &                    54,'eta_hn 2.3%',
     &                    55,'eta_hn 15.9% ',
     &                    56,'eta_hn logave',
     &                    57,'eta_hn 84.1%',
     &                    58,'eta_hn 97.7%',
     &                    59,'eta_hn 99.8%'


            endif
            write(67,'((I18,1x),58(Es18.6,1x))') iloop,gt,
     &      rho_max*unit_density,
     &      ranges(1,1:7),ranges(2,1:7),ranges(3,1:7),ranges(4,1:7),
     &      ranges(5,1:7)*unit_eta,ranges(6,1:7)*unit_eta,
     &      ranges(7,1:7)*unit_eta,ranges(8,1:7)*unit_eta

         endif ! end if kexists
      enddo  ! end kloop
      close(67)

#else
      print*, 'This must be compiled with nonideal=yes'
#endif

      END

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c     The next set of routines are to characterise data.  The first
c     is a generalised subroutine to calculate the spred of data,
c     while the latter will calculate histograms to plot the
c     distribution
c ----------------------------------------------------------------------
c     Get the range encompassing selected percentages:
      subroutine get_spread(npart,list,val,valrange,aval,findn)
      integer :: npart,list(npart),indx(npart),i,j,k
      real    :: val(npart),valrange(7),aval
      real    :: std(3)
      logical :: findn,keep_searching

      std(1) = 0.341  ! all particles within one standard deviation    (skewed)
      std(2) = 0.477  ! all particles within two standard deviations   (skewed)
      std(3) = 0.498  ! all particles within three standard deviations (skewed)
      
      CALL indexx(npart, list, val, indx)

      if (findn) then
         i = 1
         keep_searching = .true.
         do while (i < npart .and. keep_searching)
            k = list(indx(i))
            if (val(k) < 1.e10) then
               i = i + 1
            else
               keep_searching = .false.
            endif
         enddo
      else
         i = npart
      endif

      if (i==0) then
         valrange = 0
      else
         valrange(1) = val(list(indx(int(i*(0.5-std(3))))))
         valrange(2) = val(list(indx(int(i*(0.5-std(2))))))
         valrange(3) = val(list(indx(int(i*(0.5-std(1))))))
         valrange(4) = 10**(aval/i)
         valrange(5) = val(list(indx(int(i*(0.5+std(1))))))
         valrange(6) = val(list(indx(int(i*(0.5+std(2))))))
         valrange(7) = val(list(indx(int(i*(0.5+std(3))))))
      endif

      end subroutine get_spread
