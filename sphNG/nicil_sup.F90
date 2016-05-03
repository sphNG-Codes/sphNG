!+
!  MODULE: nicil_sup
!          
!  DESCRIPTION:
!  Contains supplementary routines to NICIL for use in sphNG; adapted 
!  from Phantom's evwrite & energies.
!  This is not part of the NICIL package.
!
!  AUTHOR: James Wurster
!
!----------------------------------------------------------------------!
module nicil_sup
 use nicil, only: nicil_get_eta
 use nicil, only: use_ohm,use_hall,use_ambi,ion_rays,ion_thermal, &
                  nelements,nelements_max,nlevels,meanmolmass

 implicit none
 !
 !--Indicies to determine which action to take
 integer(kind=1),   parameter :: iev0 = 0
 integer(kind=1),   parameter :: ievN = 1
 integer(kind=1),   parameter :: ievA = 2  
 integer(kind=1),   parameter :: ievX = 3
 integer,           parameter :: ievfile = 1701
 logical,           private   :: firstcall = .true.
 !
 integer,           private   :: itimeni,itX,itA,itN,ietaFX,ietaFA &
                                ,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN &
                                ,ihallX,ihallA,ihallN,iahallX,iahallA,iahallN &
                                ,ihallfX,ihallfA,ihallfN &
                                ,iahallfX,iahallfA,iahallfN &
                                ,iambiX,iambiA,iambiN,iambifX,iambifA,iambifN &
                                ,inenX,inenA,ineX,ineA,innX,innA &
                                ,izgrainX,izgrainA,izgrainN &
                                ,inihrX,inihrA,inimrX,inimrA,ingX,ingA &
                                ,inistX,inistA,inidtX,inidtA &
                                ,inhX,inhA,inheX,inheA,innaX,innaA,inmgX,inmgA &
                                ,inkX,inkA,inhedX,inhedA,innadX,innadA &
                                ,inmgdX,inmgdA,inkdX,inkdA
  integer,          private   :: ielements
  integer,          parameter :: inumev=100
  real,             private   :: rhocrit_nimhd,rhocrit2_nimhd,rhocrit3_nimhd
  real,             private   :: coef_gammaeq1,coef_gamma,coef_gam,coef_gamdh,coef_gamah
  character(len=19),private   :: ev_label(inumev)
  integer(kind=1),  public    :: ev_action(inumev)
  !--Subroutines
  public                      :: nimhd_init_gastemp,nimhd_gastemp,inform_nimhd
  private                     :: get_coef_gamma,fill_one_label,fill_xan_label

 private
!
contains
!
!-----------------------------------------------------------------------
!+
!  calculate the gas temperature coefficients & locally store critical 
!  densities for encal=="v"
!+
!-----------------------------------------------------------------------
subroutine nimhd_init_gastemp(Rg,gamma,gam,gamdh,gamah,uergg,rhocrit,rhocrit2,rhocrit3)
 real,             intent(in) :: gamma,gam,gamdh,gamah,rhocrit,rhocrit2,rhocrit3
 real(kind=8),     intent(in) :: Rg,uergg
 !
 coef_gammaeq1  = get_coef_gamma(5.0/3.0,uergg,Rg)
 coef_gamma     = get_coef_gamma(gamma,  uergg,Rg)
 coef_gam       = get_coef_gamma(gam,    uergg,Rg)
 coef_gamdh     = get_coef_gamma(gamdh,  uergg,Rg)
 coef_gamah     = get_coef_gamma(gamah,  uergg,Rg)
 !
 rhocrit_nimhd  = rhocrit
 rhocrit2_nimhd = rhocrit2
 rhocrit3_nimhd = rhocrit3
 !
end subroutine nimhd_init_gastemp
!
!--Function to get the coefficient for the gas temperature
real function get_coef_gamma(gamma,uergg,Rg)
 real,             intent(in) :: gamma
 real(kind=8),     intent(in) :: Rg,uergg
 !
 if (gamma > 1.0001) then
   get_coef_gamma = (gamma-1.0)*meanmolmass*uergg/Rg
 else
   get_coef_gamma = 2.0/3.0    *meanmolmass*uergg/Rg
 end if
 !
end function get_coef_gamma
!-----------------------------------------------------------------------
!+
!  calculate the gas temperature
!+
!-----------------------------------------------------------------------
real function nimhd_gastemp(encal,vxyzu4i,ekcle3i,rhoi,gamma)
 real,             intent(in) :: vxyzu4i,ekcle3i,rhoi,gamma
 character(len=*), intent(in) :: encal
 !
 if (encal=='r') then
    nimhd_gastemp = vxyzu4i/ekcle3i
! elseif (encal=='m') then
!!   NOTE: encal=="m" does not seem to be a current option
!    nimhd_gastemp = vxyzu4i/getcv(rhoi,vxyzu4i)
 elseif (encal=='v') then
    if (rhoi < rhocrit_nimhd) then
       nimhd_gastemp = coef_gamma*vxyzu4i
    elseif (rhoi < rhocrit2_nimhd) then
       nimhd_gastemp = coef_gam*vxyzu4i
    elseif (rhoi < rhocrit3_nimhd) then
       nimhd_gastemp = coef_gamdh*vxyzu4i
    else 
       nimhd_gastemp = coef_gamah*vxyzu4i
    endif
 elseif (gamma==1.0) then
    nimhd_gastemp = coef_gammaeq1*vxyzu4i
 else
    nimhd_gastemp = coef_gamma*vxyzu4i
 endif
 !
end function nimhd_gastemp
!----------------------------------------------------------------
!+
!  Prints useful data to file
!+
!----------------------------------------------------------------
subroutine inform_nimhd(time,npart,xyzmh,vxyzu,ekcle,Bxyz,rho,vsound,&
                        iphase,Zg_on_ag,n_electronT,alphaMM,alphamin,gamma,ifsvi,&
                        nimhdfile,encal)
 real,             intent(in) :: time,gamma
 real,             intent(in) :: xyzmh(:,:),vxyzu(:,:),ekcle(:,:),Bxyz(:,:),rho(:),vsound(:), &
                                 Zg_on_ag(:),n_electronT(:),alphaMM(:,:),alphamin(:)
 integer,          intent(in) :: npart,ifsvi
 integer(kind=1),  intent(in) :: iphase(:)
 character(len=*), intent(in) :: nimhdfile,encal
 integer                      :: i,j,ierr
 integer(kind=8)              :: np
 real                         :: dnptot,B2i,Bi,temperature,vsigi, &
                                 etaart,etaart1,etaohm,etahall,etaambi
 real                         :: data_out(13+nelements_max*nlevels-1)
 real                         :: ev_data(0:inumev),ev_data_thread(0:inumev)
 character(len=27)            :: ev_fmt
 character(len=35)            :: ev_format

 !
 !--Determine what values should be in the header
 !
 if ( firstcall ) then
    firstcall = .false.
    i = 1
    write(ev_fmt,'(a)') "(1x,'[',i2.2,1x,a11,']',2x)"
    call fill_one_label(ev_fmt,'time',        i,itimeni, iev0)
    call fill_xan_label(ev_fmt,'temp',        i,itX,     itA,     itN     )
    call fill_xan_label(ev_fmt,'eta_ar',      i,ietaFX,  ietaFA           )
    if (use_ohm) then
       call fill_xan_label(ev_fmt,'eta_o',    i,iohmX,   iohmA,   iohmN   )
       call fill_xan_label(ev_fmt,'eta_o/art',i,iohmfX,  iohmfA,  iohmfN  )
    endif
    if (use_hall) then
       call fill_xan_label(ev_fmt,'eta_h',    i,ihallX,  ihallA,  ihallN  )
       call fill_xan_label(ev_fmt,'|eta_h|',  i,iahallX, iahallA, iahallN )
       call fill_xan_label(ev_fmt,'eta_h/art',i,ihallfX, ihallfA, ihallfN )
       call fill_xan_label(ev_fmt,'|e_h|/art',i,iahallfX,iahallfA,iahallfN)
    endif
    if (use_ambi) then
      call fill_xan_label(ev_fmt,'eta_a',    i,iambiX,  iambiA,  iambiN  )
      call fill_xan_label(ev_fmt,'eta_a/art',i,iambifX, iambifA, iambifN )
    endif
    call fill_xan_label(ev_fmt,'n_e/n',       i,inenX,   inenA            )
    call fill_xan_label(ev_fmt,'n_e',         i,ineX,    ineA             )
    call fill_xan_label(ev_fmt,'n_n',         i,innX,    innA             )
    if (ion_rays) then
       call fill_xan_label(ev_fmt,'Z_Grain',  i,izgrainX,izgrainA,izgrainN)
       call fill_xan_label(ev_fmt,'n_ihR',    i,inihrX,  inihrA           )
       call fill_xan_label(ev_fmt,'n_imR',    i,inimrX,  inimrA           )
       call fill_xan_label(ev_fmt,'n_gR',     i,ingX,    ingA             )
    endif
    if (ion_thermal) then
       call fill_xan_label(ev_fmt,   'n_isT', i,inistX,  inistA           )
       if (nlevels>=2) then
          call fill_xan_label(ev_fmt,'n_idT', i,inidtX,  inidtA           )
       endif
       if (nelements>=2) then
          call fill_xan_label(ev_fmt,'n_H+',  i,inhX,    inhA             )
          call fill_xan_label(ev_fmt,'n_He+', i,inheX,   inheA            )
       endif
       if (nelements>=5) then
          call fill_xan_label(ev_fmt,'n_Na+', i,innaX,   innaA            )
          call fill_xan_label(ev_fmt,'n_Mg+', i,inmgX,   inmgA            )
          call fill_xan_label(ev_fmt,'n_K+',  i,inkX,    inkA             )
       endif
       if (nlevels>=2) then
          if (nelements>=2) then
             call fill_xan_label(ev_fmt,'n_He++',i,inhedX,   inhedA       )
          endif
          if (nelements>=5) then
             call fill_xan_label(ev_fmt,'n_Na++',i,innadX,   innadA       )
             call fill_xan_label(ev_fmt,'n_Mg++',i,inmgdX,   inmgdA       )
             call fill_xan_label(ev_fmt,'n_K++', i,inkdX,    inkdA        )
          endif
       endif
    endif
    ielements = i - 1 ! The number of columns to be calculates
    !
    !--write a header line
    open(ievfile, file=nimhdfile,status='replace')
    write(ev_fmt,'(a,I3,a)') '(',ielements+1,'a)'
    write(ievfile,ev_fmt)'#',ev_label(1:ielements)
 else
    !--open the file
    open(ievfile, file=nimhdfile,position='append')
 end if
 !
 !--Compile non-ideal MHD data
 !
 call initialise_ev_data(ielements,ev_action,ev_data)
 np = 0
!$omp parallel default(none) &
!$omp shared(xyzmh,vxyzu,ekcle,Bxyz,rho,vsound,iphase,npart,encal) &
!$omp shared(Zg_on_ag,n_electronT,alphaMM,alphamin,gamma,ifsvi) &
!$omp shared(use_ohm,use_hall,use_ambi,ion_rays,ion_thermal,nelements) &
!$omp shared(ielements,ev_data,ev_action) &
!$omp shared(itX,itA,itN,ietaFX,ietaFA,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN) &
!$omp shared(ihallX,ihallA,ihallN,iahallX,iahallA,iahallN) &
!$omp shared(ihallfX,ihallfA,ihallfN,iahallfX,iahallfA,iahallfN) &
!$omp shared(iambiX,iambiA,iambiN,iambifX,iambifA,iambifN) &
!$omp shared(inenX,inenA,ineX,ineA,innX,innA,izgrainX,izgrainA,izgrainN) &
!$omp shared(inihrX,inihrA,inimrX,inimrA,ingX,ingA,inistX,inistA,inidtX,inidtA) &
!$omp shared(inhX,inhA,inheX,inheA,innaX,innaA,inmgX,inmgA,inkX,inkA) &
!$omp shared(inhedX,inhedA,innadX,innadA,inmgdX,inmgdA,inkdX,inkdA) &
!$omp private(i,j,ierr,B2i,Bi,temperature,vsigi,etaart,etaart1,etaohm,etahall,etaambi,data_out) &
!$omp private(ev_data_thread)&
!$omp reduction(+:np) 
 call initialise_ev_data(ielements,ev_action,ev_data_thread)
!$omp do
 do i=1,npart
    if (iphase(i)==0) then
       np  = np + 1
       B2i = dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
       Bi  = sqrt(B2i)
       temperature = nimhd_gastemp(encal,vxyzu(4,i),ekcle(3,i),rho(i),gamma)
       call nicil_get_eta(etaohm,etahall,etaambi,Bi,rho(i),temperature, &
                          Zg_on_ag(i),n_electronT(i),ierr,data_out)
       vsigi = sqrt(vsound(i)**2 + B2i/rho(i))
       if (ifsvi >= 6) then
          etaart = alphaMM(2,i)*vsigi*xyzmh(5,i)
       else
          etaart = alphamin(2)*vsigi*xyzmh(5,i)
       endif
       if (etaart > 0.) then
          etaart1 = 1.0/etaart
       else
          etaart1 = 0.0
       endif
       call ev_update(ev_data_thread,temperature,      itX,   itA   ,itN   )
       call ev_update(ev_data_thread,etaart,           ietaFX,ietaFA)
       if (use_ohm) then
          call ev_update(ev_data_thread,etaohm,        iohmX, iohmA, iohmN )
          call ev_update(ev_data_thread,etaohm*etaart1,iohmfX,iohmfA,iohmfN)
       endif
       if (use_hall) then
          call ev_update(ev_data_thread,    etahall,         ihallX,  ihallA,  ihallN  )
          call ev_update(ev_data_thread,abs(etahall),        iahallX, iahallA, iahallN )
          call ev_update(ev_data_thread,    etahall*etaart1, ihallfX, ihallfA, ihallfN )
          call ev_update(ev_data_thread,abs(etahall)*etaart1,iahallfX,iahallfA,iahallfN)
       endif
       if (use_ambi) then
          call ev_update(ev_data_thread,etaambi,        iambiX, iambiA, iambiN )
          call ev_update(ev_data_thread,etaambi*etaart1,iambifX,iambifA,iambifN)
       endif
       call ev_update(ev_data_thread,data_out(7)/data_out(8),inenX,inenA)
       call ev_update(ev_data_thread,      data_out( 7), ineX,   ineA             )
       call ev_update(ev_data_thread,      data_out( 8), innX,   innA             )
       if (ion_rays) then
          call ev_update(ev_data_thread,   data_out( 4),izgrainX,izgrainA,izgrainN)
          call ev_update(ev_data_thread,   data_out( 9),inihrX,   inihrA          )
          call ev_update(ev_data_thread,   data_out(10),inimrX,   inimrA          )
          call ev_update(ev_data_thread,   data_out(13),ingX,     ingA            )
       endif
       if (ion_thermal) then
          call ev_update(ev_data_thread,   data_out(11),inistX,   inistA          )
          call ev_update(ev_data_thread,   data_out(12),inidtX,   inidtA          )
          j = 14
          if (nelements>=2) then
             call ev_update(ev_data_thread,data_out( j),inhX,    inhA             ); j=j+1
             call ev_update(ev_data_thread,data_out( j),inheX,   inheA            ); j=j+1
          endif
          if (nelements>=5) then
             call ev_update(ev_data_thread,data_out( j),innaX,   innaA            ); j=j+1
             call ev_update(ev_data_thread,data_out( j),inmgX,   inmgA            ); j=j+1
             call ev_update(ev_data_thread,data_out( j),inkX,    inkA             ); j=j+1
          endif
          if (nlevels>=2) then
             if (nelements>=2) then
                call ev_update(ev_data_thread,data_out( j),inhedX,   inhedA       ); j=j+1
             endif
             if (nelements>=5) then
                call ev_update(ev_data_thread,data_out( j),innadX,   innadA       ); j=j+1
                call ev_update(ev_data_thread,data_out( j),inmgdX,   inmgdA       ); j=j+1
                call ev_update(ev_data_thread,data_out( j),inkdX,    inkdA        ); j=j+1
             endif
          endif
       endif
    endif
 enddo
!$omp enddo
!$omp critical(collatedata)
 call collate_ev_data(ielements,ev_action,ev_data_thread,ev_data)
!$omp end critical(collatedata)
!$omp end parallel
 !
 !--Finalise the arrays
 if (np > 0.) then
    dnptot = 1./real(np)
 else
    dnptot = 0.
 endif
 call finalise_ev_data(ielements,ev_data,ev_action,dnptot)
 ev_data(itimeni) = time
 !--Print data to file
   write(ev_format,'(a,I2,a)')"(",ielements,"(1pe18.10,1x))"
   write(ievfile,ev_format) ev_data(1:ielements)
 close(ievfile)
 !
 return
end subroutine inform_nimhd
!----------------------------------------------------------------
!+
!  creates a single label for the output file
!+
!----------------------------------------------------------------
subroutine fill_one_label(ev_fmt,label,i,i0,iev)
 integer(kind=1),  intent(in)    :: iev
 integer,          intent(inout) :: i,i0
 character(len=*), intent(in)    :: label,ev_fmt
 !
 i0 = i
 i  = i + 1
 write(ev_label(i0),ev_fmt)i0,trim(label)
 ev_action(i0) = iev
 !
end subroutine fill_one_label
!----------------------------------------------------------------
!+
!  creates up to three lables per input value: Max, Ave, Min
!+
!----------------------------------------------------------------
subroutine fill_xan_label(ev_fmt,label,i,iX,iA,iN)
 integer,          intent(inout) :: i
 integer,optional, intent(inout) :: iX,iA,iN
 character(len=*), intent(in)    :: label,ev_fmt
 character(len=3)                :: lX,lA,lN
 character(len=11)               :: label0
 !
 if (len(label)<=7) then
   lX = "max"
   lA = "ave"
   lN = "min"
 else
   lX = "X"
   lA = "A"
   lN = "N"
 endif
 !
 if (present(iX)) then
    write(label0,'(a,1x,a)')trim(label),trim(lX); call fill_one_label(ev_fmt,label0,i,iX,ievX)
 endif
 if (present(iA)) then
    write(label0,'(a,1x,a)')trim(label),trim(lA); call fill_one_label(ev_fmt,label0,i,iA,ievA)
 endif
 if (present(iN)) then
    write(label0,'(a,1x,a)')trim(label),trim(lN); call fill_one_label(ev_fmt,label0,i,iN,ievN)
 endif
 !
end subroutine fill_xan_label
!----------------------------------------------------------------
!+
!  initiallised the ev_data array to zero or huge (for min)
!+
!----------------------------------------------------------------
subroutine initialise_ev_data(ielements,evaction,evdata)
 integer, intent(in)            :: ielements
 integer(kind=1), intent(in)    :: evaction(:)
 real,            intent(inout) :: evdata(0:inumev)
 integer                        :: i
 !
 evdata = 0.0
 do i = 1,ielements
   if (evaction(i)==ievN) evdata(i) = huge(evdata(i))
 enddo
 !
end subroutine initialise_ev_data
!----------------------------------------------------------------
!+
!  combines the ev_data from the various threads
!+
!----------------------------------------------------------------
subroutine collate_ev_data(ielements,evaction,evdata_thread,evdata)
 integer,         intent(in)    :: ielements
 integer(kind=1), intent(in)    :: evaction(:)
 real,            intent(in)    :: evdata_thread(0:inumev)
 real,            intent(inout) :: evdata(0:inumev)
 integer                        :: i
 !
 do i = 1,ielements
    select case(evaction(i))
    case(ievA)
      evdata(i) = evdata(i) + evdata_thread(i)
    case(ievX)
      evdata(i) = max(evdata(i),evdata_thread(i))
    case(ievN)
      evdata(i) = min(evdata(i),evdata_thread(i))
    end select
 enddo
 !
end subroutine collate_ev_data
!----------------------------------------------------------------
!+
!  update the ev_data_array
!+
!----------------------------------------------------------------
subroutine ev_update(evdata,val,iX,iA,iN)
 integer, optional, intent(in)    :: iX,iA,iN
 real,              intent(in)    :: val
 real,              intent(inout) :: evdata(0:inumev)
 !
 if (present(iX)) evdata(iX) = max(evdata(iX),val)
 if (present(iA)) evdata(iA) =     evdata(iA)+val
 if (present(iN)) evdata(iN) = min(evdata(iN),val)
 !
end subroutine ev_update
!
!----------------------------------------------------------------
!+
!  Performs final generic housekeeping on the ev_data array
!+
!----------------------------------------------------------------
subroutine finalise_ev_data(ielements,evdata,evaction,dnptot)
 integer,         intent(in)    :: ielements
 integer(kind=1), intent(in)    :: evaction(:)
 real,            intent(inout) :: evdata(0:inumev)
 real,            intent(in)    :: dnptot
 integer                        :: i
 !
 do i = 1,ielements
    if (evaction(i)==ievA) evdata(i) = evdata(i)*dnptot
    if ( evdata(i) >  0.01*huge(evdata(i)) .or. &
         evdata(i) < -0.01*huge(evdata(i))) evdata(i) = 0.0
 enddo
 !
end subroutine finalise_ev_data
!----------------------------------------------------------------------!
end module nicil_sup