!----------------------------------------------------------------------!
!+
!  MODULE: nicil_sup
!          
!  DESCRIPTION:
!  Contains supplementary routines to NICIL for use in sphNG; adapted 
!  from Phantom's evwrite & energies.
!  This is not part of the NICIL package.
!
!  AUTHOR: James Wurster
!+
!----------------------------------------------------------------------!
module nicil_sup
 use nicil, only: nicil_get_eta,nimhd_get_dudt
 use nicil, only: use_ohm,use_hall,use_ambi,ion_rays,ion_thermal, &
                  nelements,nelements_max,nlevels,meanmolmass, &
                 zeta_of_rho,n_data_out
 implicit none
 !
 !--Indicies to determine which action to take
 integer(kind=1),   parameter :: iev0 = 0
 integer(kind=1),   parameter :: ievN = 1
 integer(kind=1),   parameter :: ievA = 2  
 integer(kind=1),   parameter :: ievX = 3
 integer,           parameter :: ievfile   = 74205
 integer,    public,parameter :: inumev    = 99
 integer,    public           :: ielements = 0
 !
 integer,           private   :: itimeni,itX,itA,itN,ietaFX,ietaFA &
                                ,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN &
                                ,ihallX,ihallA,ihallN,iahallX,iahallA,iahallN &
                                ,ihallfX,ihallfA,ihallfN &
                                ,iahallfX,iahallfA,iahallfN &
                                ,iambiX,iambiA,iambiN,iambifX,iambifA,iambifN &
                                ,inenX,inenA,inenN,ineX,ineA,innX,innA &
                                ,inihrX,inihrA,inimrX,inimrA,ingnX,ingnA,ingX,ingA,ingpX,ingpA &
                                ,inistX,inistA,inidtX,inidtA,izetaA,izetaN
  real,             private   :: rhocrit_nimhd,rhocrit2_nimhd,rhocrit3_nimhd
  real,             private   :: coef_gammaeq1,coef_gamma,coef_gam,coef_gamdh,coef_gamah
  character(len=19),private   :: ev_label(inumev)
  integer(kind=1),  private   :: ev_action(inumev)
  !--Subroutines
  public                      :: nimhd_init_gastemp,nimhd_gastemp,energ_nimhd
  public                      :: nimhd_write_evheader,nimhd_write_ev,nimhd_get_stats
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
real function nimhd_gastemp(encal,vxyzu4i,ekcle3i,rhoi,gamma,n1,n2,ipart0)
 integer,          intent(in) :: n1,n2
 integer(kind=8),  intent(in) :: ipart0
 real,             intent(in) :: vxyzu4i,ekcle3i,rhoi,gamma
 character(len=*), intent(in) :: encal
 !
 if (encal=='r') then
    if (ipart0 > n1 .AND. ipart0 <= n1+n2) then
       ! for sphere-in-box, this is the medium
       nimhd_gastemp = coef_gamma*vxyzu4i
    else
       ! for sphere-in-box, this is the sphere
       if (ekcle3i > 0.0) then
          nimhd_gastemp = vxyzu4i/ekcle3i
       else
          nimhd_gastemp = 0.0  ! Should not happen if parent code is written correctly
       endif
    endif
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
!-----------------------------------------------------------------------
!+
!  wrapper routine to calculate and include du/dt
!+
!-----------------------------------------------------------------------
subroutine energ_nimhd(fxyzu4i,jcurrenti,Bxyzi,vxyzu4i,ekcle3i,rhoi,n_Ri,n_electronTi,gamma,encal,n1,n2,ipart0)
 integer,          intent(in)    :: n1,n2
 integer(kind=8),  intent(in)    :: ipart0
 real,             intent(inout) :: fxyzu4i
 real,             intent(in)    :: jcurrenti(3),Bxyzi(3),n_Ri(:)
 real,             intent(in)    :: vxyzu4i,ekcle3i,rhoi,n_electronTi,gamma
 character(len=*), intent(in)    :: encal
 integer                         :: ierr
 real                            :: Bi,temperaturei,etaohmi,etahalli,etaambii
 real                            :: dudtnonideal
 !
 Bi           = sqrt( dot_product(Bxyzi,Bxyzi) )
 temperaturei = nimhd_gastemp(encal,vxyzu4i,ekcle3i,rhoi,gamma,n1,n2,ipart0)
 call nicil_get_eta(etaohmi,etahalli,etaambii,Bi,rhoi,temperaturei,n_Ri,n_electronTi,ierr)
 call nimhd_get_dudt(dudtnonideal,etaohmi,etaambii,rhoi,jcurrenti,Bxyzi)
 fxyzu4i      = fxyzu4i + dudtnonideal
 !
end subroutine energ_nimhd
!-----------------------------------------------------------------------
!+
!  Calculates Header and initialises indicies; initialises .ev file if
!  if does not exist.
!+
!-----------------------------------------------------------------------
subroutine nimhd_write_evheader(nimhdfile,iproc)
 integer,          intent(in)  :: iproc
 character(len=*), intent(in)  :: nimhdfile
 integer                       :: i
 logical                       :: iexists
 character(len=  27)           :: ev_fmt
 !
 inquire(file=nimhdfile,exist=iexists)
 if ( iexists .and. ielements > 0 ) then
    return
 else
    ! Determine what values should be in the header
    i = 1
    write(ev_fmt,'(a)') "(1x,'[',i2.2,1x,a11,']',2x)"
    call fill_one_label(ev_fmt,'time',        i,itimeni, iev0             )
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
    call fill_xan_label(ev_fmt,'n_e/n',       i,inenX,   inenA,  inenN    )
    call fill_xan_label(ev_fmt,'n_e',         i,ineX,    ineA             )
    call fill_xan_label(ev_fmt,'n_n',         i,innX,    innA             )
    if (ion_rays) then
       call fill_xan_label(ev_fmt,'n_ihR',    i,inihrX,  inihrA           )
       call fill_xan_label(ev_fmt,'n_imR',    i,inimrX,  inimrA           )
       call fill_xan_label(ev_fmt,'n_gR',     i,ingX,    ingA             )
       call fill_xan_label(ev_fmt,'n_g(Z=-1)',i,ingnX,   ingnA            )
       call fill_xan_label(ev_fmt,'n_g(Z= 0)',i,ingX,    ingA             )
       call fill_xan_label(ev_fmt,'n_g(Z=+1)',i,ingpX,   ingpA            )
       if (zeta_of_rho) then
          call fill_xan_label(ev_fmt,'zeta',  i,iA=izetaA, iN=izetaN      )
       endif
    endif
    if (ion_thermal) then
       call fill_xan_label(ev_fmt,'n_isT',    i,inistX,  inistA           )
       call fill_xan_label(ev_fmt,'n_idT',    i,inidtX,  inidtA           )
    endif
    ielements = i - 1 ! The number of columns to be calculates
    !
    ! Write a header line if file does not exists (else this was done to initialise the indicies)
    if ( .not. iexists .and. iproc==0 ) then
       open(ievfile, file=nimhdfile,status='replace')
       write(ev_fmt,'(a,I3,a)') '(',ielements+1,'a)'
       write(ievfile,ev_fmt)'#',ev_label(1:ielements)
       close(ievfile)
    endif
 end if
 !
end subroutine nimhd_write_evheader
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD characteristics of interest
!+
!-----------------------------------------------------------------------
subroutine nimhd_get_stats(imhd2,npart,xyzmh,vxyzu,ekcle,Bxyz,rho,vsound,&
                        alphaMM,alphamin,n_R,n_electronT,iphase,&
                        gamma,ifsvi,encal,np,et,ev_data,ionfrac_eta,n1,n2,iunique,iorig)
 real,             intent(in)  :: gamma
 real,             intent(in)  :: xyzmh(:,:),vxyzu(:,:),ekcle(:,:),Bxyz(:,:),&
                                  n_R(:,:),n_electronT(:),alphamin(:)
 real(kind=4),     intent(in)  :: rho(:),vsound(:),alphaMM(:,:)
 integer,          intent(in)  :: imhd2,npart,ifsvi,n1,n2
 integer,          intent(in)  :: iorig(:)
 integer(kind=8),  intent(in)  :: iunique(:)
 integer(kind=1),  intent(in)  :: iphase(:)
 character(len=*), intent(in)  :: encal
 integer,          intent(out) :: np
 real,             intent(out) :: et,ev_data(0:inumev)
 real(kind=4),     intent(out) :: ionfrac_eta(4,imhd2)
 integer                       :: i,ierr,c0,c1,crate,cmax,iekcle
 real                          :: rhoi,B2i,Bi,temperature,vsigi, &
                                  etaart,etaart1,etaohm,etahall,etaambi,n_total
 real                          :: data_out(n_data_out)
 real                          :: ev_data_thread(0:inumev)
 !
 ! To determine the runtime in this routine; not using sphNG's routine since
 ! this file is compiled first
 call system_clock(c0,crate,cmax)
 !
 call initialise_ev_data(ielements,ev_action,ev_data)
 np     = 0
 iekcle = 1
!$omp parallel default(none) &
!$omp shared(xyzmh,vxyzu,ekcle,Bxyz,rho,vsound,iphase,npart,encal,n1,n2,iunique,iorig) &
!$omp shared(n_R,n_electronT,alphaMM,alphamin,gamma,ifsvi) &
!$omp shared(use_ohm,use_hall,use_ambi,ion_rays,ion_thermal,nelements) &
!$omp shared(ielements,ev_data,ev_action,ionfrac_eta) &
!$omp shared(itX,itA,itN,ietaFX,ietaFA,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN) &
!$omp shared(ihallX,ihallA,ihallN,iahallX,iahallA,iahallN) &
!$omp shared(ihallfX,ihallfA,ihallfN,iahallfX,iahallfA,iahallfN) &
!$omp shared(iambiX,iambiA,iambiN,iambifX,iambifA,iambifN) &
!$omp shared(inenX,inenA,ineX,ineA,innX,innA,inihrX,inihrA,inimrX,inimrA) &
!$omp shared(ingnX,ingnA,ingX,ingA,ingpX,ingpA,inistX,inistA,inidtX,inidtA,izetaA,izetaN) &
!$omp private(i,ierr,rhoi,B2i,Bi,temperature,vsigi,etaart,etaart1,etaohm,etahall,etaambi) &
!$omp private(data_out,n_total) &
!$omp firstprivate(iekcle) &
!$omp private(ev_data_thread) &
!$omp reduction(+:np) 
 call initialise_ev_data(ielements,ev_action,ev_data_thread)
!$omp do
 do i=1,npart
     if ( iphase(i)==0 ) then
       np   = np + 1
       B2i  = dot_product(Bxyz(1:3,i),Bxyz(1:3,i))
       Bi   = sqrt(B2i)
       rhoi = dble(rho(i))
       if (encal=='r') iekcle = i
       temperature = nimhd_gastemp(encal,vxyzu(4,i),ekcle(3,iekcle),rhoi,gamma,n1,n2,iunique(iorig(i)))
       call nicil_get_eta(etaohm,etahall,etaambi,Bi,rhoi,temperature, &
                          n_R(:,i),n_electronT(i),ierr,data_out)
       vsigi = sqrt( dble(vsound(i))*dble(vsound(i)) + B2i/rhoi)
       if (ifsvi >= 6) then
          etaart = dble(alphaMM(2,i))*vsigi*xyzmh(5,i)
       else
          etaart = alphamin(2)*vsigi*xyzmh(5,i)
       endif
       if (etaart > 0.0.and. etaart < huge(etaart)) then
          etaart1 = 1.0/etaart
       else
          etaart  = 0.0
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
       n_total = data_out(7) + data_out(8) + data_out(9) + data_out(10) + data_out(11)
       call ev_update(ev_data_thread,data_out(6)/n_total,inenX,inenA,inenN)
       ionfrac_eta(1,i) = real(data_out(6)/n_total,kind=4)
       ionfrac_eta(2,i) = real(etaohm, kind=4)  ! Save eta_OR for the dump file
       ionfrac_eta(3,i) = real(etahall,kind=4)  ! Save eta_HE for the dump file
       ionfrac_eta(4,i) = real(etaambi,kind=4)  ! Save eta_AD for the dump file
       call ev_update(ev_data_thread,      data_out( 6), ineX,   ineA             )
       call ev_update(ev_data_thread,      data_out( 7), innX,   innA             )
       if (ion_rays) then
          call ev_update(ev_data_thread,   data_out( 8),inihrX,   inihrA          )
          call ev_update(ev_data_thread,   data_out( 9),inimrX,   inimrA          )
          call ev_update(ev_data_thread,   data_out(12),ingnX,    ingpA           )
          call ev_update(ev_data_thread,   data_out(13),ingX,     ingA            )
          call ev_update(ev_data_thread,   data_out(14),ingpX,    ingpA           )
       endif
       if (ion_thermal) then
          call ev_update(ev_data_thread,   data_out(10),inistX,   inistA          )
          call ev_update(ev_data_thread,   data_out(11),inidtX,   inidtA          )
       endif
    endif
 enddo
!$omp enddo
!$omp critical(collatedata)
 call collate_ev_data(ielements,ev_action,ev_data_thread,ev_data)
!$omp end critical(collatedata)
!$omp end parallel
 !
 ! to prevent printing useless data
 do i = 1,ielements
    if ( ev_data(i) >  0.01*huge(ev_data(i)) .or. &
         ev_data(i) < -0.01*huge(ev_data(i))) then
        ev_data(i) = 0.0
    end if
 enddo
 !
 !--To complete the timing
 call system_clock(c1,crate,cmax)
 et = (c1-c0)/real(crate)
 !
end subroutine nimhd_get_stats
!-----------------------------------------------------------------------
!+
!  Writes non-ideal MHD data to file
!+
!-----------------------------------------------------------------------
subroutine nimhd_write_ev(time,np,iprint,et,nimhdfile,ev_data_min,ev_data_max,ev_data_sum)
 integer,          intent(in) :: np,iprint
 real,             intent(in) :: time,et
 real,             intent(in) :: ev_data_min(0:inumev),ev_data_max(0:inumev),ev_data_sum(0:inumev)
 character(len=*), intent(in) :: nimhdfile
 integer                      :: i
 real                         :: dnptot
 real                         :: ev_data(0:inumev)
 character(len=  35)          :: ev_format
 !
 ! Average the relevant elements in the ev_data array
 if (np > 0.) then
    dnptot = 1./real(np)
 else
    dnptot = 0.
 endif
 ! Return data to common array and average when necessary
 do i = 1,ielements
    select case(ev_action(i))
    case(ievA)
      ev_data(i) = ev_data_sum(i)*dnptot
    case(ievX)
      ev_data(i) = ev_data_max(i)
    case(ievN)
      ev_data(i) = ev_data_min(i)
    end select
 enddo
 ! Set the time
 ev_data(itimeni) = time
 !
 !--Open file and write data
 ! (it exist since this was previously checked in nimhd_write_evheader)
 open(ievfile, file=nimhdfile,position='append')
 write(ev_format,'(a,I2,a)')"(",ielements,"(1pe18.10,1x))"
 write(ievfile,ev_format) ev_data(1:ielements)
 call flush(ievfile)
 close(ievfile)
 !
 !--To state the time taken
 write(iprint,'(a,Es12.3,a)') & 
    ' nicil: inform_nimhd: The max elapsed time on a processor is ',et,' seconds.'
 if (et.gt.120.0) then
    write(iprint,'(a,Es12.3,a)') &
    ' nicil: inform_nimhd: The max elapsed time on a processor is ',et/60.0,' minutes.'
 end if
 !
 return
end subroutine nimhd_write_ev
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
!----------------------------------------------------------------------!
end module nicil_sup
