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
 use nicil, only: nicil_update_nimhd,nicil_get_dudt_nimhd,nicil_translate_error, &
                  nicil_get_halldrift,nicil_get_ambidrift
 use nicil, only: use_ohm,use_hall,use_ambi,n_nden,meanmolmass,zeta_of_rho,n_data_out,n_warn
 use nicil, only: nicil_initialise
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
                                ,ivelX,ivelA,ivelN,ivhallX,ivhallA,ivhallN &
                                ,ivionX,ivionA,ivionN,ivdriftX,ivdriftA,ivdriftN &
                                ,inenX,inenA,inenN,ineX,ineA,innX,innA &
                                ,ingnX,ingnA,ingX,ingA,ingpX,ingpA &
                                ,izetaA,izetaN
  real,             private   :: rhocrit_nimhd,rhocrit2_nimhd,rhocrit3_nimhd
  real,             private   :: coef_gammaeq1,coef_gamma,coef_gam,coef_gamdh,coef_gamah
  character(len=19),private   :: ev_label(inumev)
  integer(kind=1),  private   :: ev_action(inumev)

  !--Subroutines
  public                      :: nimhd_init_gastemp,nimhd_gastemp,nimhd_get_eta,energ_nimhd
  public                      :: nimhd_write_evheader,nimhd_write_ev,nimhd_get_stats
  public                      :: nicil_initialise_dust,nicil_update_nimhd_dust
  private                     :: get_coef_gamma,fill_one_label,fill_xan_label

 private

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
subroutine energ_nimhd(fxyzu4i,rhoi,eta_nimhd,jcurrenti,Bxyzi)
 real,             intent(inout) :: fxyzu4i
 real,             intent(in)    :: rhoi
 real,             intent(in)    :: eta_nimhd(3),jcurrenti(3),Bxyzi(3)
 real                            :: dudtnonideal
 !
 call nicil_get_dudt_nimhd(dudtnonideal,eta_nimhd(1),eta_nimhd(3),rhoi,jcurrenti,Bxyzi)
 fxyzu4i = fxyzu4i + dudtnonideal
 !
end subroutine energ_nimhd
!-----------------------------------------------------------------------
!+
!  A loop to calculated the non-ideal MHD coefficients
!+
!-----------------------------------------------------------------------
subroutine nimhd_get_eta(n1,n2,nlst_in,nlst_tot,vxyzu,ekcle,Bxyz,trho &
                       ,eta_nimhd,nden_nimhd,iphase,iunique,iorig,ivar,gamma,encal,gtime &
                       ,itry_array,nicil_FatalOnly,pass_dust,HY09bin_rho)
 integer,          intent(in)    :: n1,n2,nlst_in,nlst_tot
 real,             intent(in)    :: vxyzu(:,:),ekcle(:,:),Bxyz(:,:)
 real(kind=4),     intent(in)    :: trho(:)
 real,             intent(inout) :: eta_nimhd(:,:),nden_nimhd(:,:)
 integer,          intent(in)    :: iorig(:),ivar(:,:)
 integer(kind=1),  intent(in)    :: iphase(:)
 integer(kind=8),  intent(in)    :: iunique(:)
 integer,          intent(inout) :: itry_array(-1:261)
 real,             intent(in)    :: gamma,gtime
 real(kind=8),     intent(in)    :: HY09bin_rho(:,:)
 character(len=*), intent(in)    :: encal
 logical,          intent(in)    :: nicil_FatalOnly,pass_dust
 integer                         :: ierrlist(n_warn),itry_array0(-1:261)
 integer                         :: n,ipart,itry
 real                            :: Bi,temperature,ekclei

 itry_array0 = 0
!$omp parallel default(none) &
!$omp shared(n1,n2,nlst_in,nlst_tot,vxyzu,ekcle,Bxyz,trho,iphase,iunique,iorig,ivar,gtime) &
!$omp shared(eta_nimhd,nden_nimhd,gamma,encal,nicil_FatalOnly,pass_dust,HY09bin_rho) &
!$omp private(n,ipart,Bi,ekclei,temperature,itry,ierrlist) &
!$omp reduction(+:itry_array0)
!$omp do
 do n = nlst_in, nlst_tot
    ipart = ivar(3,n)
    if (iphase(ipart)==0) then
       Bi = dot_product(Bxyz(1:3,ipart),Bxyz(1:3,ipart))
       Bi = sqrt(Bi)
       if (encal=='r') then
          ekclei = ekcle(3,ipart)
       else
          ekclei = 0.0
       endif
       temperature = nimhd_gastemp(encal,vxyzu(4,ipart),ekclei,real(trho(ipart)),gamma &
                                  ,n1,n2,iunique(iorig(ipart)))
       ierrlist = 0
       if (pass_dust) then
          call nicil_update_nimhd(0,eta_nimhd(1,ipart),eta_nimhd(2,ipart),eta_nimhd(3,ipart), &
                                  Bi,real(trho(ipart)),temperature,nden_nimhd(:,ipart),ierrlist, &
                                  itry=itry,fdg_in=real(HY09bin_rho(:,ipart)))
       else
          call nicil_update_nimhd(0,eta_nimhd(1,ipart),eta_nimhd(2,ipart),eta_nimhd(3,ipart), &
                                  Bi,real(trho(ipart)),temperature,nden_nimhd(:,ipart),ierrlist,itry=itry)
       endif

       itry_array0(itry) = itry_array0(itry) + 1
       if (nicil_FatalOnly) then
          ! print only fatal messages and kill
          if ( any(ierrlist > 0) .and. gtime > 0.) then
            call nicil_translate_error(ierrlist,nicil_FatalOnly)
            call quit(1)
          endif
       else
          ! print all messages and kill if necessary
          if ( any(ierrlist /= 0)  .and. gtime > 0.) then
             call nicil_translate_error(ierrlist,nicil_FatalOnly,real(trho(ipart)),Bi,temperature, &
                                        eta_nimhd(1,ipart),eta_nimhd(2,ipart),eta_nimhd(3,ipart), &
                                        nden_nimhd(:,ipart),real(HY09bin_rho(:,ipart)))
             if ( any(ierrlist >  0) ) call quit(1)
          endif
       endif
    endif
 enddo
!$omp enddo
!$omp end parallel
itry_array = itry_array + itry_array0

end subroutine nimhd_get_eta
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
       call fill_xan_label(ev_fmt,'v_hall',   i,ivhallX, ivhallA, ivhallN )
    endif
    if (use_ambi) then
      call fill_xan_label(ev_fmt,'eta_a',    i,iambiX,  iambiA,  iambiN  )
      call fill_xan_label(ev_fmt,'eta_a/art',i,iambifX, iambifA, iambifN )
      call fill_xan_label(ev_fmt,'velocity', i,ivelX,   ivelA,   ivelN   )
      call fill_xan_label(ev_fmt,'v_ion',    i,ivionX,  ivionA,  ivionN  )
      call fill_xan_label(ev_fmt,'v_drift',  i,ivdriftX,ivdriftA,ivdriftN)
    endif
    call fill_xan_label(ev_fmt,'n_e/n',       i,inenX,   inenA,  inenN    )
    call fill_xan_label(ev_fmt,'n_e',         i,ineX,    ineA             )
    call fill_xan_label(ev_fmt,'n_n',         i,innX,    innA             )
    call fill_xan_label(ev_fmt,'n_g(Z=-1)',   i,ingnX,   ingnA            )
    call fill_xan_label(ev_fmt,'n_g(Z= 0)',   i,ingX,    ingA             )
    call fill_xan_label(ev_fmt,'n_g(Z=+1)',   i,ingpX,   ingpA            )
    if (zeta_of_rho) then
       call fill_xan_label(ev_fmt,'zeta',     i,iA=izetaA, iN=izetaN      )
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
subroutine nimhd_get_stats(imhd2,npart,xyzmh,vxyzu,ekcle,Bxyz,jcurrent,rho,vsound,&
                        alphaMM,alphamin,nden_nimhd,iphase,&
                        gamma,ifsvi,encal,np,et,ev_data,eta_nimhd,n1,n2,iunique,iorig)
 real,             intent(in)    :: gamma
 real,             intent(in)    :: xyzmh(:,:),vxyzu(:,:),ekcle(:,:),Bxyz(:,:),jcurrent(:,:),alphamin(:)
 real(kind=4),     intent(in)    :: rho(:),vsound(:),alphaMM(:,:)
 integer,          intent(in)    :: imhd2,npart,ifsvi,n1,n2
 integer,          intent(in)    :: iorig(:)
 integer(kind=8),  intent(in)    :: iunique(:)
 integer(kind=1),  intent(in)    :: iphase(:)
 character(len=*), intent(in)    :: encal
 integer,          intent(out)   :: np
 real,             intent(in)    :: eta_nimhd(4,imhd2),nden_nimhd(:,:)
 real,             intent(out)   :: et,ev_data(0:inumev)
 integer                         :: i,c0,c1,crate,cmax,iekcle
 integer                         :: ierrlist(n_warn)
 real                            :: rhoi,B2i,Bi,temperature,vsigi,v2i, &
                                    etaart,etaart1,n_total,n_total1
 real                            :: vhalli(3),vdrifti(3),curlBi(3),data_out(n_data_out)
 real                            :: ev_data_thread(0:inumev)
 !
 ! To determine the runtime in this routine; not using sphNG's routine since
 ! this file is compiled first
 call system_clock(c0,crate,cmax)
 !
 call initialise_ev_data(ielements,ev_action,ev_data)
 np     = 0
 iekcle = 1
!$omp parallel default(none) &
!$omp shared(xyzmh,vxyzu,ekcle,Bxyz,jcurrent,rho,vsound,iphase,npart,encal,n1,n2,iunique,iorig) &
!$omp shared(nden_nimhd,alphaMM,alphamin,gamma,ifsvi) &
!$omp shared(use_ohm,use_hall,use_ambi) &
!$omp shared(ielements,ev_data,ev_action,eta_nimhd) &
!$omp shared(itX,itA,itN,ietaFX,ietaFA,iohmX,iohmA,iohmN,iohmfX,iohmfA,iohmfN) &
!$omp shared(ihallX,ihallA,ihallN,iahallX,iahallA,iahallN) &
!$omp shared(ihallfX,ihallfA,ihallfN,iahallfX,iahallfA,iahallfN) &
!$omp shared(iambiX,iambiA,iambiN,iambifX,iambifA,iambifN) &
!$omp shared(ivelX,ivelA,ivelN,ivhallX,ivhallA,ivhallN,ivionX,ivionA,ivionN,ivdriftX,ivdriftA,ivdriftN) &
!$omp shared(inenX,inenA,inenN,ineX,ineA,innX,innA) &
!$omp shared(ingnX,ingnA,ingX,ingA,ingpX,ingpA,izetaA,izetaN) &
!$omp private(i,rhoi,B2i,Bi,temperature,vsigi,etaart,etaart1,ierrlist) &
!$omp private(curlBi,vhalli,vdrifti,v2i) &
!$omp private(data_out,n_total,n_total1) &
!$omp firstprivate(iekcle) &
!$omp private(ev_data_thread) &
!$omp reduction(+:np) 
 call initialise_ev_data(ielements,ev_action,ev_data_thread)
!$omp do
 do i=1,npart
     if ( iphase(i)==0 ) then
       np   = np + 1
       B2i  = dot_product(Bxyz(1:3,i), Bxyz(1:3,i) )
       v2i  = dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
       Bi   = sqrt(B2i)
       rhoi = dble(rho(i))
       if (encal=='r') iekcle = i
       curlBi = jcurrent(1:3,i)
       temperature = nimhd_gastemp(encal,vxyzu(4,i),ekcle(3,iekcle),rhoi,gamma,n1,n2,iunique(iorig(i)))
       call nicil_get_halldrift(eta_nimhd(2,i),Bxyz(1,i),Bxyz(2,i),Bxyz(3,i),curlBi,vhalli)
       call nicil_get_ambidrift(eta_nimhd(3,i),Bxyz(1,i),Bxyz(2,i),Bxyz(3,i),curlBi,vdrifti)

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
          call ev_update(ev_data_thread,eta_nimhd(1,i),        iohmX, iohmA, iohmN )
          call ev_update(ev_data_thread,eta_nimhd(1,i)*etaart1,iohmfX,iohmfA,iohmfN)
       endif
       if (use_hall) then
          call ev_update(ev_data_thread,    eta_nimhd(2,i),              ihallX,  ihallA,  ihallN  )
          call ev_update(ev_data_thread,abs(eta_nimhd(2,i)),             iahallX, iahallA, iahallN )
          call ev_update(ev_data_thread,    eta_nimhd(2,i)*etaart1,      ihallfX, ihallfA, ihallfN )
          call ev_update(ev_data_thread,abs(eta_nimhd(2,i))*etaart1,     iahallfX,iahallfA,iahallfN)
          call ev_update(ev_data_thread,sqrt(dot_product(vhalli,vhalli)),ivhallX, ivhallA, ivhallN )

       endif
       if (use_ambi) then
          call ev_update(ev_data_thread,eta_nimhd(3,i),                    iambiX,  iambiA,  iambiN  )
          call ev_update(ev_data_thread,eta_nimhd(3,i)*etaart1,            iambifX, iambifA, iambifN )
          call ev_update(ev_data_thread,sqrt(v2i),                         ivelX,   ivelA,   ivelN   )
          call ev_update(ev_data_thread,sqrt(dot_product(vdrifti,vdrifti)),ivdriftX,ivdriftA,ivdriftN)
       endif

       !--data_out is no longer calculated here; this is left here in case we want to rea-dd it at a later date
       !n_total = data_out(8) + data_out(5)  ! n_electron + n_neutral
       !if (n_total > 0.) then
       !   n_total1 = 1.0/n_total
       !else
       !   n_total1 = 0.0         ! only possible if eta_constant = .true.
       !endif
       !eta_nimhd(4,i) = data_out(8)*n_total1    ! n_electron / (n_electron + n_neutral);  Save ionisation fraction for the dump file
       !call ev_update(ev_data_thread,data_out( 8)*n_total1,inenX,inenA,inenN)
       !call ev_update(ev_data_thread,data_out( 8),         ineX,   ineA     )
       !call ev_update(ev_data_thread,data_out( 5),         innX,   innA     )
       !call ev_update(ev_data_thread,data_out(20),         ingnX,  ingpA    )
       !call ev_update(ev_data_thread,data_out(21),         ingX,   ingA     )
       !call ev_update(ev_data_thread,data_out(22),         ingpX,  ingpA    )
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
   if (evaction(i)==ievX) evdata(i) = -huge(evdata(i))
   if (evaction(i)==ievN) evdata(i) =  huge(evdata(i))
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
!----------------------------------------------------------------
!+
!  wrapper initialisation routine since Fortran77 does not include optional arguments
!+
!----------------------------------------------------------------
subroutine nicil_initialise_dust(utime,udist,umass,umagfd,HY09binsizes,HY09_ndust_bins,nicil_passdust,ierr)
 use nicil, only: na,use_fdg_in
 real,    intent(in)  :: utime,umass,udist,umagfd
 integer, intent(in)  :: HY09_ndust_bins
 integer, intent(out) :: ierr
 real,    intent(in)  :: HY09binsizes(HY09_ndust_bins)
 logical, intent(out) :: nicil_passdust
 integer              :: i

 if (HY09_ndust_bins /= na) then
    print*, "NICIL & sphNG are using a different number of grains: ",na,HY09_ndust_bins
    print*, "NICIL will use it's internal dust properties"
    call nicil_initialise(utime,udist,umass,umagfd,ierr)
 else
    nicil_passdust = .true.
    use_fdg_in     = .true.
    call nicil_initialise(utime,udist,umass,umagfd,ierr,a_grain_cgs_in=HY09binsizes)
    print*, "NICIL: Passing in dust properties from sphNG."
    open(unit=74205,file='nicil_input_dust_sizes.dat')
    do i = 1,HY09_ndust_bins
       write(74205,*) i,HY09binsizes(i)
    enddo
    close(74205)
 endif

end subroutine nicil_initialise_dust
!----------------------------------------------------------------
!+
!  wrapper routine for main call since Fortran77 does not include optional arguments
!+
!----------------------------------------------------------------
subroutine nicil_update_nimhd_dust(icall,eta_ohm,eta_hall,eta_ambi, &
                                   Bfield,rho,temperature,gtime,nden_save,fdg_in,nicil_FatalOnly,call_quit,ierrlist)
 integer, intent(in)    :: icall
 integer, intent(inout) :: ierrlist(n_warn)
 real,    intent(out)   :: eta_ohm,eta_hall,eta_ambi
 real,    intent(in)    :: Bfield,rho,temperature,gtime,fdg_in(:)
 real,    intent(inout) :: nden_save(:)
 logical, intent(in)    :: nicil_FatalOnly
 logical, intent(out)   :: call_quit

 call nicil_update_nimhd(icall,eta_ohm,eta_hall,eta_ambi,Bfield,rho,temperature,nden_save,ierrlist,fdg_in=fdg_in)
 if (nicil_FatalOnly) then
    ! print only fatal messages and kill
    if ( any(ierrlist > 0) .and. gtime > 0.) then
       call nicil_translate_error(ierrlist,nicil_FatalOnly,rho=rho,T=temperature,fdg_in=fdg_in)
       print*, 'nicilsup',fdg_in
       call_quit = .true.
    endif
 else
    ! print all messages and kill if necessary
    if ( any(ierrlist /= 0)  .and. gtime > 0.) then
       call nicil_translate_error(ierrlist,nicil_FatalOnly,rho=rho,T=temperature,fdg_in=fdg_in)
       if ( any(ierrlist >  0) ) call_quit = .true.
       print*, 'nicilsup',fdg_in
       call_quit = .true.
    endif
 endif

end subroutine nicil_update_nimhd_dust
!----------------------------------------------------------------

end module nicil_sup
