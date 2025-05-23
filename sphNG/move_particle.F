      SUBROUTINE move_particle(ifrom,ito)
c************************************************************
c                                                           *
c     Moves a particle from location ifrom to ito (within   *
c     an MPI process, not between MPI processes)            *
c                                                           *
c     NOTE: Does not alter listbins for timesteps or ireal  *
c     or hasghost for ghost particles.                      *
c     Does not alter tree.                                  *
c                                                           *
c************************************************************

#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/f2'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/accnum'
      INCLUDE 'COMMONS/accurpt'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/radsink'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/dustimplicit'
      INCLUDE 'COMMONS/grain'
      INCLUDE 'COMMONS/compv'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09rho'
      INCLUDE 'COMMONS/HY09dumrho'
      INCLUDE 'COMMONS/f1HY09'
      INCLUDE 'COMMONS/f2HY09'
      INCLUDE 'COMMONS/HY09check'
      INCLUDE 'COMMONS/HY09accel'
      INCLUDE 'COMMONS/dustfluid'
      INCLUDE 'COMMONS/dustfluidderiv'
      INCLUDE 'COMMONS/dustfluidvelu'
      INCLUDE 'COMMONS/tstopvar'

      IF (ito.EQ.ifrom) RETURN

      DO i = 1, 5
         xyzmh(i,ito) = xyzmh(i,ifrom)
         dumxyzmh(i,ito) = dumxyzmh(i,ifrom)
      END DO
      DO i = 1, 4
         vxyzu(i,ito) = vxyzu(i,ifrom)
         dumvxyzu(i,ito) = dumvxyzu(i,ifrom)
         f1vxyzu(i,ito) = f1vxyzu(i,ifrom)
         f2vxyzu(i,ito) = f2vxyzu(i,ifrom)         
      END DO
      IF (idim_grow.EQ.idim) THEN
        rgrain(ito) = rgrain(ifrom)
        dumrgrain(ito) = dumrgrain(ifrom)
        f1rgrain(ito) = f1rgrain(ifrom)
        f2rgrain(ito) = f2rgrain(ifrom)
        vturb(ito) = vturb(ifrom)
        vcomp(ito) = vcomp(ifrom)
        rhoG(ito) =  rhoG(ifrom)
      ENDIF
        
      isteps(ito) = isteps(ifrom)
      IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1) THEN
         icolumnsteps(ito) = icolumnsteps(ifrom)
         icolumnnext(ito) = icolumnnext(ifrom)
      ENDIF
      it0(ito) = it0(ifrom)
      it1(ito) = it1(ifrom)
      it2(ito) = it2(ifrom)
      iphase(ito) = iphase(ifrom)
      IF (iphase(ifrom).GE.1 .AND. iphase(ifrom).LT.10) THEN
         ipt = listrealpm(ifrom)
         IF (ipt.LE.0 .OR. ipt.GT.nptmass) THEN
            WRITE (*,*) 'ERROR - ipt.LE.0 .OR. ipt.GT.nptmass'
            CALL quit(1)
         ENDIF
         listrealpm(ito) = ipt
         listpm(ipt) = ito
         listrealpm(ifrom) = 0
      ENDIF
      iphase(ifrom) = -1
c      iorig(ito) = ito
c      isort(ito) = ito
c      iorig(ifrom) = 0
c      isort(ifrom) = 0
      iunique(iorig(ito)) = iunique(iorig(ifrom))
      iunique(iorig(ifrom)) = 0
      rho(ito) = rho(ifrom)
      dumrho(ito) = dumrho(ifrom)
      divv(ito) = divv(ifrom)
      curlv(ito) = curlv(ifrom)
      ddv(ito) = ddv(ifrom)
      pr(ito) = pr(ifrom)
      vsound(ito) = vsound(ifrom)
      poten(ito) = poten(ifrom)
      IF (idim_h2.EQ.idim) THEN
c         h2mol(ito) = h2mol(ifrom)
         h2ratio(ito) = h2ratio(ifrom)
      ENDIF
      IF (idustIMPL.GT.0) dustnorm(ito) = dustnorm(ifrom)
      IF (idustFluid.NE.0) THEN
         dustvar(:,ito) = dustvar(:,ifrom)
         dumdustvar(:,ito) = dumdustvar(:,ifrom)
         f1dust(:,ito) = f1dust(:,ifrom)
         f2dust(:,ito) = f2dust(:,ifrom)
         dustfluidvxyz(:,:,ito) = dustfluidvxyz(:,:,ifrom)
         stoppingtime(:,ito) = stoppingtime(:,ifrom)
         omega_disc(ito) = omega_disc(ifrom)
      ENDIF
      IF (idimHY09.EQ.idim) THEN
         HY09bin_rho(:,ito) = HY09bin_rho(:,ifrom)
         HY09bin_dumrho(:,ito) = HY09bin_dumrho(:,ifrom)
         f1HY09drhodt(:,ito) = f1HY09drhodt(:,ifrom)
         f2HY09drhodt(:,ito) = f2HY09drhodt(:,ifrom)
         HY09frac_h_moved(ito) = HY09frac_h_moved(ifrom)
         f1HY09vrel_h(ito) = f1HY09vrel_h(ifrom)
         gas_accel(:,ito) = gas_accel(:,ifrom)
         grav_accel(:,ito) = grav_accel(:,ifrom)
         HY09discflag(ito) = HY09discflag(ifrom)
      ENDIF
      iscurrent(ito) = iscurrent(ifrom)
      notacc(ito) = notacc(ifrom)
      DO i = 1, 2
         gradhs(i,ito) = gradhs(i,ifrom)
      END DO
      DO i = 1, isizealphaMM
         alphaMM(i,ito) = alphaMM(i,ifrom)
         dumalpha(i,ito) = dumalpha(i,ifrom)
      END DO
      DO i = 1, 1+isizealphaMM
         f1ha(i,ito) = f1ha(i,ifrom)
         f2ha(i,ito) = f2ha(i,ifrom)
      END DO
      IF (imhd.EQ.idim) THEN
         DO i = 1, imhdevol
            Bevolxyz(i,ito) = Bevolxyz(i,ifrom)
            dumBevolxyz(i,ito) = dumBevolxyz(i,ifrom)
            f1Bxyz(i,ito) = f1Bxyz(i,ifrom)
            f2Bxyz(i,ito) = f2Bxyz(i,ifrom)
         END DO
         DO i = 1, 3
            Bxyz(i,ito) = Bxyz(i,ifrom)
         END DO
         DO i = 1, 5
            divcurlB(i,ito) = divcurlB(i,ifrom)
         END DO
#ifdef NONIDEAL
         DO i = 1, 3
            jcurrent(i,ito) = jcurrent(i,ifrom)
         END DO
         DO i = 1,n_nden
            nden_nimhd(i,ito) = nden_nimhd(i,ifrom)
         END DO
         DO i = 1,4
            eta_nimhd(i,ito) = eta_nimhd(i,ifrom)
         END DO
#endif
      ENDIF
      IF (iradtrans.EQ.idim) THEN
         DO i = 1, 5
            ekcle(i,ito) = ekcle(i,ifrom)
            dumekcle(i,ito) = dumekcle(i,ifrom)
         END DO
         IF (idustRT.GT.0) THEN
            DO i = 1, 2
               dust_tk(i,ito) = dust_tk(i,ifrom)
            END DO
            DO i = 1, nheatingISR
               heatingISR(i,ito) = heatingISR(i,ifrom)
            END DO
            IF (ioptimise_column.EQ.1) THEN
               DO i = 1, 2
                  heatingISRold(i,ito) = heatingISRold(i,ifrom)
               END DO
            ENDIF
            DO i = 1, nchemistry
               chemistry(i,ito) = chemistry(i,ifrom)
            END DO
            h2frac(ito) = h2frac(ifrom)
            DO i = 1, 2
               dh2dt(i,ito) = dh2dt(i,ifrom)
            END DO
         ENDIF
      ENDIF
      IF (istellar_ionisation) THEN
         HIIion(ito) = HIIion(ifrom)
         dHIIdt(1:2,ito) = dHIIdt(1:2,ifrom)
         stellarrad(1:4,1:iptdim,ito) =stellarrad(1:4,1:iptdim,ifrom)
      ENDIF

      RETURN
      END

c===========================================================================
