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
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif
      INCLUDE 'COMMONS/stellarradiation'

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
         DO i = 1,4
            n_R(i,ito) = n_R(i,ifrom)
         END DO
         n_electronT(ito) = n_electronT(ifrom)
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
