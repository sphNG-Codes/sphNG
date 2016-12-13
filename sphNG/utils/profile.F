      PROGRAM profile
c***********************************************************
c                                                          *
c  This program reads sphNG dump files and calculates      *
c     radial profiles of various quantities.               *
c                                                          *
c***********************************************************

#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'

      REAL*8 udist, umass, utime, umagfd, udens, uvel, uergg

      REAL*4 Jx, Jy, Jz
      REAL Jmag,JcrossBx,JcrossBy,JcrossBz,JcrossB,lorentzmean,
     &     lorentzmeanvert,lorentzmeanvertout,JcrossB_dot_rhat

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/gtime'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/vargam'
      INCLUDE 'COMMONS/physcon'
#endif

      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      COMMON /gridv/ iline, isec
      COMMON /origi/ cmx, cmy, cmz
      COMMON /files/ fileout, filein(100), inname
      COMMON /optio/ what, idump, icont
      COMMON /forma/ ien
      COMMON /eost/ temp(idim), entropy_per_baryon(idim)
      COMMON /mhd/ Bmag(idim),
     &     Jx(idim),Jy(idim),Jz(idim),Jmag(idim),
     &     JcrossBx(idim),JcrossBy(idim),JcrossBz(idim),JcrossB(idim),
     &     JcrossB_dot_rhat(idim)

      REAL*4 alphaMM(idim)

#ifdef NONIDEAL
      REAL nebynmean,nzminusmean,nzzeromean,nzplusmean,etaohmmean,
     &     etahallmean,etaambimean
      REAL  data_out(17+nelements_max*nlevels-3)
#endif

      DIMENSION radius(idim), rvel(idim), rvelt(idim),
     &     vsound(idim), rxy(idim)

      CHARACTER*100 fileident
      DIMENSION indx(idim),list(idim)

      CHARACTER*7 filein
      CHARACTER*11 inname, imageout
      CHARACTER*11 fileout, fileout2, fileout3, fileout4
      CHARACTER*9  name
      CHARACTER*1  ien, icont, icor, icent, isink, itorq, itype
      CHARACTER*1  ivarbnd, inewdump, itstep

      DATA pi/3.141592654/

      varmhd = 'Brho'
c
c--read options
c
      WRITE (*,*) 'Enter number of input files:'
      READ (*, *) number_files
      
      WRITE (*,*) 'Enter names of input files:'
      DO k=1, number_files
         READ(*,90) filein(k)
 90      FORMAT(A7)
      END DO

 110  FORMAT(A1)
c
c--loop over files present
c
      WRITE (*,210)
 210  FORMAT('Do you want to centre on the centre of mass (m) or',
     &     ' on the maximum density (d) or not at all?')
      READ (*,110) icent

      WRITE (*,*) 'Undo rotating reference frame (y/n) ?'
      READ (*,110) icor

      IF (icor.EQ.'y') THEN
         WRITE (*,*) 'Enter rotation frequency (sec) ?'
         READ (*, *) rotfreq
      ENDIF

      WRITE (*,*) 'Average in sphere, or only within X degs ',
     &     'of plane(s/a)?'
      READ (*,110) itype

      IF (itype.EQ.'a') THEN
         READ (*,*) angletake
      ENDIF

      iav = 1
      WRITE (*,*) 'Average according to r or rxy (1,2)'
      READ (*,*) iav

      WRITE (*,*)'Minimum number of particles to average over (e.g. 20)'
      READ (*,*) naveragemin

      radfac = 1.259
      radmaxstart = 100.
      WRITE (*,*) 'Enter radfac (e.g. 1.259)'
      READ (*,*) radfac

      WRITE (*,*) 'Enter radmax (e.g. 7.0)'
      READ (*,*) radmax
      radmaxstart = radmax

      WRITE (*,*) 'Enter radin (e.g. 0.999)'
      READ (*,*) radin

      OPEN(19,file='interest')

c
c--Process data files
c
      DO k=1,number_files

         fileout = 'QG' // filein(k)
         fileout2 = 'QD' // filein(k)
         fileout3 = 'QA' // filein(k)
         fileout4 = 'QN' // filein(k)
         imageout = 'IP' // filein(k)
         maxrec = 100*idim
         OPEN(15,file=imageout)

         OPEN (UNIT = 11, FILE = filein(k), FORM = 'unformatted',
     &        RECL=maxrec)

         WRITE (*,*) 'File open'

         CALL rdump(11, ichkl, 0)

         udist = udisti
         umass = umassi
         utime = utimei
         umagfd = umagfdi
         udens = umass/udist**3
         uvel = udist/utime
         uergg = udist**2/utime**2
#ifdef NONIDEAL
         CALL nicil_initialise(utime,udist,umass,umagfd,ierr)
         CALL nimhd_init_gastemp(Rg,gamma,gam,gamdh,gamah,uergg
     &                          ,rhocrit,rhocrit2,rhocrit3)
#endif
c
c--process file dumps
c
         ifile = 1
         iadd = 0
         izero = 0
         WRITE(fileout, 99002) fileout, iadd, ifile
         WRITE(fileout2, 99002) fileout2, iadd, ifile
         WRITE(fileout3, 99002) fileout3, iadd, ifile
         WRITE(fileout4, 99002) fileout4, iadd, ifile
99002    FORMAT(A9, I1, I1)

   10    OPEN(21,file=fileout)
         WRITE(*,99004) fileout
99004    FORMAT(A11)

         OPEN(22,file=fileout2)
         IF (nptmass.GE.1 .AND. (isink.EQ.'y' .OR. isink.EQ.'Y')) THEN
            OPEN(23,file=fileout3)
         ENDIF
#ifdef NONIDEAL
         OPEN(24,file=fileout4)
#endif
         WRITE(*,*) npart, rhozero, nptmass,xyzmh(4,1),rho(3),iphase(1)
         WRITE(*,*) 'RK2=',RK2,' gamma=',gamma,' rhozero=',rhozero
         WRITE(*,*) 'Units=',udist, umass, utime

         IF (icor.EQ.'y') THEN
            romega = rotfreq*utime
            WRITE(*,*) romega
            DO i = 1, npart
               r = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
               th = ATAN2(xyzmh(2,i),xyzmh(1,i))
               th = th + gt*romega
               xyzmh(1,i) = r*COS(th)
               xyzmh(2,i) = r*SIN(th)

               r = SQRT(vxyzu(1,i)**2 + vxyzu(2,i)**2)
               th = ATAN2(vxyzu(2,i),vxyzu(1,i))
               th = th + gt*romega
               vxyzu(1,i) = r*COS(th)
               vxyzu(2,i) = r*SIN(th)
            END DO
         ENDIF

         tcomp = sqrt((3.0*pi)/(32.0*rhozero))
         tff = tcomp * utime
         timeff = gt/tcomp
         valjeans = 0.
         IF (tterm.NE.0) valjeans = abs(tgrav/tterm)
c         WRITE(21,*) timeff, valjeans, nptmass
c
c--reset coordinates to center of mass
c
         IF (icent.EQ.'m') CALL origin (n1, n2, 1.0)
         IF (icent.EQ.'d') THEN

            irhomax = 1
            DO i = 1, npart
               IF (rho(i).GT.rho(irhomax)) irhomax = i
            END DO
            WRITE (*,*) 'irhomax ',irhomax, rho(irhomax), 
     &           xyzmh(1,irhomax), xyzmh(2,irhomax), xyzmh(3,irhomax)
            ncen = 0
            xcen = 0.
            ycen = 0.
            zcen = 0.
            vxcen = 0.
            vycen = 0.
            vzcen = 0.
            DO i = 1, npart
               dx = xyzmh(1,i) - xyzmh(1,irhomax)
               dy = xyzmh(2,i) - xyzmh(2,irhomax)
               dz = xyzmh(3,i) - xyzmh(3,irhomax)
               dvx = vxyzu(1,i) - vxyzu(1,irhomax)
               dvy = vxyzu(2,i) - vxyzu(2,irhomax)
               dvz = vxyzu(3,i) - vxyzu(3,irhomax)
               r2 = dx*dx + dy*dy + dz*dz
               IF (r2.LT.16.0*xyzmh(5,irhomax)**2) THEN
                  ncen = ncen + 1
                  xcen = xcen + xyzmh(1,i)
                  ycen = ycen + xyzmh(2,i)
                  zcen = zcen + xyzmh(3,i)
                  vxcen = vxcen + vxyzu(1,i)
                  vycen = vycen + vxyzu(2,i)
                  vzcen = vzcen + vxyzu(3,i)
               ENDIF
            END DO
            xcen = xcen/ncen
            ycen = ycen/ncen
            zcen = zcen/ncen
            vxcen = vxcen/ncen
            vycen = vycen/ncen
            vzcen = vzcen/ncen
            WRITE (*,*) 'centre: ',ncen, xcen, ycen, zcen
            WRITE (*,*) 'vel centre: ',ncen, vxcen, vycen, vzcen
            DO i = 1, npart
               xyzmh(1,i) = xyzmh(1,i) - xcen
               xyzmh(2,i) = xyzmh(2,i) - ycen
               xyzmh(3,i) = xyzmh(3,i) - zcen
c               vxyzu(1,i) = vxyzu(1,i) - vxcen
c               vxyzu(2,i) = vxyzu(2,i) - vycen
c               vxyzu(3,i) = vxyzu(3,i) - vzcen
            END DO
         ENDIF

         ipt = listpm(1)

c      WRITE(*,*) 'Ptm ',x(ipt),y(ipt),z(ipt),vx(ipt),vy(ipt),vz(ipt)

         DO j = 1, nptmass
            i = listpm(j)
            WRITE (*,*) 'Ptmass ',i,' mass=',xyzmh(4,i), xyzmh(1,i), 
     &           xyzmh(2,i)
         END DO

         tmass1 = 0.0
         tmass2 = 0.0
         tmass3 = 0.0
         tmass4 = 0.0
         tmass5 = 0.0
         dmass = 0.0

         ikount = 0
         totmass = 0.
         totmass_dust = 0.
         distance = 1.0e+20
         denscen = 0.
         DO 500 i = 1, npart
            list(ikount+1) = i
            IF (iphase(i).NE.-1) THEN
               ikount = ikount + 1
               IF (iphase(i).EQ.0) THEN
                  totmass = totmass + xyzmh(4,i)
               ELSEIF (iphase(i).GE.11) THEN
                  totmass_dust = totmass_dust + xyzmh(4,i)
               ENDIF
c            dx = xyzmh(1,i)-xyzmh(1,ipt)
c            dy = xyzmh(2,i)-xyzmh(2,ipt)
c            dz = xyzmh(3,i)-xyzmh(3,ipt)
c            dvx = vxyzu(1,i)-vxyzu(1,ipt)
c            dvy = vxyzu(2,i)-vxyzu(2,ipt)
c            dvz = vxyzu(3,i)-vxyzu(3,ipt)
               dx = xyzmh(1,i)
               dy = xyzmh(2,i)
               dz = xyzmh(3,i)
               dvx = vxyzu(1,i)
               dvy = vxyzu(2,i)
               dvz = vxyzu(3,i)

               radius(i) = SQRT(dx*dx + dy*dy + dz*dz)

               JcrossB_dot_rhat(i) = (JcrossBx(i)*xyzmh(1,i) + 
     &              JcrossBy(i)*xyzmh(2,i) +
     &              JcrossBz(i)*xyzmh(3,i))/radius(i)
               Bmag(i) = SQRT(Bxyz(1,i)**2 +Bxyz(2,i)**2 +Bxyz(3,i)**2)

               rxy(i) = SQRT(dx*dx + dy*dy)
               rvel(i) = (dvx*dx + dvy*dy + dvz*dz)/radius(i)
               rvelt(i) = (dvy*dx - dvx*dy)/rxy(i)
               vsound(i) = SQRT(2./3.*vxyzu(4,i))
               temp(i) = vxyzu(4,i)/ekcle(3,i)

               entropy_per_baryon(i) = 2.5 - 
     &              LOG( rho(i) / 1.67E-24 / 2.0 * 
     &              5.32363E-21 / temp(i)**1.5)

#ifdef NONIDEAL
              IF (encal.EQ.'r') THEN
                  ekclei = getcv(rho(i),vxyzu(4,i))
               ELSE
                  ekclei = 0.0
               ENDIF
               temperature = nimhd_gastemp(encal,vxyzu(4,i),ekclei
     &                      ,dble(rho(i)),gamma,n1,n2,iunique(iorig(i)))
               CALL nicil_get_ion_n(dble(rho(i)),temperature
     &                      ,n_R(:,i),n_electronT(i),ierr)
#endif

c            WRITE (42,*) radius(i)*udist,entropy_per_baryon(i)

               IF (radius(i).LT.0.15 .AND. rvel(i).LT.0.1) THEN
                  tmass1 = tmass1 + xyzmh(4,i)
               ENDIF

            ELSEIF (nptmass.GE.1) THEN
               WRITE(23,88333) xyzmh(1,i),xyzmh(2,i),xyzmh(3,i)
            ENDIF
 500     CONTINUE

         WRITE (*,*) 'ikount=',ikount,' npart=',npart

         IF (iav.EQ.1) THEN
            CALL indexx(ikount, list, radius, indx)
         ELSEIF (iav.EQ.2) THEN
            CALL indexx(ikount, list, rxy, indx)
         ELSE
            WRITE (*,*) 'Error for iav'
            STOP
         ENDIF

         radmax = radmaxstart
         radmin = radmax/radfac
         radcentre = SQRT(radmax*radmin)
         radmean = 0.
         radmean_dust = 0.
         dradius = radcentre
         densitymean = 0.
         densitymean_dust = 0.
         densityvert = 0.
         densityvert_dust = 0.
         radtot = 0.
         velmean = 0.
         vertvelmean = 0.       ! -1.E30
         veltmean = 0.
         velmean_dust = 0.
         vertvelmean_dust = 0.       ! -1.E30
         veltmean_dust = 0.
         vsmean = 0.
         hmean = 0.
         hmean_dust = 0.
         number = 0
         number_dust = 0
         radvertmean = 0.
         radvertmean_dust = 0.
         numbervert = 0
         numbervert_dust = 0
         tempmean = 0.
         tempmeanvert = 0.
         fieldmean = 0.
         fieldmeanvert = 0.
         lorentzmean = 0.
         lorentzmeanvert = 0.
         pressuremean = 0.
         pressureold = 0.
         pressuremeanvert = 0.
         Emean = 0.
         Eold = 0.
         Emeanvert = 0.
         fluxmean = 0.
         fluxmeanvert = 0.
         entropymean = 0.
         entropymeanvert = 0.
         iflg = 0
         xmass = 0.
         xmass_dust = 0.
         shellmass = 0.
         shellmass_dust = 0.
#ifdef NONIDEAL
         nebynmean = 0.
         nzminusmean = 0.
         nzzeromean = 0.
         nzplusmean = 0.
         etaohmmean = 0.
         etahallmean = 0.
         etaambimean = 0.
#endif
c
c--Find quantities with spherical radius or near x-y plane
c
         DO i = ikount, 1, -1
            ipos = list(indx(i))

            IF (iphase(ipos).EQ.0) THEN
               xmass = xmass + xyzmh(4,ipos)
            ELSEIF (iphase(ipos).GE.11) THEN
               xmass_dust = xmass_dust + xyzmh(4,ipos)
            ENDIF

            angle = ABS(180./pi*ASIN(xyzmh(3,ipos)/radius(ipos)))
            IF (iav.EQ.1) THEN
               radiuswant = radius(ipos)
            ELSE
               radiuswant = rxy(ipos)
            ENDIF
            IF (iflg.EQ.0 .AND. i.LE.ikount*radin) THEN
               iflg = 1
               WRITE (*,*) 'Radius ',radin,radiuswant,
     &              radius(indx(list(ikount))),ikount*radin,i
            ENDIF

            IF (radiuswant.LT.radmax .AND. radiuswant.GE.radmin
     &           .AND. (itype.EQ.'s' .OR. itype.EQ.'a')) THEN

               IF (itype.EQ.'s' .OR. 
     &              itype.EQ.'a' .AND. angle.LT.angletake) THEN
                  IF (iphase(ipos).EQ.0) THEN
                     number = number + 1
                     shellmass = shellmass + xyzmh(4,ipos)
                     densitymean = densitymean + rho(ipos)
                     radmean = radmean + radiuswant
                     velmean = velmean + rvel(ipos)
                     veltmean = veltmean + rvelt(ipos)
                     vsmean = vsmean + vsound(ipos)
                     hmean = hmean + xyzmh(5,ipos)
                     densityvar = densityvar + rho(ipos)**2
                     velvar = velvar + rvel(ipos)**2
                     veltvar = veltvar + rvelt(ipos)**2
                     vsvar = vsvar + vsound(ipos)**2
                     hvar = hvar + xyzmh(5,ipos)**2
                     tempmean = tempmean + temp(ipos)
                     fieldmean = fieldmean + Bmag(ipos)
c               lorentzmean = lorentzmean + JcrossB_dot_rhat(ipos)
                     lorentzmean = lorentzmean + JcrossB(ipos)
c               lorentzmean = lorentzmean + JcrossBz(ipos)
                     entropymean = entropymean +entropy_per_baryon(ipos)

c               WRITE (42,*) radiuswant,JcrossB(ipos),lorentzmean

                     pressuremean = pressuremean + 2./3.*vxyzu(4,ipos)*
     &                    uergg*rho(ipos)
                     Emean = Emean + ekcle(1,ipos)*rho(ipos)
                     fluxmean = fluxmean + 
     &                    ekcle(4,ipos)*3.0E+10/ekcle(2,ipos)/rho(ipos)
#ifdef NONIDEAL
                      IF (encal.EQ.'r') THEN
                         ekclei = getcv(rho(ipos),vxyzu(4,ipos))
                      ELSE
                         ekclei = 0.0
                      ENDIF
                      IF (Bmag(ipos).GT.0.0) THEN
                         Bi1 = 1.0/Bmag(ipos)
                      ELSE
                         Bi1 = 0.0
                      ENDIF
                      temperature = nimhd_gastemp(encal,vxyzu(4,ipos)
     &                             ,ekclei,dble(rho(ipos)),gamma,n1,n2
     &                             ,iunique(iorig(ipos)))
                      CALL nicil_get_eta(eta_ohmi,eta_halli
     &                        ,eta_ambii,Bmag(ipos),dble(rho(ipos))
     &                        ,temperature,n_R(:,ipos),n_electronT(ipos)
     &                        ,ierr,data_out)
                      write(666,*) eta_halli,Bmag(ipos),temperature
                      IF (data_out(7) > 0.0) THEN
                         nebynmean = nebynmean + data_out(6)/data_out(7)
                      ENDIF
                      nzminusmean = nzminusmean + data_out(12)
                      nzzeromean  = nzzeromean  + data_out(13)
                      nzplusmean  = nzplusmean  + data_out(14)
                      etaohmmean  = etaohmmean  + eta_ohmi
                      etahallmean = etahallmean + eta_halli
                      etaambimean = etaambimean + eta_ambii
#endif
                  ELSEIF (iphase(ipos).GE.11) THEN
                     number_dust = number_dust + 1
                     shellmass_dust = shellmass_dust + xyzmh(4,ipos)
                     densitymean_dust = densitymean_dust + rho(ipos)
                     radmean_dust = radmean_dust + radiuswant
                     velmean_dust = velmean_dust + rvel(ipos)
                     veltmean_dust = veltmean_dust + rvelt(ipos)
                     hmean_dust = hmean_dust + xyzmh(5,ipos)
                  ENDIF
               ELSEIF (itype.EQ.'a' .AND. angle.GT.90.-2.0*angletake)
     &                 THEN

                  IF (iphase(ipos).EQ.0) THEN
                     numbervert = numbervert + 1
                     radvertmean = radvertmean + radius(ipos)
                     vertvelmean = vertvelmean + rvel(ipos)
c               vertvelmean = MAX(vertvelmean,rvel(ipos))
                     tempmeanvert = tempmeanvert + temp(ipos)
                     fieldmeanvert = fieldmeanvert + Bmag(ipos)
                     lorentzmeanvert = lorentzmeanvert + 
     &                 JcrossBz(ipos)*xyzmh(3,ipos)/ABS(xyzmh(3,ipos))
c     &              JcrossB(ipos)
c     &              JcrossB_dot_rhat(ipos)
c               lorentzmeanvert = MAX(lorentzmeanvert,
c     &              JcrossBz(ipos))
c     &              JcrossB(ipos))
c     &              JcrossB_dot_rhat(ipos))
                     entropymeanvert = entropymeanvert + 
     &                    entropy_per_baryon(ipos)
                     densityvert = densityvert + rho(ipos)
                     pressuremeanvert = pressuremeanvert +
     &                    2./3.*vxyzu(4,ipos)*uergg*rho(ipos)
                     Emeanvert = Emeanvert + ekcle(1,ipos)*rho(ipos)
                     fluxmeanvert = fluxmeanvert + 
     &                    ekcle(4,ipos)*3.0E+10/ekcle(2,ipos)/rho(ipos)
                  ELSEIF (iphase(ipos).GE.11) THEN
                     numbervert_dust = numbervert_dust + 1
                     radvertmean_dust = radvertmean_dust + radius(ipos)
                     vertvelmean_dust = vertvelmean_dust + rvel(ipos)
                     densityvert_dust = densityvert_dust + rho(ipos)
                  ENDIF
               ENDIF
            ELSEIF (radiuswant.LT.radmin) THEN
               IF (number.GT.naveragemin) THEN
                  surfacedensity = shellmass/(pi*(radmax**2-radmin**2))
                  densitymean = densitymean/number
                  radmean = radmean/number
                  velmean = velmean/number
                  veltmean = veltmean/number
                  vsmean = vsmean/number
                  hmean = hmean/number
                  tempmean = tempmean/number
                  fieldmean = fieldmean/number
                  lorentzmean = lorentzmean/number
                  entropymean = entropymean/number
                  pressuremean = pressuremean/number
#ifdef NONIDEAL
                  nebynmean   = nebynmean/number
                  nzminusmean = nzminusmean/number
                  nzzeromean  = nzzeromean/number
                  nzplusmean  = nzplusmean/number
                  etaohmmean  = etaohmmean/number
                  etahallmean = etahallmean/number
                  etaambimean = etaambimean/number
#endif
                  pressuregradient =ABS((pressureold-pressuremean)/
     &                 (dradius*udist))
                  pressureold = pressuremean
                  fluxmean = fluxmean/number

                  Emean = Emean/number
                  Egradient = (Emean - Eold)/(dradius*udist)
                  Eold = Emean

                  fluxmean = - fluxmean*Egradient
!     &                 *4*pi*(radmean*udist)**2

                  densityvar = SQRT(densityvar/(number-1) - 
     &                 densitymean**2)
                  velvar = SQRT(velvar/(number-1) - velmean**2)
                  veltvar = SQRT(veltvar/(number-1) - veltmean**2)
                  vsvar = SQRT(vsvar/(number-1) - vsmean**2)
                  hvar = SQRT(hvar/(number-1) - hmean**2)
                  IF (radmean.NE.0.) THEN
                     omega = veltmean/radmean
                  ELSE
                     omega = 0.
                  ENDIF

                  ndivisor = 50
                  IF (numbervert.GT.naveragemin/ndivisor) THEN
                     radvertmeanout = radvertmean/numbervert
                     vertvelmeanout = vertvelmean  /numbervert
                     tempmeanvertout = tempmeanvert/numbervert
                     fieldmeanvertout = fieldmeanvert/numbervert
                     lorentzmeanvertout = lorentzmeanvert/numbervert
c                  lorentzmeanvertout = lorentzmeanvert
                     entropymeanvertout = entropymeanvert/numbervert
                     pressuremeanvertout = pressuremeanvert/numbervert
                     pressuregradientvert = ABS((pressureoldvert-
     &                    pressuremeanvertout)/(dradius*udist))
                     pressureoldvert = pressuremeanvertout

                     fluxmeanvertout = fluxmeanvert/numbervert

                     Emeanvertout = Emeanvert/numbervert
                     Egradientvertout = (Emeanvertout-Emeanvertold)/
     &                    (dradius*udist)
                     Emeanvertold = Emeanvertout

                     fluxmeanvertout = -fluxmeanvertout*Egradientvertout
!     &                 *4*pi*(radvertmeanout*udist)**2

                     densityvertout = densityvert/numbervert
                  ELSE
                     radvertmeanout = 0.
                     vertvelmeanout = 0. ! -1.E+30
                     tempmeanvertout = 0.
                     fieldmeanvertout = 0.
                     lorentzmeanvertout = 0.
                     entropymeanvertout = 0.
                     pressuremeanvertout = 0.
                     Emeanvertout = 0.
                     Egradientvertout = 0.
                     fluxmeanvertout = 0.
                     densityvertout = 0.
                  ENDIF
c
c--Normalise dust quantities
c
                  IF (number_dust.GT.naveragemin) THEN
                     surfacedensity_dust = 
     &                    shellmass_dust/(pi*(radmax**2-radmin**2))
                     densitymean_dust = densitymean_dust/number_dust
                     radmean_dust = radmean_dust/number_dust
                     velmean_dust = velmean_dust/number_dust
                     veltmean_dust = veltmean_dust/number_dust
                     hmean_dust = hmean_dust/number_dust

                     ndivisor = 50
                     IF (numbervert_dust.GT.naveragemin/ndivisor) THEN
                        radvertmeanout_dust = radvertmean_dust/
     &                       numbervert_dust
                        vertvelmeanout_dust = vertvelmean_dust/
     &                       numbervert_dust
                        densityvertout_dust = densityvert_dust/
     &                       numbervert_dust
                     ELSE
                        radvertmeanout_dust = 0.
                        vertvelmeanout_dust = 0. ! -1.E+30
                        densityvertout_dust = 0.
                     ENDIF
                  ELSE
                     surfacedensity_dust = 0.
                     densitymean_dust = 0.
                     radmean_dust = 0.
                     velmean_dust = 0.
                     veltmean_dust = 0.
                     hmean_dust = 0.
                  ENDIF

                  gravity = 6.67E-08*(totmass-xmass)*umass/
     &                 (radmin*udist)**2 *densitymean

                  WRITE (21,88551) radmean*udist, velmean*uvel, 
     &                 veltmean*uvel, 
     &                 densitymean*udens, hmean, number, 
     &                 velvar*uvel, veltvar*uvel,
     &                 densityvar*udens, hvar, vsmean*uvel, 
     &                 vsvar*uvel, omega,
     &                 radmin*udist,totmass-xmass,surfacedensity,
     &                 tempmean,radvertmeanout*udist,
     &                 vertvelmeanout*uvel,tempmeanvertout,
     &                 densityvertout*udens,
     &                 fieldmean,fieldmeanvertout,pressuremean,
     &                 pressuremeanvertout,
     &                 pressuremean/(fieldmean**2/(2.*4.*pi)),
     &            pressuremeanvertout/(fieldmeanvertout**2/(2.*4.*pi)),
     &                 lorentzmean,lorentzmeanvertout/densityvertout,
     &                 pressuregradient,
     &                 pressuregradientvert/densityvertout,
     &                 gravity/densityvertout,
     &                 entropymean,entropymeanvertout,
     &                 Emean,Emeanvertout,Egradient,Egradientvertout,
     &                 fluxmean, fluxmeanvertout,
     &                 0.5*densitymean*velmean**3,
     &                 0.5*densityvertout*vertvelmeanout**3,
     &                 gravity/densitymean

                  WRITE (22,88515) radmean_dust*udist,velmean_dust*uvel, 
     &                 veltmean_dust*uvel, 
     &                 densitymean_dust*udens, hmean_dust, number_dust, 
     &                 totmass_dust-xmass_dust,
     &                 surfacedensity_dust,
     &                 radvertmeanout_dust*udist,
     &                 vertvelmeanout_dust*uvel,
     &                 densityvertout_dust*udens

#ifdef NONIDEAL
                  WRITE (24,88888) radmean*udist,
     &                 nebynmean,nzminusmean,nzzeromean,nzplusmean,
     &                 etaohmmean*udist*udist/utime,
     &                 etahallmean*udist*udist/utime,
     &                 etaambimean*udist*udist/utime
#endif

                  WRITE (*,*) radmean,omega
                  number = 0.
                  shellmass = 0.
                  densitymean = 0.
                  radmean = 0.
                  velmean = 0.
                  veltmean = 0.
                  vsmean = 0.
                  tempmean = 0.
                  fieldmean = 0.
                  lorentzmean = 0.
                  entropymean = 0.
                  pressuremean = 0.
                  Emean = 0.
                  Egradient = 0.
                  fluxmean = 0.
#ifdef NONIDEAL
                  nebynmean = 0.
                  nzminusmean = 0.
                  nzzeromean = 0.
                  nzplusmean = 0.
                  etaohmmean = 0.
                  etahallmean = 0.
                  etaambimean = 0.
#endif

                  number_dust = 0.
                  shellmass_dust = 0.
                  densitymean_dust = 0.
                  radmean_dust = 0.
                  velmean_dust = 0.
                  veltmean_dust = 0.
                  hmean_dust = 0.

                  IF (numbervert.GT.naveragemin/ndivisor) THEN
                     numbervert = 0
                     radvertmean = 0.
                     vertvelmean = 0. ! -1.E+30
                     tempmeanvert = 0.
                     fieldmeanvert = 0.
                     lorentzmeanvert = 0.
                     entropymeanvert = 0.
                     pressuremeanvert = 0.
                     densityvert = 0.
                     Emeanvert = 0.
                     Egradientvert = 0.
                     fluxmeanvert = 0.
                  ENDIF
                  IF (numbervert_dust.GT.naveragemin/ndivisor) THEN
                     numbervert_dust = 0
                     radvertmean_dust = 0.
                     vertvelmean_dust = 0. ! -1.E+30
                     densityvert_dust = 0.
                  ENDIF
               ENDIF
c
c--Do next radial range
c
               radmax = radmin
               radmin = radmax/radfac
               radcentreold = radcentre
               radcentre = SQRT(radmax*radmin)
               dradius = radcentre - radcentreold
            ENDIF
         END DO

88333    FORMAT(3(1PE14.7,1X))
88334    FORMAT(3(1PE14.7,1X),I6)
88444    FORMAT(4(1PE12.5,1X))
88445    FORMAT(4(1PE12.5,1X),I6)
88515    FORMAT(5(1PE12.5,1X),I7,1X,5(1PE12.5,1X))
88555    FORMAT(5(1PE12.5,1X))
88551    FORMAT(5(1PE12.5,1X),I7,1X,37(1PE12.5,1X))
88556    FORMAT(5(1PE12.5,1X),I6)
88166    FORMAT(I6, 1X, 5(1PE12.5,1X))
88666    FORMAT(6(1PE12.5,1X))
88777    FORMAT(7(1PE12.5,1X))
88888    FORMAT(8(1PE12.5,1X))
88999    FORMAT(9(1PE12.5,1X))
89000    FORMAT(10(1PE12.5,1X))
89111    FORMAT(11(1PE12.5,1X))
89222    FORMAT(12(1PE12.5,1X))
89333    FORMAT(13(1PE12.5,1X))
89444    FORMAT(14(1PE12.5,1X))

c      WRITE (*,*) 'tmass2=',tmass2,' tmass3=',tmass3

         WRITE(15,99080) gt, timeff, fileout
99080    FORMAT(1PE12.5,1X,1PE12.5,1X,A11)
         WRITE(19,88444) gt, timeff, gt*utime, tmass1

c      WRITE (*,*) 'total mass=', totmass
c
         WRITE (*,*) 'dmass = ',dmass
         WRITE (*,*) 'tmass1 = ',tmass1
         WRITE (*,*) 'tmass2 = ',tmass2
         WRITE (*,*) 'tmass3 = ',tmass3
         WRITE (*,*) 'tmass4 = ',tmass4
         WRITE (*,*) 'tmass5 = ',tmass5

         CLOSE(21)
         CLOSE(22)
         IF (isink.EQ.'y' .OR. isink.EQ.'Y') THEN
            CLOSE(23)
         ENDIF
#ifdef NONIDEAL
         CLOSE(24)
#endif
         READ(fileout, 99006) name, ifile
99006    FORMAT(A9,I2)
         ifile = ifile + 1
         IF (ifile.LE.9) THEN
            WRITE(fileout, 99002) name, izero, ifile
         ELSE
            WRITE(fileout, 99006) name, ifile
         ENDIF
         READ(fileout2, 99006) name, ifile
         ifile = ifile + 1
         IF (ifile.LE.9) THEN
            WRITE(fileout2, 99002) name, izero, ifile
         ELSE
            WRITE(fileout2, 99006) name, ifile
         ENDIF
         READ(fileout3, 99006) name, ifile
         ifile = ifile + 1
         IF (ifile.LE.9) THEN
            WRITE(fileout3, 99002) name, izero, ifile
         ELSE
            WRITE(fileout3, 99006) name, ifile
         ENDIF
#ifdef NONIDEAL
         READ(fileout4, 99006) name, ifile
         ifile = ifile + 1
         IF (ifile.LE.9) THEN
            WRITE(fileout4, 99002) name, izero, ifile
         ELSE
            WRITE(fileout4, 99006) name, ifile
         ENDIF
#endif
   20    CONTINUE

         CLOSE(15)

      END DO

      CLOSE (19)

      STOP
      END


      SUBROUTINE quit(i)
      
      STOP
      END


      SUBROUTINE endrun
      CALL quit(0)
      END

