      SUBROUTINE grid (zobs,wangle,ifile1,ifile2,xeye,idt,contour1)   
c
c*****************************************************************
c                                                                *
c  This subroutines interpolates quantitites on a given grid.    *
c  The quantity plotted is column density through the simulation *
c     calculated by integrating density through the SPH          *
c     smoothing kernel for all particles contributing to the     *
c     line of sight density.                                     *
c                                                                *
c*****************************************************************
c
      use ieee_arithmetic

      INCLUDE 'idim'
c      PARAMETER (maxline=27010)
      PARAMETER (maxline=5121)
      PARAMETER (npert=30)
      PARAMETER (npert1=npert+1)

      REAL*8 ucolumndens
      INTEGER*4 ilinex,iliney
      REAL*4 timeout, xmindensout, xmaxdensout
      CHARACTER*1 idt

      CHARACTER*13 contour1,filename

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/raddust'

      REAL xmetallicity
      COMMON /chemmetal/ xmetallicity

      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      COMMON /aux  / hmin
      COMMON /gridv/ ilinex,iliney
      COMMON /origi/ cmx, cmy, cmz
      COMMON /scale/ vmin, vmax
      COMMON /phase2/ sinksize, sinkrho, denshigh, ihigh
      COMMON /numfile/ nfile
      COMMON /columntab/ coltable(1000)
c
      CHARACTER*1 ihigh
      REAL*4 v1(maxline,maxline), xmaxdens, xmindens
      REAL*4 gx(maxline),gy(maxline)
      REAL*4 gylocal(maxline)
      REAL*4 v2(maxline,maxline), xmaxdens2, xmindens2
c
c--Arrays to store images for dustRT (radiation temperature, chemistry)
c
      REAL*4 radT(maxline*idustRT,maxline*idustRT)
      REAL*4 radD(maxline*idustRT,maxline*idustRT)
      REAL*4 Cplus(maxline*idustRT,maxline*idustRT)
      REAL*4 Catomic(maxline*idustRT,maxline*idustRT)
      REAL*4 CO(maxline*idustRT,maxline*idustRT)
      REAL*4 Hatomic(maxline*idustRT,maxline*idustRT)
      REAL*4 Hmolecular(maxline*idustRT,maxline*idustRT)

      DATA pi/3.141592654/
c
      IF (ilinex.GT.maxline) THEN
         WRITE(*,*) 'ERROR - Dimensions not big enough'
         STOP
      ENDIF
c
c--Compute steps and narrow the number of particles
c
      stepx = (2.0*wangle) / (ilinex)
      stepy = stepx

      udens = umassi/(udisti**3)
      ucolumndens = udens*udisti
      uerg = DBLE(umassi)*DBLE(udisti)**2/DBLE(utimei)**2
      uradtemp = DBLE(uerg)/(DBLE(udisti)**3)*
     &     2.99792458E+10/(4.0*5.6704E-05)

      xmaxdens = -1.0
      xmindens = 1.0E+20
      xmaxdens2 = -1.0
      xmindens2 = 1.0E+20
c
c--Convert high density gas to look like a sink if not within 10 AU of sink
c
      nlistpm = nptmass
      nlistpmadd = 0

      IF (ihigh.EQ.'Y' .OR. ihigh.EQ.'y') THEN
         DO i = 1, npart
            IF (iphase(i).EQ.0 .AND. rho(i)*udens.GT.denshigh) THEN
               DO j = 1, nlistpm
                  r2 = (xyzmh(1,i)-xyzmh(1,listpm(j)))**2 + 
     &                 (xyzmh(2,i)-xyzmh(2,listpm(j)))**2 +
     &                 (xyzmh(3,i)-xyzmh(3,listpm(j)))**2
                  IF (r2.LT.(10.0*1.5E+13/udisti)**2) GOTO 12
               END DO
               iphase(i) = 2
               nlistpmadd = nlistpmadd + 1
               nlistpm = nlistpm + 1
               listpm(nlistpm) = i
            ENDIF
 12         CONTINUE
         END DO
      ENDIF
      write(*,*)'Number of sinks (true, tot):  ',nptmass,nlistpm
      write(*,*)'Num gas convert to sink    :  ',nlistpmadd

      IF (idustRT.GT.0) print *,'Metallicity ',xmetallicity
c
c--initialize grid
c     
C$OMP PARALLEL default(none)
C$OMP& shared(ilinex,v1,v2,xyzmh,rho,stepx,stepy,x0,y0,zobs)
C$OMP& shared(ucolumndens,idt,vxyzu,ekcle)
C$OMP& shared(npart,iphase,xeye)
C$OMP& shared(sinksize,sinkrho,wangle,iliney,coltable,gx,gy,hmin)
C$OMP& shared(chemistry,h2frac,uradtemp,dust_tk,xmetallicity)
C$OMP& shared(radT,radD,Cplus,Catomic,CO,Hatomic,Hmolecular)
C$OMP& private(i,j,xi,yi,zi,pmassi,hi,hproj,ix,iy,icx,icy)
C$OMP& private(idepx,idepy,ifinx,ifiny,kx,ky,val,valuei)
C$OMP& private(temperaturei,zrel,vr,vr2,ipos1,hi12,gxlocal,gylocal)
C$OMP& private(xz,yz,xmid,ymid,xsize,ysize)
C$OMP& private(radTi,radDi,Cplusi,Catomici,COi,Hatomici,Hmoleculari)
C$OMP& reduction(MIN:xmindens,xmindens2,xminradT,xminCplus,xminCatomic)
C$OMP& reduction(MIN:xminCO,xminHatomic,xminHmol,xminradD)
C$OMP& reduction(MAX:xmaxdens,xmaxdens2,xmaxradT,xmaxCplus,xmaxCatomic)
C$OMP& reduction(MAX:xmaxCO,xmaxHatomic,xmaxHmol,xmaxradD)
C$OMP DO SCHEDULE(static)
      DO i = 1, ilinex
         gx(i) = (i - 1 -ilinex/2) * stepx
      END DO
C$OMP END DO
C$OMP DO SCHEDULE(static)
      DO j = 1, iliney
         gy(j) = (j - 1 -iliney/2) * stepy
      END DO
C$OMP END DO
C$OMP DO SCHEDULE(dynamic, 10)
      DO 15 i = 1, ilinex
         DO 15 j = 1, iliney
            IF (idustRT.EQ.1) THEN
               radT(i,j) = 0.0
               radD(i,j) = 0.0
               Cplus(i,j) = 0.0
               Catomic(i,j) = 0.0
               CO(i,j) = 0.0
               Hatomic(i,j) = 0.0
               Hmolecular(i,j) = 0.0
            ENDIF
            v2(i,j) = 0.0
 15         v1(i,j) = 0.0
C$OMP END DO
c
c--compute values on grid
c
C$OMP DO SCHEDULE(dynamic, 1000)
      DO 30 i = 1, npart
         IF (MOD(i,1000000).EQ.0) WRITE (*,*) i
         IF (iphase(i).GE.0) THEN
            IF (iphase(i).EQ.0) THEN
               hi = MAX(xyzmh(5,i),hmin)
               temperaturei = vxyzu(4,i)/ekcle(3,i)
               IF (idustRT.EQ.1) THEN
                  radTi = (ekcle(1,i)*rho(i)*uradtemp)**0.25
                  radDi = dust_tk(1,i)
                  Cplusi = chemistry(1,i)/xmetallicity
                  IF (ieee_is_nan(Cplusi)) Cplusi=0.
                  Catomici = chemistry(2,i)/xmetallicity
                  IF (ieee_is_nan(Catomici)) Catomici=0.
                  COi = chemistry(3,i)/xmetallicity
                  IF (ieee_is_nan(COi)) COi=0.
                  Hatomici = 1.0 - 2.0*h2frac(i)
                  Hmoleculari = 2.0*h2frac(i)
               ENDIF
            ELSEIF (iphase(i).GT.0) THEN
               xyzmh(5,i) = sinksize
               hi = xyzmh(5,i)
               xyzmh(4,i) = sinkrho*3.142*(2.0*sinksize)**2
               temperaturei = 100000.0
               IF (idustRT.EQ.1) THEN
                  radTi = 100000.0
                  radDi = 100000.0
                  Cplusi = 1000.0
                  Catomici = 1000.0
                  COi = 1000.0
                  Hatomici = 1000.0
                  Hmoleculari = 1000.0
               ENDIF
            ENDIF
            xi = xyzmh(1,i)
            yi = xyzmh(2,i)
            zi = xyzmh(3,i)
            zrel = zi + zobs
            valuei = xyzmh(4,i)
            hproj = 2.0*hi
c
c--find index of closest grid point
c
            IF (zrel.LT.0.0 .OR. (ABS(xi-xeye)-hproj)/zrel.GT.wangle.OR.
     &           (ABS(yi)-hproj)/zrel.GT.wangle) GOTO 30

            hi12 = 1.0/hi**2

            xz = 1.0/(zrel*stepx)
            yz = 1.0/(zrel*stepy)

            xmid = (xi - xeye ) * xz
            ymid = yi * yz
            xsize = hproj * xz
            ysize = hproj * yz
            idepx = max0(1, NINT(xmid - xsize) + 1 + ilinex/2)
            idepy = max0(1, NINT(ymid - ysize) + 1 + iliney/2)
            ifinx = min0(ilinex, NINT(xmid + xsize) + 1 + ilinex/2)
            ifiny = min0(iliney, NINT(ymid + ysize) + 1 + iliney/2)
c
c--compute particle's contribution to all grid points
c
            DO ky = idepy,ifiny
               gylocal(ky) = (yi - gy(ky)* zrel)**2 * hi12
            END DO

            DO 20 kx = idepx,ifinx
               gxlocal = (xi - (gx(kx)* zrel + xeye))**2 * hi12

               DO 21 ky = idepy, ifiny
                  vr2 = gxlocal + gylocal(ky)

                  IF (vr2.GE.4.0) THEN
                     val = 0.
                  ELSE
                     ipos1=MIN(1000, INT(vr2*250.0)+1)
                     val = coltable(ipos1)*valuei*hi12
C$OMP ATOMIC
                     v1(kx,ky) = v1(kx,ky) + val
C$OMP ATOMIC
                     v2(kx,ky) = v2(kx,ky) + temperaturei*val

                     IF (idustRT.EQ.1) THEN
C$OMP ATOMIC
                        radT(kx,ky) = radT(kx,ky) + radTi*val
C$OMP ATOMIC
                        radD(kx,ky) = radD(kx,ky) + radDi*val
C$OMP ATOMIC
                        Cplus(kx,ky) = Cplus(kx,ky) + Cplusi*val
C$OMP ATOMIC
                        Catomic(kx,ky) = Catomic(kx,ky) + Catomici*val
C$OMP ATOMIC
                        CO(kx,ky) = CO(kx,ky) + COi*val

c                        IF (kx.GT.418 .AND. kx.LT.423 .AND. 
c     &                       ky.GT.205 .AND. ky.LT.211) 
c     &                       print *,'Cp ',kx,ky,CO(kx,ky),COi,val
C$OMP ATOMIC
                        Hatomic(kx,ky) = Hatomic(kx,ky) + Hatomici*val
C$OMP ATOMIC
                        Hmolecular(kx,ky) = Hmolecular(kx,ky) + 
     &                       Hmoleculari*val
                     ENDIF
                  ENDIF

 21            CONTINUE
 20         CONTINUE
         ENDIF
 30   CONTINUE
C$OMP END DO
c
c--create output file
c
C$OMP DO SCHEDULE(dynamic, 10)
      DO j=1, iliney
         DO i=1, ilinex

c            IF (i.GT.26 .AND. i.LT.30 .OR. j.GT.26 .AND. i.LT.30)
c     &           print *,v1(i,j),v2(i,j)

            IF (v1(i,j).GT.0.) THEN
               v2(i,j) = v2(i,j)/v1(i,j)
               IF (idustRT.EQ.1) THEN
                  radT(i,j) = radT(i,j)/v1(i,j)
                  radD(i,j) = radD(i,j)/v1(i,j)
                  Cplus(i,j) = Cplus(i,j)/v1(i,j)
                  Catomic(i,j) = Catomic(i,j)/v1(i,j)
c               IF (i.GT.418 .AND. i.LT.423 .AND. j.GT.205 .AND.
c     &              j.LT.211) print *,'COp ',i,j,CO(i,j),v1(i,j)

                  CO(i,j) = CO(i,j)/v1(i,j)

c               IF (i.GT.418 .AND. i.LT.423 .AND. j.GT.205 .AND.
c     &              j.LT.211) print *,'COn ',i,j,CO(i,j),v1(i,j)

                  Hatomic(i,j) = Hatomic(i,j)/v1(i,j)
                  Hmolecular(i,j) = Hmolecular(i,j)/v1(i,j)
               ENDIF
            ENDIF
            xmaxdens2 = MAX(xmaxdens2, v2(i,j))
            xmindens2 = MIN(xmindens2, v2(i,j))

            IF (idustRT.EQ.1) THEN
               xmaxradT = MAX(xmaxradT, radT(i,j))
               xminradT = MIN(xminradT, radT(i,j))
               xmaxradD = MAX(xmaxradD, radD(i,j))
               xminradD = MIN(xminradD, radD(i,j))
               xmaxCplus = MAX(xmaxCplus, Cplus(i,j))
               xminCplus = MIN(xminCplus, Cplus(i,j))
               xmaxCatomic = MAX(xmaxCatomic, Catomic(i,j))
               xminCatomic = MIN(xminCatomic, Catomic(i,j))
               xmaxCO = MAX(xmaxCO, CO(i,j))
               xminCO = MIN(xminCO, CO(i,j))
               xmaxHatomic = MAX(xmaxHatomic, Hatomic(i,j))
               xminHatomic = MIN(xminHatomic, Hatomic(i,j))
               xmaxHmol = MAX(xmaxHmol, Hmolecular(i,j))
               xminHmol = MIN(xminHmol, Hmolecular(i,j))
            ENDIF
            v1(i,j)=v1(i,j)*ucolumndens
            xmaxdens = MAX(xmaxdens, v1(i,j))
            xmindens = MIN(xmindens, v1(i,j))
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL
      WRITE(*,*) 'Max col density (g cm^-2)=',xmaxdens
      WRITE(*,*) 'Min col density (g cm^-2)=',xmindens
      WRITE(*,*) 'Max temperature (K)      =',xmaxdens2
      WRITE(*,*) 'Min temperature (K)      =',xmindens2

      IF (idustRT.EQ.1) THEN
         WRITE (*,*) 'Max and min Rad T = ',xmaxradT,xminradT
         WRITE (*,*) 'Max and min Rad D = ',xmaxradD,xminradD
         WRITE (*,*) 'Max and min Cplus = ',xmaxCplus,xminCplus
         WRITE (*,*) 'Max and min Catomic = ',xmaxCatomic,xminCatomic
         WRITE (*,*) 'Max and min CO = ',xmaxCO,xminCO
         WRITE (*,*) 'Max and min Hatomic = ',xmaxHatomic,xminHatomic
         WRITE (*,*) 'Max and min Hmolecular = ',xmaxHmol,xminHmol
      ENDIF

c      v1(1,1)=10000.0

c
c--If xeye=0 then output both density and temperature files, but if it
c     is non-zero, then output EITHER density or temperature to different
c     files.
c
      ifile1out = ifile1
      ifile2out = ifile2
      IF (xeye.GT.0.) THEN
         ifile2out = ifile1
      ELSEIF (xeye.LT.0.) THEN
         ifile1out = ifile2
      ENDIF

      timeout = time
      IF (idt.EQ.'t' .OR. xeye.EQ.0.) THEN
         xmaxdensout = xmaxdens2
         xmindensout = xmindens2
         WRITE(ifile2out,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(ifile2out,IOSTAT=io) (v2(1:ilinex,i), i=1,iliney)
c         DO i = 1, iliney
c            DO j = 1, ilinex
c               IF (i.GT.418 .AND. i.LT.423 .AND. j.GT.205 .AND.
c     &              j.LT.211) print *,'T ',i,j,v2(i,j)
c            END DO
c         END DO
      ENDIF
      IF (idt.EQ.'d' .OR. xeye.EQ.0.) THEN
         xmaxdensout = xmaxdens
         xmindensout = xmindens
         WRITE(ifile1out,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(ifile1out,IOSTAT=io) (v1(1:ilinex,i), i=1,iliney)
      ENDIF

      IF (idustRT.EQ.1) THEN
c
c--Radiation temperature
c
         filename = contour1(1:2)//'R'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
         timeout = time
         xmaxdensout = xmaxradT
         xmindensout = xminradT
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (radT(1:ilinex,i), i=1,iliney)
         CLOSE(18)
c
c--Dust temperature
c
         filename = contour1(1:2)//'U'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
         timeout = time
         xmaxdensout = xmaxradD
         xmindensout = xminradD
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (radD(1:ilinex,i), i=1,iliney)
         CLOSE(18)
c
c--Cplus
c
         filename = contour1(1:2)//'P'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
c         OPEN (18,file=filename)
         timeout = time
         xmaxdensout = xmaxCplus
         xmindensout = xminCplus
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (Cplus(1:ilinex,i), i=1,iliney)
         CLOSE(18)
c
c--Catomic
c
         filename = contour1(1:2)//'A'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
         timeout = time
         xmaxdensout = xmaxCatomic
         xmindensout = xminCatomic
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (Catomic(1:ilinex,i), i=1,iliney)
         CLOSE(18)
c
c--CO
c
         filename = contour1(1:2)//'O'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
         timeout = time
         xmaxdensout = xmaxCO
         xmindensout = xminCO
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (CO(1:ilinex,i), i=1,iliney)
         DO i = 1, iliney
            DO j = 1, ilinex
               IF (i.GT.418 .AND. i.LT.423 .AND. j.GT.205 .AND.
     &              j.LT.211) print *,'Cp ',i,j,CO(i,j)
            END DO
         END DO
         CLOSE(18)
c
c--Hatomic
c
         filename = contour1(1:2)//'I'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
         timeout = time
         xmaxdensout = xmaxHatomic
         xmindensout = xminHatomic
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (Hatomic(1:ilinex,i), i=1,iliney)
         CLOSE(18)
c
c--Hmolecular
c
         filename = contour1(1:2)//'M'//contour1(4:13)
         print *,'Doing ',filename
         OPEN (18,file=filename, FORM = 'unformatted')
         timeout = time
         xmaxdensout = xmaxHmol
         xmindensout = xminHmol
         WRITE(18,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(18,IOSTAT=io) (Hmolecular(1:ilinex,i), i=1,iliney)
         CLOSE(18)
      ENDIF

      RETURN
      END
