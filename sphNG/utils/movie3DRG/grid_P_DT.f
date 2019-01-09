      SUBROUTINE grid (zobs, wangle, ifile1, ifile2, xeye, idt)
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
      INCLUDE 'idim'
c      PARAMETER (maxline=27010)
      PARAMETER (maxline=5121)
      PARAMETER (npert=30)
      PARAMETER (npert1=npert+1)

      REAL*8 ucolumndens
      INTEGER*4 ilinex,iliney
      REAL*4 timeout, xmindensout, xmaxdensout
      CHARACTER*1 idt

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/interstellar'

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
c
c--initialize grid
c     
C$OMP PARALLEL default(none)
C$OMP& shared(ilinex,v1,v2,xyzmh,rho,stepx,stepy,x0,y0,zobs)
C$OMP& shared(ucolumndens,idt,vxyzu,ekcle)
C$OMP& shared(npart,iphase,xeye)
C$OMP& shared(sinksize,sinkrho,wangle,iliney,coltable,gx,gy,hmin)
C$OMP& shared(chemistry)
C$OMP& private(i,j,xi,yi,zi,pmassi,hi,hproj,ix,iy,icx,icy)
C$OMP& private(idepx,idepy,ifinx,ifiny,kx,ky,val,valuei)
C$OMP& private(temperaturei,zrel,vr,vr2,ipos1,hi12,gxlocal,gylocal)
C$OMP& private(xz,yz,xmid,ymid,xsize,ysize)
C$OMP& reduction(MIN:xmindens)
C$OMP& reduction(MIN:xmindens2)
C$OMP& reduction(MAX:xmaxdens)
C$OMP& reduction(MAX:xmaxdens2)
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
            v2(i,j) = 0.0
 15         v1(i,j) = 0.0
C$OMP END DO
c
c--compute values on grid
c
C$OMP DO SCHEDULE(dynamic, 1000)
      DO 30 i = 1, npart
         IF (MOD(i,100000).EQ.0) WRITE (*,*) i
         IF (iphase(i).GE.0) THEN
            IF (iphase(i).EQ.0) THEN
               hi = MAX(xyzmh(5,i),hmin)
               temperaturei = vxyzu(4,i)/ekcle(3,i)
c               temperaturei = chemistry(3,i)
            ELSEIF (iphase(i).GT.0) THEN
               xyzmh(5,i) = sinksize
               hi = xyzmh(5,i)
               xyzmh(4,i) = sinkrho*3.142*(2.0*sinksize)**2
               temperaturei = 10000.0
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
                  ENDIF
C$OMP ATOMIC
                  v1(kx,ky) = v1(kx,ky) + val
C$OMP ATOMIC
                  v2(kx,ky) = v2(kx,ky) + temperaturei*val

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
            IF (v1(i,j).GT.0.) v2(i,j)=v2(i,j)/v1(i,j)
            v1(i,j)=v1(i,j)*ucolumndens
            xmaxdens = MAX(xmaxdens, v1(i,j))
            xmindens = MIN(xmindens, v1(i,j))
            xmaxdens2 = MAX(xmaxdens2, v2(i,j))
            xmindens2 = MIN(xmindens2, v2(i,j))
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL
      WRITE(*,*) 'Max col density (g cm^-2)=',xmaxdens
      WRITE(*,*) 'Min col density (g cm^-2)=',xmindens
      WRITE(*,*) 'Max temperature (K)      =',xmaxdens2
      WRITE(*,*) 'Min temperature (K)      =',xmindens2

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

      IF (idt.EQ.'t' .OR. xeye.EQ.0.) THEN
         timeout = time
         xmaxdensout = xmaxdens2
         xmindensout = xmindens2
         WRITE(ifile2out,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(ifile2out,IOSTAT=io) ((v2(i,j), i=1,ilinex), j=1,iliney)
      ENDIF
      IF (idt.EQ.'d' .OR. xeye.EQ.0.) THEN
         timeout = time
         xmaxdensout = xmaxdens
         xmindensout = xmindens
         WRITE(ifile1out,IOSTAT=io) ilinex,iliney,timeout,xmaxdensout,
     &        xmindensout
         WRITE(ifile1out,IOSTAT=io) ((v1(i,j), i=1,ilinex), j=1,iliney)
      ENDIF

      RETURN
      END
