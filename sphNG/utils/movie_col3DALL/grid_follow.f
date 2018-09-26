      SUBROUTINE grid(xmin, xmax, ymin, ymax, nformloop, idt, iglobal)
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
      PARAMETER (maxline=3001)
      PARAMETER (mline2=maxline*maxline)
      PARAMETER (npert=30)
      PARAMETER (npert1=npert+1)

      REAL*8 udist, umass, utime, udens, ucolumndens,uergcc
      INTEGER*4 iline
      REAL*4 timeout, xmindensout, xmaxdensout

      PARAMETER (nwidemax = 5000)
      PARAMETER (ntallmax = 3300)
      REAL*4 global
      COMMON /gimage/ global(nwidemax,ntallmax)
      COMMON /gvals / nwide, ntall
      COMMON /gpix  / ntallpix, nwidepix

      CHARACTER*1 idt, iglobal

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/interstellar'

      COMMON /unitsin/ umassi, udisti, utimei, umagfdi

      COMMON /partu/ npp
      COMMON /tempp/ temperaturep(idim)
      COMMON /gridv/ iline
      COMMON /scale/ vmin, vmax
      COMMON /rhonm/ rhozero
      COMMON /phase2/ sinksize, sinkrho, denshigh, ihigh
      COMMON /columntab/ coltable(1000)
      COMMON /obs/ ixyloc(2,idim), iphase2(idim)
c
      INTEGER*8 OMP_lock
      COMMON /locks / OMP_lock(maxline,maxline)
c
      DIMENSION xyzmhr(6,idim),listpm2(10000)
c
      CHARACTER*1 ihigh
      REAL*4 gx(maxline),gy(maxline)
      REAL*4 gxlocal(maxline)
      REAL*4 v1(maxline,maxline), xmaxdens, xmindens
      REAL*4 v2(maxline,maxline), xmaxdens2, xmindens2
c
      DATA pi/3.141592654/
      LOGICAL ifirst
      DATA ifirst/.TRUE./
c
      IF (ifirst) THEN

         print *,'Doing init lock'
         DO i = 1, maxline
            DO j = 1, maxline
               CALL OMP_INIT_LOCK(OMP_lock(j,i))
            END DO
         END DO
         ifirst = .FALSE.
      ENDIF

      IF (iline.GT.maxline) THEN
         WRITE(*,*) 'ERROR - Dimensions not big enough'
         STOP
      ENDIF
      udens = umassi/(udisti**3)
      ucolumndens = udens*udisti
      uergcc = DBLE(umassi)/(DBLE(udisti)*DBLE(utimei)**2)
c
c--compute steps and narrow the number of particles
c
      x0 = xmin
      y0 = ymin
      stepx = (xmax - xmin) / (iline - 1)
      stepy = (ymax - ymin) / (iline - 1)
      npp = 0
      nlistpm = nptmass
      nlistpmadd = 0
      DO 12 i = 1, npart
         IF (iphase(i).GE.0) THEN
            IF (xyzmh(1,i)+2.0*xyzmh(5,i).LT.xmin.OR.
     &           xyzmh(1,i)-2.0*xyzmh(5,i).GT.xmax) GOTO 12
            IF (xyzmh(2,i)+2.0*xyzmh(5,i).LT.ymin.OR.
     &           xyzmh(2,i)-2.0*xyzmh(5,i).GT.ymax) GOTO 12
            npp = npp + 1
            xyzmhr(1,npp) = xyzmh(1,i)
            xyzmhr(2,npp) = xyzmh(2,i)
            xyzmhr(3,npp) = xyzmh(3,i)
            xyzmhr(4,npp) = xyzmh(4,i)
            xyzmhr(6,npp) = rho(i)

            IF (iphase(i).GE.1) THEN
               DO j = 1, nptmass
                  IF (listpm(j).EQ.i) THEN
                     listpm2(j)=npp
                     GOTO 19
                  ENDIF
               END DO
            ENDIF
 19          CONTINUE




            IF (idt.EQ.'r' .OR. idt.EQ.'w') THEN
               temperaturep(npp) = (ekcle(1,i)*rho(i)/
     &              (7.5646e-15/uergcc))**0.25
c
c--Gives not-mass weighted
c
              IF (idt.EQ.'r') temperaturep(npp)=temperaturep(npp)/rho(i)
            ELSEIF (idt.EQ.'o') THEN
               temperaturep(npp) = ekcle(2,i)
            ELSE
               temperaturep(npp) = vxyzu(4,i)/ekcle(3,i)
            ENDIF
            IF (iphase(i).GE.1) THEN
               DO j = 1, nptmass
                  IF (listpm(j).EQ.i) THEN
                     listpm2(j)=npp
                     GOTO 9
                  ENDIF
               END DO
            ENDIF
 9          CONTINUE

c            write (*,*) i,iphase(i),rho(i),udens,denshigh

            IF (iphase(i).EQ.0 .AND. .NOT.((ihigh.EQ.'Y' .OR. 
     &           ihigh.EQ.'y') .AND. rho(npp)*udens.GT.denshigh)) THEN
               xyzmhr(5,npp) = xyzmh(5,i)
            ELSE
               IF (iphase(i).EQ.0) THEN
c
c--Convert high density gas to look like a sink if not within 10 AU of sink
c
                  DO j = 1, nlistpm
                     r2 = (xyzmh(1,i)-xyzmh(1,listpm(j)))**2 +
     &                    (xyzmh(2,i)-xyzmh(2,listpm(j)))**2 + 
     &                    (xyzmh(3,i)-xyzmh(3,listpm(j)))**2
c                     write (*,*) '   ',j,listpm(j),sqrt(r2)
                     IF (r2.LT.(10.0*1.5E+13/udist)**2) GOTO 10
                  END DO
c                  write (*,*) i,iphase(i),rho(i),udens,denshigh
                  iphase(i) = 2
                  nlistpmadd = nlistpmadd + 1
                  nlistpm = nlistpm + 1
                  listpm2(nlistpm) = npp
                  GOTO 11

 10               xyzmhr(5,npp) = xyzmh(5,i)
                  GOTO 12
               ENDIF
 11            xyzmhr(5,npp) = sinksize*(xmax-xmin)
               xyzmhr(4,npp) = sinkrho*3.142*(2.0*sinksize*
     &              (xmax-xmin))**2
            ENDIF
         ENDIF
 12   CONTINUE

      write(*,*)'Number of particles kept   :  ',npp
      write(*,*)'Number of sinks (real, new):  ',nptmass,nlistpm
      write(*,*)'Num gas convert to sink    :  ',nlistpmadd

      xmaxdens = -1.0
      xmindens = 1.0E+20
      xmaxdens2 = -1.0
      xmindens2 = 1.0E+20
c
c--initialize grid
c     
C$OMP PARALLEL default(none)
C$OMP& shared(iline,v1,v2,npp,xyzmhr,stepx,stepy,x0,y0)
C$OMP& shared(ucolumndens,temperaturep,idt,gx,gy,coltable,OMP_lock)
C$OMP& shared(ixyloc)
C$OMP& private(i,j,xi,yi,zi,pmassi,rhoi,hi,hproj,ix,iy,icx,icy)
C$OMP& private(idepx,idepy,ifinx,ifiny,kx,ky,val,valuei,hi12)
C$OMP& private(gxlocal,gylocal,vr2,ipos1)
C$OMP& reduction(MIN:xmindens)
C$OMP& reduction(MIN:xmindens2)
C$OMP& reduction(MAX:xmaxdens)
C$OMP& reduction(MAX:xmaxdens2)
C$OMP DO SCHEDULE(static)
      DO i = 1, iline
         gx(i) = x0 + (i - 1) * stepx
      END DO
C$OMP END DO
C$OMP DO SCHEDULE(static)
      DO j = 1, iline
         gy(j) = y0 + (j - 1) * stepy
      END DO
C$OMP END DO
C$OMP DO SCHEDULE(dynamic, 10)
      DO i = 1, iline
         DO j = 1, iline
            v2(i,j) = 0.0
            v1(i,j) = 0.0
         END DO
      END DO
C$OMP END DO
c
c--compute values on grid
c
C$OMP DO SCHEDULE(dynamic, 10)
      DO 30 i = 1, npp
c         IF (MOD(i,100).EQ.0) WRITE (*,*) i,xyzmhr(1,i),
c     &        xyzmhr(4,i),xyzmhr(5,i)
         xi = xyzmhr(1,i)
         yi = xyzmhr(2,i)
         zi = xyzmhr(3,i)
         valuei = xyzmhr(4,i)
         rhoi = xyzmhr(6,i)
         hi = xyzmhr(5,i)
         hproj = 2.0*hi
         hi12 = 1.0/hi**2
c
c--find index of closest grid point
c
         ix = nint( (xi - x0) / stepx) + 1
         iy = nint( (yi - y0) / stepy) + 1
         ixyloc(1,i) = ix
         ixyloc(2,i) = iy
         icx = hproj / stepx + 2
         icy = hproj / stepy + 2
         idepx = max0(1, ix - icx)
         idepy = max0(1, iy - icy)
         ifinx = min0(iline, ix + icx)
         ifiny = min0(iline, iy + icy)
c
c--compute particle's contribution to all grid points
c
         DO kx = idepx, ifinx
            gxlocal(kx) = (xi - gx(kx))**2 * hi12
         END DO
         DO 20 ky = idepy, ifiny
c            gy = y0 + (ky - 1) * stepy
            gylocal = (yi - gy(ky))**2 * hi12
            DO 20 kx = idepx, ifinx
c               gx = x0 + (kx - 1) * stepx
c               CALL value (xi, yi, zi, hi, gx, gy, val, valuei)

               vr2 = gxlocal(kx) + gylocal
               IF (vr2.GE.4.0) THEN
                  val = 0.
               ELSE
                  ipos1=MIN(1000, INT(vr2*250.0)+1)
                  val = coltable(ipos1)*valuei*hi12
               ENDIF
cC$OMP ATOMIC
               CALL OMP_SET_LOCK(OMP_lock(kx,ky))
               
               v1(kx,ky) = v1(kx,ky) + val
cC$OMP ATOMIC
               v2(kx,ky) = v2(kx,ky) + temperaturep(i)*val

               CALL OMP_UNSET_LOCK(OMP_lock(kx,ky))
 20      CONTINUE
 30   CONTINUE
C$OMP END DO
c
c--create output file
c
C$OMP DO SCHEDULE(dynamic, 10)
      DO j=1, iline
         DO i=1, iline
            IF (idt.NE.'r'.AND.idt.NE.'o') v2(i,j)=v2(i,j)/v1(i,j)
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

      IF (iglobal.EQ.'y') THEN
c
c--Move image to global picture
c
         nwidestart = MOD((nformloop-1),nwide)*(iline+1)
         ntallstart = ntallpix - ((nformloop-1)/nwide+1)*(iline+1)
         print *,nwidestart,ntallstart,nwide,ntall,iline
         DO i = 1, iline
            DO j = 1, iline
               global(nwidestart + i, ntallstart + j) = v1(i,j)
            END DO
         END DO
      ELSE
c
c--Write individual images
c
         IF (idt.NE.'d') THEN
            timeout = time
            xmaxdensout = xmaxdens2
            xmindensout = xmindens2
            write (*,*) iline,iline,timeout,xmaxdensout,xmindensout
      WRITE(16,IOSTAT=io) iline,iline,timeout,xmaxdensout,xmindensout
      WRITE(16,IOSTAT=io) ((v2(i,j), i=1,iline), j=1,iline)
         ELSE
            timeout = time
            xmaxdensout = xmaxdens
            xmindensout = xmindens
            write (*,*) iline,iline,timeout,xmaxdensout,xmindensout
      WRITE(16,IOSTAT=io) iline,iline,timeout,xmaxdensout,xmindensout
      WRITE(16,IOSTAT=io) ((v1(i,j), i=1,iline), j=1,iline)
         ENDIF
      ENDIF

      OPEN (22,FILE='contourD.dat')
      DO i = 1, iline
         DO j = 1, iline
            WRITE (22,998) v1(i,j)
 998        FORMAT(1PE12.5)
         END DO
      END DO
      CLOSE(22)
         
      OPEN (22,FILE='contourT.dat')
      DO i = 1, iline
         DO j = 1, iline
            WRITE (22,997) v2(i,j)
 997        FORMAT(1PE12.5)
         END DO
      END DO
      CLOSE(22)

      OPEN (22,FILE='stars.dat')
      DO i = 1, nptmass
         WRITE (22,*) ixyloc(1,listpm2(i)),ixyloc(2,listpm2(i))
      END DO
      CLOSE(22)

 999  RETURN
      END
