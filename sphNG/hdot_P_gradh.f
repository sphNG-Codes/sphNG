      SUBROUTINE hdot(npart, ipart, dt, xyzmh, dha)
c************************************************************
c                                                           *
c  This subroutine computes the derivative of the smoothing *
c     length.                                               *
c  DJP: no neighbour correction for grad h version          *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/soft'

      DIMENSION xyzmh(5,idim)
      REAL*4 dha(2,idim)
c
c--Compute derivative of h, try to enforce finite range
c     of neighbors
c 
      IF (iphase(ipart).GE.1) THEN
         dha(1,ipart) = 0.
         GOTO 15
      ELSE
         numneigh = nneigh(ipart)
      ENDIF

      IF (xyzmh(5,ipart).LT.hmin .AND. numneigh.GT.neimin) THEN
         dha(1,ipart) = 0.
         GOTO 15
      END IF 

      dhi = xyzmh(5,ipart) * divv(ipart) / rho(ipart) / 3.0
      dha(1,ipart) = dhi

 200  hnewg = xyzmh(5,ipart) + dha(1,ipart)*dt
      IF (hnewg.LE.0) THEN
c      IF (ABS(hnewg-xyzmh(5,ipart)).GE.ABS(xyzmh(5,ipart)/2.0)) THEN
ccc         WRITE(iprint,*) ' h < 0 ', dha(1,ipart), dhi, iorig(ipart),
ccc     &        ipart, numneigh, xyzmh(5,ipart), hnewg, divv(ipart)
         IF (dha(1,ipart).LT.0.)  THEN
            dha(1,ipart) = dha(1,ipart)/2.0
ccc            WRITE(iprint,*) ipart, iorig(ipart), dha(1,ipart)
            GOTO 200
         ENDIF
      ENDIF

ccc        dha(1,ipart) = 0.0

 15   CONTINUE

      RETURN
      END
