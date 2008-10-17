      SUBROUTINE hcalc
c************************************************************
c                                                           *
c  This subroutine computes the derivative of the smoothing *
c     length.                                               *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/Bxyz'
      
      CHARACTER*4 varmhdtemp
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine hcalc')

      third = 1./3.
      neimin = 30
      neimax = 70
c      neimin = 80
c      neimax = 120
      acc = 0.9
      dxbound = xmax - xmin
      dybound = ymax - ymin
      dzbound = zmax - zmin

      neimean = (neimax + neimin) / 2
      neirange = (neimax - neimin) / 4
      neisup = neimean + neirange
      neiinf = neimean - neirange

      ikount = 0

      WRITE(*,*) neimax, neimean, neimin, neisup, neiinf, neirange

      nlst = npart
      DO i = 1, npart
         iscurrent(i) = .TRUE.
         llist(i) = i
         isteps(i) = 0
         it0(i) = 0
         it1(i) = 0
         it2(i) = 0
      END DO

      CALL ghostp(ntot, npart, xyzmh, vxyzu, ekcle, Bevolxyz)

      IF (igrape.EQ.0) THEN
         DO i = 1, ntot
            DO j = 1, 5
               dumxyzmh(j,i) = xyzmh(j,i)
            END DO
         END DO
         DO i = 1, ntot
            DO j = 1, 4
               dumvxyzu(j,i) = vxyzu(j,i)
            END DO
         END DO
         DO i = 1, ntot
            DO j = 1, isizealphaMM
               dumalpha(j,i) = alphaMM(j,i)
            END DO
         END DO
         DO i = 1, ntot
            DO j = 1, 3
               dumBevolxyz(j,i) = Bevolxyz(j,i)
            END DO
         END DO
         DO i = 1, ntot
            DO j = 1, 5
               dumekcle(j,i) = ekcle(j,i)
            END DO
         END DO
         WRITE(*,*) ' Making tree'
         CALL insulate(1, ntot, npart, dumxyzmh, f1vxyzu)
      ENDIF

      icount = 0

      icount = icount + 1
      ikount = ikount + 1
      WRITE(*,*) ' Calculating neighbour changes'
      icall = 1

      nlst = 0
      DO i = 1, npart
         IF (iphase(i).NE.-1) THEN
            nlst = nlst + 1
            llist(nlst) = i
            iscurrent(i) = .TRUE.
         ENDIF
      END DO    
      itime = 0
c
c--Needs non-zero dt for radiative transfer in derivi
c
      dt = 1.0E-12

      nlst_in = 1
      nlst_end = nlst
      nlst0 = nlst
      imax = 1073741824
      imaxstep = imax/2
c
c--Must already know B/rho to call derivi
c  rather than do another density call here, temporarily use Bvol instead
c  this also avoids slight changes in B from setup to first dump file
c
      IF (imhd.EQ.idim) THEN
         IF (varmhd.EQ.'Brho') THEN
            varmhdtemp = varmhd
            varmhd = 'Bvol'
         ELSE
            varmhdtemp = varmhd
         ENDIF
      ENDIF
c
c--Call derivi to get div B, curl B etc initially
c
      DO i = 1, ntot
         dumekcle(3,i) = 1.0
      END DO

      CALL derivi(dt,itime,dumxyzmh,dumvxyzu,f1vxyzu,f1ha,npart,ntot,
     &            ireal,dumalpha,dumekcle,dumBevolxyz,f1Bxyz)
c
c--Now set B/rho from B (ie. now that we know rho)
c
      IF (imhd.EQ.idim) THEN
         varmhd = varmhdtemp
         IF (varmhd.EQ.'Brho') THEN
            WRITE(*,99145)
99145       FORMAT(' Calculating B/rho from B for MHD evolution...')
            DO i=1,npart
               IF (rho(i).LE.0.) STOP 'ERROR: rho = 0 setting up B/rho'
               rho1i = 1./rho(i)
               DO j=1,3
                  Bevolxyz(j,i) = Bevolxyz(j,i)*rho1i
               ENDDO
            ENDDO
         ENDIF
      ENDIF
c
c--set h to new self-consistent value
c      
      DO ipart = 1,ntot
         xyzmh(5,ipart) = dumxyzmh(5,ipart)
      ENDDO
c
c--Get h max and min and print to screen
c
      hhmin = 1.e12
      hhmax = -1.e12
      hhav = 0.
      DO ipart=1,npart
         IF (xyzmh(5,ipart).LT.hhmin) hhmin = xyzmh(5,ipart)
         IF (xyzmh(5,ipart).GT.hhmax) hhmax = xyzmh(5,ipart)
         hhav = hhav + xyzmh(5,ipart)
      ENDDO
      hhav = hhav/real(npart)
      WRITE(*,*) ' min h = ',hhmin,' max h = ',hhmax,' av h =',hhav

      DO i = 1, npart
         iscurrent(i) = .FALSE.
      END DO

      RETURN
      END
