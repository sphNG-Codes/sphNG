      SUBROUTINE smoothd
c************************************************************
c                                                           *
c  This subroutine smoothes the initial conditions.         *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/dum'
      INCLUDE 'COMMONS/tlist'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/f1'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'

      REAL*4 sm, q
      DIMENSION sm(idim), q(idim)

      EQUIVALENCE (sm, pr), (q, vsound)

      CHARACTER*1 iok
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine smooth')
c
c--Compute tables for kernel quantities
c
      CALL ktable
c
c--Make ghosts
c
      DO i = 1, npart
         iscurrent(i) = .TRUE.
      END DO

      CALL ghostp(ntot, npart, xyzmh, vxyzu, ekcle, Bevolxyz)
c
c--Build tree
c
      IF (igrape.EQ.0) THEN
         DO i = 1, ntot
            DO j = 1, 5
               dumxyzmh(j,i) = xyzmh(j,i)
            END DO
         END DO
         WRITE(*,*) ' Making tree'
         CALL insulate(1, ntot, npart, dumxyzmh, f1vxyzu)
      ENDIF
      acc = 50.
c
c--Smooth all needed quantities
c
      DO 300 ismooth = 1, 5
         IF (ismooth.EQ.1) THEN
            WRITE (*,*) ' smoothing density'
            DO i = 1, npart
               q(i) = 1.
            END DO
         ELSEIF (ismooth.EQ.2) THEN
            WRITE (*,*) ' smoothing energy'
            DO i = 1, npart
               q(i) = vxyzu(4,i)/rho(i)
            END DO
         ELSEIF (ismooth.EQ.3) THEN
            WRITE (*,*) ' smoothing velocity vx'
            DO i = 1, npart
               q(i) = vxyzu(1,i)/rho(i)
            END DO
         ELSEIF (ismooth.EQ.4) THEN
            WRITE (*,*) ' smoothing velocity vy'
            DO i = 1, npart
               q(i) = vxyzu(2,i)/rho(i)
            END DO
         ELSEIF (ismooth.EQ.5) THEN
            WRITE (*,*) ' smoothing velocity vz'
            DO i = 1, npart
               q(i) = vxyzu(3,i)/rho(i)
            END DO
         ENDIF
c
c--Initialise array
c
         DO i = 1, npart
            sm(i) = 0.
         END DO
       
         WRITE (*,*) ' do you want a minimum h for smoothing? (y/n)'
         READ (*,99010) iok
99010    FORMAT (A1)
         IF (iok.EQ.'y') THEN 
            hmn = 1000.
            hmx = 0.
            DO i = 1, npart
               IF (xyzmh(5,i).LT.hmn) hmn = xyzmh(5,i)
               IF (xyzmh(5,i).GT.hmx) hmx = xyzmh(5,i)
            END DO
            hmins = hmn
            WRITE (*,99011) hmn, hmx
99011       FORMAT ('enter minimum smoothing length', 
     &           'min h = ',1PE11.4,'  max h = ',1PE11.4)
            READ (*,*) hmins
         ENDIF
c
c--Get neighbours
c
         DO i = 1, ntot
            DO j = 1, 5
               dumxyzmh(j,i) = xyzmh(j,i)
            END DO
         ENDDO
         IF (igrape.EQ.0) THEN
            CALL insulate(3, ntot, npart, dumxyzmh, f1vxyzu)
         ELSEIF (igrape.EQ.1) THEN
            CALL insulate(4, ntot, npart, dumxyzmh, f1vxyzu)
         ENDIF
c
c--Compute smoothed variable for each particle
c
         DO 100 ipart = 1, npart
            xi = xyzmh(1,ipart)
            yi = xyzmh(2,ipart)
            zi = xyzmh(3,ipart)
            pmassi = xyzmh(4,ipart)
            pmassqi = pmassi*q(ipart)
            hi = xyzmh(5,ipart)
            IF (hi.GT.hmins) GOTO 100

            smi = 0.
            DO k = 1, nneigh(ipart)
               IF (k.GE.nlmax) THEN
                  j = neighover(k-nlmax+1,ABS(neighb(nlmax,ipart)))
               ELSE
                  j = neighb(k,ipart)
               ENDIF
c
c--Define mean h
c
               hmean = 0.5*(hi + xyzmh(5,j))
               hmean21 = 1./hmean**2
               hmean31 = hmean21/hmean
               v2 = ((xi-xyzmh(1,j))**2 +(yi-xyzmh(2,j))**2 +
     &              (zi-xyzmh(3,j))**2)*hmean21
c
c--Get kernel quantities from interpolation in table
c
               index = v2*ddvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dwdx = (wij(index1) - wij(index))*ddvtable
               wtij = (wij(index) + dwdx*dxx)*hmean31
c
c--Add contribution
c
               smi = smi + q(j)*xyzmh(4,j)*wtij
               sm(j) = sm(j) + pmassqi*wtij

            END DO
            sm(ipart) = sm(ipart) + smi
 100     CONTINUE
c
c--Normalisation
c
         DO i = 1, npart
            sm(i) = cnormk*(sm(i) + q(i)*xyzmh(4,i)/xyzmh(5,i)**3)
         END DO
c
c--Load appropriate array
c
         IF (ismooth.EQ.1) THEN
            DO i = 1, npart
               rho(i) = sm(i)
            END DO
         ELSEIF (ismooth.EQ.2) THEN
            DO i = 1, npart
               vxyzu(4,i) = sm(i)
            END DO
         ELSEIF (ismooth.EQ.3) THEN
            DO i = 1, npart
               vxyzu(1,i) = sm(i)
            END DO
         ELSEIF (ismooth.EQ.4) THEN
            DO i = 1, npart
               vxyzu(2,i) = sm(i)
            END DO
         ELSEIF (ismooth.EQ.5) THEN
            DO i = 1, npart
               vxyzu(3,i) = sm(i)
            END DO
         ENDIF
 300  CONTINUE

      DO i = 1, npart
         iscurrent(i) = .FALSE.
      END DO
c
c--If idebug='density' the densities are dumped
c
      IF (idebug.EQ.'density') THEN
         WRITE (iprint, 99002) (rho(i), i=1, npart)
99002    FORMAT (1X, 6(1PE12.6,1X))
      ENDIF

      RETURN
      END
