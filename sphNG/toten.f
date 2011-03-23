      SUBROUTINE toten
c************************************************************
c                                                           *
c  This routine computes all energies per unit total mass   *
c     kinetic, rotational, potential and internal.          *
c     It also computes the escapor mass.                    *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/Bxyz'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine toten')
c
c--Compute first mechanical energies and escapors
c
      gama1 = gamma - 1.
      tkin = 0.
      trotz = 0.
      trotx = 0.
      tgrav = 0.
      tpot = 0.
      tterm = 0.
      tmag = 0.
      trad = 0.
      escap = 0.
c
c--Scaling factors
c
      CALL scaling(gt, rscale, drdt, dlnrdt)

      DO i = 1, npart
         IF (iphase(i).GE.0) THEN
            xi = xyzmh(1,i)
            yi = xyzmh(2,i)
            zi = xyzmh(3,i)
            pmassi = xyzmh(4,i)
            vxi = vxyzu(1,i)
            vyi = vxyzu(2,i)
            vzi = vxyzu(3,i)
c
c--Kinetic energy
c
            vtot2 = vxi**2 + vyi**2 + vzi**2
            tkin = tkin + pmassi*vtot2
c
c--Rotational energy around z
c
            r2xy = xi*xi + yi*yi
            rvz = xi*vyi - yi*vxi
            IF (r2xy.NE.0) trotz = trotz + pmassi*rvz*rvz/r2xy
c
c--Rotational energy around x
c
            r2yz = yi*yi + zi*zi
            rvx = yi*vzi - zi*vyi
            IF (r2yz.NE.0) trotx = trotx + pmassi*rvx*rvx/r2yz
c
c--Potential energy
c
            poteni = pmassi*(poten(i) + dgrav(i))/rscale
            tgrav = tgrav + poteni
c
c--Magnetic energy
c
            IF (imhd.EQ.idim) THEN
               B2i = Bxyz(1,i)**2 + Bxyz(2,i)**2 + Bxyz(3,i)**2
               tmagi = pmassi*B2i/rho(i)
               tmag = tmag + tmagi
            ELSE
               tmagi = 0.
            ENDIF
c
c--Escapors
c
            IF (xi*vxi + yi*vyi + zi*vzi.LT.0.) vtot2 = 0.
            tinout = 0.5*pmassi*vtot2
            tpot = poteni
            IF (varsta.EQ.'intener') THEN
               ttherm = vxyzu(4,i)*pmassi
            ELSE
               ttherm = pmassi*pr(i)/(gama1*rho(i))
            ENDIF
            IF (encal.EQ.'r' .AND. iradtrans.EQ.idim) THEN
               tradi = ekcle(1,i)*pmassi
            ELSE
               tradi = 0.
            ENDIF
            total = tinout + tpot + ttherm + tradi + tmagi
            IF (total.GT.0.) escap = escap + pmassi
         ENDIF
      END DO
c
c--Normalisations
c
      tkin = 0.5*tkin
      trotz = 0.5*trotz
      trotx = 0.5*trotx
      tgrav = 0.5*tgrav
      tmag = 0.5*tmag
      escap = escap
c
c--Thermal energy
c
      IF (varsta.NE.'entropy') THEN
c
c--Variable of state is specific internal energy
c
         DO i = 1, npart
            IF (iphase(i).EQ.0) THEN
               tterm = tterm + xyzmh(4,i)*vxyzu(4,i)
               IF (encal.EQ.'r' .AND. iradtrans.EQ.idim) THEN
                  trad = trad + xyzmh(4,i)*ekcle(1,i)
               ENDIF
            ENDIF
         END DO
c
c--Variable of state is specific entropy
c
      ELSEIF (gama1.EQ.0.) THEN
         tterm = 1.5*vxyzu(4,1)
      ELSE
         DO i = 1, npart
            IF (iphase(i).EQ.0) tterm = 
     &           tterm + xyzmh(4,i)*vxyzu(4,i)*rho(i)**gama1/gama1
         END DO
         tterm = tterm
      ENDIF

      IF (idebug(1:5).EQ.'toten') THEN
         WRITE (iprint, 99002) tkin,trotz,trotx,tgrav,tterm,tmag,trad
99002    FORMAT (1X, 7(1PE12.5,1X))
      ENDIF

      RETURN
      END
