      SUBROUTINE cartdis(igeom,idist,np,h1)
c************************************************************
c                                                           *
c     This subroutine positions particles in cartesian      *
c              coordinate distribution                      *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/maspres'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/regionslocal'

      CHARACTER*1 iok
c
c--Allow for tracing flow
c
      third = 1./3.
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine cartdis ')
c
c--Set Condensed Areas In Cartesian Coordinates
c
      IF (idist.EQ.1 .OR. idist.EQ.2) THEN
         WRITE (*,*) 'Enter amplitude of sine wave '
         READ (*,*) ampl
c
c--Shift particles to sinusoidal density profile (from Daniel Price)
c
         itsmax = 100
         tol = 1.0e-5
         dxmax = xmax - xmin
         wk = 2.0*pi/(dxmax/1.)
         denom = dxmax + ampl/wk*(COS(wk*dxmax)-1.0)

         DO i = 1,npart

            dxi = xyzmh(1,i)-xmin
            dxprev = dxmax*2.
            xmassfrac = dxi/dxmax ! current mass fraction
                                  ! (for uniform density)
c
c--Use rootfinder on the integrated density perturbation
c  to find the new position of the particle
c    
            its = 0

            DO WHILE ((abs(dxi-dxprev).GT.tol).AND.(its.LT.itsmax))
               dxprev = dxi
               func = xmassfrac*denom - (dxi +ampl/wk*(COS(wk*dxi)-1.0))
               fderiv = -1.0 + ampl*SIN(wk*dxi)
               dxi = dxi - func/fderiv ! Newton-Raphson iteration
               its = its + 1 
c              PRINT*,'iteration',its,'dxi =',dxi,xmin,xmax,xyzmh(1,i)
            END DO

            IF (its.GE.itsmax) THEN
               WRITE (*,*) 'Error: soundwave - too many iterations'
               CALL quit(0)
            ENDIF

c            PRINT *,xyzmh(1,i),xmin + dxi,dxi

            xyzmh(1,i) = xmin + dxi
         END DO
      ELSE
         npart = np + nptmass
         rcyl2 = rcyl * rcyl 
         rmind2 = rmind * rmind
         rmax2 = rmax * rmax

 40      WRITE (*, 88002)
88002    FORMAT (' Enter number of different regions (max 10)')   
         READ (iread,*) ireg
         IF (ireg.GT.10) GOTO 40

         icount = 1
         DO 120 i = 1, ireg
 60         WRITE (*, 88004) icount
88004       FORMAT ('   Enter xmin of region ', I5)
            READ (iread,*) rxmin(icount)
            WRITE (*, 88006) icount
88006       FORMAT ('   Enter xmax of region ', I5)
            READ (iread,*) rxmax(icount)
            WRITE (*, 88008) icount
88008       FORMAT ('   Enter ymin of region ', I5)
            READ (iread,*) rymin(icount)
            WRITE (*, 88010) icount
88010       FORMAT ('   Enter ymax of region ', I5)
            READ (iread,*) rymax(icount)
            WRITE (*, 88012) icount
88012       FORMAT ('   Enter zmin of region ', I5)
            READ (iread,*) rzmin(icount)
            WRITE (*, 88014) icount
88014       FORMAT ('   Enter zmax of region ', I5)
            READ (iread,*) rzmax(icount)

 80         WRITE (*, 88016) icount
88016       FORMAT ('   Enter relative density of region ', I5, 
     &             ' (0.0 to 1.0)')
            READ (iread,*) dens(icount)
            IF ((dens(icount).LT.0.).OR.(dens(icount).GT.1.)) GOTO 80

 100        WRITE (*, 88018) icount
88018       FORMAT (' Is region ', I2,' correct (y/n)? ')
            READ (iread, 99004) iok
            IF ((iok.NE.'y').AND.(iok.NE.'n')) GOTO 100
            IF (iok.NE.'y') GOTO 60

            icount = icount + 1
 120     CONTINUE

         xmax5 = xmax * 0.5
         ymax5 = ymax * 0.5
         zmax5 = zmax * 0.5
         probavn = 0.

         DO 240 i = nptmass + 1, npart
 200        x1 = 2.*(xmax*ran1(1) - xmax5)
            y1 = 2.*(ymax*ran1(1) - ymax5)
            z1 = 2.*(zmax*ran1(1) - zmax5)

            rc2 = x1*x1 + y1*y1
            r2 = rc2 + z1*z1
            IF ((igeom.EQ.2) .AND. (rc2.GT.rcyl2)) GOTO 200
            IF ((igeom.EQ.3) .AND. (r2.GT.rmax2)) GOTO 200
            IF ((igeom.EQ.4).AND.((rc2.GT.rcyl2).OR.(rc2.LT.rmind2))) 
     &          GOTO 200

            rnd = ran1(1)
            DO 220 j = 1, ireg
               IF ((x1.GE.rxmin(j)) .AND. (x1.LE.rxmax(j)) .AND.
     &             (y1.GE.rymin(j)) .AND. (y1.LE.rymax(j)) .AND.
     &             (z1.GE.rzmin(j)) .AND. (z1.LE.rzmax(j))) THEN
                  IF (rnd.GT.dens(j)) GOTO 200
                  probavn = probavn + dens(j)
                  xyzmh(1,i) = x1
                  xyzmh(2,i) = y1
                  xyzmh(3,i) = z1
                  xyzmh(5,i) = (1.0/dens(j)) ** third
                  GOTO 240
               ENDIF
 220        CONTINUE

 240     CONTINUE 

         probavn = probavn/FLOAT(npart-nptmass)
         WRITE (*,*) probavn
         DO 300 i = nptmass + 1, npart
            xyzmh(5,i) = xyzmh(5,i) * h1 * probavn ** third
 300     CONTINUE

      ENDIF

      DO i = nptmass + 1, npart
         disfrac(i) = 1.0
      END DO

      RETURN
      END
