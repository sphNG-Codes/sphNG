      SUBROUTINE cartdis(igeom,idist,np,h1,isettletest)
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
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/astrcon'

      CHARACTER*1 iok, idirect
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
         WRITE (*, 99002)
99002    FORMAT (' Do you want to stretch the particles? (y/n)')
         READ (iread, 99004) iok
         IF (iok.EQ.'y') THEN
            WRITE (*, 99003)
99003       FORMAT (' Choose how you want to stretch the particles: ',/,
     &           '            horizontal sine wave (wave test) : (1)',/,
     &           '           vertical gaussian (settling test) : (2)')
            READ (iread,*) istretch
c
c--Shift particles to sinusoidal density profile (from Daniel Price)
c
            IF (istretch.EQ.1) THEN
               WRITE (*,*) 'Enter amplitude of sine wave '
               READ (*,*) ampl
               WRITE (*,*) 'Enter direction of sine wave (x,y,z) '
               READ (*,99004) idirect
               itsmax = 100
               tol = 1.0e-5
               IF (idirect.EQ.'x') THEN
                  distmax = xmax - xmin
                  index = 1
                  coordmin = xmin
               ELSEIF (idirect.EQ.'y') THEN
                  distmax = ymax - ymin
                  index = 2
                  coordmin = ymin
               ELSEIF (idirect.EQ.'z') THEN
                  distmax = zmax - zmin
                  index = 3
                  coordmin = zmin
               ELSE
                  WRITE (*,*) 'Invalid wave direction'
                  CALL quit(0)
               ENDIF
               wk = 2.0*pi/(distmax/1.)
               denom = distmax + ampl/wk*(COS(wk*distmax)-1.0)

               DO i = 1, npart
                  disti = xyzmh(index,i)-coordmin
                  dprev = distmax*2.
                  xmassfrac = disti/distmax ! current mass fraction
                                            ! (for uniform density)
c
c--Use rootfinder on the integrated density perturbation
c  to find the new position of the particle
c    
                  its = 0

                  DO WHILE ((ABS(disti-dprev).GT.tol)
     &                 .AND. (its.LT.itsmax))
                     dprev = disti
                     func = xmassfrac*denom - (disti + 
     &                    ampl/wk*(COS(wk*disti)-1.0))
                     fderiv = -1.0 + ampl*SIN(wk*disti)
                     disti =disti-func/fderiv ! Newton-Raphson iteration
                     its = its + 1 
                  END DO

                  IF (its.GE.itsmax) THEN
                     WRITE (*,*)'Error: soundwave - too many iterations'
                     CALL quit(0)
                  ENDIF

                  xyzmh(index,i) = coordmin + disti
               END DO
            ELSEIF (istretch.EQ.2) THEN
               isettletest = 1
c                                                                
c--Get total mass from integration of density profile             
c
               itsmax = 100
               tol    = 1.e-5
               dxmax  = xmax - xmin
               dymax  = ymax - ymin
               dzmax  = zmax - zmin

               WRITE (*,*) 'Warning! This test makes some hard-coded',
     &                     'assumptions about the disc params...:'

               WRITE (*,*) '...Turning on vertical gravity from a ',
     &                     'centrally-located, solars-mass star'
               iexf = 11 ! i.e. xmass = 1. (set in setpart.f)
               WRITE (*,*) '...Assuming zmax = 4*H where H is the gas ',
     &                     'scale-height of the disc'
               H0 = zmax/4.
               WRITE (*,*) '...Assuming a H/R = 0.05'
               honr  = 0.05
               Rdisc = (H0/honr)
               WRITE (*,*) '...which puts the simulated disc cross ',
     &                     'section at radius (au) = ',
     &                     Rdisc*udist/au
               IF (ABS(Rdisc*udist/au-50.).GT.0.01) THEN
                  zmax_suggested = (50.*au/udist)*(4.*honr)
                  WRITE (*,*) 'Oops...externf.f requires a radius of ',
     &                        '50 au for this test. Try again with ',
     &                        'a zmax =',zmax_suggested
                  CALL quit(0)
               ENDIF

               ampl = 10.
               tot_mass = dzmax + ampl*2.*SQRT(0.5*pi)*H0*
     &                            ERF(zmax/(SQRT(2.)*H0))
               DO i = 1,npart
                  dzi = xyzmh(3,i)
                  dzprev = 2.*dzmax
                  zmassfrac = (dzi-zmin)/dzmax ! current mass fraction
                                               ! (for uniform density)
c
c--Use rootfinder on the integrated density perturbation
c  to find the new position of the particle
c
                  its = 0
                  DO WHILE ((ABS(dzi-dzprev).GT.tol)
     &                         .AND.(its.LT.itsmax))
                     dzprev = dzi
                     func = zmassfrac*tot_mass - (dzi - zmin +
     &                      ampl*H0*SQRT(pi/2.)*(ERF(dzi/(SQRT(2.)*H0))
     &                      + ERF(-zmin/(SQRT(2.)*H0))))

                     fderiv = -(1.0 + ampl*EXP(-0.5*(dzi/H0)**2))
                     ! Newton-Raphson iteration
                     dzi = dzi - func/fderiv
                     its = its + 1
                  END DO

                  IF (its.GE.itsmax) THEN
                     WRITE (*,*) 'Error: too many iterations'
                     CALL quit(0)
                  ENDIF
                  xyzmh(3,i) = dzi
               ENDDO
               !--expand the z boundaries to avoid ghost particles
               zmin = 2.*zmin
               zmax = 2.*zmax
            ENDIF
         ENDIF
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
