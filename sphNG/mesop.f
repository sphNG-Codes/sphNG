      SUBROUTINE mesop
c************************************************************
c                                                           *
c  This subroutine handles the messages received from the   *
c     operator.                                             *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/stop'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/timeextra'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/secret'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/sync'
      INCLUDE 'COMMONS/polyk2'

      CHARACTER*20 string(10)
      CHARACTER*7 where
      CHARACTER*5 dummy5
      CHARACTER*3 dummy3

      DATA where/'mesop'/
c
c--Open message file first
c
      xlog2 = 0.301029996

      iwait = 0
      OPEN (idisk2, FILE='message', FORM='formatted')
c
c--Read message file (one mes. per line , max=10)
c
      nmes = 0
      DO i = 1, 10
         READ (idisk2, 99000, END=200) string(i)
99000    FORMAT (A20)
         nmes = i
      END DO
c
c--Read messages one ofter another
c
 200  DO 300 i = 1, nmes
c
c--Identify each message
c
c  0) wait before processing messages
c
         IF (string(i)(1:2).EQ.'at') THEN
            READ (string(i), *) dummy3, attime
            IF (gt.LT.attime) THEN
               iwait = 1
               WRITE (*, 99001) i
99001          FORMAT (' message ', I2, ' read. waiting to process.')
               GOTO 400
            ELSE
               GOTO 300
            ENDIF
         ENDIF
c
c  1) stop job
c
         IF (string(i)(1:4).EQ.'stop') THEN
            istop = 1
            WRITE (*, 99002) i
99002       FORMAT (' message ', I2, ' read. please wait for stop.')
            GOTO 300
         ENDIF
c
c  2) change dump frequency
c
         IF (string(i)(1:5).EQ.'nstep') THEN
            READ (string(i), 99005) nstep
99005       FORMAT (6X, I2)
            WRITE (*, 99006) i, nstep
99006       FORMAT (' message ', I2, ' read. nstep set to : ', I2)
            GOTO 300
         ENDIF
c
c  3) change syncronisation time
c
         IF (string(i)(1:5).EQ.'synct') THEN
            READ (string(i), 99007) ipower
99007       FORMAT (6X, I2)
            WRITE (*, 99008) i, ipower
99008       FORMAT (' message ', I2, ' read. synct changed by : ', I2)

            ifactor = 2**ABS(ipower)
            imaxnew = imaxstep/ifactor
            iminnew = 2*ifactor
 
            IF (ipower.LT.0) THEN
               dtmax = dtmax/FLOAT(ifactor)
               istepmin = MIN(istepmin, imaxnew)
               istepmax = MIN(istepmax, imaxnew)
               istepmin = istepmin*ifactor
               istepmax = istepmax*ifactor
 
               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     isteps(j) = MIN(isteps(j), imaxnew)
                     isteps(j) = isteps(j)*ifactor
                     it1(j) = isteps(j)/2
                     it2(j) = isteps(j)
                  ENDIF
               END DO
            ELSEIF (ipower.GT.0) THEN
               IF (istepmin/ifactor .LT. 2) THEN
                  CALL error(where, 1)
                  istop = 1
                  GOTO 300
               ENDIF
               dtmax = dtmax*FLOAT(ifactor)
               istepmin = istepmin/ifactor
               istepmax = istepmax/ifactor
 
               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     IF (isteps(j)/ifactor .LT. 2) CALL error(where, 2)
                     isteps(j) = isteps(j)/ifactor
                     it1(j) = isteps(j)/2
                     it2(j) = isteps(j)
                  ENDIF
               END DO
            ENDIF

            DO j = 1, nbinmax
               nlstbins(j) = 0
               it1bin(j) = 2**j/2
               it2bin(j) = 2**j
            END DO

            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  ibin = INT(LOG10(REAL(isteps(j)))/xlog2+0.5)
                  IF (ibin.GT.nbinmax) THEN
                     WRITE (*,*) 'ERROR - ibin.GT.nbinmax mes'
                     WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax mes'
                     CALL quit
                  ENDIF
                  nlstbins(ibin) = nlstbins(ibin) + 1
                  listbins(nlstbins(ibin),ibin) = j
                  IF (integrator.EQ.0) THEN
                     IF (it1bin(ibin).NE.it1(j)) THEN
                        WRITE (*,*) 'ERROR - it1bin mes'
                        CALL quit
                     ENDIF
                  ENDIF
                  IF (it2bin(ibin).NE.it2(j)) THEN
                     WRITE (*,*) 'ERROR - it2bin mes'
                     CALL quit
                  ENDIF
               ENDIF
            END DO
 
            GOTO 300
         ENDIF
c
c  4) change minimum time for keeping GRAPE
c
         IF (string(i)(1:5).EQ.'tkeep') THEN
            READ (string(i), *) dummy5, tkeep
            WRITE (*, 99010) i, tkeep
99010    FORMAT (' message ', I2, ' read. tkeep changed to : ',1PE12.5)
            GOTO 300
         ENDIF
c
c--Unexpected messages
c
         WRITE (*, 99090) string(i)
99090    FORMAT (' unexpected message received : ', /, A20)

 300  CONTINUE
c
c--Erase message file
c
 400  IF (iwait.EQ.0) THEN
         CLOSE (idisk2, STATUS='delete')
      ELSE
         CLOSE (idisk2)
      ENDIF
c
c--Automatic change of syncronisation time using local free-fall time
c
      OPEN (idisk2, FILE='synctcontrol', FORM='formatted')
      READ (idisk2, 99000, END=500) string(i)

      IF (string(i)(1:4).EQ.'auto') THEN
         READ (string(i)(6:20),*) fffac
         WRITE (*, 99016) fffac
         WRITE (iprint, 99016) fffac
99016    FORMAT (' Synctcontrol : ', 1PE12.5)

         ipower = 0
         rhomax = 0.
         DO j = 1, npart
            IF (rho(j).GT.rhomax) rhomax = rho(j)
         END DO
         freefalltime = SQRT((3 * pi) / (32 * rhomax))
         timemax = fffac*freefalltime
         IF (dtmax.GT.timemax) THEN
            ipower = - (INT(LOG10(dtmax/timemax)/0.30103) + 1)
            WRITE (*, 99009) ipower, freefalltime, rhomax
            WRITE (iprint, 99009) ipower, fffac, freefalltime, rhomax
99009       FORMAT (' Synct autochange : ', I2,' fffac ',1PE12.5,
     &        ' free-fall ',1PE12.5,' rhomax ',1PE12.5)

            ifactor = 2**ABS(ipower)
            imaxnew = imaxstep/ifactor
            iminnew = 2*ifactor

            IF (ipower.LT.0) THEN
               dtmax = dtmax/FLOAT(ifactor)
               istepmin = MIN(istepmin, imaxnew)
               istepmax = MIN(istepmax, imaxnew)
               istepmin = istepmin*ifactor
               istepmax = istepmax*ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     isteps(j) = MIN(isteps(j), imaxnew)
                     isteps(j) = isteps(j)*ifactor
                     it1(j) = isteps(j)/2
                     it2(j) = isteps(j)
                  ENDIF
               END DO
            ELSEIF (ipower.GT.0) THEN
               IF (istepmin/ifactor .LT. 2) THEN
                  CALL error(where, 1)
                  istop = 1
                  GOTO 500
               ENDIF
               dtmax = dtmax*FLOAT(ifactor)
               istepmin = istepmin/ifactor
               istepmax = istepmax/ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     IF (isteps(j)/ifactor .LT. 2) CALL error(where, 2)
                     isteps(j) = isteps(j)/ifactor
                     it1(j) = isteps(j)/2
                     it2(j) = isteps(j)
                  ENDIF
               END DO
            ENDIF

            DO j = 1, nbinmax
               nlstbins(j) = 0
               it1bin(j) = 2**j/2
               it2bin(j) = 2**j
            END DO

            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  ibin = INT(LOG10(REAL(isteps(j)))/xlog2+0.5)
                  IF (ibin.GT.nbinmax) THEN
                     WRITE (*,*) 'ERROR - ibin.GT.nbinmax mes'
                     WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax mes'
                     CALL quit
                  ENDIF
                  nlstbins(ibin) = nlstbins(ibin) + 1
                  listbins(nlstbins(ibin),ibin) = j
                  IF (integrator.EQ.0) THEN
                     IF (it1bin(ibin).NE.it1(j)) THEN
                        WRITE (*,*) 'ERROR - it1bin mes'
                        CALL quit
                     ENDIF
                  ENDIF
                  IF (it2bin(ibin).NE.it2(j)) THEN
                     WRITE (*,*) 'ERROR - it2bin mes'
                     CALL quit
                  ENDIF
               ENDIF
            END DO
         ENDIF
      ENDIF
c
c--Automatic change of syncronisation time using density change
c
      IF (string(i)(1:4).EQ.'dens') THEN
         READ (string(i)(6:20),*) contrast, dtmaxsync

         rhomaxsyncold = rhomaxsync
         ipower = 0
         rhomaxsync = 0.
         DO j = 1, npart
            rhomaxsync = MAX(rhomaxsync, rho(j))
         END DO

         freefalltime = SQRT((3 * pi) / (32 * rhozero))

         ratio = LOG10(rhomaxsync/rhomaxsyncold)

         contrast = LOG10(contrast)
         ipower = - ABS(INT(ratio/contrast))
         IF (ABS(ratio/contrast).LT.0.5) ipower = 1

         WRITE (*, 99017) dtmaxsync, 10**contrast, 
     &        rhomaxsync/rhomaxsyncold, ipower
         WRITE (iprint, 99017) dtmaxsync, 10**contrast, 
     &        rhomaxsync/rhomaxsyncold, ipower
99017    FORMAT (' Dtmaxsync : ', 1PE12.5,' contrast ',1PE12.5,
     &        ' ratio ',1PE12.5,' ipower ',I2)

         IF (ipower.GE.0 .AND. tstep/60.0.GT.720.0) THEN
            WRITE (iprint,99034) tstep/60.0
            WRITE (*,99034) tstep/60.0
99034       FORMAT('  Synct autochange since time ',1PE12.5)
            ipower = -1
         ENDIF
         IF (ipower.EQ.1) THEN
            tempvar = gt/(2.0*dtmax)
            diff = tempvar-INT(tempvar)
            IF (diff.GT.0.25 .AND. diff.LT.0.75) THEN
               WRITE (iprint,99035) diff
               WRITE (*,99035) diff
99035          FORMAT('  Synct autochange attempt, but sync ',1PE12.5)
               ipower = 0
            ENDIF
         ENDIF
         IF (ipower.NE.0) THEN

            ifactor = 2**ABS(ipower)
            imaxnew = imaxstep/ifactor
            iminnew = 2*ifactor

            IF (ipower.LT.0) THEN
               WRITE (*, 99021) ipower
               WRITE (iprint, 99021) ipower
99021          FORMAT ('  Synct autochange : ', I2)
               dtmax = dtmax/FLOAT(ifactor)
               istepmin = MIN(istepmin, imaxnew)
               istepmax = MIN(istepmax, imaxnew)
               istepmin = istepmin*ifactor
               istepmax = istepmax*ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     isteps(j) = MIN(isteps(j), imaxnew)
                     isteps(j) = isteps(j)*ifactor
                     it1(j) = isteps(j)/2
                     it2(j) = isteps(j)
                  ENDIF
               END DO
            ELSEIF (ipower.GT.0 .AND. 
     &              dtmax*FLOAT(ifactor).LE.dtmaxsync .AND.
     &              tstep/60.*FLOAT(ifactor).LE.720.) THEN
               IF (istepmin/ifactor .LT. 2) THEN
                  CALL error(where, 1)
                  istop = 1
                  GOTO 500
               ENDIF
               WRITE (*, 99021) ipower
               WRITE (iprint, 99021) ipower

               dtmax = dtmax*FLOAT(ifactor)
               istepmin = istepmin/ifactor
               istepmax = istepmax/ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     IF (isteps(j)/ifactor .LT. 2) CALL error(where, 2)
                     isteps(j) = isteps(j)/ifactor
                     it1(j) = isteps(j)/2
                     it2(j) = isteps(j)
                  ENDIF
               END DO
            ELSE
               IF (dtmax*FLOAT(ifactor).GT.dtmaxsync) THEN
                  WRITE (iprint,99040)
                  WRITE (*,99040)             
99040       FORMAT('  Synct autochange attempt, but dtmaxsync reached')
               ELSEIF (tstep/60.*FLOAT(ifactor).GT.720.) THEN
                  WRITE (iprint,99041)
                  WRITE (*,99041)             
99041       FORMAT('  Synct autochange attempt, but tstepmax reached')
               ENDIF
            ENDIF

            DO j = 1, nbinmax
               nlstbins(j) = 0
               it1bin(j) = 2**j/2
               it2bin(j) = 2**j
            END DO

            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  ibin = INT(LOG10(REAL(isteps(j)))/xlog2+0.5)
                  IF (ibin.GT.nbinmax) THEN
                     WRITE (*,*) 'ERROR - ibin.GT.nbinmax mes'
                     WRITE (iprint,*) 'ERROR - ibin.GT.nbinmax mes'
                     CALL quit
                  ENDIF
                  nlstbins(ibin) = nlstbins(ibin) + 1
                  listbins(nlstbins(ibin),ibin) = j
                  IF (it1(j).NE.0) THEN
                     IF (it1bin(ibin).NE.it1(j)) THEN
                        WRITE (*,*) 'ERROR - it1bin mes'
                        CALL quit
                     ENDIF
                  ENDIF
                  IF (it2bin(ibin).NE.it2(j)) THEN
                     WRITE (*,*) 'ERROR - it2bin mes'
                     CALL quit
                  ENDIF
               ENDIF
            END DO
         ENDIF
      ENDIF

 500  CLOSE(idisk2)
      RETURN
      END
