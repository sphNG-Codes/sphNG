      PROGRAM ptmass
c
c--Calculates the orbits and mass accretion of all point masses
c     using a BINARY output 'P file'
c
c     This code also re-orders sinks for MPI dumps according to iunique
c
c     The code also handles merges properly (both mergers that are
c     recorded in the P-files, and mergers that occur without being
c     recorded in the P-files.  
c
      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL (A-C, E-H, O-Z)

      PARAMETER (idim = 1300)
      PARAMETER (nbinimf = 100)
      PARAMETER (nbinbin = 40)
      PARAMETER (nmassratio = 10)

      PARAMETER (bdmass = 0.10)
      PARAMETER (bdrealmass = 0.075)

      PARAMETER (pi = 3.14159265)

      REAL*4 rho4,ptmass4
      INTEGER*2 nactotal2

      DIMENSION ipart(idim), x(idim), y(idim), 
     &     z(idim), vx(idim), 
     &     vy(idim), vz(idim), pmass(idim), 
     &     rho(idim), nactotal(idim),
     &     ptmassinner(idim), spinx(idim), 
     &     spiny(idim), spinz(idim)
      DIMENSION xn(idim), yn(idim), zn(idim), vxn(idim), 
     &     vyn(idim), vzn(idim), pmassn(idim), name(idim), lname(idim),
     &     nbincomp(idim),interest(idim),pmassmax(idim),pmassmin(idim),
     &     tform(idim),tffform(idim),vrel(idim),vrelb(idim),
     &     ejectime(idim),orbit(3,idim),period(idim),periodratio(idim),
     &     ejectimelast(idim),radcom(idim), radmax(idim)
      DIMENSION oldmasses(idim,10000), oldtimes(idim,10000), 
     &     oldvx(idim,10000),oldvy(idim,10000),oldvz(idim,10000),
     &     accrete(idim), acctime(idim), iupto(idim),
     &     radiusmin(idim),
     &     radiuscreate(idim,idim),radiusend(idim,idim),
     &     accel(idim),endacctime(idim)
      DIMENSION nimf(nbinimf),nimfa(nbinimf)
      DIMENSION seplast(idim),seplast2(idim),qratiolast(idim)
      DIMENSION nbinary(nbinbin),nbinarybd(nbinbin),
     &     nbinarybdobs(nbinbin),ntriple(nbinbin),ntriplebd(nbinbin),
     &     nquad(nbinbin),
     &     nbinary_notrealbd(nbinbin),nbinary_01_02(nbinbin),
     &     nbinary_02_05(nbinbin),nbinary_05_08(nbinbin),
     &     nbinary_08_12(nbinbin),nbinary_12(nbinbin),
     &     nbinary_01_02_notrealbd(nbinbin),
     &     nbinary_02_05_notrealbd(nbinbin),
     &     nbinary_05_08_notrealbd(nbinbin),
     &     nbinary_08_12_notrealbd(nbinbin),
     &     nbinary_12_notrealbd(nbinbin),
     &     ntriple_01_02(nbinbin),
     &     ntriple_02_05(nbinbin),
     &     ntriple_05_08(nbinbin),
     &     ntriple_08_12(nbinbin),
     &     ntriple_12(nbinbin)
      DIMENSION nbinaryq(nmassratio),nbinarybdq(nmassratio),
     &     nbinarybdobsq(nmassratio),ntripleq(nmassratio),
     &     nbinaryq_notrealbd(nbinbin),nbinaryq_01_02(nbinbin),
     &     nbinaryq_02_05(nbinbin),nbinaryq_05_08(nbinbin),
     &     nbinaryq_08_12(nbinbin),nbinaryq_12(nbinbin),
     &     nbinaryq_01_02_notrealbd(nbinbin),
     &     nbinaryq_02_05_notrealbd(nbinbin),
     &     nbinaryq_05_08_notrealbd(nbinbin),
     &     nbinaryq_08_12_notrealbd(nbinbin),
     &     nbinaryq_12_notrealbd(nbinbin)

c
c--List merge is an array of pointers that point to the currently active
c     sinks in the list of all sinks (some of which may have died in mergers)
c     In other words listmerge(i) >= i
c
      DIMENSION listmerge(idim), invert(idim)

      PARAMETER (nmultiple_max = 1000)
      DIMENSION imultiple_values(20,nmultiple_max)
      DIMENSION xmultiple_values(20,nmultiple_max)
      DIMENSION noutputline(idim)

      INTEGER*8 iunique, listunique(idim),imerged1_int8,imerged2_int8
      LOGICAL*1 ifound(idim)

      CHARACTER*21 infile, outfile
      CHARACTER*19 outbase
      CHARACTER*4 num,tempstring
      CHARACTER*1 ans,iout,ansunique
      CHARACTER*60 name

      nmultipleout = 0

      iline_value = 0

      WRITE (*,*) 'Enter input filename:' 
      READ (*,99001) infile
99001 FORMAT(A21)
      WRITE (*,*) 'Do you want individual output files?'
      READ (*,99003) iout
      IF (iout.EQ.'Y' .OR. iout.EQ.'y') THEN
         WRITE (*,*) 'Enter base name for output files:' 
         READ (*,99002) outbase
99002    FORMAT(A19)
      ENDIF

      WRITE (*,*) 'Rotating reference frame?'
      READ (*,99003) ans
99003 FORMAT(A1)
      IF (ans.EQ.'y' .OR. ans.EQ.'Y') THEN
         WRITE (*,*) '  Rotation frequency?'
         READ (*,*) rotfreq
         WRITE (*,*) '  Free Fall Time (secs)?'
         READ (*,*) tunit
         omega = rotfreq*tunit
      ELSE
         omega = 0.0
      ENDIF

      WRITE (*,*) 'Number dumps to average accretion over (max 10000)?'
      READ (*,*) ndumps
      DO i = 1, ndumps
         oldmasses(1,i) = 0.
         oldtimes(1,i) = 0.
         oldmasses(2,i) = 0.
         oldtimes(2,i) = 0.
      END DO

      WRITE (*,*) 'Enter units of length in pc?'
      READ (*,*) ulength
      runit = 206283.42*ulength

      WRITE (*,*) 'Enter free-fall time to stop at'
      READ (*,*) stoptime

      WRITE (*,*) 'Does the dump file use iunique (y/n)?'
      READ (*,99003) ansunique
c
c--Compute
c
      nptmass = 0
      nptmasscurrent = 0
      nmerge = 0
      nlines = 0

      DO i = 1, idim
         listmerge(i) = i
         interest(i) = 0
         ejectime(i) = 0
      END DO

      iflag = 0

      OPEN (22, FILE='MERGED22', STATUS='unknown')
      OPEN (15, FILE=infile, STATUS='unknown', FORM='unformatted')
 100  nlines = nlines + 1
c
c--Read time and number of point masses
c
 110  nptmassold = nptmass
 111  READ (15, END=200, ERR=200) 
     &     fftime, time, nptmass
c         write (*,*) fftime, time, nptmass

         iline_value = iline_value + 1

      IF (fftime.GT.stoptime .AND. iflag.EQ.1) GOTO 200
      IF (fftime.GT.stoptime .AND. iflag.EQ.0) iflag=1
c      IF (fftime.GT.stoptime) GOTO 200

      IF (nptmass.LT.0) THEN
         IF (ansunique.EQ.'y' .OR. ansunique.EQ.'Y') THEN
            READ (15) imerged1_int8, imerged2_int8
            imerged1 = imerged1_int8
            imerged2 = imerged2_int8
         ELSE
            READ (15) imerged1, imerged2
         ENDIF
         WRITE (22,*) 'Merged from -ve nptmass ',imerged1, imerged2
         WRITE (*,*) 'Merged file ',imerged1, imerged2
         IF (iline_value.EQ.1) THEN
            WRITE (*,*) 'Merger at start of P-file - ignoring'
            GOTO 111
         ENDIF
         IF (ansunique.EQ.'y' .OR. ansunique.EQ.'Y') THEN
c
c--Search for iunique in existing sink list
c
            DO iii = 1, nptmasscurrent
               IF (imerged2.EQ.listunique(iii)) THEN
                  i = listmerge(iii)
                  iremove_sink = iii
                  GOTO 5556
               ENDIF
            END DO
            PRINT *,'Failed to find merge ',imerged2
            STOP
 5556       nmerge = nmerge + 1
c
c--Print existing sinks
c
            DO iii = 1, nptmasscurrent
               PRINT *,'     Existing ',iii,listmerge(iii),
     &              listunique(iii),pmass(listmerge(iii))
            END DO

            nptmassold = nptmassold - 1
            nptmasscurrent = nptmasscurrent - 1
            pmass(i) = -1.0
            DO iiii = iremove_sink, idim-1
               listmerge(iiii) = listmerge(iiii+1)
               listunique(iiii) = listunique(iiii+1)
            END DO
            listmerge(idim) = listmerge(idim-1)+1

c
c--Print new sinks
c
            DO iii = 1, nptmasscurrent
               PRINT *,'     New List ',iii,listmerge(iii),
     &              listunique(iii),pmass(listmerge(iii))
            END DO

            GOTO 111
         ELSE
            nmerge = nmerge + 1
            nptmassold = nptmassold - 1
            pmass(listmerge(imerged2)) = -1.0
            DO i = imerged2, idim
               listmerge(i) = listmerge(i+1)
            END DO
            listmerge(idim) = listmerge(idim-1)+1
            GOTO 111
         ENDIF
      ENDIF
      ifindmerge = 0
      IF (nptmassold .GT. nptmass) THEN
         WRITE (*,*) 'MERGER not done via -ve nptmass ',
     &        nptmass, nptmassold, fftime, time
         WRITE (22,*) 'MERGER not done via -ve nptmass ',
     &        nptmass, nptmassold, fftime, time
         nptmassold = nptmassold - 1
         ifindmerge = 1
      ENDIF
c
c--Set ifound to false for all current sink particles
c
      DO i = 1, nptmasscurrent
         ifound(i) = .FALSE.
      END DO

      IF (nptmass+nmerge.GT.idim) THEN
      WRITE (*,*) 'GREATER THAN ',idim,' POINT MASSES ',nptmass,nmerge
         STOP
      ENDIF
c
c--Read each time dump
c
      DO k = 1, nptmass
         i = listmerge(k)
         IF (ansunique.EQ.'y' .OR. ansunique.EQ.'Y') THEN
            READ (15, END=200, ERR=200) iunique, xi, 
     &           yi, zi, vxi, 
     &           vyi, vzi, pmassi, 
     &           rho4, nactotal2,
     &           ptmass4, 
     &           spinxi, spinyi, spinzi
c
c--Search for iunique in existing sink list
c
            DO iii = 1, nptmasscurrent
               IF (iunique.EQ.listunique(iii)) THEN
                  i = listmerge(iii)
                  ifound(iii) = .TRUE.
                  GOTO 555
               ENDIF
            END DO
            nptmasscurrent = nptmasscurrent + 1
            listunique(nptmasscurrent) = iunique
            i = listmerge(nptmasscurrent)
            WRITE (*,*) 'New sink found: ',iunique,' assigned ',
     &           nptmasscurrent, i
            ifound(nptmasscurrent) = .TRUE.

            tform(i) = fftime
            tffform(i) = time
            write (*,*) 'Setting formation ',fftime, time, nptmass, i

 555        CONTINUE

            x(i) = xi
            y(i) = yi
            z(i) = zi
            vx(i) = vxi
            vy(i) = vyi
            vz(i) = vzi
            pmass(i) = pmassi
            rho(i) = rho4
            nactotal(i) = nactotal2
            ptmassinner(i) = ptmass4
            spinx(i) = spinxi
            spiny(i) = spinyi
            spinz(i) = spinzi

c            WRITE (88,99088) fftime,pmass(i),spinx(i),spiny(i),spinz(i),
c     &           ATAN2(spinz(i),spinx(i))*180./pi,
c     &           ATAN2(spinz(i),spiny(i))*180./pi
99088       FORMAT(7(1PE12.5,1X))

            IF (iunique.GT.2000000000) THEN
               WRITE (*,*) 'ERROR - iunique > 2^31'
               STOP
            ENDIF

            ipart(i) = iunique
         ELSE
            READ (15, END=200, ERR=200) ipart(i), x(i),
     &           y(i), z(i), vx(i),
     &           vy(i), vz(i), pmass(i),
     &           rho4, nactotal2,
     &           ptmass4,
     &           spinx(i), spiny(i), spinz(i)

c            print *,ipart(i),x(i),vx(i),'pmass ',pmass(i),' rho ',rho4,
c     &           nactotal2,
c     &           ptmass4,' spins ',spinx(i), spiny(i), spinz(i)
         ENDIF


c         write (*,*) ipart(i), x(i), y(i)

         IF (iupto(i).LT.ndumps) iupto(i) = iupto(i) + 1

         r = SQRT(x(i)**2 + y(i)**2)
         th = ATAN2(y(i),x(i))
         th = th + time*omega
         x(i) = r*COS(th)
         y(i) = r*SIN(th)

         r = SQRT(vx(i)**2 + vy(i)**2)
         th = ATAN2(vy(i),vx(i))
         th = th + time*omega
         vx(i) = r*COS(th)
         vy(i) = r*SIN(th)

         DO j = 1, ndumps - 1
            oldmasses(i, j) = oldmasses(i, j + 1)
            oldtimes(i, j) = oldtimes(i, j + 1)
            oldvx(i, j) = oldvx(i, j + 1)
            oldvy(i, j) = oldvy(i, j + 1)
            oldvz(i, j) = oldvz(i, j + 1)
         END DO
         oldmasses(i, ndumps) = pmass(i)
         oldtimes(i, ndumps) = time
         oldvx(i, ndumps) = vx(i)
         oldvy(i, ndumps) = vy(i)
         oldvz(i, ndumps) = vz(i)

         IF (iupto(i).GT.1) THEN
            nlast = ndumps-iupto(i)+1
            IF (time-oldtimes(i,nlast).NE.0.0) THEN
               accrete(i) = (pmass(i)-
     &              oldmasses(i,nlast))/
     &              (time-oldtimes(i,nlast))
            ELSE
               accrete(i) = 0.
            ENDIF
            acctime(i) = (time-
     &           oldtimes(i,nlast))/2.0 + 
     &           oldtimes(i,nlast)
c
c--Accretion rate of 0.0427 was used until 11/2013
c--Accretion rate of 0.04713 is 10^-7 M_sol/yr
c
         IF(endacctime(i).EQ.0.AND.accrete(i).LT.0.04713*1.)THEN
                endacctime(i) = acctime(i)
                WRITE (*,*) 'End acc time ',i, endacctime(i)
            ELSEIF (accrete(i).GE.0.04713*1.) THEN
               endacctime(i) = 0.0
            ENDIF
            IF (time-oldtimes(i,nlast).NE.0.0) THEN
               accel(i) = SQRT((vx(i)-oldvx(i,nlast))**2 +
     &              (vy(i)-oldvy(i,nlast))**2 +
     &              (vz(i)-oldvz(i,nlast))**2)/
     &              (time-oldtimes(i,nlast))
            ELSE
               accel(i) = 0.
            ENDIF
         ELSE
            accrete(i) = 0.0
            acctime(i) = 0.0
            accel(i) = 0.0
         ENDIF

         criticalvalue = 2271.6*2.0
c 2271.6 equals 1000 km/s/Myr acceleration
c         IF (accel(i).GT.2271.6) THEN
c         IF (accel(i).GT.2271.6*5.0) THEN

         IF (accel(i).GT.criticalvalue) THEN
            IF (ejectime(i).GT.0.0 .AND. 
     &           (time-ejectime(i).GT.0.00424339)) THEN
c Note that 0.00424339 is 2000 yr
               ejectimelast(i) = ejectime(i)
               WRITE (*,*) 'Saved ejtime ',i,ejectimelast(i)
            ENDIF
            ejectime(i) = time
c            WRITE (*,*) 'Eject time(s) ',i, ejectime(i),endacctime(i)
         ENDIF
 
      END DO
c
c--If merger identified through iunique
c
      nptmasscurrentnew = nptmasscurrent
c      IF (ifindmerge.EQ.1) THEN
         DO imerged = 1, nptmasscurrent
            IF (.NOT.ifound(imerged)) THEN
               WRITE (*,*) 'Merged: ifound ',imerged,listmerge(imerged),
     &              listunique(imerged),pmass(listmerge(imerged))
               WRITE (22,*) 'Merged: ifound ',imerged,
     &              listmerge(imerged),listunique(imerged),
     &              pmass(listmerge(imerged))
               nmerge = nmerge + 1
               pmass(listmerge(imerged)) = -1.0
               DO i = imerged, idim-1
                  listmerge(i) = listmerge(i+1)
                  listunique(i) = listunique(i+1)
               END DO
               nptmasscurrentnew = nptmasscurrentnew - 1

               GOTO 876
            END IF
         END DO
c      ENDIF
 876  CONTINUE
      nptmasscurrent = nptmasscurrentnew
c
c--Compute minimum radius between each two particles
c     Also, compute the centre of mass velocity of the sink particles
c
      vcmxtot = 0.
      vcmytot = 0.
      vcmztot = 0.
      DO k = 1, nptmass
         n = listmerge(k)
         radius2 = 1.0E+30
         vcmxtot = vcmxtot + pmass(n)*vx(n)
         vcmytot = vcmytot + pmass(n)*vy(n)
         vcmztot = vcmztot + pmass(n)*vz(n)
         DO kk = 1, nptmass
            n2 = listmerge(kk)
            IF (n.NE.n2) THEN
               rx = x(n)-x(n2)
               ry = y(n)-y(n2)
               rz = z(n)-z(n2)
               r2 = (rx**2 + ry**2 + rz**2)
               IF (r2.LT.radius2) THEN
                  radius2 = MIN(radius2, r2)
                  nearest = n2
                  pmassnearest = pmass(n2)
               ENDIF
               radiusend(n,n2) = SQRT(r2)
               IF (k.GT.nptmassold) THEN
                  radiuscreate(n,n2) = radiusend(n,n2)
               ENDIF
            ENDIF
         END DO
         IF (nlines .NE. 1) THEN
            IF (k.GT.nptmassold) THEN
               radiusmin(n)=1.0E+15
               write (*,*) 'Sep ',SQRT(radius2)*runit,nearest,
     &              pmassnearest
            ENDIF

            IF (radius2 .LT. radiusmin(n)**2) THEN
               radiusmin(n) = SQRT(radius2)
            ENDIF
         ELSE
            radiusmin(n) = 1.0E+15
         ENDIF
      END DO
      vcmxtot = vcmxtot/nptmass
      vcmytot = vcmytot/nptmass
      vcmztot = vcmztot/nptmass
c
c--Read next time
c
      GOTO 100
c
c--End of P-file
c
 200  CLOSE(15)
      CLOSE(22)
      nlines = nlines - 1 

      WRITE (*,*) 'Total number of lines was ',nlines

      DO i = 1, nbinimf 
         nimf(i) = 0
         nimfa(i) = 0
      END DO
      DO i = 1, nbinbin
         nbinary(i) = 0
         nbinary_notrealbd(i) = 0
         nbinary_01_02(i) = 0
         nbinary_02_05(i) = 0
         nbinary_05_08(i) = 0
         nbinary_08_12(i) = 0
         nbinary_12(i) = 0
         nbinary_01_02_notrealbd(i) = 0
         nbinary_02_05_notrealbd(i) = 0
         nbinary_05_08_notrealbd(i) = 0
         nbinary_08_12_notrealbd(i) = 0
         nbinary_12_notrealbd(i) = 0
         ntriple(i) = 0
         ntriplebd(i) = 0
         nquad(i) = 0
         ntriple_01_02(i) = 0
         ntriple_02_05(i) = 0
         ntriple_05_08(i) = 0
         ntriple_08_12(i) = 0
         ntriple_12(i) = 0
         nbinarybd(i) = 0
         nbinarybdobs(i) = 0
      END DO
      DO i = 1, nmassratio
         nbinaryq(i) = 0
         nbinaryq_notrealbd(i) = 0
         nbinaryq_01_02(i) = 0
         nbinaryq_02_05(i) = 0
         nbinaryq_05_08(i) = 0
         nbinaryq_08_12(i) = 0
         nbinaryq_12(i) = 0
         nbinaryq_01_02_notrealbd(i) = 0
         nbinaryq_02_05_notrealbd(i) = 0
         nbinaryq_05_08_notrealbd(i) = 0
         nbinaryq_08_12_notrealbd(i) = 0
         nbinaryq_12_notrealbd(i) = 0
         ntripleq(i) = 0
         nbinarybdq(i) = 0
         nbinarybdobsq(i) = 0
      END DO
      uvel = 0.20748
      write (*,*) 'Centre of mass velocity ',
     &     sqrt(vcmxtot**2+vcmytot**2+vcmztot**2)*uvel
c
c--Compute velocity of each sink particle relative to centre of mass
c     and IMF
c
      DO k = 1, nptmass
         n = listmerge(k)
c
c--If two ejections have occurred separated by more than 2000 yrs, then
c     choose the one closest to the end of accretion time (but set negative
c     if the second-to-last ejection was taken).
c
         IF (ABS(ejectimelast(n)-endacctime(n)).LT.
     &        ABS(ejectime(n)-endacctime(n))) THEN
            ejectime(n) = -ejectimelast(n)
         ENDIF

         IF (endacctime(n).EQ.0.0) endacctime(n) = time
         IF (ejectime(n).EQ.0.0) ejectime(n) = time

         vrelx = vx(n) - vcmxtot
         vrely = vy(n) - vcmytot
         vrelz = vz(n) - vcmztot
         vrel(n) = sqrt(vrelx**2+vrely**2+vrelz**2)*uvel
         vrelb(n) = vrel(n)
         write (*,*) 'Velocity ',n,': ',vrel(n)

         ipos = INT((ALOG10(pmass(n))-1)*5.0+nbinimf)
         IF (ipos.GE.1 .AND. ipos.LE.nbinimf) THEN
            nimf(ipos) = nimf(ipos) + 1
            IF (accrete(n).NE.0.0) 
     &        nimfa(ipos) = nimfa(ipos) + 1
         ENDIF
      END DO
c
c--For each point mass, output evolution file
c
c      DO n = 1, nptmass + nmerge
c         IF (n.GT.99) THEN
c            WRITE (num, 88000) n
c88000       FORMAT (I3)
c         ELSEIF (n.GT.9) THEN
c            WRITE (num, 88001) n
c88001       FORMAT ('0',I2)
c         ELSE 
c            WRITE (num, 88002) n
c88002       FORMAT ('00',I1)
c         ENDIF
c         DO junk = LEN(outbase), 1, -1
c            IF (outbase(junk:junk).NE.' ') THEN
c               iplace = junk
c               GOTO 150
c            ENDIF
c         END DO
c 150     outfile = outbase(1:junk) // num
c
c         IF (iout.EQ.'Y' .OR. iout.EQ.'y') THEN
c            OPEN (15, FILE=outfile)
c            DO i = 1, nlines
c               IF (pmass(n,i).GT.0.) THEN
c                  WRITE(15, 99005) time(i), fftime(i), 
c     &                 x(n,i), y(n,i), z(n,i), vx(n,i), 
c     &           vy(n,i), vz(n,i), pmass(n,i), rho(n,i), acctime(n,i),
c     &           accrete(n,i), spinx(n,i), spiny(n,i), spinz(n,i),
c     &           radiusmin(n,i)*runit,sqrt(vx(n,i)**2+vy(n,i)**2+
c     &           vz(n,i)**2),accel(n,i)
c99005             FORMAT (18(1PE15.8,1X))
c               ENDIF
c            END DO
c            CLOSE (15)
c         ENDIF
c         DO i = nlines, 1, -1
cc            IF (accel(n,i).GT.2271.6) THEN
cc            IF (accel(n,i).GT.2271.6*5.0) THEN
c            IF (accel(n,i).GT.2271.6*2.0) THEN
cc 2271.6 equals 1000 km/s/Myr acceleration
c               ejectime(n) = acctime(n,i)
c               WRITE (*,*) 'Eject time ',n, ejectime(n),endacctime(n)
c               GOTO 151
c            ENDIF
c         END DO
c 151     CONTINUE
c      END DO
c
c--For last dump, output endfile containing data for each point mass
c
      DO k = 1, nptmass
         n = listmerge(k)
         invert(n) = k
         WRITE (88,99555) tffform(n),tform(n),k,listunique(k)
99555    FORMAT(2(1PE12.5,1X),I5,1X,I11)
      END DO

      OPEN (15,FILE='endfile')
      WRITE (15,99100) time,fftime
99100 FORMAT('# ',2(1PE15.8,1X))
      centreofmassx = 0.
      centreofmassy = 0.
      centreofmassz = 0.
      totalmass = 0.
      xmaximummass = 0.
      DO k = 1, nptmass
         n = listmerge(k)
         IF (pmass(n).GT.xmaximummass) THEN
            xmaximummass = pmass(n)
            xmaximummassx = x(n)
            xmaximummassy = y(n)
            xmaximummassz = z(n)
         ENDIF
         centreofmassx = centreofmassx + pmass(n)*x(n)
         centreofmassy = centreofmassy + pmass(n)*y(n)
         centreofmassz = centreofmassz + pmass(n)*z(n)
         totalmass = totalmass + pmass(n)
         WRITE(15, 99012) tffform(n),tform(n),
     &        x(n), y(n), z(n), vx(n),
     &        vy(n), vz(n), pmass(n), 
     &        accrete(n), spinx(n), spiny(n), spinz(n),
     &        radiusmin(n)*runit,vrel(n),endacctime(n),
     &        endacctime(n)*tform(n)/tffform(n),
     &        endacctime(n)-tffform(n),
     &        endacctime(n)*tform(n)/tffform(n)-tform(n),
     &        ejectime(n),ejectime(n)*tform(n)/tffform(n),
     &        ABS(ejectime(n))-tffform(n),
     &        ABS(ejectime(n))*tform(n)/tffform(n)-tform(n)
99012    FORMAT (23(1PE15.8,1X))
         WRITE (33,99333) n,pmass(n),endacctime(n)*tform(n)/
     &        tffform(n)-tform(n),
     &        ABS(ejectime(n))*tform(n)/tffform(n)-tform(n)
99333    FORMAT (I4,1X,3(1PE15.8,1X))
      END DO
      CLOSE(15)
      centreofmassx = centreofmassx/totalmass
      centreofmassy = centreofmassy/totalmass
      centreofmassz = centreofmassz/totalmass

      DO k = 1, nptmass
         n = listmerge(k)
         radcom(n) = SQRT((x(n)-centreofmassx)**2 + 
     &        (y(n)-centreofmassy)**2 + (z(n)-centreofmassz)**2)
         radmax(n) = SQRT((x(n)-xmaximummassx)**2 + 
     &        (y(n)-xmaximummassy)**2 + (z(n)-xmaximummassz)**2)
      END DO

      OPEN (15,FILE='radiuscreate')
      DO k = 1, nptmass
         n = listmerge(k)
         WRITE(15, 99014) (radiuscreate(n,i)*runit,i=1,nptmass)
99014    FORMAT (51(1PE15.8,1X))
      END DO
      CLOSE(15)
      OPEN (15,FILE='radiusend')
      DO k = 1, nptmass
         n = listmerge(k)
         sepmin = 1.0E+10
         nkeep = 0
         DO kk = 1, n-1
            n2 = listmerge(kk)
            IF (radiuscreate(n,n2).LT.sepmin) THEN
               sepmin = radiuscreate(n,n2)
               nkeep = n2
            ENDIF
         END DO
         IF (nkeep.NE.0) THEN
            WRITE(15, *) sepmin*runit, radiusend(n,nkeep)*runit
         ENDIF
      END DO
      CLOSE(15)
c
c--Make IMF
c
      DO i = 1, nbinimf
         write (*,*) 10**((i-nbinimf)/5.0+1.0), nimf(i), nimfa(i)
      END DO
c
c--Calculate binaries and multiple systems
c
      numberbd = 0
      numberbdobs = 0
      numberstar = 0
      numberbd3 = 0
      numberbd4 = 0
      numberbd345 = 0
      numberbdobs3 = 0
      numberbdobs4 = 0
      numberbdobs345 = 0
      numberstar3 = 0
      numberstar4 = 0
      numberstar345 = 0
      numbermixed3 = 0
      numbermixed4 = 0
      number2 = 0
      number2_11 = 0
      number2_08 = 0
      number2_05 = 0
      number2_02 = 0
      number2_01 = 0
      number2_007 = 0
      number2_003 = 0
      number2_001 = 0
      number2_0003 = 0
      number2_11notbd = 0
      number2_08notbd = 0
      number2_05notbd = 0
      number2_02notbd = 0
      number2_01notbd = 0
      number_11 = 0
      number_08 = 0
      number_05 = 0
      number_02 = 0
      number_01 = 0
      number_007 = 0
      number_003 = 0
      number_001 = 0
      number_0003 = 0
      number_11notbd = 0
      number_08notbd = 0
      number_05notbd = 0
      number_02notbd = 0
      number_01notbd = 0
      number3 = 0
      number3_11 = 0
      number3_08 = 0
      number3_05 = 0
      number3_02 = 0
      number3_01 = 0
      number3_007 = 0
      number3_003 = 0
      number3_001 = 0
      number3_0003 = 0
      number3_11notbd = 0
      number3_08notbd = 0
      number3_05notbd = 0
      number3_02notbd = 0
      number3_01notbd = 0
      number4 = 0
      number4_11 = 0
      number4_08 = 0
      number4_05 = 0
      number4_02 = 0
      number4_01 = 0
      number4_007 = 0
      number4_003 = 0
      number4_001 = 0
      number4_0003 = 0
      number4_11notbd = 0
      number4_08notbd = 0
      number4_05notbd = 0
      number4_02notbd = 0
      number4_01notbd = 0
      number5 = 0
      number2unique = 0
      number3unique = 0
      number4unique = 0
      number5unique = 0
      number2isolated = 0
      number3isolated = 0
      number4isolated = 0
      number5isolated = 0
      number2bd = 0
      number2bdobs = 0
      number3bd = 0
      number4bd = 0
      number5bd = 0
      number2star = 0
      numberstarbd = 0
      DO k = 1, nptmass
         n = listmerge(k)
         xn(k) = x(n)
         yn(k) = y(n)
         zn(k) = z(n)
         vxn(k) = vx(n)
         vyn(k) = vy(n)
         vzn(k) = vz(n)
         spinx(k) = spinx(n)
         spiny(k) = spiny(n)
         spinz(k) = spinz(n)
         pmassn(k) = pmass(n)
         pmassmax(k) = pmassn(k)
         pmassmin(k) = pmassn(k)
         IF (endacctime(k).NE.time) interest(k)=interest(k)+10
         IF (pmassn(k).LT.bdmass) THEN
            interest(k) = interest(k) + 1000
            IF (pmassn(k).GT.0.03) THEN
               numberbdobs = numberbdobs + 1
               interest(k) = interest(k) - 900
            ENDIF
            numberbd= numberbd + 1
         ELSE
            numberstar= numberstar + 1
         ENDIF
         IF (pmassn(k).LT.0.01) THEN
            number_0003 = number_0003 + 1
         ELSEIF (pmassn(k).GE.1.2) THEN
            number_11 = number_11 + 1
            number_11notbd = number_11notbd + 1
         ELSEIF (pmassn(k).GE.0.8) THEN
            number_08 = number_08 + 1
            number_08notbd = number_08notbd + 1
         ELSEIF (pmassn(k).GE.0.5) THEN
            number_05 = number_05 + 1
            number_05notbd = number_05notbd + 1
         ELSEIF (pmassn(k).GE.0.2) THEN
            number_02 = number_02 + 1
            number_02notbd = number_02notbd + 1
         ELSEIF (pmassn(k).GE.0.1) THEN
            number_01 = number_01 + 1
            number_01notbd = number_01notbd + 1
         ELSEIF (pmassn(k).GE.0.07) THEN
            number_007 = number_007 + 1
         ELSEIF (pmassn(k).GE.0.03) THEN
            number_003 = number_003 + 1
         ELSEIF (pmassn(k).GE.0.01) THEN
            number_001 = number_001 + 1
         ENDIF
c         write (tempstring,77001) n
         write (tempstring,77001) invert(n)
         name(k) = tempstring
         lname(k) = 4
         nbincomp(k) = 1
77001    FORMAT(I4)
      END DO
      nnodes = nptmass
      numberbdsingle = numberbd
      numberbdobssingle = numberbdobs
      numberstarsingle = numberstar

      radius2 = 1.0E+30
 1500 emin = 1.0E+10
      rmin = 1.0E+10
      DO n = 1, nnodes
         DO n2 = n+1, nnodes
            rx = xn(n)-xn(n2)
            ry = yn(n)-yn(n2)
            rz = zn(n)-zn(n2)
            vxdiff = vxn(n)-vxn(n2)
            vydiff = vyn(n)-vyn(n2)
            vzdiff = vzn(n)-vzn(n2)
            r2 = (rx**2 + ry**2 + rz**2)
            radius = sqrt(r2)
            radius2 = MIN(radius2, r2)
            pmasstot = pmassn(n) + pmassn(n2)
            ekin = 0.5*(vxdiff**2 + vydiff**2 + vzdiff**2)
            epot = -pmasstot/sqrt(r2)
            etot = ekin + epot
            dangx = vydiff*rz - vzdiff*ry
            dangy = vzdiff*rx - vxdiff*rz
            dangz = vxdiff*ry - vydiff*rx
            dang2 = (dangx**2 + dangy**2 + dangz**2)
            
            semimajor = - pmasstot/2.0/etot
            eccentricity = SQRT(1.0 - dang2/pmasstot/semimajor)

cc            write (*,*) n,n2,etot,ekin,epot
c
c--Select most bound
c
c            IF (etot.LT.emin) THEN
c               emin = MIN(emin,etot)
c               nkeep1 = n
c               nkeep2 = n2
c               semimajorkeep = semimajor
c               eccentricitykeep = eccentricity
c            ENDIF
c
c--Select closest
c
            IF (etot.LT.0.0 .AND. radius.LT.rmin
     &           .AND. nbincomp(n)+nbincomp(n2).LT.5
     &           ) THEN
c
c--Select *mutual* nearest neighbours
c
               radius_test2_keep = 1.0E+30
               DO n3 = 1, nnodes
                  IF (n3.NE.n2) THEN
                     radius_test2 = (xn(n2)-xn(n3))**2 + 
     &                    (yn(n2)-yn(n3))**2 + (zn(n2)-zn(n3))**2
                     IF (radius_test2.LT.radius_test2_keep) THEN
                        radius_test2_keep = radius_test2
                        n23_keep = n3
                     ENDIF
                  ENDIF
               END DO
 
               radius_test2_keep = 1.0E+30
               DO n3 = 1, nnodes
                  IF (n3.NE.n) THEN
                     radius_test2 = (xn(n)-xn(n3))**2 + 
     &                    (yn(n)-yn(n3))**2 + (zn(n)-zn(n3))**2
                     IF (radius_test2.LT.radius_test2_keep) THEN
                        radius_test2_keep = radius_test2
                        n13_keep = n3
                     ENDIF
                  ENDIF
               END DO

c               print *,'Check is ',radius_test2_keep,n23_keep,n13_keep,
c     &              n,n2
c               print *,'Values ',pmassn(n3_keep),pmassn(n),pmassn(n2)
c               print *,'xyz-15 ',xn(15),yn(15),zn(15)
c               print *,'xyz-129 ',xn(129),yn(129),zn(129)
c               print *,' '

               IF (n23_keep.EQ.n .AND. n13_keep.EQ.n2) THEN
                  rmin = MIN(rmin,radius)
                  emin = etot
                  nkeep1 = n
                  nkeep2 = n2
                  semimajorkeep = semimajor
                  eccentricitykeep = eccentricity
c
c--Orbital plane
c
                  planex = -(vydiff*rz - vzdiff*ry)
                  planey = -(vzdiff*rx - vxdiff*rz)
                  planez = -(vxdiff*ry - vydiff*rx)
                  planelength = SQRT(planex**2 + planey**2 + planez**2)
                  planex = planex/planelength
                  planey = planey/planelength
                  planez = planez/planelength
               ENDIF
            ENDIF
            
         END DO
ccc         WRITE (99,*) x(n),y(n),z(n)
      END DO
      WRITE (*,*) 'emin ',emin
      WRITE (*,*) 'rmin ',rmin
      IF (emin.LT.0.0) THEN
         qratio = MIN(pmassn(nkeep1)/pmassn(nkeep2), 
     &        pmassn(nkeep2)/pmassn(nkeep1))
         periodkeep = 2.0*pi*SQRT(semimajorkeep**3/(6.67E-08*
     &        (pmassn(nkeep1)+pmassn(nkeep2))))*(1.496E+13*runit)**1.5/
     &        SQRT(1.991E+33)/3.1558118E+7

         xmassmax = MAX(pmassmax(nkeep1),
     &        pmassmax(nkeep2))
         xmassmin = MIN(pmassmin(nkeep1),
     &        pmassmin(nkeep2))
         WRITE (*,*) name(nkeep1), ' is bound to ',name(nkeep2)
         WRITE (*,*) '      Semi-major axis: ',semimajorkeep*runit
         WRITE (*,*) '      Period (yr)      ',periodkeep
         WRITE (*,*) '      Eccentricity:    ',eccentricitykeep
         WRITE (*,*) '      Mass Ratio:      ',qratio
         WRITE (*,*) '      Masses:          ',pmassn(nkeep1),
     &        pmassn(nkeep2)
         WRITE (*,*) '      Max Mass:        ',xmassmax
         WRITE (*,*) '      Min Mass:        ',xmassmin 
         WRITE (*,*) '      Spin1:        ',spinx(nkeep1),
     &        spiny(nkeep1),spinz(nkeep1)
         WRITE (*,*) '      Spin2:        ',spinx(nkeep2),
     &        spiny(nkeep2),spinz(nkeep2)
         WRITE (*,*) '      Orbital plane    ',planex,planey,planez
         WRITE (*,*) '      Interest values  ',interest(nkeep1),
     &        interest(nkeep2)
c
c--Merge binaries into 1 node
c
         nbincomp1old = nbincomp(nkeep1)
         nbincomp(nkeep1) = nbincomp(nkeep1) + nbincomp(nkeep2)
c
c--Keep track of maximum order of multiple (e.g. if a binary is part of quad)
c
         imultiple_values(4,noutputline(nkeep1)) = nbincomp(nkeep1)
         imultiple_values(4,noutputline(nkeep2)) = nbincomp(nkeep1)

         angleorbit1 = 0.
         angleorbit2 = 0.
         anglerelative = 0.
         IF (nbincomp(nkeep1).LE.4) THEN
            IF (nbincomp(nkeep1).EQ.3) THEN
               IF (nbincomp1old.EQ.2) THEN
                  qratio_out = pmassn(nkeep2)/pmassn(nkeep1)
                  angleorbit = ACOS(planex*orbit(1,nkeep1)+
     &                 planey*orbit(2,nkeep1)+
     &                 planez*orbit(3,nkeep1))*180.0/pi
                  WRITE (*,*) 'Triple ',planex,planey,planez
                  WRITE (*,*) 'Binary ',orbit(1,nkeep1),
     &                 orbit(2,nkeep1),orbit(3,nkeep1)
                  WRITE (*,*) 'Angle Orbit ',angleorbit
                  periodratiokeep = periodkeep/period(nkeep1)
                  WRITE (*,*) 'Period Ratio1 ',periodratiokeep,
     &                 period(nkeep1),nkeep1
                  WRITE (81,*) period(nkeep1)

                  orbit1length = SQRT(orbit(1,nkeep1)**2 + 
     &                 orbit(2,nkeep1)**2 + orbit(3,nkeep1)**2)
                  orbit2length = SQRT(planex**2 + 
     &                 planey**2 + planez**2)

                  orbit1incl = ACOS(orbit(3,nkeep1)/orbit1length)
                  orbit1omega1 = ASIN(orbit(1,nkeep1)/orbit1length/
     &                 SIN(orbit1incl))
                  orbit1omega2 = ACOS(-orbit(2,nkeep1)/orbit1length/
     &                 SIN(orbit1incl))

                  orbit2incl = ACOS(planez/orbit2length)
                  orbit2omega1 = ASIN(planex/orbit2length/
     &                 SIN(orbit2incl))
                  orbit2omega2 = ACOS(-planey/orbit2length/
     &                 SIN(orbit2incl))

                  supposedx = orbit1length*SIN(orbit1incl)*
     &                 SIN(orbit1omega1)
                  supposedy = - orbit1length*SIN(orbit1incl)*
     &                 COS(orbit1omega1)
                  supposedz = orbit1length*COS(orbit1incl)

                  IF(ABS((supposedx-orbit(1,nkeep1))/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-orbit(2,nkeep1))/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit1omega = orbit1omega1

                     WRITE (83,5) periodratiokeep,orbit(1,nkeep1),
     &                    orbit(2,nkeep1),orbit(3,nkeep1),
     &                    supposedx,supposedy,supposedz
 5                   FORMAT('1 ',7(1PE12.5,1X))
                  ELSEIF(ABS((-supposedx-orbit(1,nkeep1))/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-orbit(2,nkeep1))/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit1omega = - orbit1omega1

                     WRITE (83,11) periodratiokeep,orbit(1,nkeep1),
     &                    orbit(2,nkeep1),orbit(3,nkeep1),
     &                    supposedx,supposedy,supposedz
 11                  FORMAT('1a ',7(1PE12.5,1X))
                  ELSE

                     supposedx = orbit1length*SIN(orbit1incl)*
     &                    SIN(orbit1omega2)
                     supposedy = - orbit1length*SIN(orbit1incl)*
     &                    COS(orbit1omega2)

                     IF (orbit(1,nkeep1)*supposedx.GE.0.0) THEN
                        orbit1omega = orbit1omega2
                     ELSE
                        orbit1omega = - orbit1omega2
                     ENDIF

                     WRITE (83,6) periodratiokeep,orbit(1,nkeep1),
     &                    orbit(2,nkeep1),orbit(3,nkeep1),
     &                    supposedx,supposedy,supposedz
 6                   FORMAT('2 ',7(1PE12.5,1X))

                  ENDIF

                  supposedx = orbit2length*SIN(orbit2incl)*
     &                 SIN(orbit2omega1)
                  supposedy = - orbit2length*SIN(orbit2incl)*
     &                 COS(orbit2omega1)
                  supposedz = orbit2length*COS(orbit2incl)

                  IF(ABS((supposedx-planex)/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-planey)/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit2omega = orbit2omega1

                     WRITE (83,7) periodratiokeep,planex,
     &                    planey,planez,supposedx,supposedy,supposedz
 7                   FORMAT('3 ',7(1PE12.5,1X))
                  ELSEIF(ABS((-supposedx-planex)/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-planey)/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit2omega = -orbit2omega1

                     WRITE (83,10) periodratiokeep,planex,
     &                    planey,planez,supposedx,supposedy,supposedz
 10                  FORMAT('3a ',7(1PE12.5,1X))
                  ELSE
                     supposedx = orbit2length*SIN(orbit2incl)*
     &                    SIN(orbit2omega2)
                     supposedy = - orbit2length*SIN(orbit2incl)*
     &                    COS(orbit2omega2)

                     IF (planex*supposedx.GE.0.0) THEN
                        orbit2omega = orbit2omega2
                     ELSE
                        orbit2omega = -orbit2omega2
                     ENDIF

                     WRITE (83,8) periodratiokeep,planex,
     &                    planey,planez,supposedx,supposedy,supposedz
 8                   FORMAT('4 ',7(1PE12.5,1X))
                  ENDIF

                  ambiguous = COS(orbit1omega-orbit2omega)
                  totalangle1 = ACOS(COS(orbit1incl)*COS(orbit2incl) +
     &                 SIN(orbit1incl)*SIN(orbit2incl)*ambiguous)
                  totalangle2 = ACOS(COS(orbit1incl)*COS(orbit2incl) -
     &                 SIN(orbit1incl)*SIN(orbit2incl)*ambiguous)
                  WRITE (82,88009) periodratiokeep,orbit1incl*180./pi,
     &                 orbit1omega1*180./pi,
     &                 orbit1omega2*180./pi,orbit2incl*180./pi,
     &                 orbit2omega1*180./pi,orbit2omega2*180./pi,
     &                 orbit1omega*180./pi,orbit2omega*180./pi,
     &                 totalangle1*180/pi,totalangle2*180/pi
88009             FORMAT(11(1PE12.5,1X))
               ELSE
                  qratio_out = pmassn(nkeep1)/pmassn(nkeep2)
                  angleorbit = ACOS(planex*orbit(1,nkeep2)+
     &                 planey*orbit(2,nkeep2)+
     &                 planez*orbit(3,nkeep2))*180.0/pi
                  WRITE (*,*) 'Triple ',planex,planey,planez
                  WRITE (*,*) 'Binary ',orbit(1,nkeep2),
     &                 orbit(2,nkeep2),orbit(3,nkeep2)
                  WRITE (*,*) 'Angle Orbit ',angleorbit
                  periodratiokeep = periodkeep/period(nkeep2)
                  WRITE (*,*) 'Period Ratio2 ',periodratiokeep,
     &                 period(nkeep2),nkeep2
                  WRITE (81,*) period(nkeep2)

                  orbit1length = SQRT(orbit(1,nkeep2)**2 + 
     &                 orbit(2,nkeep2)**2 + orbit(3,nkeep2)**2)
                  orbit2length = SQRT(planex**2 + 
     &                 planey**2 + planez**2)

                  orbit1incl = ACOS(orbit(3,nkeep2)/orbit1length)
                  orbit1omega1 = ASIN(orbit(1,nkeep2)/orbit1length/
     &                 SIN(orbit1incl))
                  orbit1omega2 = ACOS(-orbit(2,nkeep2)/orbit1length/
     &                 SIN(orbit1incl))

                  orbit2incl = ACOS(planez/orbit2length)
                  orbit2omega1 = ASIN(planex/orbit2length/
     &                 SIN(orbit2incl))
                  orbit2omega2 = ACOS(-planey/orbit2length/
     &                 SIN(orbit2incl))

                  supposedx = orbit1length*SIN(orbit1incl)*
     &                 SIN(orbit1omega1)
                  supposedy = - orbit1length*SIN(orbit1incl)*
     &                 COS(orbit1omega1)
                  supposedz = orbit1length*COS(orbit1incl)

                  IF(ABS((supposedx-orbit(1,nkeep2))/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-orbit(2,nkeep2))/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit1omega = orbit1omega1

                     WRITE (83,15) periodratiokeep,orbit(1,nkeep2),
     &                    orbit(2,nkeep2),orbit(3,nkeep2),
     &                    supposedx,supposedy,supposedz
 15                   FORMAT('1b ',7(1PE12.5,1X))
                  ELSEIF(ABS((-supposedx-orbit(1,nkeep2))/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-orbit(2,nkeep2))/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit1omega = - orbit1omega1

                     WRITE (83,21) periodratiokeep,orbit(1,nkeep2),
     &                    orbit(2,nkeep2),orbit(3,nkeep2),
     &                    supposedx,supposedy,supposedz
 21                  FORMAT('1ba ',7(1PE12.5,1X))
                  ELSE

                     supposedx = orbit1length*SIN(orbit1incl)*
     &                    SIN(orbit1omega2)
                     supposedy = - orbit1length*SIN(orbit1incl)*
     &                    COS(orbit1omega2)

                     IF (orbit(1,nkeep2)*supposedx.GE.0.0) THEN
                        orbit1omega = orbit1omega2
                     ELSE
                        orbit1omega = - orbit1omega2
                     ENDIF

                     WRITE (83,16) periodratiokeep,orbit(1,nkeep2),
     &                    orbit(2,nkeep2),orbit(3,nkeep2),
     &                    supposedx,supposedy,supposedz
 16                   FORMAT('2b ',7(1PE12.5,1X))

                  ENDIF

                  supposedx = orbit2length*SIN(orbit2incl)*
     &                 SIN(orbit2omega1)
                  supposedy = - orbit2length*SIN(orbit2incl)*
     &                 COS(orbit2omega1)
                  supposedz = orbit2length*COS(orbit2incl)

                  IF(ABS((supposedx-planex)/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-planey)/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit2omega = orbit2omega1

                     WRITE (83,17) periodratiokeep,planex,
     &                    planey,planez,supposedx,supposedy,supposedz
 17                   FORMAT('3b ',7(1PE12.5,1X))
                  ELSEIF(ABS((-supposedx-planex)/supposedx).LT.
     &                 1.0e-4 .AND. 
     &                 ABS((supposedy-planey)/supposedy).LT.
     &                 1.0e-4) THEN
                     orbit2omega = -orbit2omega1

                     WRITE (83,20) periodratiokeep,planex,
     &                    planey,planez,supposedx,supposedy,supposedz
 20                  FORMAT('3ba ',7(1PE12.5,1X))
                  ELSE
                     supposedx = orbit2length*SIN(orbit2incl)*
     &                    SIN(orbit2omega2)
                     supposedy = - orbit2length*SIN(orbit2incl)*
     &                    COS(orbit2omega2)

                     IF (planex*supposedx.GE.0.0) THEN
                        orbit2omega = orbit2omega2
                     ELSE
                        orbit2omega = -orbit2omega2
                     ENDIF

                     WRITE (83,18) periodratiokeep,planex,
     &                    planey,planez,supposedx,supposedy,supposedz
 18                   FORMAT('4b ',7(1PE12.5,1X))
                  ENDIF

                  ambiguous = COS(orbit1omega-orbit2omega)
                  totalangle1 = ACOS(COS(orbit1incl)*COS(orbit2incl) +
     &                 SIN(orbit1incl)*SIN(orbit2incl)*ambiguous)
                  totalangle2 = ACOS(COS(orbit1incl)*COS(orbit2incl) -
     &                 SIN(orbit1incl)*SIN(orbit2incl)*ambiguous)
                  WRITE (82,88009) periodratiokeep,orbit1incl*180./pi,
     &                 orbit1omega1*180./pi,
     &                 orbit1omega2*180./pi,orbit2incl*180./pi,
     &                 orbit2omega1*180./pi,orbit2omega2*180./pi,
     &                 orbit1omega*180./pi,orbit2omega*180./pi,
     &                 totalangle1*180/pi,totalangle2*180/pi
               ENDIF
            ELSE
               qratio_out = qratio
               angleorbit = ACOS(planez)*180.0/pi
               WRITE (*,*) 'Angle Orbit ',angleorbit
               spinlength1 = SQRT(spinx(nkeep1)**2 + spiny(nkeep1)**2 +
     &              spinz(nkeep1)**2)
               spinlength2 = SQRT(spinx(nkeep2)**2 + spiny(nkeep2)**2 +
     &              spinz(nkeep2)**2)

               IF (nbincomp(nkeep1).EQ.2) THEN
                  angleorbit1 = ACOS((spinx(nkeep1)*planex+
     &                 spiny(nkeep1)*planey+
     &                 spinz(nkeep1)*planez)/spinlength1)*180.0/pi
                  angleorbit2 = ACOS((spinx(nkeep2)*planex+
     &                 spiny(nkeep2)*planey+
     &                 spinz(nkeep2)*planez)/spinlength2)*180.0/pi
                  anglerelative = ACOS((spinx(nkeep2)*spinx(nkeep1)+
     &                 spiny(nkeep2)*spiny(nkeep1)+
     &                 spinz(nkeep2)*spinz(nkeep1))/
     &                 spinlength1/spinlength2)*180.0/pi
                  WRITE (*,*) 'Angle Spin-Orbit 1 ',angleorbit1 
                  WRITE (*,*) 'Angle Spin-Orbit 2 ',angleorbit2
                  WRITE (*,*) 'Angle Spin-Spin ',anglerelative
                  periodratiokeep = 0.
               ELSE
                  periodratiokeep = periodkeep/MAX(period(nkeep1),
     &                 period(nkeep2))
                  WRITE (*,*) 'Period Ratio ',periodratiokeep
                  IF (nbincomp1old.EQ.2) THEN
                     WRITE (81,*) period(nkeep1)
                  ELSEIF (nbincomp(nkeep2).EQ.2) THEN
                     WRITE (81,*) period(nkeep1)
                  ENDIF
               ENDIF
            ENDIF
            
            nmultipleout = nmultipleout + 1

            imultiple_values(1,nmultipleout) = nbincomp(nkeep1)
            imultiple_values(2,nmultipleout) = interest(nkeep1)
            imultiple_values(3,nmultipleout) = interest(nkeep2)
            imultiple_values(4,nmultipleout) = nbincomp(nkeep1)
            noutputline(nkeep1) = nmultipleout

            xmultiple_values(1,nmultipleout) = xmassmax
            xmultiple_values(2,nmultipleout) = xmassmin
            xmultiple_values(3,nmultipleout) = pmassn(nkeep1)
            xmultiple_values(4,nmultipleout) = pmassn(nkeep2)
            xmultiple_values(5,nmultipleout) = qratio_out
            xmultiple_values(6,nmultipleout) = semimajorkeep*runit
            xmultiple_values(7,nmultipleout) = eccentricitykeep
            xmultiple_values(8,nmultipleout) = angleorbit
            xmultiple_values(9,nmultipleout) = angleorbit1
            xmultiple_values(10,nmultipleout) = angleorbit2
            xmultiple_values(11,nmultipleout) = anglerelative
            xmultiple_values(12,nmultipleout) = periodkeep
            xmultiple_values(13,nmultipleout) = periodratiokeep

c            WRITE (50,99444) nbincomp(nkeep1),xmassmax,xmassmin,
c     &           interest(nkeep1),interest(nkeep2),
c     &           pmassn(nkeep1),pmassn(nkeep2),
c     &           qratio_out,semimajorkeep*runit,eccentricitykeep,
c     &           angleorbit,angleorbit1,angleorbit2,anglerelative,
c     &           periodkeep,periodratiokeep
99444       FORMAT(I2,2(1X,1PE12.5),2(1X,I6),11(1X,1PE12.5),(1X,I2))
         ENDIF
         orbit(1,nkeep1) = planex
         orbit(2,nkeep1) = planey
         orbit(3,nkeep1) = planez
         period(nkeep1) = periodkeep
         periodratio(nkeep1) = periodratiokeep
         interest(nkeep1) = interest(nkeep1) + interest(nkeep2)
         IF (nbincomp(nkeep1).GT.2) THEN
            IF (nbincomp(nkeep2).GT.nbincomp1old) THEN
               seplast(nkeep1) = seplast(nkeep2)
               seplast2(nkeep1) = seplast2(nkeep2)
               qratiolast(nkeep1) = qratiolast(nkeep2)
            ENDIF
            seplastkeep = seplast(nkeep1)
            qratiolastkeep = qratiolast(nkeep1)
            seplast2keep = seplast2(nkeep1)
            seplast2(nkeep1) = seplastkeep
         ELSE
            seplastkeep = -1.
            qratiolastkeep = -1.
         ENDIF
         seplast(nkeep1) = semimajorkeep
         qratiolast(nkeep1) = qratio_out

         IF (nbincomp(nkeep1).EQ.2) THEN
            number2 = number2 + 1
            number2unique = number2unique + 1
            number2isolated = number2isolated + 1
            ipos = MAX(1,MIN(INT((ALOG10(semimajorkeep*runit))*5.0+16),
     &           nbinbin))
            IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
               nbinary(ipos) = nbinary(ipos) + 1
               IF (xmassmin.GT.bdrealmass) THEN
                  nbinary_notrealbd(ipos) = nbinary_notrealbd(ipos) + 1
               ENDIF
                     
               IF(pmassn(nkeep1).LT.bdmass.AND.pmassn(nkeep2).LT.bdmass)
     &              THEN
                  IF (pmassn(nkeep1).GT.0.03.OR.pmassn(nkeep2).GT.0.03)
     &                 THEN
                     number2bdobs = number2bdobs + 1
                     IF (pmassn(nkeep1).GT.0.03) THEN
                        numberbdobssingle = numberbdobssingle - 1
                     ENDIF
                     IF (pmassn(nkeep2).GT.0.03) THEN
                        numberbdobssingle = numberbdobssingle - 1
                     ENDIF
                     nbinarybdobs(ipos) = nbinarybdobs(ipos) + 1
                  ENDIF
                  number2bd = number2bd + 1
                  numberbdsingle = numberbdsingle - 2
                  nbinarybd(ipos) = nbinarybd(ipos) + 1
                  WRITE (*,*) 'BBD'
           ELSEIF (pmassn(nkeep1).LT.bdmass.OR.pmassn(nkeep2).LT.bdmass)
     &                 THEN
                  numberstarbd = numberstarbd + 1
                  numberbdsingle = numberbdsingle - 1
                  numberbdobssingle = numberbdobssingle - 1
               IF (pmassn(nkeep1).LT.0.03.OR.pmassn(nkeep2).LT.0.03)THEN
                     numberbdobssingle = numberbdobssingle + 1
                  ENDIF
                  numberstarsingle = numberstarsingle - 1
                  WRITE (*,*) 'S-BD'
               ELSE
                  number2star = number2star + 1
                  numberstarsingle = numberstarsingle - 2
                  WRITE (*,*) 'Binary'
               ENDIF
            ENDIF
            print *,' numberbdobssingle ',numberbdobssingle
            iposq = MAX(1,MIN(INT(qratio*nmassratio+1),nmassratio))
            nbinaryq(iposq) = nbinaryq(iposq) + 1
            IF (xmassmin.GT.bdrealmass) THEN
               nbinaryq_notrealbd(iposq) = nbinaryq_notrealbd(iposq) + 1
            ENDIF

            IF (iposq.EQ.1) WRITE (*,*) 'EXTREME Q'
            IF(pmassn(nkeep1).LT.bdmass .AND. pmassn(nkeep2).LT.bdmass)
     &           THEN
               IF(pmassn(nkeep1).GT.0.03 .OR. pmassn(nkeep2).GT.0.03)
     &              nbinarybdobsq(iposq) = nbinarybdobsq(iposq) + 1
               nbinarybdq(iposq) = nbinarybdq(iposq) + 1
            ENDIF
            IF (xmassmax.GE.0.1) THEN
               IF (xmassmax.LT.0.2) THEN
                  nbinary_01_02(ipos) = nbinary_01_02(ipos) + 1

                  print *,'01_02 inc ',nbinary_01_02(ipos),ipos

                  nbinaryq_01_02(iposq) = nbinaryq_01_02(iposq) + 1
               ELSEIF (xmassmax.LT.0.5) THEN
                  nbinary_02_05(ipos) = nbinary_02_05(ipos) + 1
                  nbinaryq_02_05(iposq) = nbinaryq_02_05(iposq) + 1

                  print *,'02_05 inc ',nbinary_02_05(ipos),ipos

               ELSEIF (xmassmax.LT.0.8) THEN
                  nbinary_05_08(ipos) = nbinary_05_08(ipos) + 1
                  nbinaryq_05_08(iposq) = nbinaryq_05_08(iposq) + 1
               ELSEIF (xmassmax.LT.1.2) THEN
                  nbinary_08_12(ipos) = nbinary_08_12(ipos) + 1
                  nbinaryq_08_12(iposq) = nbinaryq_08_12(iposq) + 1
               ELSEIF (xmassmax.GE.1.2) THEN
                  nbinary_12(ipos) = nbinary_12(ipos) + 1
                  nbinaryq_12(iposq) = nbinaryq_12(iposq) + 1
               ENDIF

               IF (xmassmin.GT.bdrealmass) THEN
                  IF (xmassmax.LT.0.2) THEN
                     nbinary_01_02_notrealbd(ipos) = 
     &                    nbinary_01_02_notrealbd(ipos) + 1

        print *,'01_02-bd inc ',nbinary_01_02_notrealbd(ipos),ipos

                     nbinaryq_01_02_notrealbd(iposq) = 
     &                    nbinaryq_01_02_notrealbd(iposq) + 1
                  ELSEIF (xmassmax.LT.0.5) THEN
                     nbinary_02_05_notrealbd(ipos) = 
     &                    nbinary_02_05_notrealbd(ipos) + 1

        print *,'02_05-bd inc ',nbinary_02_05_notrealbd(ipos),ipos

                     nbinaryq_02_05_notrealbd(iposq) = 
     &                    nbinaryq_02_05_notrealbd(iposq) + 1
                  ELSEIF (xmassmax.LT.0.8) THEN
                     nbinary_05_08_notrealbd(ipos) = 
     &                    nbinary_05_08_notrealbd(ipos) + 1

        print *,'05_08-bd inc ',nbinary_05_08_notrealbd(ipos),ipos

                     nbinaryq_05_08_notrealbd(iposq) = 
     &                    nbinaryq_05_08_notrealbd(iposq) + 1
                  ELSEIF (xmassmax.LT.1.2) THEN
                     nbinary_08_12_notrealbd(ipos) = 
     &                    nbinary_08_12_notrealbd(ipos) + 1
                     nbinaryq_08_12_notrealbd(iposq) = 
     &                    nbinaryq_08_12_notrealbd(iposq) + 1
                  ELSEIF (xmassmax.GE.1.2) THEN
                     nbinary_12_notrealbd(ipos) = 
     &                    nbinary_12_notrealbd(ipos) + 1
                     nbinaryq_12_notrealbd(iposq) = 
     &                    nbinaryq_12_notrealbd(iposq) + 1
                  ENDIF
               ENDIF
            ENDIF

            IF (pmassn(nkeep1).LT.0.01) THEN
               number_0003 = number_0003 - 1
            ELSEIF (pmassn(nkeep1).GE.1.2) THEN
               number_11 = number_11 - 1
               IF (xmassmin.GT.bdrealmass) number_11notbd = 
     &              number_11notbd - 1
            ELSEIF (pmassn(nkeep1).GE.0.8) THEN
               number_08 = number_08 - 1
               IF (xmassmin.GT.bdrealmass) number_08notbd = 
     &              number_08notbd - 1
            ELSEIF (pmassn(nkeep1).GE.0.5) THEN
               number_05 = number_05 - 1
               IF (xmassmin.GT.bdrealmass) number_05notbd = 
     &              number_05notbd - 1
            ELSEIF (pmassn(nkeep1).GE.0.2) THEN
               number_02 = number_02 - 1
               IF (xmassmin.GT.bdrealmass) number_02notbd = 
     &              number_02notbd - 1
            ELSEIF (pmassn(nkeep1).GE.0.1) THEN
               number_01 = number_01 - 1
               IF (xmassmin.GT.bdrealmass) number_01notbd = 
     &              number_01notbd - 1
            ELSEIF (pmassn(nkeep1).GE.0.07) THEN
               number_007 = number_007 - 1
            ELSEIF (pmassn(nkeep1).GE.0.03) THEN
               number_003 = number_003 - 1
            ELSEIF (pmassn(nkeep1).GE.0.01) THEN
               number_001 = number_001 - 1
            ENDIF

            IF (pmassn(nkeep2).LT.0.01) THEN
               number_0003 = number_0003 - 1
            ELSEIF (pmassn(nkeep2).GE.1.2) THEN
               number_11 = number_11 - 1
               IF (xmassmin.GT.bdrealmass) number_11notbd =
     &              number_11notbd - 1
            ELSEIF (pmassn(nkeep2).GE.0.8) THEN
               number_08 = number_08 - 1
               IF (xmassmin.GT.bdrealmass) number_08notbd =
     &              number_08notbd - 1
            ELSEIF (pmassn(nkeep2).GE.0.5) THEN
               number_05 = number_05 - 1
               IF (xmassmin.GT.bdrealmass) number_05notbd =
     &              number_05notbd - 1
            ELSEIF (pmassn(nkeep2).GE.0.2) THEN
               number_02 = number_02 - 1
               IF (xmassmin.GT.bdrealmass) number_02notbd =
     &              number_02notbd - 1
            ELSEIF (pmassn(nkeep2).GE.0.1) THEN
               number_01 = number_01 - 1
               IF (xmassmin.GT.bdrealmass) number_01notbd =
     &              number_01notbd - 1
            ELSEIF (pmassn(nkeep2).GE.0.07) THEN
               number_007 = number_007 - 1
            ELSEIF (pmassn(nkeep2).GE.0.03) THEN
               number_003 = number_003 - 1
            ELSEIF (pmassn(nkeep2).GE.0.01) THEN
               number_001 = number_001 - 1
            ENDIF

            IF (MAX(pmassn(nkeep1),pmassn(nkeep2)).LT.0.01) THEN
               number2_0003 = number2_0003 + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.1.2) THEN
               number2_11 = number2_11 + 1
               IF (xmassmin.GT.bdrealmass) 
     &              number2_11notbd = number2_11notbd + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.8) THEN
               number2_08 = number2_08 + 1
               IF (xmassmin.GT.bdrealmass) 
     &              number2_08notbd = number2_08notbd + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.5) THEN
               number2_05 = number2_05 + 1

               print *,'Inc number2_05 ',number2_05

               IF (xmassmin.GT.bdrealmass) 
     &              number2_05notbd = number2_05notbd + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.2) THEN
               number2_02 = number2_02 + 1

               print *,'Inc number2_02 ',number2_02

               IF (xmassmin.GT.bdrealmass) 
     &              number2_02notbd = number2_02notbd + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.1) THEN
               number2_01 = number2_01 + 1
               IF (xmassmin.GT.bdrealmass) THEN
                  number2_01notbd = number2_01notbd + 1
                  print *,'number2_01notbd ',number2_01notbd
               ENDIF
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.07) THEN
               number2_007 = number2_007 + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.03) THEN
               number2_003 = number2_003 + 1
            ELSEIF (MAX(pmassn(nkeep1),pmassn(nkeep2)).GE.0.01) THEN
               number2_001 = number2_001 + 1
            ENDIF

         ELSEIF (nbincomp(nkeep1).EQ.3) THEN
            number3 = number3 + 1
            number3unique = number3unique + 1
            number3isolated = number3isolated + 1

            IF (nbincomp1old.EQ.2 .OR. nbincomp(nkeep2).EQ.2) THEN
               number2unique = number2unique - 1
               number2isolated = number2isolated - 1
            ENDIF

            IF (seplastkeep.LT.0.0) THEN
               print *,'ERROR - seplastkeep.LT.0.0 ',seplastkeep,
     &              nkeep1,nkeep2
               STOP
            ENDIF

            ipos = INT((ALOG10(seplastkeep*runit))*5.0+16)
            IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
               nbinary(ipos) = nbinary(ipos) - 1

c               IF (nbinary(ipos).LT.0) THEN
               print *,'Dec Bin T ',nbincomp(nkeep1),nbincomp(nkeep2),
     &              seplastkeep*runit,seplast(nkeep1)*runit,
     &              seplast(nkeep2)*runit
c               ENDIF

               IF (nbincomp1old.EQ.2) THEN
                  xmaxhere = pmassmax(nkeep1)
                  xminhere = pmassmin(nkeep1)

                  xmaxother = pmassmax(nkeep2)

                  iposq = MAX(1,MIN(INT(qratiolastkeep*
     &                 nmassratio+1),nmassratio))
               ELSEIF (nbincomp(nkeep2).EQ.2) THEN
                  xmaxhere = pmassmax(nkeep2)
                  xminhere = pmassmin(nkeep2)

                  xmaxother = pmassmax(nkeep1)

                  iposq = MAX(1,MIN(INT(qratiolast(nkeep2)*
     &                 nmassratio+1),nmassratio))
               ELSE
                  print *,'ERROR xxx ',nbincomp1old,nkeep1,nkeep2
                  STOP
               ENDIF

               IF (xmaxhere.LT.0.03) THEN
                  nbinarybd(ipos) = nbinarybd(ipos) - 1
c                  nbinarybdq(iposq) = nbinarybdq(iposq) - 1
               ELSEIF (xmaxhere.LT.bdmass) THEN
                  nbinarybd(ipos) = nbinarybd(ipos) - 1
c                  nbinarybdq(iposq) = nbinarybdq(iposq) - 1
                  nbinarybdobs(ipos) = nbinarybdobs(ipos) - 1
c                  nbinarybdobsq(iposq) = nbinarybdobsq(iposq) - 1
               ELSEIF (xmaxhere.LT.0.2) THEN
                  nbinary_01_02(ipos) = nbinary_01_02(ipos) - 1

                  print *,'01_02 dec ',nbinary_01_02(ipos),ipos

c                  nbinaryq_01_02(iposq) = nbinaryq_01_02(iposq) - 1
                  IF (xminhere.GT.bdrealmass .AND. 
     &                 xmaxother.GT.bdrealmass) THEN
                     nbinary_01_02_notrealbd(ipos) = 
     &                    nbinary_01_02_notrealbd(ipos) - 1

       print *,'01_02-bd dec ',nbinary_01_02_notrealbd(ipos),ipos

c                     nbinaryq_01_02_notrealbd(iposq) = 
c     &                    nbinaryq_01_02_notrealbd(iposq) - 1
                  ENDIF
               ELSEIF (xmaxhere.LT.0.5) THEN
                  nbinary_02_05(ipos) = nbinary_02_05(ipos) - 1

                  print *,'02_05 dec ',nbinary_02_05(ipos),ipos

                  IF (xminhere.GT.bdrealmass .AND. 
     &                 xmaxother.GT.bdrealmass) THEN
                     nbinary_02_05_notrealbd(ipos) = 
     &                    nbinary_02_05_notrealbd(ipos) - 1

       print *,'02_05-bd dec ',nbinary_02_05_notrealbd(ipos),ipos

c                     nbinaryq_02_05_notrealbd(iposq) = 
c     &                    nbinaryq_02_05_notrealbd(iposq) - 1
                  ENDIF
               ELSEIF (xmaxhere.LT.0.8) THEN
                  nbinary_05_08(ipos) = nbinary_05_08(ipos) - 1
                  IF (xminhere.GT.bdrealmass .AND. 
     &                 xmaxother.GT.bdrealmass) THEN
                     nbinary_05_08_notrealbd(ipos) = 
     &                    nbinary_05_08_notrealbd(ipos) - 1

       print *,'05_08-bd dec ',nbinary_05_08_notrealbd(ipos),ipos

c                     nbinaryq_05_08_notrealbd(iposq) = 
c     &                    nbinaryq_05_08_notrealbd(iposq) - 1
                  ENDIF
               ELSEIF (xmaxhere.LT.1.2) THEN
                  nbinary_08_12(ipos) = nbinary_08_12(ipos) - 1
                  IF (xminhere.GT.bdrealmass .AND. 
     &                 xmaxother.GT.bdrealmass) THEN
                     nbinary_08_12_notrealbd(ipos) = 
     &                    nbinary_08_12_notrealbd(ipos) - 1
c                     nbinaryq_08_12_notrealbd(iposq) = 
c     &                    nbinaryq_08_12_notrealbd(iposq) - 1
                  ENDIF
               ELSEIF (xmaxhere.GE.1.2) THEN
                  nbinary_12(ipos) = nbinary_12(ipos) - 1
                  IF (xminhere.GT.bdrealmass .AND. 
     &                 xmaxother.GT.bdrealmass) THEN
                     nbinary_12_notrealbd(ipos) = 
     &                    nbinary_12_notrealbd(ipos) - 1
c                     nbinaryq_12_notrealbd(iposq) = 
c     &                    nbinaryq_12_notrealbd(iposq) - 1
                  ENDIF
               ENDIF

               ntriple(ipos) = ntriple(ipos) + 1

               print *,'Dec-Mov T ',nbincomp(nkeep1),nbincomp(nkeep2),
     &              seplastkeep*runit,seplast(nkeep1)*runit,
     &              seplast(nkeep2)*runit

               IF (xmassmax.LT.bdmass) THEN
                  ntriplebd(ipos) = ntriplebd(ipos) + 1
               ELSEIF (xmassmax.LT.0.2) THEN
                  ntriple_01_02(ipos) = ntriple_01_02(ipos) + 1
               ELSEIF (xmassmax.LT.0.5) THEN
                  ntriple_02_05(ipos) = ntriple_02_05(ipos) + 1
               ELSEIF (xmassmax.LT.0.8) THEN
                  ntriple_05_08(ipos) = ntriple_05_08(ipos) + 1
               ELSEIF (xmassmax.LT.1.2) THEN
                  ntriple_08_12(ipos) = ntriple_08_12(ipos) + 1
               ELSE
                  ntriple_12(ipos) = ntriple_12(ipos) + 1
               ENDIF
            ENDIF

            ipos = INT((ALOG10(semimajorkeep*runit))*5.0+16)
            IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
               ntriple(ipos) = ntriple(ipos) + 1

               print *,'Dec-Inc T ',nbincomp(nkeep1),nbincomp(nkeep2),
     &              seplastkeep*runit,seplast(nkeep1)*runit,
     &              seplast(nkeep2)*runit,semimajorkeep*runit

               IF (xmassmax.LT.bdmass) THEN
c                  ntriplebd(ipos) = ntriplebd(ipos) + 1
c  This is done below
               ELSEIF (xmassmax.LT.0.2) THEN
                  ntriple_01_02(ipos) = ntriple_01_02(ipos) + 1
               ELSEIF (xmassmax.LT.0.5) THEN
                  ntriple_02_05(ipos) = ntriple_02_05(ipos) + 1
               ELSEIF (xmassmax.LT.0.8) THEN
                  ntriple_05_08(ipos) = ntriple_05_08(ipos) + 1
               ELSEIF (xmassmax.LT.1.2) THEN
                  ntriple_08_12(ipos) = ntriple_08_12(ipos) + 1
               ELSE
                  ntriple_12(ipos) = ntriple_12(ipos) + 1
               ENDIF

            ENDIF
c
c--Remove from counts of individual objects, those that are companions
c
            IF (nbincomp1old.EQ.1) THEN
               IF (pmassn(nkeep1).LT.0.01) THEN
                  number_0003 = number_0003 - 1
               ELSEIF (pmassn(nkeep1).GE.1.2) THEN
                  number_11 = number_11 - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08 = number_08 - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05 = number_05 - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02 = number_02 - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01 = number_01 - 1
               ELSEIF (pmassn(nkeep1).GE.0.07) THEN
                  number_007 = number_007 - 1
               ELSEIF (pmassn(nkeep1).GE.0.03) THEN
                  number_003 = number_003 - 1
               ELSEIF (pmassn(nkeep1).GE.0.01) THEN
                  number_001 = number_001 - 1
               ENDIF

               IF (pmassmin(nkeep2).GE.bdrealmass .AND. 
     &              pmassn(nkeep1).GE.bdrealmass) THEN

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number3_11notbd = number3_11notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number3_08notbd = number3_08notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number3_05notbd = number3_05notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number3_02notbd = number3_02notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number3_01notbd = number3_01notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ENDIF

c
c--Need to increment number of apparent stellar binaries (ignoring BD object)
c
               ELSEIF (pmassmax(nkeep2).GE.bdrealmass .AND.
     &              pmassn(nkeep1).GE.bdrealmass) THEN

                  ipos = INT((ALOG10(semimajorkeep*runit))*5.0+16)
                  IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN

                     print *,'App-Inc T ',nbincomp(nkeep1),
     &                    nbincomp(nkeep2),
     &                    seplastkeep*runit,seplast(nkeep1)*runit,
     &                    seplast(nkeep2)*runit,semimajorkeep*runit

                  ELSE
                     print *,'ERROR AA'
                     STOP
                  ENDIF
                  

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number2_11notbd = number2_11notbd + 1

               nbinary_12_notrealbd(ipos) = 
     &              nbinary_12_notrealbd(ipos) + 1

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number2_08notbd = number2_08notbd + 1

               nbinary_08_12_notrealbd(ipos) = 
     &              nbinary_08_12_notrealbd(ipos) + 1

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number2_05notbd = number2_05notbd + 1

               nbinary_05_08_notrealbd(ipos) = 
     &              nbinary_05_08_notrealbd(ipos) + 1

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number2_02notbd = number2_02notbd + 1

               nbinary_02_05_notrealbd(ipos) = 
     &              nbinary_02_05_notrealbd(ipos) + 1

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number2_01notbd = number2_01notbd + 1
               print *,'number2_01notbd B ',number2_01notbd

               nbinary_01_02_notrealbd(ipos) = 
     &              nbinary_01_02_notrealbd(ipos) + 1

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ENDIF

               ENDIF
            ENDIF

            IF (nbincomp(nkeep2).EQ.1) THEN
               IF (pmassn(nkeep2).LT.0.01) THEN
                  number_0003 = number_0003 - 1
               ELSEIF (pmassn(nkeep2).GE.1.2) THEN
                  number_11 = number_11 - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08 = number_08 - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05 = number_05 - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02 = number_02 - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01 = number_01 - 1
               ELSEIF (pmassn(nkeep2).GE.0.07) THEN
                  number_007 = number_007 - 1
               ELSEIF (pmassn(nkeep2).GE.0.03) THEN
                  number_003 = number_003 - 1
               ELSEIF (pmassn(nkeep2).GE.0.01) THEN
                  number_001 = number_001 - 1
               ENDIF

               IF (pmassmin(nkeep1).GE.bdrealmass .AND.
     &              pmassn(nkeep2).GE.bdrealmass) THEN

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number3_11notbd = number3_11notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number3_08notbd = number3_08notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number3_05notbd = number3_05notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number3_02notbd = number3_02notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number3_01notbd = number3_01notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number2_11notbd = number2_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number2_08notbd = number2_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number2_05notbd = number2_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number2_02notbd = number2_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ENDIF

               ELSEIF (pmassmax(nkeep1).GE.bdrealmass .AND.
     &              pmassn(nkeep2).GE.bdrealmass) THEN
c
c--Need to add in apparent binaries
c
                  ipos = INT((ALOG10(semimajorkeep*runit))*5.0+16)
                  IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN

                     print *,'App-Inc T2 ',nbincomp(nkeep1),
     &                    nbincomp(nkeep2),
     &                    seplastkeep*runit,seplast(nkeep1)*runit,
     &                    seplast(nkeep2)*runit,semimajorkeep*runit

                  ELSE
                     print *,'ERROR AA2'
                     STOP
                  ENDIF

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number2_11notbd = number2_11notbd + 1

               nbinary_12_notrealbd(ipos) = 
     &              nbinary_12_notrealbd(ipos) + 1

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number2_08notbd = number2_08notbd + 1

               nbinary_08_12_notrealbd(ipos) = 
     &              nbinary_08_12_notrealbd(ipos) + 1

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number2_05notbd = number2_05notbd + 1

               nbinary_05_08_notrealbd(ipos) = 
     &              nbinary_05_08_notrealbd(ipos) + 1

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number2_02notbd = number2_02notbd + 1

               nbinary_02_05_notrealbd(ipos) = 
     &              nbinary_02_05_notrealbd(ipos) + 1

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number2_01notbd = number2_01notbd + 1

               print *,'number2_01notbd C ',number2_01notbd

               nbinary_01_02_notrealbd(ipos) = 
     &              nbinary_01_02_notrealbd(ipos) + 1

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ENDIF

               ENDIF
            ENDIF
c
c--Increment count of triples
c
            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).LT.0.001) THEN
               number3_0003 = number3_0003 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number3_11 = number3_11 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number3_08 = number3_08 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number3_05 = number3_05 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number3_02 = number3_02 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number3_01 = number3_01 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.07) THEN
               number3_007 = number3_007 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.03) THEN
               number3_003 = number3_003 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.01) THEN
               number3_001 = number3_001 + 1
            ENDIF
            IF (nbincomp(nkeep2).EQ.2) THEN
               nkeephere = nkeep2
            ELSE
               nkeephere = nkeep1
            ENDIF
c
c--Remove the binaries that are components of triples
c
            IF (pmassmax(nkeephere).LT.0.01) THEN
               number2_0003 = number2_0003 - 1
            ELSEIF (pmassmax(nkeephere).GE.1.2) THEN
               number2_11 = number2_11 - 1
            ELSEIF (pmassmax(nkeephere).GE.0.8) THEN
               number2_08 = number2_08 - 1
            ELSEIF (pmassmax(nkeephere).GE.0.5) THEN
               number2_05 = number2_05 - 1

               print *,'Dec number2_05 ',number2_05

            ELSEIF (pmassmax(nkeephere).GE.0.2) THEN
               number2_02 = number2_02 - 1

               print *,'Dec number2_02 ',number2_02

            ELSEIF (pmassmax(nkeephere).GE.0.1) THEN
               number2_01 = number2_01 - 1
            ELSEIF (pmassmax(nkeephere).GE.0.07) THEN
               number2_007 = number2_007 - 1
            ELSEIF (pmassmax(nkeephere).GE.0.03) THEN
               number2_003 = number2_003 - 1
            ELSEIF (pmassmax(nkeephere).GE.0.01) THEN
               number2_001 = number2_001 - 1
            ENDIF
c
c--Remove the VLM binaries that are components of triples
c
            IF (pmassmax(nkeephere).LT.bdmass) THEN
               number2bd = number2bd - 1
               IF (pmassmax(nkeephere).GT.0.03) THEN
                  number2bdobs = number2bdobs - 1
               ENDIF
            ENDIF
c
c--Remove the stellar binaries that are components of triples
c
            IF (pmassmin(nkeephere).GE.bdmass) THEN
               number2star = number2star - 1
            ELSEIF (pmassmax(nkeephere).GE.bdmass) THEN
               numberstarbd = numberstarbd - 1
            ENDIF

            print *,interest(nkeep1),interest(nkeep1)/1000,
     &           interest(nkeep1)/1000,INT(interest(nkeep1)-
     &              1000*(interest(nkeep1)/1000))/100 +
     &              interest(nkeep1)/1000

            IF (interest(nkeep1).LT.100) THEN
c
c--All components of the triple are stars
c
               numberstar3 = numberstar3 + 1
            ELSEIF (INT(interest(nkeep1)-
     &              1000*(interest(nkeep1)/1000))/100 + 
     &              interest(nkeep1)/1000 .EQ. 3) THEN
c
c--All components of the triple are VLM (brown dwarfs)
c
               numberbd3 = numberbd3 + 1
c
c--Observable if any of the BDs are observable
c
               IF (INT(interest(nkeep1)-
     &              1000*(interest(nkeep1)/1000))/100.GE.1) THEN
                  numberbdobs3 = numberbdobs3 + 1
               ENDIF
               IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
                  ntriplebd(ipos) = ntriplebd(ipos) + 1
               ENDIF
            ELSE
               numbermixed3 = numbermixed3 + 1
            ENDIF

            interest(nkeep1) = interest(nkeep1) + 10000
            WRITE (*,*) 'Triple'

            IF (nbincomp1old.EQ.1) THEN
               IF (pmassn(nkeep1).LT.bdmass) THEN
                  numberbd345 = numberbd345 + 1
                  numberbdsingle = numberbdsingle - 1
                  IF (pmassn(nkeep1).GT.0.03) THEN
                     numberbdobs345 = numberbdobs345 + 1
                     numberbdobssingle = numberbdobssingle - 1
                  ENDIF
               ELSE
                  numberstar345 = numberstar345 + 1
                  numberstarsingle = numberstarsingle - 1
               ENDIF
            ELSEIF (nbincomp(nkeep2).EQ.1) THEN
               IF (pmassn(nkeep2).LT.bdmass) THEN
                  numberbd345 = numberbd345 + 1
                  numberbdsingle = numberbdsingle - 1
                  IF (pmassn(nkeep2).GT.0.03) THEN
                     numberbdobs345 = numberbdobs345 + 1
                     numberbdobssingle = numberbdobssingle - 1
                  ENDIF
               ELSE
                  numberstar345 = numberstar345 + 1
                  numberstarsingle = numberstarsingle - 1
               ENDIF
            ENDIF
c
c--Quadruples
c
         ELSEIF (nbincomp(nkeep1).EQ.4) THEN
            number4 = number4 + 1
            number4unique = number4unique + 1
            number4isolated = number4isolated + 1

            ipos = INT((ALOG10(semimajorkeep*runit))*5.0+16)
            IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
               nquad(ipos) = nquad(ipos) + 1
            ENDIF

            IF (nbincomp1old.EQ.2) THEN
               number2unique = number2unique - 1
               number2isolated = number2isolated - 1

               IF (seplastkeep.LT.0.0) THEN
                  print *,'ERROR - seplastkeep.LT.0.0 -Q ',seplastkeep,
     &                 nkeep1,nkeep2
                  STOP
               ENDIF

               ipos = INT((ALOG10(seplastkeep*runit))*5.0+16)
               IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
                  nbinary(ipos) = nbinary(ipos) - 1
                  nquad(ipos) = nquad(ipos) + 1

c                  IF (nbinary(ipos).LT.0) THEN
                  print *,'Dec Q ',nbincomp1old,
     &                 nbincomp(nkeep1),nbincomp(nkeep2),
     &                 seplastkeep*runit,seplast(nkeep1)*runit,
     &                 seplast(nkeep2)*runit,pmassmax(nkeep1),
     &                 pmassmin(nkeep1)
c                  ENDIF

                  xmaxhere = pmassmax(nkeep1)
                  xminhere = pmassmin(nkeep1)

                  iposq = MAX(1,MIN(INT(qratiolastkeep*
     &                 nmassratio+1),nmassratio))

                  IF (xmaxhere.LT.0.03) THEN
                     nbinarybd(ipos) = nbinarybd(ipos) - 1
c                     nbinarybdq(iposq) = nbinarybdq(iposq) - 1
                  ELSEIF (xmaxhere.LT.bdmass) THEN
                     nbinarybd(ipos) = nbinarybd(ipos) - 1
c                     nbinarybdq(iposq) = nbinarybdq(iposq) - 1
                     nbinarybdobs(ipos) = nbinarybdobs(ipos) - 1
c                     nbinarybdobsq(iposq) = nbinarybdobsq(iposq) - 1
                  ELSEIF (xmaxhere.LT.0.2) THEN
                     nbinary_01_02(ipos) = nbinary_01_02(ipos) - 1

                  print *,'01_02 dec Q1 ',nbinary_01_02(ipos),ipos

c                     nbinaryq_01_02(iposq) = nbinaryq_01_02(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_01_02_notrealbd(ipos) =
     &                       nbinary_01_02_notrealbd(ipos) - 1

       print *,'01_02-bd dec Q1 ',nbinary_01_02_notrealbd(ipos),ipos

c                        nbinaryq_01_02_notrealbd(iposq) =
c     &                       nbinaryq_01_02_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.LT.0.5) THEN
                     nbinary_02_05(ipos) = nbinary_02_05(ipos) - 1

                  print *,'02_05 dec Q1 ',nbinary_02_05(ipos),ipos

c                     nbinaryq_02_05(iposq) = nbinaryq_02_05(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_02_05_notrealbd(ipos) =
     &                       nbinary_02_05_notrealbd(ipos) - 1

       print *,'02_05-bd dec Q1 ',nbinary_02_05_notrealbd(ipos),ipos

c                        nbinaryq_02_05_notrealbd(iposq) =
c     &                       nbinaryq_02_05_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.LT.0.8) THEN
                     nbinary_05_08(ipos) = nbinary_05_08(ipos) - 1
c                     nbinaryq_05_08(iposq) = nbinaryq_05_08(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_05_08_notrealbd(ipos) =
     &                       nbinary_05_08_notrealbd(ipos) - 1

       print *,'05_08-bd dec Q1 ',nbinary_05_08_notrealbd(ipos),ipos

c                        nbinaryq_05_08_notrealbd(iposq) =
c     &                       nbinaryq_05_08_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.LT.1.2) THEN
                     nbinary_08_12(ipos) = nbinary_08_12(ipos) - 1
c                     nbinaryq_08_12(iposq) = nbinaryq_08_12(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_08_12_notrealbd(ipos) =
     &                       nbinary_08_12_notrealbd(ipos) - 1
c                        nbinaryq_08_12_notrealbd(iposq) =
c     &                       nbinaryq_08_12_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.GE.1.2) THEN
                     nbinary_12(ipos) = nbinary_12(ipos) - 1
c                     nbinaryq_12(iposq) = nbinaryq_12(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_12_notrealbd(ipos) =
     &                       nbinary_12_notrealbd(ipos) - 1
c                        nbinaryq_12_notrealbd(iposq) =
c     &                       nbinaryq_12_notrealbd(iposq) - 1
                     ENDIF
                  ENDIF

               ENDIF

            ENDIF
            IF (nbincomp(nkeep2).EQ.2) THEN
               number2unique = number2unique - 1
               number2isolated = number2isolated - 1

               IF (seplast(nkeep2).LT.0.0) THEN
                  print *,'ERROR - seplastkeep.LT.0.0-Q2 ',
     &                 seplast(nkeep2),nkeep1,nkeep2
                  STOP
               ENDIF

               ipos = INT((ALOG10(seplast(nkeep2)*runit))*5.0+16)
               IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
                  nbinary(ipos) = nbinary(ipos) - 1
                  nquad(ipos) = nquad(ipos) + 1

                  print *,'Dec Q2 ',nbincomp1old,
     &                 nbincomp(nkeep1),nbincomp(nkeep2),
     &                 seplastkeep*runit,seplast(nkeep1)*runit,
     &                 seplast(nkeep2)*runit

               IF (nbinary(ipos).LT.0) THEN
                  print *,'Change42 ',nbincomp1old,
     &                 nbincomp(nkeep1),nbincomp(nkeep2),
     &                 seplastkeep,seplast(nkeep1),seplast(nkeep2)
               ENDIF

                  xmaxhere = pmassmax(nkeep2)
                  xminhere = pmassmin(nkeep2)

                  iposq = MAX(1,MIN(INT(qratiolast(nkeep2)*
     &                 nmassratio+1),nmassratio))

                  IF (xmaxhere.LT.0.03) THEN
                     nbinarybd(ipos) = nbinarybd(ipos) - 1
c                     nbinarybdq(iposq) = nbinarybdq(iposq) - 1
                  ELSEIF (xmaxhere.LT.bdmass) THEN
                     nbinarybd(ipos) = nbinarybd(ipos) - 1
c                     nbinarybdq(iposq) = nbinarybdq(iposq) - 1
                     nbinarybdobs(ipos) = nbinarybdobs(ipos) - 1
c                     nbinarybdobsq(iposq) = nbinarybdobsq(iposq) - 1
                  ELSEIF (xmaxhere.LT.0.2) THEN
                     nbinary_01_02(ipos) = nbinary_01_02(ipos) - 1

                  print *,'01_02 dec Q2 ',nbinary_01_02(ipos),ipos

c                     nbinaryq_01_02(iposq) = nbinaryq_01_02(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_01_02_notrealbd(ipos) =
     &                       nbinary_01_02_notrealbd(ipos) - 1

       print *,'01_02-bg dec Q2 ',nbinary_01_02_notrealbd(ipos),ipos

c                        nbinaryq_01_02_notrealbd(iposq) =
c     &                       nbinaryq_01_02_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.LT.0.5) THEN
                     nbinary_02_05(ipos) = nbinary_02_05(ipos) - 1

                  print *,'02_05 dec Q2 ',nbinary_02_05(ipos),ipos

c                     nbinaryq_02_05(iposq) = nbinaryq_02_05(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_02_05_notrealbd(ipos) =
     &                       nbinary_02_05_notrealbd(ipos) - 1

       print *,'02_05-bg dec Q2 ',nbinary_02_05_notrealbd(ipos),ipos

c                        nbinaryq_02_05_notrealbd(iposq) =
c     &                       nbinaryq_02_05_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.LT.0.8) THEN
                     nbinary_05_08(ipos) = nbinary_05_08(ipos) - 1
c                     nbinaryq_05_08(iposq) = nbinaryq_05_08(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
c                        nbinary_05_08_notrealbd(ipos) =
c     &                       nbinary_05_08_notrealbd(ipos) - 1

       print *,'05_08-bg dec Q2 ',nbinary_05_08_notrealbd(ipos),ipos

c                        nbinaryq_05_08_notrealbd(iposq) =
c     &                       nbinaryq_05_08_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.LT.1.2) THEN
                     nbinary_08_12(ipos) = nbinary_08_12(ipos) - 1
c                     nbinaryq_08_12(iposq) = nbinaryq_08_12(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_08_12_notrealbd(ipos) =
     &                       nbinary_08_12_notrealbd(ipos) - 1
c                        nbinaryq_08_12_notrealbd(iposq) =
c     &                       nbinaryq_08_12_notrealbd(iposq) - 1
                     ENDIF
                  ELSEIF (xmaxhere.GE.1.2) THEN
                     nbinary_12(ipos) = nbinary_12(ipos) - 1
c                     nbinaryq_12(iposq) = nbinaryq_12(iposq) - 1
                     IF (xminhere.GT.bdrealmass) THEN
                        nbinary_12_notrealbd(ipos) =
     &                       nbinary_12_notrealbd(ipos) - 1
c                        nbinaryq_12_notrealbd(iposq) =
c     &                       nbinaryq_12_notrealbd(iposq) - 1
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
            IF (nbincomp1old.EQ.3 .OR. nbincomp(nkeep2).EQ.3) THEN
               number3unique = number3unique - 1
               number3isolated = number3isolated - 1

               IF (nbincomp1old.EQ.3) THEN
                  seplasthere = seplastkeep
                  seplast2here = seplast2keep
               ELSE
                  seplasthere = seplast(nkeep2)
                  seplast2here = seplast2(nkeep2)
               ENDIF

               IF (seplasthere.LT.0.0) THEN
                  print *,'ERROR - seplasthere.LT.0.0-Q3 ',seplasthere,
     &                 nkeep1,nkeep2
                  STOP
               ENDIF

               IF (seplast2here.LT.0.0) THEN
                  print *,'ERROR - seplast2here.LT.0.0-Q ',seplast2here,
     &                 nkeep1,nkeep2
                  STOP
               ENDIF

               ipos = INT((ALOG10(seplasthere*runit))*5.0+16)
               IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
                  ntriple(ipos) = ntriple(ipos) - 1
                  nquad(ipos) = nquad(ipos) + 1

                  print *,'Dec QT ',nbincomp1old,
     &                 nbincomp(nkeep1),nbincomp(nkeep2),
     &                 seplasthere*runit,seplast(nkeep1)*runit,
     &                 seplast(nkeep2)*runit

                  IF (nbincomp1old.EQ.3) THEN
                     IF (pmassmax(nkeep1).LT.bdmass) THEN
                        ntriplebd(ipos) = ntriplebd(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.0.2) THEN
                        ntriple_01_02(ipos) = ntriple_01_02(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.0.5) THEN
                        ntriple_02_05(ipos) = ntriple_02_05(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.0.8) THEN
                        ntriple_05_08(ipos) = ntriple_05_08(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.1.2) THEN
                        ntriple_08_12(ipos) = ntriple_08_12(ipos) - 1
                     ELSE
                        ntriple_12(ipos) = ntriple_12(ipos) - 1
                     ENDIF

                  ELSE
                     IF (pmassmax(nkeep2).LT.bdmass) THEN
                        ntriplebd(ipos) = ntriplebd(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.0.2) THEN
                        ntriple_01_02(ipos) = ntriple_01_02(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.0.5) THEN
                        ntriple_02_05(ipos) = ntriple_02_05(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.0.8) THEN
                        ntriple_05_08(ipos) = ntriple_05_08(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.1.2) THEN
                        ntriple_08_12(ipos) = ntriple_08_12(ipos) - 1
                     ELSE
                        ntriple_12(ipos) = ntriple_12(ipos) - 1
                     ENDIF
                  ENDIF
               ENDIF

               ipos = INT((ALOG10(seplast2here*runit))*5.0+16)
               IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
                  ntriple(ipos) = ntriple(ipos) - 1
                  nquad(ipos) = nquad(ipos) + 1

                  print *,'Dec QT2 ',nbincomp1old,
     &                 nbincomp(nkeep1),nbincomp(nkeep2),
     &                 seplast2here*runit,seplast(nkeep1)*runit,
     &                 seplast(nkeep2)*runit

                  IF (nbincomp1old.EQ.3) THEN
                     IF (pmassmax(nkeep1).LT.bdmass) THEN
                        ntriplebd(ipos) = ntriplebd(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.0.2) THEN
                        ntriple_01_02(ipos) = ntriple_01_02(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.0.5) THEN
                        ntriple_02_05(ipos) = ntriple_02_05(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.0.8) THEN
                        ntriple_05_08(ipos) = ntriple_05_08(ipos) - 1
                     ELSEIF (pmassmax(nkeep1).LT.1.2) THEN
                        ntriple_08_12(ipos) = ntriple_08_12(ipos) - 1
                     ELSE
                        ntriple_12(ipos) = ntriple_12(ipos) - 1
                     ENDIF

                  ELSE
                     IF (pmassmax(nkeep2).LT.bdmass) THEN
                        ntriplebd(ipos) = ntriplebd(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.0.2) THEN
                        ntriple_01_02(ipos) = ntriple_01_02(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.0.5) THEN
                        ntriple_02_05(ipos) = ntriple_02_05(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.0.8) THEN
                        ntriple_05_08(ipos) = ntriple_05_08(ipos) - 1
                     ELSEIF (pmassmax(nkeep2).LT.1.2) THEN
                        ntriple_08_12(ipos) = ntriple_08_12(ipos) - 1
                     ELSE
                        ntriple_12(ipos) = ntriple_12(ipos) - 1
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
c
c--Remove from counts of individual objects, those that are companions
c
            IF (nbincomp1old.EQ.1) THEN
               IF (pmassn(nkeep1).LT.0.01) THEN
                  number_0003 = number_0003 - 1
               ELSEIF (pmassn(nkeep1).GE.1.2) THEN
                  number_11 = number_11 - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08 = number_08 - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05 = number_05 - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02 = number_02 - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01 = number_01 - 1
               ELSEIF (pmassn(nkeep1).GE.0.07) THEN
                  number_007 = number_007 - 1
               ELSEIF (pmassn(nkeep1).GE.0.03) THEN
                  number_003 = number_003 - 1
               ELSEIF (pmassn(nkeep1).GE.0.01) THEN
                  number_001 = number_001 - 1
               ENDIF

               IF (pmassmin(nkeep2).GE.bdrealmass .AND.
     &              pmassn(nkeep1).GE.bdrealmass) THEN

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number4_11notbd = number4_11notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number4_08notbd = number4_08notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number4_05notbd = number4_05notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number4_02notbd = number4_02notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number4_01notbd = number4_01notbd + 1

               IF (pmassmax(nkeep2).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep2).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep1).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep1).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ENDIF

               ELSEIF (pmassmax(nkeep2).GE.bdrealmass .AND.
     &              pmassn(nkeep1).GE.bdrealmass) THEN
c
c--Can't be sure whether the triple is 1 star and 2 BDs or 2 stars and 1 BD.
c

               ENDIF
            ENDIF

            IF (nbincomp(nkeep2).EQ.1) THEN
               IF (pmassn(nkeep2).LT.0.01) THEN
                  number_0003 = number_0003 - 1
               ELSEIF (pmassn(nkeep2).GE.1.2) THEN
                  number_11 = number_11 - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08 = number_08 - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05 = number_05 - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02 = number_02 - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01 = number_01 - 1
               ELSEIF (pmassn(nkeep2).GE.0.07) THEN
                  number_007 = number_007 - 1
               ELSEIF (pmassn(nkeep2).GE.0.03) THEN
                  number_003 = number_003 - 1
               ELSEIF (pmassn(nkeep2).GE.0.01) THEN
                  number_001 = number_001 - 1
               ENDIF

               IF (pmassmin(nkeep1).GE.bdrealmass .AND.
     &              pmassn(nkeep2).GE.bdrealmass) THEN

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number4_11notbd = number4_11notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number4_08notbd = number4_08notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number4_05notbd = number4_05notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number4_02notbd = number4_02notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number4_01notbd = number4_01notbd + 1

               IF (pmassmax(nkeep1).GE.1.2) THEN
                  number3_11notbd = number3_11notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.8) THEN
                  number3_08notbd = number3_08notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.5) THEN
                  number3_05notbd = number3_05notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.2) THEN
                  number3_02notbd = number3_02notbd - 1
               ELSEIF (pmassmax(nkeep1).GE.0.1) THEN
                  number3_01notbd = number3_01notbd - 1
               ENDIF

               IF (pmassn(nkeep2).GE.1.2) THEN
                  number_11notbd = number_11notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.8) THEN
                  number_08notbd = number_08notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.5) THEN
                  number_05notbd = number_05notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.2) THEN
                  number_02notbd = number_02notbd - 1
               ELSEIF (pmassn(nkeep2).GE.0.1) THEN
                  number_01notbd = number_01notbd - 1
               ENDIF
            ENDIF

               ELSEIF (pmassmax(nkeep1).GE.bdrealmass .AND.
     &              pmassn(nkeep2).GE.bdrealmass) THEN
c
c--Can't be sure whether the triple is 1 star and 2 BDs or 2 stars and 1 BD.
c

               ENDIF
            ENDIF
c
c--Increment count of quadruples
c
            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).LT.0.01) THEN
               number4_0003 = number4_0003 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number4_11 = number4_11 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number4_08 = number4_08 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number4_05 = number4_05 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number4_02 = number4_02 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number4_01 = number4_01 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.07) THEN
               number4_007 = number4_007 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.03) THEN
               number4_003 = number4_003 + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.01) THEN
               number4_001 = number4_001 + 1
            ENDIF

            IF (nbincomp(nkeep2).EQ.2) THEN
c
c--Quad is composed of two binaries
c
               nkeephere = nkeep1
               IF (pmassmax(nkeephere).LT.0.01) THEN
                  number2_0003 = number2_0003 - 1
               ELSEIF (pmassmax(nkeephere).GE.1.2) THEN
                  number2_11 = number2_11 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.8) THEN
                  number2_08 = number2_08 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.5) THEN
                  number2_05 = number2_05 - 1

               print *,'Dec number2_05 ',number2_05

               ELSEIF (pmassmax(nkeephere).GE.0.2) THEN
                  number2_02 = number2_02 - 1

               print *,'Dec number2_02 ',number2_02

               ELSEIF (pmassmax(nkeephere).GE.0.1) THEN
                  number2_01 = number2_01 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.07) THEN
                  number2_007 = number2_007 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.03) THEN
                  number2_003 = number2_003 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.01) THEN
                  number2_001 = number2_001 - 1
               ENDIF
c
c--Remove the VLM triples that are components of quads
c
               IF (pmassmax(nkeephere).LT.bdmass) THEN
                  number2bd = number2bd - 1
                  IF (pmassmax(nkeephere).GT.0.03) THEN
                     number2bdobs = number2bdobs - 1
                  ENDIF
               ENDIF
c
c--Remove the stellar binaries that are components of quads
c
               IF (pmassmin(nkeephere).GE.bdmass) THEN
                  number2star = number2star - 1
               ELSEIF (pmassmax(nkeephere).GE.bdmass) THEN
                  numberstarbd = numberstarbd - 1
               ENDIF
c
c--Count for statistics without real BDs
c
               IF (pmassmax(nkeep1).GE.bdrealmass .AND.
     &              pmassmax(nkeep2).GE.bdrealmass) THEN
                  IF (pmassmin(nkeephere).GE.bdrealmass) THEN
c
c----Single is star, then new object is
c
                     IF (pmassmin(nkeephere).GE.1.2) THEN
                        number2_11notbd = number2_11notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.8) THEN
                        number2_08notbd = number2_08notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.5) THEN
                        number2_05notbd = number2_05notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.2) THEN
                        number2_02notbd = number2_02notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.1) THEN
                        number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
                     ENDIF
                  ELSE
                     IF (pmassmin(nkeephere).GE.1.2) THEN
                        number_11notbd = number_11notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.8) THEN
                        number_08notbd = number_08notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.5) THEN
                        number_05notbd = number_05notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.2) THEN
                        number_02notbd = number_02notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.1) THEN
                        number_01notbd = number_01notbd - 1
                     ENDIF
                  ENDIF
               ENDIF
c
c--Do other binary
c
               nkeephere = nkeep2
               IF (pmassmax(nkeephere).LT.0.01) THEN
                  number2_0003 = number2_0003 - 1
               ELSEIF (pmassmax(nkeephere).GE.1.2) THEN
                  number2_11 = number2_11 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.8) THEN
                  number2_08 = number2_08 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.5) THEN
                  number2_05 = number2_05 - 1

               print *,'Dec number2_05 ',number2_05

               ELSEIF (pmassmax(nkeephere).GE.0.2) THEN
                  number2_02 = number2_02 - 1

               print *,'Dec number2_02 ',number2_02

               ELSEIF (pmassmax(nkeephere).GE.0.1) THEN
                  number2_01 = number2_01 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.07) THEN
                  number2_007 = number2_007 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.03) THEN
                  number2_003 = number2_003 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.01) THEN
                  number2_001 = number2_001 - 1
               ENDIF
c
c--Remove the VLM triples that are components of quads
c
               IF (pmassmax(nkeephere).LT.bdmass) THEN
                  number2bd = number2bd - 1
                  IF (pmassmax(nkeephere).GT.0.03) THEN
                     number2bdobs = number2bdobs - 1
                  ENDIF
               ENDIF
c
c--Remove the stellar binaries that are components of quads
c
               IF (pmassmin(nkeephere).GE.bdmass) THEN
                  number2star = number2star - 1
               ELSEIF (pmassmax(nkeephere).GE.bdmass) THEN
                  numberstarbd = numberstarbd - 1
               ENDIF
c
c--Count for statistics without real BDs
c
               IF (pmassmax(nkeep1).GE.bdrealmass .AND.
     &              pmassmax(nkeep2).GE.bdrealmass) THEN
                  IF (pmassmin(nkeephere).GE.bdrealmass) THEN
                     IF (pmassmin(nkeephere).GE.1.2) THEN
                        number2_11notbd = number2_11notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.8) THEN
                        number2_08notbd = number2_08notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.5) THEN
                        number2_05notbd = number2_05notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.2) THEN
                        number2_02notbd = number2_02notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.1) THEN
                        number2_01notbd = number2_01notbd - 1
                  print *,'number2_01notbd Dec ',number2_01notbd
                     ENDIF
                  ELSE
                     IF (pmassmin(nkeephere).GE.1.2) THEN
                        number_11notbd = number_11notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.8) THEN
                        number_08notbd = number_08notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.5) THEN
                        number_05notbd = number_05notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.2) THEN
                        number_02notbd = number_02notbd - 1
                     ELSEIF (pmassmin(nkeephere).GE.0.1) THEN
                        number_01notbd = number_01notbd - 1
                     ENDIF
                  ENDIF
               ENDIF
c
c--For statistics without BDs need to decide whether this new system is
c     unchanged, a new binary, a new triple or a new quadruple
c
               IF (pmassmax(nkeep1).GE.bdrealmass .AND.
     &              pmassmax(nkeep2).GE.bdrealmass) THEN

               IF (pmassmin(nkeep1).GE.bdrealmass .AND.
     &              pmassmin(nkeep2).GE.bdrealmass) THEN

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number4_11notbd = number4_11notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number4_08notbd = number4_08notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number4_05notbd = number4_05notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number4_02notbd = number4_02notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number4_01notbd = number4_01notbd + 1
            ENDIF

               ELSEIF (pmassmin(nkeep1).GE.bdrealmass .AND.
     &              pmassmin(nkeep2).LT.bdrealmass .OR.
     &              pmassmin(nkeep2).GE.bdrealmass .AND.
     &              pmassmin(nkeep1).LT.bdrealmass) THEN

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number3_11notbd = number3_11notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number3_08notbd = number3_08notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number3_05notbd = number3_05notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number3_02notbd = number3_02notbd + 1
            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number3_01notbd = number3_01notbd + 1
            ENDIF

               ELSEIF (pmassmin(nkeep1).LT.bdrealmass .AND.
     &              pmassmin(nkeep2).LT.bdrealmass) THEN
c
c--Need to add in apparent binaries
c
                  ipos = INT((ALOG10(semimajorkeep*runit))*5.0+16)
                  IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN

                     print *,'App-Inc Q ',nbincomp(nkeep1),
     &                    nbincomp(nkeep2),
     &                    seplastkeep*runit,seplast(nkeep1)*runit,
     &                    seplast(nkeep2)*runit,semimajorkeep*runit

                  ELSE
                     print *,'ERROR AA3'
                     STOP
                  ENDIF

            IF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.1.2) THEN
               number2_11notbd = number2_11notbd + 1

               nbinary_12_notrealbd(ipos) = 
     &              nbinary_12_notrealbd(ipos) + 1

            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.8) THEN
               number2_08notbd = number2_08notbd + 1

               nbinary_08_12_notrealbd(ipos) = 
     &              nbinary_08_12_notrealbd(ipos) + 1

            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.5) THEN
               number2_05notbd = number2_05notbd + 1

               nbinary_05_08_notrealbd(ipos) = 
     &              nbinary_05_08_notrealbd(ipos) + 1

            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.2) THEN
               number2_02notbd = number2_02notbd + 1

               nbinary_02_05_notrealbd(ipos) = 
     &              nbinary_02_05_notrealbd(ipos) + 1

            ELSEIF (MAX(pmassmax(nkeep1),pmassmax(nkeep2)).GE.0.1) THEN
               number2_01notbd = number2_01notbd + 1
               
               print *,'number2_01notbd D ',number2_01notbd

               nbinary_01_02_notrealbd(ipos) = 
     &              nbinary_01_02_notrealbd(ipos) + 1

            ENDIF
                  
               ENDIF
               ENDIF
            ELSE
c
c--Quad is composed of a triple and a single
c
               IF (nbincomp(nkeep2).EQ.3) THEN
                  nkeephere = nkeep2
               ELSE
                  nkeephere = nkeep1
               ENDIF

               IF (pmassmax(nkeephere).LT.0.01) THEN
                  number3_0003 = number3_0003 - 1
               ELSEIF (pmassmax(nkeephere).GE.1.2) THEN
                  number3_11 = number3_11 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.8) THEN
                  number3_08 = number3_08 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.5) THEN
                  number3_05 = number3_05 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.2) THEN
                  number3_02 = number3_02 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.1) THEN
                  number3_01 = number3_01 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.07) THEN
                  number3_007 = number3_007 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.03) THEN
                  number3_003 = number3_003 - 1
               ELSEIF (pmassmax(nkeephere).GE.0.01) THEN
                  number3_001 = number3_001 - 1
               ENDIF
c
c--Remove the VLM triples that are components of quads
c
               IF (pmassmax(nkeephere).LT.bdmass) THEN
                  numberbd3 = numberbd3 - 1
                  IF (pmassmax(nkeephere).GT.0.03) THEN
                     numberbdobs3 = numberbdobs3 - 1
                  ENDIF
               ENDIF
c
c--Remove the stellar triples that are components of quads
c
               IF (pmassmax(nkeephere).GE.bdmass) THEN
                  IF (pmassmin(nkeephere).GE.bdmass) THEN
                     numberstar3 = numberstar3 - 1
                  ELSE
                     numbermixed3 = numbermixed3 - 1
                  ENDIF
               ENDIF
            ENDIF

            interesttest = interest(nkeep1) - 
     &           (interest(nkeep1)/10000)*10000
            IF (interesttest.LT.100) THEN
c
c--All components of the quadruple are stars
c
               numberstar4 = numberstar4 + 1
            ELSEIF (INT(interesttest-
     &              1000*(interesttest/1000))/100 + 
     &              interesttest/1000 .EQ. 4) THEN
c
c--All components of the quadruple are VLM (brown dwarfs)
c
               numberbd4 = numberbd4 + 1
               IF (INT(interesttest-
     &              1000*(interesttest/1000))/100.GE.1) THEN
                  numberbdobs4 = numberbdobs4 + 1
               ENDIF
            ELSE
               numbermixed4 = numbermixed4 + 1
            ENDIF

            interest(nkeep1) = interest(nkeep1) + 100000
            WRITE (*,*) 'Quadruple'
            IF (nbincomp1old.EQ.1) THEN
               IF (pmassn(nkeep1).LT.bdmass) THEN
                  numberbd345 = numberbd345 + 1
                  numberbdsingle = numberbdsingle - 1
                  IF (pmassn(nkeep1).GT.0.03) THEN
                     numberbdobs345 = numberbdobs345 + 1
                     numberbdobssingle = numberbdobssingle - 1
                  ENDIF
               ELSE
                  numberstar345 = numberstar345 + 1
                  numberstarsingle = numberstarsingle - 1
               ENDIF
            ELSEIF (nbincomp(nkeep2).EQ.1) THEN
               IF (pmassn(nkeep2).LT.bdmass) THEN
                  numberbd345 = numberbd345 + 1
                  numberbdsingle = numberbdsingle - 1
                  IF (pmassn(nkeep2).GT.0.03) THEN
                     numberbdobs345 = numberbdobs345 + 1
                     numberbdobssingle = numberbdobssingle - 1
                  ENDIF
               ELSE
                  numberstar345 = numberstar345 + 1
                  numberstarsingle = numberstarsingle - 1
               ENDIF
            ENDIF
c
c--Quintuples
c
         ELSEIF (nbincomp(nkeep1).EQ.5) THEN
            number5 = number5 + 1
            number5unique = number5unique + 1
            number5isolated = number5isolated + 1

            IF (nbincomp1old.EQ.2) THEN
               number2unique = number2unique - 1
               number2isolated = number2isolated - 1
            ENDIF
            IF (nbincomp(nkeep2).EQ.2) THEN
               number2unique = number2unique - 1
               number2isolated = number2isolated - 1
            ENDIF
            IF (nbincomp1old.EQ.3) THEN
               number3unique = number3unique - 1
               number3isolated = number3isolated - 1
            ENDIF
            IF (nbincomp(nkeep2).EQ.3) THEN
               number3unique = number3unique - 1
               number3isolated = number3isolated - 1
            ENDIF
            IF (nbincomp1old.EQ.4 .OR. nbincomp(nkeep2).EQ.4) THEN
               number4unique = number4unique - 1
               number4isolated = number4isolated - 1
            ENDIF

            interest(nkeep1) = interest(nkeep1) + 1000000            
            WRITE (*,*) 'Quintuple'
            IF (nbincomp1old.EQ.1) THEN
               IF (pmassn(nkeep1).LT.bdmass) THEN
                  numberbd345 = numberbd345 + 1
                  numberbdsingle = numberbdsingle - 1
                  IF (pmassn(nkeep1).GT.0.03) THEN
                     numberbdobs345 = numberbdobs345 + 1
                     numberbdobssingle = numberbdobssingle - 1
                  ENDIF
               ELSE
                  numberstar345 = numberstar345 + 1
                  numberstarsingle = numberstarsingle - 1
               ENDIF
            ELSEIF (nbincomp(nkeep2).EQ.1) THEN
               IF (pmassn(nkeep2).LT.bdmass) THEN
                  numberbd345 = numberbd345 + 1
                  numberbdsingle = numberbdsingle - 1
                  IF (pmassn(nkeep2).GT.0.03) THEN
                     numberbdobs345 = numberbdobs345 + 1
                     numberbdobssingle = numberbdobssingle - 1
                  ENDIF
               ELSE
                  numberstar345 = numberstar345 + 1
                  numberstarsingle = numberstarsingle - 1
               ENDIF
            ENDIF
         ELSEIF (nbincomp(nkeep1).GE.6) THEN
            IF (nbincomp1old.EQ.2) THEN
               number2isolated = number2isolated - 1
            ENDIF
            IF (nbincomp(nkeep2).EQ.2) THEN
               number2isolated = number2isolated - 1
            ENDIF
            IF (nbincomp1old.EQ.3) THEN
               number3isolated = number3isolated - 1
            ENDIF
            IF (nbincomp(nkeep2).EQ.3) THEN
               number3isolated = number3isolated - 1
            ENDIF
            IF (nbincomp1old.EQ.4) THEN
               number4isolated = number4isolated - 1
            ENDIF
            IF (nbincomp(nkeep2).EQ.4) THEN
               number4isolated = number4isolated - 1
            ENDIF
            IF (nbincomp1old.EQ.5 .OR. nbincomp(nkeep2).EQ.5) THEN
               number5isolated = number5isolated - 1
            ENDIF
         ENDIF

         cmx = pmassn(nkeep1)*xn(nkeep1) + pmassn(nkeep2)*xn(nkeep2)
         cmy = pmassn(nkeep1)*yn(nkeep1) + pmassn(nkeep2)*yn(nkeep2)
         cmz = pmassn(nkeep1)*zn(nkeep1) + pmassn(nkeep2)*zn(nkeep2)
         cmvx = pmassn(nkeep1)*vxn(nkeep1) + pmassn(nkeep2)*vxn(nkeep2)
         cmvy = pmassn(nkeep1)*vyn(nkeep1) + pmassn(nkeep2)*vyn(nkeep2)
         cmvz = pmassn(nkeep1)*vzn(nkeep1) + pmassn(nkeep2)*vzn(nkeep2)
         pmassmax(nkeep1) = MAX(pmassmax(nkeep1),pmassmax(nkeep2))
         pmassmin(nkeep1) = MIN(pmassmin(nkeep1),pmassmin(nkeep2))
         pmassn(nkeep1) = pmassn(nkeep1) + pmassn(nkeep2) 
         xn(nkeep1) = cmx/pmassn(nkeep1)
         yn(nkeep1) = cmy/pmassn(nkeep1)
         zn(nkeep1) = cmz/pmassn(nkeep1)
         vxn(nkeep1) = cmvx/pmassn(nkeep1)
         vyn(nkeep1) = cmvy/pmassn(nkeep1)
         vzn(nkeep1) = cmvz/pmassn(nkeep1)
         velcentre = SQRT(vxn(nkeep1)**2+vyn(nkeep1)**2+
     &        vzn(nkeep1)**2)*uvel
         WRITE (*,*) '      CM Velocity:     ',velcentre
         IF (nbincomp(nkeep1).EQ.2) THEN
            tempstring = name(nkeep1)(1:4)
            print *,'written ',tempstring
            READ (tempstring,*) iname1
            print *,'read'
            vrelb(iname1) = velcentre
            tempstring = name(nkeep2)(1:4)
            READ (tempstring,*) iname2
            vrelb(iname2) = velcentre
            WRITE (*,*) '      Gave com vel to  ',iname1,iname2
         ENDIF

c         write (tempstring,77001) nkeep2
         write (tempstring,77001) invert(nkeep2)
         name(nkeep1) = name(nkeep1)(1:lname(nkeep1)) // ',' // 
     &       name(nkeep2)(1:lname(nkeep2))
         lname(nkeep1) = lname(nkeep1) + 1 + lname(nkeep2)
c
c--Compact list
c
         ipos = 0
         DO n = 1, nnodes
            IF (n.NE.nkeep2) THEN
               ipos = ipos + 1
               xn(ipos) = xn(n)
               yn(ipos) = yn(n)
               zn(ipos) = zn(n)
               vxn(ipos) = vxn(n)
               vyn(ipos) = vyn(n)
               vzn(ipos) = vzn(n)
               spinx(ipos) = spinx(n)
               spiny(ipos) = spiny(n)
               spinz(ipos) = spinz(n)
               pmassn(ipos) = pmassn(n)
               pmassmax(ipos) = pmassmax(n)
               pmassmin(ipos) = pmassmin(n)
               name(ipos) = name(n)
               lname(ipos) = lname(n)
               nbincomp(ipos) = nbincomp(n)

               noutputline(ipos) = noutputline(n)

               interest(ipos) = interest(n)
               DO kkkk = 1, 3
                  orbit(kkkk,ipos) = orbit(kkkk,n)
               END DO
               seplast(ipos) = seplast(n)
               seplast2(ipos) = seplast2(n)
               qratiolast(ipos) = qratiolast(n)
               period(ipos) = period(n)
               periodratio(ipos) = period(n)
            ENDIF
         END DO
         nnodes = ipos
c
c--If last pair bound, get next most bound pair
c
         GOTO 1500
      ENDIF
c
c--Output systems
c
      DO n = 1, nnodes
         WRITE (51,98444) nbincomp(n),pmassn(n),pmassmax(n),pmassmin(n)
98444    FORMAT(I2,3(1X,1PE12.5))
      END DO
      WRITE (*,*) 'Smallest final sink separation: ',
     &     SQRT(radius2)*runit

      WRITE (*,*)
      WRITE (*,*) 'B  T  Q  Q'
      WRITE (*,*) number2, number3, number4, number5
      WRITE (*,*) number2unique, number3unique, number4unique,
     &     number5unique
      WRITE (*,*) number2isolated, number3isolated, number4isolated,
     &     number5isolated
      WRITE (*,*) 'BD  BBD SinBD BD3 BD4 BDin345'
      WRITE (*,*) numberbd, number2bd, numberbdsingle, numberbd3, 
     &     numberbd4, numberbd345
      WRITE (*,*) numberbdobs, number2bdobs, numberbdobssingle,
     &     numberbdobs3, numberbdobs4, numberbdobs345
      WRITE (*,*) '*  Bin* Sin* *+BD Trip* Quad* *in345 Mix3 Mix4'
      WRITE (*,*) numberstar, number2star, numberstarsingle,
     &     numberstarbd,numberstar3,numberstar4,numberstar345,
     &      numbermixed3,numbermixed4
      WRITE (*,*)

      WRITE (*,*)
      WRITE (*,*) 'B<0.01: ',number_0003,number2_0003,number3_0003,
     &     number4_0003
      WRITE (*,*) 'B<0.03: ',number_001,number2_001,number3_001,
     &     number4_001
      WRITE (*,*) 'B<0.07: ',number_003,number2_003,number3_003,
     &     number4_003
      WRITE (*,*) 'B<0.10: ',number_007,number2_007,number3_007,
     &     number4_007
      WRITE (*,*) 'B<0.20: ',number_01,number2_01,number3_01,
     &     number4_01
      WRITE (*,*) 'B<0.50: ',number_02,number2_02,number3_02,
     &     number4_02
      WRITE (*,*) 'B<0.80: ',number_05,number2_05,number3_05,
     &     number4_05
      WRITE (*,*) 'B<1.2: ',number_08,number2_08,number3_08,number4_08
      WRITE (*,*) 'B>1.2: ',number_11,number2_11,number3_11,number4_11
      WRITE (*,*)
      WRITE (*,*) 'B<0.20 (-BD): ',number_01notbd,number2_01notbd,
     &     number3_01notbd,number4_01notbd
      WRITE (*,*) 'B<0.50 (-BD): ',number_02notbd,number2_02notbd,
     &     number3_02notbd,number4_02notbd
      WRITE (*,*) 'B<0.80 (-BD): ',number_05notbd,number2_05notbd,
     &     number3_05notbd,number4_05notbd
      WRITE (*,*) 'B<1.2 (-BD): ',number_08notbd,number2_08notbd,
     &     number3_08notbd,number4_08notbd
      WRITE (*,*) 'B>1.2 (-BD): ',number_11notbd,number2_11notbd,
     &     number3_11notbd,number4_11notbd
      WRITE (*,*)

      WRITE (*,*) 'Binary separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary(i)
      END DO

      WRITE (*,*) 'Binary separation distribution (not BD comp)'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary 0.1-0.2 separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_01_02(i)
      END DO

      WRITE (*,*) 'Binary 0.2-0.5 separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_02_05(i)
      END DO

      WRITE (*,*) 'Binary 0.5-0.8 separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_05_08(i)
      END DO

      WRITE (*,*) 'Binary 0.8-1.2 separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_08_12(i)
      END DO

      WRITE (*,*) 'Binary >1.2 separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_12(i)
      END DO

      WRITE (*,*) 'Binary 0.1-0.2 separation distribution (not BD comp)'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_01_02_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary 0.2-0.5 separation distribution (not BD comp)'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_02_05_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary 0.5-0.8 separation distribution (not BD comp)'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_05_08_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary 0.8-1.2 separation distribution (not BD comp)'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_08_12_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary >1.2 separation distribution (not BD comp)'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary_12_notrealbd(i)
      END DO

      WRITE (*,*) 'Triple separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),ntriple(i)
      END DO

      WRITE (*,*) 'Quad separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nquad(i)
      END DO

      WRITE (*,*) 'Brown dwarf/VLM separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinarybd(i)
      END DO

      WRITE (*,*) bdmass,'-0.03 Brown dwarf/VLM separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinarybdobs(i)
      END DO

      WRITE (*,*) 'Triple brown dwarf/VLM separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),ntriplebd(i)
      END DO
      WRITE (*,*)

      WRITE (*,*) 'Binary mass ratio distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq(i)
      END DO

      WRITE (*,*) 'Binary mass ratio distribution (not BD comp)'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.1-0.2 distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_01_02(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.2-0.5 distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_02_05(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.5-0.8 distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_05_08(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.8-1.2 distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_08_12(i)
      END DO

      WRITE (*,*) 'Binary mass ratio >1.2 distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_12(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.1-0.2 distribution (not BD comp)'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_01_02_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.2-0.5 distribution (not BD comp)'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_02_05_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.5-0.8 distribution (not BD comp)'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_05_08_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary mass ratio 0.8-1.2 distribution (not BD comp)'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_08_12_notrealbd(i)
      END DO

      WRITE (*,*) 'Binary mass ratio >1.2 distribution (not BD comp)'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq_12_notrealbd(i)
      END DO

      WRITE (*,*) 'Triple mass ratio distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),ntripleq(i)
      END DO

      WRITE (*,*) 'Brown dwarf/VLM mass ratio distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinarybdq(i)
      END DO

      WRITE (*,*) bdmass,'-0.03 Brown dwarf/VLM mass ratio distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinarybdobsq(i)
      END DO

      OPEN (15,FILE='discveldisp')
      WRITE (15,99101)
99101 FORMAT('# Num Time           Free-Fall Time  Mass            ',
     &   'Truncation      Vel Disp        Vel Disp        ')
      WRITE (15,99100) time,fftime
      i = nlines
      DO k = 1, nptmass
         n = listmerge(k)
         WRITE(15, 99013) n,tffform(n),tform(n),pmass(n), 
     &        radiusmin(n)*runit,vrel(n),vrelb(n),radcom(n)*runit,
     &        radmax(n)*runit
99013    FORMAT (I4,1X,8(1PE15.8,1X))
      END DO
      CLOSE(15)

      DO i = 1, nmultipleout
      WRITE (50,99444) imultiple_values(1,i),xmultiple_values(1,i),
     &        xmultiple_values(2,i),imultiple_values(2,i),
     &        imultiple_values(3,i),(xmultiple_values(j,i),j=3,13),
     &        imultiple_values(4,i)
c      WRITE (50,99444) nbincomp(nkeep1),xmassmax,xmassmin,
c     &           interest(nkeep1),interest(nkeep2),
c     &           pmassn(nkeep1),pmassn(nkeep2),
c     &           qratio_out,semimajorkeep*runit,eccentricitykeep,
c     &           angleorbit,angleorbit1,angleorbit2,anglerelative,
c     &           periodkeep,periodratiokeep
      END DO

      STOP
      END
