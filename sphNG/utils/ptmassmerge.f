      PROGRAM ptmass
c
c--Calculates the orbits and mass accretion of all point masses
c     using a BINARY output 'P file'
c
      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL (A-C, E-H, O-Z)

      PARAMETER (idim = 1000)
      PARAMETER (ilines = 100000)
      PARAMETER (nbinimf = 100)
      PARAMETER (nbinbin = 40)
      PARAMETER (nmassratio = 10)
      DIMENSION ipart(idim,ilines), x(idim,ilines), y(idim,ilines), 
     &     z(idim,ilines), vx(idim,ilines), 
     &     vy(idim,ilines), vz(idim,ilines), pmass(idim,ilines), 
     &     rho(idim,ilines), nactotal(idim,ilines),
     &     ptmassinner(idim,ilines), spinx(idim,ilines), 
     &     spiny(idim,ilines), spinz(idim,ilines)
      DIMENSION xn(idim), yn(idim), zn(idim), vxn(idim), 
     &     vyn(idim), vzn(idim), pmassn(idim), name(idim), lname(idim),
     &     nbincomp(idim),
     &     tform(idim),tffform(idim),vrel(idim)
      DIMENSION oldmasses(idim,10000), oldtimes(idim,10000), 
     &     accrete(idim,ilines), acctime(idim,ilines), iupto(idim),
     &     time(ilines),fftime(ilines),radiusmin(idim,ilines),
     &     radiuscreate(idim,idim),radiusend(idim,idim)
      DIMENSION nimf(nbinimf),nimfa(nbinimf)
      DIMENSION nbinary(nbinbin),nbinarybd(nbinbin)
      DIMENSION nbinaryq(nmassratio),nbinarybdq(nmassratio)
      DIMENSION listmerge(idim)

      CHARACTER*21 infile, outfile
      CHARACTER*19 outbase
      CHARACTER*3 num,tempstring
      CHARACTER*1 ans,iout
      CHARACTER*60 name

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
c
c--Compute
c
      nlines = 0
      nptmass = 0
      nmerge = 0

      DO i = 1, idim
         listmerge(i) = i
      END DO

      OPEN (22, FILE='MERGED23', STATUS='unknown')
      OPEN (15, FILE=infile, STATUS='unknown', FORM='unformatted')
 100  nlines = nlines + 1
      IF (nlines.GT.ilines) THEN
         WRITE (*,*) 'MAX NUMBER OF LINES REACHED'
         GOTO 200
      ENDIF
c
c--Read time and number of point masses
c
 110  nptmassold = nptmass
 111  READ (15, END=200, ERR=200) 
     &     fftime(nlines), time(nlines), nptmass
c         write (*,*) fftime(nlines), time(nlines), nptmass
      IF (nptmass.LT.0) THEN
         READ (15) imerged1, imerged2
         WRITE (*,*) 'Merged ',imerged1, imerged2
c         READ (22,*,END=112) imerged1, imerged2
         WRITE (*,*) 'Merged file ',imerged1, imerged2
         nmerge = nmerge + 1
         nptmassold = nptmassold - 1
         pmass(listmerge(imerged2),nlines) = -1.0
         DO i = imerged2, idim
            listmerge(i) = listmerge(i+1)
         END DO
         GOTO 111
 112     nptmass = nptmassold
         GOTO 200
      ENDIF
      IF (nptmassold .NE. nptmass) THEN
         DO i=listmerge(nptmassold+1),listmerge(nptmass)
            tform(i) = fftime(nlines)
            tffform(i) = time(nlines)
         END DO
         write (*,*) fftime(nlines), time(nlines), nptmass
      ENDIF

      IF (nptmass+nmerge.GT.idim) THEN
      WRITE (*,*) 'GREATER THAN ',idim,' POINT MASSES ',nptmass,nmerge
         STOP
      ENDIF
c
c--Read each time dump
c
      DO k = 1, nptmass
         i = listmerge(k)
         READ (15, END=200, ERR=200) ipart(i,nlines), x(i,nlines), 
     &        y(i,nlines), z(i,nlines), vx(i,nlines), 
     &        vy(i,nlines), vz(i,nlines), pmass(i,nlines), 
     &        rho(i,nlines), nactotal(i,nlines),
     &        ptmassinner(i,nlines), 
     &        spinx(i,nlines), spiny(i,nlines), spinz(i,nlines)

c         write (*,*) ipart(i,nlines), x(i,nlines), y(i,nlines)

         IF (iupto(i).LT.ndumps) iupto(i) = iupto(i) + 1

         r = SQRT(x(i,nlines)**2 + y(i,nlines)**2)
         th = ATAN2(y(i,nlines),x(i,nlines))
         th = th + time(nlines)*omega
         x(i,nlines) = r*COS(th)
         y(i,nlines) = r*SIN(th)

         r = SQRT(vx(i,nlines)**2 + vy(i,nlines)**2)
         th = ATAN2(vy(i,nlines),vx(i,nlines))
         th = th + time(nlines)*omega
         vx(i,nlines) = r*COS(th)
         vy(i,nlines) = r*SIN(th)

         DO j = 1, ndumps - 1
            oldmasses(i, j) = oldmasses(i, j + 1)
            oldtimes(i, j) = oldtimes(i, j + 1)
         END DO
         oldmasses(i, ndumps) = pmass(i,nlines)
         oldtimes(i, ndumps) = time(nlines)

         IF (iupto(i).GT.1) THEN
            accrete(i,nlines) = (pmass(i,nlines)-
     &           oldmasses(i,ndumps-iupto(i)+1))/
     &           (time(nlines)-oldtimes(i,ndumps-iupto(i)+1))
            acctime(i,nlines) = (time(nlines)-
     &           oldtimes(i,ndumps-iupto(i)+1))/2.0 + 
     &           oldtimes(i,ndumps-iupto(i)+1)
         ELSE
            accrete(i,nlines) = 0.0
            acctime(i,nlines) = 0.0
         ENDIF
 
      END DO
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
         vcmxtot = vcmxtot + pmass(n,nlines)*vx(n,nlines)
         vcmytot = vcmytot + pmass(n,nlines)*vy(n,nlines)
         vcmztot = vcmztot + pmass(n,nlines)*vz(n,nlines)
         DO kk = 1, nptmass
            n2 = listmerge(kk)
            IF (n.NE.n2) THEN
               rx = x(n,nlines)-x(n2,nlines)
               ry = y(n,nlines)-y(n2,nlines)
               rz = z(n,nlines)-z(n2,nlines)
               r2 = (rx**2 + ry**2 + rz**2)
               radius2 = MIN(radius2, r2)
               radiusend(n,n2) = SQRT(r2)
               IF (k.GT.nptmassold) THEN
                  radiuscreate(n,n2) = radiusend(n,n2)
               ENDIF
            ENDIF
         END DO         
         IF (nlines .NE. 1) THEN
            IF (k.GT.nptmassold) THEN
               radiusmin(n,nlines-1)=1.0E+15
               write (*,*) 'Sep ',SQRT(radius2)*runit
            ENDIF

            IF (radius2 .LT. radiusmin(n,nlines-1)**2) THEN
               radiusmin(n,nlines) = SQRT(radius2)
            ELSE
               radiusmin(n,nlines) = radiusmin(n,nlines-1)
            ENDIF
         ELSE
            radiusmin(n,nlines) = 1.0E+15
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

      DO i = 1, nbinimf 
         nimf(i) = 0
         nimfa(i) = 0
      END DO
      DO i = 1, nbinbin
         nbinary(i) = 0
         nbinarybd(i) = 0
      END DO
      DO i = 1, nmassratio
         nbinaryq(i) = 0
         nbinarybdq(i) = 0
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
         vrelx = vx(n,nlines) - vcmxtot
         vrely = vy(n,nlines) - vcmytot
         vrelz = vz(n,nlines) - vcmztot
         vrel(n) = sqrt(vrelx**2+vrely**2+vrelz**2)*uvel
         write (*,*) 'Velocity ',n,': ',vrel(n)

         ipos = INT((ALOG10(pmass(n,nlines))-1)*5.0+nbinimf)
         IF (ipos.GE.1 .AND. ipos.LE.nbinimf) THEN
            nimf(ipos) = nimf(ipos) + 1
            IF (accrete(n,nlines).NE.0.0) 
     &        nimfa(ipos) = nimfa(ipos) + 1
         ENDIF
      END DO
c
c--For each point mass, output evolution file
c
      IF (iout.EQ.'Y' .OR. iout.EQ.'y') THEN
         DO n = 1, nptmass + nmerge
            IF (n.GT.99) THEN
               WRITE (num, 88000) n
88000          FORMAT (I3)
            ELSEIF (n.GT.9) THEN
               WRITE (num, 88001) n
88001          FORMAT ('0',I2)
            ELSE 
               WRITE (num, 88002) n
88002          FORMAT ('00',I1)
            ENDIF
            DO junk = LEN(outbase), 1, -1
               IF (outbase(junk:junk).NE.' ') THEN
                  iplace = junk
                  GOTO 150
               ENDIF
            END DO
 150        outfile = outbase(1:junk) // num

            OPEN (15, FILE=outfile)
            DO i = 1, nlines
               IF (pmass(n,i).GT.0.) THEN
                  WRITE(15, 99005) time(i), fftime(i), 
     &                 x(n,i), y(n,i), z(n,i), vx(n,i), 
     &           vy(n,i), vz(n,i), pmass(n,i), rho(n,i), acctime(n,i),
     &           accrete(n,i), spinx(n,i), spiny(n,i), spinz(n,i),
     &           radiusmin(n,i)*runit
99005             FORMAT (16(1PE15.8,1X))
               ENDIF
            END DO
            CLOSE (15)
         END DO
      ENDIF
c
c--For last dump, output endfile containing data for each point mass
c
      OPEN (15,FILE='endfile')
      WRITE (15,99100) time(nlines),fftime(nlines)
99100 FORMAT('# ',2(1PE15.8,1X))
      i = nlines
      DO k = 1, nptmass
         n = listmerge(k)
         WRITE(15, 99012) tffform(n),tform(n),
     &        x(n,i), y(n,i), z(n,i), vx(n,i),
     &        vy(n,i), vz(n,i), pmass(n,i), 
     &        accrete(n,i), spinx(n,i), spiny(n,i), spinz(n,i),
     &        radiusmin(n,i)*runit,vrel(n)
99012    FORMAT (15(1PE15.8,1X))
      END DO
      CLOSE(15)
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
         DO kk = 1, n-1
            n2 = listmerge(kk)
            IF (radiuscreate(n,n2).LT.sepmin) THEN
               sepmin = radiuscreate(n,n2)
               nkeep = n2
            ENDIF
         END DO
         WRITE(15, *) sepmin*runit, radiusend(n,nkeep)*runit
      END DO
      CLOSE(15)
      OPEN (15,FILE='discveldisp')
      WRITE (15,99101)
99101 FORMAT('# Num Time           Free-Fall Time  Mass            ',
     &   'Truncation      Vel Disp        Vel Disp        ')
      WRITE (15,99100) time(nlines),fftime(nlines)
      i = nlines
      DO k = 1, nptmass
         n = listmerge(k)
         WRITE(15, 99013) n,tffform(n),tform(n),pmass(n,i), 
     &        radiusmin(n,i)*runit,vrel(n),vrel(n)
99013    FORMAT (I3,1X,6(1PE15.8,1X))
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
      numberstar = 0
      number2 = 0
      number3 = 0
      number4 = 0
      number5 = 0
      number2bd = 0
      number3bd = 0
      number4bd = 0
      number5bd = 0
      number2star = 0
      numberstarbd = 0
      DO k = 1, nptmass
         n = listmerge(k)
         xn(k) = x(n,nlines)
         yn(k) = y(n,nlines)
         zn(k) = z(n,nlines)
         vxn(k) = vx(n,nlines)
         vyn(k) = vy(n,nlines)
         vzn(k) = vz(n,nlines)
         pmassn(k) = pmass(n,nlines)
         IF (pmassn(k).LT.0.09) THEN
            numberbd= numberbd + 1
         ELSE
            numberstar= numberstar + 1
         ENDIF
         write (tempstring,77001) n
         name(k) = tempstring
         lname(k) = 3
         nbincomp(k) = 1
77001    FORMAT(I3)
      END DO
      nnodes = nptmass
      numberbdsingle = numberbd
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
            IF (etot.LT.0.0 .AND. radius.LT.rmin) THEN
               rmin = MIN(rmin,radius)
               emin = etot
               nkeep1 = n
               nkeep2 = n2
               semimajorkeep = semimajor
               eccentricitykeep = eccentricity
            ENDIF
            
         END DO
ccc         WRITE (99,*) x(n,nlines),y(n,nlines),z(n,nlines)
      END DO
      WRITE (*,*) 'emin ',emin
      WRITE (*,*) 'rmin ',rmin
      IF (emin.LT.0.0) THEN
         qratio = MIN(pmassn(nkeep1)/pmassn(nkeep2), 
     &        pmassn(nkeep2)/pmassn(nkeep1))

         WRITE (*,*) name(nkeep1), ' is bound to ',name(nkeep2)
         WRITE (*,*) '      Semi-major axis: ',semimajorkeep*runit
         WRITE (*,*) '      Eccentricity:    ',eccentricitykeep
         WRITE (*,*) '      Mass Ratio:      ',qratio
         WRITE (*,*) '      Masses:          ',pmassn(nkeep1),
     &        pmassn(nkeep2)
c
c--Merge binaries into 1 node
c
         nbincomp(nkeep1) = nbincomp(nkeep1) + nbincomp(nkeep2)
         IF (nbincomp(nkeep1).EQ.2) THEN
            number2 = number2 + 1
            ipos = INT((ALOG10(semimajorkeep*runit))*5.0+16)
            IF (ipos.GE.1 .AND. ipos.LE.nbinbin) THEN
               nbinary(ipos) = nbinary(ipos) + 1
               IF(pmassn(nkeep1).LT.0.09.AND.pmassn(nkeep2).LT.0.09)THEN
                  number2bd = number2bd + 1
                  numberbdsingle = numberbdsingle - 2
                  nbinarybd(ipos) = nbinarybd(ipos) + 1
                  WRITE (*,*) 'BBD'
           ELSEIF (pmassn(nkeep1).LT.0.09.OR.pmassn(nkeep2).LT.0.09)THEN
                  numberstarbd = numberstarbd + 1
                  numberbdsingle = numberbdsingle - 1
                  numberstarsingle = numberstarsingle - 1
                  WRITE (*,*) 'S-BD'
               ELSE
                  number2star = number2star + 1
                  numberstarsingle = numberstarsingle - 2
               ENDIF
            ENDIF
            ipos = MAX(1,MIN(INT(qratio*nmassratio+1),nmassratio))
            nbinaryq(ipos) = nbinaryq(ipos) + 1
            IF (ipos.EQ.1) WRITE (*,*) 'EXTREME Q'
            IF(pmassn(nkeep1).LT.0.09 .AND. pmassn(nkeep2).LT.0.09)THEN
               nbinarybdq(ipos) = nbinarybdq(ipos) + 1
            ENDIF
         ELSEIF (nbincomp(nkeep1).EQ.3) THEN
            number3 = number3 + 1
         ELSEIF (nbincomp(nkeep1).EQ.4) THEN
            number4 = number4 + 1
         ELSEIF (nbincomp(nkeep1).EQ.5) THEN
            number5 = number5 + 1
         ENDIF

         cmx = pmassn(nkeep1)*xn(nkeep1) + pmassn(nkeep2)*xn(nkeep2)
         cmy = pmassn(nkeep1)*yn(nkeep1) + pmassn(nkeep2)*yn(nkeep2)
         cmz = pmassn(nkeep1)*zn(nkeep1) + pmassn(nkeep2)*zn(nkeep2)
         cmvx = pmassn(nkeep1)*vxn(nkeep1) + pmassn(nkeep2)*vxn(nkeep2)
         cmvy = pmassn(nkeep1)*vyn(nkeep1) + pmassn(nkeep2)*vyn(nkeep2)
         cmvz = pmassn(nkeep1)*vzn(nkeep1) + pmassn(nkeep2)*vzn(nkeep2)
         pmassn(nkeep1) = pmassn(nkeep1) + pmassn(nkeep2) 
         xn(nkeep1) = cmx/pmassn(nkeep1)
         yn(nkeep1) = cmy/pmassn(nkeep1)
         zn(nkeep1) = cmz/pmassn(nkeep1)
         vxn(nkeep1) = cmvx/pmassn(nkeep1)
         vyn(nkeep1) = cmvy/pmassn(nkeep1)
         vzn(nkeep1) = cmvz/pmassn(nkeep1)
         WRITE (*,*) '      CM Velocity:     ',
     &      SQRT(vxn(nkeep1)**2+vyn(nkeep1)**2+vzn(nkeep1)**2)*uvel
         write (tempstring,77001) nkeep2
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
               pmassn(ipos) = pmassn(n)
               name(ipos) = name(n)
               lname(ipos) = lname(n)
               nbincomp(ipos) = nbincomp(n)
            ENDIF
         END DO
         nnodes = ipos
c
c--If last pair bound, get next most bound pair
c
         GOTO 1500
      ENDIF
      WRITE (*,*) 'Smallest final sink separation: ',
     &     SQRT(radius2)*runit

      WRITE (*,*)
      WRITE (*,*) 'B  T  Q  Q'
      WRITE (*,*) number2, number3, number4, number5
      WRITE (*,*) 'BD  BBD SinBD'
      WRITE (*,*) numberbd, number2bd, numberbdsingle
      WRITE (*,*) '*  Bin* Sin* *+BD'
      WRITE (*,*) numberstar, number2star, numberstarsingle,numberstarbd
      WRITE (*,*)

      WRITE (*,*) 'Binary separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinary(i)
      END DO
      WRITE (*,*) 'Brown dwarf/VLM separation distribution'
      DO i = 1, nbinbin
         WRITE (*,*) 10**((i-16)/5.0),nbinarybd(i)
      END DO
      WRITE (*,*)

      WRITE (*,*) 'Binary mass ratio distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinaryq(i)
      END DO
      WRITE (*,*) 'Brown dwarf/VLM separation distribution'
      DO i = 1, nmassratio
         WRITE (*,*) (i-1)/REAL(nmassratio),nbinarybdq(i)
      END DO

      STOP
      END
