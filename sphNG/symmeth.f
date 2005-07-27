      PROGRAM symmeth

      INCLUDE 'idim'

      REAL*8 udist, umass, utime

      COMMON /part / npart, x(idim), y(idim), z(idim), vx(idim),
     &               vy(idim), vz(idim), u(idim), time
      COMMON /dens / rho(idim), pmass(idim), dgrav(idim)
      COMMON /kerne/ h(idim)
      COMMON /aux  / gamma
      COMMON /gridv/ iline, isec
      COMMON /bodys/ n1, n2
      COMMON /origi/ cmx, cmy, cmz
      COMMON /files/ fileout, filein, inname
      COMMON /optio/ what, idump, icont
      COMMON /forma/ ien
      COMMON /rhonm/ rhozero
      COMMON /ptmass/ nptmass,spinx(iptdim),spiny(iptdim),
     &           spinz(iptdim),listpm(iptdim),iptmass,hacc,ptmcrit,radcrit
      COMMON /phase/ iphase(idim)
      COMMON /gho  / ireal(idim)
      COMMON /timei/ isteps(idim)
      COMMON /tor/   angaddx(iptdim), angaddy(iptdim),angaddz(iptdim),
     &     spinadx(iptdim),spinady(iptdim),spinadz(iptdim),torqt(idim),
     &     torqg(idim),torqp(idim),torqv(idim),torqc(idim)

      DIMENSION hnew(idim)

      CHARACTER*7 filein
      CHARACTER*11 inname, imageout
      CHARACTER*11 fileout, fileout2, fileout3
      CHARACTER*9  name
      CHARACTER*1  ien, icont, icor, icent, isink, itorq

      WRITE (*,*) 'Enter name of input file:'
         READ(*,90) filein
 90      FORMAT(A7)

      OPEN (UNIT = 11, FILE = filein, FORM = 'unformatted',
     &     RECL=maxrec)

      WRITE (*,*) 'Reading input file'

      READ (11, END = 20) udist, umass, utime, npart, n1, n2,
     &        time,gamma,rhozero,K2,(h(i),i=1,npart),escap,tkin,tgrav,
     &        tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,npart),
     &        (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     &        (u(i),i=1,npart),(pmass(i),i=1,npart),(rho(i),i=1,
     &        npart),(dgrav(i),i=1,npart),dtmaxdp,
     &        (isteps(i), i=1, npart),
     &        (iphase(i),i=1,npart),nptmass,
     &        (listpm(i),i=1,nptmass),(spinx(i),i=1,nptmass),
     &        (spiny(i),i=1,nptmass),(spinz(i),i=1,nptmass),
     &        (angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass),
     &        (angaddz(i),i=1,nptmass),
     &        anglostx, anglosty, anglostz,
     &        nreassign, naccrete, nkill, specang, ptmassin,
     &        (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass),
     &        (spinadz(i),i=1,nptmass),
     &        (torqt(i), i=1, npart), (torqg(i), i=1, npart),
     &        (torqp(i), i=1, npart),(torqv(i), i=1, npart)

 20   CLOSE(11)

      WRITE (*,*) 'Calculating h(i)s'

      iout = 0
      DO i = 1, npart
         IF (MOD(i,50).EQ.0) THEN
            iout = 1
            WRITE (*,*) 'P i=',i,' h: ',h(i)
         ENDIF
         x1 = ABS(x(i))
         y1 = ABS(y(i))
         z1 = ABS(z(i))
         hnew(i) = 0.
         icount = 0
         DO j = 1, npart
            x2 = ABS(x(j))
            y2 = ABS(y(j))
            z2 = ABS(z(j))
            IF (ABS(x1-x2).LT.0.00001 .AND. ABS(y1-y2).LT.0.00001 
     &           .AND. ABS(z1-z2).LT.0.00001) THEN
               icount = icount + 1
               hnew(i) = hnew(i) + h(j)
               IF (iout.EQ.1) WRITE (*,*) '  j=',j,' h: ',h(j)
            ENDIF
         END DO
         iout = 0
         IF (icount.NE.8) THEN
            WRITE (*,*) icount
            STOP
         ENDIF
         hnew(i) = hnew(i)/icount
      END DO

      DO i = 1, npart
         h(i) = hnew(i)
      END DO

      OPEN (UNIT = 11, FILE = 'OUTPUT1', FORM = 'unformatted',
     &     RECL=maxrec)

      WRITE (11, ERR = 30) udist, umass, utime, npart, n1, n2,
     &        time,gamma,rhozero,K2,(h(i),i=1,npart),escap,tkin,tgrav,
     &        tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,npart),
     &        (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     &        (u(i),i=1,npart),(pmass(i),i=1,npart),(rho(i),i=1,
     &        npart),(dgrav(i),i=1,npart),dtmaxdp,
     &        (isteps(i), i=1, npart),
     &        (iphase(i),i=1,npart),nptmass,
     &        (listpm(i),i=1,nptmass),(spinx(i),i=1,nptmass),
     &        (spiny(i),i=1,nptmass),(spinz(i),i=1,nptmass),
     &        (angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass),
     &        (angaddz(i),i=1,nptmass),
     &        anglostx, anglosty, anglostz,
     &        nreassign, naccrete, nkill, specang, ptmassin,
     &        (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass),
     &        (spinadz(i),i=1,nptmass),
     &        (torqt(i), i=1, npart), (torqg(i), i=1, npart),
     &        (torqp(i), i=1, npart),(torqv(i), i=1, npart),
     &        (torqc(i), i=1, npart)

 30   CLOSE(11)

      END
