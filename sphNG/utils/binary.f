      PROGRAM binary
c
c--Calculates the orbital parameters and mass accretion of two point masses
c     using a BINARY output 'P file'
c
      IMPLICIT DOUBLE PRECISION (D)
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL (A-C, E-H, O-Z)

      PARAMETER (idim = 4)
      DIMENSION ipart(idim), x(idim), y(idim), z(idim), vx(idim), 
     &     vy(idim), vz(idim), pmass(idim), spinx(idim), spiny(idim), 
     &     spinz(idim)
      DIMENSION oldmasses(idim,10000), oldtimes(idim,10000), 
     &     accrete(idim), acctime(idim)

      CHARACTER*21 infile, outfile

      WRITE (*,*) 'Enter input filename:' 
      READ (*,99001) infile
99001 FORMAT(A21)
      WRITE (*,*) 'Enter output filename:' 
      READ (*,99001) outfile

      WRITE (*,*) 'Number dumps to average accretion over (max 10000)?'
      READ (*,*) ndumps
      DO i = 1, ndumps
         oldmasses(1,i) = 0.
         oldtimes(1,i) = 0.
         oldmasses(2,i) = 0.
         oldtimes(2,i) = 0.
      END DO

      iupto = 0
      imorethantwo = 0

      OPEN (15, FILE=infile, STATUS='unknown', FORM='unformatted')
      OPEN (16, FILE=outfile)

 100  READ (15, END=200) fftime, time, nptmass
      IF (imorethantwo.EQ.0 .AND. nptmass.GT.2) THEN
         WRITE (*,*) 'MORE THAN TWO POINT MASSES'
         imorethantwo = 1
      ENDIF
      IF (nptmass.GT.idim) THEN
         WRITE (*,*) 'GREATER THAN ',idim,' POINT MASSES'
         STOP
      ENDIF

      DO i = 1, nptmass
         READ (15, END=200) ipart(i), x(i), y(i), z(i), vx(i), vy(i), 
     &        vz(i), pmass(i), rho, nactotal, ptmassinner, spinx(i),
     &        spiny(i), spinz(i)
      END DO
      
      IF (nptmass.GE.2) THEN
         IF (iupto.LT.ndumps) iupto = iupto + 1
         DO i = 1, ndumps - 1
            DO j = 1, 2
               oldmasses(j, i) = oldmasses(j, i + 1)
               oldtimes(j, i) = oldtimes(j, i + 1)
            END DO
         END DO
         DO j = 1, 2
            oldmasses(j, ndumps) = pmass(j)
            oldtimes(j, ndumps) = time

            IF (iupto.GT.1) THEN
               accrete(j) = (pmass(j)-oldmasses(j,ndumps-iupto+1))/
     &              (time-oldtimes(j,ndumps-iupto+1))
               acctime(j) = (time-oldtimes(j,ndumps-iupto+1))/2.0 + 
     &              oldtimes(j,ndumps-iupto+1)
            ELSE
               accrete(j) = 0.0
               acctime(j) = 0.0
            ENDIF
         END DO

         dq = pmass(2)/pmass(1)
         dpmasstot = pmass(1) + pmass(2)
         dx = x(1) - x(2)
         dy = y(1) - y(2)
         dz = z(1) - z(2)
         dvx = vx(1) - vx(2)
         dvy = vy(1) - vy(2)
         dvz = vz(1) - vz(2)
         dr = SQRT(dx*dx + dy*dy + dz*dz)

         denergy = - dpmasstot/dr + (dvx*dvx + dvy*dvy + dvz*dvz)/2.0
         dangx = dvy*dz - dvz*dy
         dangy = dvz*dx - dvx*dz
         dangz = dvx*dy - dvy*dx
         dang2 = (dangx*dangx + dangy*dangy + dangz*dangz)

         da = - dpmasstot/2.0/denergy

         decc = 1.0D+0 - dang2/dpmasstot/da
         IF (decc.LT.-0.000001) THEN
            WRITE (*,*) 'ERROR: ecc^2 = ',decc
         ELSEIF (decc.LT.0.0) THEN
            decc = 0.0
         ENDIF
         decc = SQRT(decc)

         IF (accrete(1).GT.0.0) THEN
            dm2dm1 = accrete(2)/accrete(1)
         ELSE
            dm2dm1 = 0.0
         ENDIF
         
         WRITE(16,99005) time, fftime, pmass(1)+pmass(2), 
     &        pmass(1), pmass(2), dq, da, decc,
     &        acctime(1), accrete(1), accrete(2), dm2dm1, dr, spinx(1),
     &        spiny(1), spinz(1), spinx(2), spiny(2), spinz(2)
99005    FORMAT (19(1PE15.8,1X))
      ENDIF

      GOTO 100

 200  CONTINUE

      WRITE(*,*) 'Times: ',fftime, time, ' Number ptmass: ',nptmass

      STOP
      END
