      SUBROUTINE sphmas
c************************************************************
c                                                           *
c  This subroutine gives particles masses in a spherical    *
c     coordinate distribution                               *
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

      CHARACTER*1 iok,prof
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine cartmas ')
99004 FORMAT (A1)
c
c--Set Condensed Areas In Spherical Coordinates
c
      fractot = 0.
      rcyl2 = rcyl*rcyl
      zmax2 = zmax*zmax
      deltaz = zmax - zmin

      WRITE (*, 99025)
99025 FORMAT (' Do you want sinusoidal mass perturbations ? (s)',/,
     &        '                       or a r mass profile ? (p)',/,
     &        '                or a linear z mass profile ? (z)',/,
     &        '             or a cos(mtheta) perturbation ? (c)',/,
     &        ' or a cos(mtheta)exp(-(k(r-r1)/r1)^2 pert. ? (e)')
      READ (*, 99004) iok

      IF (iok.EQ.'s') THEN
         npert = 0
         densc = 1.
         WRITE (*, 99026)
99026    FORMAT (' Enter number of perturbations',
     &        ' and the density contrast')   
         READ (*, *) npert, densc
         densc1 = densc - 1
         DO i = nptmass + 1, npart
            r1 = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2)
            z1 = (xyzmh(3,i) - zmin) / deltaz
            disfrac(i) = disfrac(i)*(1.0 + densc1*(SIN(pi*npert*z1))**2
     &           * ((1.0 - r1/rcyl)**2))
            fractot = fractot + disfrac(i)
         END DO

      ELSEIF (iok.EQ.'p') THEN
         kprof = 1
         WRITE (*, 99029)
99029    FORMAT(' Enter profile along r ',/,
     &          '                  r^-0 : 0',/, 
     &          '                  r^-1 : 1',/,
     &          '                  r^-2 : 2',/,
     &          '           exponential : e')   
         READ (*, 99004) prof

         IF (prof.NE.'e') THEN
            iprofr = 0
            IF (prof.EQ.'1') iprofr = 1
            IF (prof.EQ.'2') iprofr = 2
            profr = iprofr/2.
            DO i = nptmass + 1, npart
               ri2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
               disfrac(i) = disfrac(i)*((rcyl2/ri2)**profr) 
               fractot = fractot + disfrac(i)
            END DO
         ELSE
            rcyl20 = rcyl2 * 4.0 / 9.0
            DO i = nptmass + 1, npart
               ri2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
               disfrac(i) =  disfrac(i)*(10.0 * EXP(-1.0*ri2/rcyl20))
               fractot = fractot + disfrac(i)
            END DO
         ENDIF

      ELSEIF (iok.EQ.'z') THEN
         WRITE (*,99031)
99031    FORMAT(' Enter density contrast ')
         READ (*,*) densc
         densc = densc - 1
         DO i = nptmass + 1, npart
            zi = (xyzmh(3,i) - zmin) / deltaz
            disfrac(i) = disfrac(i)*(1.0 + zi*densc)
            fractot = fractot + disfrac(i)
         END DO

      ELSEIF (iok.EQ.'c') THEN
         WRITE (*,99034)
99034    FORMAT(' Enter m, and density contrast ')
         READ(*,*) m, densc
         DO i = nptmass + 1, npart
            theta = ATAN2(xyzmh(2,i), xyzmh(1,i))
            disfrac(i) = disfrac(i)*(1.0 + densc*COS(m*theta))
            fractot = fractot + disfrac(i)
         END DO

      ELSEIF (iok.EQ.'e') THEN
         WRITE (*,99035)
99035    FORMAT(' Enter m, and density contrast ')
         READ(*,*) m, densc
         WRITE (*,99036)
99036    FORMAT(' Enter r1 and k ')
         READ(*,*) rad1, width
         DO i = nptmass + 1, npart
            theta = ATAN2(xyzmh(2,i), xyzmh(1,i))
            r1 = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2)
            disfrac(i) = disfrac(i)*(1.0 + densc*COS(m*theta)*
     &           EXP(-(width*(r1-rad1)/rad1)**2.0))
            fractot = fractot + disfrac(i)
         END DO
      ENDIF

      fractot = fractot / FLOAT(npart - nptmass)
      WRITE (*,*) fractot
      DO i = nptmass + 1, npart
         disfrac(i) = disfrac(i)/fractot
      END DO

      RETURN
      END
