      SUBROUTINE ktable
c************************************************************
c                                                           *
c  This subroutine builds a table for the various values    *
c     of the kernel, the gradient of the kernel, the mass   *
c     fraction and the potential energy.                    *
c     The entry is v**2.                                    *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      REAL*8 sum, v, v2, v3, v4, v5, v6, v7

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ktable')
c
c--Maximum interaction length and step size
c
      radkernel = 2.5
      part1kernel = 0.5
      part2kernel = 1.5
      v2max = radkernel*radkernel
      dvtable = v2max/itable
      i1 = part1kernel/dvtable
      i2 = part2kernel/dvtable
c
c--Build tables
c
c  a) v less than 1
c
      DO i = 0, i1
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         v7 = v*v6
         sum = 1./40.*(575./8. - 105.*v2 + 54.*v4)
         wij(i) = sum
         sum = 1./40.*(-210.*v + 216.*v3)
         grwij(i) = sum
         sum = 1./10.*(575./24.*v3 - 21.*v5 + 54./7.*v7)
         fmass(i) = sum
         sum = 1./10.*(575./48.*v2 - 21./4.*v4 + 9./7.*v6 - 1199./64.)
         fpoten(i) = sum
         sum = 1./10.*(1199./64. - 575./16.*v2 + 105./4.*v4 - 9.*v6)
         dphidh(i) = sum
      END DO
c
c  b) v greater than 0.5
c
      DO i = i1 + 1, i2
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         v7 = v*v6
         sum = 1./40.*(275./4. + 30.*v - 210.*v2 + 160.*v3 - 36.*v4)
         wij(i) = sum
         sum = 1./40.*(30. - 420.*v + 480.*v2 - 144.*v3)
         grwij(i) = sum
         sum = 1./10.*(275./12.*v3 + 15./2.*v4 - 42.*v5 + 80./3.*v6 
     &        - 36./7.*v7 + 1./672.)
         fmass(i) = sum
         sum = 1./10.*(275./24.*v2 + 5./2.*v3 - 21./2.*v4 + 16./3.*v5 
     &        - 6./7.*v6 - 599./32.)
         fpoten(i) = sum
         sum = 1./10.*(599./32. - 275./8.*v2 - 
     &        10.*v3 + 105./2.*v4 - 32.*v5 + 6.*v6)
         dphidh(i) = sum
      END DO
c
c  c) v greater than 1.5
c
      DO i = i2 + 1, itable
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         v7 = v*v6
         sum = 1./40.*(3125./16. - 375.*v + 525./2.*v2 - 
     &        80.*v3 + 9.*v4)
         wij(i) = sum
         sum = 1./40.*(-375. + 525.*v - 240.*v2 + 36.*v3)
         grwij(i) = sum
         sum = 1./10.*(3125./48.*v3 - 375./4.*v4 + 105./2.*v5 - 
     &        40./3.*v6 + 9./7.*v7 - 2185./1344.)
         fmass(i) = sum
         sum = 1./10.*(3125./96.*v2 - 125./4.*v3 + 105./8.*v4 - 
     &        8./3.*v5 + 3./14.*v6 - 3125./128.)
         fpoten(i) = sum
         sum = 1./10.*(3125./128. - 3125./32.*v2 + 
     &        125.*v3 - 525./8.*v4 + 16.*v5 - 3./2.*v6)
         dphidh(i) = sum
      END DO
c
c--Normalisation constant
c
      cnormk = 1.0/pi
      selfnormkernel = 115./64.
      part1potenkernel = -1.0/6720.0
      part2potenkernel = 437.0/2688.0

      RETURN
      END
