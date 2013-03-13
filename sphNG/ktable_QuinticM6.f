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
      radkernel = 3.0
      part1kernel = 1.0
      part2kernel = 2.0
      v2max = radkernel*radkernel
      dvtable = v2max/itable
      ddvtable = itable/v2max
      i1 = part1kernel*part1kernel*ddvtable
      i2 = part2kernel*part2kernel*ddvtable
c
c--Build tables
c
c  a) v less than 1.0
c
      DO i = 0, i1
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         v7 = v*v6
         sum = (3.0-v)**5 - 6.0*(2.0-v)**5 + 15.0*(1.0-v)**5
         wij(i) = sum
         sum = - 5.0*(3.0-v)**4 + 30.0*(2.0-v)**4 - 75.0*(1.0-v)**4 
         grwij(i) = sum
         sum = 1.3333333333*v3 - 1.2*v5 + 0.5*v6
         fmass(i) = sum
         sum = 0.66666666666*v2 - 0.3*v4 + 0.1*v5 - 1.4
         fpoten(i) = sum
         sum = -1.4 + 2.*v2 - 1.5*v4 + 0.6*v5
         dphidh(i) = sum
      END DO
c
c  b) v greater than 1.0
c
      DO i = i1 + 1, i2
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         v7 = v*v6
         sum = (3.0-v)**5 - 6.0*(2.0-v)**5
         wij(i) = sum
         sum = - 5.0*(3.0-v)**4 + 30.0*(2.0-v)**4
         grwij(i) = sum
         sum = -0.16666666666*v6 + 1.2*v5 - 3.*v4 + 2.66666666666*v3 -
     &         0.0666666666666
         fmass(i) = sum
         sum = -0.033333333333*v5 + 0.3*v4 - v3 + 1.3333333333*v2 - 1.6
         fpoten(i) = sum
         sum = -1.6 + 4.*v2 - 4.*v3 + 1.5*v4 - 0.2*v5
         dphidh(i) = sum
      END DO
c
c  c) v greater than 2.0
c
      DO i = i2 + 1, itable
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         v7 = v*v6
         sum = (3.0-v)**5
         wij(i) = sum
         sum = - 5.0*(3.0-v)**4
         grwij(i) = sum
         sum = 1.
         fmass(i) = sum
         sum = 1.
         fpoten(i) = sum
         sum = 0.
         dphidh(i) = sum
      END DO
c
c--Normalisation constant
c
      cnormk = 1.0/(120.0*pi)
      selfnormkernel = 66.0
      part1potenkernel = 1.0/15.0
      part2potenkernel = -1.0
c
c--For dust/gas drag, need double humped kernel
c
      doublehumpnormk = 1./(168.*pi)

      RETURN
      END
