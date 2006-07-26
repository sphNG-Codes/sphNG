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

      REAL*8 sum, v, v2, v3, v4, v5, v6

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
      radkernel = 2.0
      part1kernel = 1.0
      part2kernel = 2.0
      v2max = radkernel*radkernel
      dvtable = v2max/itable
      ddvtable = itable/v2max
      i1 = part1kernel*ddvtable
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
ccc         sum = 1. - 1.5*v2 + 0.75*v3
         sum = 3.0/4.0*(10./3. - 7.*v2 + 4.*v3)
         wij(i) = sum
ccc         sum = -3.*v + 2.25*v2
         sum = 3.0/4.0*(-14.*v + 12.*v2)
         grwij(i) = sum
ccc         sum = 1.3333333333*v3 - 1.2*v5 + 0.5*v6
         sum = 10./3.*v3 - 4.2*v5 + 2.*v6
         fmass(i) = sum
ccc         sum = 0.66666666666*v2 - 0.3*v4 + 0.1*v5 - 1.4
         sum = 5./3.*v2 - 21./20.*v4 + 0.4*v5 - 2.1
         fpoten(i) = sum
ccc         sum = -1.4 + 2.*v2 - 1.5*v4 + 0.6*v5
         sum = -2.1 + 5.*v2 - 5.25*v4 + 2.4*v5
         dphidh(i) = sum
      END DO
c
c  b) v greater than 1
c
      DO i = i1 + 1, itable
         v2 = i*dvtable
         v = DSQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         dif2 = 2. - v
ccc         sum = 0.25*dif2*dif2*dif2
         sum = 3.0/4.0*(dif2*dif2*(5. - 4.*v)/3.)
         wij(i) = sum
ccc         sum = -0.75*v2 + 3*v - 3.
         sum = 3.0/4.0*(-4.*v2 + 14.*v - 12.)
         grwij(i) = sum
ccc         sum = -0.16666666666*v6 + 1.2*v5 - 3.*v4 + 2.66666666666*v3 -
ccc     &         0.0666666666666
         sum = -2./3.*v6 + 4.2*v5 - 9.*v4 + 20./3.*v3 -
     &         0.0666666666666
         fmass(i) = sum
ccc         sum = -0.033333333333*v5 + 0.3*v4 - v3 + 1.3333333333*v2 - 1.6
         sum = -2./15.*v5 + 21./20.*v4 - 3.*v3 + 10./3.*v2 - 2.4
         fpoten(i) = sum
ccc         sum = -1.6 + 4.*v2 - 4.*v3 + 1.5*v4 - 0.2*v5
         sum = -2.4 + 10.*v2 - 12.*v3 + 5.25*v4 - 12./15.*v5
         dphidh(i) = sum
      END DO
c
c--Normalisation constant
c
      cnormk = 1.0/pi
      selfnormkernel = 10.0/4.0
      part1potenkernel = 1.0/15.0
      part2potenkernel = 0.0

      RETURN
      END
