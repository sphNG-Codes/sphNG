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

      WRITE (iprint, 99002)
99002 FORMAT (' KTABLE: M4 cubic ')

c
c--Maximum interaction length and step size
c
      radkernel = 2.0
      part1kernel = 1.0
      part2kernel = 2.0
      v2max = radkernel*radkernel
      dvtable = v2max/itable
      ddvtable = itable/v2max
      i1 = part1kernel*part1kernel*ddvtable
c
c--Build tables
c
c  a) v less than 1
c
      DO i = 0, i1
         v2 = i*dvtable
         v = SQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         sum = 1. - 1.5*v2 + 0.75*v3
         wij(i) = sum
         sum = -3.*v + 2.25*v2
         grwij(i) = sum
         sum = 1.3333333333*v3 - 1.2*v5 + 0.5*v6
         fmass(i) = sum
         sum = 0.66666666666*v2 - 0.3*v4 + 0.1*v5 - 1.4
         fpoten(i) = sum
         sum = -1.4 + 2.*v2 - 1.5*v4 + 0.6*v5
         dphidh(i) = sum
      END DO
c
c  b) v greater than 1
c
      DO i = i1 + 1, itable
         v2 = i*dvtable
         v = SQRT(v2)
         v3 = v*v2
         v4 = v*v3
         v5 = v*v4
         v6 = v*v5
         dif2 = 2. - v
         sum = 0.25*dif2*dif2*dif2
         wij(i) = sum
         sum = -0.75*v2 + 3*v - 3.
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
c--Normalisation constant
c
      cnormk = 1.0/pi
      selfnormkernel = 1.0
      part1potenkernel = 1.0/15.0
      part2potenkernel = 0.0
c
c--For dust/gas drag, need double humped kernel
c
      doublehumpnormk = 10./(9.*pi)

      RETURN
      END
