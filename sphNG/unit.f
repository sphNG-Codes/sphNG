      SUBROUTINE unit
c************************************************************
c                                                           *
c  This routine computes the transformation between the     *
c     physical units (cgs) and the units used in the code.  *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/vargam'
      INCLUDE 'COMMONS/physeos'

      REAL*8 uerg
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine units')
c
c--Specify mass unit (g)
c
      fnbtot = 1.
c      fnbtot = 1000.
c      fnbtot = 0.1
c      fnbtot = 0.2
      umass = DBLE(fnbtot)*solarm
c
c--Specify distance unit (cm)
c
      udist = 1.e16
c      udist = 1.e15
c      udist = 0.1*pc
c      udist = pc
c      udist = 1.e14
c      udist = 1.e13
c      udist = 1000.0*pc
c
c--------------------
c--VARIABLE GAMMA EOS
c--------------------
c
c--Critical density for changing gamma from 1 to 1.4 for variable eq. of state
c
c      rhocrit = 1.e-14
c      rhocrit = 2.e-13
      rhocrit = 1.e-13
c      rhocrit = 1.e-12
c      rhocrit = 1.e-16
      gam = 1.4
c      gam = 5./3.
c
c--Critical density for changing gamma from 1.4 to 1.1 for variable e.o.s.
c     i.e. at 2000K assuming trans to 1.4 at 10K
c
c-- gam = 1.4
c
c      rhocrit2 = rhocrit * (200.**2.5)
c
c-- gam = 5/3
c
c      rhocrit2 = rhocrit * (200.**1.5)
c
      rhocrit2 = 1.0e-10
c      rhocrit2 = 1.0e-11
c      rhocrit2 = 1.0e-12
c      gamdh = 1.10
c      gamdh = 1.15
c      gamdh = 1.05
      gamdh = 1.0
c
c--Critical density for changing gamma from 1.1 to 5/3 for variable e.o.s.
c
c      rhocrit3 = 1.0e-0
      rhocrit3 = 1.0e-3
c      rhocrit3 = 1.0e-10
c      rhocrit3 = 1.0e-11
c      rhocrit3 = 1.0e-12
      gamah = 5./3.
c
c--****** For turning off 2nd collapse phase! ******
c
c      rhocrit2 = rhocrit3
c
c
c
c-------------------------------------------
c--PHYSICAL EOS (Bodenheimer, Bate, Burkert)
c-------------------------------------------
c
c--Changing gamma from 1 to 5/3 for physical eq. of state
c
      gamphys1 = 5./3.
      rhophys1 = 2.0e-13
      rhoref1 = rhophys1
      rhochange1 = rhophys1*27.
c
c--Changing gamma from 5/3 to 1.4 for physical e.o.s.
c
      gamphys2 = 1.4
      rhoref2 = rhophys1/9.
      rhochange2 = 7.5e-11
c
c--Changing gamma from 1.4 to 2 for physical e.o.s.
c
      gamphys3 = 2.0
      rhoref3 = 2.9e-12
c
c
c
c
c
c--Transformation factor for :
c
c  a) density
c
      udens = DBLE(umass)/DBLE(udist)**3

      rhocrit = rhocrit / udens
      rhocrit2 = rhocrit2 / udens
      rhocrit3 = rhocrit3 / udens

      rhophys1 = rhophys1 / udens
      rhoref1 = rhoref1 / udens
      rhoref2 = rhoref2 / udens
      rhoref3 = rhoref3 / udens
      rhochange1 = rhochange1 / udens
      rhochange2 = rhochange2 / udens
c
c  b) time
c
      utime = DSQRT(DBLE(udist)**3/(DBLE(gg)*DBLE(umass)))
c
c  c) ergs
c
      uerg = DBLE(umass)*DBLE(udist)**2/DBLE(utime)**2
c
c  c) ergs per gram
c
      uergg = DBLE(udist)**2/DBLE(utime)**2
c
c  d) ergs per cc
c
      uergcc = DBLE(umass)/(DBLE(udist)*DBLE(utime)**2)
c
c  e) magnetic flux density
c
c     (specify charge unit in esu)
c
      ucharge = DSQRT(DBLE(umass)*DBLE(udist)/cgsmu0)
c
c     (set units for magnetic field)
c
      umagfd = DBLE(umass)/(DBLE(utime)*DBLE(ucharge))

      IF (idebug(1:4).EQ.'unit') THEN
         WRITE (iprint, 99002) umass, udist, udens, utime, uergg, uergcc
99002    FORMAT (1X, 5(1PE12.5,1X), /, 1X, 2(1PE12.5,1X))
      ENDIF

      RETURN
      END
