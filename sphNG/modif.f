      SUBROUTINE modif
c************************************************************
c                                                           *
c  This subroutine allows modifications to be made during   *
c     the transfer of dumps.                                *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/trans'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine modif')
c
c--Find particles to change
c
      ncount = 0
      step = 0.005
      rsup = 0.000
      need = npart*frac + 1
c
c--Find radius containing particle fraction
c
 100  DO i = 1, npart
         d = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2)
         IF (d.LT.rsup) THEN
            IF (ncount + 1.GT.need) THEN
               rsup = rsup - step
               step = step/1.5
               IF (step.GE.1.E-4) GOTO 300
               GOTO 400
            ENDIF
            ncount = ncount + 1
         ENDIF
      END DO
      IF (ncount.EQ.need) GOTO 400
 300  rsup = rsup + step
      ncount = 0
      GOTO 100
c
c--Make the change
c
 400  rsup2 = rsup**2
      WRITE (*, *) 'rsup=', rsup
      totm = 0.
      totmin = 0.
      DO 500 i = 1, npart
         d2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
         IF (ichang.EQ.1) THEN
            IF (d2.LE.rsup2) THEN
               vxyzu(4,i) = vxyzu(4,i) + energc
               totmin = totmin + xyzmh(4,i)
            ENDIF
            totm = totm + xyzmh(4,i)
            GOTO 500
         ENDIF
         IF (ichang.EQ.2) THEN
            IF (d2.LE.rsup2) THEN
               vxyzu(4,i) = vxyzu(4,i) + energc
               totmin = totmin + xyzmh(4,i)
            ENDIF
            totm = totm + xyzmh(4,i)
            rx = SQRT(xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2)
            vrad = vexpan*( - 13.0*rx**10 - 7.0*rx**5)/20.0
            vxyzu(1,i) = vrad*xyzmh(1,i)/rnorm - omeg0*xyzmh(2,i)
            vxyzu(2,i) = vrad*xyzmh(2,i)/rnorm + omeg0*xyzmh(1,i)
            vxyzu(3,i) = vrad*xyzmh(3,i)/rnorm
         ENDIF
 500  CONTINUE
      WRITE (*, *) 'mass ratio ', totmin/totm

      RETURN
      END
