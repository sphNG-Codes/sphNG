      SUBROUTINE angmom
c************************************************************
c                                                           *
c  Computes the total angular momentum of the system        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/angm'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/phase'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine angmom')
c
c--Initialisation
c
      angx = 0.
      angy = 0.
      angz = 0.
c
c--Compute total angular momentum
c
      DO i = 1, npart
         IF (iphase(i).GE.0) THEN
            angx = angx + xyzmh(4,i)*(xyzmh(2,i)*vxyzu(3,i) - 
     &           vxyzu(2,i)*xyzmh(3,i))
            angy = angy + xyzmh(4,i)*(vxyzu(1,i)*xyzmh(3,i) - 
     &           xyzmh(1,i)*vxyzu(3,i))
            angz = angz + xyzmh(4,i)*(xyzmh(1,i)*vxyzu(2,i) - 
     &           vxyzu(1,i)*xyzmh(2,i))
         ENDIF
      END DO
c
c--Add spin angular momentum of point masses
c
      DO i = 1, nptmass
         angx = angx + spinx(i)
         angy = angy + spiny(i)
         angz = angz + spinz(i)
      END DO

      angto = SQRT(angx**2 + angy**2 + angz**2)

      IF (idebug(1:6).EQ.'angmom') THEN
         WRITE (iprint, 99002) angx, angy, angz, angto
99002    FORMAT (1X, 4(1PE12.5,1X))
      ENDIF

      RETURN
      END
