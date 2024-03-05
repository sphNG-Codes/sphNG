      SUBROUTINE quadratic(a1,a2,a3,old,soln,moresweep,ipartin)
c************************************************************
c                                                           *
c  Solves a quadratic equation.                             *
c    Used for implicit one-fluid dust (Elsender & Bate 2024)*
c                                                           *
c     Code written by DE (2024)                             *
c                                                           *
c************************************************************

      IMPLICIT NONE

      REAL tiny
      PARAMETER (tiny = 1.0E-30)

      REAL a1,a2,a3,old,soln
      REAL dis,z1,z2
      REAL biggest_term

      INTEGER ipartin

      LOGICAL moresweep

      IF (a1.NE.0.0) THEN
c
c--First check the discriminant 
c
         dis = a2**2 - 4*a1*a3
         biggest_term = MAX(ABS(a2**2),ABS(-4*a1*a3))
         IF (dis.LT.0.AND.ABS(dis)/biggest_term.LT.1d-12) THEN
            dis = 0.0
         ELSE IF (dis.LT.0.0) THEN
            PRINT*, 'quadratic: discriminant is negative ',ipartin
            PRINT*, 'discriminant: ',dis
            moresweep=.TRUE.
            RETURN
         ENDIF ! dis.LT.0

         z1 = 0.5*((-a2)+SQRT(dis))/a1
         z2 = 0.5*((-a2)-SQRT(dis))/a1
      
         IF (z1.GT.0.0 .AND. z2.LE.0.0) THEN
            soln = z1
         ELSEIF (z1.LE.0.0 .AND. z2.GT.0.0) THEN
            soln = z2
         ELSEIF (z1.LE.0.0 .AND. z2.LE.0.0) THEN
            PRINT*, 'Failed 1 ',z1,z2,old,ipartin,a1,a2,a3
c         PRINT*, a1,a2,a3
c         PRINT*,'    ',dis
            CALL quit(1)
         ELSEIF (ABS(z1-old).GT.ABS(z2-old)) THEN
            soln = z2
            IF (ABS((soln-old)/old).GT.1.0) THEN
               PRINT*, 'Change big a',z1,z2,old,ipartin
               PRINT*, a1,a2,a3
               PRINT*, dis
               CALL quit(1)
            ENDIF
         ELSE
            soln = z1
            IF (ABS((soln-old)/old).GT.1.0) THEN
               PRINT*, 'Change big b',z1,z2,old,ipartin
               PRINT*, a1,a2,a3
               PRINT*, dis
               CALL quit(1)
            ENDIF
         ENDIF
      ELSE
         soln = -a3/a2
      ENDIF

      RETURN

      END
