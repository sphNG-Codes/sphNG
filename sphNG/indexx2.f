      SUBROUTINE indexx2(n, arrin, indx)
c************************************************************
c                                                           *
c  Subroutine for sorting see W. Press.                     *
c                                                           *
c************************************************************

      DIMENSION arrin(n), indx(n)

      IF (n .EQ. 1) THEN
         indx(n) = n
      ELSEIF (n .GT. 1) THEN

         DO i = 1, n
            indx(i) = i
         END DO

         l = n/2 + 1
         ir = n
  200    IF (l.GT.1) THEN
            l = l - 1
            indxt = indx(l)
            q = arrin(indxt)
         ELSE
            indxt = indx(ir)
            q = arrin(indxt)
            indx(ir) = indx(1)
            ir = ir - 1
            IF (ir.EQ.1) THEN
               indx(1) = indxt
               RETURN
            ENDIF
         ENDIF

         i = l
         j = l + l
  300    IF (j.LE.ir) THEN
            IF (j.LT.ir) THEN
               IF (arrin(indx(j)).LT.arrin(indx(j+1))) j = j + 1
            ENDIF
            IF (q.LT.arrin(indx(j))) THEN
               indx(i) = indx(j)
               i = j
               j = j + j
            ELSE
               j = ir + 1
            ENDIF
            GOTO 300
         ENDIF
         indx(i) = indxt

         GOTO 200
      ENDIF
      END
