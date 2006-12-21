      SUBROUTINE indexx3unique(n, list, arr, ipos, indx)
c************************************************************
c                                                           *
c  This is INDEXX using the quicksort algorithm.            *
c                                                           *
c************************************************************

      PARAMETER (m=7, nstack=50)

      DIMENSION list(n), istack(nstack), indx(n), arr(5,*)

      CHARACTER*7 where

      DATA where/'indexx'/

      DO j = 1, n
         indx(j) = j
      END DO
      jstack = 0
      l = 1
      ir = n

 1    IF (ir - l.LT.m) THEN
         DO j = l + 1, ir
            indxt = indx(j)
            listindxt = list(indxt)
            a = arr(ipos,listindxt)
            DO i = j - 1, 1, -1
               IF (arr(ipos,list(indx(i))).LT.a) GOTO 2
               IF (arr(ipos,list(indx(i))).EQ.a .AND.
     &              list(indx(i)).LT.listindxt) GOTO 2
               indx(i + 1) = indx(i)
            END DO
            i = 0
 2          indx(i + 1) = indxt
         END DO
         IF (jstack.EQ.0) RETURN
         ir = istack(jstack)
         l = istack(jstack - 1)
         jstack = jstack - 2
      ELSE
         k = (l + ir)/2
         itemp = indx(k)
         indx(k) = indx(l + 1)
         indx(l + 1) = itemp
         IF (arr(ipos,list(indx(l+1))).GT.arr(ipos,list(indx(ir)))) THEN
            itemp = indx(l + 1)
            indx(l + 1) = indx(ir)
            indx(ir) = itemp
         ELSEIF (arr(ipos,list(indx(l+1))).EQ.arr(ipos,list(indx(ir)))
     &           .AND. list(indx(l + 1)).GT.list(indx(ir))) THEN
            itemp = indx(l + 1)
            indx(l + 1) = indx(ir)
            indx(ir) = itemp
         ENDIF
         IF (arr(ipos,list(indx(l))).GT.arr(ipos,list(indx(ir)))) THEN
            itemp = indx(l)
            indx(l) = indx(ir)
            indx(ir) = itemp
         ELSEIF (arr(ipos,list(indx(l))).EQ.arr(ipos,list(indx(ir)))
     &           .AND. list(indx(l)).GT.list(indx(ir))) THEN
            itemp = indx(l)
            indx(l) = indx(ir)
            indx(ir) = itemp
         ENDIF
         IF (arr(ipos,list(indx(l+1))).GT.arr(ipos,list(indx(l)))) THEN
            itemp = indx(l + 1)
            indx(l + 1) = indx(l)
            indx(l) = itemp
         ELSEIF (arr(ipos,list(indx(l+1))).EQ.arr(ipos,list(indx(l)))
     &           .AND. list(indx(l + 1)).GT.list(indx(l))) THEN
            itemp = indx(l + 1)
            indx(l + 1) = indx(l)
            indx(l) = itemp
         ENDIF
         i = l + 1
         j = ir
         indxt = indx(l)
         listindxt = list(indxt)
         a = arr(ipos,listindxt)

 3       CONTINUE
         i = i + 1
         IF (arr(ipos,list(indx(i))).LT.a) GOTO 3
         IF (arr(ipos,list(indx(i))).EQ.a .AND.
     &        list(indx(i)).LT.listindxt) GOTO 3
 4       CONTINUE
         j = j - 1
         IF (arr(ipos,list(indx(j))).GT.a) GOTO 4
         IF (arr(ipos,list(indx(j))).EQ.a .AND.
     &        list(indx(j)).GT.listindxt) GOTO 4
         IF (j.LT.i) GOTO 5
         itemp = indx(i)
         indx(i) = indx(j)
         indx(j) = itemp
         GOTO 3

 5       indx(l) = indx(j)
         indx(j) = indxt
         jstack = jstack + 2
         IF (jstack.GT.nstack) CALL error(where, 1)
         IF (ir - i + 1.GE.j - l) THEN
            istack(jstack) = ir
            istack(jstack - 1) = i
            ir = j - 1
         ELSE
            istack(jstack) = j - 1
            istack(jstack - 1) = l
            l = i
         ENDIF
      ENDIF

      GOTO 1
      END
