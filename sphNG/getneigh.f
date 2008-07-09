      SUBROUTINE getneigh(ipart,ntot,hi,xyzmh,limits,nlist,
     &     nlistmax,list)

      INCLUDE 'idim'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'

      DIMENSION xyzmh(5,limits), list(iptneigh)
c
c--Not setting nlist=0 allows particles to be added onto an existing list.
c     However, it also means that nlist MUST BE SET BEFORE CALLING getneigh!
c
ccc      nlist = 0

      DO i = 1, ntot
         IF (iphase(i).EQ.0 .OR. iphase(i).GE.10) THEN
            dx = xyzmh(1,ipart) - xyzmh(1,i)
            dy = xyzmh(2,ipart) - xyzmh(2,i)
            dz = xyzmh(3,ipart) - xyzmh(3,i)
            r2 = dx**2 + dy**2 + dz**2
            IF (r2.LT.radkernel*radkernel*hi*hi .AND. i.NE.ipart) THEN
               nlist = nlist + 1
               IF (nlist.GT.nlistmax) THEN
                  WRITE (iprint,*) 'ERROR - getneigh list overflow'
                  CALL quit
               ENDIF
               list(nlist) = i
            ENDIF
         ENDIF
      END DO

      RETURN
      END
