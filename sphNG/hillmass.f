      FUNCTION hillmass(planetmass, sswitch, xmass)

c--Added sswitch such that this calculation now works for sink
c  particles as well as fixed potentials.

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/pxpy'

      EQUIVALENCE (ilist, next1)
      EQUIVALENCE (rr, tempsort)

      REAL rr(idim), hillr, hillm, hillmass, planetmass, xmass
      INTEGER ilist(idim), sswitch

      hillm = planetmass
      hillr = (hillm/(3.*xmass))**(1.0/3.0)

      IF (irotpot.EQ.1 .AND. (px.EQ.0.0 .AND. py.EQ.0.0)) THEN
         print *, 'px and py = 0 in externf.'
         STOP
      ENDIF

      DO i = 1, npart
         IF (sswitch.EQ.0) THEN
            IF (irotpot.EQ.0) THEN
               rr(i) = SQRT((xyzmh(1,i)-1.0)**2 + xyzmh(2,i)**2 +
     &              xyzmh(3,i)**2)
            ELSE
               rr(i) = SQRT((xyzmh(1,i)-px)**2 + (xyzmh(2,i)-py)**2 +
     &              xyzmh(3,i)**2)
            ENDIF
         ELSEIF (sswitch.EQ.1) THEN
            rr(i) = SQRT((xyzmh(1,i)-xyzmh(1,listpm(1)))**2 +
     &           (xyzmh(2,i) - xyzmh(2,listpm(1)))**2 +
     &           (xyzmh(3,i) - xyzmh(3,listpm(1)))**2)
         ENDIF
      END DO

      CALL indexx2(npart, rr, ilist)

      DO j = 1, npart
         i = ilist(j)
         IF (iphase(i).EQ.0) THEN
            IF(rr(i).LE.hillr) hillm = hillm + xyzmh(4,i)
            hillr = (hillm/(3.*xmass))**(1.0/3.0)
         ENDIF
      END DO

      hillmass = hillm

      END FUNCTION hillmass
