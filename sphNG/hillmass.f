      FUNCTION hillmass(planetmass, sswitch)

c--Added sswitch such that this calculation now works for sink
c  particles as well as fixed potentials.

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/ptmass'

      EQUIVALENCE (ilist, next1)
      EQUIVALENCE (rr, tempsort)

      REAL rr(idim), hillr, hillm, hillmass, planetmass
      INTEGER ilist(idim), sswitch

      hillm = planetmass
      hillr = (planetmass/3.0)**(1.0/3.0)

      DO i = 1, npart
         IF (sswitch.EQ.0) THEN
            rr(i) = SQRT((xyzmh(1,i)-1.0)**2 + xyzmh(2,i)**2 +
     &           xyzmh(3,i)**2)
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
            hillr = (hillm/3.0)**(1.0/3.0)
         ENDIF
      END DO

      hillmass = hillm
c      hillmass = 0.0001

      END FUNCTION hillmass
