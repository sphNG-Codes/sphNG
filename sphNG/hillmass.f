      FUNCTION hillmass(planetmass)

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/sort'

      EQUIVALENCE (ilist, next1)
      EQUIVALENCE (rr, tempsort)

      REAL rr(idim), hillr, hillm, hillmass, planetmass
      INTEGER ilist(idim)

      hillm = planetmass
      hillr = (planetmass/3.0)**(1.0/3.0)

      DO i = 1, npart
         rr(i) = SQRT((xyzmh(1,i)-1.0)**2 + xyzmh(2,i)**2 +
     &        xyzmh(3,i)**2)
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

      END FUNCTION hillmass
