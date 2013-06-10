      FUNCTION hillmass(planetmass, sswitch, xmass, isink)

c--Added sswitch such that this calculation now works for sink (1)
c  particles as well as fixed potentials (0).

c--Added isink, so it no longer assumes listpm(1), allowing multiple
c  sink particles to be used.

      INCLUDE 'idim'

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/pxpy'

c      EQUIVALENCE (ilist, next1)
      EQUIVALENCE (rr, tempsort)

      REAL rr(idim), hillr, hillm, hillmass, planetmass, xmass
      REAL rorbit, xsub, ysub, zsub
      INTEGER ilist(idim), sswitch, isink

      hillm = planetmass

      IF (irotpot.EQ.1 .AND. (px.EQ.0.0 .AND. py.EQ.0.0)) THEN
         print *, 'px and py = 0 in externf.'
         STOP
      ENDIF

      IF (sswitch.EQ.0) THEN
         IF (irotpot.EQ.0) THEN
            rorbit = 1.0
            xsub = 1.0
            ysub = 0.0
            zsub = 0.0
         ELSE
            rorbit = sqrt(px**2 + py**2 + pz**2)
            xsub = px
            ysub = py
            zsub = pz
         ENDIF
      ELSEIF (sswitch.EQ.1) THEN
         IF (isink.EQ.0) THEN
            print *, 'ERROR: hillmass, isink = 0'
            STOP
         ENDIF
         rorbit = sqrt(xyzmh(1,listpm(isink))**2 +
     &        xyzmh(2,listpm(isink))**2 +
     &        xyzmh(3,listpm(isink))**2)
         xsub = xyzmh(1,listpm(isink))
         ysub = xyzmh(2,listpm(isink))
         zsub = xyzmh(3,listpm(isink))
      ENDIF

      DO i = 1, npart
         rr(i) = SQRT((xyzmh(1,i) - xsub)**2 +
     &                (xyzmh(2,i) - ysub)**2 +
     &                (xyzmh(3,i) - zsub)**2)
      END DO

      CALL indexx2(npart, rr, ilist)

      hillr = rorbit*(hillm/(3.*xmass))**(1.0/3.0)

      DO j = 1, npart
         i = ilist(j)
         IF (iphase(i).EQ.0) THEN
            IF(rr(i).LE.hillr) hillm = hillm + xyzmh(4,i)
            hillr = (hillm/(3.*xmass))**(1.0/3.0)
         ENDIF
      END DO

      hillmass = hillm

      END FUNCTION hillmass
