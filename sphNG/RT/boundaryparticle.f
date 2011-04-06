      FUNCTION boundaryparticle(i,xyzmh,rhoi)

      INCLUDE 'idim'
      INCLUDE 'igrape'

      LOGICAL boundaryparticle
      INTEGER i
      DIMENSION xyzmh(5,mmax2)
c
c--NOTE: rhoi here is REAL not REAL*4
c
      REAL rhoi

      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/boundheight'

      REAL mu
c
c--If there is a second object, this object acts as the boundary
c
      IF (n2.GT.0) THEN
         boundaryparticle = iunique(iorig(i)).GT.n1
c
c--Disc boundary
c
      ELSEIF (ibound.EQ.100) THEN
         boundtest = ABS(xyzmh(3,i)/(hoverr*SQRT(
     &        xyzmh(1,i)**2 + xyzmh(2,i)**2)))
         boundaryparticle = boundtest.GT.bounddens
c
c--Disc + planet boundary
c
      ELSEIF (ibound/10.EQ.10) THEN
         radu = sqrt(xyzmh(1,i)**2 + xyzmh(2,i)**2)
         izoverh = INT((radu - rmind)/deltar + 1)
         IF (izoverh.LT.1) izoverh = 1
         IF (izoverh.GT.1000) izoverh = 1000
         mu = (radu - (rmind + (izoverh-1.0)*deltar))/deltar
         bounddenslocal = (zoverh(izoverh+1) -
     &        zoverh(izoverh))*mu + zoverh(izoverh)
         IF (use_tprof) THEN
            boundtest = ABS(xyzmh(3,i)/(hoverr*radu**
     &           (0.5*(tprof+3))))
         ELSE
            boundtest = ABS(xyzmh(3,i)/(hoverr*radu))
         ENDIF
         boundaryparticle = boundtest.GT.bounddenslocal
         IF (radu.LT.0.385 .AND. use_tprof)
     &        boundaryparticle = .TRUE.
c
c--Low density gas boundary (e.g. for stellar cluster formation)
c
      ELSE
         boundaryparticle = rhoi.LT.bounddens
      ENDIF

      RETURN
      END
