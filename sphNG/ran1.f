      FUNCTION ran1(iseed)

c************************************************************
c                                                           *
c  Uniform random deviate (from Numerical Recipes 2)        *
c  This rountine must be called with a negative seed        *
c     and then a POSITIVE number from then on               *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMONS/savernd'

      INTEGER iseed, IA, IM, IQ, IR, NTAB, NDIV
      REAL ran1, AM, EPS, RNMX

      PARAMETER (IA=16807, IM=2147483647, AM=1./IM, IQ=127773, IR=2836,
     &     NTAB=32, NDIV=1+(IM-1)/NTAB, EPS=1.2e-7, RNMX=1.-EPS)

      INTEGER j,k,iv(NTAB),iy
      SAVE iv, iy
      DATA iv /NTAB*0/, iy /0/

      IF (iseed.LE.0) THEN
         idum = iseed
      ENDIF

      IF (idum.LE.0 .OR. iy.EQ.0) THEN
         idum=MAX(-idum, 1)
         DO j = NTAB + 8, 1, -1
            k = idum/IQ
            idum = IA*(idum - k*IQ) - IR*k
            IF (idum.LT.0) idum = idum + IM
            IF (j.LE.NTAB) iv(j) = idum
         END DO
         iy = iv(1)
      ENDIF
      k = idum/IQ
      idum = IA*(idum - k*IQ) - IR*k
      IF (idum.LT.0) idum = idum + IM
      j = 1 + iy/NDIV
      iy = iv(j)
      iv(j) = idum
      ran1 = MIN(AM*iy, RNMX)
 
      RETURN
      END
