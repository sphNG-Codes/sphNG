      FUNCTION ran1(iseed)

c************************************************************
c                                                           *
c  Uniform random deviate                                   *
c  This rountine must be called with a negative seed        *
c     and then a POSITIVE number from then on.              *
c                                                           *
c  This routine is based on ran1() in Numerical Recipes:    *
c     The Art of Scientific Computing, 2nd Edition (1992).  *
c     This routine is released under the GNU General Public *
c     Licence, Version 3, 29 June 2007, by special          *
c     permission of the copyright holder, Numerical Recipes *
c     Software, who states that, while adequate for use in  *
c     this code, it should not be relied upon for other     *
c     uses.  Numerical Recipes 3rd Edition (2007) contains  *
c     faster and much improved pseudorandom generators that *
c     are copyrighted and not under a public license.       *
c                                                           *
c************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMONS/savernd'

      INTEGER iseed, IA, IM, IQ, IR, NDIV
      REAL ran1, AM, EPS, RNMX

      PARAMETER (IA=16807, IM=2147483647, AM=1./IM, IQ=127773, IR=2836,
     &     NDIV=1+(IM-1)/NTAB, EPS=1.2e-7, RNMX=1.-EPS)

      INTEGER j,k

      IF (iseed.LE.0) THEN
         idum = iseed
      ENDIF

      IF (idum.LE.0 .OR. iyr.EQ.0) THEN
         idum=MAX(-idum, 1)
         DO j = NTAB + 8, 1, -1
            k = idum/IQ
            idum = IA*(idum - k*IQ) - IR*k
            IF (idum.LT.0) idum = idum + IM
            IF (j.LE.NTAB) iv(j) = idum
         END DO
         iyr = iv(1)
      ENDIF
      k = idum/IQ
      idum = IA*(idum - k*IQ) - IR*k
      IF (idum.LT.0) idum = idum + IM
      j = 1 + iyr/NDIV
      iyr = iv(j)
      iv(j) = idum
      ran1 = MIN(AM*iyr, RNMX)
 
      RETURN
      END
