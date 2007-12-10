      SUBROUTINE gptfromgas(ipart, npart, ntot, xyzmh, fxyzu)
c************************************************************
c                                                           *
c  Subroutine by IAB, MRB 1994.  Evaluates forces on point  *
c     mass due to all other particles.                      *
c     Returns list of nearest neighbours to point mass !    *
c     THIS ROUTINE VECTORIZABLE. (?)                        *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,mmax)
      DIMENSION fxyzu(4,idim)

      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/gravi'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nearmpt'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/mpidebug'

      CHARACTER*7 where

      DATA where/'gptfrom'/

      CALL getused(tgforpt1)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine gptfromgas')

      IF (iphase(ipart).LE.0) THEN
         WRITE (*,*) 'ERROR - gptfromgas called by non-sink'
         CALL quit
      ENDIF
c
c--Needed for MPI code
c
      IF (ipart.GT.npart) THEN
         iparttree = ipart + ntot + 2
      ELSE
         iparttree = ipart
      ENDIF

      rrx = xyzmh(1,iparttree)
      rry = xyzmh(2,iparttree)
      rrz = xyzmh(3,iparttree)
      DO 101 j = 1, npart
         IF(iphase(j).NE.0) GOTO 101

         pmassj = xyzmh(4,j)
     
         difx = xyzmh(1,j) - rrx
         dify = xyzmh(2,j) - rry
         difz = xyzmh(3,j) - rrz     

         rr = difx**2 + dify**2 + difz**2 + tiny
c
c--Add forces
c
c--The force definition:
c
         IF(iphase(j).EQ.0 .OR. iptsoft.EQ.0) THEN 
            rr05 = SQRT(rr)
            fff = pmassj/(rr*rr05)
            potn = pmassj/rr05
         ELSEIF (iphase(j).GE.1) THEN
c
c--Non-pointmass force
c
c
c--Softened potential 1 (spline kernel softening)
c
            IF (iptsoft.EQ.1) THEN
               hptmass = 0.5*ptsoft
               v2 = rr/hptmass**2
               IF (v2.LT.4.0) THEN
                  v1 = SQRT(v2)
                  rr05 = v1*hptmass
                  third = 1.0/3.0
                  IF (v2.LT.1.0) THEN
                     v4 = v2*v2
                     fff = pmassj/(rr*rr05)*(4.0*third*v2*v1-
     &                    1.2*v4*v1+0.5*v4*v2)
                     potn = -(pmassj/hptmass*(2.0*third*v2-0.3*v4+
     &                    0.1*v4*v1-1.4))
                  ELSE
                     v4 = v2*v2
                     fff = pmassj/(rr*rr05)*(8.0*third*v2*v1-
     &                    3.0*v4+1.2*v4*v1-0.5*third*v4*v2-0.2*third)
                     potn = -(pmassj/hptmass*(4.0*third*v2-v2*v1+
     &                    0.3*v4-0.1*third*v4*v1-1.6) +
     &                       0.2*third*pmassj/rr05)
                  ENDIF
               ELSE
                  rr05 = SQRT(rr)
                  fff = pmassj/(rr*rr05)
                  potn = pmassj/rr05
               ENDIF
c     
c--Softened potential 2
c
            ELSEIF(iptsoft.EQ.2) THEN 
               rr = rr + ptsoft*ptsoft
               rrs05 = SQRT(rr)
               rr32 = rr*rrs05
               potn = pmassj / rrs05
               fff = pmassj / rr32
c
c--Softened potential 3
c
            ELSEIF (iptsoft.EQ.3) THEN
               rr4 = rr*rr + ptsoft*ptsoft*ptsoft*ptsoft
               rr4s025 = (rr4)**0.25
               rr54 = rr4*rr4s025
               potn = pmassj / rr4s025
               fff = rr*pmassj / rr54
            ENDIF
         ENDIF
c
c--End force definition
c
         fxyzu(1,ipart) = fxyzu(1,ipart) + fff*difx
         fxyzu(2,ipart) = fxyzu(2,ipart) + fff*dify
         fxyzu(3,ipart) = fxyzu(3,ipart) + fff*difz
         fxyzu(4,ipart) = fxyzu(4,ipart) - potn
 101  CONTINUE

      CALL getused(tgforpt2)
      tgforpt = tgforpt + (tgforpt2 - tgforpt1)

      RETURN
      END
