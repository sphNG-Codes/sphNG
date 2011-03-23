      SUBROUTINE gpti(ipart, npart, ntot, xyzmh, fxyzu)
c************************************************************
c                                                           *
c  Subroutine by MRB 2007.  Evaluates forces on a particle  *
c     due to all sink particles.                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,mmax2)
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
      INCLUDE 'COMMONS/rbnd'

      LOGICAL icalc
      CHARACTER*7 where

      DATA where/'gpti'/

      CALL getused(tgforpt1)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine gpti')

c
c--Needed for MPI code
c
      IF (ipart.GT.npart) THEN
         iparttree = ipart + ntot + 2
      ELSE
         iparttree = ipart
      ENDIF

      xi = xyzmh(1,iparttree)
      yi = xyzmh(2,iparttree)
      zi = xyzmh(3,iparttree)
      pmassi = xyzmh(4,iparttree)

      DO iptn = 1, nptmass
         ipt = listpm(iptn)
c
c--Don't include self interaction
c
         IF (ipt.EQ.ipart) GOTO 101

         difx = xyzmh(1,ipt) - xi
         dify = xyzmh(2,ipt) - yi
         difz = xyzmh(3,ipt) - zi
         pmassipt = xyzmh(4,ipt)
         hipt = xyzmh(5,ipt)

         rr = difx**2 + dify**2 + difz**2 + tiny
c
c--Add forces
c
c--The force definition:
c
         IF(iphase(ipart).EQ.0 .OR. iphase(ipart).GE.10
     &        .OR. iptsoft.EQ.0) THEN 
            IF (iphase(ipt).NE.5) THEN   
               rr05 = SQRT(rr)
               fff = pmassipt/(rr*rr05)
               potn = pmassipt/rr05
            ELSE
               rr05 = SQRT(rr)
               rsurface = rplanet*pradfac
               IF (rr05.LE.(2.*rsurface)) THEN
                  fsurface = (((2.*rsurface)-rr05)/
     &                 (rsurface))**4
               ELSE
                  fsurface = 0.0
               ENDIF
               fff = (1.0-fsurface)*pmassipt/(rr*rr05)
               potn = pmassipt/rr05
            ENDIF
         ELSEIF (iphase(ipart).GE.1 .AND. iphase(ipart).LT.10) THEN
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
                     fff = pmassipt/(rr*rr05)*(4.0*third*v2*v1-
     &                    1.2*v4*v1+0.5*v4*v2)
                     potn = -(pmassipt/hptmass*(2.0*third*v2-0.3*v4+
     &                    0.1*v4*v1-1.4))
                  ELSE
                     v4 = v2*v2
                     fff = pmassipt/(rr*rr05)*(8.0*third*v2*v1-
     &                    3.0*v4+1.2*v4*v1-0.5*third*v4*v2-0.2*third)
                     potn = -(pmassipt/hptmass*(4.0*third*v2-v2*v1+
     &                    0.3*v4-0.1*third*v4*v1-1.6) +
     &                    0.2*third*pmassi/rr05)
                  ENDIF
               ELSE
                  rr05 = SQRT(rr)
                  fff = pmassipt/(rr*rr05)
                  potn = pmassipt/rr05
               ENDIF
c
c--Softened potential 2
c
            ELSEIF(iptsoft.EQ.2) THEN 
               rr = rr + ptsoft*ptsoft
               rrs05 = SQRT(rr)
               rr32 = rr*rrs05
               potn = pmassipt / rrs05
               fff = pmassipt / rr32
c
c--Softened potential 3
c
            ELSEIF (iptsoft.EQ.3) THEN
               rr4 = rr*rr + ptsoft*ptsoft*ptsoft*ptsoft
               rr4s025 = (rr4)**0.25
               rr54 = rr4*rr4s025
               potn = pmassipt / rr4s025
               fff = rr*pmassipt / rr54
            ELSE
c
c--ERROR! unknown softening option
c           
               print*,' ERROR: unknown softening option in gpti'
               fff = 0.
               potn = 0.
            ENDIF
         ENDIF
c
c--End force definition
c
         fxyzu(1,ipart) = fxyzu(1,ipart) + fff*difx
         fxyzu(2,ipart) = fxyzu(2,ipart) + fff*dify
         fxyzu(3,ipart) = fxyzu(3,ipart) + fff*difz
         fxyzu(4,ipart) = fxyzu(4,ipart) - potn
 101     CONTINUE
      END DO

      CALL getused(tgforpt2)
      tgforpt = tgforpt + (tgforpt2 - tgforpt1)

      RETURN
      END
