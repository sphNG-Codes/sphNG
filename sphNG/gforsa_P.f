      SUBROUTINE gforsa(m,ntot,nlistga,listga,xyzmh,fsx,fsy,fsz,epot)
c************************************************************
c                                                           *
c  Subroutine by W. Press (11/21/86).  Evaluates force on   *
c     particle M due to a list of other particles.          *
c     THIS ROUTINE VECTORIZABLE.                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptsoft'

      DIMENSION listga(idim), xyzmh(5,mmax2)

c
c--Need for MPI code
c
      IF (m.GT.ntot) THEN
         mtot2 = m + ntot + 2
         rrx = xyzmh(1,mtot2)
         rry = xyzmh(2,mtot2)
         rrz = xyzmh(3,mtot2)
      ELSE
         rrx = xyzmh(1,m)
         rry = xyzmh(2,m)
         rrz = xyzmh(3,m)
      ENDIF

      DO 101 j = 1, nlistga
         n = listga(j)
         difx = xyzmh(1,n) - rrx
         dify = xyzmh(2,n) - rry
         difz = xyzmh(3,n) - rrz

         rr = difx**2 + dify**2 + difz**2 + tiny

         IF (n.LE.natom) THEN
            pmassn = imfac(n)*xyzmh(4,n)
         ELSE
            pmassn = xyzmh(4,n)
         ENDIF
c
c--If one of the particles is a point mass, then soften interaction via ptsoft
c
         IF (iphase(m).GE.1 .AND. iphase(m).LT.10 .OR. 
     &        iphase(n).GE.1 .AND. iphase(n).LT.10) THEN
            IF (iptsoft.EQ.0) THEN 
               rr05 = SQRT(rr)
               fff = pmassn/(rr*rr05)
               potn = pmassn/rr05
c
c--Softened potential 1 (spline kernel softening)
c
            ELSEIF (iptsoft.EQ.1) THEN
               hptmass = 0.5*ptsoft
               v2 = rr/hptmass**2
               IF (v2.LT.4.0) THEN
                  v1 = SQRT(v2)
                  rr05 = v1*hptmass
                  third = 1.0/3.0
                  IF (v2.LT.1.0) THEN
                     v4 = v2*v2
                     fff = pmassn/(rr*rr05)*(4.0*third*v2*v1-
     &                    1.2*v4*v1+0.5*v4*v2)
                     potn = -(pmassn/hptmass*(2.0*third*v2-0.3*v4+
     &                    0.1*v4*v1-1.4))
                  ELSE
                     v4 = v2*v2
                     fff = pmassn/(rr*rr05)*(8.0*third*v2*v1-
     &                    3.0*v4+1.2*v4*v1-0.5*third*v4*v2-0.2*third)
                     potn = -(pmassn/hptmass*(4.0*third*v2-v2*v1+
     &                    0.3*v4-0.1*third*v4*v1-1.6) +
     &                    0.2*third*pmassn/rr05)
                  ENDIF
               ELSE
                  rr05 = SQRT(rr)
                  fff = pmassn/(rr*rr05)
                  potn = pmassn/rr05
               ENDIF
c
c--Softened potential 2
c
            ELSEIF (iptsoft.EQ.2) THEN 
               rr = rr + ptsoft*ptsoft
               rrs05 = SQRT(rr)
               rr32 = rr*rrs05
               potn = pmassn / rrs05
               fff = pmassn / rr32
c
c--Softened potential 3
c
            ELSEIF (iptsoft.EQ.3) THEN
               rr4 = rr*rr + ptsoft*ptsoft*ptsoft*ptsoft
               rr4s025 = (rr4)**0.25
               rr54 = rr4*rr4s025
               potn = pmassn / rr4s025
               fff = rr*pmassn / rr54
            ENDIF            
c
c--Otherwise, both particles are gas and are not neighbours
c
         ELSE
            IF (isoft.EQ.1) rr = rr + psoft**2

            rr05 = SQRT(rr)
            potn = pmassn/rr05
            fff = potn/rr
         ENDIF

         epot = epot + potn
         fsx = fsx + fff*difx
         fsy = fsy + fff*dify
         fsz = fsz + fff*difz
 101  CONTINUE

      RETURN
      END
