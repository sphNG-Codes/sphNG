      SUBROUTINE gforspt(nlst, llist, npart, ntot, xyzmh, fxyzu)
c************************************************************
c                                                           *
c  Subroutine by MRB 2000.  Evaluates forces on point       *
c     masses due ONLY to other point masses (NOT GAS).      *
c     ALSO evaluates forces on current gas particles from   *
c     point masses (PT-GAS).  The GAS-GAS and GAS-PT are    *
c     done in the TREE.                                     *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,mmax)
      DIMENSION fxyzu(4,idim)
      DIMENSION llist(idim)

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

      DATA where/'gforspt'/

      CALL getused(tgforpt1)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine gforspt')

C$OMP PARALLEL default(none)
C$OMP& shared(nptmass,listpm,fxyzu,poten,xyzmh,npart,ntot)
C$OMP& shared(pmass,iptsoft,ptsoft,nlst,llist,iphase)
C$OMP& private(i,ipart,iparttree,hipt,rrx,rry,rrz)
C$OMP& private(jptn,jpt,pmassj,difx,dify,difz,rr,rr05,fff,potn)
C$OMP& private(rrs05,rr32,rr4,rr4s025,rr54)
C$OMP& private(hptmass,v2,v1,third,v4)
C$OMP DO SCHEDULE(runtime)
      DO i = 1, nlst
         ipart = llist(i)
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
         hipt = xyzmh(5,iparttree)

         DO jptn = 1, nptmass
            jpt = listpm(jptn)
            IF (jpt.NE.ipart) THEN
               pmassj = xyzmh(4,jpt)

               difx = xyzmh(1,jpt) - rrx
               dify = xyzmh(2,jpt) - rry
               difz = xyzmh(3,jpt) - rrz     

               rr = difx**2 + dify**2 + difz**2 + tiny
c
c--Add forces
c
c--NOTE: Interactions between two sinks and force on gas due to sink are
c     softened using iptsoft and ptsoft now (28/09/2005)
c
c--The force definition:
c
c               IF (iphase(ipart).EQ.0 .OR. iptsoft.EQ.0) THEN 
               IF (iptsoft.EQ.0) THEN 
                  rr05 = SQRT(rr)
                  fff = pmassj/(rr*rr05)
                  potn = pmassj/rr05
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
                        fff = pmassj/(rr*rr05)*(4.0*third*v2*v1-
     &                                           1.2*v4*v1+0.5*v4*v2)
                        potn = -(pmassj/hptmass*(2.0*third*v2-0.3*v4+
     &                                                0.1*v4*v1-1.4))
                     ELSE
                        v4 = v2*v2
                        fff = pmassj/(rr*rr05)*(8.0*third*v2*v1-
     &                    3.0*v4+1.2*v4*v1-0.5*third*v4*v2-0.2*third)
                        potn = -(pmassj/hptmass*(4.0*third*v2-v2*v1+
     &                       0.3*v4-0.1*third*v4*v1-1.6) +
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
               ELSEIF (iptsoft.EQ.2) THEN 
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
c
c--End force definition
c
               fxyzu(1,ipart) = fxyzu(1,ipart) + fff*difx
               fxyzu(2,ipart) = fxyzu(2,ipart) + fff*dify
               fxyzu(3,ipart) = fxyzu(3,ipart) + fff*difz
               fxyzu(4,ipart) = fxyzu(4,ipart) - potn
               poten(ipart) = poten(ipart) - potn
            ENDIF
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL

      CALL getused(tgforpt2)
      tgforpt = tgforpt + (tgforpt2 - tgforpt1)

      RETURN
      END
