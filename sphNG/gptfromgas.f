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

      DIMENSION xyzmh(5,mmax2)
      DIMENSION fxyzu(4,idim3)

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
      INCLUDE 'COMMONS/rbnd'

      CALL getused(tgforpt1)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine gptfromgas')

      IF (iphase(ipart).LE.0 .OR. iphase(ipart).GE.10) THEN
         WRITE (*,*) 'ERROR - gptfromgas called by non-sink'
         CALL quit(1)
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
         IF(iphase(j).LT.0 .OR. 
     &        iphase(j).GE.1 .AND. iphase(j).LT.10) GOTO 101

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
         IF (iphase(ipart).NE.5) THEN
            rr05 = SQRT(rr)
            fff = pmassj/(rr*rr05)
            potn = pmassj/rr05
         ELSE
            rr05 = SQRT(rr)
            rsurface = xyzmh(5,ipart)*pradfac(listrealpm(ipart))
            IF (rr05.LE.(2.*rsurface)) THEN
               fsurface = (((2.*rsurface)-rr05)/
     &              (rsurface))**4
            ELSE
               fsurface = 0.0
            ENDIF
            fff = (1.0-fsurface)*pmassj/(rr*rr05)
            potn = pmassj/rr05
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
