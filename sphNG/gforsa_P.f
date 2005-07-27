      SUBROUTINE gforsa(m, nlistga, listga, xyzmh, fsx, fsy, fsz, epot)
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

      DIMENSION listga(idim), xyzmh(5,idim)

      rrx = xyzmh(1,m)
      rry = xyzmh(2,m)
      rrz = xyzmh(3,m)

      DO 101 j = 1, nlistga
         n = listga(j)
         difx = xyzmh(1,n) - rrx
         dify = xyzmh(2,n) - rry
         difz = xyzmh(3,n) - rrz

         IF (isoft.EQ.1) THEN
            rr = difx**2 + dify**2 + difz**2 + psoft**2
         ELSE
            rr = difx**2 + dify**2 + difz**2 + tiny
         ENDIF

         rr05 = SQRT(rr)
c
c--The force definition
c
         IF (n.LE.natom) THEN
            pmassn = imfac(n)*xyzmh(4,n)
         ELSE
            pmassn = xyzmh(4,n)
         ENDIF
         potn = pmassn/rr05
         fff = potn/rr

         epot = epot + potn
         fsx = fsx + fff*difx
         fsy = fsy + fff*dify
         fsz = fsz + fff*difz
 101  CONTINUE

      RETURN
      END
