      SUBROUTINE gforsn(m,ntot,nlistgn,listgn,xyzmh,fsx,fsy,fsz,epot)
c************************************************************
c                                                           *
c  Subroutine by W. Press (11/21/86).  Evaluates force on   *
c     particle M due to a list of composite nodes.          *
c     THIS ROUTINE VECTORIZABLE.                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/phase'

      DIMENSION listgn(idim), xyzmh(5,mmax)

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

      DO j = 1, nlistgn
         n = listgn(j)
         difx = xyzmh(1,n) - rrx
         dify = xyzmh(2,n) - rry
         difz = xyzmh(3,n) - rrz

         IF (isoft.EQ.1) THEN
            rr = difx**2 + dify**2 + difz**2 + psoft**2
         ELSE 
            rr = difx**2 + dify**2 + difz**2 + tiny
         ENDIF

         IF (iphase(m).GE.1 .AND. rr.LT.ptsoft**2) THEN
            WRITE (*,*) 'ptsoft gforsn ',m,rr,ptsoft
            STOP
         ENDIF

c
c--The force definition
c         
         rri = 1./rr
         rri05 = SQRT(rri)
         f = rri*rri05
         fpr = ( - 3.0)*f*rri
         fpprr = ( - 4.0)*fpr*rri

         tx = qrad(2,n)*difx + qrad(5,n)*dify + qrad(7,n)*difz
         ty = qrad(5,n)*difx + qrad(3,n)*dify + qrad(6,n)*difz
         tz = qrad(7,n)*difx + qrad(6,n)*dify + qrad(4,n)*difz
         rqr = difx*tx + dify*ty + difz*tz

         IF (n.LE.natom) THEN
            pmassn = imfac(n)*xyzmh(4,n)
         ELSE
            pmassn = xyzmh(4,n)
         ENDIF
         fac = pmassn*f
     &              + 0.5*(rqr*(fpprr-fpr/rr)+
     &               (qrad(2,n)+qrad(3,n)+qrad(4,n))*fpr)
         phim = pmassn*rri05
         phiq = 0.5*rqr*rri*f
         fsx = fsx + fac*difx 
     &             + fpr*tx
         fsy = fsy + fac*dify
     &             + fpr*ty
         fsz = fsz + fac*difz
     &             + fpr*tz
         epot = epot + phim
     &             + phiq
      END DO

      RETURN
      END
