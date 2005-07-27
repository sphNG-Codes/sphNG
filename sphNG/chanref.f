      SUBROUTINE chanref(idir)
c************************************************************
c                                                           *
c  This subroutine transforms the various quantities        *
c     between two frame of references.                      *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/typef'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine chanref')
c
c--Scaling factors
c
      CALL scaling(gt, rscale, drdt, dlnrdt)
      r2 = drdt/rscale**2
      sdens = 1./rscale**3
      sergg = 1./rscale
      or = omeg0*rscale
      or1 = omeg0/rscale
c
c--Transform into original unit system
c
      IF ( idir.NE.1 ) THEN
         IF ( ifcor.EQ.1 ) THEN
            DO 50 i = 1, npart
               vxyzu(1,i) = vxyzu(1,i)/rscale - xyzmh(1,i)*r2 + 
     &              or1*xyzmh(2,i)
               vxyzu(2,i) = vxyzu(2,i)/rscale - xyzmh(2,i)*r2 - 
     &              or1*xyzmh(1,i)
               vxyzu(3,i) = vxyzu(3,i)/rscale - xyzmh(3,i)*r2

               xyzmh(1,i) = xyzmh(1,i)/rscale
               xyzmh(2,i) = xyzmh(2,i)/rscale
               xyzmh(3,i) = xyzmh(3,i)/rscale
               xyzmh(5,i) = xyzmh(5,i)/rscale

               rho(i) = rho(i)/sdens
               vxyzu(4,i) = vxyzu(4,i)/sergg
 50         CONTINUE
         ELSE IF ( ifcor.EQ.2 ) THEN
            DO 80 i = 1, npart
               vxyzu(1,i) = vxyzu(1,i)/rscale - xyzmh(1,i)*r2 
               vxyzu(2,i) = vxyzu(2,i)/rscale - xyzmh(2,i)*r2 + 
     &              or1*xyzmh(3,i)
               vxyzu(3,i) = vxyzu(3,i)/rscale - xyzmh(3,i)*r2 - 
     &              or1*xyzmh(2,i)

               xyzmh(1,i) = xyzmh(1,i)/rscale
               xyzmh(2,i) = xyzmh(2,i)/rscale
               xyzmh(3,i) = xyzmh(3,i)/rscale
               xyzmh(5,i) = xyzmh(5,i)/rscale

               rho(i) = rho(i)/sdens
               vxyzu(4,i) = vxyzu(4,i)/sergg
 80         CONTINUE
         ENDIF
         RETURN
c
c--Transform into expanding frame units
c
      ELSE
         IF ( ifcor.EQ.1 ) THEN
            DO 100 i = 1, npart
               vxyzu(1,i) = rscale*vxyzu(1,i) + xyzmh(1,i)*drdt - 
     &              or*xyzmh(2,i)
               vxyzu(2,i) = rscale*vxyzu(2,i) + xyzmh(2,i)*drdt + 
     &              or*xyzmh(1,i)
               vxyzu(3,i) = rscale*vxyzu(3,i) + xyzmh(3,i)*drdt

               xyzmh(1,i) = xyzmh(1,i)*rscale
               xyzmh(2,i) = xyzmh(2,i)*rscale
               xyzmh(3,i) = xyzmh(3,i)*rscale
               xyzmh(5,i) = xyzmh(5,i)*rscale

               rho(i) = rho(i)*sdens
               vxyzu(4,i) = vxyzu(4,i)*sergg
               poten(i) = poten(i)*sergg
 100        CONTINUE
         ELSE IF ( ifcor.EQ.2 ) THEN
            DO 120 i = 1, npart
               vxyzu(1,i) = rscale*vxyzu(1,i) + xyzmh(1,i)*drdt 
               vxyzu(2,i) = rscale*vxyzu(2,i) + xyzmh(2,i)*drdt - 
     &              or*xyzmh(3,i)
               vxyzu(3,i) = rscale*vxyzu(3,i) + xyzmh(3,i)*drdt + 
     &              or*xyzmh(2,i)

               xyzmh(1,i) = xyzmh(1,i)*rscale
               xyzmh(2,i) = xyzmh(2,i)*rscale
               xyzmh(3,i) = xyzmh(3,i)*rscale
               xyzmh(5,i) = xyzmh(5,i)*rscale

               rho(i) = rho(i)*sdens
               vxyzu(4,i) = vxyzu(4,i)*sergg
               poten(i) = poten(i)*sergg
 120        CONTINUE
         ENDIF
         RETURN
      ENDIF

      END
