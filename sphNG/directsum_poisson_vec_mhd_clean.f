      SUBROUTINE directsum_poisson_vec(i,ntot,
     &           xyzmh,sourcexyz,
     &           curlAx,curlAy,curlAz)
c************************************************************
c                                                           *
c  Calculates solution to the vector poisson equation       *
c                                                           *
c   del^2 A = 4*pi*J                                        *
c                                                           *
c  (where A and J are both vector quantities)               * 
c                                                           *
c   by direct summation over the particles                  *
c  Returns curl A for a particular particle i     *
c							    *
c  Subroutine by D. Price (2003)                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/soft'
      
      DIMENSION xyzmh(5,idim)
      DIMENSION sourcexyz(3,idim)

      xi = xyzmh(1,i)
      yi = xyzmh(2,i)
      zi = xyzmh(3,i)
      hi2 = 0.25*xyzmh(5,i)**2 ! use h/2
      curlAx = 0.
      curlAy = 0.
      curlAz = 0.
!      Ax = 0.
!      Ay = 0.
!      Az = 0.

      DO j = 1, ntot
	 IF (j.NE.i) THEN
	    dx = xyzmh(1,j) - xi
            dy = xyzmh(2,j) - yi
            dz = xyzmh(3,j) - zi
c
c--define radius (with softenings etc)
c
c            IF (isoft.EQ.1) THEN
c               rr = dx**2 + dy**2 + dz**2 + psoft**2
c            ELSE 
               rr = dx**2 + dy**2 + dz**2 + hi2
c            ENDIF
c	    rr = dx**2. + dy**2. + dz**2.

            rr05 = SQRT(rr)
c
c--calculate curl A
c	 
	    denom = 1./(rr*rr05)
	    
	    curlAx = curlAx - (sourcexyz(2,j)*dz 
     &                       - sourcexyz(3,j)*dy)*denom
	    curlAy = curlAy - (sourcexyz(3,j)*dx
     &                       - sourcexyz(1,j)*dz)*denom
	    curlAz = curlAz - (sourcexyz(1,j)*dy
     &                       - sourcexyz(2,j)*dx)*denom
c
c--calculate A
c	    
!            drr05 = 1./rr05
!	    Ax = Ax + sourcex(j)*drr05
!	    Ay = Ay + sourcey(j)*drr05
!	    Az = Az + sourcez(j)*drr05
	    
         ENDIF
      ENDDO

      RETURN
      END
