      SUBROUTINE divBclean(nlst_in,nlst_end,npart,ntot,list,xyzmh,rho,
     &                     Bevolxyz)
c************************************************************
c                                                           *
c   This subroutine cleans up the divergence of the         *
c   magnetic field using a current projection               *
c                                                           *
c   given a field Bx, By, Bz, it calculates the current     *
c   using SPH operators, then solves the vector             *
c   Poisson equation:                                       *
c                                                           *
c   del^2 A* = curl B*                                      *
c                                                           *
c   using the tree, and corrects the magnetic field using:  *
c                                                           *
c   B = curl A*                                             *
c                                                           *
c   This gives a divergence-free projection    *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      DIMENSION xyzmh(5,idim)
      REAL*4 rho(idim)
      DIMENSION Bevolxyz(3,imhd)
      DIMENSION curlBmrhoxyz(3,imhd)
      DIMENSION list(idim)

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/varmhd'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine divBclean')
c
c--get B from B/rho if necessary
c
c      IF (varmhd.eq.'Bvol') THEN
c         DO i=1,npart
c            Bxyz(1,i) = Bevolxyz(1,i)
c            Bxyz(2,i) = Bevolxyz(2,i)
c            Bxyz(3,i) = Bevolxyz(3,i)
c         ENDDO
c      ELSEIF (varmhd.EQ.'Brho') THEN
c         DO i=1,npart
c            Bxyz(1,i) = Bevolxyz(1,i)*rho(ipart)
c            Bxyz(2,i) = Bevolxyz(2,i)*rho(ipart)
c            Bxyz(3,i) = Bevolxyz(3,i)*rho(ipart)
c         ENDDO
c      ELSEIF (varmhd.NE.'eulr') THEN
c         STOP 'unknown MHD variable in derivi'
c      ENDIF
c
c--check that the current is non-zero (calculated in last force call)
c
      IF (abs(divcurlB(2,1)).LT.tiny) THEN
         WRITE(*,*) 'WARNING: current = 0 in call to divBclean'
         DO i=1,10
            WRITE(iprint,*) divcurlB(2,i),divcurlB(3,i),divcurlB(4,i)
         ENDDO
         RETURN
      ENDIF
     
      WRITE(iprint,99005)
99005 FORMAT (' INVOKING VECTOR PROJECTION ON DIV B...')

c
c--first of all calculate the curl of B for *every* particle
c  (returns the source term m*curlB/(4*pi*rho) required in the
c   Poisson equation)
c
      DO i=1,npart
         curlBmrhoxyz(1,i) = xyzmh(4,i)*divcurlB(2,i)/(4.*pi*rho(i))
         curlBmrhoxyz(2,i) = xyzmh(4,i)*divcurlB(3,i)/(4.*pi*rho(i))
         curlBmrhoxyz(3,i) = xyzmh(4,i)*divcurlB(4,i)/(4.*pi*rho(i))
      ENDDO       
c
c--then calculate the correction to B by solution of a Poisson equation
c
c--make the tree  
c
c          WRITE(iprint,*) 'making tree...',ntot,npart        
c          CALL mtree_vec(ntot,npart,x,y,z,pmass,
c     &                   curlBmrhox,curlBmrhoy,curlBmrhoz,h)
c 
c--climb the tree to get forces (ie. grad psi) - these come back as gravx,y,z from insulate
c
c          WRITE(iprint,*) 'climbing tree...'
c          accdivB = 0.5
c
c--Compute the neighbour indexes & gravitational forces of the distant 
c     particles for all the particles in the list
c
cC$OMP PARALLEL default(none), shared(nlst,npart,accdivB)
cC$OMP& shared(h,Bx,By,Bz,poten,dphit,iphase)
cC$OMP& private(ipart,fsx,fsy,fsz,potx,poty,potz)
cC$OMP DO SCHEDULE(runtime)
c         DO ipart = 1, npart
         !!print*,ipart,' Bprev = ',Bx(ipart),By(ipart),Bz(ipart)
c
c--Walk through the tree to get the neighbours, the acceleration
c     and the potential due to outside 2h particles
c
c            CALL treef_vec(ipart,npart,h,accdivB,1,
c     &           Bcorrxi,Bcorryi,Bcorrzi,potx,poty,potz)
c
c            Bevolxyz(1,i) = Bcorrxi + Bextx
c            Bevolxyz(2,i) = Bcorryi + Bexty
c            Bevolxyz(3,i) = Bcorrzi + Bextz
c            IF (varmhd.EQ.'Brho') THEN
c               Bevolxyz(1,i) = Bevolxyz(1,i)/rho(i)
c               Bevolxyz(2,i) = Bevolxyz(2,i)/rho(i)
c               Bevolxyz(3,i) = Bevolxyz(3,i)/rho(i)
c            ENDIF
c         ENDDO
cC$OMP END DO
cC$OMP END PARALLEL

       WRITE(iprint,*) 'getting corrected field by direct sum...',
     &      nlst_end-nlst_in+1
       DO i=nlst_in,nlst_end
          ipart = list(i)
          CALL directsum_poisson_vec(ipart,npart,xyzmh,curlBmrhoxyz,
     &            Bcorrxi,Bcorryi,Bcorrzi)
c          if (ipart.le.10) then
c          print*,ipart,'curlB=  ',divcurlB(:,ipart)
c          IF (varmhd.EQ.'Brho') THEN
c             print*,'old=  ',Bevolxyz(:,ipart)*rho(ipart)
c          ELSE
c             print*,'old = ',Bevolxyz(:,ipart)   
c          ENDIF
c          print*,'new = ',Bcorrxi + Bextx,Bcorryi+Bexty,Bcorrzi+Bextz
c          read*
c          endif
          Bevolxyz(1,ipart) = Bcorrxi + Bextx
          Bevolxyz(2,ipart) = Bcorryi + Bexty
          Bevolxyz(3,ipart) = Bcorrzi + Bextz
          IF (varmhd.EQ.'Brho') THEN
             Bevolxyz(1,ipart) = Bevolxyz(1,ipart)/rho(ipart)
             Bevolxyz(2,ipart) = Bevolxyz(2,ipart)/rho(ipart)
             Bevolxyz(3,ipart) = Bevolxyz(3,ipart)/rho(ipart)
          ENDIF
       ENDDO
       WRITE(iprint,*) '...done, Bext = ',Bextx,Bexty,Bextz
c
c--update these on ghosts
c
       DO i=npart+1,ntot
          j = ireal(i)
          Bevolxyz(1,i) = Bevolxyz(1,j)
          Bevolxyz(2,i) = Bevolxyz(2,j)
          Bevolxyz(3,i) = Bevolxyz(3,j)
       ENDDO
c
c--remake the tree with mass so as not to screw up the rest of the code
c      
c      CALL insulate(1,ntot,npart,x,y,z,pmass,h)
c      WRITE(*,99002)

      RETURN
      END
