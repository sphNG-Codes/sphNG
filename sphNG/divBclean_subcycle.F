      SUBROUTINE divBclean_subcycle(dtmax,npart,ntot,xyzmh,rho,Bevolxyz)
c************************************************************
c                                                           *
c   This subroutine cleans up the divergence of the         *
c   magnetic field using subcycling of the SPH version      *
c   of the Dedner et al. (2002) method.                     *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/varmhd'

      DIMENSION xyzmh(5,mmax2)
      REAL*4 rho(idim2)
      DIMENSION Bevolxyz(imhdevol,imhd2)
      DIMENSION Bxyzporig(imhdevol,imhd)
c     DIMENSION dBevolxyz(imhdevol,imhd2)

      totalmagenergy = 0.

C$OMP DO SCHEDULE(static)
      DO n = 1, ncompact
         i = ivar(3,n)
c
c--NOTE: Operates on Bevolxyz, not Bxyz
c
         Bxyzporig(1,i) = Bevolxyz(1,i)
         Bxyzporig(2,i) = Bevolxyz(2,i)
         Bxyzporig(3,i) = Bevolxyz(3,i)
         Bxyzporig(4,i) = Bevolxyz(4,i)

c         totalmagenergy = totalmagenergy + xyzmh(4,i)*(Bxyz(1,i)**2 +
c     &        Bxyz(2,i)**2 + Bxyz(3,i)**2)/rho(i)
      END DO
C$OMP END DO

      divBmax = 0.

      CALL divBsubcycle(ndivBsubcycles,dtmax,npart,ntot,xyzmh,rho,
     &     Bevolxyz,divBmax)

c sanity check that substepping hasn't massively screwed up the field
      IF (divBmax.GT.1.0e+10) THEN
C$OMP DO SCHEDULE(static)
            ! reset field to pre substepping state
            DO n = 1, ncompact
               i = ivar(3,n)

               Bevolxyz(1,i) = Bxyzporig(1,i)
               Bevolxyz(2,i) = Bxyzporig(2,i)
               Bevolxyz(3,i) = Bxyzporig(3,i)
               Bevolxyz(4,i) = Bxyzporig(4,i)
            END DO
C$OMP END DO
            WRITE (*,*) 'divB clean FAILED '
      ENDIF

c      WRITE (*,*) 'divB clean end ',divBmax,ndivBsubcycles
c
c--DO not need to replace divcurlB() since it has not yet been
c         calculated (calculated in forcei)
c
      RETURN
      END
