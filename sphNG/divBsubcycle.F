      SUBROUTINE divBsubcycle(nsubsteps,dtmax,npart,ntot,xyzmh,rho,
     &           Bxyzpnew,divBmax)
 
#ifdef MPIALL
#include "mpi_sup.h"
#endif
 
      INCLUDE 'idim'
#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      REAL dtmax,divBmax
      INTEGER npart
      DIMENSION xyzmh(5,mmax2)
      REAL*4 rho(idim2)
      DIMENSION Bxyzpnew(imhdevol,imhd2)
c      DIMENSION dBevolxyz(imhdevol,imhd2)

      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/call'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/eosq'

      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/divve'
      INCLUDE 'COMMONS/gradhterms'

      REAL Bxyzphalf(imhdevol,imhd2)
      REAL dBxyzp(imhdevol,imhd2)

      INTEGER nosweep
      REAL dr
      REAL dx,dy,dz

      REAL subcycle_dtimin, subcycle_dtimax
      INTEGER subcycle_istepmin, subcycle_istepmax

      REAL dtsize
      COMMON /implicit3/ vari2(3,idim2)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine divBsubcycle')

#ifdef MPIALL
      CALL MPI_TYPE_CONTIGUOUS(4, MPI_REAL8, i4REAL8, ierr)
      CALL MPI_TYPE_COMMIT(i4REAL8, ierr)
#endif
c
c--Set up constants in Code units
c
      tolerance = 1.0e-3
      dtimax = dtmax/imaxstep

c      print *, dtmax, imaxstep
c
c--Set errors to zero for iteration start
c
      numoscillations = 0
      numequal = 0
      numcomp = 0
      ipos = 1

      ihasghostcount = 0
      ihasghost = 0
      icompact = 0
      Bxyzmax = 0.0

      subcycle_dtimax = 0.0
      subcycle_istepmax = 0.0
      subcycle_dtimin = huge(subcycle_dtimin)
      subcycle_istepmin = huge(subcycle_istepmin)

      dtsize = huge(dtsize)
      vsigi = 0.0
      hi = huge(hi)

      ! Future optimization would be to calculate time step only once
      ! rather than on each substep

C$OMP PARALLEL default(none)
C$OMP& shared(icall,dtimax,dtmax,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij,varij2,vari2)
C$OMP& shared(xyzmh,dvtable,ddvtable,grwij,cnormk)
C$OMP& shared(nlst0,ihasghost,radkernel,gradhs,ntot)
C$OMP& shared(Bxyzpnew,Bxyzp,Bxyzphalf,dBxyzp,vsound,ncompactlocal)
C$OMP& private(n,i,iparttree)
C$OMP& private(rhoi, rho1i, B2i ,vs2i)
C$OMP& reduction(MAX:vsigi)
C$OMP& reduction(MIN:dtsize)
C$OMP& reduction(MIN:hi)
c
C$OMP DO SCHEDULE(static)
      DO n = 1, ncompactlocal
         i = ivar(3,n)
c
c--Needed for MPI code
c
         IF (i.GT.ntot) THEN
            iparttree = i + ntot + 2
         ELSE
            iparttree = i
         ENDIF

         rhoi = rho(i)
         rho1i = 1.0/rhoi

         B2i = Bxyzpnew(1,i)**2 + Bxyzpnew(2,i)**2 + Bxyzpnew(3,i)**2
         vs2i = vsound(i)**2 + B2i*rhoi

         vsigi = MAX(vsigi, sqrt(vs2i))

         IF (icall.EQ.1) THEN
            IF (isteps(i).EQ.0) THEN
               dtsize= MIN(dtsize, dtmax*1.0d-12)
            ELSE
               dtsize = MIN(dtsize, dtimax*isteps(i))
            ENDIF
         ELSEIF (n.LE.nlst0) THEN
            dtsize = MIN(dtsize, dtimax*isteps(i))
         ELSE
            dtsize = MIN(dtsize, dtimax/2.0*isteps(i))
         ENDIF

         hi = MIN(hi, xyzmh(5,iparttree))
      END DO  
C$OMP END DO
C$OMP END PARALLEL

#ifdef MPI
c
c--Need to get values of hi, dtsize and vsigi across all MPI processes
c
      CALL MPI_ALLREDUCE(dtsize,dtsizetot,1,MPI_REAL8,MPI_MIN,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(hi,hitot,1,MPI_REAL8,MPI_MIN,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(vsigi,vsigitot,1,MPI_REAL8,MPI_MAX,
     &     MPI_COMM_WORLD,ierr)
      dtsize = dtsizetot
      hi = hitot
      vsigi = vsigitot

c      print *,iproc,': dtsize,hi,vsigi ',dtsize,hi,vsigi
#endif

      ! make sure that for this wave speed, the Courant condition is
      ! upheld on all particles

      vs2i = vsigi * vsigi
      dtsize = MIN(dtsize, 0.5 * hi / vsigi)         
c      print *, '  dtsize=', dtsize, '  vsigi=', vsigi      

C$OMP PARALLEL default(none)
C$OMP& shared(icall,dtimax,dtmax,isteps,npart,hasghost,ireal,nghost)
C$OMP& shared(rho,vari,ivar,ijvar,varij,varij2,vari2)
C$OMP& shared(xyzmh,dvtable,ddvtable,grwij,cnormk)
C$OMP& shared(nlst0,ihasghost,radkernel,gradhs)
C$OMP& shared(Bxyzpnew,Bxyzp,Bxyzphalf,dBxyzp,vsound,ncompact)
C$OMP& shared(dtsize, vsigi, vs2i, ntot)
C$OMP& shared(iunique,iorig, iproc)
C$OMP& private(n,i,j,k,icompact,pmi,hj,hj21,hj41,grkerni,grkernj)
C$OMP& private(rxyi,ryzi,rxzi,rxxi,ryyi,rzzi)
C$OMP& private(dti,dx,dy,dz,iparttree)
C$OMP& private(rij2,rij,rij1,dr,pmj,rhoi,rhoj,hi,hi21,hi41)
C$OMP& private(v2,v,index,dxx,index1,dgrwdx,grwtij)
C$OMP& private(pmjdWrij1rhoj,runix,runiy,runiz,denom1)
C$OMP& private(dBxi,dByi,dBzi,projdB,rho1i)
C$OMP& reduction(+:ihasghostcount)
C$OMP& reduction(MAX:Bxyzmax)
C$OMP& reduction(MIN:subcycle_dtimin, subcycle_istepmin)
C$OMP& reduction(MAX:subcycle_dtimax, subcycle_istepmax)
c
c--Set up values that don't change during sweeps
c
C$OMP DO SCHEDULE(static)
      DO n = 1, ncompact
         i = ivar(3,n)
c
c--Needed for MPI code
c
         IF (i.GT.ntot) THEN
            iparttree = i + ntot + 2
         ELSE
            iparttree = i
         ENDIF

         IF(hasghost(i)) THEN
            ihasghostcount = ihasghostcount + 1
         ENDIF

         DO j = 1, imhdevol
            dBxyzp(j,i) = 0.
         END DO

         dBxi = 0.
         dByi = 0.
         dBzi = 0.

         pmi = xyzmh(4,iparttree)
         hi = xyzmh(5,iparttree)
         hi21 = 1./(hi*hi)
         hi41 = hi21*hi21
         rhoi = rho(i)
         rho1i = 1.0/rhoi

         DO k = 1, ivar(1,n)
            icompact = ivar(2,n) + k
            j = ijvar(icompact)

            dx = xyzmh(1,iparttree) - xyzmh(1,j)
            dy = xyzmh(2,iparttree) - xyzmh(2,j)
            dz = xyzmh(3,iparttree) - xyzmh(3,j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            rij = SQRT(rij2)
            rij1 = 1./rij
            dr = rij

            pmj = xyzmh(4,j)
            hj = xyzmh(5,j)
            hj21 = 1.0/(hj*hj)
            hj41 = hj21*hj21
            rhoj = rho(j)

            v2 = rij2*hi21
            v = rij/hi
            IF (v.LT.radkernel) THEN
               index = v2*ddvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dgrwdx = (grwij(index1) - grwij(index))*ddvtable
               grkerni = (grwij(index) + dgrwdx*dxx)*cnormk*hi41*
     &              gradhs(1,i)
            ELSE
               grkerni = 0.
            ENDIF

            v2 = rij2*hj21
            v = rij/hj
            IF (v.LT.radkernel) THEN
               index = v2*ddvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dgrwdx = (grwij(index1) - grwij(index))*ddvtable
               grkernj = (grwij(index) + dgrwdx*dxx)*cnormk*hj41*
     &              gradhs(1,j)
            ELSE
               grkernj = 0.
            ENDIF

            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1

c            B2j = Bxyz(1,j)**2 + Bxyz(2,j)**2 + Bxyz(3,j)**2
c--rough signal velocity only (no div v term)
c            vsigij = 0.5*(vsigi + sqrt(vsound(j)**2 + B2j/rhoj))
c            pmjdWrij1rhoj = pmj*dW*vsigij/rhoj

            varij(1,icompact) = grkerni * pmj * rho1i * rho1i
            varij(2,icompact) = grkerni * pmi * rho1i * rho1i
            varij(3,icompact) = grkernj * pmj / (rhoj * rhoj)
            varij(4,icompact) = rhoj

            varij2(1,icompact) = runix
            varij2(2,icompact) = runiy
            varij2(3,icompact) = runiz
         END DO

         vari(1,n) = rhoi

c don't need 0.5 divv term during substepping, only when part of full MHD equations
         vari2(1,n) = vs2i
         vari2(2,n) = 0.8*vsigi/hi !+ 0.5*divv(i)*rho1i
         vari2(3,n) = dtsize
      END DO
C$OMP END DO
C$OMP END PARALLEL

      IF (ihasghostcount.GE.1) ihasghost = 1
cC$OMP END SINGLE
cC$OMP DO SCHEDULE(static)
c      DO i = npart + 1, npart + nghost*ihasghost
c         j = ireal(i)
c         DO k = 1, 7
c            dvdx(k,i) = dvdx(k,j)
c         END DO
c      END DO
cC$OMP END DO
c
c--Begin iterating
c
      divBinitial = 0.

      DO nosweep = 1, nsubsteps
c
c--Calculate div B and dBpsi
c
         CALL divB_dBpsi(Bxyzpnew,ntot)
c
c--Can now calculate half-timestep values of B and psi
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompactlocal,ivar,vari2,Bxyzphalf,Bxyzpnew,varinew)
C$OMP& shared(dBxyzp,nosweep)
C$OMP& shared(dtsize)
C$OMP& private(n,i,dti,dpsidt)
C$OMP& reduction(MAX:divBinitial)
         DO n = 1, ncompactlocal
            i = ivar(3,n)
c
c--Do 1st bit of full timestep integration
c
            dti = dtsize /256.0

            Bxyzphalf(1,i) = Bxyzpnew(1,i) + varinew(1,i)*dti

            Bxyzphalf(2,i) = Bxyzpnew(2,i) + varinew(2,i)*dti
            Bxyzphalf(3,i) = Bxyzpnew(3,i) + varinew(3,i)*dti
            dpsidt = -vari2(1,n)*varinew(4,i) - vari2(2,n)*Bxyzpnew(4,i)
            Bxyzphalf(4,i) = Bxyzpnew(4,i) + dpsidt*dti

            DO j = 1,3
               dBxyzp(j,i) = dBxyzp(j,i) + varinew(j,i)*dti
            END DO
            dBxyzp(4,i) = dBxyzp(4,i) + dpsidt*dti


            IF (nosweep.EQ.1) THEN
               divBinitial = MAX(divBinitial,ABS(varinew(4,i)))
            ENDIF
c
c--Integrate to 1/2 timestep
c
            dti = 0.5 * dtsize

            Bxyzpnew(1,i) = Bxyzpnew(1,i) + varinew(1,i)*dti
            Bxyzpnew(2,i) = Bxyzpnew(2,i) + varinew(2,i)*dti
            Bxyzpnew(3,i) = Bxyzpnew(3,i) + varinew(3,i)*dti
            Bxyzpnew(4,i) = Bxyzpnew(4,i) + dpsidt*dti
         END DO
C$OMP END PARALLEL DO
c
c--Modify values on ghosts for half-timestep
c
         DO i = npart + 1, npart + nghost*ihasghost
            j = ireal(i)
            DO k = 1, 4
               Bxyzpnew(k,i) = Bxyzpnew(k,j)
            END DO
         END DO

#ifdef MPI
c
c--Update values of Bxyzpnew for copies of particles on MPI processes
c
      inumbertotallocal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
c     is because here it is assumed that the total number of remote particles
c     on which forces need to be calculated by this MPI process is less
c     than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
c
c--Only sends anything at all to another MPI process if it needs to (the
c     other MPI process knows when it doesn't need to expect anything)
c
                     IF (nneighsentany(j+1)) THEN
#ifdef MPIDEBUG
                        print *,iproc,': sending Bxyzpnew to ',j,
     &                     ' starting at ',llistsentback(1,j+1),'+1 ',
     &                       nneighsentback(j+1),' list '
c     &             ,(llistsentback(kkk,j+1),kkk=1,nneighsentback(j+1))
#endif
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i4REAL8,indexB,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i4REAL8,indexB,ierr)
                CALL MPI_TYPE_COMMIT(indexB,ierr)

                        CALL MPI_SEND(Bxyzpnew,1,indexB,j,154,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexB,ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent Bxyzpnew'
#endif
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
               print *,iproc,': receiving Bxyzpnew data ',
     &              inumberreturned,inumofreturns
#endif
               istart = ntot + inumbertotallocal + 1
               CALL MPI_RECV(Bxyzpnew(1,istart), idim, MPI_REAL8,
     &              i, 154, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
               print *,iproc,': got Bxyzpnew from ',iprocrec,inumber
#endif
               inumbertotallocal = inumbertotallocal + inumber/4
               IF (inumberindiv(inumberreturned)*4.NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumbercumm(inumberreturned).NE.istart) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm ',istart,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberproc ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit(1)
      ENDIF

#endif
c
c--Then calculate derivatives of B and psi at half-timestep
c
         CALL divB_dBpsi(Bxyzpnew,ntot)
c
c--Them calculate full-timestep values of B and psi
c
         f22 = 255./256.

         divBmax = 0.
         divBmean = 0.
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(ncompactlocal,ivar,vari2,Bxyzphalf,Bxyzpnew,dBxyzp)
C$OMP& shared(dtsize,varinew,f22)
C$OMP& private(n,i,dti,dpsidt)
C$OMP& reduction(MAX:divBmax)
C$OMP& reduction(+:divBmean)
         DO n = 1, ncompactlocal
            i = ivar(3,n)

            divBmax = MAX(divBmax,ABS(varinew(4,i)))
            divBmean = divBmean + ABS(varinew(4,i))
c
c--Do the 2nd bit of full timestep integration
c
            dti = f22 * dtsize

            Bxyzpnew(1,i) = Bxyzphalf(1,i) + varinew(1,i)*dti
            Bxyzpnew(2,i) = Bxyzphalf(2,i) + varinew(2,i)*dti
            Bxyzpnew(3,i) = Bxyzphalf(3,i) + varinew(3,i)*dti
            dpsidt = -vari2(1,n)*varinew(4,i) -vari2(2,n)*Bxyzpnew(4,i)
            Bxyzpnew(4,i) = Bxyzphalf(4,i) + dpsidt*dti

            DO j = 1,3
               dBxyzp(j,i) = dBxyzp(j,i) + varinew(j,i)*dti
            END DO
            dBxyzp(4,i) = dBxyzp(4,i) + dpsidt*dti
         END DO
C$OMP END PARALLEL DO
         divBmean = divBmean / ncompactlocal

c         DO n=1, ncompactlocal
c            i = ivar(3,n)
c            IF (ABS(varinew(4,i).EQ.divBmax)) THEN
c               print *,iunique(iorig(i)),divBmax
c               STOP
c            ENDIF
c         END DO

c
c--Modify values on ghosts for next timestep
c
         IF (nosweep.LT.nsubsteps) THEN
            DO i = npart + 1, npart + nghost*ihasghost
               j = ireal(i)
               DO k = 1, 4
                  Bxyzpnew(k,i) = Bxyzpnew(k,j)
               END DO
            END DO
         ENDIF

#ifdef MPI
c
c--Update values of Bxyzpnew for copies of particles on MPI processes
c
      inumbertotallocal = 0
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
            IF (nneighsentanyatall) THEN
c
c--Otherwise this process does not need to send any particle back!
c
c--Send active node data to be processed by other processes.  Data from all
c     other processes is received before any processing is done.  This
c     is because here it is assumed that the total number of remote particles
c     on which forces need to be calculated by this MPI process is less
c     than idim.
c
               DO j = 0, numproc - 1
                  IF (j.NE.iproc) THEN
c
c--Only sends anything at all to another MPI process if it needs to (the
c     other MPI process knows when it doesn't need to expect anything)
c
                     IF (nneighsentany(j+1)) THEN
#ifdef MPIDEBUG
                        print *,iproc,': sending Bxyzpnew to ',j,
     &                     ' starting at ',llistsentback(1,j+1),'+1 ',
     &                       nneighsentback(j+1),' list '
c     &             ,(llistsentback(kkk,j+1),kkk=1,nneighsentback(j+1))
#endif
c               CALL MPI_TYPE_CREATE_INDEXED_BLOCK(nneighsentback(j+1),
c     &             1,llistsentback(1,j+1),i4REAL8,indexB,ierr)
                CALL MPI_TYPE_INDEXED(nneighsentback(j+1),lblocklengths,
     &                llistsentback(1,j+1),i4REAL8,indexB,ierr)
                CALL MPI_TYPE_COMMIT(indexB,ierr)

                        CALL MPI_SEND(Bxyzpnew,1,indexB,j,154,
     &                       MPI_COMM_WORLD, ierr)

                        CALL MPI_TYPE_FREE(indexB,ierr)
#ifdef MPIDEBUG
                        print *,iproc,' sent Bxyzpnew'
#endif
                     ENDIF
                  ENDIF
               END DO
            ENDIF
c
c--Other processes receive the particles being sent
c
         ELSE
            IF (inumberreturned.LT.inumofreturns .AND.
     &           nneightogetback(i+1).GT.0) THEN
               inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
               print *,iproc,': receiving Bxyzpnew data ',
     &              inumberreturned,inumofreturns
#endif
               istart = ntot + inumbertotallocal + 1
               CALL MPI_RECV(Bxyzpnew(1,istart), idim, MPI_REAL8,
     &              i, 154, MPI_COMM_WORLD, istatus, ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, inumber, ierr)
               iprocrec = istatus(MPI_SOURCE)

#ifdef MPIDEBUG
               print *,iproc,': got Bxyzpnew from ',iprocrec,inumber
#endif
               inumbertotallocal = inumbertotallocal + inumber/4
               IF (inumberindiv(inumberreturned)*4.NE.inumber) THEN
                  WRITE (*,*) iproc,': ERR - inumberindiv ',inumber,
     &                 inumberindiv(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumbercumm(inumberreturned).NE.istart) THEN
                  WRITE (*,*) iproc,': ERR - inumbercumm ',istart,
     &                 inumbercumm(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
               IF (inumberproc(inumberreturned).NE.iprocrec) THEN
                  WRITE (*,*) iproc,': ERR - inumberproc ',iprocrec,
     &                 inumberproc(inumberreturned),inumberreturned
                  CALL quit(1)
               ENDIF
            ENDIF
         ENDIF
      END DO
      IF (inumbertotallocal.NE.inumbertotal) THEN
         WRITE (*,*) iproc,': inumbertotallocal.NE.inumbertotal ',
     &        inumbertotallocal,inumbertotal
         CALL quit(1)
      ENDIF

#endif

c         print *,'it ',nosweep,divBmax,divBmean/ncompact,ncompact

c      print *,'Mean/Max divergence is ',nosweep,divBmax,divBmean,
c     &        divBinitial
      END DO ! End loop over nsubsteps

c
c--Then calculate derivatives of B and psi at the end
c
      CALL divB_dBpsi(Bxyzpnew,ntot)

#ifdef MPIALL
      CALL MPI_TYPE_FREE(i4REAL8,ierr)
#endif

      RETURN

      END





      SUBROUTINE divB_dBpsi(Bxyzpnew,ntot)

      INCLUDE 'idim'
#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/implicit'
      INCLUDE 'COMMONS/implicit2'
      INCLUDE 'COMMONS/current'

      INCLUDE 'COMMONS/sort'

      REAL Bxyzpnew(imhdevol,imhd2)

c
c--Zero arrays for storing div B and dBpsi
c
C$OMP PARALLEL default(none)
C$OMP& shared(ivar,ijvar,vari,varij,varij2)
C$OMP& shared(Bxyzpnew,iscurrent,varinew,ncompact)
C$OMP& shared(iunique,iorig,iproc,ikeep)
C$OMP& private(i,j,k,n,divBi,dBpsix,dBpsiy,dBpsiz,rhoi,rhoj)
C$OMP& private(icompact,pmjdWirho2i,pmidWirho2i,pmjdWjrho2j)
C$OMP& private(runix,runiy,runiz,dBx,dBy,dBz,projdB,tempval)
C$OMP& private(pmjdWirho1i)
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)
         DO k = 1, 4
            varinew(k,i) = 0.
         END DO
      END DO
C$OMP END DO
c
c--Calculate divB values without using separate subroutine
c
c--Particle I loop
c
C$OMP DO SCHEDULE(runtime)
      DO n = 1, ncompact
         i = ivar(3,n)

         divBi = 0.
         dBpsix = 0.
         dBpsiy = 0.
         dBpsiz = 0.

         rhoi = vari(1,n)
c
c--All the neighbours loop
c
         DO k = 1, ivar(1,n)
            icompact = ivar(2,n) + k
            j = ijvar(icompact)

            pmjdWirho2i = varij(1,icompact)
            pmjdWirho1i = pmjdWirho2i * rhoi
            pmidWirho2i = varij(2,icompact)
            pmjdWjrho2j = varij(3,icompact)
            rhoj = varij(4,icompact)

            runix = varij2(1,icompact)
            runiy = varij2(2,icompact)
            runiz = varij2(3,icompact)
c
c--Calculate div B (only contributions from within hi)
c
            dBx = Bxyzpnew(1,i)*rhoi - Bxyzpnew(1,j)*rhoj
            dBy = Bxyzpnew(2,i)*rhoi - Bxyzpnew(2,j)*rhoj
            dBz = Bxyzpnew(3,i)*rhoi - Bxyzpnew(3,j)*rhoj
            projdB = dBx*runix + dBy*runiy + dBz*runiz

            divBi = divBi - pmjdWirho1i*projdB
c
c--Calculate dBpsi
c
c--Add terms from particle j within hi
c
            tempval = pmjdWirho2i*Bxyzpnew(4,i)
            dBpsix = dBpsix - tempval*runix
            dBpsiy = dBpsiy - tempval*runiy
            dBpsiz = dBpsiz - tempval*runiz
c
c--Add terms to particle j (if j is active)
c
            IF (iscurrent(j)) THEN
               tempval = pmidWirho2i*Bxyzpnew(4,i)
C$OMP ATOMIC
               varinew(1,j) = varinew(1,j) + tempval*runix
C$OMP ATOMIC
               varinew(2,j) = varinew(2,j) + tempval*runiy
C$OMP ATOMIC
               varinew(3,j) = varinew(3,j) + tempval*runiz
c
c--Add terms to particle i from inactive j particle
c
            ELSE
               tempval = pmjdWjrho2j*Bxyzpnew(4,j)
               dBpsix = dBpsix - tempval*runix
               dBpsiy = dBpsiy - tempval*runiy
               dBpsiz = dBpsiz - tempval*runiz
            ENDIF
         END DO                 !J-loop
c
c--Store div B and dBpsi
c
C$OMP ATOMIC
         varinew(1,i) = varinew(1,i) + dBpsix
C$OMP ATOMIC
         varinew(2,i) = varinew(2,i) + dBpsiy
C$OMP ATOMIC
         varinew(3,i) = varinew(3,i) + dBpsiz

         varinew(4,i) = divBi

      END DO
C$OMP END DO
C$OMP END PARALLEL

#ifdef MPI
c
c--Need to add contributions to varinew(1-4) from other MPI processes
c
      istart = ntot + 1
      istartrec = istart + inumbertotal
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back quantities to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' varinew ',j
#endif
               CALL MPI_RECV(varinew(1,istartrec),
     &              nelements_varinew*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,153,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.
     &              nelements_varinew*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nelements',
     &                 '*nnsentback varinw'
                  CALL quit(1)
               ENDIF

               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1

                  DO k = 1, 4
                     varinew(k,ipos) = varinew(k,ipos) + varinew(k,jpos)
                  END DO
               END DO

#ifdef MPIDEBUG
               print *,iproc,': got varinew from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
#endif

            END DO
c
c--Other processes send the values back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                  print *,iproc,': sending varinew to ',iprocsend,
     &                 ' istart ',istart,' returned ',inumberreturned,
     &                 ' of ',inumofreturns
#endif
                  CALL MPI_SEND(varinew(1,istart),
     &                 nelements_varinew*inumber,MPI_REAL8,
     &                 iprocsend,153,MPI_COMM_WORLD,ierr)
c
c--Only need to loop over jjj until found appropriate index for process 'i'
c     then can jump out of loop
c
                  GOTO 7703
               ENDIF
            END DO
 7703       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO
#endif

      RETURN
      END
