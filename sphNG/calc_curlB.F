      SUBROUTINE calc_curlB(nlst_in,nlst_tot,ntot,npart,xyzmh,Bxyz,trho)
c*****************************************************************
c                                                                *
c     Calculate J = curl B for use in non-ideal MHD              *
c                                                                *
c*****************************************************************

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

      DIMENSION xyzmh(5,mmax2), Bxyz(3,imhd2)
      REAL*4 trho(idim2)

      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/table'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
#endif
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/compact'

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,ntot,nlst_in,nlst_tot,xyzmh,Bxyz,jcurrent,iphase)
C$OMP& shared(radkernel,dvtable,ddvtable,grwij)
C$OMP& shared(gradhs,cnormk,trho,ivar,ijvar)
C$OMP& private(n,ipart,iparttree,j,k,xi,yi,zi,pmassi)
C$OMP& private(hi,hi1,hi21,hi31,hi41,Bxi,Byi,Bzi)
C$OMP& private(dx,dy,dz,pmassj,rij2,v2,rij1,dBx,dBy,dBz,index,index1)
C$OMP& private(dxx,dgrwdx,grwtij,numneigh,ioffset)
C$OMP& private(gradhi,grkerntable,grkerni,grpmi,rhoi,rhoi1)
      DO n = nlst_in, nlst_tot
         numneigh = ivar(1,n)
         ioffset  = ivar(2,n)
         ipart    = ivar(3,n)

         DO k = 1, 3
            jcurrent(k,ipart) = 0.
         END DO

         IF (ipart.GT.npart) THEN
            iparttree = ipart + ntot + 2
         ELSE
            iparttree = ipart
         ENDIF
         IF (iphase(ipart).EQ.0) THEN
            xi     = xyzmh(1,iparttree)
            yi     = xyzmh(2,iparttree)
            zi     = xyzmh(3,iparttree)
            pmassi = xyzmh(4,iparttree)
            hi     = xyzmh(5,iparttree)

            hi1  = 1./hi
            hi21 = hi1*hi1
            hi31 = hi21*hi1
            hi41 = hi21*hi21

            Bxi = Bxyz(1,ipart)
            Byi = Bxyz(2,ipart)
            Bzi = Bxyz(3,ipart)

            gradhi = gradhs(1,ipart)
            rhoi   = trho(ipart)
            rhoi1  = 1.0/rhoi

            DO k = 1, ivar(1,n)
               j = ijvar(ioffset + k)
               IF (iphase(j).EQ.0) THEN
                  dx     = xi - xyzmh(1,j)
                  dy     = yi - xyzmh(2,j)
                  dz     = zi - xyzmh(3,j)
                  pmassj =      xyzmh(4,j)
                  rij2   = dx*dx + dy*dy + dz*dz + tiny
                  v2     = rij2*hi21

                  IF (v2.LT.radkernel**2) THEN
                     rij1 = SQRT(rij2)
                     dBx  = Bxi - Bxyz(1,j)
                     dBy  = Byi - Bxyz(2,j)
                     dBz  = Bzi - Bxyz(3,j)
c
c--Get kernel quantities from interpolation in table
c
                     index  = v2*ddvtable
                     dxx    = v2 - index*dvtable
                     index1 = index + 1
                     IF (index.GE.itable) THEN
                        index  = itable
                        index1 = itable
                     ENDIF
                     dgrwdx = (grwij(index1) - grwij(index))*ddvtable
                     grkerntable = (grwij(index)+ dgrwdx*dxx)*cnormk
                     grkerni     = grkerntable*hi41*gradhi
                     grpmi       = grkerni*pmassj/rij1
                     jcurrent(1,ipart) = jcurrent(1,ipart) 
     &                    + grpmi*(dBy*dz - dBz*dy)
                     jcurrent(2,ipart) = jcurrent(2,ipart) 
     &                    + grpmi*(dBz*dx - dBx*dz)
                     jcurrent(3,ipart) = jcurrent(3,ipart) 
     &                    + grpmi*(dBx*dy - dBy*dx)
                  ENDIF
               ENDIF
            END DO
            jcurrent(1,ipart) = jcurrent(1,ipart)*rhoi1
            jcurrent(2,ipart) = jcurrent(2,ipart)*rhoi1
            jcurrent(3,ipart) = jcurrent(3,ipart)*rhoi1
         ENDIF
      END DO
C$OMP END PARALLEL DO
#ifdef MPI
c
c--For MPI job, need to add contributions to jcurrent 
c     from other MPI processes
c
c
c--Now need to return contributions from neighbours on remote nodes.
c
      istart = ntot + 1
      istartrec = istart + inumbertotal
      IF (istartrec + maxnneighsentback.GT.2*idim) THEN
         WRITE (*,*) 'ERROR - istartrec + maxnneighsentback.GT.2*idim'
         CALL quit(1)
      ENDIF
      inumberreturned = 0
      DO i = 0, numproc - 1
         IF (iproc.EQ.i) THEN
c
c--Receive back jcurrent to be added on to local values
c
            DO j = 1, inumofsends
#ifdef MPIDEBUG
               print *,iproc,': expecting to get ',inumofsends,
     &              ' jcurrent contributions ',j,
     &              ' with maxnneighsentback ',maxnneighsentback
#endif
               CALL MPI_RECV(jcurrent(1,istartrec),3*maxnneighsentback,
     &              MPI_REAL8,MPI_ANY_SOURCE,215,MPI_COMM_WORLD,istatus,
     &              ierr)
               CALL MPI_GET_COUNT(istatus, MPI_REAL8, ireturned, ierr)
               IF (istartrec+ireturned.GT.idim2) THEN
                  WRITE (*,*) 'ERROR - istartrec+ireturned.GT.idim2 ',
     &                 istartrec,ireturned
                  CALL quit(1)
               ENDIF
               iprocrec = istatus(MPI_SOURCE)
               IF (ireturned.NE.3*nneighsentback(iprocrec+1)) THEN
                  WRITE (*,*) 'ERROR - ireturned.NE.nnsentback jcurr ',
     &                 ireturned,nneighsentback(iprocrec+1)
                  CALL quit(1)
               ENDIF
#ifdef MPIDEBUG
               print *,iproc,': got jcurrent from ',iprocrec,
     &              ' put into ',llistsentback(1,iprocrec+1)+1
               print *,iproc,': received values ',
     &              (jcurrent(1,j3),j3=istartrec,
     &              istartrec+1)
#endif

C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nneighsentback,iprocrec,llistsentback,jcurrent)
C$OMP& shared(istartrec,iphase)
C$OMP& private(l,ipos,jpos,k)
               DO l = 1, nneighsentback(iprocrec+1)
                  ipos = llistsentback(l,iprocrec+1)+1
                  jpos = istartrec + l - 1
                  DO k = 1, 3
                     jcurrent(k,ipos) = jcurrent(k,ipos) + 
     &                    jcurrent(k,jpos)
                  END DO
               END DO
C$OMP END PARALLEL DO
            END DO
c
c--Other processes send the jcurrent values back
c
         ELSE
            DO jjj = 1, inumofreturns
               IF (inumberproc(jjj).EQ.i) THEN
                  istart = inumbercumm(jjj)
                  inumber = inumberindiv(jjj)
                  iprocsend = inumberproc(jjj)
                  inumberreturned = inumberreturned + 1
#ifdef MPIDEBUG
                  print *,iproc,': sending jcurrent to ',iprocsend,
     &                 ' istart ',istart,' number ',inumber,
     &                 ' returned ',inumberreturned,
     &                 ' of ',inumofreturns
                  print *,iproc,': sending values ',
     &                 (jcurrent(1,j2),j2=istart,istart+1)
#endif
                  CALL MPI_SEND(jcurrent(1,istart),3*inumber,MPI_REAL8,
     &                 iprocsend,215,MPI_COMM_WORLD, ierr)

#ifdef MPIDEBUG
                  print *,iproc,': sent jcurrent to ',
     &                 iprocsend,' and has sent ',inumber,' values'
#endif
                  GOTO 8765
               ENDIF
            END DO
 8765       CONTINUE
         ENDIF
c
c--If there are no particles sent, don't need to do anything (no MPI_SEND
c     because none will be expected)
c
      END DO

      IF (inumberreturned.NE.inumofreturns) THEN
         WRITE (*,*) 'ERRORJ - inumberreturned.NE.inumofreturns ',iproc,
     &        inumberreturned,inumofreturns
         CALL quit(1)
      ENDIF

#endif
c
c--Update J for inactive particles that are neighbours of active particles
c
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(npart,ntot,ireal,jcurrent)
C$OMP& private(i,j,k)
      DO i = npart + 1, ntot
         j = ireal(i)
         DO k = 1, 3
            jcurrent(k,i) = jcurrent(k,j)
         END DO
      END DO
C$OMP END PARALLEL DO

      RETURN
      END
