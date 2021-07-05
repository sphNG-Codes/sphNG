      SUBROUTINE dust_growth(nlst_in,nlst_end,npart,ntot,
     &     xyzmh,ilenx,vxyzu,ilenv,drgrain,ilendr)
c***********************************************************************
c                                                                      *
c     This subroutine calculates the dust growth in a poly disperse    *
c     population. We need the double loop structure over types         *
c     11 (growing) and types 12 (non-growing) because types 11 and 12  *
c     are not necessarily mutual neighbours.                           *
c                                                                      *
c***********************************************************************

#ifdef MPIALL
#define IMPLICIT NONE
#include "mpi_sup.h"
#else
      IMPLICIT NONE
#endif

      INCLUDE 'idim'
c 
c--MPI headers
c
#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif
c
c--I/O variables
c
      INTEGER   nlst_in, nlst_end, npart, ntot, ilenx, ilenv, ilendr
      REAL      xyzmh(5,ilenx), vxyzu(4,ilenv), drgrain(ilendr)
c
c--Common variables
c
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/dustimplicit'
      INCLUDE 'COMMONS/current'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/grain'
      INCLUDE 'COMMONS/compv'
      INCLUDE 'COMMONS/sort'
c
c--Local variables
c
      INTEGER   ipart, j, k, n, ioffset, index, index1, iparttot2, 
     &          background
      REAL      dr2, drx, dry, drz, m_b, h_b, h2_b, h3_b, u2_b,
     &          wk_sym, wk_sym2, dwdx, dxx, rho_planetesimal, 
     &          ddvx, ddvy, ddvz, scalar, term, vterm, ratio
c
c--If we use a polydisperse model, we want to calculate dust growth 
c     using background type 12 particles that are within the background
c     smoothing length of growing dust particles (types 11).  
c     Cannot simply loop over growing dust particles to look at its
c     neighbours, because this will not include *active* background
c     dust particles for which h_background > h_growing.
c
c--If we use a monodisperse model, we will only have a single type of
c     particle (type 11) and so we will be able to calculate growth
c     using a single loop. To do that, we will simply use a poly/mono 
c     flag a use type 11 as "background" particles in the first loop, 
c     and skip the second.
c
      IF (igrowthmodel.EQ.1) THEN
        background = 11
      ELSEIF (igrowthmodel.EQ.2) THEN
        background = 12
      ELSE
         WRITE (*,*) 'ERROR - dust growth model ',igrowthmodel,
     &        ' is not defined'
         CALL quit(0)
      ENDIF
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,ivar,iphase,ijvar,xyzmh,vxyzu)
C$OMP& shared(dumrho,radkernel,ddvtable,dvtable)
C$OMP& shared(wij,cnormk,drgrain,vfrag)
C$OMP& shared(doublehumpnormk,npart,ntot,iscurrent)
C$OMP& shared(udens,vturb,term,igrowthmodel,background)
C$OMP& shared(vcomp,iunique,iorig)
C$OMP& private(j,k,n,ipart,ioffset,drx,dry,drz,dr2)
C$OMP& private(ddvx,ddvy,ddvz,scalar)
C$OMP& private(m_b,h_b,h2_b,h3_b,u2_b)
C$OMP& private(wk_sym,wk_sym2)
C$OMP& private(dwdx,dxx,index,index1)
C$OMP& private(iparttot2,rho_planetesimal)
C$OMP& private(vterm,ratio)
c
c--Normalization calculation from looping over dust particles
c     (and zeroes the values initially too)
c
C$OMP DO SCHEDULE(runtime)
      DO 20 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
c
c--Loop only through growing dust particles (types 11)
c
         IF (iphase(ipart).NE.11) GOTO 20
c
c--Check for non-MPI code that the particle number makes sense.
c     For MPI code, need to define correct index to the xyzmh() array,
c     which is different to most arrays because the tree structure is
c     stored in the array.
c
         iparttot2 = ipart
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#else
         IF (ipart.GT.npart) THEN
            WRITE (*,*) 'ERROR - ipart.GT.npart'
            CALL quit(0)
         ENDIF
#endif
c
c--Loop over neighbouring particles. They will be type 11 if we are in a
c  monodisperse model, and type 12 if we are in a polydisperse one
c
         ioffset = ivar(2,n)
         DO k = 1, ivar(1,n)
c
c--Get neighbouring particle index from the list
c
            j = ijvar(ioffset + k)
c
c--Loop only through background particles. They will be type 11 in a 
c     monodisperse grow model and type 12 in a polydisperse one.
c
            IF (iphase(j).EQ.background) THEN

               IF (dumrho(j).EQ.0.0) THEN
                  WRITE (*,*) 'ERROR - dumrho(j).EQ.0 - dust_growth',j
                  CALL quit(0)
               ENDIF
c
c--Dust-dust pair distances
c
               drx = xyzmh(1,iparttot2) - xyzmh(1,j)
               dry = xyzmh(2,iparttot2) - xyzmh(2,j)
               drz = xyzmh(3,iparttot2) - xyzmh(3,j)
               dr2  = drx*drx + dry*dry + drz*drz
c
c--Relative velocities needed for dust growth
c
               IF (dr2.NE.0.0) THEN
                  ddvx = vxyzu(1,ipart) - vxyzu(1,j)
                  ddvy = vxyzu(2,ipart) - vxyzu(2,j)
                  ddvz = vxyzu(3,ipart) - vxyzu(3,j)
                  scalar = (ddvx*drx + ddvy*dry + ddvz*drz)/SQRT(dr2)
               ELSE
                  scalar = 0.
               ENDIF
c
c--Mass and smoothing length coefficients
c
               m_b      = xyzmh(4,j)
               h_b      = xyzmh(5,j)
c
c--Kernel (double hump) calculations
c
               h2_b     = h_b*h_b
               h3_b     = h2_b*h_b
               u2_b     = dr2/h2_b
               IF (u2_b.LT.radkernel**2) THEN
                  index  = u2_b*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx    = u2_b - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  IF (igrowthmodel.EQ.1) THEN
                     wk_sym = cnormk*(wij(index) + dwdx*dxx)/h3_b
                  ELSEIF (igrowthmodel.EQ.2) THEN
                     wk_sym  = doublehumpnormk*(wij(index) + dwdx*dxx)*
     &                 u2_b/h3_b
                  ENDIF
c
c--Grain growth term (driven by turbulent velocity by default). 
c     If fragmentation is on, growth is reversed when the fragmentation
c     threshold is crossed
c
                  CALL rho_grain(rho_planetesimal)
                  term = vturb(j)/rho_planetesimal
                  IF (vfrag.GT.0.0 .AND. vturb(j).GT.vfrag) THEN
                     term = -term
                  ENDIF
c
c--If we are in a polydisperse model, we additionally choose between 
c     the turbulent velocity and the ballistic velocity as a function 
c     of who is larger. If fragmentation is on, growth is reversed 
c     when the fragmentation threshold is crossed as well
c
                  IF (igrowthmodel.EQ.2) THEN
                     IF (scalar.LT.0.0.AND.ABS(scalar).GT.vturb(j)) THEN
                        term = ABS(scalar)/rho_planetesimal
                        IF (vfrag.GT.0.0.AND.ABS(scalar).GT.vfrag) THEN
                           term = -term
                        ENDIF
                     ENDIF
                  ENDIF

                  drgrain(ipart) = drgrain(ipart) + m_b*term*wk_sym
c
c--Average of the turbulent velocity vs ballistic velocity ratio
c
                  IF (scalar.LT.0.0) THEN
                     vterm = ABS(scalar)/vturb(j)
                     vcomp(ipart) = vcomp(ipart)  +
     &                    (m_b/dumrho(j))*vterm*wk_sym
                  ENDIF
               ENDIF
            ENDIF
         END DO   ! End background neighbours loop
 20   END DO
C$OMP END DO

      IF (igrowthmodel.EQ.1) GOTO 40
c
c--Now, add contributions from looping over background dust particles.
c     This loop is only necessary in a polydisperse growth model.
c
C$OMP DO SCHEDULE(runtime)
      DO 30 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
c
c--Loop only through background particles
c
         IF (iphase(ipart).NE.background) GOTO 30
c
c--Check for non-MPI code that the particle number makes sense.
c     For MPI code, need to define correct index to the xyzmh() array,
c     which is different to most arrays because the tree structure is
c     stored in the array.
c
         iparttot2 = ipart
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#else
         IF (ipart.GT.npart) THEN
            WRITE (*,*) 'ERROR - ipart.GT.npart'
            CALL quit(0)
         ENDIF
#endif
         m_b = xyzmh(4,iparttot2)
         h_b = xyzmh(5,iparttot2)
c
c--Loop over neighbouring particles
c
         ioffset = ivar(2,n)
         DO 50 k = 1, ivar(1,n)
c
c--Get neighbouring particle index from the list
c
            j = ijvar(ioffset + k)
c
c--Dust pair distances
c
            drx = xyzmh(1,iparttot2) - xyzmh(1,j)
            dry = xyzmh(2,iparttot2) - xyzmh(2,j)
            drz = xyzmh(3,iparttot2) - xyzmh(3,j)
            dr2  = drx*drx + dry*dry + drz*drz
c
c--If the growing dust particle is active, only include dust particles for 
c     which h_background > h_growing (otherwise their contribution has already 
c     been accounted for).  Alternately, if the growing dust particle is 
c     inactive, need to add the contribution here regardless of how
c     h_growing and h_background compare.
c
            IF (iphase(j).EQ.11 .AND. 
     &           (dr2.GT.(radkernel*xyzmh(5,j))**2 .OR. 
     &           .NOT.iscurrent(j))) THEN
c
c--Relative velocities
c
               IF (dr2.NE.0.0) THEN
                  ddvx = vxyzu(1,ipart) - vxyzu(1,j)
                  ddvy = vxyzu(2,ipart) - vxyzu(2,j)
                  ddvz = vxyzu(3,ipart) - vxyzu(3,j)
                  scalar = (ddvx*drx + ddvy*dry + ddvz*drz)/SQRT(dr2)
               ELSE
                  scalar = 0.
               ENDIF
c
c--Kernel (double hump) calculations
c
               h2_b     = h_b*h_b
               h3_b     = h2_b*h_b
               u2_b     = dr2/h2_b
               IF (u2_b.LT.radkernel**2) THEN
                  index  = u2_b*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx    = u2_b - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  IF (igrowthmodel.EQ.1) THEN
                    wk_sym = cnormk*(wij(index) + dwdx*dxx)/h3_b
                  ELSEIF (igrowthmodel.EQ.2) THEN
                    wk_sym = doublehumpnormk*(wij(index) + dwdx*dxx)*
     &                 u2_b/h3_b
                  ENDIF  
c
c--Grain growth term (driven by turbulent velocity by default). If fragmentation
c  is on, growth is reversed when the fragmentation threshold is crossed
c
                  CALL rho_grain(rho_planetesimal)
                  term = vturb(ipart)/rho_planetesimal
                  IF (vfrag.GT.0.0 .AND. vturb(ipart).GT.vfrag) THEN
                     term = -term
                  ENDIF
c
c--If we are in a polydisperse model, we additionally choose between the turbulent velocity and the
c  ballistic velocity as a function of who is larger. If fragmentation
c  is on, growth is reversed when the fragmentation threshold is crossed as well
c
                  IF (igrowthmodel.EQ.2) THEN
                     IF (scalar.LT.0.0 .AND.
     &                    ABS(scalar).GT.vturb(ipart)) THEN
                        term = ABS(scalar)/rho_planetesimal
                        IF (vfrag.GT.0.0.AND.ABS(scalar).GT.vfrag) THEN
                           term = -term
                        ENDIF
                     ENDIF
                  ENDIF
C$OMP ATOMIC
                  drgrain(j) = drgrain(j) + m_b*term*wk_sym
               ENDIF
c
c--Average of turbulent vs ballistic velocity ratio
c
               IF (scalar.LT.0.0) THEN
                  vterm = ABS(scalar)/vturb(ipart)
C$OMP ATOMIC
                  vcomp(j) = vcomp(j)  +
     &                 (m_b/dumrho(ipart))*vterm*wk_sym
               ENDIF
            ENDIF
 50      END DO                 ! End of growing dust neighbour loop
 30   END DO
C$OMP END DO

 40   CONTINUE
C$OMP END PARALLEL

      END SUBROUTINE dust_growth
