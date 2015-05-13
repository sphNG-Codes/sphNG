      SUBROUTINE dust_norm_rhoG(nlst_in,nlst_end,npart,ntot,xyzmh,ilenx)
c***********************************************************************
c                                                                      *
c     This subroutine calculates the dustnorm() which essentially      *
c     is a measure of how well the dust particle kernel is normalised  *
c     (the values should be unity if the kernel is well sampled).      *
c     It also allows calculation of the gas density in the vicinity of *
c     each dust particle (though this is not currently used anywhere   *
c     else in the code).                                               *
c     The quantity dustnorm() is used in the implementation of         *
c     the method of Loren-Aguilar & Bate (2014) which treats drag      *
c     forces between gas and dust imposing momentum conservation +     *
c     energy dissipation.                                              *
c                                                                      *
c***********************************************************************
c
      IMPLICIT NONE
c
c--Common variables
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/dustimplicit'
c 
c--MPI headers
c
#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif
c
c--I/O variables
c
      INTEGER   nlst_in, nlst_end, npart, ntot, ilenx
      REAL      xyzmh(5,ilenx)
c
c--Local variables
c
      INTEGER   ipart, i, j, k, n, ioffset, index, index1, iparttot2
      REAL      m_j, 
     &          dr2, drx, dry, drz, h_j, h2_j, h3_j,      
     &          wk_j, wk_sym, wk_sym2, dwdx, dxx, u2_j

c      print *,'Entered dust_norm_rhoG ',nlst_end

C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,ivar,iphase,ijvar,xyzmh)
C$OMP& shared(dumrho,radkernel,ddvtable,dvtable)
C$OMP& shared(wij,dustnorm,cnormk,rhoG)
C$OMP& shared(doublehumpnormk,npart,ntot)
C$OMP& private(i,j,k,n,ipart,ioffset,drx,dry,drz,dr2)
C$OMP& private(m_j,h_j,h2_j,h3_j,u2_j)
C$OMP& private(wk_sym,wk_sym2,wk_j)
C$OMP& private(dwdx,dxx,index,index1)
C$OMP& private(iparttot2)
c
c--Normalization calculation
c
C$OMP DO SCHEDULE(runtime)
      DO 20 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
         iparttot2 = ipart
         dustnorm(ipart) = 0.0
         rhoG(ipart) = 0.0
c
c--Check for non-MPI code that the particle number makes sense.
c     For MPI code, need to define correct index to the xyzmh() array,
c     which is different to most arrays because the tree structure is
c     stored in the array.
c
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#else
         IF (ipart.GT.npart) THEN
            WRITE (*,*) 'ERROR - ipart.GT.npart'
            CALL quit(0)
         ENDIF
#endif
c
c--Loop only through real dust particles
c
         IF (iphase(ipart).LT.11) GOTO 20
c
c--Loop over neighbouring particles
c
         ioffset = ivar(2,n)
         DO k = 1, ivar(1,n)
c
c--Get gas neighbouring particle index from the list
c
            j = ijvar(ioffset + k)
c
c--Loop only through gas particles
c
            IF (iphase(j).EQ.0) THEN

               IF (dumrho(j).EQ.0.0) THEN
                  WRITE (*,*) 'ERROR - dumrho(j).EQ.0.0 ',j
                  CALL quit(0)
               ENDIF
c
c--Dust-gas pair distances
c
               drx = xyzmh(1,iparttot2) - xyzmh(1,j)
               dry = xyzmh(2,iparttot2) - xyzmh(2,j)
               drz = xyzmh(3,iparttot2) - xyzmh(3,j)
               dr2  = drx*drx + dry*dry + drz*drz
c
c--Mass coefficients
c
               m_j     = xyzmh(4,j)
c
c--Kernel (double hump) calculations
c
               h_j      = xyzmh(5,j)
               h2_j     = h_j*h_j
               h3_j     = h2_j*h_j
               u2_j     = dr2/h2_j
               IF (u2_j.LT.radkernel**2) THEN
                  index  = u2_j*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx    = u2_j - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  wk_j   = doublehumpnormk*(wij(index) + dwdx*dxx)*u2_j/
     &                 h3_j
                  wk_sym   = wk_j
                  wk_sym2  = cnormk*(wij(index) + dwdx*dxx)/h3_j
               ELSE
                  wk_sym   = 0.0
                  wk_sym2   = 0.0
               ENDIF
c
c--SPH normalization equation
c
               dustnorm(ipart) = dustnorm(ipart) +(m_j/dumrho(j))*wk_sym
c
c--Gas density at dust particle location
c
               rhoG(ipart) = rhoG(ipart) + m_j*wk_sym2
            ENDIF
         END DO   ! End neighbour loop

        IF (dustnorm(ipart).GT.1.5) THEN
           print *,'ERROR -- dustnorm seems too big ',ipart,
     &       dustnorm(ipart),iphase(ipart)
           CALL quit(1)
        ENDIF

20    END DO
C$OMP END DO
C$OMP END PARALLEL

      END SUBROUTINE dust_norm_rhoG
