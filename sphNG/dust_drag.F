      SUBROUTINE dust_drag(nlst_in,nlst_end,npart,ntot,
     &     xyzmh,ilenx,vxyzu,ilenv,itime,dt)
c***********************************************************************
c                                                                      *
c     This subroutine calculates the drag forces between gas and dust  *
c     using the method of Loren-Aguilar & Bate (2014), which           *
c     imposes momentum conservation + energy dissipation               *
c                                                                      *
c***********************************************************************
c
c--Common variables
c
      INCLUDE 'idim'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/compact'
      INCLUDE 'COMMONS/dumderivi'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/eosq'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/units'
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
c--Uncomment this line if you want to check momentum conservation 
c
c#define debug_drag
c
c--I/O variables
c
      INTEGER   nlst_in, nlst_end, npart, ntot, ilenx, ilenv, itime
      REAL      xyzmh(5,ilenx), vxyzu(4,ilenv), dt
c
c--Local variables
c
      INTEGER   ipart, i, j, k, n, ioffset, index, index1
      REAL      m_ipart, m_j, scalar, ddvx, ddvy, ddvz, dr, 
     &          dr2, drx, dry, drz, chi, h_j, h2_j, h3_j, 
     &          wk_j, wk_sym, wk_sym2, nu,
     &          mD_Ks, dwdx, dxx, vthermal, mG, ts
      PARAMETER (nu=3.0)
#ifdef debug_drag
      REAL      px1, py1, pz1, px0, py0, pz0
      REAL      Lx1, Ly1, Lz1, Lx0, Ly0, Lz0, L0, L1
#endif

c      print *,'Entered dust_drag ',nlst,nlst0
c
c--Store present velocities for drag interaction
c
C$OMP PARALLEL default(none)
C$OMP& shared(nlst_in,nlst_end,ivar,iphase,ijvar,xyzmh,vxyzu)
C$OMP& shared(dumrho,dt,radkernel,ddvtable,dvtable,nlst0)
C$OMP& shared(wij,dustnorm,r_planetesimals,rho_planetesimals)
C$OMP& shared(utime,udist,cnormk)
C$OMP& shared(doublehumpnormk,npart,ntot,itime,imaxstep,isteps)
C$OMP& shared(n_planetesimal_types,vxyzu_change)
C$OMP& private(i,j,k,n,ipart,ioffset,drx,dry,drz,dr,ddvx,ddvy,ddvz)
C$OMP& private(m_ipart,m_j,dr2)
C$OMP& private(h_j,h2_j,h3_j,chi,u2_j)
C$OMP& private(wk_sym,wk_sym2,wk_j,mG,ts)
C$OMP& private(dwdx,dxx,index,index1,beta,scalar,deltat)
C$OMP& private(r_planetesimal,rho_planetesimal)
C$OMP& private(vthermal,mD_Ks,radius,r2,iparttot2)
c
c--Drag calculation
c
C$OMP DO SCHEDULE(runtime)
      DO 40 n = nlst_in, nlst_end
c
c--Particle in the list
c
         ipart = ivar(3,n)
         iparttot2 = ipart
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#endif
c
c--Loop only through real dust particles
c
         IF (iphase(ipart).LT.11) GOTO 40
c
c--Check that dumrho() is defined (to avoid divide by zero)
c
         IF (dumrho(ipart).EQ.0.0) THEN
            WRITE (*,*) 'ERROR - idumrho(ipart).EQ.0.0 ',ipart
            CALL quit(0)
         ENDIF
c
c--Dust particle time-step
c
         deltat = dt*isteps(ipart)/imaxstep
c
c--Drag coefficient -- in code units (Epstein regime)
c
         IF (iphase(ipart).GE.11 .AND. 
     &        iphase(ipart)-10.LE.n_planetesimal_types) THEN
            r_planetesimal =  r_planetesimals(iphase(ipart)-10)
            rho_planetesimal =  rho_planetesimals(iphase(ipart)-10)
         ELSE
            WRITE (*,*) 'ERROR - iphase-10.GT.n_planetesimal_types'
            CALL quit(0)
         ENDIF

c         r_planetesimal =  1.0E-6
c         r_planetesimal = r_planetesimal*1.E5/udist
c         rho_planetesimal = 3.0

         vthermal = 1.0E5*utime/udist
         mD_Ks = r_planetesimal*rho_planetesimal/vthermal
c
c--Loop over neighbouring particles
c
         ioffset = ivar(2,n)
         DO k = 1,ivar(1,n)
c
c--Get gas neighbouring particle index from the list
c
            j = ijvar(ioffset + k)  
c
c--Loop only through gas type particles
c
            IF (iphase(j).EQ.0) THEN

               IF (dumrho(j).EQ.0.0) THEN
                  WRITE (*,*) 'ERROR - dumrho(j).EQ.0.0 ',j
                  CALL quit(0)
               ENDIF
c  
c--Dust-gas pair distances
c
               drx = xyzmh(1,iparttot2)-xyzmh(1,j)
               dry = xyzmh(2,iparttot2)-xyzmh(2,j)
               drz = xyzmh(3,iparttot2)-xyzmh(3,j)
               dr  = SQRT(drx*drx + dry*dry + drz*drz)
c
c--Calculate unit vectors
c
               IF (dr.NE.0.0) THEN
                  drx = drx/dr
                  dry = dry/dr
                  drz = drz/dr
               ELSE
                  drx = 0.0
                  dry = 0.0
                  drz = 0.0
               ENDIF
c
c--Velocity differences
c
               ddvx = vxyzu(1,ipart)-vxyzu(1,j)
               ddvy = vxyzu(2,ipart)-vxyzu(2,j)
               ddvz = vxyzu(3,ipart)-vxyzu(3,j)
               scalar = ddvx*drx + ddvy*dry + ddvz*drz
c
c--Mass coefficients
c
               m_ipart = xyzmh(4,iparttot2)
               m_j     = xyzmh(4,j)
c
c--Kernel (double hump) calculations
c
               dr2      = dr*dr
               h_j      = xyzmh(5,j)
               h2_j     = h_j*h_j
               h3_j     = h2_j*h_j
               u2_j     = dr2/h2_j

               IF (u2_j.LT.radkernel**2) THEN
                  index  = u2_j*ddvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dxx  = u2_j - index*dvtable
                  dwdx = (wij(index1) - wij(index))*ddvtable
                  wk_j = doublehumpnormk*(wij(index) + dwdx*dxx)*u2_j/
     &                 h3_j
                  wk_sym   = wk_j
               ELSE
                  wk_sym   = 0.0
               ENDIF
c
c--Constant drag
c
c           ts   = mD_Ks/(1.0+m_ipart/mG(ipart))
c
c--Epstein profile
c
               ts = mD_Ks/(dumrho(ipart)+dumrho(j))
c
c--Non-linear profile
c
c           ts   = mD_Ks/((dumrho(ipart)+dumrho(j))*vDG(ipart))
c
c--Chi calculation
c
               IF (ts.GT.0.0) THEN
                  beta = 1.0/(dumrho(ipart)+ dumrho(j))
                  IF (n.GT.nlst0) THEN
                     chi = beta*(1.0-dexp(-0.5*deltat/ts))
                  ELSE
                     chi = beta*(1.0-dexp(-deltat/ts))
                  ENDIF
               ELSE
                  WRITE (*,*) 'ERROR - ts.LE.0.0 ',ts,ipart
                  CALL quit(1)
               ENDIF
c
c--Careful with norm!!!
c
               IF (dustnorm(ipart).EQ.0.0) THEN
                  WRITE (*,*) 'ERROR - dustnorm(ipart).EQ.0.0 ',ipart,
     &                 dustnorm(ipart)
                  CALL quit(1)
               ENDIF
c
c--Change dust velocity
c
               vxyzu_change(1,ipart) = vxyzu_change(1,ipart) - nu*m_j*
     &              chi*scalar*drx*wk_sym/dustnorm(ipart)
               vxyzu_change(2,ipart) = vxyzu_change(2,ipart) - nu*m_j*
     &              chi*scalar*dry*wk_sym/dustnorm(ipart)
               vxyzu_change(3,ipart) = vxyzu_change(3,ipart) - nu*m_j*
     &              chi*scalar*drz*wk_sym/dustnorm(ipart)
c
c--Change gas velocity and internal energy
c
C$OMP ATOMIC
               vxyzu_change(1,j) = vxyzu_change(1,j) + nu*m_ipart*
     &              chi*scalar*drx*wk_sym/dustnorm(ipart)
C$OMP ATOMIC
               vxyzu_change(2,j) = vxyzu_change(2,j) + nu*m_ipart*
     &              chi*scalar*dry*wk_sym/dustnorm(ipart)
C$OMP ATOMIC
               vxyzu_change(3,j) = vxyzu_change(3,j) + nu*m_ipart*
     &              chi*scalar*drz*wk_sym/dustnorm(ipart)
cC$OMP ATOMIC
c               vxyzu_change(4,j)=vxyzu_change(4,j) + nu*m_ipart*
c     &              chi*scalar*scalar*wk_sym/dustnorm(ipart)

            ENDIF
         END DO                 ! End of GAS neighbors loop
 40   END DO                    ! End of DUST particles loop
C$OMP END DO 
C$OMP END PARALLEL
c
c--Check momentum after dust-gas drag interaction (NOTE: Only works for
c     global timesteps, with no ghost boundaries).
c
#ifdef debug_drag
      px1 = 0.0
      py1 = 0.0
      pz1 = 0.0
      Lx1 = 0.0
      Ly1 = 0.0
      Lz1 = 0.0
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(nlst_in,nlst_end,npart,xyzmh,vxyzu_change,ivar)
C$OMP& private(n,ipart,iparttot2)
C$OMP& reduction(+:px1)
C$OMP& reduction(+:py1)
C$OMP& reduction(+:pz1)
C$OMP& reduction(+:Lx1)
C$OMP& reduction(+:Ly1)
C$OMP& reduction(+:Lz1)
      DO 60 n = nlst_in, nlst_end
         ipart = ivar(3,n)
         iparttot2 = ipart
#ifdef MPI
         IF (ipart.GT.npart) iparttot2 = ipart + ntot + 2
#endif
         px1 = px1 + xyzmh(4,iparttot2)*vxyzu_change(1,ipart)
         py1 = py1 + xyzmh(4,iparttot2)*vxyzu_change(2,ipart)
         pz1 = pz1 + xyzmh(4,iparttot2)*vxyzu_change(3,ipart)
         Lx1 = Lx1 + xyzmh(4,iparttot2)*
     &        (xyzmh(2,iparttot2)*vxyzu_change(3,ipart) -   
     &        vxyzu_change(2,ipart)*xyzmh(3,iparttot2))
         Ly1 = Ly1 - xyzmh(4,iparttot2)*
     &        (xyzmh(1,iparttot2)*vxyzu_change(3,ipart) -   
     &        vxyzu_change(1,ipart)*xyzmh(3,iparttot2))
         Lz1 = Lz1 + xyzmh(4,iparttot2)*
     &        (xyzmh(1,iparttot2)*vxyzu_change(2,ipart) -   
     &        vxyzu_change(1,ipart)*xyzmh(2,iparttot2))
 65      CONTINUE
 60   END DO
C$OMP END PARALLEL DO
      L1 = SQRT(Lx1*Lx1 + Ly1*Ly1 + Lz1*Lz1)
      PRINT *, 'Momentum loss in dust-gas drag =',px1,py1,pz1
      PRINT *, 'Angular momentum loss in dust-gas drag =',L1
#endif

      END SUBROUTINE dust_drag
