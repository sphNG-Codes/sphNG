      SUBROUTINE column_los(m,itime,ntot,xyzmh)
c************************************************************
c                                                           *
c  Subroutine by M.R Bate (15/02/16).  Determines an        *
c     estimate of the line of sight column density to a     *
c     particle by using tree nodes and leaves.              *
c                                                           *
c     This routine was a derivative of gforsa for           *
c     calculating gravity from the tree leaves.             *
c                                                           *
c     THIS ROUTINE VECTORIZABLE.                            *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif

      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/sightlines'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/table'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/stellarradiation'

      INCLUDE 'COMMONS/dumderivi'
c
c--Do not need columndens(), unless want to access from elsewhere
c     (e.g. writing out from step_P.F before calculation is run).
c
c      COMMON /stuff/ columndens(idim)

      DIMENSION xyzmh(5,mmax2)
c
c--These arrays needs to be thread private because treef_P.F
c     which calls this subroutine is called by multiple OpenMP threads
c     
      DIMENSION listnodes(idim), indx(idim)
      DIMENSION distance2(idim), dimpact2(idim)
      DIMENSION column_to_sinks(4,iptdim)
c      DIMENSION sinksightline(3)
      SAVE listnodes, indx, distance2, dimpact2
      SAVE column_to_sinks
C$OMP THREADPRIVATE(column_to_sinks,listnodes,indx)
C$OMP THREADPRIVATE(distance2,dimpact2)

      PARAMETER (ntable=1000)
      COMMON /columnkernel/ coltable(ntable)
      COMMON /columnkerneln2/ coltable_n2(ntable), coltable_n2y2(ntable)
c     
c--Stop for iptintree=2 because then sink particles are put in the tree
c     with their full mass and the column densities will include the 
c     sink particle masses
c     
      IF (iptintree.LT.0 .OR. iptintree.GT.1) THEN
         WRITE (*,*) 'ERROR - cannot use column_los.f when sink',
     &        ' particles are stored in the tree with their mass'
         CALL quit(0)
      ENDIF

c      IF (iphase(m).LT.1 .OR. iphase(m).GE.10) THEN
c         WRITE (*,*) 'ERROR - column_los is for sink particles only'
c         CALL quit(0)
c      ENDIF
c
c--Need for MPI code
c
      IF (m.GT.ntot) THEN
         mtot2 = m + ntot + 2
         rrx = xyzmh(1,mtot2)
         rry = xyzmh(2,mtot2)
         rrz = xyzmh(3,mtot2)
         pmassm = xyzmh(4,mtot2)
         hm = xyzmh(5,mtot2)
      ELSE
         rrx = xyzmh(1,m)
         rry = xyzmh(2,m)
         rrz = xyzmh(3,m)
         pmassm = xyzmh(4,m)
         hm = xyzmh(5,m)
      ENDIF
c
c--For each massive sink particle, need to calculate direct radiation
c
      DO isink = 1, nptmasstot
         sinkx = stellar_xyzmhrti(1,isink)
         sinky = stellar_xyzmhrti(2,isink)
         sinkz = stellar_xyzmhrti(3,isink)

c         sinkh = 7.50000E-04
c         sinkh = 0.
c         sinkh = stellar_xyzmhrti(5,isink)

         column_to_sinks(1:4,isink) = 0.
c
c--Add self-contribution of particle
c
         IF (m.LE.ntot) THEN
c            contribution = 0.5*pmassm/hm**2
            contribution = 0.
            column_to_sinks(1,isink) = contribution
            column_to_sinks(2,isink) = contribution*2.0*h2frac(m)
            column_to_sinks(3,isink) = 0.
            column_to_sinks(4,isink) = 0.
         ENDIF
c
c--Find particles that interest line of sight between sink and particle
c
         CALL getneigh_los(m,sinkx,sinky,sinkz,ntot,
     &        numneigh,listnodes,dimpact2,xyzmh)

c         IF (ipt.GT.ntot .AND. numneigh.GT.0) THEN
c            print *,iproc,': foreign sink ',numneigh, sinkx, sinky,
c     &           dimpact2(1)
c         ENDIF
c
c--Sort into increasing distance
c
         DO i = 1, numneigh
            difx = xyzmh(1,listnodes(i)) - sinkx
            dify = xyzmh(2,listnodes(i)) - sinky
            difz = xyzmh(3,listnodes(i)) - sinkz

            distance2(i) = difx**2 + dify**2 + difz**2
         END DO

         CALL indexx2(numneigh,distance2,indx)
c
c--Unit vector for line of sight between particle and sink particle
c
c         sinksightline(1) = sinkx - rrx
c         sinksightline(2) = sinky - rry
c         sinksightline(3) = sinkz - rrz
c         rr = 0.
c         DO k = 1, 3
c            rr = rr + sinksightline(k)**2
c         END DO
c         radius_to_sink = SQRT(rr)
c         DO k = 1, 3
c            sinksightline(k) = sinksightline(k)/radius_to_sink
c         END DO
c
c--This adds some column because of finite sink particle size.
c     May or may not want to do this.
c
c         IF (ipt.LE.ntot) THEN
c            column_to_sinks(1,isink) = column_to_sinks(1,isink) +
c     &           sinkh*dumrho(m)
c            column_to_sinks(2,isink) = column_to_sinks(2,isink) +
c     &           sinkh*dumrho(m)
c         ENDIF
c
c--Determine which nodes and particles contribute to line of sight
c     
         DO j = 1, numneigh
            n = listnodes(indx(j))

            IF (n.GT.natom) THEN
               WRITE (*,*) 'ERROR - column_los', n, natom
               CALL quit(1)
            ENDIF
c
c--For particles, use radkernel times the smoothing length
c
            pmassn = xyzmh(4,n)
            hn = xyzmh(5,n)

            fraction = 0.
            v2 = dimpact2(indx(j))

            IF (v2.LT.radkernel**2) THEN
c
c--For line of sight column density
c
               ipos1 = MIN(ntable-1, INT(v2*(ntable/4.0))+1)
               r1=(ipos1-1)*(4.0/ntable)
               r2=r1+(4.0/ntable)
               valuemean = coltable(ipos1) + (v2-r1)*(ntable/4.0)*
     &              (coltable(ipos1+1)-coltable(ipos1))
               fraction = valuemean
c
c--For integral of number density squared
c
               quantity2 = coltable_n2(ipos1) + (v2-r1)*(ntable/4.0)*
     &              (coltable_n2(ipos1+1)-coltable_n2(ipos1))
               quantity3 = coltable_n2y2(ipos1) + (v2-r1)*(ntable/4.0)*
     &              (coltable_n2y2(ipos1+1)-coltable_n2y2(ipos1))

c               fraction = 1.

               contribution = pmassn/(hn)**2 * fraction
c
c--These quantities are for:
c     1: Total column density
c     2: Column density of H_2
c     3: Integral of density^2 (for all forms of H)
c     4: Column density excluding HII
c
               column_to_sinks(1,isink) = column_to_sinks(1,isink) + 
     &              contribution !cR    *       h2frac(n)

               column_to_sinks(2,isink) = column_to_sinks(2,isink) + 
     &              contribution*2.0*h2frac(n)

               column_to_sinks(3,isink) = column_to_sinks(3,isink) +
     &              (distance2(indx(j)) * quantity2 * 
     &              pmassn**2/(hn)**5
     &              + quantity3 * pmassn**2/(hn)**3)

               IF (istellar_ionisation) THEN
                  column_to_sinks(4,isink) = column_to_sinks(4,isink) +
     &              distance2(indx(j)) * contribution * (1.0-HIIion(n))
               ENDIF

            ENDIF

         END DO   ! numneigh
      END DO   ! nptmasstot
c
c--Finish off line of sight integrals
c
      DO isink = 1, nptmasstot
         selfshield = 0.

         columnnumberdensity = 
     &        column_to_sinks(1,isink)*
     &        umass/udist**2 /(gmw*1.67E-24)
         column_n2 = column_to_sinks(3,isink)* hfact**3 *
     &        (umass**2/udist**3) ! /(gmw*1.67E-24)**2

         ext = columnnumberdensity * Qv(c/0.0000550,Qabs550)
c         ext = 0.
c
c--Compute the extinction of the ISR for the dust, the photoelectric
c     heating of the gas, and the value of the optical extinction (A_V)
c
         IF (m.GT.ntot) THEN
            stellarrad(1,isink,m) = ext
            stellarrad(2,isink,m) = columnnumberdensity
            stellarrad(3,isink,m) = column_n2

            stellarrad(4,isink,m) = column_to_sinks(4,isink)
     &        *umass !           /(gmw*1.67E-24)/(pmassm*umass/2.0)
c     &           /(gmw*1.67E-24)/(pmassm/2.0)
         ELSE
            stellarrad(1,isink,m) = stellarrad(1,isink,m) + ext
            stellarrad(2,isink,m) = stellarrad(2,isink,m) + 
     &           columnnumberdensity
            stellarrad(3,isink,m) = stellarrad(3,isink,m) + column_n2

            stellarrad(4,isink,m) = stellarrad(4,isink,m) + 
     &           column_to_sinks(4,isink)
     &        *umass !           /(gmw*1.67E-24)/(pmassm*umass/2.0)
c     &           /(gmw*1.67E-24)/(pmassm/2.0)
         ENDIF
      END DO

      RETURN
      END
