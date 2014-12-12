      SUBROUTINE column_density(m,ntot,nlistga,listga,nlistgn,listgn,
     &     xyzmh)
c************************************************************
c                                                           *
c  Subroutine by M.R Bate (06/11/13).  Determines an        *
c     estimate of the mean column density to a particle     *
c     by using lines of sights and tree nodes and leaves    *
c     This is similar to TREECOL by Clark et al. (2011)     *
c     except that they consider the contribution of nodes   *
c     to HEALpix areas rather than just using lines of      *
c     sight.  The method used here is described by Bate     *
c     & Keto (2015).                                        *
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
c
c--Do not need columndens(), unless want to access from elsewhere
c     (e.g. writing out from step_P.F before calculation is run).
c
c      COMMON /stuff/ columndens(idim)

c      REAL get_heatingISR
      DIMENSION listga(idim), listgn(idim), xyzmh(5,mmax2)
c
c--The columndensity array needs to be thread private because treef_P.F
c     which calls this subroutine is called by multiple OpenMP threads
c
      DIMENSION columndensity(nsightlinesmax)
      SAVE columndensity
C$OMP THREADPRIVATE(columndensity)

c
c--Stop for iptintree=2 because then sink particles are put in the tree
c     with their full mass and the column densities will include the 
c     sink particle masses
c
      IF (iptintree.LT.0 .OR. iptintree.GT.1) THEN
         WRITE (*,*) 'ERROR - cannot use column_density.f when sink',
     &        ' particles are stored in the tree with their mass'
         CALL quit(0)
      ENDIF
      IF (iphase(m).NE.0) GOTO 100
c
c--Need for MPI code
c
      IF (m.GT.ntot) THEN
         mtot2 = m + ntot + 2
         rrx = xyzmh(1,mtot2)
         rry = xyzmh(2,mtot2)
         rrz = xyzmh(3,mtot2)
      ELSE
         rrx = xyzmh(1,m)
         rry = xyzmh(2,m)
         rrz = xyzmh(3,m)
      ENDIF

      columndensitymean = 0.

      DO k = 1, nsightlines
         columndensity(k) = 0.
      END DO
c
c--Don't distinguish between atoms and nodes
c
      DO 101 j = 1, nlistga + nlistgn
         IF (j.GT.nlistga) THEN
            n = listgn(j-nlistga)
c
c--From empirical testing, 0.5*qrad(1,n) seems to work best for nodes.
c
c            sizenode = qrad(1,n) + 2.0*xyzmh(5,n)
            sizenode = 0.5*qrad(1,n)
         ELSE
c
c--For particles, use twice the smoothing length
c
            n = listga(j)
            sizenode = 2.0*xyzmh(5,n)
         ENDIF
c
c--If one of the particles is a point mass or a star particle, then ignore
c
         IF (j.GT. nlistga .OR. j.LE.nlistga .AND. .NOT. 
     &        (iphase(n).GE.1 .AND. iphase(n).LE.10)) THEN

            difx = xyzmh(1,n) - rrx
            dify = xyzmh(2,n) - rry
            difz = xyzmh(3,n) - rrz

            rr = difx**2 + dify**2 + difz**2 + tiny
            rad = SQRT(rr)
c
c--Unit vector pointing from particle to node
c
            difx = difx/rad
            dify = dify/rad
            difz = difz/rad
            angle_node = sizenode/rad

            IF (n.LE.natom) THEN
               pmassn = imfac(n)*xyzmh(4,n)
            ELSE
               pmassn = xyzmh(4,n)
            ENDIF

            DO k = 1, nsightlines
c
c--HEALpix theta and phi are converted to Cartesian unit vector coordinates
c     Angular distance between HEALpix vector and node distance is then
c     cos(angle) = dot product of normalised vectors
c
               cos_angle = difx*sightline(1,k) + dify*sightline(2,k) + 
     &              difz*sightline(3,k)

               v = 2.0*ABS(ACOS(cos_angle)/ATAN(angle_node))

               IF (v.LT.2.0 .AND. cos_angle.GT.0.0) THEN
                  columndensity(k) = columndensity(k) + 
     &                 pmassn/(pi*sizenode**2)
c
c--If we only wanted columndensity could use the expression below, but
c     with MPI we need to add the contributions from other MPI processes
c     and then calculate the extinction which is a non-linear function 
c     of the column density (i.e. exp(-tau) ).
c
c                  columndensitymean = columndensitymean + 
c     &                 pmassn/(pi*sizenode**2)
c
c--Also tried other possible expressions (for example specifically 
c     integrating the column density through a particle based on the 
c     impact parameter of the ray).  But this was not found to produce
c     a better result and is more complicated.
c
c                  impact_parameter = MIN(1000,INT(ncolumntable/4.0*
c     &                 v**2)+1)
c                  columndensity(k) = columndensity(k) + 
c     &                 pmassn/((0.5*sizenode)**2)
c     &                 *coltable(impact_parameter)
c     &                 /(2.0*pi)
c
               ENDIF
            END DO
         ENDIF
 101  CONTINUE

c
c--NOTE: If this was a straight mean, it could be done without storing 
c     individual values along each sight line.  However, the 
c     attenuation is given by an exponential and it is these that must
c     be averaged.
c
c     If not MPI, then can calculate the mean attenuation right here.
c     However, with MPI each sightline must be stored for each particle
c     because the sightlines will, in general, have contributions from
c     other MPI processes.  Unfortunately, this gives a significant
c     memory overhead.
c
#ifdef MPI
c
c--These are zeroed in densityiterate_gradh_MPI.F (treated the same as 
c     gravitational forces which also come from the tree).
c
      IF (m.GT.ntot) THEN
         DO k = 1, nsightlines
            sightcolumns(k,m) = columndensity(k)
         END DO
      ELSE
         DO k = 1, nsightlines
            sightcolumns(k,m) = sightcolumns(k,m) + columndensity(k)
         END DO
      ENDIF
#else
      heatingISRi = 0.
      Gphotoelectrici = 0.
      exp_Av = 0.
      DO k = 1, nsightlines
         columndensitymean = columndensitymean + columndensity(k)
c
c--Based on Zucconi et al. (2001)
c
         ext = columndensity(k)*umass/udist**2 *
     &     Qv(c/0.0000550)/(gmw*1.67E-24)
c         ext = 0.

c
c--Compute the extinction of the ISR for the dust, the photoelectric
c     heating of the gas, and the value of the optical extinction (A_V)
c
         heatingISRi = heatingISRi + get_heatingISR(ext)
         Gphotoelectrici = Gphotoelectrici + get_Gphotoelectric(ext)
         exp_Av = exp_Av + EXP(-ext)
      END DO
      columndensitymean = columndensitymean/nsightlines
c      columndens(m) = columndensitymean
      heatingISR(1,m) = heatingISRi/nsightlines
      heatingISR(2,m) = Gphotoelectrici/nsightlines
      heatingISR(3,m) = exp_Av/nsightlines
c
c--Print out the progress (if required)
c
c      IF (MOD(m,100000).EQ.0) print *,m,heatingISRi,columndensitymean,
c     &     Qv(c/0.0000550)/(gmw*1.67E-24)
#endif

 100  RETURN
      END
