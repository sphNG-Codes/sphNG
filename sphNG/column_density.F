      SUBROUTINE column_density(m,itime,ntot,nlistga,listga,nlistgn,
     &     listgn,xyzmh)
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

      DIMENSION listga(idim), listgn(idim), xyzmh(5,mmax2)
      LOGICAL itest
c
c--The columndensity array needs to be thread private because treef_P.F
c     which calls this subroutine is called by multiple OpenMP threads
c
      DIMENSION columndensity(nsightlinesmax,2)
      SAVE columndensity
C$OMP THREADPRIVATE(columndensity)

c
c--Stop for iptintree=2 because then sink particles are put in the tree
c     with their full mass and the column densities will include the 
c     sink particle masses
c
      IF (idustRT.NE.1) THEN
         WRITE (*,*) 'ERROR - idustRT.NE.1'
         CALL quit(0)
      ENDIF

      IF (iptintree.LT.0 .OR. iptintree.GT.1) THEN
         WRITE (*,*) 'ERROR - cannot use column_density.f when sink',
     &        ' particles are stored in the tree with their mass'
         CALL quit(0)
      ENDIF
      IF (iphase(m).NE.0) GOTO 100
c
c--Optimisation so that if column density doesn't change much from step to
c     step, don't recalculate it every particle time step.
c
      IF (ioptimise_column.EQ.1) THEN
         IF (itime.LT.icolumnnext(m)) GOTO 100
         IF (itime.GT.icolumnnext(m)) THEN
#ifdef MPI
            WRITE (*,*) 'ERROR - itime.GT.icolnext(m) ',itime,
     &           icolumnnext(m), m, ntot, iproc
#else
            WRITE (*,*) 'ERROR - itime.GT.icolnext(m) ',itime,
     &           icolumnnext(m), m, ntot
#endif
         ENDIF
      ENDIF
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
      h2columnmean = 0.

      DO k = 1, nsightlines
         columndensity(k,1) = 0.
         columndensity(k,2) = 0.
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
c--Excludes the mass of ghost particles & sink particles (skips 'leaves'
c     that have imfac=0)
c
         IF (n.LE.natom) THEN
            IF (imfac(n).EQ.0) GOTO 101
         ENDIF
         pmassn = xyzmh(4,n)
         
c
c--If one of the particles is a point mass or a star particle, then ignore
c
         IF (j.LE.nlistga) THEN
            itest = .NOT. (iphase(n).GE.1 .AND. iphase(n).LE.10)
         ELSE
            itest = .FALSE.
         ENDIF
         IF (j.GT. nlistga .OR. itest) THEN

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
            angle_node = ATAN(sizenode/rad)
            angle_node2 = angle_node**2
            angle_pixel2 = 4.0/nsightlines   ! pi's cancel
            angle_pixel = SQRT(angle_pixel2)

            DO k = 1, nsightlines
c
c--HEALpix theta and phi are converted to Cartesian unit vector coordinates
c     Angular distance between HEALpix vector and node distance is then
c     cos(angle) = dot product of normalised vectors
c
               cos_angle = difx*sightline(1,k) + dify*sightline(2,k) + 
     &              difz*sightline(3,k)
c
c--Ensure that cos_angle lies between -1 and 1 (rounding error can cause
c     trouble because of the finite precision of the HEALpix vectors)
c
               cos_angle = MIN(1.0,MAX(-1.0,cos_angle))
               centres = ACOS(cos_angle)

c
c--Calculate overlap of two circles on the 'sky' -- the area of the HEALPix
c     pixel (4*pi/nsightlines) and the area of the tree node
c
               fraction = 0.
c
c--No overlap
c
               IF (angle_node+angle_pixel.LE.centres
     &              .OR. cos_angle.LT.0) THEN
                  GOTO 99
c
c--Node entirely within pixel
c
               ELSEIF (centres+angle_node.LE.angle_pixel) THEN
                  fraction = angle_node2/angle_pixel2
c
c--Pixel entirely within node
c
               ELSEIF (centres+angle_pixel.LE.angle_node) THEN
                  fraction = 1.0
c
c--Non-trivial overlap calculation
c
               ELSE
                  centres2 = centres**2

                  CBD = (angle_node2 + centres2 -angle_pixel2)/
     &                 (2.0*angle_node*centres)
                  IF (ABS(CBD).GT.1.) THEN 
                     print *,'ERROR - column_density.F 1 ',CBD
                     CALL quit(1)
                  ENDIF
                  CBD = ACOS(CBD)
                  CAD = (angle_pixel2 + centres2 -angle_node2)/
     &                 (2.0*angle_pixel*centres)
                  IF (ABS(CAD).GT.1.) THEN 
                     print *,'ERROR - column_density.F 2 ',CBD
                     CALL quit(1)
                  ENDIF
                  CAD = ACOS(CAD)
                  value = (-centres+angle_node+angle_pixel)*
     &                 (centres+angle_node-angle_pixel)*
     &                 (centres-angle_node+angle_pixel)*
     &                 (+centres+angle_node+angle_pixel)
                  IF (value.LT.0.) THEN
                     print *,'ERROR - column_density.F 3 ',value
                     CALL quit(1)
                  ENDIF
                  fraction = (angle_node2*CBD + angle_pixel2*CAD -
     &                 0.5*SQRT(value))
     &                 / (pi*angle_pixel2)
               ENDIF

               contribution = pmassn/(pi*sizenode**2)*fraction
               columndensity(k,1) = columndensity(k,1) + contribution
               columndensity(k,2) = columndensity(k,2) + 
     &              contribution*h2frac(n)

 99            CONTINUE
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
            sightcolumns(k,1,m) = columndensity(k,1)
            sightcolumns(k,2,m) = columndensity(k,2)
         END DO
      ELSE
         DO k = 1, nsightlines
            sightcolumns(k,1,m) = sightcolumns(k,1,m) + 
     &           columndensity(k,1)
            sightcolumns(k,2,m) = sightcolumns(k,2,m) + 
     &           columndensity(k,2)
         END DO
      ENDIF
#else
      heatingISRi = 0.
      Gphotoelectrici = 0.
      exp_Av = 0.
      selfshield = 0.
      DO k = 1, nsightlines
         columndensitymean = columndensitymean + columndensity(k,1)
c
c--Based on Zucconi et al. (2001)
c
         columnnumberdensity = columndensity(k,1)*
     &        umass/udist**2 /(gmw*1.67E-24) + bounddens
         h2column = columndensity(k,2)*
     &        umass/udist**2 /(gmw*1.67E-24) + bounddens

         ext = columnnumberdensity * Qv(c/0.0000550)
c         ext = 0.

c
c--Compute the extinction of the ISR for the dust, the photoelectric
c     heating of the gas, and the value of the optical extinction (A_V)
c
         heatingISRi = heatingISRi + get_heatingISR(ext)
         Gphotoelectrici = Gphotoelectrici + get_Gphotoelectric(ext)
         expext = EXP(-ext)
         exp_Av = exp_Av + expext

         selfval = SQRT(1.0+h2column/(5.0E+14))
         selffac = 0.965/(1.0+h2column/(5.0E+14))**2 +
     &        0.035/selfval*EXP(-8.5E-4*selfval)
         selfshield = selfshield + expext**3.74 * selffac
      END DO
      columndensitymean = columndensitymean/nsightlines
c      columndens(m) = columndensitymean
      heatingISR(1,m) = heatingISRi/nsightlines
      heatingISR(2,m) = Gphotoelectrici/nsightlines
      heatingISR(3,m) = exp_Av/nsightlines
      heatingISR(4,m) = selfshield/nsightlines
c
c--Print out the progress (if required)
c
c      IF (MOD(m,100000).EQ.0) print *,m,heatingISRi,columndensitymean,
c     &     Qv(c/0.0000550)/(gmw*1.67E-24)
#endif

 100  RETURN
      END
