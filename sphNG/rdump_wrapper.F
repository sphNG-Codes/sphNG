      SUBROUTINE rdump_wrapper(idisk1, ichkl, iskip)
c************************************************************
c                                                           *
c  This routine reads a dump into memory                    *
c                                                           *
c************************************************************


#ifdef MPIALL
#include "mpi_sup.h"
#endif
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/gtdble'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
c      INCLUDE 'COMMONS/torq'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/tokamak'
      INCLUDE 'COMMONS/xforce'      
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/savernd'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/dustfluid'
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/grain'
      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09rho'
      INCLUDE 'COMMONS/f2HY09'
      INCLUDE 'COMMONS/HY09check'

      DIMENSION itempsort(idim)
c      EQUIVALENCE (itempsort, next1)

      CHARACTER*7 where

      DATA icall/2/
      DATA where/'rdump_w'/
c     
c--Read
c
      IF (itrace.EQ.'all') WRITE (*, 99001)
99001 FORMAT (' entry subroutine rdump_wrapper')

      CALL rdump(idisk1, ichkl, iskip)
c
c--Make sink particle accretion radii consistent with i-file
c
      DO i = 1, nptmass
         write (*,*) 'Setting ',i,listpm(i),xyzmh(5,listpm(i)),hacc,
     &        iphase(listpm(i))
         xyzmh(5,listpm(i)) = hacc
      END DO
c
c--For initial conditions for protoplanet accretion
c
      IF (gt.EQ.0.0.AND.ibound.EQ.100) THEN
         DO i = 1, npart
            xyzmh(4,i) = xyzmh(4,i)/gapfac
         END DO
      ENDIF

      gtdouble = DBLE(gt)
c
c--Sort particles to ensure most efficient running.  Note that this 
c     should not be visible to the outside observer.  In other words,
c     an array must be kept of original list of particles and this
c     must be used to index *ANY* value from an array which is written 
c     to the outside.  This requires modification to almost every output
c     line in the code.  Done 21 Nov 2000.
c
      IF (itiming) CALL getused(tsort1)
c
c--Avoid sorting (or not)
c
      xminimum = 1.0E+30
      xmaximum = -1.0E+30
      DO i = 1, npart
         xmaximum = MAX(xmaximum, xyzmh(1,i))
         xminimum = MIN(xminimum, xyzmh(1,i))
      END DO
      xrange = xmaximum-xminimum
      istepmin = imax
      istepmax = 0
      DO i = 1, npart
         llist(i) = i
         IF (iphase(i).EQ.-1) THEN
            tempsort(i) = LOG(REAL(imax))/LOG(2.0)
         ELSE
            IF (gt.EQ.0.0 .OR. isteps(i).EQ.0) THEN
               tempsort(i) = (xyzmh(1,i)-xminimum)/xrange
            ELSE
               tempsort(i) = LOG(REAL(isteps(i)))/LOG(2.0) +
     &           (xyzmh(1,i)-xminimum)/xrange
            ENDIF
            istepmin = MIN(istepmin, isteps(i))
            istepmax = MAX(istepmax, isteps(i))
         ENDIF
      END DO
#ifdef MPI
      CALL MPI_ALLREDUCE(istepmin, istepmintot, 1, MPI_INTEGER, MPI_MIN,
     &     MPI_COMM_WORLD,ierr)
      istepmin = istepmintot
      CALL MPI_ALLREDUCE(istepmax, istepmaxtot, 1, MPI_INTEGER, MPI_MAX,
     &     MPI_COMM_WORLD,ierr)
      istepmax = istepmaxtot
#endif
c
c--Initialise timesteps to be consistent with integration method
c
      IF (individualtimesteps.EQ.0) THEN
         DO i = 1, npart
            isteps(i) = istepmin
         END DO
      ELSEIF (individualtimesteps.EQ.1) THEN
         DO i = 1, npart
            IF (iphase(i).GE.1 .AND. iphase(i).LT.10 .AND. 
     &           isteps(i).GT.istepmin) 
     &           isteps(i) = istepmin
         END DO
      ENDIF
c
c--Sort particles based on their individual timesteps and x
c

      CALL indexx(npart, llist, tempsort, iorig)

      DO i = 1, npart
         isort(iorig(i)) = i
      ENDDO
      DO k = 1, 5
         DO i = 1, npart
            tempsort(i) = xyzmh(k,iorig(i))
         END DO
         DO i = 1, npart
            xyzmh(k,i) = tempsort(i)
         END DO
      ENDDO
      DO k = 1, 4
         DO i = 1, npart
            tempsort(i) = vxyzu(k,iorig(i))
         END DO
         DO i = 1, npart
            vxyzu(k,i) = tempsort(i)
         END DO
      ENDDO
      IF (idim_grow.EQ.idim) THEN
         DO i = 1, npart
            tempsort(i) = rgrain(iorig(i))
         END DO
         DO i = 1, npart
            rgrain(i) = tempsort(i)
         END DO
      ENDIF
      DO i = 1, npart
         tempsort(i) = rho(iorig(i))
      END DO
      DO i = 1, npart
         rho(i) = tempsort(i)
      END DO
      DO i = 1, npart
         tempsort(i) = dgrav(iorig(i))
      END DO
      DO i = 1, npart
         dgrav(i) = tempsort(i)
      END DO
      IF (iener.EQ.4) THEN
         DO i = 1, npart
            tempsort(i) = h2ratio(iorig(i))
         END DO
         DO i = 1, npart
            h2ratio(i) = tempsort(i)
         END DO
         DO i = 1, npart
            tempsort(i) = abHIq(iorig(i))
         END DO
         DO i = 1, npart
            abHIq(i) = tempsort(i)
         END DO
         DO i = 1, npart
            tempsort(i) = abhpq(iorig(i))
         END DO
         DO i = 1, npart
            abhpq(i) = tempsort(i)
         END DO
         DO i = 1, npart
            tempsort(i) = abeq(iorig(i))
         END DO
         DO i = 1, npart
            abeq(i) = tempsort(i)
         END DO
         DO i = 1, npart
            tempsort(i) = abco(iorig(i))
         END DO
         DO i = 1, npart
            abco(i) = tempsort(i)
         END DO
      END IF
      DO k = 1, isizealphaMM
         DO i = 1, npart
            tempsort(i) = alphaMM(k,iorig(i))
         END DO
         DO i = 1, npart
            alphaMM(k,i) = tempsort(i)
         END DO
      ENDDO
      
      IF (encal.EQ.'r') THEN
         DO k = 1, 5
            DO i = 1, npart
               tempsort(i) = ekcle(k,iorig(i))
            END DO
            DO i = 1, npart
               ekcle(k,i) = tempsort(i)
            END DO
         END DO
         IF (idustRT.GT.0) THEN
            DO k = 1, 2
               DO i = 1, npart
                  tempsort(i) = dust_tk(k,iorig(i))
               END DO
               DO i = 1, npart
                  dust_tk(k,i) = tempsort(i)
               END DO
            END DO
            DO k = 1, nchemistry
               DO i = 1, npart
                  tempsort(i) = chemistry(k,iorig(i))
               END DO
               DO i = 1, npart
                  chemistry(k,i) = tempsort(i)
               END DO
            END DO
            DO k = 1, nheatingISR
               DO i = 1, npart
                  tempsort(i) = heatingISR(k,iorig(i))
               END DO
               DO i = 1, npart
                  heatingISR(k,i) = tempsort(i)
               END DO
            END DO
            DO i = 1, npart
               tempsort(i) = h2frac(iorig(i))
            END DO
            DO i = 1, npart
               h2frac(i) = tempsort(i)
            END DO
            DO k = 1, 2
               DO i = 1, npart
                  tempsort(i) = dh2dt(k,iorig(i))
               END DO
               DO i = 1, npart
                  dh2dt(k,i) = tempsort(i)
               END DO
            END DO
         ENDIF
      ENDIF
      IF (imhd.EQ.idim) THEN
         DO k = 1, imhdevol
            DO i = 1, npart
               tempsort(i) = Bevolxyz(k,iorig(i))
            END DO
            DO i = 1, npart
               Bevolxyz(k,i) = tempsort(i)
            END DO
         ENDDO
      ENDIF
      IF (idustFluid.NE.0) THEN
         DO i = 1, npart
            tempsort(i) = dustvar(iorig(i))
         END DO
         DO i = 1, npart
            dustvar(i) = tempsort(i)
         END DO
      ENDIF
      IF (idimHY09.EQ.idim) THEN
         DO i = 1, npart
            f2HY09drhodt(:,i) = HY09bin_rho(:,iorig(i))
         END DO
         DO i = 1, npart
            HY09bin_rho(:,i) = f2HY09drhodt(:,i)
         END DO
         DO i = 1, npart
            tempsort(i) = HY09frac_h_moved(iorig(i))
         END DO
         DO i = 1, npart
            HY09frac_h_moved(i) = tempsort(i)
         END DO         
      ENDIF
      IF (istellar_ionisation) THEN
         DO i = 1, npart
            tempsort(i) = HIIion(iorig(i))
         END DO
         DO i = 1, npart
            HIIion(i) = tempsort(i)
         END DO
         DO i = 1, npart
            tempsort(i) = dHIIdt(1,iorig(i))
         END DO
         DO i = 1, npart
            dHIIdt(1,i) = tempsort(i)
         END DO
      ENDIF

      DO i = 1, npart
         itempsort(i) = isteps(iorig(i))
      END DO
      DO i = 1, npart
         isteps(i) = itempsort(i)
      END DO
      DO i = 1, npart
         itempsort(i) = iphase(iorig(i))
      END DO
      DO i = 1, npart
         iphase(i) = itempsort(i)
      END DO
      DO i = 1, nptmass
         listpm(i) = isort(listpm(i))
      END DO
      DO i = 1, nlistinactive
         itempsort(i) = listinactive(i)
      END DO

      DO i = 1, nlistinactive
         listinactive(i) = isort(itempsort(i))
         IF (iphase(listinactive(i)).NE.-1) THEN
            WRITE (*,*) 'ERROR - sorting listinactive ',i,
     &           iphase(listinactive(i)),
     &           listinactive(i)
            CALL quit(0) 
         ENDIF
      END DO

      IF (itiming) THEN
         CALL getused(tsort2)
         tsort = tsort + (tsort2 - tsort1)
         tsortmtree = 0.
      ENDIF
c
c--Zero torques
c
c      DO i = 1, idim
c         torqt(i) = 0.
c         torqg(i) = 0.
c         torqp(i) = 0.
c         torqv(i) = 0.
c         torqc(i) = 0.
c      END DO
c
c--Check units in file the same as in the code!
c
      IF (udisti.LT.0.99999*udist .OR. udisti.GT.1.00001*udist) THEN
         print *,'ERROR - dump file distance units:',udisti,udist
         CALL error(where,1)
      ELSEIF (umassi.LT.0.99999*umass .OR.umassi.GT.1.00001*umass) THEN
         CALL error(where,2)
      ELSEIF (imhd.EQ.idim) THEN
         IF (umagfdi.LT.0.9999*umagfd 
     &      .OR.umagfdi.GT.1.00001*umagfd) THEN
            CALL error(where,4)
         ENDIF
      ENDIF
      IF (npart.GT.idim) THEN
         CALL error(where,3)
      ENDIF
c
c--Check that dtmax times are the same.  If not, modify isteps(i) as in mesop.f
c
ccc      GOTO 50

      IF (gt.NE.0.0 .AND. 
     &     (dtmaxdp.LT.0.9999*dtmax .OR. dtmaxdp.GT.1.0001*dtmax)) THEN
         ipower = INT(LOG10(dtmax/dtmaxdp)/LOG10(2.0))

         ifactor = 2**ABS(ipower)
         imaxnew = imaxstep/ifactor
         iminnew = 2*ifactor

         IF (ipower.LT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  isteps(j) = MIN(isteps(j), imaxnew)
                  isteps(j) = isteps(j)*ifactor
               ENDIF
            END DO
         ELSEIF (ipower.GT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  IF (isteps(j)/ifactor .LT. 2) CALL error(where, 4)
                  isteps(j) = isteps(j)/ifactor
               ENDIF
            END DO
         ENDIF
      ENDIF
c
c--Change reference frame
c
c      IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2.AND.gt.NE.0.0)) THEN
      IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2.AND.ibound.NE.100)
     &     .OR. (ifcor.GT.0.AND.ifcor.LE.2.AND.ibound.EQ.100.AND.
     &     gt.NE.0.0)) THEN
         CALL chanref(icall)
      ELSEIF (ifcor.GT.2) THEN
         ifcor = ifcor - 2
      ENDIF

#ifdef MPI
#ifdef MPIDEBUG
      xxmin = 1.0E+30
      DO i = 1, npart
         xxnew = xyzmh(1,i)**2 + xyzmh(2,i)**2 + (xyzmh(3,i)-0.0)**2
         IF (xxnew.LT.xxmin) THEN
            xxmin = xxnew
            ipos1 = i
         ENDIF
      END DO
      print *,iproc,': CLOSEST z=+0.0 ',ipos1,xxmin
      xxmin = 1.0E+30
      DO i = 1, npart
         xxnew = xyzmh(1,i)**2 + xyzmh(2,i)**2 + (xyzmh(3,i)+0.0)**2
         IF (xxnew.LT.xxmin) THEN
            xxmin = xxnew
            ipos2 = i
         ENDIF
      END DO
      print *,iproc,': CLOSEST z=-0.0 ',ipos2,xxmin
#endif
#endif

#ifdef MPIALL
      print *,iproc,' EXITED RDUMP_WRAPPER'
#else
      print *,' EXITED RDUMP_WRAPPER'
#endif

      IF (itrace.EQ.'all') WRITE (*, 99002)
99002 FORMAT (' exit subroutine rdump_wrapper')
      RETURN

      RETURN
      END
