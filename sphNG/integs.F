      SUBROUTINE integs
c************************************************************
c                                                           *
c  This subroutine integrates the systeme of differential   *
c     equations over one timestep and measures time neded.  *
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
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/ghost'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/gtdble'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/outneigh'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/neighbor_P'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/gaptbl'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/timecrit'

      CHARACTER*80 outputstring

      DIMENSION nsteplist(2,30),nsteplisttot(2,30)
      INTEGER nstepphase(2,30)
      PARAMETER (numtimestepcrit = 15)
      INTEGER nstepcriterion(numtimestepcrit,30)
      INTEGER nstepcriteriontot(numtimestepcrit,30)
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE(iprint,250)
 250  FORMAT(' entry subroutine integs')

      ifail = 0
      ioutinf = 0
      ioutsup = 0
      ioutmin = 0
      ioutmax = 0
      inmin = 10000
      inmax = -1
      inminsy = 10000
      inmaxsy = 0
      imaxit = 0
      inumit = 0
      numparticlesdone = 0
      inumfixed = 0
      inumrecalc = 0
c
c--Get time
c
      CALL getused(tbefor)
c
c--Advance system one hydro time-step
c
      dt = dtmax

      IF (itiming) CALL getused(tsteponly1)

      CALL step(dt)

      IF (itiming) THEN
         CALL getused(tsteponly2)
         tsteponly = tsteponly2 - tsteponly1
      ENDIF

      gtdouble = gtdouble + DBLE(dtmax)
      gt = REAL(gtdouble)

c
c--If modelling a planet in a section, rebuild ZEUS table 
c  for planet including gas mass

      IF (iexf.EQ.7 .AND. ibound.EQ.100) THEN
         rhillmasstemp = hillmass(planetmass(1), 0, xmass, 0)
         
         IF(gt.GT.(2.*pi).AND.rhillmasstemp.GT.(1.1*hmass)) THEN
            hmass = rhillmasstemp !Store new Hillmass
            CALL zeustable
            gapfac = 0.0
            CALL gapfactor(variation, hmass, gapfac)
         ENDIF
      ENDIF

#ifdef MPIALL
      IF (iproc.EQ.0) THEN
#endif
      IF (ioutinf.NE.0) THEN
         WRITE (iprint, *) 'h too small', ioutinf, ' times'
         WRITE (iprint, *) 'minimum no. of neigh ', inmin
      ENDIF
      WRITE (iprint, *) 'minimum no. neigh for sync step', inminsy
      WRITE (iprint, *) 'maximum no. neigh for sync step', inmaxsy
      IF (ioutsup.NE.0) THEN
         WRITE (iprint, *) 'h too big  ', ioutsup, ' times'
         WRITE (iprint, *) 'maximum no. of neigh (not hmin) ', inmax
      ENDIF
      IF (numoverflowmax.NE.0) THEN
         WRITE (iprint, *) 'numoverflowmax ', numoverflowmax
         numoverflowmax = 0
      ENDIF
      IF (ioutmin.NE.0)
     &     WRITE (iprint, *) 'h less than ', hmin, ioutmin,' times'
      IF (hmaximum.NE.0.0 .AND. ioutmax.NE.0)
     &     WRITE (iprint, *) 'h greater than ', hmaximum, ioutmax,
     &        ' times'
      IF (numparticlesdone.GT.0) THEN
         WRITE (iprint, 99200) REAL(inumit)/REAL(numparticlesdone)
         WRITE (iprint, *) ' maximum no. density iterations ',imaxit
         WRITE (iprint, 99201) REAL(inumrecalc)/REAL(numparticlesdone)
         WRITE (iprint, 99202) REAL(inumfixed)/REAL(numparticlesdone)
99200    FORMAT(' mean    no. density iterations ',F8.4)
99201    FORMAT(' mean    no. neighbour finds    ',F8.4)
99202    FORMAT(' mean    no. neighbour fixed pt ',1PE11.4)
      ENDIF
      IF (nactive-nptmass.GT.0) THEN
         WRITE (iprint, 99018) REAL(nneightot)/REAL(nactive-nptmass)
99018    FORMAT('Mean number of nneigh(i) = ',F8.4)
      ENDIF
#ifdef MPIALL
      ENDIF
#endif
c
c--Report number of particles with each timestep
c
      nsteplist(:,:) = 0
      nstepphase(:,:) = 0
      nstepcriterion(:,:) = 0

      DO j = 1, npart
         IF (iupdated(j)) THEN
            WRITE (iprint,*) 'ERROR - iupdated true on full step ',j
         ENDIF
         IF (iphase(j).GE.0) THEN
            dttime = dtmax*isteps(j)/imaxstep
            ipostime = MAX(1,MIN(30,
     &           INT((LOG10(dtmax/dttime) * 3.3219281) + 1.01) ))

            nsteplist(1,ipostime) = nsteplist(1,ipostime) + 1
            IF (iphase(j).EQ.0)
     &           nstepphase(1,ipostime) = nstepphase(1,ipostime) + 1
            IF (iphase(j).GT.0)
     &           nstepphase(2,ipostime) = nstepphase(2,ipostime) + 1
            IF (timestepcriterion(j).GT.numtimestepcrit) THEN
               WRITE (iprint,*) 'ERROR - numtimestepcrit too small'
               timestepcriterion(j) = numtimestepcrit
            ENDIF
            IF (timestepcriterion(j).GT.0) 
     &           nstepcriterion(timestepcriterion(j),ipostime) =
     &           nstepcriterion(timestepcriterion(j),ipostime) + 1
            IF (integrator.EQ.0) THEN
               IF (timestepcriterion(j).LE.0) print *,'TS ',j,
     &              timestepcriterion(j),iphase(j)
            ENDIF

            IF (idustRT.GT.0 .AND. ioptimise_column.EQ.1 
     &           .AND. iphase(j).EQ.0) THEN
               dttime = dtmax*icolumnsteps(j)/imaxstep
               ipostime = INT((LOG10(dtmax/dttime) * 3.3219281) + 1.01)
               IF ((ipostime.GE.1) .AND. (ipostime.LE.29)) THEN
                  nsteplist(2,ipostime) = nsteplist(2,ipostime)+1
               ELSE
                  nsteplist(2,30) = nsteplist(2,30) + 1
               ENDIF
            ENDIF
         ENDIF
      END DO

      DO j = 1, 30
         nsteplisttot(:,j) = nsteplist(:,j)
         nstepcriteriontot(:,j) = nstepcriterion(:,j)
      END DO

#ifdef MPI
      CALL MPI_REDUCE(nsteplist,nsteplisttot,2*30,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nstepcriterion,nstepcriteriontot,
     &     numtimestepcrit*30,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (iproc.EQ.0) THEN
         WRITE (iprint,*)
         WRITE (iprint,99100) nparttot, nactivetot
#else
         WRITE (iprint,99100) npart, nactive
#endif
c
c--Write out timestep distributions
c
99100 FORMAT ('Timestep distribution ',I10,' number active ',I10)

      IF (ioptimise_column.NE.1) THEN
         WRITE (iprint,99101)
99101    FORMAT ('Bin Timestep  All Types  Local Gas  Local Non-Gas')
      ELSE
         WRITE (iprint,99102)
99102    FORMAT ('Bin Timestep  All Types  Local Gas  Local Non-Gas',
     &     ' Optimise Col')
      ENDIF

      DO j = 0, 29
         IF (ioptimise_column.NE.1) THEN
            WRITE (iprint,99110) j, dtmax/2**j, nsteplisttot(1,j+1),
     &           nstepphase(1,j+1),nstepphase(2,j+1)
         ELSE
            WRITE (iprint,99111) j, dtmax/2**j, nsteplisttot(1,j+1),
     &           nstepphase(1,j+1),nstepphase(2,j+1),nsteplisttot(2,j+1)
         ENDIF
      END DO
99110 FORMAT (I2,1PE12.4,1X,I10,1X,I10,1X,I10)
99111 FORMAT (I2,1PE12.4,1X,I10,1X,I10,1X,I10,1X,I10)
c
c--Write out how many particles have timesteps set by which criterion
c
      WRITE (iprint,*)
      WRITE (iprint,99120)
99120 FORMAT ('Distributions of timestep criteria - numbers',
     &     ' are LOG10(n)')
      WRITE (iprint,99121)
99121 FORMAT ('Bin RK-x RK-v RK-u RK-h CR-2 CR-1 Cool Drag DuDt',
     &     ' Ohm  Hall Ambi CCol etaB 1Dust')
      DO j = 0, 29
         WRITE (outputstring,'(I2)') j
         DO k = 1, 15
            IF (nstepcriteriontot(k,j+1).GT.0) THEN
               WRITE (outputstring(4+(k-1)*5:4+k*5),'(1X,F4.1)')
     &              LOG10(REAL(nstepcriteriontot(k,j+1)))
c         WRITE (iprint,99106) j, LOG10(REAL(nstepcriteriontot(1:15,j+1)))
            ELSE
               WRITE (outputstring(4+(k-1)*5:4+k*5),*) '  -  '
            ENDIF
         END DO
         WRITE (iprint,'(A80)') outputstring
      END DO
c
c--Other MPI processes (iproc>0), write out local timestep distributions
c
#ifdef MPI
      ELSE
         WRITE (iprint,99100) npart, nactive
         DO j = 0, 29
            WRITE (iprint,99130) j, dtmax/2**j, nsteplist(1,j+1)
         END DO
      ENDIF
99130 FORMAT (I2,1PE12.4,I10)
#endif

#ifdef MPI
c
c--Load balance MPI processes
c
c      IF (iresort.EQ.0) THEN
      CALL balance(1)
c      ENDIF
#endif
c
c--Get time and compute time needed for this timestep
c
      CALL getused(tafter)
      tstep = (tafter - tbefor)
c
c--Check for messages
c
      CALL mesop
c
c--Write time, and number of ghosts
c
      nghosttot = nghost
#ifdef MPI
      CALL MPI_REDUCE(nghost,nghosttot,1,MPI_INTEGER,MPI_SUM,0,
     &   MPI_COMM_WORLD,ierr)
#endif

      WRITE (iprint, 98001) tstep/60., nghosttot
98001 FORMAT (' time used  : ', F8.3, ' min. , ghosts : ', I10)

      IF (itiming) THEN
98022    FORMAT (' step only    time: ', F8.3, ' min. ',F6.2,' %')
98008    FORMAT (' densityi     time: ', F8.3, ' min. ',F6.2,' %')
98009    FORMAT (' forcei       time: ', F8.3, ' min. ',F6.2,' %')
98016    FORMAT (' gforspt      time: ', F8.3, ' min. ',F6.2,' %')
98003    FORMAT (' insulate     time: ', F8.3, ' min. ',F6.2,' %')
98017    FORMAT ('  TRE mtree   time: ', F8.3, ' min. ',F6.2,' %')
98018    FORMAT ('  TRE treef   time: ', F8.3, ' min. ',F6.2,' %')
98019    FORMAT ('  TRE revtree time: ', F8.3, ' min. ',F6.2,' %')
98119    FORMAT ('  TRE    rt-1 time: ', F8.3, ' min. ',F6.2,' %')
98120    FORMAT ('  TRE    rt-2 time: ', F8.3, ' min. ',F6.2,' %')
98121    FORMAT ('  TRE    rt-3 time: ', F8.3, ' min. ',F6.2,' %')
98122    FORMAT ('  TRE    rt-4 time: ', F8.3, ' min. ',F6.2,' %')
98123    FORMAT ('  TRE    rt-5 time: ', F8.3, ' min. ',F6.2,' %')
98021    FORMAT (' accrete      time: ', F8.3, ' min. ',F6.2,' %')
98023    FORMAT (' step-1       time: ', F8.3, ' min. ',F6.2,' %')
98024    FORMAT (' step-2       time: ', F8.3, ' min. ',F6.2,' %')
98025    FORMAT (' step-3       time: ', F8.3, ' min. ',F6.2,' %')
98026    FORMAT (' step-4       time: ', F8.3, ' min. ',F6.2,' %')
98027    FORMAT (' step-5       time: ', F8.3, ' min. ',F6.2,' %')
98028    FORMAT (' step-6       time: ', F8.3, ' min. ',F6.2,' %')
98029    FORMAT (' step-7       time: ', F8.3, ' min. ',F6.2,' %')
98030    FORMAT (' step-8       time: ', F8.3, ' min. ',F6.2,' %')
98037    FORMAT (' step-13      time: ', F8.3, ' min. ',F6.2,' %')
98038    FORMAT (' step-14      time: ', F8.3, ' min. ',F6.2,' %')
98039    FORMAT (' step-15      time: ', F8.3, ' min. ',F6.2,' %')
98041    FORMAT (' ass rad tran time: ', F8.3, ' min. ',F6.2,' %')
98303    FORMAT (' sort         time: ', F8.3, ' min. ')
98304    FORMAT (' sort - mtree time: ', F8.3, ' min. ')
98305    FORMAT (' divBsubcycle time: ', F8.3, ' min. ',F6.2,' %')
98332    FORMAT ('  densit-1    time: ', F8.3, ' min. ',F6.2,' %')
98333    FORMAT ('  densit-2    time: ', F8.3, ' min. ',F6.2,' %')
98342    FORMAT ('  divv        time: ', F8.3, ' min. ',F6.2,' %')
         WRITE (iprint, 98020)
         WRITE (iprint, 98022) tsteponly/60., tsteponly/tstep*100.
         IF (encal.EQ.'r' .OR. encal.EQ.'m') 
     &        WRITE (iprint, 98041) tass/60.,   tass/tstep*100.
         WRITE (iprint, 98008) tdens/60.,   tdens/tstep*100.
         WRITE (iprint, 98332) td1/60.,     td1/tdens*100.
         WRITE (iprint, 98333) td2/60.,     td2/tdens*100.
         WRITE (iprint, 98342) td3/60.,     td3/tdens*100.
         WRITE (iprint, 98009) tforce/60.,  tforce/tstep*100.
         WRITE (iprint, 98016) tgforpt/60., tgforpt/tstep*100.
         WRITE (iprint, 98021) taccrete/60., taccrete/tstep*100.
         WRITE (iprint, 98305) tdivBsub/60., tdivBsub/tstep*100.
         WRITE (iprint, 98023) ts1/60.,     ts1/tstep*100.
         WRITE (iprint, 98024) ts2/60.,     ts2/tstep*100.
         WRITE (iprint, 98025) ts3/60.,     ts3/tstep*100.
         WRITE (iprint, 98026) ts4/60.,     ts4/tstep*100.
         WRITE (iprint, 98027) ts5/60.,     ts5/tstep*100.
         WRITE (iprint, 98028) ts6/60.,     ts6/tstep*100.
         WRITE (iprint, 98029) ts7/60.,     ts7/tstep*100.
         WRITE (iprint, 98030) ts8/60.,     ts8/tstep*100.
         WRITE (iprint, 98037) ts13/60.,    ts13/tstep*100.
         WRITE (iprint, 98038) ts14/60.,    ts14/tstep*100.
         WRITE (iprint, 98039) ts15/60.,    ts15/tstep*100.
         WRITE (iprint, 98003) tins/60.,    tins/tstep*100.
         IF (tins.NE.0.0) THEN
            WRITE (iprint, 98017) tmtree/60.,  tmtree/tins*100.
            WRITE (iprint, 98018) ttreef/60.,  ttreef/tins*100.
            WRITE (iprint, 98019) trevt/60.,   trevt/tins*100.
         ELSE
            WRITE (iprint, 98017) tmtree/60.,  tmtree
            WRITE (iprint, 98018) ttreef/60.,  ttreef
            WRITE (iprint, 98019) trevt/60.,   trevt
         ENDIF
         IF (trevt.NE.0.0) THEN
           WRITE (iprint, 98119) revtreep1/60., revtreep1/trevt*100.
           WRITE (iprint, 98120) revtreep2/60., revtreep2/trevt*100.
           IF (revtreep2.GT.0.0) THEN
              WRITE (iprint, 98121) revtreep3/60., revtreep3/
     &              revtreep2*100.
              WRITE (iprint, 98122) revtreep4/60., revtreep4/
     &              revtreep2*100.
              WRITE (iprint, 98123) revtreep5/60., revtreep5/
     &              revtreep2*100.
           ENDIF
         ENDIF
         WRITE (iprint, 98303) tsort/60.
         IF (tsort.NE.0.0) THEN
            WRITE (iprint, 98304) tsortmtree/60.
         ENDIF
         tdens = 0.
         tforce = 0.
         tgforpt = 0.
         taccrete = 0.
         td1 = 0.
         td2 = 0.
         td3 = 0.
         ts1 = 0.
         ts2 = 0.
         ts3 = 0.
         ts4 = 0.
         ts5 = 0.
         ts6 = 0.
         ts7 = 0.
         ts8 = 0.
         ts9 = 0.
         ts10 = 0.
         ts11 = 0.
         ts12 = 0.
         ts13 = 0.
         ts14 = 0.
         ts15 = 0.
         ts16 = 0.
         ts17 = 0.
         ts18 = 0.
         ts19 = 0.
         tass = 0.
         tsort = 0.
         tsortmtree = 0.
         tins = 0.
         tins = 0.
         tmtree = 0.
         ttreef = 0.
         trevt = 0.
         revtreep1 = 0.
         revtreep2 = 0.
         revtreep3 = 0.
         revtreep4 = 0.
         revtreep5 = 0.
         tdivBsub = 0.
98020    FORMAT (' TREE Timing:')
      ENDIF

      RETURN
      END
