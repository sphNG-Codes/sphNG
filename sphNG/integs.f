      SUBROUTINE integs
c************************************************************
c                                                           *
c  This subroutine integrates the systeme of differential   *
c     equations over one timestep and measures time neded.  *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

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
      INCLUDE 'COMMONS/secret'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/neighbor_P'

      DIMENSION nsteplist(30)

      CHARACTER*7 where

      DATA where/'integs'/

      ifail = 0
      ioutinf = 0
      ioutsup = 0
      ioutmin = 0
      ioutmax = 0
      inmin = 10000
      inmax = -1
      inminsy = 10000
      inmaxsy = 0
      IF (igrape.EQ.1) THEN
         searchfacmax = 0.
         searchfacmin = 1.0E+10
         neighovergrape = 0
         neighretry = 0
      ENDIF
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
      WRITE (iprint, *) 'numoverflowmax ', numoverflowmax
      numoverflowmax = 0
      IF (ioutmin.NE.0)
     &     WRITE (iprint, *) 'h less than ', hmin, ioutmin,' times'
      IF (hmaximum.NE.0.0 .AND. ioutmax.NE.0)
     &     WRITE (iprint, *) 'h greater than ', hmaximum, ioutmax,
     &        ' times'
c
c--If GRAPE, report number of neighbour lists overflows and retries
c
      IF (igrape.EQ.1) THEN
         WRITE (iprint, 99008) searchfacmax
         WRITE (iprint, 99009) searchfacmin
         WRITE (iprint, 99010) neighovergrape
         WRITE (iprint, 99011) neighretry
         WRITE (iprint, 99014) REAL(nneightot)/REAL(nactive-nptmass)
         WRITE (iprint, 99015) tkeep

99008    FORMAT('Maximum value of the search factor       = ',F8.4)
99009    FORMAT('Minimum value of the search factor       = ',F8.4)
99010    FORMAT('Number of GRAPE-neighbour-list overflows = ', I6)
99011    FORMAT('Number of GRAPE-neighbour-list retries   = ', I6)
99014    FORMAT('Mean number of nneigh(i) = ',F8.4)
99015    FORMAT('Keep GRAPE for ',F6.1,' seconds')
      ELSE
         WRITE (iprint, 99018) REAL(nneightot)/REAL(nactive-nptmass)

99018    FORMAT('Mean number of nneigh(i) = ',F8.4)
      ENDIF
c
c--Report number of particles with each timestep
c
      DO j = 1, 30
         nsteplist(j) = 0
      END DO

      DO j = 1, npart
         IF (iphase(j).GE.0) THEN
            dttime = dtmax*isteps(j)/imaxstep
            ipos = INT((LOG10(dtmax/dttime) * 3.3219281) + 1.01)
            IF ((ipos.GE.1) .AND. (ipos.LE.29)) THEN
               nsteplist(ipos) = nsteplist(ipos) + 1
            ELSE
               nsteplist(30) = nsteplist(30) + 1
            ENDIF
         ENDIF
      END DO

      WRITE (iprint,99100) npart, nactive
99100 FORMAT ('Timestep distribution ',I8,' number active ',I8)

      DO j = 0, 29
         WRITE (iprint,99101) j, dtmax/2**j, nsteplist(j+1)
      END DO
99101 FORMAT (I2,1PE12.4,I8)
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
c--Check for errors during integration
c
      IF (ifail.NE.0) CALL error(where, ifail)
c
c--Write time, and number of ghosts
c
      WRITE (iprint, 98001) tstep/60., nghost
98001 FORMAT (' time used  : ', F8.3, ' min. , ghosts : ', I8)

      IF (itiming) THEN
         IF (igrape.EQ.1) THEN
            WRITE (iprint, 98002)
            WRITE (iprint, 98022) tsteponly/60., tsteponly/tstep*100.
            WRITE (iprint, 98008) tdens/60.,   tdens/tstep*100.
            WRITE (iprint, 98009) tforce/60.,  tforce/tstep*100.
            WRITE (iprint, 98016) tgforpt/60., tgforpt/tstep*100.
            WRITE (iprint, 98021) taccrete/60., taccrete/tstep*100.
            WRITE (iprint, 98032) td1/60.,     td1/tstep*100.
            WRITE (iprint, 98033) td2/60.,     td2/tstep*100.
            WRITE (iprint, 98023) ts1/60.,     ts1/tstep*100.
            WRITE (iprint, 98024) ts2/60.,     ts2/tstep*100.
            WRITE (iprint, 98025) ts3/60.,     ts3/tstep*100.
            WRITE (iprint, 98026) ts4/60.,     ts4/tstep*100.
            WRITE (iprint, 98027) ts5/60.,     ts5/tstep*100.
            WRITE (iprint, 98028) ts6/60.,     ts6/tstep*100.
            WRITE (iprint, 98029) ts7/60.,     ts7/tstep*100.
            WRITE (iprint, 98030) ts8/60.,     ts8/tstep*100.
            WRITE (iprint, 98031) ts9/60.,     ts9/tstep*100.
            WRITE (iprint, 98003) tins/60.,    tins/tstep*100.
            WRITE (iprint, 98004) tginit/60.,  tginit/tins*100.
            WRITE (iprint, 98005) tgsort/60.,  tgsort/tins*100.
            WRITE (iprint, 98006) tgload/60.,  tgload/tins*100.
            WRITE (iprint, 98007) tgcall1/60., tgcall1/tins*100.
            WRITE (iprint, 98010) tggrav/60.,  tggrav/tins*100.
            WRITE (iprint, 98012) tgcall2/60., tgcall2/tins*100.
            WRITE (iprint, 98011) tgnei/60.,   tgnei/tins*100.
            WRITE (iprint, 98015) tgnmisc/60., tgnmisc/tins*100.
            WRITE (iprint, 98014) tgnstor/60., tgnstor/tins*100.
            WRITE (iprint, 98013) tgcall3/60., tgcall3/tins*100.
            WRITE (iprint,*) 'ninit = ',ninit1,ninit2,ninit3,
     &                       ninit4,ninit5
c            WRITE (iprint,*) 'ttest = ',ttest/60., ttest/tins*100.
            tdens = 0.
            tforce = 0.
            tgforpt = 0.
            taccrete = 0.
            td1 = 0.
            td2 = 0.
            ts1 = 0.
            ts2 = 0.
            ts3 = 0.
            ts4 = 0.
            ts5 = 0.
            ts6 = 0.
            ts7 = 0.
            ts8 = 0.
            ts9 = 0.
            tins = 0.
            tginit = 0.
            tgsort = 0.
            tgload = 0.
            tgcall1 = 0.
            tggrav = 0.
            tgcall2 = 0.
            tgnei = 0.
            tgnmisc = 0.
            tgnstor = 0.
            tgcall3 = 0.

            ninit1 = 0
            ninit2 = 0
            ninit3 = 0
            ninit4 = 0
            ninit5 = 0
            ttest = 0.

98002       FORMAT (' GRAPE Timing:')
98022       FORMAT (' step only    time: ', F8.3, ' min. ',F6.2,' %')
98008       FORMAT (' densityi     time: ', F8.3, ' min. ',F6.2,' %')
98009       FORMAT (' forcei       time: ', F8.3, ' min. ',F6.2,' %')
98016       FORMAT (' gforspt      time: ', F8.3, ' min. ',F6.2,' %')
98003       FORMAT (' insulate     time: ', F8.3, ' min. ',F6.2,' %')
98004       FORMAT ('  GRP init    time: ', F8.3, ' min. ',F6.2,' %')
98005       FORMAT ('  GRP sort    time: ', F8.3, ' min. ',F6.2,' %')
98006       FORMAT ('  GRP load    time: ', F8.3, ' min. ',F6.2,' %')
98007       FORMAT ('  GRP call1   time: ', F8.3, ' min. ',F6.2,' %')
98010       FORMAT ('  GRP gravity time: ', F8.3, ' min. ',F6.2,' %')
98011       FORMAT ('  GRP nei get time: ', F8.3, ' min. ',F6.2,' %')
98012       FORMAT ('  GRP call2   time: ', F8.3, ' min. ',F6.2,' %')
98013       FORMAT ('  GRP call3   time: ', F8.3, ' min. ',F6.2,' %')
98014       FORMAT ('  GRP nstore  time: ', F8.3, ' min. ',F6.2,' %')
98015       FORMAT ('  GRP nmisc   time: ', F8.3, ' min. ',F6.2,' %')
98017       FORMAT ('  TRE mtree   time: ', F8.3, ' min. ',F6.2,' %')
98018       FORMAT ('  TRE treef   time: ', F8.3, ' min. ',F6.2,' %')
98019       FORMAT ('  TRE revtree time: ', F8.3, ' min. ',F6.2,' %')
98119       FORMAT ('  TRE revtree1time: ', F8.3, ' min. ',F6.2,' %')
98120       FORMAT ('  TRE revtree2time: ', F8.3, ' min. ',F6.2,' %')
98121       FORMAT ('  TRE revtree3time: ', F8.3, ' min. ',F6.2,' %')
98122       FORMAT ('  TRE revtree4time: ', F8.3, ' min. ',F6.2,' %')
98123       FORMAT ('  TRE revtree5time: ', F8.3, ' min. ',F6.2,' %')
98021       FORMAT (' accrete      time: ', F8.3, ' min. ',F6.2,' %')
98023       FORMAT (' step-1       time: ', F8.3, ' min. ',F6.2,' %')
98024       FORMAT (' step-2       time: ', F8.3, ' min. ',F6.2,' %')
98025       FORMAT (' step-3       time: ', F8.3, ' min. ',F6.2,' %')
98026       FORMAT (' step-4       time: ', F8.3, ' min. ',F6.2,' %')
98027       FORMAT (' step-5       time: ', F8.3, ' min. ',F6.2,' %')
98028       FORMAT (' step-6       time: ', F8.3, ' min. ',F6.2,' %')
98029       FORMAT (' step-7       time: ', F8.3, ' min. ',F6.2,' %')
98030       FORMAT (' step-8       time: ', F8.3, ' min. ',F6.2,' %')
98031       FORMAT (' step-9       time: ', F8.3, ' min. ',F6.2,' %')
98032       FORMAT (' derivi-1     time: ', F8.3, ' min. ',F6.2,' %')
98033       FORMAT (' derivi-2     time: ', F8.3, ' min. ',F6.2,' %')
98034       FORMAT (' derivi-tot   time: ', F8.3, ' min. ',F6.2,' %')
98035       FORMAT (' step-kick1/2 time: ', F8.3, ' min. ',F6.2,' %')
98036       FORMAT (' timestep     time: ', F8.3, ' min. ',F6.2,' %')
98037       FORMAT (' step-13      time: ', F8.3, ' min. ',F6.2,' %')
98038       FORMAT (' step-14      time: ', F8.3, ' min. ',F6.2,' %')
98039       FORMAT (' step-15      time: ', F8.3, ' min. ',F6.2,' %')
98040       FORMAT (' step-16      time: ', F8.3, ' min. ',F6.2,' %')
         ELSE
            WRITE (iprint, 98020)
            WRITE (iprint, 98022) tsteponly/60., tsteponly/tstep*100.
            WRITE (iprint, 98008) tdens/60.,   tdens/tstep*100.
            WRITE (iprint, 98009) tforce/60.,  tforce/tstep*100.
            WRITE (iprint, 98016) tgforpt/60., tgforpt/tstep*100.
            WRITE (iprint, 98021) taccrete/60., taccrete/tstep*100.
            WRITE (iprint, 98034) ts10/60.,    ts10/tstep*100.
            WRITE (iprint, 98032) td1/60.,     td1/tstep*100.
            WRITE (iprint, 98033) td2/60.,     td2/tstep*100.
            WRITE (iprint, 98023) ts1/60.,     ts1/tstep*100.
            WRITE (iprint, 98024) ts2/60.,     ts2/tstep*100.
            WRITE (iprint, 98025) ts3/60.,     ts3/tstep*100.
            WRITE (iprint, 98026) ts4/60.,     ts4/tstep*100.
            WRITE (iprint, 98027) ts5/60.,     ts5/tstep*100.
            WRITE (iprint, 98028) ts6/60.,     ts6/tstep*100.
            WRITE (iprint, 98029) ts7/60.,     ts7/tstep*100.
            WRITE (iprint, 98030) ts8/60.,     ts8/tstep*100.
            WRITE (iprint, 98031) ts9/60.,     ts9/tstep*100.
            WRITE (iprint, 98035) ts11/60.,    ts11/tstep*100.
            WRITE (iprint, 98036) ts12/60.,    ts12/tstep*100.
            WRITE (iprint, 98037) ts13/60.,    ts13/tstep*100.
            WRITE (iprint, 98038) ts14/60.,    ts14/tstep*100.
            WRITE (iprint, 98039) ts15/60.,    ts15/tstep*100.
            WRITE (iprint, 98040) ts16/60.,    ts16/tstep*100.
            WRITE (iprint, 98003) tins/60.,    tins/tstep*100.
            WRITE (iprint, 98017) tmtree/60.,  tmtree/tins*100.
            WRITE (iprint, 98018) ttreef/60.,  ttreef/tins*100.
            WRITE (iprint, 98019) trevt/60.,   trevt/tins*100.
            IF (trevt.NE.0.0) THEN
              WRITE (iprint, 98119) revtreep1/60., revtreep1/trevt*100.
              WRITE (iprint, 98120) revtreep2/60., revtreep2/trevt*100.
              WRITE (iprint, 98121) revtreep3/60., revtreep3/
     &             revtreep2*100.
              WRITE (iprint, 98122) revtreep4/60., revtreep4/
     &             revtreep2*100.
              WRITE (iprint, 98123) revtreep5/60., revtreep5/
     &             revtreep2*100.
            ENDIF
            tdens = 0.
            tforce = 0.
            tgforpt = 0.
            taccrete = 0.
            td1 = 0.
            td2 = 0.
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

98020       FORMAT (' TREE Timing:')
         ENDIF
      ENDIF

      RETURN
      END
