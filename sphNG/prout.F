      SUBROUTINE prout(where)
c************************************************************
c                                                           *
c  This routine prints out all interesting quantities at    *
c     the present time                                      *
c                                                           *
c************************************************************
      INCLUDE 'idim'

#ifdef MPI
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/angm'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/latti'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/trans'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/out1'
      INCLUDE 'COMMONS/out2'
      INCLUDE 'COMMONS/new'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/task'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/active'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/physeos'
      INCLUDE 'COMMONS/outmhd'

#ifdef MPI
      REAL*4 romax1tot, romax2tot, valphamin1tot, valphamax1tot,
     &     valphamin2tot, valphamax2tot
#endif

      CHARACTER*24 sentenc
      CHARACTER*7 where
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine prout')
c
c--Check for output type
c
      IF (where(1:6).EQ.'inform') THEN

         tkintot = tkin
         tgravtot = tgrav
         ttermtot = tterm
         tmagtot = tmag
         tradtot = trad
         trotztot = trotz
         trotxtot = trotx
         angtotot = angto
         totmomtot = totmom
         nactivetot = nactive
         nptmasstot = nptmass

         nreassigntot = nreassign
         naccretetot = naccrete
         nkilltot = nkill

         n1tot = n1
         cmx1tot = cmx1
         cmy1tot = cmy1
         cmz1tot = cmz1
         vcmx1tot = vcmx1
         vcmy1tot = vcmy1
         vcmz1tot = vcmz1
         hmi1tot = hmi1
         hma1tot = hma1
         dmax1tot = dmax1
         zmax1tot = zmax1
         romean1tot = romean1
         romax1tot = romax1
         rocen1tot = rocen1
         valphamin1tot = valphamin1
         valphamax1tot = valphamax1
         tgmean1tot = tgmean1
         tgmax1tot = tgmax1
         tgcen1tot = tgcen1
         trmean1tot = trmean1
         trmax1tot = trmax1
         trcen1tot = trcen1

         n2tot = n2
         cmx2tot = cmx2
         cmy2tot = cmy2
         cmz2tot = cmz2
         vcmx2tot = vcmx2
         vcmy2tot = vcmy2
         vcmz2tot = vcmz2
         hmi2tot = hmi2
         hma2tot = hma2
         dmax2tot = dmax2
         zmax2tot = zmax2
         romean2tot = romean2
         romax2tot = romax2
         rocen2tot = rocen2
         valphamin2tot = valphamin2
         valphamax2tot = valphamax2
         tgmean2tot = tgmean2
         tgmax2tot = tgmax2
         tgcen2tot = tgcen2
         trmean2tot = trmean2
         trmax2tot = trmax2
         trcen2tot = trcen2

#ifdef MPI
      CALL MPI_REDUCE(tkin,tkintot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tgrav,tgravtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tterm,ttermtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmag,tmagtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(trad,tradtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(trotz,trotztot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(trotx,trotxtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(angto,angtotot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(totmom,totmomtot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(escap,escaptot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nptmass,nptmasstot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nreassign,nreassigntot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(naccrete,naccretetot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(nkill,nkilltot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)

      CALL MPI_REDUCE(n1,n1tot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(cmx1,cmx1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(cmy1,cmy1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(cmz1,cmz1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(vcmx1,vcmx1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(vcmy1,vcmy1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(vcmz1,vcmz1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(hmi1,hmi1tot,1,MPI_REAL8,MPI_MIN,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(hma1,hma1tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(dmax1,dmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(zmax1,zmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(romean1,romean1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      romean1tot = romean1tot/numproc
      CALL MPI_REDUCE(romax1,romax1tot,1,MPI_REAL4,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(rocen1,rocen1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      rocen1tot = rocen1tot/numproc
      CALL MPI_REDUCE(valphamin1,valphamin1tot,1,MPI_REAL4,MPI_MIN,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(valphamax1,valphamax1tot,1,MPI_REAL4,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tgmean1,tgmean1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      tgmean1tot = tgmean1tot/numproc
      CALL MPI_REDUCE(tgmax1,tgmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tgcen1,tgcen1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      tgcen1tot = tgcen1tot/numproc
      CALL MPI_REDUCE(trmean1,trmean1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      trmean1tot = trmean1tot/numproc
      CALL MPI_REDUCE(trmax1,trmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(trcen1,trcen1tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      trcen1tot = trcen1tot/numproc

      CALL MPI_REDUCE(n2,n2tot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(cmx2,cmx2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(cmy2,cmy2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(cmz2,cmz2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(vcmx2,vcmx2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(vcmy2,vcmy2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(vcmz2,vcmz2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(hmi2,hmi2tot,1,MPI_REAL8,MPI_MIN,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(hma2,hma2tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(dmax2,dmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(zmax2,zmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(romean2,romean2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      romean2tot = romean2tot/numproc
      CALL MPI_REDUCE(romax2,romax2tot,1,MPI_REAL4,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(rocen2,rocen2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      rocen2tot = rocen2tot/numproc
      CALL MPI_REDUCE(valphamin2,valphamin2tot,1,MPI_REAL4,MPI_MIN,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(valphamax2,valphamax2tot,1,MPI_REAL4,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tgmean2,tgmean2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      tgmean2tot = tgmean2tot/numproc
      CALL MPI_REDUCE(tgmax2,tgmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tgcen2,tgcen2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      tgcen2tot = tgcen2tot/numproc
      CALL MPI_REDUCE(trmean2,trmean2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      trmean2tot = trmean2tot/numproc
      CALL MPI_REDUCE(trmax2,trmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(trcen2,trcen2tot,1,MPI_REAL8,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      trcen2tot = trcen2tot/numproc

      IF (iproc.EQ.0) THEN
#else
      escaptot = escap
#endif

c
c--Write integration output :
c
c--Number of cycle
c
c         IF (nbuild.EQ.1 .OR. MOD(ncount, nstep).EQ.0 .OR. 
         IF (MOD(ncount, nstep).EQ.0 .OR. 
     &                                        iptcreat.EQ.1) THEN
            WRITE (iprint, 99002, ERR=100) irec, MOD(ncount, nstep)+1
         ELSE
            WRITE (iprint, 99202, ERR=100) irec, MOD(ncount, nstep)+1
         ENDIF
99002    FORMAT (//, ' ------> C Y C L E   N O  : ', I4,'/', I2,
     &           '  W R I T T E N   O N   D I S K <------', /)
99202    FORMAT (//, ' ------> C Y C L E   N O  : ', I4,'/', I2,
     &           ' <------', /)
c
c--Write present time
c
         WRITE (iprint, 99003, ERR=100) gt
99003    FORMAT (1X, 'TIME  : ', 1PE16.10)
         tcomp = SQRT((3 * pi) / (32 * rhozero))
         WRITE (iprint, 99203, ERR=100) gt/tcomp
99203    FORMAT (1X, 'Free Fall Time  : ', 1PE16.10)
c
c--Energies + total angular momentum
c
         total = tkintot + tgravtot
         IF (igrp.NE.0) total = total + ttermtot + tmagtot + tradtot
         WRITE (iprint, 99004, ERR=100) total,tkintot,trotztot,trotxtot,
     &        tgravtot,ttermtot,tmagtot,tradtot,angtotot,totmomtot
99004    FORMAT (' General properties of system : ', /,
     &           ' total energy                       : ', 1PE14.5, /,
     &           ' kinetic energy                     : ', 1PE14.5, /,
     &           ' rotational energy around z         : ', 1PE14.5, /,
     &           ' rotational energy around x         : ', 1PE14.5, /,
     &           ' potential energy                   : ', 1PE14.5, /,
     &           ' internal energy                    : ', 1PE14.5, /,
     &           ' magnetic energy                    : ', 1PE14.5, /,
     &           ' radiation energy                   : ', 1PE14.5, /,
     &           ' total angular momentum             : ', 1PE14.5, /,
     &           ' total linear momentum              : ', 1PE14.5)
         
         IF (ABS(tgravtot).GT.tiny) THEN
            alph = ttermtot / ABS(tgravtot)
            betl = trotztot / ABS(tgravtot)
            betr = trotxtot / ABS(tgravtot)
            rmag2grav = tmagtot / ABS(tgravtot)
            ajeans = 1 / alph
         ELSE
            alph = 0.
            betl = 0.
            betr = 0.
            ajeans = 0.
            rmag2grav = 0.
         ENDIF
         IF (ABS(tterm).GT.tiny) THEN
            rmag2therm = tmagtot / ABS(ttermtot)
         ELSE
            rmag2therm = 0.
         ENDIF
         WRITE (iprint, 88001, ERR=100) alph, betl, betr, ajeans
88001    FORMAT (' evolutionary energy parmeters : ', /,
     &           '                        alpha  : ', 1PE14.5, /,
     &           '       (z)       beta parallel : ', 1PE14.5, /,
     &           '       (x)       beta perpend  : ', 1PE14.5, /,
     &           '                 Jeans number  : ', 1PE14.5)
         IF (imhd.EQ.idim) THEN
            WRITE(iprint,89001,ERR=100) rmag2grav,rmag2therm
89001       FORMAT ('         magnetic/grav energy  : ', 1PE14.5, /,
     &              '             magnetic/thermal  : ', 1PE14.5)
         ENDIF
c
c--Escapors
c
         WRITE (iprint, 99005, ERR=100) escaptot
99005    FORMAT (' escapors (total mass)         : ', 1PE14.5)
c
c--Accretion
c
         WRITE (iprint, 88002, ERR=100) nactivetot, nptmasstot
88002    FORMAT (' number of active particles    : ', I8, I8)
         IF (ibound.EQ.8 .OR. ibound/10.EQ.9 .OR. ibound.EQ.100) THEN
            WRITE (iprint, 88003, ERR=100) nreassigntot
88003       FORMAT (' number of reassigned part.    : ', I8)
            WRITE (iprint, 88004, ERR=100) naccretetot
88004       FORMAT (' number of accreted part.      : ', I8)
            WRITE (iprint, 88005, ERR=100) nkilltot
88005       FORMAT (' number of killed particles    : ', I8)
         ENDIF
c         IF (iaccevol.EQ.'v' .OR. iaccevol.EQ.'s') THEN
c            WRITE (iprint, 88006, ERR=100)
c88006       FORMAT (' point mass data:')
c            DO i = 1, nptmass
c               ipt = listpm(i)
c               WRITE (iprint, 88007, ERR=100) i, iphase(ipt), 
c     &              xyzmh(5,ipt)
c88007          FORMAT ('  point mass: ',I2,' type: ',I1,
c     &                 ' hacc: ',1PE12.3)
c            END DO
c            WRITE (iprint,88008) rmax
c88008       FORMAT (' new boundary radius: ',1PE12.5)
c         ENDIF
c
c--Magnetic field parameters	(computed in mhdparams)
c
         IF (imhd.EQ.idim) THEN    
	    WRITE (iprint, 66001, ERR=100) betamhdmin,betamhdmax,
     &           betamhdav,divBmax,divBav,curlBmax,curlBav,div2curlBmax,
     &           div2curlBav,omegamhdmax,omegamhdav,omegtol,fracdivBok,
     &           valphaBmin,valphaBmax,fluxtot,crosshel
66001       FORMAT (/,' Magnetic field parameters : ', /,
     &   '          plasma beta  min : ', 1PE14.5,'  max : ', 1PE14.5,
     &   ' mean : ', 1PE14.5, /,     
     &   '                 div B max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '                curl B max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '      (div B)/(curl B) max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '              divB*h/B max : ', 1PE14.5,
     &   '  mean: ', 1PE14.5,' frac < ',1PE9.2,' : ',0PF6.2,'%',/,
     &   '   resistivity switch  min : ', 1PE14.5, '  max : ',1PE14.5,//,
     &   ' total magnetic flux  (int div B dV)   : ', 1PE14.5,/,
     &   ' total cross helicity (int v.B dV)     : ', 1PE14.5)
         ENDIF
c
c--Object no 1
c
         WRITE (iprint, 99006, ERR=100) n1tot,cmx1tot,cmy1tot,cmz1tot,
     &        vcmx1tot,vcmy1tot, vcmz1tot, hmi1tot, hma1tot, dmax1tot,
     &        zmax1tot, romean1tot, romax1tot, rocen1tot,
     &        valphamin1tot,valphamax1tot,tgmean1tot,tgmax1tot,tgcen1tot
99006    FORMAT (/, ' Object number 1 (', I8, ' particles ) : ', /,
     &           ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,
     &           '  z :', 1PE14.5, /, ' velocity cm    vx :', 1PE14.5,
     &           ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,
     &           ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' max. dist. cm   r :', 1PE14.5, '  z :', 1PE14.5, /,
     &           ' density      mean :', 1PE14.5, ' max:', 1PE14.5, 
     &           ' cen:', 1PE14.5, /,
     &           ' visc. switch  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' temperature  mean :', 1PE14.5, ' max:', 1PE14.5,
     &           ' cen:', 1PE14.5, /)
         IF (imhd.EQ.idim) THEN
            WRITE (iprint, 99206) Bmin, Bmean, Bmax
99206       FORMAT (' mag field    min  :', 1PE14.5, 'mean:', 1PE14.5,
     &              ' max:', 1PE14.5,/)         
         ENDIF
         IF (encal.EQ.'r' .OR. encal.EQ.'m') THEN
            WRITE (iprint, 99106, ERR=100)trmean1tot,trmax1tot,trcen1tot
99106       FORMAT (' temper. rad. mean :', 1PE14.5, ' max:', 
     &           1PE14.5, ' cen:', 1PE14.5, /)
         ENDIF
c
c--Object no 2
c
         IF (n2.NE.0) THEN
            WRITE (iprint, 99007, ERR=100) n2tot,cmx2tot,cmy2tot,
     &           cmz2tot,vcmx2tot,vcmy2tot,vcmz2tot,hmi2tot,hma2tot,
     &           dmax2tot,zmax2tot,romean2tot,
     &           romax2tot,rocen2tot,valphamin2tot,valphamax2tot,
     &           tgmean2, tgmax2, tgcen2
99007       FORMAT (/, ' Object number 2 (', I8, ' particles ) : ', /,
     &           ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,
     &           '  z :', 1PE14.5, /, ' velocity cm    vx :', 1PE14.5,
     &           ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,
     &           ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' max. dist. cm   r :', 1PE14.5, '  z :', 1PE14.5, /,
     &           ' density      mean :', 1PE14.5, ' max:', 1PE14.5, 
     &           ' cen:', 1PE14.5, /,
     &           ' visc. switch  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' temperature  mean :', 1PE14.5, ' max:', 1PE14.5,
     &           ' cen:', 1PE14.5, /)
            IF (encal.EQ.'r' .OR. encal.EQ.'m') THEN
               WRITE (iprint, 99106, ERR=100) trmean2, trmax2, trcen2
            ENDIF
         ENDIF
#ifdef MPI
      ENDIF
#endif
c
c--Write transfer output
c

      ELSEIF (where(1:5).EQ.'trans') THEN
         WRITE (iprint, 99008) ibegin, iend, file1, file2
99008    FORMAT (//, ' dumps no ', I4, ' to ', I4, ' copied from file ',
     &           A7, ' to file ', A7, /)
         IF (ichang.EQ.0) WRITE (iprint, 99009)
99009    FORMAT (' no change done. ')
         IF (ichang.EQ.1) WRITE (iprint, 99010) frac, energc
99010    FORMAT (' the inner ', F5.3,
     &           ' of the total number of particles',
     &           ' have had their internal energy multiplied by: ',
     &           1PE12.5)
         IF (ichang.EQ.2) THEN
            WRITE (iprint, 99010) frac, energc
            WRITE (iprint, 99011) vexpan*udist/utime, rnorm*udist
99011       FORMAT (' a homologous expansion of ', 1PE12.5, ' cm/s',
     &              ' at ', 1PE12.5, ' cm has been added')
         ENDIF
c
c--Write output for initialisation
c
      ELSEIF (where(1:6).EQ.'newrun') THEN
         IF (varsta.EQ.'entropy') sentenc = 'specific entropy'
         IF (varsta.EQ.'intener') sentenc = 'specific internal energy' 

         IF (what.EQ.'scratch' .OR. what.EQ.'s') THEN
            WRITE (iprint, 99012)
99012       FORMAT (' New initial conditions defined from scratch')
            WRITE (iprint, 99013) idist
99013       FORMAT (
     &           ' Particles distributed according to distribution no :'
     &           , I2)
            WRITE (iprint, 99014) sentenc, thermal
99014       FORMAT (' Value of ', A24, ' : ', 1PE12.5)

            IF (encal.EQ.'i') WRITE (iprint, 99100)
99100       FORMAT (' Isothermal equation of state')
            IF (encal.EQ.'t') WRITE (iprint, 99777)
99777       FORMAT (' Thermal instability equation of state')
            IF (encal.EQ.'a') WRITE (iprint, 99103)
99103       FORMAT (' Adiabatic equation of state')
            IF (encal.EQ.'c') WRITE (iprint, 99123)
99123       FORMAT (' Adiabatic equation of state with cooling')
            IF (encal.EQ.'p') WRITE (iprint, 99101) gamma
99101       FORMAT (' Polytropic equation of state gamma is ',1PE12.5)
            IF (encal.EQ.'v') WRITE (iprint, 99102)
99102       FORMAT (' Gamma variable equation of state')
            IF (encal.EQ.'x') WRITE (iprint, 99104)
99104       FORMAT (' Physical equation of state')
            IF (encal.EQ.'r') WRITE (iprint, 99105)
99105       FORMAT (' Implicit Radiation/Adiabatic equation of state')

            IF (idist.GE.1) WRITE (iprint, 99015) xlmax, xlmin, ylmax,
     &                               ylmin, zlmax, zlmin
99015       FORMAT (/,' Particles are in volume given by :', /,
     &              ' xmax : ', 1PE12.5, ' xmin : ', 1PE12.5, /,
     &              ' ymax : ', 1PE12.5, ' ymin : ', 1PE12.5, /,
     &              ' zmax : ', 1PE12.5, ' zmin : ', 1PE12.5)
         ENDIF

         IF (what(1:5).EQ.'exist' .OR. what.EQ.'e') THEN
            WRITE (iprint, 99016) file1, file2
99016       FORMAT (//, ' new initial conditions made from file ', A7,
     &              ' and ', A7)
            WRITE (iprint, 99017) n1, xx1, yy1, zz1, vvx1, vvy1, vvz1
99017       FORMAT (/, ' first   object is made of ', I4,
     &              ' particles : ', /, '  x : ', 1PE12.5, '   y : ',
     &              1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,
     &              '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
            WRITE (iprint, 99018) n2, xx2, yy2, zz2, vvx2, vvy2, vvz2
99018       FORMAT (/, ' second  object is made of ', I4,
     &              ' particles : ', /, '  x : ', 1PE12.5, '   y : ',
     &              1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,
     &              '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
            WRITE (iprint, 99019) angto
99019       FORMAT (/, ' total angular momentum of the system : ',
     &              1PE12.5)
         ENDIF
      ENDIF
      GOTO 200
c
c--Handle errors during writing
c
 100  where = 'prout'
      CALL error(where, 1)

 200  CALL FLUSH (iprint)
      RETURN
      END
