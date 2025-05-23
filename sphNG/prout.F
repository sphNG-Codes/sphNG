      SUBROUTINE prout(where)
c************************************************************
c                                                           *
c  This routine prints out all interesting quantities at    *
c     the present time                                      *
c                                                           *
c************************************************************


#ifdef MPIALL
#include "mpi_sup.h"
#endif
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'
#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
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
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/accstat'
      INCLUDE 'COMMONS/rtstat'
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/HY09dustprops'
      INCLUDE 'COMMONS/HY09rho'
      INCLUDE 'COMMONS/dustimplstat'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'

      INTEGER nptot
      REAL evdatamintot(0:inumev)
      REAL evdatamaxtot(0:inumev)
      REAL evdatasumtot(0:inumev)
      REAL et
#endif

#ifdef MPI
      REAL*4 romax1tot, romax2tot, valphamin1tot, valphamax1tot,
     &     valphamin2tot, valphamax2tot, valphaBmaxtot, valphaBmintot
#endif

      CHARACTER*24 sentenc
      CHARACTER*7 where
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine prout')
#ifdef MPIDEBUG
      print *,iproc,': Entered prout'
#endif
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

         tkintot_o1  = tkin_o1
         tgravtot_o1 = tgrav_o1
         ttermtot_o1 = tterm_o1
         tmagtot_o1  = tmag_o1
         trotztot_o1 = trotz_o1
         trotxtot_o1 = trotx_o1

         nreassigntot = nreassign
         naccretetot = naccrete
         nkilltot = nkill

c         n1tot = n1
         cmx1tot = cmx1
         cmy1tot = cmy1
         cmz1tot = cmz1
         vcmx1tot = vcmx1
         vcmy1tot = vcmy1
         vcmz1tot = vcmz1
         hmi1tot = hmi1
         hma1tot = hma1
         dmax1tot = dmax1
         rmin1tot = rmin1
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
         IF (idustRT.GT.0) THEN
            tdmean1tot = tdmean1
            tdmax1tot = tdmax1
            tdmin1tot = tdmin1
            h2mass1tot = h2mass1
         ENDIF
         IF (idustFluid.NE.0) dustmasstot = dustmass1

c         n2tot = n2
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
         IF (idustRT.GT.0) THEN
            tdmean2tot = tdmean2
            tdmax2tot = tdmax2
            tdmin2tot = tdmin2
            h2mass2tot = h2mass2
         ENDIF

#ifdef MPI
         CALL MPI_REDUCE(npart,nparttot,1,MPI_INTEGER,
     &        MPI_SUM,0,MPI_COMM_WORLD,ierr)

         CALL MPI_REDUCE(tkin,tkintot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tgrav,tgravtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tterm,ttermtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tmag,tmagtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trad,tradtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trotz,trotztot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trotx,trotxtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(angx,angxtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(angy,angytot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(angz,angztot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         angtotot = SQRT(angxtot**2 + angytot**2 + angztot**2)
         CALL MPI_REDUCE(xmom,xmomtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(ymom,ymomtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(zmom,zmomtot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         totmomtot = SQRT(xmomtot**2 + ymomtot**2 + zmomtot**2)
         CALL MPI_REDUCE(escap,escaptot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(nactive,nactivetot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)

         CALL MPI_ALLREDUCE(nptmass,nptmasstot,1,MPI_INTEGER,MPI_SUM,
     &        MPI_COMM_WORLD,ierr)

         CALL MPI_REDUCE(tkin_o1,tkintot_o1,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tgrav_o1,tgravtot_o1,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tterm_o1,ttermtot_o1,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tmag_o1,tmagtot_o1,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trotz_o1,trotztot_o1,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trotx_o1,trotxtot_o1,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)

         CALL MPI_REDUCE(nreassign,nreassigntot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(naccrete,naccretetot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(nkill,nkilltot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
c
c--n1 and n2 are kept as global constant values
c
c      CALL MPI_REDUCE(n1,n1tot,1,MPI_INTEGER,MPI_SUM,0,
c     &     MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(n1new,n1newtot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(fmas1,fmas1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(cmx1*fmas1,cmx1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(cmy1*fmas1,cmy1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(cmz1*fmas1,cmz1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(vcmx1*fmas1,vcmx1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(vcmy1*fmas1,vcmy1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(vcmz1*fmas1,vcmz1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (fmas1tot.NE.0.0) THEN
            cmx1tot = cmx1tot/fmas1tot
            cmy1tot = cmy1tot/fmas1tot
            cmz1tot = cmz1tot/fmas1tot
            vcmx1tot = vcmx1tot/fmas1tot
            vcmy1tot = vcmy1tot/fmas1tot
            vcmz1tot = vcmz1tot/fmas1tot
         ENDIF
         CALL MPI_REDUCE(hmi1,hmi1tot,1,MPI_REAL8,MPI_MIN,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(hma1,hma1tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(dmax1,dmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(rmin1,rmin1tot,1,MPI_REAL8,MPI_MIN,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(zmax1,zmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(romean1*n1new,romean1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (n1newtot.NE.0) THEN
            romean1tot = romean1tot/n1newtot
         ENDIF
         CALL MPI_REDUCE(romax1,romax1tot,1,MPI_REAL4,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(rocen1,rocen1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         rocen1tot = rocen1tot/numproc
         CALL MPI_REDUCE(valphamin1,valphamin1tot,1,MPI_REAL4,MPI_MIN,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(valphamax1,valphamax1tot,1,MPI_REAL4,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tgmean1*n1new,tgmean1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (n1newtot.NE.0) THEN
            tgmean1tot = tgmean1tot/n1newtot
         ENDIF
         CALL MPI_REDUCE(tgmax1,tgmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tgcen1,tgcen1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         tgcen1tot = tgcen1tot/numproc
         CALL MPI_REDUCE(trmean1*n1new,trmean1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (n1newtot.NE.0) THEN
            trmean1tot = trmean1tot/n1newtot
         ENDIF
         CALL MPI_REDUCE(trmax1,trmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trcen1,trcen1tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         trcen1tot = trcen1tot/numproc

         IF (idustRT.GT.0) THEN
            CALL MPI_REDUCE(tdmean1*n1new,tdmean1tot,1,MPI_REAL8,
     &           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (n1newtot.NE.0) THEN
               tdmean1tot = tdmean1tot/n1newtot
            ENDIF
            CALL MPI_REDUCE(tdmax1,tdmax1tot,1,MPI_REAL8,MPI_MAX,0,
     &           MPI_COMM_WORLD,ierr)
            CALL MPI_REDUCE(tdmin1,tdmin1tot,1,MPI_REAL8,MPI_MIN,0,
     &           MPI_COMM_WORLD,ierr)
            CALL MPI_REDUCE(h2mass1,h2mass1tot,1,MPI_REAL8,
     &           MPI_SUM,0,MPI_COMM_WORLD,ierr)            
         ENDIF
         IF (idustFluid.NE.0) THEN
            CALL MPI_REDUCE(dustmass1,dustmasstot,1,MPI_REAL8,
     &           MPI_SUM,0,MPI_COMM_WORLD,ierr)
         ENDIF
c
c--n1 and n2 are kept as global constant values
c
c      CALL MPI_REDUCE(n2,n2tot,1,MPI_INTEGER,MPI_SUM,0,
c     &     MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(n2new,n2newtot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(fmas2,fmas2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(cmx2*fmas2,cmx2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(cmy2*fmas2,cmy2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(cmz2*fmas2,cmz2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(vcmx2*fmas2,vcmx2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(vcmy2*fmas2,vcmy2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(vcmz2*fmas2,vcmz2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (fmas2tot.NE.0.0) THEN
            cmx2tot = cmx2tot/fmas2tot
            cmy2tot = cmy2tot/fmas2tot
            cmz2tot = cmz2tot/fmas2tot
            vcmx2tot = vcmx2tot/fmas2tot
            vcmy2tot = vcmy2tot/fmas2tot
            vcmz2tot = vcmz2tot/fmas2tot
         ENDIF
         CALL MPI_REDUCE(hmi2,hmi2tot,1,MPI_REAL8,MPI_MIN,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(hma2,hma2tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(dmax2,dmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(zmax2,zmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(romean2*n2new,romean2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (n2newtot.NE.0) THEN
            romean2tot = romean2tot/n2newtot
         ENDIF
         CALL MPI_REDUCE(romax2,romax2tot,1,MPI_REAL4,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(rocen2,rocen2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         rocen2tot = rocen2tot/numproc
         CALL MPI_REDUCE(valphamin2,valphamin2tot,1,MPI_REAL4,MPI_MIN,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(valphamax2,valphamax2tot,1,MPI_REAL4,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tgmean2*n2new,tgmean2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (n2newtot.NE.0) THEN
            tgmean2tot = tgmean2tot/n2newtot
         ENDIF
         CALL MPI_REDUCE(tgmax2,tgmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(tgcen2,tgcen2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         tgcen2tot = tgcen2tot/numproc
         CALL MPI_REDUCE(trmean2*n2new,trmean2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         IF (n2newtot.NE.0) THEN
            trmean2tot = trmean2tot/n2newtot
         ENDIF
         CALL MPI_REDUCE(trmax2,trmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(trcen2,trcen2tot,1,MPI_REAL8,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         trcen2tot = trcen2tot/numproc

         IF (idustRT.GT.0) THEN
            CALL MPI_REDUCE(tdmean2*n2new,tdmean2tot,1,MPI_REAL8,
     &           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (n2newtot.NE.0) THEN
               tdmean2tot = tdmean2tot/n2newtot
            ENDIF
            CALL MPI_REDUCE(tdmax2,tdmax2tot,1,MPI_REAL8,MPI_MAX,0,
     &           MPI_COMM_WORLD,ierr)
            CALL MPI_REDUCE(tdmin2,tdmin2tot,1,MPI_REAL8,MPI_MIN,0,
     &           MPI_COMM_WORLD,ierr)
            CALL MPI_REDUCE(h2mass2,h2mass2tot,1,MPI_REAL8,
     &           MPI_SUM,0,MPI_COMM_WORLD,ierr)            
         ENDIF

         IF (imhd.EQ.idim) THEN
            CALL MPI_REDUCE(betamhdmin,betamhdmintot,1,MPI_REAL8,
     &           MPI_MIN,0,MPI_COMM_WORLD,ierr)
            betamhdmin = betamhdmintot
            CALL MPI_REDUCE(betamhdmax,betamhdmaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            betamhdmax = betamhdmaxtot
            CALL MPI_REDUCE(betamhdav*npart,betamhdavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) betamhdav = betamhdavtot/nparttot
            
            CALL MPI_REDUCE(divBmax,divBmaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            divBmax = divBmaxtot
            CALL MPI_REDUCE(divBav*npart,divBavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) divBav = divBavtot/nparttot
            
            CALL MPI_REDUCE(curlBmax,curlBmaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            curlBmax = curlBmaxtot
            CALL MPI_REDUCE(curlBav*npart,curlBavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) curlBav = curlBavtot/nparttot
            
            CALL MPI_REDUCE(div2curlBmax,div2curlBmaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            div2curlBmax = div2curlBmaxtot
            CALL MPI_REDUCE(div2curlBav*npart,div2curlBavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) div2curlBav = div2curlBavtot/nparttot

            CALL MPI_REDUCE(omegamhdmax,omegamhdmaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            omegamhdmax = omegamhdmaxtot
            CALL MPI_REDUCE(omegamhdav*npart,omegamhdavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) omegamhdav = omegamhdavtot/nparttot

            CALL MPI_REDUCE(fracdivBok*npart,fracdivBoktot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) fracdivBok = fracdivBoktot/nparttot
            
            CALL MPI_REDUCE(valphaBmin,valphaBmintot,1,MPI_REAL4,
     &           MPI_MIN,0,MPI_COMM_WORLD,ierr)
            valphaBmin = valphaBmintot
            CALL MPI_REDUCE(valphaBmax,valphaBmaxtot,1,MPI_REAL4,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            valphaBmax = valphaBmaxtot

            CALL MPI_REDUCE(Bmin,Bmintot,1,MPI_REAL8,
     &           MPI_MIN,0,MPI_COMM_WORLD,ierr)
            Bmin = Bmintot
            CALL MPI_REDUCE(Bmax,Bmaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            Bmax = Bmaxtot
            CALL MPI_REDUCE(Bmean*npart,Bmeantot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) Bmean = Bmeantot/nparttot
            
            CALL MPI_REDUCE(etamin,etamintot,1,MPI_REAL8,
     &           MPI_MIN,0,MPI_COMM_WORLD,ierr)
            etamin = etamintot
            CALL MPI_REDUCE(etamax,etamaxtot,1,MPI_REAL8,
     &           MPI_MAX,0,MPI_COMM_WORLD,ierr)
            etamax = etamaxtot
            CALL MPI_REDUCE(etaav*npart,etaavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) etaav = etaavtot/nparttot
            
            CALL MPI_REDUCE(etarealoverartmax,etarealoverartmaxtot,
     &           1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
            etarealoverartmax = etarealoverartmaxtot
            CALL MPI_REDUCE(etarealoverartmin,etarealoverartmintot,
     &           1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
            etarealoverartmin = etarealoverartmintot
            CALL MPI_REDUCE(etarealoverartav*npart,etarealavtot,1,
     &           MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            IF (nparttot.NE.0) etarealoverartav = etarealavtot/nparttot

            CALL MPI_REDUCE(fluxtot,fluxtot_tot,1,MPI_REAL8,MPI_SUM,
     &           0,MPI_COMM_WORLD,ierr)
            fluxtot = fluxtot_tot
            CALL MPI_REDUCE(crosshel,crossheltot,1,MPI_REAL8,MPI_SUM,
     &           0,MPI_COMM_WORLD,ierr)
            crosshel = crossheltot
         ENDIF
#ifdef NONIDEAL
         CALL MPI_REDUCE(np_nimhd,nptot,1,MPI_INTEGER,MPI_SUM,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(et_nimhd,et,1,MPI_REAL8,MPI_MAX,0,
     &        MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(evdata(1:ielements),evdatamintot(1:ielements),
     &        ielements,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(evdata(1:ielements),evdatamaxtot(1:ielements),
     &        ielements,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(evdata(1:ielements),evdatasumtot(1:ielements),
     &        ielements,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif

         IF (iproc.EQ.0) THEN
#else
         escaptot = escap
#ifdef NONIDEAL
         et           = et_nimhd
         nptot        = np_nimhd
         evdatamintot(1:ielements) = evdata(1:ielements)
         evdatamaxtot(1:ielements) = evdata(1:ielements)
         evdatasumtot(1:ielements) = evdata(1:ielements)
#endif
#endif
c
c--Write integration output :
c
c--Number of cycle
c
c         IF (nbuild.EQ.1 .OR. MOD(ncount, nstep).EQ.0 .OR. 
         IF (MOD(ncount, nstep).EQ.0) THEN
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
99003    FORMAT (1X, 'TIME  : ', 1PE17.10)
         tcomp = SQRT((3 * pi) / (32 * rhozero))
         WRITE (iprint, 99203, ERR=100) gt/tcomp
99203    FORMAT (1X, 'Free Fall Time  : ', 1PE17.10)
c
c--Energies + total angular momentum
c
         total = tkintot + tgravtot + ttermtot + tmagtot + tradtot
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
            IF (alph.GT.tiny) THEN
               ajeans = 1. / alph
            ELSE
               ajeans = 0.0
            ENDIF
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
88002    FORMAT (' number of active particles    : ', I10, I8)
         IF (ibound.EQ.8 .OR. ibound/10.EQ.9 .OR. ibound.EQ.100) THEN
            WRITE (iprint, 88003, ERR=100) nreassigntot
88003       FORMAT (' number of reassigned part.    : ', I8)
            WRITE (iprint, 88004, ERR=100) naccretetot
88004       FORMAT (' number of accreted part.      : ', I8)
            WRITE (iprint, 88005, ERR=100) nkilltot
88005       FORMAT (' number of killed particles    : ', I8)
         ENDIF

         IF (iexf.EQ.7) THEN
            WRITE (iprint, 99407, ERR=100) hmass
99407       FORMAT (' mass for zeus table           : ', 1PE14.5)
            WRITE (iprint, 99408, ERR=100) planetradius(1)*pradfac(1,gt)
99408       FORMAT (' radius of planet surface      : ', 1PE14.5)
            WRITE (iprint, 99409, ERR=100) gapfac
99409       FORMAT (' injection adjustment factor   : ', 1PE14.5)
         ENDIF
         
         IF (ibound.EQ.102 .OR. ibound.EQ.103) THEN
            IF (initialptm.EQ.5) THEN
               WRITE (iprint, 11407, ERR=100) nptmass
11407          FORMAT(/,' Number of planets :',I2)
               DO i = 1, nptmass
                  WRITE (iprint, 11408, ERR=100) i, xyzmh(5,listpm(i))*
     &                 pradfac(0,gt)
               END DO
11408          FORMAT('  radius of planet ',I2,' surface       : ',
     &              1PE14.5)
               DO i = 1, nptmass
                  WRITE (iprint, 99410, ERR=100) i,
     &                 hillmass(xyzmh(4,listpm(i)), 1, xmass, i)
               END DO
            ELSEIF (irotpot.EQ.1) THEN
               WRITE (iprint, 11407, ERR=100) numplanet
               DO i = 1, numplanet
                  WRITE (iprint, 11408) i,planetradius(i)*pradfac(i,gt)
               END DO
               DO i = 1, numplanet
                  WRITE (iprint, 99410) i,
     &                 hillmass(planetmass(i), 0, xmass, 0)
               END DO
            ENDIF
         ELSEIF (ibound.EQ.100 .AND. iexf.EQ.7) THEN
            WRITE (iprint, 99410, ERR=100) 1, hillmass(planetmass(1), 0,
     &           xmass, 0)
99410    FORMAT ('  mass within planet ',I2,' Hill radius : ', 1PE14.5)
         ENDIF

         IF (imigrate.EQ.1) THEN
            IF (irotpot.EQ.1) THEN
               WRITE (iprint, 11001) sqrt(px**2 + py**2),
     &              planetradius(1)
11001          FORMAT (' migrated to ', 1PE14.5, ' from ', 1PE14.5)
            ENDIF
            IF (nptmass.GE.1) THEN
               is = listpm(1)
               rmass1 = xmass
               rmass2 = xyzmh(4,is)
               rmass2 = 0.0

               rcradius = SQRT(xyzmh(1,is)**2 + xyzmh(2,is)**2)
               rmradius = SQRT(xyzmh(1,is)**2 + xyzmh(2,is)**2 +
     &              xyzmh(3,is)**2)

               amx = xyzmh(2,is)*vxyzu(3,is) - xyzmh(3,is)*vxyzu(2,is)
               amy = xyzmh(3,is)*vxyzu(1,is) - xyzmh(1,is)*vxyzu(3,is)
               amz = xyzmh(1,is)*vxyzu(2,is) - xyzmh(2,is)*vxyzu(1,is)
               
               angmom = amx + amy + amz
               angmom2 = amx**2 + amy**2 + amz**2
               
               vel2 = vxyzu(1,is)**2 + vxyzu(2,is)**2 + vxyzu(3,is)**2
               energy = vel2/2. - (rmass1+rmass2)/rmradius

               e = SQRT(1.0 + (2.*energy*angmom2/
     &              (rmass1+rmass2)**2))
               a = angmom2/((rmass1+rmass2)*(1.-e**2))

               WRITE (iprint, 11001) rcradius, planetradius(1)
            ENDIF
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
c--Magnetic field parameters      (computed in mhdparams)
c
         IF (imhd.EQ.idim) THEN    
          WRITE (iprint, 66001, ERR=100) betamhdmin,betamhdmax,
     &           betamhdav,divBmax,divBav,curlBmax,curlBav,div2curlBmax,
     &           div2curlBav,omegamhdmax,omegamhdav,omegtol,fracdivBok,
     &           valphaBmin,valphaBmax
66001       FORMAT (/,' Magnetic field parameters : ', /,
     &   '          plasma beta  min : ', 1PE14.5,'  max : ', 1PE14.5,
     &   ' mean : ', 1PE14.5, /,     
     &   '                 div B max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '                curl B max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '      (div B)/(curl B) max : ', 1PE14.5,'  mean: ', 1PE14.5,/,
     &   '              divB*h/B max : ', 1PE14.5,
     &   '  mean: ', 1PE14.5,' frac < ',1PE9.2,' : ',0PF6.2,'%',/,
     &   '   resistivity switch  min : ', 1PE14.5, '  max : ',1PE14.5)
            IF (ABS(iresist).GE.2) THEN
               WRITE (iprint, 66002, ERR=100) etarealoverartmin,
     &               etarealoverartmax,etarealoverartav
            ENDIF
66002       FORMAT (
     &   '     eta real/eta art. min : ', 1PE14.5,'  max : ', 1PE14.5,
     &   ' mean : ', 1PE14.5)

            WRITE (iprint, 66003, ERR=100) fluxtot,crosshel
66003       FORMAT (/,
     &   ' total magnetic flux  (int div B dV)   : ', 1PE14.5,/,
     &   ' total cross helicity (int v.B dV)     : ', 1PE14.5)
         ENDIF
c
c--Stellar radiation
c
         IF (istellar_radiation .AND. nptmasstot.GT.0) THEN
            WRITE (iprint, *, ERR=100)
            WRITE (iprint, *, ERR=100) 'Stellar radiation :'
            DO i = 1, nptmasstot
               WRITE (iprint, 67000, ERR=100) i,stellar_xyzmhrti(4,i),
     &              stellar_xyzmhrti(6:8,i)
67000          FORMAT (' Star: ',I4,'   M*: ',1PE8.2,' R*: ',1PE8.2,
     &              ' T*: ',1PE10.4,' Ionising flux: ',1PE8.2)
            END DO
         ENDIF
c
c--Object no 1
c
c         WRITE (96,*) gt, rmin1tot

         WRITE (iprint, 99006, ERR=100) n1,cmx1tot,cmy1tot,cmz1tot,
     &        vcmx1tot,vcmy1tot, vcmz1tot, hmi1tot, hma1tot, dmax1tot,
     &        zmax1tot, romean1tot, romax1tot, rocen1tot,
     &        valphamin1tot,valphamax1tot,tgmean1tot,tgmax1tot,tgcen1tot
99006    FORMAT (/, ' Object number 1 (', I10, ' particles ) : ', /,
     &           ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,
     &           '  z :', 1PE14.5, /, ' velocity cm    vx :', 1PE14.5,
     &           ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,
     &           ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' max. dist. cm   r :', 1PE14.5, '  z :', 1PE14.5, /,
     &           ' density      mean :', 1PE14.5, ' max:', 1PE14.5, 
     &           ' cen:', 1PE14.5, /,
     &           ' visc. switch  min :', 1PE14.5, ' max:', 1PE14.5, /,
     &           ' temperature  mean :', 1PE14.5, ' max:', 1PE14.5,
     &           ' cen:', 1PE14.5)
         IF (imhd.EQ.idim) THEN
            WRITE (iprint, 99206) Bmin, Bmean, Bmax
99206       FORMAT (' mag field    min  :', 1PE14.5, '  av:', 1PE14.5,
     &              ' max:', 1PE14.5)         
            IF (ABS(iresist).GE.2) THEN
               WRITE (iprint, 99207) etamin, etaav, etamax
99207          FORMAT (' eta resist   min  :', 1PE14.5, '  av:',1PE14.5,
     &              ' max:', 1PE14.5)
            ENDIF  
         ENDIF
         IF (encal.EQ.'r' .OR. encal.EQ.'m') THEN
            WRITE (iprint, 99106, ERR=100)trmean1tot,trmax1tot,trcen1tot
99106       FORMAT (' temper. rad. mean :', 1PE14.5, ' max:', 
     &           1PE14.5, ' cen:', 1PE14.5)
            IF (idustRT.GT.0) THEN
               WRITE (iprint, 99107, ERR=100) tdmean1tot,tdmax1tot,
     &              tdmin1tot
99107          FORMAT (' dust temper. mean :', 1PE14.5, ' max:', 
     &              1PE14.5, ' min:', 1PE14.5)
               WRITE (iprint, 99108, ERR=100) h2mass1tot
99108          FORMAT (' H_2 total mass    :', 1PE14.5)
            ENDIF
         ENDIF
c
c--Output mass distribution of dust in bins for central particle
c     Need to divide by the logarithmic bin spacing because without
c     this the dust mass fraction per bin changes as the number of bins
c     being used or as the max/min dust size range changes.
c     Note that this still includes the dust to gas ratio (i.e. this
c     gives the dust mass distribution as a fraction of the total (gas
c     plus dust) mass so the numbers will be smaller if the dust-to-gas 
c     ratio is smaller).
c
         IF (idimHY09.EQ.idim) THEN
            WRITE (iprint,99120, ERR=100) HY09_ndust_bins
            dlog_binsize = LOG(HY09_size_max/HY09_size_min)/
     &           HY09_ndust_bins
            DO i = 1, HY09_ndust_bins
               WRITE (iprint,99121, ERR=100) HY09binsizes(i)*1.0E+04,
     &              HY09cen(i), HY09cen(i)/dlog_binsize
c#ifndef MPI
c               WRITE (100+irec,*) HY09binsizes(i)*1.0E+04,
c     &              HY09cen(i), HY09cen(i)/dlog_binsize
c#endif
            END DO

99120       FORMAT (' Central dust properties (',I3,' bins):')
99121       FORMAT ('   Dust size: ',F9.4,' micron: ep_j:',1PE11.3,
     &           ' ep_j/dln(binsize):',1PE11.3)
            WRITE (iprint,99122, ERR=100) 
     &           HY09cen_totaldens/(1.0-HY09cen_totaldens),
     &           HY09_dustgas_ratio
99122       FORMAT ('   Total dust/gas ratio:   ',1PE14.5,' initially:',
     &           1PE14.5)
         ENDIF
         IF (idustFluid.NE.0) THEN
            WRITE (iprint,99124, ERR=100) dustmasstot
99124       FORMAT ('   Total 1-fluid dust mass:',1PE14.5)
         ENDIF
c        Write energies for object 1; this is the same as total energies
c        above if only one object present.
c        Using short-hand notation to not confuse grepping with total
c        values
         IF (n2.NE.0) THEN
            IF (ABS(tgravtot_o1).GT.tiny) THEN
               alph = ttermtot_o1 / ABS(tgravtot_o1)
               betl = trotztot_o1 / ABS(tgravtot_o1)
               betr = trotxtot_o1 / ABS(tgravtot_o1)
               betm = tmagtot_o1  / ABS(tgravtot_o1)
            ELSE
               alph = 0.
               betl = 0.
               betr = 0.
            ENDIF
            WRITE (iprint, 99300, ERR=100)tkintot_o1,
     &                                    trotztot_o1,trotxtot_o1,
     &                                    tgravtot_o1,ttermtot_o1,
     &                                    tmagtot_o1,alph,betl,betr,betm
99300       FORMAT (' E_kinetic              : ', 1PE14.5, /,
     &              ' E_rotational_about_z   : ', 1PE14.5, /,
     &              ' E_rotational_about_x   : ', 1PE14.5, /,
     &              ' E_potential            : ', 1PE14.5, /,
     &              ' E_internal             : ', 1PE14.5, /,
     &              ' E_magnetic             : ', 1PE14.5, /,
     &              ' alpha                  : ', 1PE14.5, /,
     &              ' beta_para (z)          : ', 1PE14.5, /,
     &              ' beta_perp (x)          : ', 1PE14.5, /,
     &              ' beta_mag               : ', 1PE14.5, /)
         ENDIF
c
c--Object no 2
c
         IF (n2.NE.0) THEN
            WRITE (iprint, 99007, ERR=100) n2,cmx2tot,cmy2tot,
     &           cmz2tot,vcmx2tot,vcmy2tot,vcmz2tot,hmi2tot,hma2tot,
     &           dmax2tot,zmax2tot,romean2tot,
     &           romax2tot,rocen2tot,valphamin2tot,valphamax2tot,
     &           tgmean2, tgmax2, tgcen2
99007       FORMAT (/, ' Object number 2 (', I10, ' particles ) : ', /,
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
               IF (idustRT.GT.0) THEN
                  WRITE (iprint, 99107, ERR=100) tdmean1tot,tdmax1tot,
     &                 tdmin1tot
                  WRITE (iprint, 99108, ERR=100) h2mass1tot
               ENDIF
            ENDIF
         ENDIF

#ifdef NONIDEAL
c
c--non-ideal MHD
c
         CALL nimhd_write_ev(gt,nptot,iprint,et,nimhdfile
     &                       ,evdatamintot,evdatamaxtot,evdatasumtot)

         WRITE(iprint, 98010) itry_array
98010 FORMAT(/, ' Non-ideal MHD iterations of ion densities:', 263I10)
         itry_array = 0
#endif

#ifdef MPI
c Note: the following two blocks of info (ptmasses and RT) are produced
c       for all MPI processes because no MPI reduction has been
c       implemented. Output goes to their respective ascii files not to
c       the main iprint (which contains iproc=1's data)
      ENDIF
#endif

#ifdef RT
c
c--RT stats
c
         IF (iradtrans.EQ.idim.AND.encal.EQ.'r') THEN
            WRITE(iprint, 98020)nintegsuc,nintegfail,
     &           REAL(nintegtot)/(REAL(nintegsuc)+tiny),nrtnca,
     &           nrtnc1,nrtnc2,nrtmaxit
            nintegfail = 0
            nintegsuc  = 0 
            nrtnc1     = 0
            nrtnc2     = 0
            nrtnca     = 0
            nrtmaxit   = 0
            nintegtot  = 0
         ENDIF
98020    FORMAT(/, ' Implicit R-T (GSIMPL) :', /,
     &    ' RT-Succeeded -',I8,' Failed    -',I8,' Av. iter. - ',F8.4,/,
     &    ' RT-Failures:',11X,'Almost N-C -', I8,' N-C 1     - ',I8, /, 
     &    22X,'N-C 1-2    - ', I8, ' Max iter. - ', I8)
#endif
c
c--Implicit one-fluid dust stats
c
         IF (idustFluid.EQ.-2) THEN
            WRITE(iprint, 98030)ndustimpl_suc,ndustimpl_fail,
     &           REAL(ndustimpl_tot)/(REAL(ndustimpl_suc)+tiny)
            ndustimpl_fail = 0
            ndustimpl_suc  = 0 
            ndustimpl_tot  = 0
            WRITE(iprint, 98040)ndustdiff_suc,ndustdiff_fail,
     &           REAL(ndustdiff_tot)/(REAL(ndustdiff_suc)+tiny)
            ndustdiff_fail = 0
            ndustdiff_suc  = 0
            ndustdiff_tot  = 0
         ENDIF
98030    FORMAT(/, ' Implicit one-fluid dust dynamics        :', /,
     &    ' Dy-Succeeded -',I8,' Failed    -',I8,' Av. iter. - ',F8.4)
98040    FORMAT(/, ' Implicit one-fluid dust alpha diffusion :', /,
     &    ' Df-Succeeded -',I8,' Failed    -',I8,' Av. iter. - ',F8.4)
c
c--ptmass stats
c
         ido = iaccf1+iaccf2+iaccf3+iaccf4+iaccf5+iaccf6
     &       +iaccph1+iaccph2+iaccsuc+iaccmg 
         IF (ido .GT. 0) THEN
           WRITE(iprint,98002)iaccsuc,iaccmg,iaccph1,iaccph2,
     &           iaccf6,iaccf5,iaccf1,iaccf2,iaccf3,iaccf4
         ENDIF
         iaccsuc = 0
         iaccmg  = 0
         iaccph1 = 0
         iaccph2 = 0
         iaccf1  = 0
         iaccf2  = 0
         iaccf3  = 0
         iaccf4  = 0
         iaccf5  = 0
         iaccf6  = 0
98002 FORMAT(/, ' ptmass - created    : ', I8, ' merged  : ', I8, /,
     &           ' passed - h          : ', I8, ' and dt  : ', I8, /,
     &           ' failed - h          : ', I8, ' dt      : ', I8, /,
     &           '          divai      : ', I8, ' alpha   : ', I8, /,
     &           '          alpha+beta : ', I8, ' total e : ', I8)

c
c-- Begin next dtmax
c
        WRITE(iprint,98003)
98003 FORMAT(/, ' ------------------------------------------------',//)
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
99010    FORMAT (' the inner ', F12.3,
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

         IF (what.EQ.'scratch' .OR. what.EQ.'s') THEN
            WRITE (iprint, 99012)
99012       FORMAT (' New initial conditions defined from scratch')
            WRITE (iprint, 99013) idist
99013       FORMAT (
     &           ' Particles distributed according to distribution no :'
     &           , I2)

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
99017       FORMAT (/, ' first   object is made of ', I8,
     &              ' particles with offsets : ', /,
     &              '  x : ', 1PE12.5, '   y : ',
     &              1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,
     &              '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
            WRITE (iprint, 99018) n2, xx2, yy2, zz2, vvx2, vvy2, vvz2
99018       FORMAT (/, ' second  object is made of ', I8,
     &              ' particles with offsets : ', /,
     &              '  x : ', 1PE12.5, '   y : ',
     &              1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,
     &              '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
            WRITE (iprint, 99019) angto
99019       FORMAT (/, ' total angular momentum of the system : ',
     &              1PE12.5)
         ENDIF
      ENDIF
c
      GOTO 200
c
c--Handle errors during writing
c
 100  where = 'prout'
      CALL error(where, 1)

 200  CALL FLUSH (iprint)

#ifdef MPIDEBUG
      print *,iproc,': Exited prout'
#endif

      RETURN
      END
