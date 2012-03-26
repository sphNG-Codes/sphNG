      SUBROUTINE options
c************************************************************
c                                                           *
c  This subroutine defines all options desired for the run  *
c                                                           *
c************************************************************

      INCLUDE 'idim'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/dissi'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/integ'
      INCLUDE 'COMMONS/varet'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/crpart'
      INCLUDE 'COMMONS/ptbin'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/prdrag'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine options')
c
c--Open input file
c
      OPEN (iterm, FILE=inname)
c
c--Determine options for evolution run
c
c--Read name of run
c
      READ (iterm, 99002) namerun
      WRITE (namenextrun,99002) namerun
      namelength = LEN(namenextrun)
99002 FORMAT (A20)
      DO i = namelength, 1, -1
         IF (namenextrun(i:i).GE.'0' .AND. namenextrun(i:i).LE.'9') THEN
            READ (namenextrun(i:i),88001) number
88001       FORMAT(I1)
            icarry = 0
            number = number + 1
            IF (number.GT.9) THEN
               icarry = 1
               number = 0
            ENDIF
            WRITE (namenextrun(i:i),88001) number
            j = i-1
            IF ((j.GT.0) .AND. (icarry.EQ.1) .AND.
     &           (namenextrun(j:j).GE.'0' .AND. 
     &           namenextrun(j:j).LE.'9')) THEN
               READ (namenextrun(j:j),88001) number
               number = number + icarry
               WRITE (namenextrun(j:j),88001) number
            ENDIF
            GOTO 100
         ELSEIF (namenextrun(i:i).NE.' ') THEN
            GOTO 100
         ENDIF
      ENDDO
c
c--Open output file
c 
 100  CONTINUE
#ifdef MPIALL
      IF (iproc.EQ.0) THEN
#endif
      IF(iprint.NE.6) OPEN (iprint, FILE=namerun)
#ifdef MPIALL
      ENDIF
#endif

      CALL labrun
c
c--Read name of file containing physical input
c
      READ (iterm, 99003) file1
99003 FORMAT (A7)
      READ (iterm, 99003) varsta
      IF (imhd.EQ.idim) THEN
         READ (iterm, *) varmhd
      ELSE
         READ (iterm, *)
      ENDIF
c
c--Read options
c
      READ (iterm, 99006) encal
99006 FORMAT (A1)
      IF (encal.EQ.'r') THEN
         IF (iradtrans.NE.idim) THEN
            WRITE (iprint,*) 'ERROR - iradtrans.NE.idim'
            CALL quit
         ENDIF
         READ (iterm, *) tolerance, boundtemp, bounddens, opdenom, 
     &        iopmodel, metallicity
      ENDIF

      print *,'Opacity model ',iopmodel

#ifdef MPI
      READ (iterm, *) mpitype, xlocationmpi, ylocationmpi, 
     &     zlocationmpi,densmpi
      print *,xlocationmpi, ylocationmpi,
     &     zlocationmpi,densmpi
#else
      READ (iterm, *)
#endif

      READ (iterm, *) initialptm
      READ (iterm, 99006) iaccevol
      IF (iaccevol.EQ.'v' .OR. iaccevol.EQ.'s') READ (iterm, *) accfac
      READ (iterm, *) iptmass
      READ (iterm, *) igrp
      READ (iterm, *) igphi
      IF (imhd.EQ.idim) THEN
         READ (iterm, *) ifsvi, alpha, beta, alphamin(2), alphamax(2)
         READ (iterm, *) iresist,etamhd
      ELSE
         READ (iterm, *) ifsvi, alpha, beta
         READ (iterm, *)
      ENDIF
c      IF (ifsvi.EQ.6) THEN
         alphamin(1) = alpha
         alphamax(1) = beta
c      ENDIF
      READ (iterm, *) ifcor
      READ (iterm, *) ichoc
      READ (iterm, *) iener
      READ (iterm, *) damp
      READ (iterm, *) ibound
      READ (iterm, *) iexf
      READ (iterm, *) iexpan
      READ (iterm, *) nstep
      IF (nstep.LT.1) nstep = 1
      READ (iterm, *) nfullstep
      IF (nstep.LT.1) nfullstep = 1

      READ (iterm, *) iptoutnum
      xlog2 = 0.30103 + 0.00001
      ibin = INT(LOG10(FLOAT(iptoutnum))/xlog2) + 1
      iptoutnum = 2**ibin

      READ (iterm, *) tol, tolptm, tolh
      READ (iterm, *) ipos
      READ (iterm, *) tmax
      READ (iterm, *) tstop
      READ (iterm, *) dtmax
      READ (iterm, *) dtini

      omeg0 = 0.
      IF (ifcor.NE.0) THEN
         READ (iterm, *) omeg0
         omeg0 = omeg0*utime
      ENDIF

      vexpan = 0.
      IF (iexpan.NE.0) THEN
         READ (iterm, *) vexpan
         vexpan = vexpan*utime/udist
      ENDIF

      pext = 0.
      hmaximum = 0.
      IF (ibound.EQ.7) THEN
         READ (iterm, *) hmaximum
         READ (iterm, *) pext
      ENDIF

      IF (ibound.EQ.8) THEN
         READ (iterm, *) deadbound
         READ (iterm, *) fractan, fracradial, nstop, nfastd
      ENDIF
      IF (ibound/10.EQ.9) THEN
         READ (iterm, *) deadbound
         READ (iterm, *) fractan, fracradial, nshell, rshell
      ENDIF

      imigrate = 0
      r_planetesimal = 0.0
      rho_planetesimal = 0.0
      use_tprof = .FALSE.
      prswitch = .FALSE.
      gasdrag = .FALSE.
      IF (ibound.EQ.100) THEN
         READ (iterm, *) gasdrag
         IF (gasdrag) THEN
            READ (iterm, *) r_planetesimal
            r_planetesimal = r_planetesimal*1.E5/udist
            READ (iterm, *) rho_planetesimal
            rho_planetesimal = rho_planetesimal/udens
            print *, 'r_pla ', r_planetesimal, rho_planetesimal
         ENDIF
         READ (iterm, *) variation,phibound,flowrate,signorm,
     &        hoverr,partm
      ENDIF

      IF (ibound.EQ.102) THEN
         READ (iterm, *) use_tprof
         IF (use_tprof) READ (iterm, *) tprof
         READ (iterm, *) prswitch
         IF (prswitch) READ (iterm, *) pbeta
         READ (iterm, *) gasdrag
         IF (gasdrag) THEN
            READ (iterm, *) r_planetesimal
            r_planetesimal = r_planetesimal*1.E5/udist
            READ (iterm, *) rho_planetesimal
            rho_planetesimal = rho_planetesimal/udens
            print *, 'r_pla ', r_planetesimal, rho_planetesimal
         ENDIF
         READ (iterm, *) signorm,hoverr,partm
         READ (iterm, *) imigrate
         IF (imigrate.EQ.1) THEN
            READ (iterm, *) pmrate
            pmrate = pmrate*(1.49598E13/3.1536E13)*(utime/udist)
            print *,'pmrate in options = ', pmrate
            READ (iterm, *) rorbitmax
         ENDIF
      ENDIF

      xmass = 0.
      iandrea = 0
      IF (iexf.EQ.5 .OR. iexf.EQ.6 .OR. iexf.EQ.7) THEN
         READ (iterm, *) xmass
      ENDIF
      IF (iexf.EQ.7 .OR. (ibound.EQ.102 .AND. initialptm.EQ.0)) THEN
         READ (iterm, *) planetmass
         READ (iterm, *) rplanet
         IF (iexf.EQ.7) THEN
            READ (iterm, *) iandrea
            READ (iterm, *) eccentricity_p
            READ (iterm, *) inclination_p
         ENDIF
      ENDIF

      IF (ibound.EQ.102 .AND. initialptm.EQ.5) THEN
            READ (iterm, *) uniformtslim
            uniformtslim2 = uniformtslim*uniformtslim
      ENDIF

      IF (iptmass.NE.0 .OR. initialptm.NE.0) THEN
         READ (iterm, *) hacc
         READ (iterm, *) haccall
      ENDIF

      IF (iptmass.NE.0) THEN
         READ (iterm, *) radcrit
         READ (iterm, *) ptmcrit
      ENDIF

      IF (encal.EQ.'c') THEN
         READ (iterm, *) coolingrate
      ENDIF
c
c--Read boundaries
c
      rmax = 0.
      rcyl = 0.
      rmind = 0.
      IF (ibound.EQ.1 .OR. ibound.EQ.3 .OR. ibound.EQ.8 .OR. 
     &                                           ibound/10.EQ.9) THEN 
          READ (iterm, *) rmind, rmax, 
     &                   xmin, xmax, ymin, ymax, zmin, zmax
      ELSEIF (ibound.EQ.2) THEN
          READ (iterm, *) rmind, rcyl, 
     &                   xmin, xmax, ymin, ymax, zmin, zmax
      ELSE
          READ (iterm, *) rd1, rd2,
     &                   xmin, xmax, ymin, ymax, zmin, zmax
          IF (rd1.GE.0.) THEN
             rcyl = rd2
             rmind = rd1
          ELSE
             rmax = rd2
          ENDIF
      ENDIF
c
c--Check for consistency
c
      CALL chekopt

      IF (idebug(1:7).EQ.'options') THEN
         WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc, iener,
     &                         ibound, damp, varsta
99004    FORMAT (1X, 7(I2,1X), E12.5, 1X, A7)
         WRITE (iprint, 99005) file1, ipos, nstep

         IF (ibound.EQ.1 .OR. ibound.EQ.3) 
     &   WRITE (iprint, *) rmax, xmin, xmax, ymin, ymax, zmin, zmax

         IF (ibound.EQ.2) WRITE (iprint, *) rmind, rcyl, 
     &                   xmin, xmax, ymin, ymax, zmin, zmax
99005    FORMAT (1X, A7, 1X, I4, 1X, I4)
      ENDIF

      CLOSE(iterm)

      RETURN
      END
