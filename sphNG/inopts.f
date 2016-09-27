      SUBROUTINE inopts
c************************************************************
c                                                           *
c  This subroutine defines all options desired for the run  *
c                                                           *
c************************************************************

      INCLUDE 'idim'

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
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/init'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/varmhd'
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine inopts')
c
c--Open input file
c
      OPEN (iterm, FILE='inspho', STATUS='old', ACTION='read')
c
c--Determine options for insph
c
c--Read options
c
      READ (iterm, *) igrp
      IF (imhd.EQ.idim) THEN
         READ (iterm, *) ifsvi,alpha,beta,alphamin(2),alphamax(2)
      ELSE
         READ (iterm, *) ifsvi,alpha,beta
      ENDIF
c--NB: we always read MHD options from inspho, even
c  if doing a hydro run, so that all inspho's are compatible
c  (and should contain this line)      
      READ (iterm, *) ndivBsubcycles,overclean_fac,iresist,etamhd
      READ (iterm, *) ifcor
      READ (iterm, *) ichoc
      READ (iterm, *) iener
      READ (iterm, *) damp
      READ (iterm, *) iexf
      READ (iterm, *) iexpan
      READ (iterm, *) nstep
      READ (iterm, *) nfullstep
      READ (iterm, *) iptoutnum
      READ (iterm, *) tol, tolptm, tolh
      READ (iterm, *) irec
      READ (iterm, *) tmax
      READ (iterm, *) tstop
      READ (iterm, *) dtmax
      READ (iterm, *) dtini
      omeg0 = 0.
      IF (ifcor.NE.0) THEN
         READ (iterm, *) omeg0
      ENDIF
      vexpan = 0.
      IF (iexpan.NE.0) THEN
         READ (iterm, *) vexpan
      ENDIF

      IF (idebug(1:7).EQ.'inopts') THEN
         WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc, iener,
     &        ibound, damp, varsta
99004    FORMAT (1X, 7(I2,1X), 2(E12.5,1X), 1X, A7)
         WRITE (iprint, 99005) file1, irec, nstep
99005    FORMAT (1X, A7, 1X, I4, 1X, I4)
      ENDIF

      CLOSE (iterm)

      RETURN
      END
