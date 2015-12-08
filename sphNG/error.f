      SUBROUTINE error(where, iwhat)
c************************************************************
c                                                           *
c  This routine identifies and handle the errors detected   *
c     by the program during execution                       *
c                                                           *
c************************************************************

      INCLUDE 'idim'

      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/active'

      CHARACTER*7 where
c
c--Write error message
c
      WRITE (iprint, 99001) where
99001 FORMAT (' ERROR DETECTED IN ROUTINE : ', A7)
c
c--Read error
c
      IF (where(1:5).EQ.'place') THEN
         WRITE (iprint, 99002)
99002    FORMAT (' EOF on disk, check record number !')
         CALL quit(0)
      ENDIF
c
c--Write error
c
      IF (where(1:5).EQ.'wdump') THEN
         WRITE (iprint, 99003) iwhat
99003    FORMAT (' error number : ', I4)
         CALL quit(0)
      ENDIF
c
c--Integration error
c
      IF (where(1:5).EQ.'integ') THEN
         IF (iwhat.EQ.1) THEN
            WRITE (iprint, 99004)
99004       FORMAT (' time step is smaller than 1e-10, check data!')
            CALL quit(1)
         ENDIF
         IF (iwhat.EQ.2) THEN
            WRITE (iprint, 99005)
99005       FORMAT (' number of rejections greater than 10, check data')
            CALL quit(1)
         ENDIF
         IF (iwhat.EQ.3) THEN
            WRITE (iprint, 99006)
99006       FORMAT (' integration is not converging, nothing can help!')
            CALL quit(1)
         ENDIF
         IF (iwhat.EQ.4) THEN
            WRITE (iprint, 99007)
99007       FORMAT (' tolerance was increased temporarily, too many ',
     &              'rejections were counted')
            RETURN
         ENDIF
         IF (iwhat.EQ.5) THEN
            WRITE (iprint, 99008)
99008       FORMAT (' tolerance was increased temporarily, integration',
     &              ' was no longer converging')
            RETURN
         ENDIF
      ENDIF
c
c--Error in transfer
c
      IF (where(1:5).EQ.'trans') THEN
         WRITE (iprint, 99009)
99009    FORMAT (' not so many dumps on input file !')
         CALL quit(0)
      ENDIF
c
c--Error in adding dumps
c
      IF (where(1:6).EQ.'addump') THEN
         IF (iwhat.EQ.1) WRITE (iprint, 99002)
         IF (iwhat.EQ.2) THEN
            WRITE (iprint, 99010)
99010       FORMAT (' smoothing length of both object not equal !')
            CALL quit(0)
         ENDIF
      ENDIF
c
c--Error in building tree
c
      IF (where(1:5).EQ.'mtree') THEN
         WRITE (iprint, 99011) iwhat
99011    FORMAT (' maximum number of nodes overflow !',I2)
         CALL quit(1)
      ENDIF
      IF (where(1:6).EQ.'naybor') THEN
         IF (iwhat.EQ.1) WRITE (iprint, 99012) iwhat
99012    FORMAT (' hash table too short !')
         IF (iwhat.EQ.2) WRITE (iprint, 99013) iwhat
99013    FORMAT (' map tables too short !')
         CALL quit(0)
      ENDIF
      IF (where(1:6).EQ.'indexx') THEN
         WRITE (iprint, 99014)
99014    FORMAT (' nstack must be made larger !')
         CALL quit(0)
      ENDIF
c
c--Error in prout
c
      IF (where(1:5).EQ.'prout') THEN
         WRITE (iprint, 99015)
99015    FORMAT (' error in writing dump on output file !')
         CALL endrun
      ENDIF
c
c--Error in options
c
      IF (where(1:7).EQ.'chekopt') THEN
         WRITE (iprint, 99016)
99016    FORMAT (' incompatible options defined in input !')
         IF (iwhat.EQ.2) WRITE (iprint, 99101)
99101    FORMAT (' general damping incompatible with total',
     &           ' energy conservation')
         IF (iwhat.EQ.3) WRITE (iprint, 99102)
99102    FORMAT (' isothermal eos cannot have shock or work',
     &           ' contributions')
         IF (iwhat.EQ.4) WRITE (iprint, 99103)
99103    FORMAT (' adiabatic eos must have shock and work',
     &           ' contributions')
         IF (iwhat.EQ.5) WRITE (iprint, 99104)
99104    FORMAT (' polytropic eos cannot have shock or work',
     &           ' contributions')
         IF (iwhat.EQ.6) WRITE (iprint, 99105)
99105    FORMAT (' iptintree = 2 cannot be used with point masses',
     &           ' without self-gravity')
         IF (iwhat.EQ.7) WRITE (iprint, 99106)
99106    FORMAT (' using real resistivity from function (iresist=2)',/,
     &           ' but multiplying factor etamhd is <= 0,',
     &           ' usually should be 1.0',/,
     &           ' (to use etamhd = 0.0, select iresist=1) ')

         CALL quit(0)
      ENDIF
c
      IF (where(1:7).EQ.'density') THEN
         WRITE (iprint, 99018) iwhat
99018    FORMAT (' dimension for neighbor list too short, iover=', I8)
         CALL endrun
      ENDIF
c
c--Error in ghost particle number
c
      IF (where(1:4).EQ.'ghos' .OR. 
     &     (where(1:5).EQ.'rdump' .AND. iwhat.EQ.3)) THEN
         WRITE (iprint, 99019) iwhat
99019    FORMAT (' number of particles+ghosts exceeds dimensions!', /,
     &           ' dimensions needed at least :', I8)
         CALL endrun
      ENDIF
c
c--Neighbours list dimensions too small
c
      IF ((where(1:5).EQ.'treef').OR.
     &     (where(1:5).EQ.'insul' .AND. iwhat.EQ.1)) THEN
         WRITE (iprint, 99020)
99020    FORMAT (' number of neighbours exceeds dimensions')
         CALL endrun
      ENDIF
c
c--Point mass neighbours list dimensions too small
c
      IF ((where(1:7).EQ.'accrete' .AND. iwhat.EQ.2) .OR. 
     &     (where(1:7).EQ.'gforspt')) THEN
         WRITE (iprint, 99021)
99021    FORMAT (' number of point mass neighbours exceeds dimensions')
         CALL endrun
      ENDIF
c
c--Units the code was compiled with do not match binary dump file units
c
      IF (where(1:5).EQ.'rdump' .AND. iwhat.NE.3) THEN
         WRITE (iprint, 99022)
99022    FORMAT (' incompatible units between code and binary file')
         IF (iwhat.EQ.1) WRITE (iprint,*)'   Units: udist'
         IF (iwhat.EQ.2) WRITE (iprint,*)'   Units: umass'
         IF (iwhat.EQ.4) WRITE (iprint,*)'   Units: umagfd'
         CALL quit(0)
      ENDIF
c
c--Number of gas particles less than minimum value
c
      IF ((where(1:4).EQ.'step') .AND. iwhat.EQ.1) THEN
         WRITE (iprint, 99023) nactive - nptmass
99023    FORMAT (' insufficient gas particles remaining ',I4)
         CALL endrun
      ENDIF
c
c--Point mass accretion error
c
      IF (where(1:7).EQ.'accrete' .AND. iwhat.EQ.1) THEN
         WRITE (iprint, 99024)
99024    FORMAT (' accretion error')
         CALL endrun
      ENDIF
c
c--Equation of state not defined
c
      IF (where(1:5).EQ.'eospg') THEN
         WRITE (iprint, 99025)
99025    FORMAT (' equation of state not defined')
         CALL endrun
      ENDIF
c
c--Accretion radius not defined
c
      IF (where(1:6).EQ.'header') THEN
         IF (iwhat.EQ.1) THEN
            WRITE (iprint, 99026)
         ELSE
            WRITE (iprint, 99027)
         ENDIF
99026    FORMAT ('Outer Accretion radius not defined')
99027    FORMAT ('Inner Accretion radius not defined')
         CALL endrun
      ENDIF
c
c--Point masses present without telling program they exist
c
      IF (where(1:6).EQ.'preset' .AND. iwhat.EQ.2) THEN
         WRITE (iprint, 99028)
99028    FORMAT (' Point Masses without initialptm set')
         CALL endrun
      ENDIF
c
c--Problem with treef
c
      IF (where(1:5).EQ.'treef') THEN
         IF (iwhat.EQ.1) THEN
            WRITE (iprint, 99029)
99029       FORMAT (' Called treef non-existant particle!!')
         ELSEIF (iwhat.EQ.2) THEN
            WRITE (iprint, 99030)
99030       FORMAT (' Called treef for point mass!!')
         ENDIF
         CALL endrun
      ENDIF
c
c--Infinite loop to set h for reassigned particle in 'step'
c
      IF ((where(1:4).EQ.'step') .AND. iwhat.EQ.2) THEN
         WRITE (iprint, 99031) 
99031    FORMAT (' Infinite loop to set h for reassigned particle')
         CALL endrun
      ENDIF
c
c--Too many new particles created
c
      IF ((where(1:4).EQ.'step') .AND. iwhat.EQ.3) THEN
         WRITE (iprint, 99032) 
99032    FORMAT (' Total number of particles > idim')
         CALL endrun
      ENDIF
c
c--Attempt to increase synchronisation time beyond possible timestep range
c
      IF ((where(1:5).EQ.'mesop') .OR. 
     &     (where(1:5).EQ.'rdump'.AND.iwhat.EQ.4)) THEN
         WRITE (iprint, 99033) iwhat
99033    FORMAT (' Attempt to increase synchron time beyond range ',I2)
         IF (iwhat.GE.2) CALL endrun
      ENDIF
c
c--GRAPE flag igrape not defined
c
      IF (where(1:6).EQ.'derivi' .AND. iwhat.EQ.1) THEN
         WRITE (iprint, 99040)
99040    FORMAT (' GRAPE flag igrape not defined')
         CALL quit(0)
      ENDIF
      IF (where(1:5).EQ.'insul' .AND. iwhat.EQ.2) THEN
         WRITE (iprint, 99045)
99045    FORMAT (' Code compiled with insulate_TREE rather',
     &        ' than insulate_GRAPE')
         CALL quit(0)
      ENDIF
      IF (where(1:5).EQ.'insul' .AND. iwhat.EQ.3) THEN
         WRITE (iprint, 99050)
99050    FORMAT (' Code compiled with insulate_GRAPE rather',
     &        ' than insulate_TREE')
         CALL quit(0)
      ENDIF
c
c--Number of point masses exceeds dimensions
c
      IF (where(1:7).EQ.'accrete' .AND. iwhat.EQ.3) THEN
         WRITE (iprint, 99055)
99055    FORMAT (' Number of point masses exceeds dimensions')
         CALL endrun
      ENDIF

c
c--quit by default on unknown errors
c
      WRITE(iprint, 99000)
99000 FORMAT(' ERROR: Unknown error.')
      CALL quit(0)

      RETURN
      END
