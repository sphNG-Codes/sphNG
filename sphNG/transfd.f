      SUBROUTINE transfd
c************************************************************
c                                                           *
c  This subroutine allows the transfer  of one given dump   *
c     of portion of file to another file with possibility   *
c     of making changes                                     *
c                                                           *
c************************************************************

      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/trans'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/diskbd'
      INCLUDE 'COMMONS/expan'
      INCLUDE 'COMMONS/rotat'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/tming'

      CHARACTER*7 where
      CHARACTER*1 iok, smooth

      DATA where/'transfd'/
c
c--Allow for tracing flow
c
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine transfd')

      WRITE (*, 99002)
99002 FORMAT (' BEGIN TRANSFER', //)
      WRITE (*, 99003)
99003 FORMAT (' give name of this run :')
      READ (*, 99004) namerun
99004 FORMAT (A20)
      WRITE (*, 99005)
99005 FORMAT (' give name of file containing data')
      READ (*, 99006) file1
99006 FORMAT (A7)
      WRITE (*, 99007)
99007 FORMAT (' give name of file where to write transfer')
      READ (*, 99006) file2

      OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)
      OPEN (idisk2, FILE=file2, STATUS='unknown', FORM='unformatted',
     &        RECL=imaxrec)

      WRITE (*, 99008)
99008 FORMAT (' first and last reccord to be transfered')
      READ (*, *) ibegin, iend
      WRITE (*, 99009)
99009 FORMAT (' modification of data desired? (y/n)')
      READ (*, 99010) iok
99010 FORMAT (A1)
      ichang = 0

      IF (iok.EQ.'y') THEN
         WRITE (*, 99011)
99011    FORMAT (' what type of modifications do yo want?', /,
     &           ' change energy in an inner sphere            : 1', /,
     &           ' same as 2 and change reference frame        : 2')
         READ (*, *) ichang
         IF (ichang.GE.1) THEN
            WRITE (*, 99012)
99012       FORMAT (' fraction of theparticles to change')
            READ (*, *) frac
            WRITE (*, 99013)
99013       FORMAT (' new value of internal energy: ')
            READ (*, *) energc
            IF (ichang.GE.2) THEN
               WRITE (*, 99014)
99014          FORMAT (' omega (1/s), expansion velocity (cm/s)',
     &              ' and distance (cm)')
               READ (*, *) omeg0, vexpan, rnorm
               omeg0 = omeg0*utime
               vexpan = utime*vexpan/udist
               rnorm = rnorm/udist
            ENDIF
         ENDIF
      ENDIF
      WRITE (*, 99015)
99015 FORMAT (' smoothing data? (y/n)')
      READ (*, 99010) smooth
c
c--If boundaries take them into account for smoothing
c
      rmax = 1.E30
      xmin = -1.E30
      xmax = 1.E30
      ymin = -1.E30
      ymax = 1.E30
      zmin = -1.E30
      zmax = 1.E30

      IF (smooth.EQ.'y') THEN
         WRITE (*, 99016)
99016    FORMAT (' give boundary type : 0 - no boundaries', /,
     &           '                      1 - cartesian    ', /,
     &           '                      2 - cylindrical  ', /,
     &           '                      3 - spherical    ', /)
         READ (*, *) ibound
         IF (ibound.EQ.1) THEN
            WRITE (*, 99017)
99017       FORMAT (' give xmin,xmax,ymin,ymax,zmin,zmax')
            READ (*, *) xmin, xmax, ymin, ymax, zmin, zmax
         ENDIF
         IF (ibound.EQ.2) THEN
            WRITE (*, 99018)
99018       FORMAT (' give rcyl, zmin, zmax')
            READ (*, *) rcyl, zmin, zmax
            xmax = rcyl
            ymax = rcyl
            xmin = - xmax
            ymin = - ymax
            WRITE (*, 98018)
98018       FORMAT (' do you want inner boundary? (y/n)')
            READ (*, 99010) iok
            IF (iok.EQ.'y') THEN
               WRITE (*, 98019)
98019          FORMAT(' enter inner boundary ')
               READ (*, *) rmind
            ENDIF
         ENDIF
         IF (ibound.EQ.3) THEN
            WRITE (*, 99019)
99019       FORMAT ('give rmax ')
            READ (*,*) rmax
            zmax = rmax
            xmax = rmax
            ymax = rmax
            zmin = - zmax
            xmin = - xmax
            ymin = - ymax
         ENDIF
      ENDIF

      WRITE (*, 99020)
99020 FORMAT (' ALL OPTIONS SET, PROCESSING DATA')
c
c--Write label
c
      CALL labrun
c
c--Go to begining of first reccord
c
      CALL place(idisk1, ibegin, irec, 1)
c
c--Start transfering
c
      DO i = ibegin, iend
         CALL rdump(idisk1, ichkl, 1)
         IF (ichkl.EQ.1) CALL error(where, ichkl)
c
c--Allow for modifications
c
         IF (ichang.NE.0) CALL modif
c
c--Smooth data
c
         IF (smooth.EQ.'y') CALL smoothd

         nfullstep = 1
         CALL wdump(idisk2)

      END DO
c
c--Write output
c
      CALL prout(where)

      IF (idebug(1:5).EQ.'trans') THEN
         WRITE (iprint, 99021) file1, file2
99021    FORMAT (1X, A7, 1X, A7)
         WRITE (iprint, 99022) ibegin, iend, ichang, frac
99022    FORMAT (1X, 2(I4,1X), I2, 1X, F5.3, 1X, I2)
      ENDIF

      CALL endrun

      RETURN
      END
