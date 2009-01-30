      PROGRAM convert_full_small
c
c--Converts full dump to small dump
c
      INCLUDE '../idim'

      PARAMETER (maxnumber = 10000)

      INCLUDE '../COMMONS/tming'
      INCLUDE '../COMMONS/units'
      INCLUDE '../COMMONS/sort'
      INCLUDE '../COMMONS/stepopt'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      CHARACTER*1 ianswer
      CHARACTER*7 fullname(maxnumber)
      CHARACTER*9 smallname

      WRITE (*,*) 'Do you wish to overwrite full with small dumps?'
      WRITE (*,*) '    (If you answer anything but y, the small dumps'
      WRITE (*,*) '    will have the original name appended with .s)'
      READ (*,88001) ianswer
88001 FORMAT(A1)

      WRITE (*,*) 'How many full dumps?'
      READ (*,*) numfulldumps
      IF (numfulldumps.GT.maxnumber) THEN
         WRITE (*,*) 'Filename array too small'
         STOP
      ENDIF

      WRITE (*,*) 'Enter names of dump files'
      DO i = 1, numfulldumps
         READ (*,*) fullname(i)
      END DO

      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO

      DO iloop = 1, numfulldumps

c      WRITE (*,*) 'Enter name of (full) file to convert'
c      READ (*,*) fullname
         IF (ianswer.EQ.'y') THEN
            smallname = fullname(iloop)
         ELSE
            smallname = fullname(iloop) // '.s'
         ENDIF

         OPEN (11,FILE=fullname,FORM='unformatted')

         CALL rdump(11,ichkl,1)

         CLOSE (11)
         IF (ichkl.EQ.1) THEN
            WRITE (*,*) 'ERROR reading input file: ichkl ',ichkl
            STOP
         ENDIF

         umass = umassi
         udist = udisti
         utime = utimei
         umagfd = umagfdi
         dtmax = dtmaxdp

         OPEN (11,FILE=smallname,FORM='unformatted')

         nfullstep = 100
         CALL wdump(11)
      
         CLOSE (11)

      END DO

      END


      SUBROUTINE quit
      STOP
      END


      SUBROUTINE endrun
      CALL quit
      END
