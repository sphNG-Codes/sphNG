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
      INCLUDE '../COMMONS/phase'
      INCLUDE '../COMMONS/part'
      INCLUDE '../COMMONS/ptmass'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      CHARACTER*1 ianswer,ianswer2
      CHARACTER*7 fullname(maxnumber)
      CHARACTER*9 smallname

      WRITE (*,*) 'Do you wish to overwrite full with small dumps?'
      WRITE (*,*) '    (If you answer anything but y, the small dumps'
      WRITE (*,*) '    will have the original name appended with .s)'
      READ (*,88001) ianswer
88001 FORMAT(A1)

      WRITE (*,*) 'Do you wish to remove inactive particles ',
     &     '(iphase=-1)?'
      READ (*,88001) ianswer2

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

         IF (ianswer2.EQ.'y') THEN
            DO i = 1, nptmass
               listrealpm(listpm(i)) = i
            END DO

            icountparticles = 0
            DO i = 1, npart
               IF (iphase(i).NE.-1) THEN
                  icountparticles = icountparticles+1
                  IF (i.NE.icountparticles) THEN
                     CALL move_particle(i,icountparticles)
                  ENDIF
               ENDIF
            END DO
            WRITE (*,*) 'New number of particles: ',icountparticles,
     &           npart
            npart = icountparticles
            n1 = npart
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

      SUBROUTINE ghostp(ntot,npart,xyzmh,vxyzu,ekcle,Bevolxyz)

      STOP
      END

      SUBROUTINE insulate(ival, ntot, npart, dumxyzmh, f1vxyzu)

      STOP
      END


