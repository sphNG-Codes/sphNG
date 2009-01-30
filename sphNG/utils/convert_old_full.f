      PROGRAM convert_old_full
c
c--Converts old (pre-Bate-Price dump format) to new full dump format
c
      INCLUDE '../idim'

      INCLUDE '../COMMONS/tming'
      INCLUDE '../COMMONS/units'
      INCLUDE '../COMMONS/sort'
      INCLUDE '../COMMONS/stepopt'
      INCLUDE '../COMMONS/part'
      INCLUDE '../COMMONS/densi'
      INCLUDE '../COMMONS/ener2'
      INCLUDE '../COMMONS/ener3'
      INCLUDE '../COMMONS/timei'
      INCLUDE '../COMMONS/phase'
      INCLUDE '../COMMONS/ptmass'
      INCLUDE '../COMMONS/fracg'
      INCLUDE '../COMMONS/bodys'
      INCLUDE '../COMMONS/gtime'
      INCLUDE '../COMMONS/cgas'
      INCLUDE '../COMMONS/polyk2'
      INCLUDE '../COMMONS/binary'

      REAL*8 umassi, udisti, utimei, umagfdi
      COMMON /unitsin/ umassi, udisti, utimei, umagfdi
      COMMON /dtmaxin/ dtmaxdp

      DIMENSION rhoin(idim), dgravin(idim), iphasein(idim)

      CHARACTER*7 fullname
      CHARACTER*9 smallname

      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO

      WRITE (*,*) 'Enter name of (full) file to convert'
      READ (*,*) fullname
      smallname = fullname // '.s'

      OPEN (11,FILE=fullname,FORM='unformatted')

      READ (11, END=100) udisti, umassi, utimei,
     &     npart, n1, n2, gt, gamma, rhozero, RK2
     &     ,(xyzmh(5,i), i=1, npart), escap, tkin, tgrav, tterm,
     &     (xyzmh(1,i), i=1, npart), (xyzmh(2,i), i=1, npart),
     &     (xyzmh(3,i), i=1, npart), (vxyzu(1,i), i=1, npart),
     &     (vxyzu(2,i), i=1, npart), (vxyzu(3,i), i=1, npart),
     &     (vxyzu(4,i), i=1, npart), (xyzmh(4,i), i=1, npart)
     &     ,(rhoin(i), i=1, npart), (dgravin(i), i=1, npart),
     &     dtmaxdp, (isteps(i), i=1, npart),
     &     (iphasein(i), i=1, npart),
     &     nptmass, (listpm(i), i=1, nptmass)
c     &     ,(spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass),
c     &     (spinz(i),i=1,nptmass)
c     &     ,(angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass),
c     &     (angaddz(i),i=1,nptmass),
c     &     anglostx, anglosty, anglostz,
c     &     nreassign, naccrete, nkill, specang, ptmassin,
c     &     (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass),
c     &     (spinadz(i),i=1,nptmass)

      DO i = 1, npart
         iphase(i) = iphasein(i)
         rho(i) = rhoin(i)
      END DO

      print *,udisti, umassi, utimei,
     &     npart, n1, n2, gt, gamma, rhozero, RK2

      print *,dtmaxdp, nptmass

c      STOP

      CLOSE (11)
      IF (ichkl.EQ.1) THEN
         WRITE (*,*) 'ERROR reading input file: ichkl ',ichkl
         STOP
      ENDIF

      umass = umassi
      udist = udisti
      utime = utimei
      umagfd = 0.
      dtmax = dtmaxdp

      OPEN (11,FILE=smallname,FORM='unformatted')

      nfullstep = 1
      CALL wdump(11)
      
 100  CLOSE (11)

      END


      SUBROUTINE quit
      STOP
      END


      SUBROUTINE endrun
      CALL quit
      END
