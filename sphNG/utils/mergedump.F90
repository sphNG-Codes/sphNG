      PROGRAM mergedump
!************************************************
!  Written for cluster chem runs.
!  Merge together n dumps for analysis
!  Needs dumpmols.txt - same one as used for chem run
!  AY Feb 2019
!************************************************
#ifdef USEKROME
      INCLUDE 'COMMONS/krome_mods'
#endif

!#ifdef KROMEMPI
!      INCLUDE mpi
!#endif
!      INCLUDE 'omp_lib.h'                                                    
      INCLUDE 'idim'
#ifdef USEKROME
      INCLUDE 'COMMONS/kromevar'
#endif

      INCLUDE 'igrape'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/ener1'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/interstellar'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/xforce'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/savernd'
      INCLUDE 'COMMONS/abundances'
      INCLUDE 'COMMONS/perform'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/makeplt'
      INCLUDE 'COMMONS/radsink'

      INTEGER :: iindump,ichkl,i,j,n,istart,ptstart
      INTEGER :: numin,nptmasstmp,nparttmp
      INTEGER :: iuniquetmp(idim2),istepstmp(idim2)
      INTEGER :: iphasetmp(idim2),iout
      REAL :: xyzmhtmp(5,idim),vxyzutmp(4,idim)
      REAL :: rhotmp(idim2),gradhstmp(2,idim2)
      REAL :: alphaMMtmp(isizealphaMM,idim),ekcletmp(5,idim)
      REAL :: dust_tktmp(2,iradtrans2*idustRT+1)
      REAL :: chemistrytmp(nchemistry,iradtrans*idustRT+1)
      REAL :: heatingISRtmp(nheatingISR,iradtrans2*idustRT+1)
      REAL :: h2fractmp(nh2frac)
      REAL :: dh2dttmp(2,iradtrans*idustRT+1),Bevolxyztmp(imhdevol,imhd)
      REAL :: divcurlBtmp(5,imhd2),etarealtmp(imhd)
      REAL :: etaartificialtmp(imhd),spinxtmp(iptdim),spinytmp(iptdim)
      REAL :: spinztmp(iptdim),angaddxtmp(iptdim),angaddytmp(iptdim)
      REAL :: angaddztmp(iptdim),spinadxtmp(iptdim),spinadytmp(iptdim)
      REAL :: spinadztmp(iptdim)
      REAL,ALLOCATABLE :: mfracstmp(:,:)
      REAL :: tlastconvtmp(idim)
      INTEGER*4 :: numfailstmp(idim)
      INTEGER :: tmp_listpm(iptdim),dumpmols
      INTEGER,PARAMETER :: maxmols=35
      CHARACTER(len=30),ALLOCATABLE :: infiles(:)
      CHARACTER(len=30) :: outfile
      character*16 :: mollist(maxmols)
      character*20 :: molfile
      integer :: dead, sink, other

#IFDEF USEKROME
         usekrome = 1
! usekrome = 1 to use dumpmols not full dump
#ENDIF
      iindump = 1
      iout = 3
      nptmasstmp = 0
      nparttmp = 0
      istart = 1
      ptstart = 1
      ifulldump = 0
      imax = 1073741824
      nfullstep = 1
      iphasetmp = 0
#IFDEF USEKROME
      molnames = get_names()
      CALL readmolfile(maxmols,mollist,dumpmols)
#ENDIF
!  Find out how many files to read
      PRINT *, "Enter number of files"
      READ (*,*) numin
      ALLOCATE(infiles(numin))
      PRINT *, "Now enter", numin, "filenames"
      DO n = 1, numin
         READ (*,*) infiles(n)
      END DO
      ALLOCATE(mfracstmp(dumpmols,idim))
! Read files and append to temp
!      add nptmass to nptmasstmp

      DO n=1,numin
         PRINT *, "opening", infiles(n), 'istart=', istart
         OPEN (UNIT=iindump,FILE=infiles(n),FORM ='unformatted',
     &            ACTION='READ')
!     RECL=maxrec,         
         CALL rdump(iindump, ichkl, 0)
         CLOSE(iindump)
         PRINT *, "READ ", infiles(n)
! count running total of npart and nptmass as dumps are added
         nptmasstmp = nptmasstmp + nptmass
         nparttmp = nparttmp + npart
         PRINT *, "running totals: nptmass=", nptmasstmp,
     &          "npart=",nparttmp 
         DO i=1, npart
            IF (iphase(i) .EQ. -1) THEN
               dead = dead + 1
            ELSE IF (iphase(i) .GT. 0 .AND. iphase(i) .LE.6) THEN
               sink = sink + 1
               print *, i, " is a sink"
            ELSE IF (iphase(i) .NE. 0) THEN
               other = other+1
            END IF
         END DO
         PRINT *, "Dead =", dead, " sink=", sink, "nptmass=", nptmass, 
     &   "other =", other
         PRINT *, "sink list", (listpm(j), j=1, nptmass)

         
         iuniquetmp(istart:istart+npart) = iunique(1:npart) 
         istepstmp(istart:istart+npart) = isteps(1:npart)
         iphasetmp(istart:istart+npart) = iphase(1:npart)
         PRINT *, "done iphase"
         xyzmhtmp(:,istart:istart+npart) = xyzmh(:,1:npart)
         vxyzutmp(:,istart:istart+npart) = vxyzu(:,1:npart)
         rhotmp(istart:istart+npart) = rho(1:npart)
         gradhstmp(:,istart:istart+npart) = gradhs(:,1:npart)
         PRINT *, "done position,vel, rho"
         IF (iradtrans.EQ.idim) THEN
            ekcletmp(:,istart:istart+npart) = ekcle(:,1:npart)
            IF (idustRT.GT.0) THEN
               dust_tktmp(:,istart:istart+npart) = dust_tk(:,1:npart)
               chemistrytmp(:,istart:istart+npart) = 
     &                         chemistry(:,1:npart)
               heatingISRtmp(:,istart:istart+npart) = 
     &                         heatingISR(:,1:npart)
               h2fractmp(istart:istart+npart) = h2frac(1:npart)
               dh2dttmp(:,istart:istart+npart) = dh2dt(:,1:npart)
            END IF
         END IF
         PRINT *, "done RT"
!       MHD stuff 
         IF (imhd.EQ.idim) THEN
            Bevolxyztmp(:,istart:istart+npart) = Bevolxyz(:,1:npart)
            divcurlBtmp(:,istart:istart+npart) = divcurlB(:,1:npart)
            alphaMMtmp(:,istart:istart+npart) = alphaMM(:,1:npart)
            etarealtmp(istart:istart+npart) = etareal(1:npart)
            etaartificialtmp(istart:istart+npart) = 
     &            etaartificial(1:npart)
         END IF
         IF (nptmass.GT.0) THEN
!       sink particles, nptmasstmp
!            tmp_listpm(ptstart:ptstart+nptmass) = 
!     &            listpm(1:nptmass) + istart
            DO i=1,nptmass
               tmp_listpm(ptstart-1+i) = 
     &            listpm(i) + istart
            
               if (iphasetmp(tmp_listpm(ptstart-1+i)) .LT. 1) THEN
                  PRINT *, i, "iphasetmp < 1"
               end if
               if (iphase(listpm(i)) .LT. 1) THEN
                  PRINT *, i, "iphase < 1"
               end if
            END DO
            PRINT *, "listpm=", (tmp_listpm(ptstart-1+i),i=1,nptmass)
            spinxtmp(ptstart:ptstart+nptmass) = spinx(1:nptmass)
            spinytmp(ptstart:ptstart+nptmass) = spiny(1:nptmass)
            spinztmp(ptstart:ptstart+nptmass) = spinz(1:nptmass)
            angaddxtmp(ptstart:ptstart+nptmass) = angaddx(1:nptmass)
            angaddytmp(ptstart:ptstart+nptmass) = angaddy(1:nptmass)
            angaddztmp(ptstart:ptstart+nptmass) = angaddz(1:nptmass)
            spinadxtmp(ptstart:ptstart+nptmass) = spinadx(1:nptmass)
            spinadytmp(ptstart:ptstart+nptmass) = spinady(1:nptmass)
            spinadztmp(ptstart:ptstart+nptmass) = spinadz(1:nptmass)
            PRINT *, "Done sinks"
         END  IF

#ifdef USEKROME
         PRINT *, "starting chemistry"
!       krome chemistry
         mfracstmp(1:dumpmols,istart:istart+npart) = 
     &      mfracs(1:dumpmols,1:npart)
         PRINT *, "done mfracs"
         tlastconvtmp(istart:istart+npart) = tlastconv(1:npart)
         numfailstmp(istart:istart+npart) = numfails(1:npart)
         PRINT *, "Done chemistry"
#endif

      istart = istart + npart + 1
      ptstart = ptstart + nptmass + 1
      END DO
! Overwrite arrays from temp and call wdump
      nptmass = nptmasstmp
      npart = nparttmp
      iuniquemax = npart

      iunique = iuniquetmp
      isteps = istepstmp
      iphase(:) = 0
      iphase = iphasetmp

      xyzmh = xyzmhtmp
      vxyzu = vxyzutmp
      rho = rhotmp
      gradhs = gradhstmp

      IF (iradtrans.EQ.idim) THEN
            ekcle = ekcletmp
            IF (idustRT.GT.0) THEN
               dust_tk = dust_tktmp
               chemistry = chemistrytmp
               heatingISR = heatingISRtmp
               h2frac = h2fractmp
               dh2dt = dh2dttmp
            END IF
         END IF
!       MHD stuff 
         IF (imhd.EQ.idim) THEN
            Bevolxyz = Bevolxyztmp
            divcurlB = divcurlBtmp
            alphaMM = alphaMMtmp
            etareal = etarealtmp
            etaartificial = etaartificialtmp
         END IF
         IF (nptmass.GT.0) THEN
!       sink particles, nptmasstmp
            listpm = tmp_listpm
            spinx = spinxtmp
            spiny = spinytmp
            spinz = spinztmp
            angaddx = angaddxtmp
            angaddy = angaddytmp
            angaddz = angaddztmp
            spinadx = spinadxtmp
            spinady = spinadytmp
            spinadz = spinadztmp
            PRINT *, "Done sinks"
         END  IF

#ifdef USEKROME
         PRINT *, "starting chemistry"
!       krome chemistry
         mfracs = 0d0
         DO j=1,dumpmols
            index = krome_get_index(adjustl(mollist(j)))
            mfracs(index,:) = mfracstmp(j,:)
         END DO
         PRINT *, "done mfracs"
         tlastconv = tlastconvtmp
         numfails = numfailstmp
         PRINT *, "Done chemistry"
#endif
         PRINT *, (rho(i), i=1002,1050)
      CALL getoutfile(infiles(1),outfile)
      DO i=1,npart
         isort(i) = i
      END DO
      OPEN(UNIT=iout,FILE=outfile,FORM ='unformatted',
     &            ACTION='WRITE',STATUS='REPLACE')
      PRINT *, "usekrome=", usekrome, iradtrans, idustRT

      usekrome = 1
      CALL wdump(iout)
      dead = 0
      sink = 0
      other = 0
      DO i=1, npart
         IF (iphase(i) .EQ. -1) THEN
            dead = dead + 1
         ELSE IF (iphase(i) .GT. 0 .AND. iphase(i) .LE.6) THEN
            sink = sink + 1
         ELSE IF (iphase(i) .NE. 0) THEN
            other = other+1
         END IF
      END DO
      PRINT *, "Dead =", dead, " sink=", sink, "nptmass=", nptmass, 
     &  "other =", other

      CLOSE(iout)

      PRINT *, "deallocating arrays"
      DEALLOCATE(infiles,mfracstmp)
      END PROGRAM mergedump
!******************************************

      SUBROUTINE quit(i)
       STOP
      END SUBROUTINE quit

      SUBROUTINE endrun
       CALL quit(0)
      END SUBROUTINE endrun

      SUBROUTINE getoutfile(inname,outfile)
        implicit none
        CHARACTER(len=30) :: inname,outfile
        INTEGER :: ind
        ind = index(inname,"_")
        outfile = inname(:ind) //"total"//inname(ind+7:)
        PRINT *, "WRITING TO ", outfile

      END SUBROUTINE getoutfile
