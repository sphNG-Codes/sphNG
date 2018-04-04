      PROGRAM krome_sph
c*******************************************************************
c     calls Krome for an sphNG dump file                           *
c     use "params" text file for list of dumps:                    *
c       eg. N/R New run/Restart from chem dump                     *
c           infofilename                                           *
c           no. input dumps  1                                     *
c           list of files  FC00186../"                             *
c           type of chem timestep: (F)ixed/read from(D)ump         *
c           fixed only [timestep (years)]                          *   
c                      [ number of timesteps]                      *
c           D only     [init timestep (yrs) ]                      *
c           small dump frequency (every X timesteps)               *
c           full dump frequency (every X timesteps)                *
c                                                                  *
c     To input initial abundances (mass fractions), the file       *
c       "abundances.dat" must be present in format "H2 0.5 ... "   *
c    OR "init_mfracs.dat" containing mass fractions "H2 0.73 ..."  * 
c                                                                  *
c     N.B. reactions_verbatim.dat  must be copied                  *
c     from krome/build to your running directory and add           *
c     path to krome/build to LD_LIBRARY_PATH                       *
c                                                                  *
c*******************************************************************
      INCLUDE 'COMMONS/krome_mods'
      implicit none
#ifdef KROMEMPI
      INCLUDE 'mpif.h'
#endif
      INCLUDE 'idim'
      INCLUDE 'COMMONS/kromevar'
      INCLUDE 'igrape'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/raddust'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/interstellar'
      
      INTEGER  maxdumps,istart,iend
      PARAMETER (maxdumps=400)
      REAL denscut
      PARAMETER (denscut=5.d-5)
      REAL*8 dust2gas
      PARAMETER (dust2gas=0.01d0)
      INTEGER i,maxrec,iout,ichkl,parunit,numdumps,idump,j 
      REAL*8 dt,spy,Tgas,rhoi,x(krome_nmols),local_dt,Av
      REAL*8 init_x(krome_nmols),ktime,init_dt !ktime=krometime
      CHARACTER*40 infile,filelist(maxdumps),infofile
      CHARACTER*40 outfilename,tempchar,tmpofilename
      INTEGER n,k,niter,sdumpfreq,nskipped,fdumpfreq,inew(idim2)
      INTEGER ioutinfo,iindump,itime,kstart,iunique_orig(idim2)
      LOGICAL init_from_dump,havefile,restart,fixed_dt,found
      CHARACTER*6 ctime
      CHARACTER*40 fdumpname
      REAL*8 mfracs_orig(krome_nmols,idim)
      CHARACTER*1 cdump
      PARAMETER (spy = 3.1536d7) !seconds per year
      INTEGER iphasecount
#ifdef KROMEMPI
      INTEGER master,status(MPI_STATUS_SIZE),bufflen,inumber
      PARAMETER (bufflen = krome_nmols*idim/2)
      PARAMETER (MASTER=0)
      INTEGER ierr,nproc,rank,lchunksize,chunksize,buffnum,count
      INTEGER totskipped,totipcount
      REAL*8 buffer(bufflen)
      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
#endif 
     
      WRITE(*,*) "Beginning with:",krome_nmols, "molecules"
      imax = 1073741824
      ifulldump = 0
      nfullstep = 1
      parunit = 1
      iout = 12
      iindump = 13
      ioutinfo = 14
      maxrec = 100 * idim
      ktime = 0d0 !seconds
c     usekrome is an INT looked at by rdump and wdump
c     usekrome = 0 read normal non-chem dump
c     usekrome = 1 write chemistry to dump from dumpmols (no mhd)
c     usekrome = 2 read chemistry from dumpfile/ write full chem dump
c     usekrome = 3 read chemistry from mhd without mhd data /
c                  write chemistry without mhd data

c     write small dump every sdumpfreq iterations
      init_from_dump = .False.
      nskipped=0
      
      WRITE(*,*) "Getting molecule names from krome"
      molnames = get_names()
      WRITE(*,*) "Found molecule names"
      CALL unit
c-- Read list of dumpfiles
      CALL readparams(parunit,maxdumps,numdumps,filelist,dt,niter,
     &     sdumpfreq,fdumpfreq,restart,infofile,fixed_dt,init_dt)
      WRITE(*,*) "Finished reading parameters"
      CALL checkfilespresent(maxdumps,filelist,infofile,restart)

      IF (restart) THEN
               OPEN (UNIT=ioutinfo,FILE=infofile,FORM ='formatted',
     &        STATUS="old", POSITION="append")         
               WRITE(ioutinfo,"(/,A,/)") "RESTARTED"
               CALL readtime(ktime,k,filelist(1),infofile,spy)
      ELSE
               OPEN (UNIT=ioutinfo,FILE=infofile,FORM ='formatted',
     &        STATUS="replace")         
               write(ioutinfo,*) "PARAMETERS:"
               IF (imhd.EQ.idim) then
                  write(ioutinfo,*) "MHD run"
               ELSE IF (imhd.EQ.1) then
                  write(ioutinfo,*) "Not MHD"
               ELSE
                  write(ioutinfo,*) "Not MHD but imhd .NE. 1"
               ENDIF
               IF (iradtrans.EQ.idim) then
                  write(ioutinfo,*) "RT run"
               ENDIF
               IF (idustRT.GT.0) then
                  write(ioutinfo,*) "Dust RT"
               ENDIF
               write(ioutinfo,*) "udist", udist
       END IF

      IF (fixed_dt) THEN
         WRITE(ioutinfo,*) "Files:", (filelist(i),i=1,numdumps) 
         WRITE(ioutinfo,"(A,F9.1,X,A,I0)") "dt=",dt,"niter=",niter
      ELSE
         WRITE(ioutinfo,"(A)") "dt set from dumps"
         WRITE(ioutinfo,"(A,F9.1)") "Initial timestep:", init_dt
      END IF
      varmhd = 'Brho'
      WRITE(*,*) "Initialising krome"
      CALL krome_init()

c-----
c----- Beginning of loop over filelist
      DO idump=1, numdumps 
        WRITE(*,"(A,X,A)") "FILE:", filelist(idump)
c---  initialise whole mfracs array to 0      
        mfracs = 0d0

c        WRITE(*,*) (mfracs_orig(krome_idx_H2,i),i=50,150)
        usekrome = 0
c ---- Determine kind of dump to read and set chem dumpfile name
        tempchar = filelist(idump)
        n = index(tempchar," ")
        IF (index(filelist(idump),"full") .GT. 0) THEN
           WRITE(*,*) "Reading full chem dump file"
           init_from_dump = .True.
           usekrome = 2
c          usekrome = 3
           n = index(tempchar,"_t")
           outfilename = tempchar(1:n+1)
        ELSE IF ( n .GT. 0) THEN
           outfilename = tempchar(1:n-1)//"_chem_t"
        ELSE
           outfilename = tempchar//"_chem_t"
        ENDIF
      
c ---- Read dumpfile      
        OPEN (UNIT=iindump,FILE=filelist(idump),FORM ='unformatted',
     &        RECL=maxrec,SHARED,ACTION='read')
        WRITE (*,*) 'File open; usekrome =', usekrome
        CALL rdump(iindump, ichkl, 0)
        CLOSE(iindump)
        WRITE(*,"(A,X,I0,X,A)") "Read", npart,"particles"
        WRITE(*,"(A7,X,A,X,E9.2,A)") filelist(idump),"sphNG time= ",
     &       gt*utime/spy,"years"

        DO i=1,npart
           isort(i) = i
c     --- Now wipe dead particles' iuniques so live particles aren't
c     --- matched up with dead ones by accident
           IF (iphase(i) .NE. 0) THEN
              iunique(i) = -1
           END IF
        END DO

#ifdef KROMEMPI
c           npart = 200000 ! for testing
           chunksize = npart / nproc
           istart = rank * chunksize + 1
           IF (rank .EQ. nproc-1) THEN
              iend = npart
              lchunksize = chunksize + MOD(npart,nproc)
           ELSE
              iend = istart + chunksize - 1
              lchunksize = chunksize
           END IF
           PRINT *, "This is rank ", rank, "doing",istart," to", iend
#else 
           istart = 1
           iend = npart
#endif

        IF (idump .EQ. 1) THEN
c        if idump > 1 iunique_orig is previous dump      
c     ---  i.e. iunique_orig(i) = iunique(i) for all i
           iunique_orig = iunique
           IF (restart .AND. .NOT. fixed_dt) THEN
              mfracs_orig = mfracs
             CYCLE
c     ---    Don't need to run chem as read from dump
           END IF
        ELSE
c ---   find order of current particles in previous dump
           PRINT *, "Matching particles with previous dump..."
C$OMP PARALLEL DO SCHEDULE(static,1) default(none)
C$OMP& PRIVATE(i,j,found) SHARED(npart,iunique,iunique_orig,inew,
C$OMP& istart,iend,iphase) 
           DO i=1,npart
              IF (iphase(i) .NE. 0 ) THEN
                 inew(i) = -1 
                 CYCLE
              END IF
              found = .FALSE.
              DO j=1,npart
                 IF (iunique(i) .EQ. iunique_orig(j) ) THEN
                    inew(i) = j ! tells you position in old list
                    found = .TRUE.
                 END IF
              END DO
              IF (.NOT. found) THEN
                 PRINT *, "ERROR: iunique not found for particle", i
              END IF
           END DO
C$OMP END PARALLEL DO

           write(cdump,'(I1)') idump
           print *, "done iuniques for idump=", idump
          print *, (mfracs_orig(krome_idx_H2,inew(i)),i=1,50)
           OPEN(UNIT=32,FILE="newunique"//cdump,FORM='formatted',
     &        STATUS='replace')
           WRITE(32,'(2I0)') (iunique(i),iphase(i),i=1,40)
           CLOSE(32)
           OPEN(UNIT=31,FILE="oldunique"//cdump,FORM='formatted',
     &          STATUS='replace')
           WRITE(31,'(I0)') (iunique_orig(inew(i)), i=1,40)
           CLOSE(31)
           
           PRINT *, "OLD H2:", mfracs_orig(krome_idx_H2,inew(23)),
     &          iunique_orig(inew(23)),iunique(23),iphase(23)

           PRINT *, "particle 500000", mfracs_orig(krome_idx_H2,
     &          inew(500000))
        END IF

c Test order
c      DO i=1,npart
c        iunique_orig(i) = iunique(i)
c      END DO
cc
c **** TEST
c      outfilename = "rhocheckdump"
c      WRITE(*,*) "writing to file: ", outfilename
c      OPEN (UNIT=iout, FILE = outfilename, FORM= 'unformatted', 
c     &     RECL=maxrec, STATUS='replace')
c      usekrome = 1
c      CALL wdump(iout)
c      STOP
c*****
c      CYCLE

c      write(*,*) "user crate=", krome_get_user_crate()
        IF (restart .AND. fixed_dt) THEN
           kstart = k + 1
        ELSE
           kstart = 1
        END IF
c---- get init_x = init abundance rel to H nuclei     
        IF (idump.EQ.1 .AND..NOT. init_from_dump ) THEN
           CALL get_init_abunds(krome_nmols,init_x)
           WRITE(ioutinfo,*) "Initial Mass Fractions:"
           DO j=1,krome_nmols
              WRITE(ioutinfo,"(A12,E16.5)") molnames(j),init_x(j) 
           END DO
        ELSE IF (idump.EQ.1 .AND. init_from_dump) THEN
           WRITE(ioutinfo,*) "Initialising mass fractions from dump"
        ENDIF

        WRITE(ioutinfo,"(/,A,6X,A,6X,A/)") "File", "Time", "Iter No."
        PRINT *, "idump= ", idump

c---  Set timestep local_dt used to call KROME
c---  Local_dt in SECONDS
        IF (fixed_dt) THEN
           local_dt = dt * spy
        ELSE if (idump.EQ.1 .AND. .NOT. restart) THEN
           local_dt = init_dt * spy
        ELSE
           local_dt = (gt * utime) - ktime 
        END IF

c--   Loop to run chemistry more than once on dumpfile, dumping
c--   after each run
        DO k=kstart, niter
           nskipped = 0
           WRITE (*,104) "Iteration number: ",k,"of",niter,"for file",
     &          filelist(idump)
 104       FORMAT (A,X,I0,X,A,X,I0,A,X,A7)
           
           WRITE (*,*) "Beginning at krometime t=", ktime/spy,
     &          "timestep = ", local_dt/spy

           iphasecount = 0
c     Main loop over particles to initialise mass fracs and call krome
c     
C$OMP PARALLEL DO SCHEDULE(runtime) default(none)
C$OMP& shared(rho,ekcle,mfracs,vxyzu,udens,npart,iunique,iunique_orig,
C$OMP& iphase,k,init_from_dump,dust_tk,heatingISR,init_x,nskipped,
C$OMP& idump,mfracs_orig,inew,local_dt,istart,iend,iphasecount)
C$OMP& private(i,j,Tgas,rhoi,x,Av)
 
           DO i=istart, iend
              IF (MOD(i,10000) .EQ. 0) THEN
                 PRINT *, "particle", i,"of",npart,"nskipped=",nskipped
              END IF

              IF (iphase(i) .NE. 0) THEN
c                 PRINT *, "skipping iunique",iunique(i),"due to iphase",
c     &                iphase(i)
c$OMP ATOMIC
                 iphasecount = iphasecount + 1
                 CYCLE
              END IF

              IF (idump.GT.1) THEN
c ----     get abundances calculated for prev dump
c           print *, "mforig",mfracs_orig(krome_idx_H2,inew(i)) 
                 DO j=1,krome_nmols
                    x(j) = mfracs_orig(j,inew(i))
                 END DO
c           PRINT *, "from prev calc, xH2=", x(krome_idx_H2)
              ELSE IF (init_from_dump .OR. k.GT.1) THEN
c     ----     Use mass fracs from prev krome run(same dump) or full_dump
                 DO j=1, krome_nmols
                    x(j) = mfracs(j,i)
                 END DO
c           PRINT *, "x = mfracs"
              ELSE IF ( k.EQ.1 ) THEN   
                 DO j=1, krome_nmols 
                    x(j) = init_x(j)
                 END DO
c           print *, "init,xh2=", x(krome_idx_H2) 
              ELSE
                 print *, "ERROR- don't know init mfracs!"
                 STOP
              END IF
         
              IF (imhd.EQ.idim .AND. rho(i).LT.denscut) THEN
                 DO j=1, krome_nmols
                    mfracs(j,i) = 0d0
                 END DO
                 print *, "skipping due to imhd"
                 CYCLE
              END IF
              
              IF (ekcle(3,i) .LT. tiny) then
                 write(*,*) ekcle(3,i), "skipping U/T = 0"
                 CYCLE
              END IF
                  
              Tgas = DBLE(vxyzu(4,i)/ekcle(3,i))

              rhoi = DBLE(rho(i) * udens)
              IF (rhoi .GT. 4d-11) THEN
                 DO j=1, krome_nmols
                    mfracs(j,i) = 0d0
                 END DO
c            print *, "chem dens cut" 
c$OMP ATOMIC
                 nskipped = nskipped + 1
                 CYCLE
              ELSE IF (x(krome_idx_H2) .LT. tiny) THEN
c             ELSE IF (ALL(x.eq.0.d0)) THEN
c-------- If H2 abundance=zero then it must have been set so            
c$OMP ATOMIC
                 nskipped = nskipped + 1

                 PRINT *, "skipping H2 = 0,iphase=",iphase(i),"i=", i,
     &                "i prev=",inew(i),iunique_orig(inew(i)),iunique(i)
                 CYCLE
              END IF

              IF (idustRT .GT. 0 ) then
                 CALL krome_set_user_Tdust(dust_tk(1,i))
              ELSE 
                 CALL krome_set_user_Tdust(Tgas)
              ENDIF

              IF (heatingISR(3,i) .GT. 0.d0) THEN
                 Av = -LOG(heatingISR(3,i))
              ELSE
                 Av = 86d0
              END IF
              CALL krome_set_user_Av(Av) 
              CALL krome_set_user_crate(1.3d-17)         
              CALL krome_set_user_rhogas(rhoi)

c         WRITE(*,*) "calling krome:"
c         WRITE(*,*) "rho, T, AV,dt", rhoi,Tgas,Av, local_dt/spy
c------- CALLING KROME
              CALL krome(x,rhoi,Tgas,local_dt)
cfor testing
c              DO j=1, krome_nmols
c                 x(j) = i*10.d0
c              END DO
c-------Update mfracs array
              DO j=1, krome_nmols
                 mfracs(j,i) = x(j)
              END DO

c         print *, "Graintot", 
c     &        x(krome_idx_GRAIN0)+x(krome_idx_GRAINk),Tgas,rhoi
           END DO
C$OMP END PARALLEL DO
c     --- End of loop over i particles

c     ---  Update time
           IF ( (idump.GT.1.AND..NOT.restart) .OR.restart ) THEN
              ktime = ktime + local_dt
           ELSE
              ktime = 0d0
           END IF

#ifdef KROMEMPI           
           IF (rank .EQ. MASTER) THEN
              DO n=1, nproc-1
                 PRINT *, "MASTER waiting to receive from process", n
                 istart = n * chunksize + 1
                 IF (n .LT. (nproc-1)) THEN
                    count = chunksize * krome_nmols
                    iend = istart + chunksize -1
                 ELSE IF (n .EQ. (nproc-1)) THEN 
                    iend = npart
                    count = (chunksize + MOD(npart,nproc)) * krome_nmols
                 END IF
                DO i=istart,iend
                   CALL MPI_RECV(mfracs(1,i),krome_nmols,MPI_REAL8,n,1,
     &                  MPI_COMM_WORLD,status,ierr)
                   CALL MPI_GET_COUNT(status,MPI_REAL8,inumber,ierr)
                   IF (inumber .ne. krome_nmols) THEN
                      WRITE (*,*) n,':ERROR partial chem transferred',
     &                     inumber
                      CALL quit(1)
                   ENDIF
                END DO
                 PRINT *, "Received", istart," to",iend,"from rank",n
                 PRINT *, "CO particle at iend",iend,
     &                mfracs(krome_idx_CO,iend)
              END DO
          ELSE

             DO i=istart, iend
                CALL MPI_SEND(mfracs(1,i),krome_nmols,
     &             MPI_REAL8,MASTER,1,MPI_COMM_WORLD,ierr)
             END DO
             PRINT *, "Rank", rank, " sending",istart,
     &             mfracs(krome_idx_CO,istart)
          END IF
           
          CALL MPI_Reduce(nskipped,totskipped,1,MPI_INTEGER,MPI_SUM,
     &          MASTER,MPI_COMM_WORLD,ierr)

          CALL MPI_Reduce(iphasecount,totipcount,1,MPI_INTEGER,
     &         MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)

          IF (rank .EQ. MASTER) THEN
             PRINT *, "Total skipped = ", totskipped,"iphase skipped=",
     &            totipcount
#endif

c     ---  Create dumpfile name
           itime = INT(ktime/spy)
           write(ctime,"(I6.6)") itime
           print *, "itime=", itime
           n= index(outfilename," ")
           tmpofilename = outfilename(1:n-1)//ctime
           print *, "tmpname ", tmpofilename
c     ---  Do small dump every <sdumpfreq> iterations
           IF (MOD(k,sdumpfreq).EQ.0 .OR. 
     &          MOD(idump,sdumpfreq).EQ.0) THEN
              WRITE(*,*) "writing to file: ", tmpofilename
              OPEN (UNIT=iout, FILE = tmpofilename, FORM= 'unformatted', 
     &             RECL=maxrec, STATUS='replace')
              usekrome = 1
              CALL wdump(iout)
              WRITE(*,*) "usekrome", usekrome
              CLOSE(iout)
              WRITE(ioutinfo,"(A,X,F14.1,X,I0)") tmpofilename,ktime,k 
           END IF

c     --- Do full dump every <fdumpfreq> iterations
           IF ((idump.EQ.1).OR.(niter.GT.1 .AND. k.EQ.niter)
     &          .OR.(MOD(k,fdumpfreq).EQ.0)
     &          .OR. (MOD(idump,fdumpfreq).EQ.0)
     &          .OR. (idump.EQ.numdumps) ) THEN
              IF (imhd.EQ. idim) THEN
                 usekrome = 2
c      usekrome=2 dumps mhd info where present
              ELSE
                 usekrome = 2
              ENDIF
              n = index(tmpofilename," ") 
              fdumpname = trim(tmpofilename)//"full"
              OPEN (UNIT=iout, FILE=fdumpname, FORM='unformatted'
     &             ,RECL=maxrec,STATUS='replace')
              CALL wdump(iout)
              PRINT *, "Written full chemistry dump: ", fdumpname
              CLOSE(iout)
              WRITE(ioutinfo,"(A,X,F14.1,X,I0)") fdumpname,ktime/spy,k 
              usekrome = 1
           END IF
#ifdef KROMEMPI
           END IF
#endif
           PRINT *, "End of timestep", k, "Nskipped=,", nskipped
        END DO
c-- End of multiple krome runs loop

        IF (idump.LT.numdumps) THEN
#ifdef KROMEMPI
           IF (rank .EQ. MASTER) THEN
#endif
              Print *, "updating mfracs_orig"
              PRINT *, "particle 23", mfracs(krome_idx_H2,238126)     
c ---- update iunique_orig and mfracs_orig
              iunique_orig(:) = iunique(:)
              mfracs_orig(:,:) = mfracs(:,:)
#ifdef KROMEMPI
              PRINT *, "orig,rank",rank,
     &             mfracs_orig(krome_idx_H2,npart-1)
           ELSE IF (rank .GT. MASTER) THEN
              iunique_orig(:) = iunique(:)
              mfracs_orig = 0d0
           END IF
c     --- all ranks onwards           
           PRINT *, "NPART=",npart
           DO i=1,npart
              IF ((i.eq.npart) .or.
     &             (MOD(i,10000) .eq. 0) ) THEN
                 PRINT *,"Before: mfracs This is rank",rank,"particle:",
     &                i, "H2:",mfracs(krome_idx_H2,i)
              END IF
              CALL MPI_Bcast(mfracs_orig(1,i),krome_nmols,
     &             MPI_REAL8,MASTER,MPI_COMM_WORLD,
     &             ierr)
c     PRINT *, "rank", rank,"completed broadcast"
              IF ((i .eq. npart) .OR.
     &             ((rank .eq. 1) .AND. (MOD(i,10000).eq.0))) THEN
                 Print *, "After: orig This is rank",rank,"particle:",
     &                i,"H2:",mfracs_orig(krome_idx_H2,i)
              END IF
           END DO
#endif
        ELSE 
           WRITE(*,*) "END"
        END IF
      END DO
c-- End of loop over filelist

      CLOSE(ioutinfo)
      
#ifdef KROMEMPI
      CALL mpi_finalize(ierr)
#endif
      CALL quit(0)
      END PROGRAM krome_sph

c ****************************************************************

      SUBROUTINE quit(i)
      STOP
      END 

c ****************************************************************

      SUBROUTINE endrun
      CALL quit(0)
      END

c*****************************************************************
c-----Reads list of dumpfiles and other parameters
      SUBROUTINE readparams(parunit,maxdumps,numdumps,filelist,dt,niter,
     & sdumpfreq,fdumpfreq,restart,infofile,fixed_dt,init_dt)
      INTEGER, INTENT(IN) :: parunit,maxdumps
      INTEGER,INTENT(OUT) ::  numdumps,niter,sdumpfreq,fdumpfreq
      CHARACTER*40, INTENT(OUT) :: filelist(maxdumps),infofile
      REAL*8, INTENT(OUT) :: dt,init_dt
      LOGICAL,INTENT(OUT) :: restart,fixed_dt
      CHARACTER(len=1) :: tmp

      OPEN(UNIT=parunit, FILE="params", FORM="FORMATTED",
     &     STATUS="OLD") 
      READ(parunit,*) tmp
      IF (tmp.EQ."N" .OR. tmp.EQ."S") THEN
         restart = .False.
         PRINT *, "New run"
      ELSE IF (tmp.EQ."R") THEN
         restart = .True.
         PRINT *, "Restarting"
      ELSE
         PRINT *, "ERROR no restart option specified in params"
         CALL QUIT(1)
      END IF
      READ(parunit,*) infofile
      READ(parunit,*) numdumps
      PRINT *, "expecting", numdumps, "dump files"
      IF (numdumps.GT.maxdumps) THEN
         PRINT *, "Error numdumps > maxdumps"
         CALL QUIT(1)
      END IF
      DO i=1, numdumps
         PRINT *, i
        READ(parunit,*) filelist(i)
        PRINT *, filelist(i)
      END DO
      IF ( i < numdumps) THEN
         PRINT *, "Error filenames < numdumps"
         CALL QUIT(1)
      END IF
      READ(parunit,*) tmp
      IF (tmp.EQ."F") THEN
         fixed_dt = .True.
      ELSE IF (tmp.EQ."D") THEN
         fixed_dt = .False.
      ELSE
         PRINT *, "ERROR chem timestep option not found"
         CALL QUIT(1)
      END IF
      IF (fixed_dt) THEN
         READ(parunit,*) dt
         READ(parunit,*) niter
         init_dt = 0.
      ELSE
         READ(parunit,*) init_dt
         dt = 0.
         niter = 1
      END IF
      READ(parunit,*) sdumpfreq
      READ(parunit,*) fdumpfreq
      CLOSE(parunit)

      END SUBROUTINE 
c ****************************************************************
      
      SUBROUTINE checkfilespresent(maxdumps,filelist,infofile,restart)
      INTEGER maxdumps,i
      CHARACTER*40 filelist(maxdumps),infofile
      LOGICAL filepresent,anymissing,restart
      
      anymissing = .FALSE.
      CALL ispresent("dumpmols.txt",filepresent,anymissing)
      CALL ispresent("reactions_verbatim.dat",filepresent,anymissing)
      CALL ispresent("init_mfracs.dat",filepresent,anymissing)
      IF (restart) THEN
         CALL ispresent(infofile,filepresent,anymissing)
      END IF
      DO i=1,maxdumps
         IF (filelist(i).EQ."") EXIT
         CALL ispresent(filelist(i),filepresent,anymissing)
      END DO
      CALL ispresent("dumpmols.txt",filepresent,anymissing)
      IF (anymissing) THEN
         WRITE (*,*) "MISSING FILES ERROR. EXITING."
         CALL QUIT(0)
      ELSE
         WRITE (*,*) "All files present"
      END IF
      END SUBROUTINE

c ****************************************************************
      SUBROUTINE ispresent(filename,filepresent,anymissing)
      CHARACTER*40 filename
      LOGICAL filepresent, anymissing

      INQUIRE(FILE=filename,EXIST=filepresent)
      IF (.NOT. filepresent) THEN
         anymissing = .True.
         WRITE(*,*), filename, "MISSING"
      END IF
      END SUBROUTINE

c ****************************************************************
      SUBROUTINE readmols(molfile)
      INCLUDE 'COMMONS/krome_mods'
      implicit none
      INCLUDE 'idim'
      INCLUDE 'COMMONS/kromevar'
      
      CHARACTER*11 molfile
      INTEGER unit, molcount,endind,j, maxsp
      CHARACTER line*100, molid*25
      unit = 23
      molcount = 0
      maxsp = 200

      WRITE(*,*) "opening"
      OPEN(UNIT=unit,FILE=molfile, FORM='formatted',
     &     STATUS='old')

      DO j=1,maxsp
         READ(unit,*) line
         IF (line(1:6).EQ."krome_") THEN
          molcount = molcount + 1
          molnames(molcount) = line(11:)
        ENDIF
c        write(*,*) "molcount=", molcount,krome_nmols
         IF (molcount.EQ.krome_nmols) THEN
            EXIT
         ENDIF
      END DO
      CLOSE(unit)

      WRITE(*,*) "found", krome_nmols, molcount,"species"
      END

c *********************************************
      SUBROUTINE ofilename(outfilename,k,n)
      implicit none
      CHARACTER*40 outfilename, tempname
      INTEGER k,n,m,i
      CHARACTER*3 tempchar,newtempchar

      write(tempchar,'(I3)')  k
      write(*,*) "tempchar=", tempchar
      DO i = 1, 3
         m = INDEX(tempchar," ")
         IF (m.EQ.1) THEN
            newtempchar = "0"//tempchar(m+1:3)
            write(*,*) "newtchar=", newtempchar
            tempchar = newtempchar
         ELSE IF (m.GT.1) THEN
            newtempchar = tempchar(1:m-1)//"0"//tempchar(m+1:3)
            tempchar = newtempchar
         ENDIF
      END DO
      write (*,*) "k=", k, "tempchar=", tempchar
      write(*,*) "tempchar=", tempchar
      n = index(outfilename,"m") 
      tempname = outfilename
      IF (n .GT. 0) then
        outfilename = tempname(1:n)//tempchar
      ELSE
        outfilename = tempname//tempchar
      ENDIF

      END 

c**********************************************
c abunds rel to H nuclei
      SUBROUTINE get_init_abunds(nmols,init_x)
        INCLUDE "COMMONS/krome_mods"
        IMPLICIT NONE
        INTEGER nmols,i
        REAL*8 init_x(nmols)
        CHARACTER*40 abundfile,mfracfile
        LOGICAL haveAbundfile,haveMfracfile
        abundfile = "abundances.dat"
        mfracfile = "init_mfracs.dat"
        DO i=1,nmols
           init_x(i) = 0.d0
        END DO
        INQUIRE(FILE=abundfile,EXIST=haveAbundfile)
        INQUIRE(FILE=mfracfile,EXIST=haveMfracfile)
        IF (haveMfracfile) THEN
           WRITE (*,*) "init_mfracs.dat found. Reading..."
           CALL read_mfracs(nmols,init_x,mfracfile)
        ELSE IF (haveAbundfile) THEN
           WRITE(*,*) "BEWARE! abunds read needs updating"
           WRITE(*,*) "abundances.dat found. Reading..."
           CALL read_abunds(nmols,init_x,abundfile)
        ELSE
           WRITE(*,*) "Using default abundances"
          init_x(krome_idx_H) = 0.73
          init_x(krome_idx_He) = 0.279
          init_x(krome_idx_Cj) = 7.13d-4
          init_x(krome_idx_O) = 1.42d-3
c          init_x(krome_idx_Si) = 6.78d-5
        ENDIF
      END SUBROUTINE get_init_abunds

c *********************************************
c **** Read abundances relative to TOTAL H nuclei
      SUBROUTINE read_abunds(nmols,init_x,abundfile)
         INCLUDE "COMMONS/krome_mods"
         IMPLICIT NONE
         INTEGER nmols
         REAL*8 init_x(nmols),abunds(nmols),abund,masses(nmols)
         CHARACTER*20 abundfile,name,molnames(nmols)
         CHARACTER*30 kromename
         INTEGER iunit,iostat,index,j
c
c     ADD GRAIn and E initialization to this subrtn.
c
         iunit=16
         iostat=0
         molnames = krome_get_names()
         masses = krome_get_mass()
         WRITE(*,*) "No mols:", nmols
c         WRITE(*,*) molnames
         OPEN(UNIT=iunit,FILE=abundfile,FORM='formatted',
     &     STATUS='old')
         DO
            READ(iunit,*,IOSTAT=iostat) name, abund  
            IF (name(1:1) == "!") THEN
               CYCLE
            END IF
c            WRITE(*,*) "reading:", name,abund
            IF (iostat > 0) THEN
               WRITE(*,*) "Error reading abundances"
               CALL quit(1)
            ELSE IF (iostat < 0) THEN
               EXIT
            ELSE
c               WRITE(*,*) name, abund
               index = krome_get_index(name)
               abunds(index) = abund
c               init_x(index) = abund
            END IF
         END DO
         CLOSE(iunit)
c       Convert rel abundances to mass fracs
         DO j=1,nmols
            init_x(j) = abunds(j) * masses(j)
         END DO
         DO j=1,nmols
            init_x(j) = init_x(j)/SUM(init_x)
         END DO
      END SUBROUTINE read_abunds

c-------------------------------------------------
      SUBROUTINE read_mfracs(nmols,init_x,mfracfile)
         INCLUDE "COMMONS/krome_mods"
         IMPLICIT NONE
         INTEGER nmols
         REAL*8 init_x(nmols),frac
         CHARACTER*20 mfracfile,name,molnames(nmols)
         CHARACTER*30 kromename
         INTEGER iunit,iostat,index,j

         iunit=17
         iostat=0
         molnames = krome_get_names()
         
         WRITE(*,*) "No mols:", nmols
c         WRITE(*,*) molnames
         OPEN(UNIT=iunit,FILE=mfracfile,FORM='formatted',
     &     STATUS='old')
         DO
            READ(iunit,*,IOSTAT=iostat) name, frac  
            IF (name(1:1) == "!") THEN
               CYCLE
            END IF
c            WRITE(*,*) "reading:", name,abund
            IF (iostat > 0) THEN
               WRITE(*,*) "Error reading abundances"
               CALL quit(1)
            ELSE IF (iostat < 0) THEN
               EXIT
            ELSE
c               WRITE(*,*) name, abund
               index = krome_get_index(name)
               init_x(index) = frac
            END IF
         END DO
         CLOSE(iunit)
       END SUBROUTINE read_mfracs

c************************************************************
c ---- reads ktime and iteration number reached previously
      SUBROUTINE readtime(ktime,k,filename,infofile,spy)
        IMPLICIT NONE
        REAL*8, INTENT(OUT) :: ktime
        REAL*8 :: spy
        INTEGER :: k,iread
        INTEGER :: reason,n
        CHARACTER*120 :: line
        CHARACTER*40 :: filename,infofile,name
        LOGICAL :: found

        iread = 16
        found = .FALSE.
        n = INDEX(filename," ")
        OPEN(UNIT=iread,FILE=infofile,FORM='formatted',
     &     STATUS='old')
        PRINT *, "Reading krome time from info file..."
        DO
           READ(iread,'(A)',IOSTAT=reason) line
           IF (reason .GT. 0) THEN
              print *, "ERROR reading time"
           ELSE IF (reason .LT.0) THEN
             EXIT
           ELSE
              IF (INDEX(line,filename(1:n-1)) .GT. 0) THEN
                 found = .TRUE.
                 EXIT
              END IF
           END IF
        END DO
        IF (.NOT.found) THEN
           PRINT *, "ERROR can't find time for: ", filename
           CALL quit(1)
        END IF
        READ(line,*) name,ktime,k
c        print *, ctime
c        print *, "converting time to real"
c        READ(ctime,"(E16.8)") time
        
        print *, ":) krome time=", ktime, "k=", k
        CLOSE(iread)
        ktime = ktime * spy

      END SUBROUTINE readtime

      INTEGER FUNCTION countmols()
      implicit none
      character*20 molfile,rmolfile
      integer dumpmols,unit,i,iostat
      logical haveRfile
c     NB rdumpmols.txt is list to read if different to write
      rmolfile = "rdumpmols.txt"
      unit = 336
      INQUIRE(FILE=rmolfile,EXIST=haveRfile)
      IF (haveRfile) THEN
         molfile = rmolfile
      ELSE
         molfile = "dumpmols.txt"
      END IF
      print *, "Counting molecules in", molfile
      OPEN (UNIT=unit,file=molfile,FORM='formatted',IOSTAT=iostat)
      i = 0
      DO
         READ(unit,*,IOSTAT=iostat)
         IF (iostat > 0) THEN
            print *, "Error reading mollist"
            STOP
         ELSE IF (iostat < 0) THEN
            EXIT
         ELSE
            i = i + 1
         ENDIF
      END DO
      countmols = i
      print *, "Found", countmols, "molecules to read."
      RETURN
      END FUNCTION countmols

c --------------

      SUBROUTINE readmolfile(maxmols,mollist,dumpmols)
      implicit none
      integer maxmols
      character*16 mollist(maxmols),tempmol
      character*20 molfile
      integer dumpmols,unit,i,iostat

      molfile = "dumpmols.txt"
      unit= 335
      print *, "Reading molecule file"
      OPEN (UNIT=unit,file=molfile,FORM='formatted',IOSTAT=iostat)
      DO i=1,maxmols
         READ(unit,*,IOSTAT=iostat) tempmol
         IF (iostat > 0) THEN
            print *, "Error reading mollist"
            STOP
         ELSE IF (iostat < 0) THEN
            EXIT
         ELSE
            mollist(i) = tempmol
         ENDIF
      END DO
      dumpmols = i - 1
      print *, "Found", dumpmols, "molecules to dump."
      CLOSE(unit)
      END SUBROUTINE readmolfile

c ***********************************************
c Used in profile.F
      SUBROUTINE get_mol_inds(mollist,maxmols,dumpmols,molinds)
      INCLUDE 'COMMONS/krome_mods'
      INCLUDE 'idim'
      INCLUDE 'COMMONS/kromevar'
      INTEGER maxmols,dumpmols,molinds(maxmols),j
      CHARACTER*16 mollist(maxmols)
      
      DO j=1,dumpmols
         molinds(j) = krome_get_index(adjustl(mollist(j)))
      END DO
      
      END SUBROUTINE get_mol_inds
