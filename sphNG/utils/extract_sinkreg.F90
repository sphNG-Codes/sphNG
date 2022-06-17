      PROGRAM extract_sinkreg
!******************************************************************
!  Extract particles within a region around a sink particle       *
!         AKY 04/2022
!******************************************************************

        implicit none
        INCLUDE 'omp_lib.h'
        INCLUDE 'idim'
        INCLUDE 'igrape'
        INCLUDE 'COMMONS/actio'
        INCLUDE 'COMMONS/binary'
        INCLUDE 'COMMONS/bodys'
        INCLUDE 'COMMONS/densi'
        INCLUDE 'COMMONS/gtime'
        INCLUDE 'COMMONS/part'
        INCLUDE 'COMMONS/phase'
        INCLUDE 'COMMONS/ptmass'
        INCLUDE 'COMMONS/raddust'
        INCLUDE 'COMMONS/radtrans'
        INCLUDE 'COMMONS/rbnd'
        INCLUDE 'COMMONS/recor'
        INCLUDE 'COMMONS/savernd'
        INCLUDE 'COMMONS/stepopt'
        INCLUDE 'COMMONS/timei'
        INCLUDE 'COMMONS/tming'
        INCLUDE 'COMMONS/varmhd'
        INCLUDE 'COMMONS/units'
        INCLUDE 'COMMONS/sort'
        INCLUDE 'COMMONS/interstellar'

        INTEGER,PARAMETER :: maxwant = 800000
        REAL,PARAMETER :: spy=31536000d0 
        INTEGER :: iindump,ichkl,i,nbin,lower,ibin
        INTEGER :: sink_ID,nselect,iout,dead,sink,other
        INTEGER :: gas,inp,ind1,ind2
        CHARACTER(len=15) :: charsinkno
        CHARACTER(len=25) :: infile
        CHARACTER(len=45) :: outfile,idlistfile
        INTEGER,DIMENSION(MAXWANT) :: wantedparts
        LOGICAL :: fromsink
        REAL*8 ::dtmaxdp
        COMMON /dtmaxin/ dtmaxdp
        iindump = 1
        iout  = 2
        gas = 0;dead=0;sink=0;other=0

        call unit
        
        PRINT *, "Enter name of input dump file"
        READ  (*,*) infile
        PRINT *, "Find particles 1) near sink or 2) from iunique?"
        READ (*,*) inp
        IF (inp .EQ. 1) THEN
           fromsink = .TRUE.
           PRINT *, "Enter iunique of desired sink"
           READ  (*,*) sink_ID
        ELSE IF (inp .EQ. 2) THEN
           fromsink = .FALSE.
           PRINT *, "Enter name of iunique list file"
           READ (*,*) idlistfile
        ELSE
           PRINT *, "Error invalid option"
           CALL quit(1)
        END IF
        
! ---- Read dumpfile               
        WRITE(*,*) "Reading", infile
        OPEN (UNIT=iindump,FILE=infile,FORM ='unformatted',
     & ACTION='READ')
        CALL rdump(iindump, ichkl, 0)
        CLOSE(iindump)
        dtmax = dtmaxdp
        nfullstep = 1
        WRITE(*,284) "sphNG time= ",
     &       gt*utime/spy," years. Read", npart,"particles."
        WRITE(*,*) "(utime=",utime," gt=", gt,"udist=",udist,
     &       "dtmax=", dtmax
284     FORMAT(A,X,E9.2,A,I0,A)

        WRITE (*,285) n1,n2,nreassign,naccrete,nkill,iyr,idum
285     FORMAT("n1=",I0,/,"n2=",I0,/,"nreassign=",I0,/,"naccrete=",I0,
     &        /,"nkill=",I0,/, "iyr=",I0,/,"idum=",I0)
        IF (fromsink) THEN 
           CALL selectfromregion(sink_ID,wantedparts,maxwant,nselect)
           CALL extract_wanted(wantedparts,maxwant,nselect)
           WRITE (*,*) "Extracted particle data"
!           WRITE (*,*) xyzmh(4,wantedparts(12)),rho(wantedparts(12)),
!     &    iphase(wantedparts(12))
           WRITE (charsinkno,'(I8)') sink_ID
           PRINT *, charsinkno, sink_ID
           charsinkno = adjustl(charsinkno)
           outfile = trim(infile) // "_sink" // charsinkno
        END IF

        IF (.NOT. fromsink) THEN
           CALL get_particles_fromID(idlistfile,maxwant,nselect)
           ind1 = INDEX(idlistfile,"_")
           ind2 = INDEX(idlistfile,".dat")
           outfile = trim(infile) // idlistfile(ind1:ind2-1) //"parts"
        END IF
        
! ----- Write to new dump file           
        PRINT *, "writing to:", outfile
        print *, "npart=", npart
        npart = nselect
        print *, "npart=", npart

        DO i=1, npart
         IF (iphase(i) .EQ. 0) THEN
            gas = gas + 1
         ELSEIF (iphase(i) .EQ. -1) THEN
            dead = dead + 1
         ELSE IF (iphase(i) .EQ. 1) THEN
            sink = sink + 1
            print *, "type 1"
            iphase(i) = -1
         ELSE IF (iphase(i) .EQ. 2 ) THEN
            sink = sink+1
            print *, "type 2"
         ELSE IF (iphase(i) .GT. 2 .AND. iphase(i) .LE.6) THEN
            sink = sink + 1
            print *, "other sink"
         ELSE
            other = other+1
            print *, "weird type", iphase(i)
         END IF
       END DO
       PRINT *, "Dead =", dead, " sink=", sink, "nptmass=", nptmass, 
     &  "other =", other, "gas=", gas
       
!----- test listinactive
       DO i=1,nlistinactive
          IF (iphase(listinactive(i)) .NE. -1) THEN
             PRINT *, "Error:iphase of",i,"is",iphase(listinactive(i))
          END IF
       END DO

!----- recentre coords and velocities on heaviest sink
       CALL recentre()
       PRINT *, "resetting time to zero"
       gt = 0.              
       DO i=1,npart
          isort(i) = i
          iorig(i) = i
       END DO
       PRINT *, "checking listpm..."
       PRINT *, (listpm(i),i=1,nptmass)
       PRINT *, "checking density array"
       DO i=1,npart
          IF (rho(i) .EQ. 0.0) THEN
             WRITE (*,*) "rho=0", i,iunique(i),iphase(i),rho(i)
          END IF
       END DO
       OPEN(UNIT=iout,FILE=trim(outfile),FORM='unformatted',
     &        STATUS='replace')
        
       CALL wdump(iout)
       CLOSE(iout)
       IF (fromsink) CALL write_IDfile(outfile)
       WRITE (*,*) "done"
      END PROGRAM extract_sinkreg
      
!------------------------------------------------------

       SUBROUTINE quit(i)
         STOP
       END SUBROUTINE quit

!------------------------------------------------------
       
       SUBROUTINE selectfromregion(sink_ID,wantedparts,maxwant,nselect)
         implicit none
         INCLUDE 'idim'
         INCLUDE 'COMMONS/astrcon'
         INCLUDE 'COMMONS/part'
         INCLUDE 'COMMONS/phase'
         INCLUDE 'COMMONS/ptmass'
         INCLUDE 'COMMONS/densi'
         INCLUDE 'COMMONS/sort'
         INCLUDE 'COMMONS/units'
         INTEGER ::  sink_ID,maxwant
         REAL :: rwant,sinkx,sinky,sinkz,rwant2,r2
         INTEGER :: i,nselect,ipart
         INTEGER,dimension(maxwant) :: wantedparts
         LOGICAL :: stoploop
         
         WRITE(*,*) "Enter radius from sink in AU"
         READ (*,*) rwant
         WRITE(*,*) "Finding particles r <",rwant,"au of sink",sink_ID
         rwant = rwant * au / udist
         rwant2 = rwant * rwant
         WRITE(*,*) "==", rwant,"code units."
! ----- Find index of sink
         stoploop = .FALSE.
         DO i=1,npart
            IF (iunique(i) .EQ. sink_ID) THEN
               ipart = i
               stoploop = .TRUE.
               EXIT
            END IF
         END DO

         IF (.NOT. stoploop) THEN
            WRITE (*,*) "Sink ID not found in list", sink_ID
            CALL quit(0)
         ELSE IF (iphase(ipart) .NE. 2) THEN
            WRITE (*,*) sink_ID, "is not a sink. iphase=",iphase(ipart)
            CALL quit(0)
         END IF
         sinkx = xyzmh(1,ipart)
         sinky = xyzmh(2,ipart)
         sinkz = xyzmh(3,ipart)
!        loop over all particles to find particles within r<rwant

         nselect = 0
         DO i=1,npart
            r2 = (xyzmh(1,i)-sinkx)**2 + (xyzmh(2,i)-sinky)**2
     &      + (xyzmh(3,i)-sinkz)**2
     
            IF (r2 .LT. rwant2) THEN
!               print *, rwant2, r2
               nselect = nselect + 1
               wantedparts(nselect) = i
               IF (iunique(i) .eq. iunique(ipart)) then
                 write(*,*) "found the sink",iunique(i),iphase(i)
               END IF
            END IF
         END DO
         WRITE (*,*) "Selected", nselect, "of",npart, "particles."
         
       END SUBROUTINE selectfromregion

       SUBROUTINE get_particles_fromID(idlistfile,maxwant,nselect)
         IMPLICIT NONE
         INCLUDE 'idim' 
         INCLUDE 'omp_lib.h'
         INCLUDE 'COMMONS/part'
         INCLUDE 'COMMONS/phase'
         INCLUDE 'COMMONS/ptmass'
         INCLUDE 'COMMONS/densi'
         INCLUDE 'COMMONS/sort'
         INCLUDE 'COMMONS/units'
         CHARACTER(len=45) :: idlistfile,junk
         INTEGER ::  maxwant,iostatus
         INTEGER :: i,nselect,ipart,iunit=19,j
         INTEGER,dimension(maxwant) :: wantedparts
         INTEGER*8, dimension(maxwant) :: wantIDs
         REAL :: t_start, t_end

         iostatus = 0
         wantIDs(:) = 0
         wantedparts(:) = 0
         WRITE (*,*) "Reading ID file:", idlistfile
         OPEN(UNIT=iunit,FILE=trim(idlistfile),FORM='formatted',
     &        STATUS='old')
         DO i=1, maxwant
            READ (iunit,*,iostat=iostatus) wantIDs(i)
            IF (iostatus .GT. 0) THEN
               WRITE (*,*) "Error in read:", iostatus, wantIDs(i)
               CALL quit(1)
            ELSE IF (iostatus .LT. 0) THEN
               WRITE (*,*) "End of file read"
               EXIT
            END IF
         END DO
         CLOSE(iunit)
         PRINT *, "iunique", i-1, "is", wantIDs(i-1)
         nselect = i-1
         WRITE (*,*) "Found", nselect, "wanted particles"
! -----  Make list of indices of wanted particles
         ipart = 0
!$OMP PARALLEL
#ifdef _OPENMP
         WRITE(*,*) OMP_GET_NUM_THREADS(), "OMP threads"
#else
         WRITE(*,*) "Running with single thread"
#endif
!$OMP END PARALLEL
         CALL CPU_TIME(t_start)
!$OMP PARALLEL DO PRIVATE(j) SHARED(npart,nselect,wantIDs,
!$OMP&   iunique,wantedparts,ipart,rho) SCHEDULE(static,100)
         DO i=1, npart
            IF (iunique(i) .EQ. 0 ) CYCLE
            DO j=1, nselect
               IF (wantIDs(j) .EQ. 0) CYCLE
               IF (iunique(i) .EQ. wantIDs(j) ) THEN
                  IF (iphase(i) .GE. 0   .AND.
     &            rho(i) .GT. 0.0) THEN
!                 get rid of dead particles                  
!$OMP ATOMIC UPDATE
                     ipart = ipart + 1  
                     wantedparts(ipart) = i
                  END IF
! wipe wantIDs(j)
                  wantIDs(j) = 0
                  EXIT
               END IF
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL CPU_TIME(t_end)
         PRINT *, "nfound =",ipart, "nselect=",nselect
         nselect = ipart
         PRINT *, "loop time = ", t_end - t_start
         WRITE (*,*) "Continue?"
         READ (*,*) junk
         CALL extract_wanted(wantedparts,maxwant,nselect)
         
       END SUBROUTINE get_particles_fromID

! -------------------------------------------------------------
       SUBROUTINE extract_wanted(wantedparts,maxwant,nselect)
         implicit none
         INCLUDE 'idim'
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
         
         INTEGER,intent(in) :: maxwant,nselect
         INTEGER,DIMENSION(maxwant),intent(in) :: wantedparts
         INTEGER :: j,nonzero=0,i
         
!       extract sinks first
         CALL extract_sinks(maxwant,wantedparts,nselect)
         DO i=1, maxwant
            if (wantedparts(i) .GT. 0.0) then
               nonzero=nonzero + 1
            end if
         END DO
         WRITE (*,*) nonzero, "nonzero values",nselect
         CALL extract_I8(iunique,maxwant,wantedparts,nselect)
         CALL extract_ID(isteps,maxwant,wantedparts,nselect)
         CALL extract_I1(iphase,maxwant,wantedparts,nselect)
         CALL rebuild_inactive(nselect)
         write (*,*) "listpm", listpm(1:2)
         write (*,*) "iphase", iphase(listpm(1)),iphase(listpm(2))

         DO j=1,5
            CALL extract_RD(xyzmh(j,:),maxwant,wantedparts,nselect)
         END DO
         PRINT *, "mass", xyzmh(4,listpm(1)), xyzmh(4,listpm(2))
         DO j=1,4
            CALL extract_RD(vxyzu(j,:),maxwant,wantedparts,nselect)
         END DO
         PRINT *, "DONE vxyzu"
         print *, "old rho", rho(wantedparts(1:10))
         CALL extract_R4(rho,maxwant,wantedparts,nselect)
         PRINT *, "new RHO", rho(1:10),rho(nselect-10), rho(nselect+10)
         DO j=1,2
            CALL extract_R4(gradhs(j,:),maxwant,wantedparts,nselect)
         END DO
         PRINT *, "DONE gradhs"
         PRINT *, "Wanted parts", wantedparts(1:10)
         CALL extract_R4(alphaMM(1,:),maxwant,wantedparts,nselect)
         DO j=1,5
            CALL extract_RD(ekcle(j,:),maxwant,wantedparts,nselect)
         END DO
         PRINT *, "DONE ekcle"
         DO j=1,2
            CALL extract_RD(dust_tk(j,:),maxwant,wantedparts,nselect)
            PRINT *, "DONE dust"
         END DO
         IF (idustRT.GT.0) THEN
            DO j=1, nchemistry
               CALL extract_R4(chemistry(j,:),maxwant,wantedparts,
     &                   nselect)
            END DO
            DO j=1, nheatingISR
               CALL extract_R4(heatingISR(j,:),maxwant,wantedparts
     &                  ,nselect)
               PRINT *, "Done heating ISR"
            END DO
            CALL extract_R4(h2frac,maxwant,wantedparts,nselect)
            CALL extract_R4(dh2dt,maxwant,wantedparts,nselect)

         END IF

         n1 = nselect
         n2 = 0
         naccrete = nlistinactive
         print *, "Done extraction"

       END SUBROUTINE extract_wanted

       !----------------------------------
       SUBROUTINE extract_I1(array,maxwant,wanted,nselect)
         implicit none
         INCLUDE 'idim'
         integer*1,dimension(:) :: array(idim)
         integer,intent(IN) :: maxwant,nselect
         integer,intent(IN) :: wanted(maxwant)
         integer*1,dimension(:) :: tmp_arr(maxwant)
         integer :: i
         print *, "extracting I1 array"
         tmp_arr(:) = -2
         do i=1, nselect
            tmp_arr(i) = array(wanted(i))
         end do
         array(:) = -2
         do i=1, nselect
            array(i) = tmp_arr(i)
         end do
       end SUBROUTINE extract_I1

!----------------------------------
       SUBROUTINE extract_I4(array,maxwant,wanted,nselect)
         implicit none
         INCLUDE 'idim'
         integer*4,dimension(:) :: array(idim)
         integer,intent(IN) :: maxwant,nselect
         integer,intent(IN) :: wanted(maxwant)
         integer*4,dimension(:) :: tmp_arr(maxwant)
         integer :: i
         print *, "extracting I4 array"
         tmp_arr(:) = 0
         do i=1, nselect
            tmp_arr(i) = array(wanted(i))
         end do
         array(:) = 0
         do i=1, nselect
            array(i) = tmp_arr(i)
         end do
       end SUBROUTINE extract_I4

!----------------------------------------------------------
       SUBROUTINE extract_I8(array,maxwant,wantedparts,nselect)
         implicit none
         INCLUDE 'idim'
         integer*8,dimension(:) :: array(idim)
         integer,intent(IN) :: maxwant,nselect
         integer,intent(IN) :: wantedparts(maxwant)
         integer*8,dimension(:) :: tmp_arr(maxwant)
         integer :: i
         print *, "extracting I8 array"
         tmp_arr(:) = 0
         do i=1, nselect
            tmp_arr(i) = array(wantedparts(i))
         end do
         array(:) = 0
         do i=1, nselect
            array(i) = tmp_arr(i)
         end do
       end SUBROUTINE extract_I8

!----------------------------------------------------------
       SUBROUTINE extract_ID(array,maxwant,wantedparts,nselect)
         implicit none
         INCLUDE 'idim'
         integer,dimension(:) :: array(idim)
         integer,intent(IN) :: maxwant,nselect
         integer,intent(IN) :: wantedparts(maxwant)
         integer,dimension(:) :: tmp_arr(maxwant)
         integer :: i
         print *, "extracting Default int array"
         tmp_arr(:) = 0
         do i=1, nselect
            tmp_arr(i) = array(wantedparts(i))
         end do
         array(:) = 0
         do i=1, nselect
            array(i) = tmp_arr(i)
         end do
       end SUBROUTINE extract_ID

         
!       extract batch of default REALs
      subroutine extract_RD(array,nmax,tmpwanted,nselect)
        implicit none
        INCLUDE 'idim'
        real,dimension(:) :: array(idim)
        integer,intent(IN) :: nmax,nselect
        integer,intent(IN) :: tmpwanted(nmax)
        real,dimension(:) :: tmp_arr(nmax)
        integer :: i
        print *, "extracting Default Real array"

        tmp_arr(:) = 0
        do i=1, nselect
           tmp_arr(i) = array(tmpwanted(i))
        end do
        array(:) = 0
        do i=1, nselect
           array(i) = tmp_arr(i)
        end do
!       print *, "array reset"
!        print *, "RD:", tmpwanted(1:10)
      end subroutine extract_RD


!-----------------------------------------------------------------

      !       extract batch of default REALs
      subroutine extract_R4(array,maxwant,wantedparts,nselect)
        implicit none
        INCLUDE 'idim'
        real*4,dimension(:) :: array(idim)
        integer,intent(IN) :: maxwant,nselect
        integer,intent(IN) :: wantedparts(maxwant)
        real*4,dimension(:) :: tmp_arr(maxwant)
        integer :: i
        print *, "extracting R4 array"

        tmp_arr(:) = 0
        do i=1, nselect
           tmp_arr(i) = array(wantedparts(i))
        end do
        array(:) = 0
        do i=1, nselect
           array(i) = tmp_arr(i)
        end do
!       print *, "array reset"
      end subroutine extract_R4

!       extract batch of default REALs
      subroutine extract_RDsink(array,nmax,inewlist,nselect,iptwant)
        implicit none
        INCLUDE 'idim'
        real,dimension(:) :: array(iptdim)
        integer,intent(IN) :: nmax,nselect
        integer,intent(IN) :: inewlist(nmax),iptwant(nmax)
        real,dimension(:) :: tmp_arr(nmax)
        integer :: i
        print *, "extracting Default Real Sink array"

        tmp_arr(:) = 0
        do i=1, nselect
           tmp_arr(i) = array(iptwant(i))
           print *, tmp_arr(i),iptwant(i),i
        end do
        array(:) = 0
        do i=1, nselect
           array(i) = tmp_arr(inewlist(i))
           print *, array(i),inewlist(i),i
        end do
      end subroutine extract_RDsink

!-----------------------------------------------------------------
! handles all the sink variables
      subroutine extract_sinks(nmax,wanted,nselect)
        implicit none
        INCLUDE 'idim'
        INCLUDE 'COMMONS/part'
        INCLUDE 'COMMONS/phase'
        INCLUDE 'COMMONS/ptmass'
        INCLUDE 'COMMONS/sort'
        integer,intent(IN) :: nmax
        integer,dimension(nmax),intent(IN) :: wanted
        INTEGER,INTENT(IN) :: nselect
        integer,dimension(:) ::inewlist(nptmass),tmp(nptmass)
        integer :: i,j,npmwant,imassive,iwantmassive,n
        INTEGER :: iwant(nptmass),iptwant(nptmass)
        REAL :: mtmp

! tmp(nptmass)= index of sink in wanted old arrays
! iwant(nptmass)= index of sink in wanted list
! inewlist(nptmass)= index of sink in new ptmass arrays
! iptwant(nptmass)= index of sink in old ptmass arrays
        tmp(:) = 0
        inewlist(:) = 0
        npmwant = 0

        iwant(:) = 0
        iptwant(:) = 0
! ----- Loop over wantedparts
        DO i=1,nselect
           IF (iphase(wanted(i)) .EQ. 2 ) then
              npmwant = npmwant + 1
              tmp(npmwant) = wanted(i)
              iwant(npmwant) = i
! find corresponding index in old ptmass arrays
              DO j=1,nptmass
                 IF (iunique(wanted(i)) .EQ. iunique(listpm(j))) THEN
                    iptwant(npmwant) = j
                    EXIT
                 END IF
              END DO
           END IF
        END DO
        IF (npmwant .LT. 1) THEN
           WRITE(*,*) "No sinks!"
           nptmass = 0
           return
        ELSE
           PRINT *, "listpms wanted:", (iptwant(j),j=1,npmwant)
        END IF

! ----- Find most massive sink and put first in list
! ----- so it can be centered.
        mtmp = 0.
        DO i=1, npmwant
           IF (xyzmh(4,tmp(i)).GT. mtmp) THEN
              mtmp = xyzmh(4,tmp(i))
              imassive = i
              iwantmassive = tmp(i)
           END IF
        END DO
        listpm(:) = 0
! making new listpm
        listpm(1) = iwant(imassive)
        n=1
        DO i=1, npmwant
           IF (i .EQ. imassive) THEN
              inewlist(1) = i
              CYCLE              
           END IF
           n = n + 1
           inewlist(n) = i
           listpm(n) = iwant(i)
        END DO
        IF (n .NE. npmwant) THEN
           WRITE(*,*) "Error: something wrong with ptm lists"
           CALL quit(0)
        END IF
        WRITE (*,*) "Sink masses:"
        WRITE (*,*) (tmp(i),xyzmh(4,tmp(i)),i=1,npmwant)
        WRITE (*,*) "biggest is ", wanted(imassive),mtmp,
     & xyzmh(4,wanted(imassive))
        
        nptmass = npmwant
        WRITE (*,*) "nptmass want", npmwant, nptmass
        print *, "listpm=",listpm(1:nptmass)
        print *, "inewlist=", inewlist(1:nptmass)
        call extract_RDsink(spinx,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(spiny,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(spinz,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(angaddx,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(angaddy,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(angaddz,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(spinadx,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(spinady,iptdim,inewlist,npmwant,iptwant)
        call extract_RDsink(spinadz,iptdim,inewlist,npmwant,iptwant)
        DO i=1,npmwant
           print *, "spins", spinx(i),spiny(i),spinz(i)
        END DO
        print *, "extracted sink data"
      end subroutine extract_sinks

! -------------------------------------------------
      SUBROUTINE write_IDfile(outfilename)
        INCLUDE 'idim'
        INCLUDE 'COMMONS/part'
        INCLUDE 'COMMONS/sort'
        CHARACTER(len=45),INTENT(IN) :: outfilename
        INTEGER :: i,iunit=11
        CHARACTER(len=45) :: filename

        filename = trim(outfilename) // "_IDs.dat"
        OPEN(UNIT=iunit,FILE=trim(filename),FORM='formatted',
     &        STATUS='replace')
        DO i=i, npart
           WRITE(iunit,'(I8)') iunique(i)
        END DO
        CLOSE(iunit)
        WRITE (*,*) "written iuniques to file:", filename
      END SUBROUTINE

! ------------
      
      SUBROUTINE rebuild_inactive(nselect)
        implicit none
        INCLUDE 'idim'
        INCLUDE 'COMMONS/rbnd'
        INCLUDE 'COMMONS/phase'

        INTEGER,INTENT(IN) :: nselect
        INTEGER :: i
        
        WRITE (*,*) "Rebuilding listinactive"
        
        listinactive(:) = 0
        nlistinactive = 0
        DO i=1,nselect
           IF (iphase(i) .LT. 0) THEN
              nlistinactive = nlistinactive + 1
              listinactive(nlistinactive) = i
           END IF
        END DO
        WRITE (*,*) "(extract) New nlistinactive=", nlistinactive
           
      END SUBROUTINE rebuild_inactive

! -----
! ----- recentre coords and vels relative to heaviest sink
! ----- (this is the first in listpm)
       SUBROUTINE recentre()
        IMPLICIT NONE
        INCLUDE 'idim'
        INCLUDE 'omp_lib.h'
        INCLUDE 'COMMONS/part'
        INCLUDE 'COMMONS/ptmass'

        INTEGER :: i,j
        REAL,DIMENSION(3) :: sink_xyz,sink_vxyz
        sink_xyz(:) = xyzmh(1:3,listpm(1))
        sink_vxyz(:) = vxyzu(1:3, listpm(1))
        PRINT *, "Sink loc:", sink_xyz
        PRINT *, sink_vxyz
!$OMP PARALLEL DO SCHEDULE(static,100) SHARED(npart,xyzmh,
!$OMP& listpm,vxyzu,sink_xyz,sink_vxyz) PRIVATE(i,j)
        DO i=1,npart
           DO j=1,3
              xyzmh(j,i) = xyzmh(j,i) - sink_xyz(j) 
              vxyzu(j,i) = vxyzu(j,i) - sink_vxyz(j)
           END DO
        END DO
!$OMP END PARALLEL DO

        PRINT *, "Sink loc 2:", xyzmh(:,listpm(1))
        PRINT *, vxyzu(:,listpm(1))
        WRITE(*,*) "Recentered particles"

      END SUBROUTINE recentre
