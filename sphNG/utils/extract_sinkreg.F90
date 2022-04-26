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
        INCLUDE 'COMMONS/densi'
        INCLUDE 'COMMONS/gtime'
        INCLUDE 'COMMONS/part'
        INCLUDE 'COMMONS/phase'
        INCLUDE 'COMMONS/ptmass'
        INCLUDE 'COMMONS/raddust'
        INCLUDE 'COMMONS/radtrans'
        INCLUDE 'COMMONS/recor'
        INCLUDE 'COMMONS/timei'
        INCLUDE 'COMMONS/tming'
        INCLUDE 'COMMONS/varmhd'
        INCLUDE 'COMMONS/units'
        INCLUDE 'COMMONS/sort'
        INCLUDE 'COMMONS/interstellar'

        INTEGER,PARAMETER :: maxwant = 800000
        REAL,PARAMETER :: spy=31536000d0 
        INTEGER :: iindump,ichkl,i,nbin,lower,ibin,sinkno
        INTEGER :: isink,ipart,nselect,iout,dead,sink,other
        INTEGER :: gas,inp
        CHARACTER(len=15) :: charsinkno
        CHARACTER(len=25) :: infile
        CHARACTER(len=45) :: outfile
        INTEGER,DIMENSION(MAXWANT) :: wantedparts
        iindump = 1
        iout  = 2
        gas = 0;dead=0;sink=0;other=0

        call unit
        
        PRINT *, "Enter name of input dump file"
        READ  (*,*) infile
        PRINT *, "Find particles 1) near sink or 2) from iunique?"
        PRINT *, "Enter number of desired sink"
        READ  (*,*) sinkno
        
        WRITE(*,*) "Reading", infile
! ---- Read dumpfile       
      OPEN (UNIT=iindump,FILE=infile,FORM ='unformatted',
     & ACTION='READ')
      CALL rdump(iindump, ichkl, 0)
      CLOSE(iindump)
      nfullstep = 1
      WRITE(*,284) "sphNG time= ",
     &       gt*utime/spy," years. Read", npart,"particles."
        WRITE(*,*) "(utime=",utime," gt=", gt,"udist=",udist
 284    FORMAT(A,X,E9.2,A,I0,A)

        ipart = listpm(sinkno)
        WRITE (*,*) "looking for particle", ipart,iphase(ipart)
! ----- Check sink exists
        IF (listpm(sinkno) .LT. 1) THEN
           WRITE (*,*) "Sink particle doesn't exist"
           CALL quit(1)
        END IF
        CALL selectfromregion(ipart,wantedparts,maxwant,nselect)
        WRITE (*,*) wantedparts(300)
        CALL extract_wanted(wantedparts,maxwant,nselect)
        WRITE (*,*) "Extracted particle data"

        WRITE (*,*) xyzmh(4,wantedparts(12)),rho(wantedparts(12)),
     &    iphase(wantedparts(12))
        WRITE (charsinkno,'(I3)') sinkno
        PRINT *, charsinkno, sinkno
        charsinkno = adjustl(charsinkno)
        outfile = trim(infile) // "_sink" // charsinkno
        PRINT *, "writing to:", outfile
        print *, "npart=", npart
        npart = nselect
        print *, "npart=", npart

        DO i=1, npart
         IF (iphase(i) .EQ. -1) THEN
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
         ELSE IF (iphase(i) .NE. 0) THEN
            other = other+1
            print *, "weird type", iphase(i)
         ELSE IF (iphase(i) .EQ. 0) THEN
            gas = gas + 1
         END IF
       END DO
       PRINT *, "Dead =", dead, " sink=", sink, "nptmass=", nptmass, 
     &  "other =", other, "gas=", gas
       
       CALL write_IDfile(outfile)
       DO i=1,npart
          isort(i) = i
       END DO
       OPEN(UNIT=iout,FILE=trim(outfile),FORM='unformatted',
     &        STATUS='replace')
        
       CALL wdump(iout)
       CLOSE(iout)

       WRITE (*,*) "done"
      END PROGRAM extract_sinkreg
      
!------------------------------------------------------

       SUBROUTINE quit(i)
         STOP
       END SUBROUTINE quit

!------------------------------------------------------
       
       SUBROUTINE selectfromregion(ipart,wantedparts,maxwant,nselect)
         implicit none
         INCLUDE 'idim'
         INCLUDE 'COMMONS/astrcon'
         INCLUDE 'COMMONS/part'
         INCLUDE 'COMMONS/phase'
         INCLUDE 'COMMONS/ptmass'
         INCLUDE 'COMMONS/densi'
         INCLUDE 'COMMONS/sort'
         INCLUDE 'COMMONS/units'
         INTEGER ::  maxwant
         REAL :: rwant,sinkx,sinky,sinkz,rwant2,r2
         INTEGER :: i,nselect,ipart
         INTEGER,dimension(maxwant) :: wantedparts
         
         WRITE(*,*) "Enter radius from sink in AU"
         READ (*,*) rwant
         WRITE(*,*) "Finding particles within", rwant, "au of sink"
         rwant = rwant * au / udist
         rwant2 = rwant * rwant
         WRITE(*,*) "==", rwant,"code units."
         sinkx = xyzmh(1,ipart)
         sinky = xyzmh(2,ipart)
         sinkz = xyzmh(3,ipart)
!        loop over all particles to find particles within r<rwant
!        OMP this!

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
         
!       extract sinks first!!
         PRINT *, "Wanted parts1", wantedparts(1:10)
         CALL extract_sinks(maxwant,wantedparts,nselect)
         PRINT *, "Wanted parts, after sinks:", wantedparts(1:10) 
         DO i=1, maxwant
            if (wantedparts(i) .GT. 0.0) then
               nonzero=nonzero + 1
            end if
         END DO
         WRITE (*,*) nonzero, "nonzero values"
         CALL extract_I8(iunique,maxwant,wantedparts,nselect)
         CALL extract_ID(isteps,maxwant,wantedparts,nselect)
         CALL extract_I1(iphase,maxwant,wantedparts,nselect)
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
      subroutine extract_RDsink(array,nmax,tmpwanted,nselect)
        implicit none
        INCLUDE 'idim'
        real,dimension(:) :: array(iptdim)
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
      end subroutine extract_RDsink

!-----------------------------------------------------------------
! handles all the sink variables
      subroutine extract_sinks(nmax,wanted,nselect)
        implicit none
        INCLUDE 'idim'
        INCLUDE 'COMMONS/phase'
        INCLUDE 'COMMONS/ptmass'
        integer,intent(IN) :: nmax
        integer,dimension(nmax),intent(IN) :: wanted
        INTEGER,INTENT(IN) :: nselect
        integer,dimension(:) :: tmplist(nptmass)
        integer :: i,j,npmwant

        tmplist(:) = 0
        npmwant = 0
        listpm(:) = 0
        !       Loop over wantedparts
        DO i=1,nselect
           IF (iphase(wanted(i)) .EQ. 2 ) then
              npmwant = npmwant + 1
              tmplist(npmwant) = wanted(i)
              listpm(npmwant) = i
           END IF
        END DO
        IF (npmwant .LT. 1) THEN
           WRITE(*,*) "Error, no sinks!"
           CALL quit(1)
        END IF
        nptmass = npmwant
        WRITE (*,*) "nptmass want", npmwant, nptmass
        print *, "listpm=",listpm(1:nptmass)
        call extract_RDsink(spinx,iptdim,tmplist,npmwant)
        call extract_RDsink(spiny,iptdim,tmplist,npmwant)
        call extract_RDsink(spinz,iptdim,tmplist,npmwant)
        call extract_RDsink(angaddx,iptdim,tmplist,npmwant)
        call extract_RDsink(angaddy,iptdim,tmplist,npmwant)
        call extract_RDsink(angaddz,iptdim,tmplist,npmwant)
        call extract_RDsink(spinadx,iptdim,tmplist,npmwant)
        call extract_RDsink(spinady,iptdim,tmplist,npmwant)
        call extract_RDsink(spinadz,iptdim,tmplist,npmwant)
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

        
        
