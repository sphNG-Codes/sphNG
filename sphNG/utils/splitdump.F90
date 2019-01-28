      PROGRAM splitdump

! ****************************************************
!  Splits an MPI dump into n smaller dumps by        *
!  iuniques. ./splitdump ifile                                          *
!  AKY Jan 2019                                      *
! ****************************************************
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
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/densi'
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
      INCLUDE 'COMMONS/divcurlB'
      INCLUDE 'COMMONS/Bxyz'
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


      INTEGER :: iindump,ichkl,i,nbin,lower,ibin
      INTEGER :: binsize,average,j,sinkcheck
      INTEGER,ALLOCATABLE :: bins(:,:),bintots(:),iout(:),countsink(:)
      INTEGER,ALLOCATABLE :: tmp_listpm(:,:)
      CHARACTER(len=7) :: infile
      CHARACTER(len=9) :: outfile
      iindump = 1
      iterm  = 84
  
      PRINT *, "Enter name of file to split"
      READ  (*,'(A7)') infile
      PRINT *, "How many files to split it into?"
      READ (*,*) nbin
      IF (nbin .GT. 9) THEN
         PRINT *, "ERROR. No. files cannot exceed 9"
         CALL QUIT(1)
      END IF
      PRINT *, "Enter name of ifile"
      READ (*,*) inname
      CALL options
! ---- Read dumpfile       
      OPEN (UNIT=iindump,FILE=infile,FORM ='unformatted',
     &            ACTION='READ')
!     RECL=maxrec,         
      CALL rdump(iindump, ichkl, 0)
      CLOSE(iindump)

      WRITE(*,284) infile,"Read ", npart,"particles."
 284        FORMAT(A7,X,A,I0,X,A)
      
!     Allocate bin arrays
      binsize = (iuniquemax/nbin) + MOD(iuniquemax,nbin)
!      binsize = get_binsize(iuniquemax,nbin)
      allocate(bins(binsize,nbin))
      allocate(tmp_listpm(iptdim,nbin))
      allocate(bintots(nbin))
      allocate(iout(nbin))
      allocate(countsink(nbin))
      DO ibin=1, nbin
         bintots(ibin) = 0
         countsink(ibin) = 0
      END DO
      DO j=1, nbin
         iout(j) = j + 3
      END DO

!     Loop over particles
      DO i=1, npart
         lower = 0
! loop over bins
        DO ibin=1, nbin 
           IF (iphase(i).LT.0) CYCLE
           IF ((iunique(i).GT.lower) .AND.
     &          (iunique(i).LE.(ibin*binsize))) THEN
            bintots(ibin) = bintots(ibin) + 1
            bins(bintots(ibin),ibin) = i
            IF (iphase(i).GT.1) THEN
               countsink(ibin) = countsink(ibin) + 1
               PRINT *, "SINK: iphase=",iphase(i)
               tmp_listpm(countsink(ibin),ibin) = i
            END IF
            EXIT
           END IF
           lower = ibin*binsize
        END DO
      END DO
      PRINT*, "Bin totals"
      DO ibin=1, nbin
         print *, bintots(ibin)
      END DO
      PRINT *, "Bin 1:", binsize
      average = 0
      sinkcheck = 0
      DO ibin=1, nbin
         average = average + bintots(ibin) 
         sinkcheck = sinkcheck + countsink(ibin)
      END DO
  
      PRINT *, "iuniquemax: ", iuniquemax
      PRINT *, "average bin occupation=", average/nbin
      PRINT *, "bin total particles= ", average
      IF (sinkcheck .NE. nptmass) THEN
         PRINT *, "ERROR: sinkcheck NE nptmass"
         STOP
      END IF

      DO ibin=1, nbin
         IF (ibin .GT. 1) THEN
            OPEN (UNIT=iindump,FILE=infile,FORM ='unformatted',
     &            ACTION='READ')
            CALL rdump(iindump, ichkl, 0)
            CLOSE(iindump)
         END IF
         PRINT *, "idustRT", idustRT, nchemistry, nheatingISR
         PRINT *, "rdumptests", npart, rho(ibin*binsize),
     &          ekcle(3,ibin*binsize)
         PRINT *, "nums 1", naccretetot, nblocks,xmax,xmin
         npart = bintots(ibin)
         PRINT *, "ibin= ", ibin,"npart= ", npart
         CALL extract_I8(iunique,binsize,bins(:,ibin),bintots(ibin))
         CALL extract_ID(isteps,binsize,bins(:,ibin),bintots(ibin))
         CALL extract_I1(iphase,binsize,bins(:,ibin),bintots(ibin))
         DO j=1,5
            CALL extract_RD(xyzmh(j,:),binsize,bins(:,ibin),
     &           bintots(ibin))
         END DO
         DO j=1,4
            CALL extract_RD(vxyzu(j,:),binsize,bins(:,ibin),
     &           bintots(ibin))
         END DO
         CALL extract_R4(rho,binsize,bins(:,ibin),bintots(ibin))
         DO j=1,2
            CALL extract_R4(gradhs(j,:),binsize,bins(:,ibin),
     &           bintots(ibin))
         END DO
         CALL extract_R4(alphaMM(1,:),binsize,bins(:,ibin),
     &          bintots(ibin))
         DO j=1,5
            CALL extract_RD(ekcle(j,:),binsize,bins(:,ibin),
     &          bintots(ibin))
         END DO
         DO j=1,2
            CALL extract_RD(dust_tk(j,:),binsize,bins(:,ibin),
     &          bintots(ibin))
         END DO
         IF (idustRT.GT.0) THEN
            DO j=1, nchemistry
               CALL extract_R4(chemistry(j,:),binsize,
     &              bins(:,ibin),bintots(ibin))
            END DO
            DO j=1, nheatingISR
               CALL extract_R4(heatingISR(j,:),binsize,bins(:,ibin),
     &              bintots(ibin))
            END DO
            CALL extract_R4(h2frac,binsize,bins(:,ibin),
     &              bintots(ibin))
            CALL extract_R4(dh2dt,binsize,bins(:,ibin),
     &              bintots(ibin))

         END IF
         IF (countsink(ibin) .GT. 0) THEN
            CALL extract_sinks(countsink(ibin),tmp_listpm(:,ibin))
         END IF
         IF (ibin .LT. nbin) THEN
            iuniquemax = ibin * binsize
         END IF
         PRINT *, "TESTS", rho(binsize),ekcle(3,binsize)
         PRINT *, "nums 2", naccretetot, nblocks,xmax,xmin
         WRITE(outfile,12) infile, ibin
12       FORMAT(A7,"_",I1)
!         OPEN(UNIT=iout(ibin),FILE=outfile,FORM="formatted", 
!     &        STATUS="replace", ACTION="write",IOSTAT=iostat)
 !        WRITE(iout(ibin),'(I6)') (iunique(i), i=1,npart)
         DO i=1,npart
            isort(i) = i
         END DO
         OPEN(UNIT=iout(ibin),FILE=outfile,FORM='unformatted',
     &        STATUS='replace')
         CALL wdump(iout(ibin))
         CLOSE(iout(ibin))
         print *, "NEW"
         print *, iunique(bintots(nbin)-50:bintots(nbin))      
      END DO
!      CALL wdump(iout(j))
!      CALL wdump(iout(1))

      deallocate(bins,bintots,iout,countsink,tmp_listpm)
      END PROGRAM
!****************************************************************
!
      SUBROUTINE quit(i)
      
       STOP
      END SUBROUTINE

!
      SUBROUTINE endrun
       CALL quit(0)
      END SUBROUTINE

!
!     extract batch of default REALs
      subroutine extract_RD(array,binsize,binlist,bintotal)
        INCLUDE 'idim'
        real,dimension(:),intent(INOUT) :: array(idim)
        integer,intent(IN) :: binsize
        integer,intent(IN) :: binlist(binsize),bintotal
        real,dimension(:) :: tmp_arr(binsize)
        integer :: i
        print *, "extracting Default Real array"
!        print *, "Doing particles", binlist(1), binlist(bintotal)
        tmp_arr(:) = 0
        do i=1, bintotal
!           print *, i,array(binlist(i))
           tmp_arr(i) = array(binlist(i))
        end do
!        print *, "wiping array"
        array(:) = 0
!        print *, "array wiped"
        do i=1, bintotal
           array(i) = tmp_arr(i)
        end do
!       print *, "array reset"
      end subroutine extract_RD

!     extract batch of REAL*4
      subroutine extract_R4(array,binsize,binlist,bintotal)
        INCLUDE 'idim'
        real*4,dimension(:),intent(INOUT) :: array(idim)
        integer,intent(IN) :: binsize
        integer,intent(IN) :: binlist(binsize),bintotal
        real*4,dimension(:) :: tmp_arr(binsize)
        integer :: i
        print *, "extracting Real*4 array"
!        print *, "Doing particles", binlist(1), binlist(bintotal)
        tmp_arr(:) = 0
        do i=1, bintotal
!           print *, i,array(binlist(i))
           tmp_arr(i) = array(binlist(i))
        end do
!        print *, "wiping array"
        array(:) = 0
!        print *, "array wiped"
        do i=1, bintotal
           array(i) = tmp_arr(i)
        end do
!       print *, "array reset"
      end subroutine extract_R4

! Extract INT*8 array
      subroutine extract_I8(array,binsize,binlist,bintotal)
        INCLUDE 'idim'
        integer*8,dimension(:),intent(INOUT) :: array(idim)
        integer,intent(IN) :: binsize
        integer,intent(IN) :: binlist(binsize),bintotal
        integer*8,dimension(:) :: tmp_arr(binsize)
        integer :: i
        print *, "extracting I8 array"
        print *, "bintotal= ", bintotal
        tmp_arr(:) = 0
        do i=1, bintotal
!           print *, i,array(binlist(i))
           tmp_arr(i) = array(binlist(i))
        end do
        print *, "wiping array"
        array(:) = 0
        print *, "array wiped"
        do i=1, bintotal
           array(i) = tmp_arr(i)
        end do
        print *, "array reset"

      end subroutine extract_I8
        
!!
! Extract INTEGER*4 array
      subroutine extract_I4(array,binsize,binlist,bintotal)
        INCLUDE 'idim'
        integer*4,dimension(:),intent(INOUT) :: array(idim)
        integer,intent(IN) :: binsize
        integer,intent(IN) :: binlist(binsize),bintotal
        integer*4,dimension(:) :: tmp_arr(binsize)
        integer :: i
        print *, "extracting I4 array"
        tmp_arr(:) = 0
        do i=1, bintotal
!           print *, i,array(binlist(i))
           tmp_arr(i) = array(binlist(i))
        end do
!        print *, "wiping array"
        array(:) = 0
!        print *, "array wiped"
        do i=1, bintotal
           array(i) = tmp_arr(i)
        end do
!        print *, "array reset"
      end subroutine extract_I4

!
! Extract INTEGER*4 array
      subroutine extract_I1(array,binsize,binlist,bintotal)
        INCLUDE 'idim'
        integer*1,dimension(:),intent(INOUT) :: array(idim)
        integer,intent(IN) :: binsize
        integer,intent(IN) :: binlist(binsize),bintotal
        integer*1,dimension(:) :: tmp_arr(binsize)
        integer :: i
        print *, "extracting I1 array"
        tmp_arr(:) = 0
        do i=1, bintotal
!           print *, i,array(binlist(i))
           tmp_arr(i) = array(binlist(i))
        end do
!        print *, "wiping array"
        array(:) = 0
!        print *, "array wiped"
        do i=1, bintotal
           array(i) = tmp_arr(i)
        end do
!        print *, "array reset"
      end subroutine extract_I1

!
! Extract DEFAULT INTEGER array
      subroutine extract_ID(array,binsize,binlist,bintotal)
        INCLUDE 'idim'
        integer,dimension(:),intent(INOUT) :: array(idim)
        integer,intent(IN) :: binsize
        integer,intent(IN) :: binlist(binsize),bintotal
        integer,dimension(:) :: tmp_arr(binsize)
        integer :: i
        print *, "extracting ID array"
        tmp_arr(:) = 0
        do i=1, bintotal
!           print *, i,array(binlist(i))
           tmp_arr(i) = array(binlist(i))
        end do
!        print *, "wiping array"
        array(:) = 0
!        print *, "array wiped"
        do i=1, bintotal
           array(i) = tmp_arr(i)
        end do
!        print *, "array reset"
      end subroutine extract_ID

!
! handles all the sink variables
      subroutine extract_sinks(totsinks,ptm_inds)
        INCLUDE 'idim'
        INCLUDE 'COMMONS/ptmass'
        integer,intent(IN) :: totsinks
        integer,dimension(:) :: ptm_inds(iptdim)
        integer :: i
        listpm(:) = 0
        do i=1,totsinks
           listpm(i) = ptm_inds(i)
        end do
        nptmass = totsinks
        call extract_RD(spinx,iptdim,ptm_inds,totsinks)
        call extract_RD(spiny,iptdim,ptm_inds,totsinks)
        call extract_RD(spinz,iptdim,ptm_inds,totsinks)
        call extract_RD(angaddx,iptdim,ptm_inds,totsinks)
        call extract_RD(angaddy,iptdim,ptm_inds,totsinks)
        call extract_RD(angaddz,iptdim,ptm_inds,totsinks)
        call extract_RD(spinadx,iptdim,ptm_inds,totsinks)
        call extract_RD(spinady,iptdim,ptm_inds,totsinks)
        call extract_RD(spinadz,iptdim,ptm_inds,totsinks)
      end subroutine extract_sinks
