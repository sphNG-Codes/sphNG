
!###################################################
! Calls krome for Tgas, rho, timestep,Av entered by user
!  !Use abundances relative to H!
!  
program krome_once
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  use krome_getphys
  implicit none
  integer,parameter::nsp=krome_nmols !number of species (common)
  real*8,parameter :: dust2gas=0.01d0
  real*8::Tgas,dt,x(nsp),rho,spy,av,time,n(nsp),masses(nsp)
  real*8 :: dtk,maxtime,time_yrs,scaleconst,newtime,rhoj(nsp)
  character*128 :: outfilename,foutfilename
  character*24 :: molnames(krome_nspec),name
  character*24 :: Tchar, rhochar,dtchar,FMT
  integer :: j,niter,k,Hnuc,Tint,i,ind
  real*8 :: crate,xprint(nsp)
  logical :: exptstep
  character(len=1) :: steptype 

  spy = 3.65d2 * 2.4d1 * 3.6d3 !seconds per year
  
  call krome_init() !init krome (mandatory)
 
  masses(:) = krome_get_mass()

  print *, "Enter Tgas"
  read (*,*) Tgas
  print *, "Enter rho"
  read (*,*) rho
  print *, "linear or exponential time steps?(l/e)?"
  read (*,*) steptype
  if (steptype.EQ."l") THEN
     exptstep = .False.
  else if (steptype.EQ."e") THEN
     exptstep = .True.
  else
     Print *, "steptype not recognised. Stop."
     STOP
  end if
  print *, "Enter Max time (yrs)"
  read (*,*) maxtime
  print *, "Enter no. time steps"
  read (*,*) niter
  print *, "Enter Av"
  read (*,*) av
  print *, "Enter CR ionization rate"
  read (*,*) crate
  call krome_set_user_crate(crate)
!  Tgas = 1d1 !gas temperature (K)
!  dt = 2d4 * spy !time-step (s)
!  rho = 1d-18 !gas density (g/cm3)
!  scaleconst = maxtime/10d0

!  x(:) = x(:) / sum(x) !normalize
! Reboussin 2014 abundances are rel to H nuclei no. density
!  do j=1, nsp
!     rhoj(j) = (rho/masses(krome_idx_H)) * masses(j) * x(j)
!  end do
!  do j=1,nsp
!     x(j) = rhoj(j)/sum(rhoj)
!  end do

  print *, "Using", Tgas, rho
  
  call get_init_abunds(krome_nmols,x(:),dust2gas,rho)
! ** don't SPECIFY GRAIN ABUNDANCE IN INPUT ABUNDANCES **
  PRINT *, "x(GRAIN-)=", x(krome_idx_GRAINk)

  Tint = INT(Tgas)
  write(Tchar,"(I0)") Tint
  write(rhochar,"(e8.2)") rho
  print *, rhochar
  write(dtchar,"(I4)") niter
  
  write ( *,*) "strings:", Tchar, rhochar,dtchar
!  outfilename= "n"//trim(Tchar)//"K_"//trim(rhochar)//"cm3_"//trim(adjustl(dtchar))
  outfilename= "n_"//trim(Tchar)//"K_"//trim(rhochar)//"cm3"
  foutfilename = "x_"//trim(Tchar)//"K_"//trim(rhochar)//"cm3"

  call krome_set_user_Av(av)
  CALL krome_set_user_Tdust(Tgas)
  call krome_set_user_rhogas(rho)
  !call the solver
  molnames = get_names()
  OPEN(UNIT=1,FILE=outfilename,FORM='formatted', &
       STATUS='replace')
  OPEN(UNIT=3,FILE=foutfilename,FORM='formatted', &
       STATUS='replace')
  OPEN(UNIT=11,FILE="molnumbers.dat",FORM='formatted', &
       STATUS='replace')

  DO j=1,krome_nmols
     WRITE(11,101) j, molnames(j)
  END DO
  CLOSE(11)
101 FORMAT (I5," ",a16)
  WRITE(FMT,'("(",I7,"e24.8)")') krome_nmols+1
  time=0.
  time_yrs = 0
  ! write initial abundances
  do j=1, nsp
     n(j) = rho *x(j)/masses(j)
     print *, molnames(j),x(j),masses(j)
  end do
  
  WRITE (1,FMT) time_yrs,(n(j),j=1,krome_nmols)
  DO k=1, niter
     IF (exptstep) THEN
        dt = spy*100d0* 10d0**(k * 60d0/niter)
     ELSE
        dt = maxtime * spy / niter
     END IF
     newtime = time + dt
     print *, "dt=", dt/spy
     PRINT *, "before:", "H:", x(krome_idx_H), "H2:", x(krome_idx_H2)
!     Print *, "initial x(:)"
 !    print *, x
     WRITE(*,*) "calling krome", x(krome_idx_H),rho,Tgas,dt/spy,Av
     call krome(x(:),rho,Tgas, dt) !call KROME
     time = newtime
     time_yrs = time/spy
     print *, "Time = ", time_yrs

     do j=1, nsp
        if (x(j) .lt. 1d-40) then
           xprint(j) = 0d0
        else
           xprint(j) = x(j)
        end if
        n(j) = rho *xprint(j)/masses(j)
     end do
     print *, "H:", x(krome_idx_H), "H2:", x(krome_idx_H2)
!   dt , mol1, mol2, mol3 ...
     WRITE (1,FMT) time_yrs,(n(j),j=1,krome_nmols)
     WRITE (3,FMT) time_yrs,(xprint(j),j=1,krome_nmols)
     IF (time_yrs > maxtime) then
        EXIT
     END IF
  END DO
!  DO j=1,krome_nmols
!     WRITE(1,100) molnames(j), x(j) 
!  END DO
!100 FORMAT (a16 e15.8)

  CLOSE(1)
  CLOSE(3)
  print *,"Test OK!"

end program krome_once

      SUBROUTINE get_init_abunds(nmols,x,dust2gas,rho)
        INCLUDE "COMMONS/krome_mods"
        IMPLICIT NONE
        INTEGER :: nmols,i
        REAL*8 :: x(nmols)
        REAL*8,intent(in) :: dust2gas,rho
        CHARACTER*20 :: abundfile,mfracfile
        LOGICAL :: haveAbundfile,haveMfracfile

        abundfile = "abundances.dat"
        mfracfile = "init_mfracs.dat"
        DO i=1,nmols
           x(i) = 0d0
        END DO
        INQUIRE(FILE=abundfile,EXIST=haveAbundfile)
        INQUIRE(FILE=mfracfile,EXIST=haveMFracfile)
        IF (haveMfracfile) THEN
           WRITE (*,*) "init_mfracs.dat found. Reading..."
           CALL read_mfracs(nmols,x,mfracfile)
        ELSE IF (haveAbundfile) THEN
           WRITE(*,*) "abundances.dat found. Reading..."
           CALL read_abunds(nmols,x,abundfile,rho,dust2gas)
        ELSE
           WRITE(*,*) "WARNING:Using default abundances"
          x(krome_idx_H) = 0.73
          x(krome_idx_He) = 0.279
          x(krome_idx_Cj) = 7.13d-4
          x(krome_idx_O) = 1.42d-3
!         init_x(krome_idx_Si) = 6.78d-5      
        ENDIF
         
      END SUBROUTINE get_init_abunds

!--------------------------------------------------------
!       Read abundances relative to TOTAL H nuclei

      SUBROUTINE read_abunds(nmols,x,abundfile,rho,dust2gas)
         INCLUDE "COMMONS/krome_mods"
         IMPLICIT NONE
         INTEGER :: nmols
         REAL*8,INTENT(IN) :: rho,dust2gas
         REAL*8 :: x(nmols),abunds(nmols),abund,masses(nmols)
         REAL*8 :: xtot, cations
         CHARACTER*20 :: abundfile,name,molnames(nmols)
         CHARACTER*30 :: kromename
         INTEGER :: iunit,iostat,ind,j,isgrain,iscation

         iunit=16
         iostat=0
         molnames = krome_get_names()
         masses = krome_get_mass()
         WRITE(*,*) "No mols:", nmols

         abunds(1:nmols) = (/ (0d0,j=1,nmols) /)
         OPEN(UNIT=iunit,FILE=abundfile,FORM='formatted', &
          STATUS='old')
         DO
            READ(iunit,*,IOSTAT=iostat) name, abund
            IF (name(1:1) == "!") THEN
               CYCLE
            END IF
            WRITE(*,*) "reading:", name,abund
            IF (iostat > 0) THEN
               WRITE(*,*) "Error reading abundances"
               STOP
            ELSE IF (iostat < 0) THEN
               EXIT
            ELSE
               WRITE(*,*) name, abund
               ind = krome_get_index(trim(name))
               abunds(ind) = abund
!               init_x(index) = abund
            END IF
         END DO
         CLOSE(iunit)

         ! SET GRAIN MASS FRACTIONS
         IF ((x(krome_idx_GRAIN0) .EQ. 0d0) .AND. (x(krome_idx_GRAINk) &
              .EQ. 0d0) ) THEN
            PRINT *, "GRAIN ABUNDANCES = 0, setting through dust2gas"
            x(krome_idx_GRAINk) = dust2gas
            x(krome_idx_GRAIN0) = 0d0
            abunds(krome_idx_GRAINk) = x(krome_idx_GRAINk) * &
                 masses(krome_idx_H) /(masses(krome_idx_GRAINk) *0.73d0)
            PRINT *, "GRAIN- x, n:", x(krome_idx_GRAINk), &
                 abunds(krome_idx_GRAINk)
         ELSE
            PRINT *, "GRAIN- fraction = ", x(krome_idx_GRAINk)
         END IF
         
         xtot = 0.d0
         cations = 0.d0

         !SET OTHER MASS FRACTIONS
         PRINT *, "Calculating mass fractions from abundances"
         DO j=1,nmols
            IF (index(molnames(j),"GRAIN") > 0) THEN
               CYCLE
            END IF
            x(j) = abunds(j) * masses(j)
           ! print *, "xj=", x(j)
            xtot = xtot + x(j)
            iscation = index(molnames(j),"+")
            IF (iscation > 0 .AND. abunds(j) >0d0) THEN
              PRINT *, "Initial cation:", molnames(j)
               cations = cations + abunds(j)
            END IF
         END DO
         PRINT *, "xtot, cations", xtot, cations

         !Calculate initial electron abundance
         abunds(krome_idx_E) = cations - abunds(krome_idx_GRAINk)
         IF ( abunds(krome_idx_E) < 0d0) THEN
           PRINT *, "ERROR n(E) is negative!!"
           STOP
         END IF
         PRINT *, "n(E)=", abunds(krome_idx_E)
         x(krome_idx_E) = abunds(krome_idx_E) * masses(krome_idx_E)
         PRINT *, "x(E)=", abunds(krome_idx_E)
         xtot = xtot + x(krome_idx_E)

         print *, "sum x = ", xtot
         DO j=1,nmols
            IF (index(molnames(j),"GRAIN") > 0 ) THEN
               PRINT *, "skipping grain", molnames(j) 
               CYCLE
            END IF
            x(j) = 0.99d0 * x(j) / xtot
         END DO
         !test loop
         xtot = 0d0
         DO j=1, nmols
            xtot = xtot + x(j)
         END DO
         PRINT *, "TOTAL should == 1:", xtot
      END SUBROUTINE read_abunds

!---------------------------------------------------
      SUBROUTINE read_mfracs(nmols,x,mfracfile)
         INCLUDE "COMMONS/krome_mods"
         IMPLICIT NONE
         INTEGER :: nmols
         REAL*8 :: x(nmols),frac
         CHARACTER*20 :: mfracfile,name,molnames(nmols)
         CHARACTER*30 :: kromename
         INTEGER :: iunit,iostat,ind,j

         iunit=17
         iostat=0
         molnames = krome_get_names()

         WRITE(*,*) "No mols:", nmols
!         WRITE(*,*) molnames                                                        
         OPEN(UNIT=iunit,FILE=mfracfile,FORM='formatted', &
          STATUS='old')
         DO
            READ(iunit,*,IOSTAT=iostat) name, frac
            IF (name(1:1) == "!") THEN
               CYCLE
            END IF
!            WRITE(*,*) "reading:", name,abund                                       
            IF (iostat > 0) THEN
               WRITE(*,*) "Error reading abundances"
               STOP
            ELSE IF (iostat < 0) THEN
               EXIT
            ELSE
!               WRITE(*,*) name, abund  
               ind = krome_get_index(name)
               x(ind) = frac
            END IF
         END DO
         CLOSE(iunit)
       END SUBROUTINE read_mfracs
