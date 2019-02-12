c Some subroutines that profile.F needs for krome=yes
c For use in profile.F
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
c makes array of indices of locations of selected molecules in
c mfracs array   
      SUBROUTINE get_mol_inds(mollist,maxmols,dumpmols,molinds)
      INCLUDE 'COMMONS/krome_mods'
      INCLUDE 'idim'
      INCLUDE 'COMMONS/kromevar'
      INTEGER maxmols,dumpmols,molinds(maxmols),j
      CHARACTER*16 mollist(maxmols)
      molinds = 0
      DO j=1,dumpmols
         print *, j, adjustl(mollist(j))
         molinds(j) = krome_get_index(adjustl(mollist(j)))
      END DO

      END SUBROUTINE get_mol_inds

c************************************************** 
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
      print *, "Counting molecules in ", molfile
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
      CLOSE(unit)
      print *, "Found", countmols, "molecules to read."
      RETURN
      END FUNCTION countmols



c****************************************************************               
      SUBROUTINE writekey(unit,nwrite,maxmols,molinds)
      INCLUDE "COMMONS/krome_mods"
      INCLUDE "idim"
      INCLUDE "COMMONS/kromevar"

      INTEGER unit,nwrite,i,maxmols,molinds(maxmols)
      CHARACTER*15 list(5)
      list = (/"r", "rho","rhovert","temp","tempvert"/)

      IF (usekrome .eq. 1) THEN
         OPEN(unit,FILE="krome_skey.txt")
         WRITE(unit,12) (i,list(i),i=1,5)
         WRITE(unit,13) (i+5,molnames(molinds(i)),i=1,nwrite)
         WRITE(unit,"(I2,1X,A15)") nwrite+6,"rad_vert"
         WRITE(unit,14) (i+6+nwrite,molnames(molinds(i)),i=1,nwrite)
 12            FORMAT(5(I2,1X,A15,/))
 13                  FORMAT(<nwrite>(I2,1X,A15,/))
 14                        FORMAT(<nwrite>(I2,1X,"vert",A15,/))
         CLOSE(unit)
      ELSE IF (usekrome .eq. 2) THEN
         OPEN(unit+1,FILE="krome_fkey.txt")
         WRITE(unit+1,12) (i,list(i),i=1,5)
         WRITE(unit+1,15) (i+5,molnames(i),i=1,krome_nmols)
         WRITE(unit+1,"(I3,1X,A15)") krome_nmols+5, "rad_vert"
         WRITE(unit+1,16) (i+5+krome_nmols,molnames(i),i=1,krome_nmols)
 15            FORMAT(<krome_nmols>(I3,1X,A15,/))
 16                  FORMAT(<krome_nmols>(I3,1X,"vert",A15,/))
         CLOSE(unit+1)
      END IF
      END SUBROUTINE
