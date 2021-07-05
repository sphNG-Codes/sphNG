      SUBROUTINE preset(ihopin)
c************************************************************
c                                                           *
c  This routine makes sure that everything not entered as   *
c     option is defined before integration begins.          *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nicil_subs'
#endif
      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'COMMONS/mpiall'
#endif

      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/nlim'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/btree'
      INCLUDE 'COMMONS/stop'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/recor'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/nextmpt'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptdump'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/ptsoft'
      INCLUDE 'COMMONS/soft'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/binfile'
      INCLUDE 'COMMONS/delay'
      INCLUDE 'COMMONS/useles'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/optbl'
      INCLUDE 'COMMONS/eostbl'
      INCLUDE 'COMMONS/dusttbl'
      INCLUDE 'COMMONS/gaptbl'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/updated'
      INCLUDE 'COMMONS/rbnd' ! for periodic_no_ghosts
      INCLUDE 'COMMONS/boundheight'
      INCLUDE 'COMMONS/andrea'
      INCLUDE 'COMMONS/astrcon'
      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/radtran3'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/pxpy'
      INCLUDE 'COMMONS/physcon'
#ifdef NONIDEAL
      INCLUDE 'COMMONS/nonideal'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/vargam'
#endif
      INCLUDE 'COMMONS/stellarradiation'
      INCLUDE 'COMMONS/lyonfiles'
      INCLUDE 'COMMONS/hosokawa'

      CHARACTER*7 where
      CHARACTER*6 dummy

      CHARACTER*100 metalfilenames, metalbase, homedir, filetable
      DIMENSION metalfilenames(opmxMet_ferguson)

      REAL*4 rho_real4

      REAL optable_line(opmxrh)

      DATA where/'preset'/

      INTEGER ihopin
#ifdef DGRID
      INTEGER dfilenum
      LOGICAL existingdragfile
      CHARACTER*4 dfilename
#endif

      IF(ihopin.NE.1) THEN
c
c--Open file1
c
#ifdef MPI
c      CALL appendiproc(7,file1,iproc)
         write (*,*) 'About to open ',file1
#endif

         OPEN (idisk1, FILE=file1, STATUS='unknown', 
     &        FORM='unformatted', RECL=imaxrec)
c
c--Read dump
c
         CALL rdump_wrapper(idisk1, ichkl, 0)
         CLOSE (idisk1)
         IF (ichkl.EQ.1) CALL error(where, ichkl)
c
c--For stellar evolution code, need to make file name to store stellar
c     evolution model that is appropriate for the *next* dump file.
c
         IF (istellarfeedback.EQ.6) THEN
            DO i = 1, nptmass
               file_stellarmodel_in(i) = file1
               file_stellarmodel_out(i) = file1
               CALL increment_dump_filename(file_stellarmodel_out(i))
            END DO 
            file_stellarmodel_new = file1
            CALL increment_dump_filename(file_stellarmodel_new)
         ENDIF
c
c--Open point mass data output file
c
         WRITE (ptfile,99990) namerun
99990    FORMAT ('P', A20)
         WRITE (accfile,99991) namerun
99991    FORMAT ('A', A20)
#ifdef MPI
c         CALL appendiproc(21,ptfile,iproc)
c         CALL appendiproc(21,accfile,iproc)
#endif
         IF (nptmasstot.NE.0 .OR. iptmass.NE.0) THEN
#ifdef MPIALL
            print *,iproc,': Preset ',ptfile
            IF (iproc.EQ.0) THEN
#endif
               OPEN (iptprint, FILE=ptfile,STATUS='replace',
     &              ACTION='write', FORM='unformatted')
               CLOSE (iptprint)
#ifdef MPIALL
            ENDIF
#endif
            OPEN (iaccpr,FILE=accfile,STATUS='replace',
     &           ACTION='write', FORM='unformatted')
         ENDIF
c
c--Open files for killing and reassignment of particles
c
         WRITE (killfile,99992) namerun
99992    FORMAT ('K', A20)
         WRITE (reassfile,99993) namerun
99993    FORMAT ('R', A20)
#ifdef MPI
c         CALL appendiproc(21,killfile,iproc)
c         CALL appendiproc(21,reassfile,iproc)
#endif
         IF (ibound.EQ.8 .OR. ibound/10.EQ.9 .OR. 
     &        ibound.GE.100 .AND. ibound.LT.109) THEN
            OPEN (ikillpr, FILE=killfile, STATUS='replace',
     &           ACTION='write', FORM='unformatted')
            OPEN (ireasspr, FILE=reassfile, STATUS='replace',
     &           ACTION='write', FORM='unformatted')
         ENDIF
         IF (irotpot.EQ.1)
     &        OPEN (ikillpr, FILE=killfile, STATUS='replace',
     &        ACTION='write', FORM='unformatted')
c
c--Open notify file
c
#ifdef MPIALL
         IF (iproc.EQ.0) THEN
#endif
            OPEN (inotify,FILE='notify')
#ifdef MPIALL
         ENDIF
#endif

      ENDIF

#ifdef MPIALL
      print *,iproc,': Preset B'
#endif

c
c--check ibound = 11 for periodic no ghosts
c
#ifdef PERIODIC_NO_GHOSTS
      IF (ibound.NE.11) THEN 
99901    FORMAT(' PERIODIC_NO_GHOSTS REQUIRES ibound=11')
         WRITE(*,99901)
         CALL quit(0)
      ENDIF
      dxbound = xmax - xmin
      dybound = ymax - ymin
      dzbound = zmax - zmin
#endif

c
c--Initialize viscosity switch
c
      IF (ifsvi.GE.6) THEN
         alphamin(1) = alpha
         alphamax(1) = beta
      ENDIF
      
      DO i = 1, npart
         IF (ifsvi.GE.6 .AND. abs(gt).LT.tiny .OR. 
     &        alphaMM(1,i).EQ.0.0) THEN
            DO k = 1, isizealphaMM
               alphaMM(k,i) = alphamin(k)
            END DO
c
c--always initialise resistivity (assumes constant value of alphamin
c  if ifsvi.lt.6)
c
         ELSEIF (imhd.EQ.idim .and. abs(gt).LT.tiny) THEN
            alphaMM(2,i) = alphamin(2)
         ENDIF
      END DO
#ifdef MPIALL
      print *,iproc,': Preset C'
#endif
c
c--Set constant for artificial viscosity
c
c      IF (alpha.EQ.0. .AND. beta.EQ.0.) THEN
c         alpha = 1.0
c         beta = 2.0
c      ENDIF
c
c--Set accuracy parameter for tree force calculation
c     theoretical limit for 3D tree is 0.57 = 1/SQRT(3)
c
c      acc = 0.7
      acc = 0.5
c      acc = 0.3
c      acc = 0.0
c
c--Set stop flag
c
      istop = 0
c
c--Set min and max limit of neighbours the code tries to inforce
c
      neimin = 30
      neimax = 70
      nrange = 12
c
c--Set mass for particle being partially accreted at which it is 
c     completely accreted
c
      pmassmin = 1.0E+10
      fmas1 = 0.
      fmas2 = 0.
      DO i = 1, npart
         IF (iphase(i).GE.0) THEN
            IF (iunique(iorig(i)).LE.n1) THEN
               fmas1 = fmas1 + xyzmh(4,i)
            ELSE
               fmas2 = fmas2 + xyzmh(4,i)
            ENDIF
            IF (pmassmin.GT.xyzmh(4,i)) pmassmin = xyzmh(4,i)
         ENDIF
      END DO
      pmassleast = pmassmin/100.
c
c--Point Mass Presets
c
c--Set critical density for point mass creation
c
      rhocrea = rhozero*ptmcrit
c
c--Standardise point mass types and accretion radii
c
#ifdef MPIALL
      print *,iproc,': Preset D ',nptmass,initialptm
#endif
      DO i = 1, nptmass

#ifdef MPIALL
         print *,iproc,': Preset Ds ',nptmass,i
         print *,iproc,': Preset Ds ',listpm(i)
#endif
         iptcur = listpm(i)
         IF (initialptm.LT.1 .OR. initialptm.GT.6) CALL error(where, 2)
         iphase(iptcur) = initialptm
cccc         xyzmh(5,iptcur) = hacc
      END DO
c
c--For Massive Accretion, only accrete mass and change ptmass properties
c     for certain particles - exclude these:
c
#ifdef MPIALL
      print *,iproc,': Preset Da '
#endif
      DO i = 1, idim
         notacc(i) = .FALSE.
         iupdated(i) = .FALSE.
      END DO
      nnotacc = 0
c      WRITE (iprint,*) '********* IDELAYACC = 1 *********'
c      OPEN (22, FILE='WANTRUN')
c 50   READ (22, *, END=100) inum
c      nnotacc = nnotacc + 1
c      notacc(inum) = .TRUE.
c      GOTO 50
c 100  CLOSE(22)
c
c--Gravitational softening for ptmass-ptmass interactions
c
      iptsoft = 1
      ptsoft = 0.02*au/udist
c
c--Define the merger radius (squared) for the sink particles (may be 
c     related to the gravitational softening radius defined above)
c
      pt_merger_radius2 = (6.0*solarr/udist)**2


#ifdef MPIALL
      print *,iproc,': Preset Db '
#endif
c
c--Set psoft for softening the gravitational potential when the 1/(r+psoft)
c     potential law is used (i.e. when igrape=1, or isoft=1).
c
      psoft = 0.01
c
c--Value of minimum h in order to save computing time (if hmin is different
c     from 0 then program does not follow high density regions accurately)
c
      hmin = 0.
c      hmin = 0.01

      IF (isoft.EQ.1 .AND. hmin.LT.psoft) hmin = psoft
c
c--Compute tables for kernel quantities
c
#ifdef MPIALL
      print *,iproc,': Preset E'
#endif
c Setpart has already called ktable for initial jobs
      IF (ihopin.NE.1) CALL ktable
c
c--Read in tables of quantities
c
      print *, 'Reading tables in preset'
      CALL GET_ENVIRONMENT_VARIABLE('SPH_HOME',homedir)
      IF (TRIM(homedir).EQ.'') THEN
        WRITE(*,99901)
      ENDIF
99901 FORMAT(' WARNING: SPH_HOME not set.')
c
c--Load in tables of opacity
c
      IF (encal.EQ.'r') THEN
c
c--Read in low temperature opacity model (used if last digit = 0)
c
         imodop = MOD(iopmodel,10)
         IF (imodop.EQ.0 .OR. imodop.EQ.1 .OR. imodop.EQ.2) THEN
c
c--This is also the default opacity model (for both high and low temperature)
c     The "opacitytbl" file here is to take the full opacity table
c     from the original table (Whitehouse & Bate 2006).  This table was later
c     split into the Pollack only bit and the Alexander only bit.  The
c     max values from each of these two bits gives back the original table.
c     But the two are split so that can allow different low-temperature and
c     high-temperature tables to be used.
c
c            OPEN(UNIT=8,FILE='/home/mrbate/tables/opacitytbl')
            filetable = TRIM(homedir)//
     &           'tables/opacitytbl_Pollacketal1985'
            OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
            DO i=1, opmxtg
               READ(itablerd,*) (optable(i,j), j=1, opmxrh)
            END DO
            CLOSE(itablerd)
         ELSE
            WRITE (*,*) 'ERROR - unknown low-T opacity model'
            CALL quit(0)
         ENDIF
c
c--If using imodop.EQ.2, then also read in opacities from other tables
c
         IF (imodop.EQ.2) THEN
            IF (idustRT.EQ.0) THEN
               WRITE (*,*) 'ERROR - need idustRT.NE.0 for non-grey RT'
               CALL quit(0)
            ENDIF
c            filetable = TRIM(homedir)//
c     &           'tables/opacitytbl_DL1984_012mu_albedo.txt'
            filetable = TRIM(homedir)//
     &           'tables/opacitytbl_DL1984_1mu_albedo.dat'
            OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
            READ(itablerd,*)
            READ(itablerd,*)
            DO i = 1, nopacity_lambda
               READ(itablerd,*) opacity_vs_lambda(1,i),
     &              opacity_vs_lambda(2,i), opacity_vs_lambda(3,i)
               opacity_vs_lambda(1,i) = LOG10(opacity_vs_lambda(1,i))
               opacity_vs_lambda(2,i) = LOG10(opacity_vs_lambda(2,i))
               opacity_vs_lambda(3,i) = LOG10(opacity_vs_lambda(3,i))
            END DO
            CLOSE(itablerd)
c
c--Test opacities
c
c            DO i = 1, 2000
c               v = c/(1.0E-05*(1.01)**i)
c               vv = Qv(v,Qabs)
c               WRITE (200,*) 1.0E-05*(1.01)**i,vv,Qabs
c            END DO
c            CALL quit(0)
         ENDIF
c
c--Read in high temperature opacity models
c
         IF (iopmodel/10.EQ.0) THEN
c
c--Read in Alexander (1975) models (iopmodel="0X").  
c     These are placed into the optable so that the full optable 
c     is the same as in the original version of
c     the code (i.e. a mix of Pollack et al. 1985 and Alexander 1975)
c
            ialexander_only = 34
            filetable = TRIM(homedir)//
     &           'tables/opacitytbl_Alexander1975_Model_KingIVa'
            OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
c
c--Skip 1st line
c
            READ(itablerd,*)
c
c--Read Alexander data
c
            DO i=opmxtg-ialexander_only+1, opmxtg
               READ(itablerd,*) (optable_line(j), j=1, opmxrh)
c
c--Take the maximum of the values already in the table (from 
c     Pollack et al. 1985) and those taken from the Alexander (1975) table
c
               DO j = 1, opmxrh
                  optable(i,j) = MAX(optable(i,j),optable_line(j))
               END DO
            END DO
            CLOSE(itablerd)
c
c--End of reading in Alexander data
c
         ELSEIF (iopmodel/10.EQ.1) THEN
c
c--Read in Ferguson et al. (2005), ApJ, 623, 585 models (iopmodel="1X").
c     Note these are provided as a tar file (f05.g93.tar).  
c     The last three lines of the file "g7.0.tron" need to be 
c     edited because the numbers are not separated by spaces, 
c     but other than that the data files are used
c     as provided by the tarfile.
c
            metaltable(1) = 0.0
            metaltable(2) = 0.00001
            metaltable(3) = 0.00003
            metaltable(4) = 0.0001
            metaltable(5) = 0.0003
            metaltable(6) = 0.001
            metaltable(7) = 0.002
            metaltable(8) = 0.004
            metaltable(9) = 0.01
            metaltable(10) = 0.02
            metaltable(11) = 0.03
            metaltable(12) = 0.04
            metaltable(13) = 0.05
            metaltable(14) = 0.06
            metaltable(15) = 0.08
            metaltable(16) = 0.1

            metalbase = TRIM(homedir)//'tables/Fergusonetal2005/'
            metalfilenames(1) = TRIM(metalbase) // 'g7.0.tron'
            metalfilenames(2) = TRIM(metalbase) // 'g7.00001.tron'
            metalfilenames(3) = TRIM(metalbase) // 'g7.00003.tron'
            metalfilenames(4) = TRIM(metalbase) // 'g7.0001.tron'
            metalfilenames(5) = TRIM(metalbase) // 'g7.0003.tron'
            metalfilenames(6) = TRIM(metalbase) // 'g7.001.tron'
            metalfilenames(7) = TRIM(metalbase) // 'g7.002.tron'
            metalfilenames(8) = TRIM(metalbase) // 'g7.004.tron'
            metalfilenames(9) = TRIM(metalbase) // 'g7.01.tron'
            metalfilenames(10) = TRIM(metalbase) // 'g7.02.tron'
            metalfilenames(11) = TRIM(metalbase) // 'g7.03.tron'
            metalfilenames(12) = TRIM(metalbase) // 'g7.04.tron'
            metalfilenames(13) = TRIM(metalbase) // 'g7.05.tron'
            metalfilenames(14) = TRIM(metalbase) // 'g7.06.tron'
            metalfilenames(15) = TRIM(metalbase) // 'g7.08.tron'
            metalfilenames(16) = TRIM(metalbase) // 'g7.1.tron'
c
c--Set Z (metal contribution) from metallicity in solar units (Z=0.02)
c
            xlogZ = MIN(-1.0,LOG10(MAX(metallicity*0.02,1.0E-10)))
            DO k = 1, opmxMet_ferguson
               metaltable(k) = LOG10(MAX(metaltable(k),1.0E-10))
               IF (xlogZ.GE.metaltable(k)) metaltable1 = k
            END DO
            metaltable1 = MIN(metaltable1,opmxMet_ferguson-1)

            DO k = 1, opmxMet_ferguson
            OPEN(itablerd,FILE=metalfilenames(k),
     &           STATUS='old',ACTION='read')
c
c--Skip initial lines, the read in temperature scale
c
            DO i=1,3
               READ(itablerd,*)
            END DO
            READ(itablerd,55001) dummy, (optable_ferguson(1,i,k),i=2,
     &           opmxR_ferguson)
55001       FORMAT(A6,19(1F7.3))
c55002       FORMAT(20(1PE12.3))
            Rmax_ferguson = optable_ferguson(1,opmxR_ferguson,k)
            Rmin_ferguson = optable_ferguson(1,2,k)
            Rdelta_ferguson = optable_ferguson(1,3,k)-
     &           optable_ferguson(1,2,k)
c
c--Read in data - ***NOTE*** that it stores the table BACKWARDS since the
c     table as provided is in order of decreasing temperature
c
            DO i=opmxtg_ferguson, 2, -1
               READ(itablerd,*) 
     &            (optable_ferguson(i,j,k), j=1, opmxR_ferguson)
            END DO
            CLOSE(itablerd)
            Tmax_ferguson = optable_ferguson(opmxtg_ferguson,1,k)
            Tmin_ferguson = optable_ferguson(2,1,k)
            Tdelta_ferguson = optable_ferguson(3,1,k)-
     &           optable_ferguson(2,1,k)

            END DO
         ELSEIF (iopmodel/10.EQ.0) THEN
            WRITE (*,*) 'Going with Alexander 1975 high-T model'
         ELSE
            WRITE (*,*) 'ERROR - unknown low-H opacity model'
            CALL quit(0)
         ENDIF
c
c--Read in dusttable, which is required to make the interface between
c     dust opacity and gas opacity smooth when using reduced dust opacities
c     using the "opdenom" parameter, and when using the Alexander and
c     Pollack tables (i.e. opacity model iopmodel=0).
c
         filetable = TRIM(homedir)//'tables/dusttbl2'
         OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
         DO i=1, duslen
            READ(itablerd,*) (dusttable(i,j), j=1, duswid)
         END DO
         CLOSE(itablerd)

      ENDIF
c
c--End reading of opacity tables
c
c--Read in equation of state tables
c      
      IF (encal.EQ.'r' .OR. 
     &     (encal.EQ.'i' .AND. (ibound.EQ.102 .OR. ibound.EQ.103))
     &     .OR. ABS(iresist).EQ.2 .OR. ABS(iresist).EQ.3) THEN

#ifdef ENDIANBIG
c         filetable = TRIM(homedir)//'tables/gasttbl'
         filetable = TRIM(homedir)//'tables/gasttbl_lowdens'
#else
      filetable = TRIM(homedir)//'tables/gasttbl_lowdens_little_endian'
#endif
         OPEN(itablerd,FILE=filetable,FORM='unformatted',
     &        STATUS='old',ACTION='read')
         DO i=1, tgmxu
            READ(itablerd) (tgtable(i,j), j=1, tgmxrh)
         END DO
         CLOSE(itablerd)

#ifdef ENDIANBIG
c         filetable = TRIM(homedir)//'tables/utbl'
         filetable = TRIM(homedir)//'tables/utbl_lowdens'
#else
         filetable = TRIM(homedir)//'tables/utbl_lowdens_little_endian'
#endif
         OPEN(itablerd,FILE=filetable,FORM='unformatted',
     &        STATUS='old',ACTION='read')
         DO i=1, umxrh
            READ(itablerd) (utable(j,i), j=1, umxt)
         END DO
         CLOSE(itablerd)
 
#ifdef ENDIANBIG
c         filetable = TRIM(homedir)//'tables/molmasstbl'
         filetable = TRIM(homedir)//'tables/molmasstbl_lowdens'
#else
         filetable = TRIM(homedir)//
     &        'tables/molmasstbl_lowdens_little_endian'
#endif
         OPEN(itablerd,FILE=filetable,FORM='unformatted',
     &        STATUS='old',ACTION='read')
         DO i=1, mumxu
            READ(itablerd) (mutable(i,j), j=1, mumxrh)
         END DO
         CLOSE(itablerd)
      ENDIF

      IF (ibound.EQ.100) THEN
         filetable = TRIM(homedir)//'tables/gaptbl'
         OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
         DO i=1, gaplen
            READ(itablerd,*) (gaptable(i,j), j=1, gapwid)
         END DO
         CLOSE(itablerd)
      ENDIF
c
c--Make sure planet position is set from the beginning
c
      IF (irotpot.EQ.1) THEN
         CALL planetpotential (px, py, pz, imigrate, rorbitmax,
     &        pmrate, rorbit_orig, gt)
         print *, 'PX PY PZ'
         print *, px, py, pz
         print *, imigrate, rorbitmax
         print *, pmrate, rorbit_orig
         print *, gt
      ENDIF
c
c--SECTION OF CODE ONLY FOR TESTING TABLES
c
c      IF (.TRUE.) THEN
      IF (.FALSE.) THEN
         temp_val = 1.0
         u_val_init = 8.314E+07/2.381*10.0
         u_fac = 1.01
         rho_val = 1.0e-16/udens
         rho_fac = 3.0
         DO i = 1, 30
            rho_val = rho_val*rho_fac
            rho_real4 = rho_val
            cv_val = getcv(rho_real4,u_val_init)
            u_val = cv_val*temp_val
            DO j = 1, 1600
               u_val = u_val*u_fac
               cv_val = getcv(rho_real4,u_val)
               op_val = getkappa(u_val,cv_val,rho_real4)
c              WRITE (80+i,55004) rho_val*udens,cv_val*uergg,
c     &           u_val/cv_val,op_val/umass*udist**2
c55004          FORMAT(4(1PE12.5))
            END DO
         END DO

         CALL quit(0)
      ENDIF
c
c--END SECTION OF CODE FOR TESTING TABLES
c
c--Andrea Fortier provided tables to introduce planetesimal heating
c  to planet calculations.
c
      e_planetesimal(:) = 0.0
c
c--Zero planetesimal-planet collision energy array (not needed)
c
      pcolen(:) = 0.0

#ifdef DGRID
c
c--Setup table for recording planetesimal energy deposition
c
      print *, 'IPLANETESIMALS = ', iplanetesimals
      IF (gasdrag.AND.iplanetesimals.GT.0) THEN
         CALL dragtable

         READ (file1, 99303) dfilename, dfilenum
99303    FORMAT (A3, I4)
99304    FORMAT (A4, I4.4)
         dfilename = 'drag'
         WRITE (dragfile, 99304) dfilename, dfilenum

         INQUIRE (FILE=dragfile, exist=existingdragfile)
         IF (.NOT.existingdragfile .AND. dfilenum.GT.0) THEN
            WRITE (iprint,*) 'ERROR READING DRAGFILE'
            WRITE (iprint,*) 'No corresponding drag file for dump,'
            WRITE (iprint,*) 'but dump not No. 0000, i.e. first'
            CALL quit(0) 
         ENDIF

         IF (existingdragfile) THEN
            OPEN (idragfile, FILE=dragfile, STATUS='unknown', 
     &           FORM='unformatted', RECL=imaxrec)
            print *, 'READING DRAG FILE'
            CALL read_dragtable
         ENDIF
      ENDIF
#endif

      IF (iandrea.EQ.1) THEN
         filetable = TRIM(homedir)//'tables/planetesimals'
         OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
         READ (itablerd,*) prmin, princ
         READ (itablerd,*) pmmin, pminc
         pmmin = pmmin*5.9742E27/umass
         pminc = pminc*5.9742E27/umass
         DO i = 1, mbins
            READ (itablerd,*) (ptable(i,j), j = 1, rbins)
         END DO
         CLOSE(itablerd)
         DO i = 1, mbins
            DO j = 1, rbins
c--This is column 2 from Andrea's work. Despite here original notes
c  stating units of erg/g, they are actually erg/(g.s) (see e-mail
c  correspondence, 4th March 2011). They are also not logged as
c  originally stated.
               ptable(i,j) = (ptable(i,j))*utime/uergg
            END DO
         END DO
      ENDIF

      IF ((ibound.EQ.101 .OR. ibound.EQ.102 .OR. ibound.EQ.103) 
     &     .AND. encal.EQ.'r') THEN
         CALL discboundary(ibound)
         print *, 'Disc optically thin boundary set'
      ENDIF

      IF (ihopin.NE.1) THEN
c
c--Build table for choosing inflow into planet/disc simulation
c 
         IF (ibound.EQ.100 .AND. iexf.EQ.7) THEN
            CALL cumradtable
c
c--Build table for choosing inflow into planet/disc simulation from ZEUS
c     output
c 
            CALL zeusread
            print *, 'Building ZEUS table for mass', hmass
            CALL zeustable
         ENDIF
      ENDIF

#ifdef RT
c
c--Build tables for diffuse interstellar medium model (Bate & Keto 2015)
c
      IF (idustRT.GT.0) THEN
c
c--Build table of sightlines on a sphere
c
         CALL buildsightlines
c
c--Build tables for interstellar radiation field, photoelectric heating,
c     and the Planck mean dust opacity
c
         CALL build_interstellar_radiation_table
         CALL build_photoelectric_table
         CALL build_planck_mean_table
      ENDIF
#else
c
c--Allow ionisation to be used without radiative transfer
c
      IF (istellar_ionisation) THEN
         metallicity = 1.
      ENDIF
#endif
c
c--Build tables for direct stellar radiation from sink particles
c
      IF (istellar_radiation) THEN
         IF (istellar_ionisation .AND. idim_ionisation.NE.idim) THEN
            WRITE (*,*) 'ERROR - stellar ionisation TRUE but ',
     &           'idim_ionisation.NE.idim'
            CALL quit(0)
         ENDIF
         CALL build_kernel_column
         CALL build_kernel_n2
         CALL build_stellar_radiation_table
         IF (istellar_ionisation)
     &        CALL build_photoionisation_heating_rate_table
c         CALL build_radpres_table
c
c--Initialise *arrays*
c
         stellarsum = 0.
         istellar_location = 1
c
c--Read Hosokawa stellar evolution tables for accreting protostars
c
         filetable = TRIM(homedir)//
     &        'tables/Hosokawa/md3.dat'
         OPEN(itablerd,FILE=filetable,STATUS='old',ACTION='read')
         READ(itablerd,*)
         nhosokawa = 0
         DO i = 1, nhosokawa_max
c
c--Read stellar mass, radius, luminosity & total luminosity (+accretion)
c
            READ(itablerd,*,END=240) dummy,hosokawa(1,i),
     &           hosokawa(2,i),dummy,hosokawa(3,i),hosokawa(4,i)
            nhosokawa = nhosokawa + 1
         END DO
 240     CLOSE(itablerd)
         IF (nhosokawa.EQ.0) THEN
            WRITE (*,*) 'ERROR - nhosokawa.EQ.0'
            CALL quit(0)
         ENDIF
      ENDIF
c
c--Initialise non-ideal MHD
c
#ifdef NONIDEAL
      CALL nicil_initialise(utime,udist,umass,umagfd,ierr)
      IF (ierr/=0) CALL quit(0)
      CALL nimhd_init_gastemp(Rg,gamma,gam,gamdh,gamah,uergg
     &                       ,rhocrit,rhocrit2,rhocrit3)
      eta_nimhd  = 0.0
      nden_nimhd = 0.0
      WRITE (nimhdfile,99994) namerun
99994 FORMAT ('N', A20)
#ifdef MPI
      iproc_in = iproc
#else
      iproc_in = 0
#endif
      CALL nimhd_write_evheader(nimhdfile,iproc_in)
#endif

      IF (itrace.EQ.'all') WRITE(iprint,250)
 250  FORMAT(' exit subroutine preset')
#ifdef MPIALL
      print *,iproc,': Preset F'
#endif

      RETURN
      END

      SUBROUTINE appendiproc(numchar,string,iproc)

      CHARACTER*25 string

      IF (iproc.LT.10) THEN
         WRITE (string,99001) string(1:numchar),iproc
99001    FORMAT(A7,'.00',I1)
      ELSEIF (iproc.LT.100) THEN
         WRITE (string,99002) string(1:numchar),iproc
99002    FORMAT(A7,'.0',I2)
      ELSEIF (iproc.LT.1000) THEN
         WRITE (string,99003) string(1:numchar),iproc
99003    FORMAT(A7,'.',I3)
      ELSE
         WRITE (*,*) 'ERROR: string ',iproc
         CALL quit(0)
      ENDIF

      RETURN
      END
