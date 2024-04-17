      SUBROUTINE makedust
c************************************************************
c                                                           *
c  This subroutine makes dust particles for 2-fluid,        *
c     gas and dust particle, simulations.                   *
c                                                           *
c     It can read in an existing dump file and bring        *
c     dead particles to life as dust particles.             *
c     Originally used to add dust to a settled disc.        *
c                                                           *
c     Alternately:                                          * 
c     Also allows two dump files, one containing gas and    *
c     the other containing dust particles to be combined.   *
c     The was used for Bate & Loren-Aguilar (2017), but     *
c     not added to the repository until 19 Oct 2022.        *
c                                                           *
c     Following Bate & Loren-Aguilar (2017), in 2023 the    *
c     routine was also enhanced to allow dust particles     *
c     to be set up at their terminal velocity initially     *
c     (instead of with zero velocity relative to the gas).  *
c     This requires two 'gas' particle dump files to be     *
c     created with idustIMPL = 1 and imakedust=.TRUE.       *
c     in 'idim'.  This is so that the dump file that will   *
c     be converted into dust particles has the pressure     *
c     gradients for all gas particles computed and written  *
c     to the dump file.  These are then used when the two   *
c     dump files are combined (this routine) to set the     *
c     the dust at its terminal velocity relative to the     *
c     gas.  The 'gas' and 'dust' dump files should be       *
c     identical, except that the numbers of particles can   *
c     differ (would usually expect Ngas >> Ndust to avoid   *
c     artificial dust particle clumping).  In particular,
c     the total masses of the two clouds should be the      *
c     same -- the dust-to-gas ratio is set in this makedust *
c     routine when the two dump files are combined.         *
c     Although using the standard M4 kernel works, the      *
c     pressure gradients are best calculated using the      *
c     Quintic SPH kernel rather than the standard M4 kernel.*
c     It's also best to run two-fluid dust with the         *
c     Quintic kernel anyway, so this is consistent.         *
c     Note that imakedust=.TRUE. should be set back to      *
c     .FALSE. for the executable that is used to actually   *
c     run the calculation (there is a check for this).      *
c                                                           *
c************************************************************

#ifdef MPIALL
#include "mpi_sup.h"
#endif
      INCLUDE 'idim'

#ifdef MPI
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/units'
      INCLUDE 'COMMONS/part'
      INCLUDE 'COMMONS/densi'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/gtime'
      INCLUDE 'COMMONS/bodys'
      INCLUDE 'COMMONS/ener2'
      INCLUDE 'COMMONS/ener3'
      INCLUDE 'COMMONS/fracg'
      INCLUDE 'COMMONS/polyk2'
      INCLUDE 'COMMONS/phase'
      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/sort'
      INCLUDE 'COMMONS/tming'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/gradhterms'
      INCLUDE 'COMMONS/varmhd'
      INCLUDE 'COMMONS/presb'
      INCLUDE 'COMMONS/Bxyz'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/files'
      INCLUDE 'COMMONS/initpt'
      INCLUDE 'COMMONS/rbnd'
      INCLUDE 'COMMONS/planetesimal'
      INCLUDE 'COMMONS/logun'
      INCLUDE 'COMMONS/grain'
      INCLUDE 'COMMONS/physcon'
      INCLUDE 'COMMONS/dustfluidvelu'

      CHARACTER*11 ifile, ifile_dust, ofile
      CHARACTER*20 name_in
      CHARACTER*1 igeom, icombine, iref, ireset_time, iterminal
      
1000  FORMAT (A7)
1001  FORMAT (A1)
1002  FORMAT (A11)
1020  FORMAT (A20)

#ifdef MPI
c
c--IF MPI setup, then read set up input from file
c
c      iread = 10
      OPEN (iread, FILE='makedust.txt', STATUS='unknown')
#endif
c
c--assume that MHD variable is B
c
      IF (imhd.EQ.idim) THEN
         varmhd = 'Bvol'
      ENDIF

      imax = 1073741824
      imaxstep = imax/2

      nactive = 0
      DO i = 1, idim
         isort(i) = i
         iorig(i) = i
      END DO
c
c============================
c--Determine how to make dust
c============================
c
      PRINT *, 'Do you want to combine separate gas and dust files',
     &     ' or add dust to a gas file (c/a) ?'
      READ (iread, 1001) icombine

      IF (icombine.EQ.'c' .OR. icombine.EQ.'C') THEN
c
c=====================================================================
c--Combine two existing dump files -- one for the gas and one for dust
c=====================================================================
c
      PRINT *,' Enter name of ASCII log file for running code'
      READ (iread, 1020) name_in

      PRINT *, 'Name of binary dump file containing gas?'
      PRINT *, 'NOTE: Should contain total (gas+dust) mass'
      READ (iread, 1000) ifile

      PRINT *, 'Name of binary dump file containing dust?'
      PRINT *, 'NOTE: Should contain total (gas+dust) mass'
      PRINT *, '      (i.e. may be different number of SPH'
      PRINT *, '      particles than the gas file, but should'
      PRINT *, '      be identical otherwise'
      READ (iread, 1000) ifile_dust

      PRINT *,'Enter dust-to-gas ratio (usually <1)'
      READ (iread,*) dust_to_gas_ratio

      PRINT *, 'Name of output binary dump file ?'
      READ (iread, 1000) ofile

      print *,iread, ifile, ifile_dust, ofile, name_in

      PRINT *, 'Name of corresponding ASCII ifile ?'
      READ (iread, 1002) inname

      print *,iread,inname

      PRINT *, 'About to open ',ifile

      OPEN (idisk1, FILE = ifile, FORM = 'unformatted')
      PRINT *, 'reading file ', ifile
c
c--Read dump file
c
      CALL options

      file1 = ofile
      namenextrun = name_in
      print *, 'Output files =', file1, namenextrun

      CALL rdump(idisk1,ichkl,1)
      CLOSE (idisk1)
      IF (ichkl.EQ.1) THEN
         PRINT*, 'ERROR READING DUMP FILE'
         CALL quit(0)
      ENDIF

      npart_gas = npart

      PRINT *, 'Got ',npart_gas,' gas particles'

      PRINT *, 'About to open ',ifile_dust

      OPEN (idisk1, FILE = ifile_dust, FORM = 'unformatted')
      PRINT *, 'reading file ', ifile_dust

      CALL rdump(idisk1,ichkl,npart)
      CLOSE (idisk1)
      IF (ichkl.EQ.1) THEN
         PRINT*, 'ERROR READING DUMP FILE'
         CALL quit(0)
      ENDIF

      PRINT *,'Got ',npart-npart_gas,' dust particles'
      PRINT *,'Total particles = ',npart
c
c--End reading of dump files
c---------------------------
c
c--Modify masses of gas to allow for dust mass
c
      rho_gas_fac = 1.0/(1.0+dust_to_gas_ratio)
      xyzmh(4,1:npart_gas) =xyzmh(4,1:npart_gas)*rho_gas_fac
c
c--Turn particles from second file to dust
c
      iphase(npart_gas+1:npart) = 11
      rho_dust_fac = (1.0-1.0/(1.0+dust_to_gas_ratio))
      xyzmh(4,npart_gas+1:npart) = xyzmh(4,npart_gas+1:npart)*
     &     rho_dust_fac

      min_rplan = 1.0
      max_rplan = 3.0

      PRINT *,'Enter the minimum and maximum radii for the dust',
     &     ' (e.g. ',min_rplan,max_rplan,' )'
      READ (iread,*) min_rplan, max_rplan

      PRINT *,'Do you want cylindrical or spherical radius (c/s)?'
      READ (iread,1001) igeom

      gas_mass_inside = 0.
      dust_mas_inside = 0.
      ngas_inside = 0
      ndust_inside = 0
      DO i = 1, npart
         IF (igeom.EQ.'s') THEN
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
         ELSE
            r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
         ENDIF
         IF (iphase(i).EQ.0) THEN
            IF (SQRT(r2).GT.min_rplan.AND.SQRT(r2).LT.max_rplan) THEN
               gas_mass_inside = gas_mass_inside + xyzmh(4,i)
               ngas_inside = ngas_inside + 1
            ENDIF
         ELSEIF (iphase(i).EQ.11) THEN
            IF (SQRT(r2).GT.min_rplan.AND.SQRT(r2).LT.max_rplan) THEN
               dust_mass_inside = dust_mass_inside + xyzmh(4,i)
               ndust_inside = ndust_inside + 1
            ELSE
               iphase(i) = -1
            ENDIF
         ENDIF
      END DO
      PRINT *,''
      PRINT *,'Total gas mass inside dust radii ',gas_mass_inside
      PRINT *,'Total gas particles inside dust radii ',ngas_inside
      PRINT *,'Total dust mass inside dust radii ',dust_mass_inside
      PRINT *,'Total dust particles inside dust radii ',ndust_inside
      PRINT *,''
c
c--Set terminal velocities if required
c
      PRINT *, 'Do you want to set dust to initial terminal velocity?'
      READ (iread, 1001) iterminal
      print *,'iterminal = ',iterminal
      IF (iterminal.EQ.'Y' .OR. iterminal.EQ.'y') THEN
         IF (encal.EQ.'i' .OR. encal.EQ.'t') THEN
            vthermal_coef = SQRT(8.0/pi)
         ELSE
            vthermal_coef = SQRT(8.0/pi/gamma)
         ENDIF
         print *,'terminal ',vthermal_coef
c
c--Store prgrad/rho (acceleration due to pressure gradient) values from
c     dust dump file and thermal velocities
c
c         DO i=1, npart
c            IF (iphase(i).GT.10) THEN
c               dustfluidvxyz_store(1:3,i) = dustfluidvxyz(1:3,i)
c            ENDIF
c         END DO
c
c--Need to call hcalc to get thermal velocities for stopping time
c
c         CALL hcalc
c
c
c
         IF (n_planetesimal_types.EQ.0) THEN
            WRITE(*,*) 'ERROR - you are using fixed-size ',        
     &           'implicit dust-drag but have not defined ',
     &           'any dust grain size in your ASCII file'
            CALL quit(0)
         ENDIF
c
c--Set dust terminal velocities (ASSUMES dist and gas have zero relative
c     velocity prior to this)
c
         DO i = 1, npart
            IF (iphase(i).GT.10) THEN
               IF (iphase(i)-10.LE.n_planetesimal_types) THEN
                  print *,'planet ',i,iphase(i)
                  rho_planetesimal = rho_grain(iphase(i))
                  print *,'planet rho ',rho_planetesimal
                  r_planetesimal = r_planetesimals(iphase(i)-10)
                  print *,'planet rad ',r_planetesimal
               ELSE
                  WRITE(*,*)'ERROR - Not enough n_planetesimal_types'
                  CALL quit(0)
               ENDIF
c
c--Stopping time -- ASSUMES that the DUST particle input file is gas
c     to start with so that dustfluidvxyz() holds the acceleration due
c     to the pressure gradient, and vsound() is defined so that it gives
c     the gas temperature at the dust particle's location
c
               stopt = r_planetesimal*rho_planetesimal/
     &              (thermal_velocity(i)*rho(i))
         print *,'terminal ',stopt,r_planetesimal,rho_planetesimal
               vxyzu(1:3,i) = vxyzu(1:3,i) +dustfluidvxyz(1:3,1,i)*stopt
               rho(i) = rho(i)*rho_dust_fac
            ELSEIF (iphase(i).EQ.0) THEN
               rho(i) = rho(i)*rho_gas_fac
            ENDIF
         END DO
      ENDIF
c
c=======================================
c--Add dust particles to a gas-only file
c=======================================
c
      ELSE

      PRINT *, 'Name of binary dump file to add dust to ?'
      READ (iread, 1000) ifile

      print *,iread,ifile

      PRINT *, 'Name of output binary dump file ?'
      READ (iread, 1000) ofile

      print *,iread,ofile

c      OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

      PRINT *, 'Name of corresponding ASCII ifile ?'
      READ (iread, 1002) inname

      PRINT *, 'Do you want to reset time to (almost) zero ?'
      READ (iread, 1001) ireset_time

      print *,iread,inname

      PRINT *, 'About to open ',ifile

      OPEN (idisk1, FILE = ifile, FORM = 'unformatted')
      PRINT *, 'reading file ', ifile
c
c--Read dump file
c
      CALL options

      file1 = ofile
      print *, 'Output file=', file1

      CALL rdump_wrapper(idisk1,ichkl,1)
      CLOSE (idisk1)
      IF (ichkl.EQ.1) THEN
         PRINT*, 'ERROR READING DUMP FILE'
         CALL quit(0)
      ENDIF
c
c--End reading of dump file
c--------------------------
c
      DO i = 1, npart
         IF (iphase(i).GE.0) nactive = nactive + 1
      END DO

      min_rplan = 1.0
      max_rplan = 3.0

      PRINT *,'Enter the minimum and maximum radii for the dust',
     &     ' (e.g. ',min_rplan,max_rplan,' )'
      READ (iread,*) min_rplan, max_rplan

      PRINT *,'Do you want cylindrical or spherical radius (c/s)?'
      READ (iread,1001) igeom
c
c--Select gas particles within radial limits for the dust distribution
c
      gas_mass_inside = 0.
      ngas_inside = 0.
      DO i = 1, npart
         IF (iphase(i).EQ.0) THEN
            IF (igeom.EQ.'s') THEN
               r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
            ELSE
               r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
            ENDIF
            IF (SQRT(r2).GT.min_rplan.AND.SQRT(r2).LT.max_rplan) THEN
               gas_mass_inside = gas_mass_inside + xyzmh(4,i)
               ngas_inside = ngas_inside + 1
            ENDIF
         ENDIF
      END DO
      PRINT *,''
      PRINT *,'Total gas mass inside dust radii ',gas_mass_inside
      PRINT *,'Total gas particles inside dust radii ',ngas_inside

      PRINT *,''
      PRINT *,'There are ',nactive,' active particles ',
     &     'out of ',npart,' total'
c
c--Now create dust particles
c
c
c--If dust growth is to be used with implicit dust drag algorithm, set:
c
      IF (idim_grow.EQ.idim) THEN
 30      PRINT *,'What type of growth model do you want?'
         PRINT*, '1: Monodisperse, 2: Polydisperse'
         READ (iread,*) igrowthmodel
         IF (igrowthmodel.LT.1 .OR. igrowthmodel.GT.2) THEN
            PRINT *,'ERROR - maximum of 2 dust types for dust growth'
            GOTO 30
         ENDIF

         PRINT *,'Enter fragmentation velocity threshold (cm/s)',
     &        '(set to -1 to switch fragmentation off)'
         READ (iread,*) vfrag
         vfrag = vfrag*utime/udist

         PRINT *,'Enter minimum grain size with fragmentation (cm)?',
     &        '(set to -1 to avoid iminimum size)'
         READ (iread,*) rgrain_min
         rgrain_min = rgrain_min/udist
      ENDIF
c
c--This section creates standard dust particles (not background grains).
c     The dust particles have fixed sizes if growth is not active.
c     Or they grow if growth is active.
c
      IF (idim_grow.EQ.idim) THEN
c
c--The variable n_planetesimal_types is not used when running the
c     monodisperse or polydisperse dust growth models, it is just
c     set here to set up the 1 or 2 dust populations.
c
         n_planetesimal_types = igrowthmodel
         IF (n_planetesimal_types.GT.nplanetesimal_maxtype) THEN
            PRINT *,'ERROR - nplanetesimal_maxtype not large enough:',
     &           nplanetesimal_maxtype
            CALL quit(0)
         ENDIF
      ELSE
 40      PRINT *,'How many fixed-size dust species you want to create?'
         READ(iread,*) n_planetesimal_types
         IF (n_planetesimal_types.GT.nplanetesimal_maxtype) THEN
            PRINT *,'Max types is: ',nplanetesimal_maxtype
            GOTO 40
         ENDIF
      ENDIF
c
c--Create standard dust particles -- npart_gas sets part of array to
c     search for gas particles
c
      npart_gas = npart
      DO k = 1, n_planetesimal_types
 50      CONTINUE
         IF (idim_grow.EQ.idim) THEN
            IF (k.EQ.1) THEN
               PRINT *,'How many growing dust grains?'
            ELSEIF (k.EQ.2) THEN
               PRINT *,'How many fixed-size background dust grains ',
     &              'do you want?'
            ELSE
               PRINT *,'ERROR - max of 2 dust types for dust growth'
               CALL quit(0)
            ENDIF
         ELSE
            PRINT *,'How many particles for population',k,'?'
         ENDIF
         READ (iread,*) ndustpart
c
c--For MPI, split dust between processes -- ASSUMES that equal amounts
c     of dust are assigned to each MPI domain !
c
#ifdef MPI
         ndustpart = ndustpart/numproc
#endif
         IF (nactive + ndustpart.GT.idim) THEN
            PRINT *,'ERROR - number must be <= to ',idim-nactive
            GOTO 50
         ENDIF

         PRINT *,'Enter desired dust to gas ratio (normally < 1)'
         READ (iread,*) dust_gas_ratio
         dust_particle_mass=dust_gas_ratio*gas_mass_inside/ndustpart

         PRINT *,'Gas particle mass, Dust particle mass ',
     &        gas_mass_inside/ngas_inside,dust_particle_mass

         IF (idim_grow.EQ.idim) THEN
            IF (k.EQ.1) THEN
               PRINT *,'Enter initial size (in cm) of growing grains'
            ELSEIF (k.EQ.2) THEN
               PRINT *,' Enter size (in cm) of the background grains'
            ENDIF
         ELSE
            PRINT *,' Enter size (in cm) of dust population',k
         ENDIF
         READ (iread,*) r_planetesimals(k)
         r_planetesimals(k) = r_planetesimals(k)/udist

         ndust = 0
         DO j = 1, idim
c
c--Re-use any dead particles for new dust particles, otherwise add to
c     the end of the current particle list
c
            IF (iphase(j).EQ.-1 .OR. j.GT.npart) THEN
c
c--Distribute dust particles based on random sampling of gas 
c     particle distribution
c
 1237          i = MIN(INT(ran1(1)*npart_gas)+1,npart_gas)
               IF (iphase(i).NE.0) GOTO 1237

               IF (igeom.EQ.'s') THEN
                  r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2 + xyzmh(3,i)**2
               ELSE
                  r2 = xyzmh(1,i)**2 + xyzmh(2,i)**2
               ENDIF
               IF (SQRT(r2).LT.min_rplan.OR.SQRT(r2).GT.max_rplan)
     &              GOTO 1237

               ndust = ndust + 1

               xi = xyzmh(1,i) + 0.01*(2.*ran1(1)-1.)*xyzmh(5,i)
               yi = xyzmh(2,i) + 0.01*(2.*ran1(1)-1.)*xyzmh(5,i)
               zi = xyzmh(3,i) + 0.01*(2.*ran1(1)-1.)*xyzmh(5,i)

               xyzmh(1,j) = xi
               xyzmh(2,j) = yi
               xyzmh(3,j) = zi
               xyzmh(4,j) = dust_particle_mass
               xyzmh(5,j) = xyzmh(5,i)*
     &              (REAL(ngas_inside)/ndustpart)**(1./3.)

               vxyzu(1,j) = vxyzu(1,i)
               vxyzu(2,j) = vxyzu(2,i)
               vxyzu(3,j) = vxyzu(3,i)
               vxyzu(4,j) = vxyzu(4,i)

               iphase(j)  = 10 + k
               rho(j)     = dust_gas_ratio*rho(i)
               isteps(j)  = isteps(i)
c
c--If growth is active, the grain size is stored in the rgrain array.
c     Otherwise, grain sizes and densities are read from ifile.
c
               IF (idim_grow.EQ.idim) rgrain(j) = r_planetesimals(k)

               IF (ndust.EQ.ndustpart) THEN
                  npart = j
                  GOTO 100
               ENDIF
            ENDIF
         END DO

         PRINT *,'FAILED - only managed to create ',ndust,' particles ',
     &        'before running out of space'
         CALL quit(0)

 100     CONTINUE

         nactive = nactive + ndustpart
      END DO
c
c==END CHOICE
c
      ENDIF
c
c--Update nactive and the list of any remaining inactive particles
c
      nlstacc = 0
      nlistinactive = 0
      nkill = 0
      n1 = 0
      nactive = 0
      DO i=1, npart
         IF (iphase(i).EQ.-1) THEN
            nlistinactive = nlistinactive + 1
            IF (nlistinactive.GT.idim) THEN
               WRITE (iprint,*) 'ERROR step nlistinactive'
               CALL quit(1)
            ENDIF
            listinactive(nlistinactive) = i

         ELSE
            n1 = i
            nactive = nactive + 1
         ENDIF
      ENDDO
      PRINT *, 'Final Numbers ',npart, nactive, n1
c
c--Dump new file
c
c
c--Zero time (but not quite, so that the timesteps are still used when
c     starting the code).
c
      IF (gt.NE.0. .AND. ireset_time.EQ.'y') gt = 1.0E-20

      CALL wrinsph

      ifulldump = 0
      nfullstep = 1
      PRINT *,'writing dump file'
c
c--Transform into original frame of reference
c
      PRINT *,'Do you want to change back to original reference frame?'
      READ (iread,1001) iref

      IF (iref.EQ.'y' .OR. iref.EQ.'Y') CALL chanref(1)
c
c--Write dump file
c
#ifdef MPI
      IF (iproc.EQ.0) THEN
#endif
         CALL file
#ifdef MPI
      ENDIF
#endif

      CALL wdump_wrapper(idisk1)
c
c--End writing of full dump file
c-------------------------------
c

#ifdef MPI
      CLOSE(iread)
#endif

      PRINT 2000, file1
 2000 FORMAT ('file ', A10, 'has been created')
 
      CALL endrun

      END
