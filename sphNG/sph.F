      PROGRAM sph
c************************************************************
c************************************************************
c************************************************************
c******************                    **********************
c******************   S  P (M, R) H    **********************
c******************                    **********************
c************************************************************
c************************************************************
c************************************************************
c                                                           *
c  This is a three-dimensional hydro, MHD and radiation     *
c  hydro code based on the Smoothed Particle Hydrodynamics  *
c  method.                                                  *
c                                                           *
c************************************************************
c                                                           *
c  W. Benz        Los Alamos National Laboratory 31/08/86   *
c                 Harvard College Observatory    10/03/87   *
c                                                           *
c  I. A. Bonnell  University de Montreal         XX/XX/91   *
c                 University of Cambridge, UK    01/10/92   *
c                                                           *
c  M. R. Bate     University of Cambridge, UK    01/10/92   *
c                 University of Exeter, UK                  *
c                                                           *
c  D. J. Price    University of Cambridge, UK    XX/XX/03   *
c                 University of Exeter, UK       XX/XX/05   *
c                 Monash University, Australia   XX/XX/11   *
c                                                           *
c************************************************************
c                                                           *
c  Modifications by: I. A. Bonnell      UdeM    28/05/92    *
c                                                           *
c  (1) random number generator from Numerical Recipes: ran1 *
c  (2) cylindrical boundaries (ibound = 2): ghosp2          *
c  (3) gforsa and gforsn are parallelized                   *
c  (4) generator of axial density perturbations (ibound = 2)*
c                                                           *
c  Major Modifications by: M. R. Bate   UofCam  31/01/95    *
c                                                           *
c  (1) individual particle timesteps                        *
c  (2) inclusion of accreting, massive point masses (sink   *
c         particles) with boundaries, identified by iphase. *
c  (3) corrected ghost boundaries                           *
c  (4) constant pressure boundaries                         *
c                                                           *
c  Major Modifications by: M. R. Bate and R. Klessen        *
c                                         MPIA  02/02/96    *
c                                                           *
c  (1) switchable between BINARY TREE and GRAPE board       *
c         computation of gravity and neighbours lists       *
c                                                           *
c  Modifications by: D. J. Price       UofExeter 2003-2006  *
c                                                           *
c  (1) grad h terms, grad soft terms                        *
c  (2) MHD                                                  *
c                                                           *
c************************************************************
c  General Description                                      *
c                                                           *
c  This code allows three different type of tasks to be     *
c  done according to variable "job" . Variable job is read  *
c  from the command line.                                   *
c                                                           *
c  1) job=initial                                           *
c     -----------                                           *
c     this allows to compute new initial conditions. This   *
c     can be done by either starting from scratch or by     *
c     adding two existing situations together               *
c     (see subroutine newrun)                               *
c                                                           *
c  2) job=transfer                                          *
c     -------------                                         *
c     this allows to copy one or several dumps from one     *
c     file to another one. During this transfer  there is   *
c     a possibility of making changes in the various        *
c     quantities                                            *
c     (see subroutine transfd)                              *
c                                                           *
c  3) job=evolution                                         *
c     -------------                                         *
c     this allows to compute the actual evolution of the    *
c     system by specifying the various forces to be         *
c     included                                              *
c     (see subroutine evol)                                 *
c                                                           *
c************************************************************
c                                                           *
c  Running SPH                                              *
c  ===========                                              *
c                                                           *
c  The main option of the code, JOB, is read from the       *
c  command line. If JOB=INITIAL or TRANSFER, the code goes  *
c  into an interactive mode. Just answer the questions      *
c  after having read the informations below. If EVOLUTION   *
c  is chosen, the code is looking for a file called INSPH   *
c  where the various options have to be defined.            *
c                                                           *
c  files used and written by the code:                      *
c  -----------------------------------                      *
c                                                           *
c  1) If JOB=INITIAL and the SCRATCH option is chosen, then *
c     the code does not need any files but will produce     *
c     a binary file containing all the quantities needed    *
c     for the evolution part. Furthemore, the code produces *
c     an ascii file callled OUTSPH summarizing the code's   *
c     actions.                                              *
c  2) If JOB=INITIAL and the EXIST option is chosen, then   *
c     the code will require two files containing data to    *
c     combine them into a third file. An ascii file OUTSPH  *
c     is also produced.                                     *
c  3) If JOB=TRANSFER the code needs an existing binary file*
c     and upon completion it will create another one. Again *
c     an ascii file OUTSPH is created                       *
c  4) If JOB=EVOLUTION the code needs a binary file with    *
c     initial conditions produces by the INITIAL evolution  *
c     or the file resulting from a prior evolution run.     *
c     Again the file OUTSPH is produced and contains all    *
c     informations regarding the status of the evolution.   *
c     It is important to check the file periodically and    *
c     monitor the total energy and angular momentum for     *
c     example since both should be conserved (within 2% say)*
c                                                           *
c  Various important details                                *
c  -------------------------                                *
c  1) The code checks the size of the binary file. If it    *
c     grows beyond 10Mb it opens a new one. To keep track   *
c     of all these files, the code assumes that the name of *
c     the binary file is of the form XXXXXYY. XXXXX being   *
c     characters and YY numbers (typically 01 at the start).*
c     When a new file is created, its name will be XXXXXYY+1*
c     (see subroutine FILE)                                 *
c  2) The number of neighbors for each particles is  stored *
c     using the parameter ILIST. If dimensions too small    *
c     the code will issue a warning message but do the job. *
c     The dimensions should however be changed.             *
c  3) The code checks periodically for a file called        *
c     MESSAGE. This file permits to change some of the      *
c     options defined in the input deck during the code's   *
c     execution (see MESOP)                                 *
c  4) The debug option IDEBUG should be used with care since*
c     it produces a lont of output! ITRACE can be set to all*
c     to trace the code through the various subroutines in  *
c     case of mysterious errors!                            *
c  5) All forces can be turned off or on at will except for *
c     pressure gradients than should always be on.          *
c  6) The code in its present version assumes a perfect gas *
c     equation of state (see EOSPG). This can be complicated*
c     at will!.                                             *
c  7) The number of particles that can be used is determined*
c     by the parameter IDIM. The number of particles NPART  *
c     has to be less (or equal) to IDIM.                    *
c  8) Before anything can be done with the code, the proper *
c     units have to be defined (see UNIT). The number to    *
c     input are the total mass and a length scale. These    *
c     numbers are actually inside the code and not options  *
c     since for any particular choice of units, they should *
c     not be changed anymore! Once both numbers have been   *
c     coded in, the code determines all other units.        *
c  9) The sum of the mass of all particles should be equal 1*
c     unless, the total mass given in UNIT is not the one   *
c     represented by the particles.                         *
c                                                           *
c                                                           *
c-----------------------------------------------------------*
c                                                           *
c  job=evolution input deck (see OPTIONS)                   *
c  --------------------------------------                   *
c                                                           *
c  namerun  = name of the run                               *
c  file1    = name of file containing starting dump and on  *
c             which new dumps are added                     *
c  varsta   = entropy : use entropy as variable of state    *
c           = intener : use specific internal energy as     *
c                       variable of state                   *
c                                                           *
c  --------------                                           *
c  igrp     = 0       : no pressure gradients               *
c             1       : with pressure gradients             *
c  igphi    = 0       : no self-gravity                     *
c             1       : with self-gravity                   *
c  ifsvi    = 0       : no artificial viscosity             *
c             1       : with artificial viscosity           *
c  ifcor    = 0       : no coriolis forces                  *
c             1       : with coriolis forces                *
c  ichoc    = 0       : no heating due to artif. viscosity  *
c             1       : with heating due to artif. viscosity*
c  iener    = 0       : no total energy conservation        *
c             1       : with total energy conservation      *
c  damp     = value of general damping coefficient          *
c  ibound   = 0       : no boundaries                       *
c             1       : with reflective boundaries          *
c  iexf     = 0       : no external forces                  *
c             1       : allow for external forces           *
c  iexpan   = 0       : no general expansion                *
c             1       : general expansion superposed        *
c                                                           *
c  the external force, the boundaries and the general       *
c  damping are NOT considered in the energy equation.       *
c                                                           *
c  --------------                                           *
c  nstep    = frequency at which timesteps are written on   *
c             disk                                          *
c  tdump    = time between dumps (alternative to number     *
c             of steps between dumps)                       *
c  tnext    = time of next dump                             *
c  tol      = tolerance used by the Runge-Kutta integrator  *
c  ipos     = 9999    : use last dump as starting dump      *
c           = number  : use number as starting dump (all    *
c                       dumps behond number are destroyed!) *
c  tmax     = max time in minute allowed for the job        *
c  tstop    = time (computer unit) at which the program     *
c             will stop.                                    *
c                                                           *
c  only if iexpan= 1                                        *
c  -----------------                                        *
c  vexpan    =  expansion velocity of the frame of ref. in  *
c               cm/s.                                       *
c                                                           *
c  only if ifcor=1                                          *
c  ---------------                                          *
c  omeg0     =  angular velocity of the reference frame in  *
c               1/s                                         *
c                                                           *
c  only if ibound=1                                         *
c  ----------------                                         *
c  (all boundaries are used, if only one type is desired    *
c   it is necessary to set the other far away ... )         *
c                                                           *
c  rmax      =  max. radius for spherical boundary          *
c  xmin,xmax =  min and max x coord. for boundary           *
c  ymin,ymax =  min and max y coord. for boundary           *
c  zmin,zmax =  min and max z coord. for boundary           *
c                                                           *
c                                                           *
c************************************************************
c                                                           *
c  alphabetical list of common variables                    *
c  =====================================                    *
c                                                           *
c  variable  common block   description                     *
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*
c                                                           *
c  alpha        (shock) : coefficient for bulk viscosity    *
c  angto        (angm ) : total angular momentum            *
c  angx         (angm ) : x componant of total angular mom. *
c  angy         (angm ) : y componant of total angular mom. *
c  angz         (angm ) : z componant of total angular mom. *
c  beta         (shock) : coefficient for Neumann-Richtmeyer*
c                         viscosity                         *
c  cmx1         (out1 ) : x coord. of c.m. body 1           *
c  cmy1         (out1 ) : y coord. of c.m. body 1           *
c  cmz1         (out1 ) : z coord. of c.m. body 1           *
c  cmx2         (out2 ) : x coord. of c.m. body 2           *
c  cmy2         (out2 ) : y coord. of c.m. body 2           *
c  cmz2         (out2 ) : z coord. of c.m. body 2           *
c  cnormk       (kerne) : normalisation constant of kernel  *
c  damp         (dissi) : general damping coefficient       *
c  dq(idim)     (ener1) : heating due to shock dissipation  *
c  dt           (time ) : current timestep                  *
c  dvtable      (table) : spacing between table values      *
c  encal        (cgas ) : eos/energy calculation method flag*
c  escap        (fracg) : mass of escapors                  *
c  file1        (file ) : name of file number 1             *
c  file2        (file ) : name of file number 2             *
c  file3        (file ) : name of file number 3             *
c  fmass(itable)(table) : mass of particle as a funct. of r *
c  fmas1        (bodys) : fraction of total mass in body 1  *
c  fmas2        (bodys) : fraction of total mass in body 2  *
c  fnbtot       (units) : number of total mass              *
c  fpoten(itable)(table): potential energy as a funct. of r *
c  fx(idim)     (force) : x componant of total force        *
c  fy(idim)     (force) : y componant of total force        *
c  fz(idim)     (force) : z componant of total force        *
c  gamma        (cgas ) : value of gamma                    *
c  grwij(itable)(table) : gradient of kernel table          *
c  gt           (gtime) : current global time in code units *
c  h(idim)      (kerne) : smoothing length                  *
c  hma1         (out1 ) : max h for body 1                  *
c  hma2         (out2 ) : max h for body 2                  *
c  hmi1         (out1 ) : min h for body 1                  *
c  hmi2         (out2 ) : min h for body 2                  *
c  ibound       (typef) : flag to allow boundaries          *
c  icell(i,j,k) (nlist) : chaining mesh cells               *
c  ichoc        (typef) : flag to allow shock heating       *
c  idebug       (debug) : debug option                      *
c  idisk1       (logun) : logical unit number for disk file *
c  idisk2       (logun) : logical unit number for disk file *
c  idisk3       (logun) : logical unit number for disk file *
c  idist        (new  ) : type of distribution for setting  *
c                         the initial particles             *
c  idumy        (random): number to initialize generator    *
c  iener        (typef) : flag to allow total energy cons.  *
c  iexf         (typef) : flag to allow external forces     *
c  iexpan       (typef) : flag to allow expansion           *
c  ifcor        (typef) : flag to allow coriolis forces     *
c  ifsvi        (typef) : flag to allow arti. viscosity     *
c  igrp         (typef) : flag to allow pressure gradients  *
c  igphi        (typef) : flag to allow body forces         *
c  ipos         (recor) : record number wanted on disk      *
c  iprint       (logun) : logical unit number for print out *
c  irec         (recor) : record number on disk             *
c  iterm        (logun) : logical unit number for read      *
c  itrace       (debug) : flag to allow flow tracing        *
c  job          (actio) : main option that tells the code   *
c                         what option it has to run         *
c  n1           (bodys) : number of particles in body 1     *
c  n2           (bodys) : number of particles in body 2     *
c  namerun      (actio) : name of the run                   *
c  ncount       (tming) : number of timestep since last dump*
c                         was written on disk               *
c  np           (new  ) : number of particles at initial.   *
c  npart        (part ) : total number of particles         *
c  nstep        (tming) : number of  timestep between disk  *
c                         dumps                             *
c  tdump        (tming) : time between disk dumps           *
c                         (alternative to number of steps)  *
c  tnext        (tming) : time of next dump                 *
c  omeg0        (rotat) : angular velocity of frame of ref. *
c  partm        (new  ) : particles mass                    *
c  pdv(idim)    (ener1) : pdv work on particle              *
c  pmass(idim)  (carac) : mass of particles                 *
c  poten(idim)  (ener1) : potential energy of particle      *
c  pr(idim)     (eosq ) : pressure at particle locations    *
c  qphys(idco)  (tcond) : physical thermal conductivity     *
c  rho(idim)    (densi) : density at particle locations     *
c  rmax         (rbnd ) : max radius of boundary            *
c  romax1       (out1 ) : max density in body 1             *
c  romax2       (out2 ) : max density in body 2             *
c  romean1      (out1 ) : mean density body 1               *
c  romean2      (out2 ) : mean density body 2               *
c  v0           (space) : frac. of h of init. lattice cell  *
c  tgrav        (ener2) : total gravitational energy        *
c  thermal      (new  ) : either entropy or int. energy at  *
c                         initialisation                    *
c  thick        (new  ) : thickness ratio of initial plate  *
c  tkin         (ener2) : total kinetic energy              *
c  tmax         (tming) : maximum time allowed for the job  *
c  totm         (new  ) : total mass represented by part.   *
c  trot         (ener2) : total rotational energy           *
c  tstep        (tming) : cpu time needed for one timestep  *
c  tstop        (tming) : time at which simulation has to   *
c                         stop (code units)                 *
c  tterm        (ener2) : total thermal energy              *
c  u(idim)      (part ) : total specific internal energy    *
c  udens        (units) : unit of density                   *
c  udist        (units) : unit of length                    *
c  uergg        (units) : unit of energy per unit mass      *
c  uergcc       (units) : unit of pressure                  *
c  umass        (units) : unit of mass                      *
c  utime        (units) : unit of time                      *
c  varsta       (varet) : variable of state beside density  *
c  vcmx1        (out1 ) : x comp. of velocity of cm body 1  *
c  vcmy1        (out1 ) : y comp. of velocity of cm body 1  *
c  vcmz1        (out1 ) : z comp. of velocity of cm body 1  *
c  vcmx2        (out2 ) : x comp. of velocity of cm body 2  *
c  vcmy2        (out2 ) : y comp. of velocity of cm body 2  *
c  vcmz2        (out2 ) : z comp. of velocity of cm body 2  *
c  version      (infor) : version of the code               *
c  vexpan       (expan) : expansion velocity                *
c  vsmax        (eosq ) : maximum sound speed               *
c  vsound(idim) (eosq ) : sound speed at particle locations *
c  vvx1         (new  ) : x componant of velocity of center *
c                         of mass body 1                    *
c  vvx2         (new  ) : x componant of velocity of center *
c                         of mass body 2                    *
c  vvy1         (new  ) : y componant of velocity of center *
c                         of mass body 1                    *
c  vvy2         (new  ) : y componant of velocity of center *
c                         of mass body 2                    *
c  vvz1         (new  ) : z componant of velocity of center *
c                         of mass body 1                    *
c  vvz2         (new  ) : z componant of velocity of center *
c                         of mass body 2                    *
c  vx(idim)     (part ) : x componant of the part. velocity *
c  vy(idim)     (part ) : y componant of the part. velocity *
c  vz(idim)     (part ) : z componant of the part. velocity *
c  what         (new  ) : gives type of initialisation      *
c  wij(itable)  (table) : table of kernel values            *
c  x(idim)      (part ) : x coordinate of the particles     *
c  xmax         (rbnd ) : position of upper reflec. bound.  *
c  xmin         (rbnd ) : position of lower reflec. bound.  *
c  xx1          (new  ) : x coordinate of center of mass    *
c                         body 1                            *
c  xx2          (new  ) : x coordinate of center of mass    *
c                         body 2                            *
c  y(idim)      (part ) : y coordinate of the particles     *
c  ymax         (rbnd ) : position of upper reflec. bound.  *
c  ymin         (rbnd ) : position of lower reflec. bound.  *
c  yy1          (new  ) : y coordinate of center of mass    *
c                         body 1                            *
c  yy2          (new  ) : y coordinate of center of mass    *
c                         body 2                            *
c  z(idim)      (part ) : z coordinate of the particles     *
c  zmax         (rbnd ) : position of upper reflec. bound.  *
c  zmin         (rbnd ) : position of lower reflec. bound.  *
c  zz1          (new  ) : z coordinate of center of mass    *
c                         body 1                            *
c  zz2          (new  ) : z coordinate of center of mass    *
c                         body 2                            *
c                                                           *
c************************************************************
c                                                           *
c  timing and memory requirements                           *
c  ==============================                           *
c                                                           *
c  Both vary enormously with the various options chosen.    *
c  Running time scales like NlogN.                          *
c                                                           *
c************************************************************

      INCLUDE 'idim'
      INCLUDE 'igrape'

#ifdef MPIALL
      INCLUDE 'mpif.h'
      INCLUDE 'COMMONS/mpiall'
      INCLUDE 'COMMONS/mpidebug'
#endif
#ifdef MPI
      INCLUDE 'COMMONS/mpi'
#endif

      INCLUDE 'COMMONS/infor'
      INCLUDE 'COMMONS/actio'
      INCLUDE 'COMMONS/debug'
      INCLUDE 'COMMONS/savernd'
      INCLUDE 'COMMONS/rbnd'

#ifdef ZEN
      INTEGER (4)        :: clist(8)
      INTEGER :: mpi_ppn, omp_get_thread_num,omp_get_num_threads
#endif
#ifdef _OPENMP
      INTEGER*8 forcei_lock,revtree_lock
      COMMON /locks / forcei_lock(idim2/10+1),revtree_lock
#endif
#ifdef ZEN
      clist(1:8) = (/ 0, 1, 4, 5, 2, 3, 6, 7 /)
#endif
#ifdef MPIALL
      iloop = 0
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)
#ifdef ZEN
!$OMP PARALLEL
!$OMP.  DEFAULT (none)
!$OMP.  SHARED  (iproc, clist)
!$OMP.  PRIVATE (j,mpi_ppn)
      mpi_ppn=8/omp_get_num_threads() 
      j = MOD (iproc, mpi_ppn) + omp_get_thread_num() * mpi_ppn 
      if (iproc==0 .or. iproc == 1) print*,'iproc j mpi_ppn',iproc,j,
     .mpi_ppn
      CALL place_me (j, clist)
!$OMP END PARALLEL        
#endif
      IF (numproc.GT.nummaxproc) THEN
         WRITE (*,*) 'ERROR - numproc.GT.nummaxproc'
         CALL quit
      ENDIF
      DO i = 1, idim
         lblocklengths(i) = 1
      END DO
#else
#ifdef ZEN
!$OMP PARALLEL
!$OMP.  DEFAULT (none)
!$OMP.  SHARED  (clist)
!$OMP.  PRIVATE (j)
      j = omp_get_thread_num()
      CALL place_me (j, clist)
!$OMP END PARALLEL      
#endif
#endif
#ifdef MPI
      IF (idim_MPI.NE.idim) THEN
         PRINT *,'ERROR - idim_MPI.NE.idim'
         CALL quit
      ENDIF
#endif
c
c--Version - date is date last updated with other code version changes
c
      version = 'SPRMHD-MPI-1.0-MRB-DJP-16-05-06'
c
c   Derived from version = 'SPRMHD-1.1-MRB-DJP-12-01-06'
c     but with MPI/OpenMP hybrid capability.
c
c      version = 'SPRMHD-1.1-MRB-DJP-12-01-06'
c
c   Derived from version = 'SPRMHD-1.0-MRB-DJP-06-09-05'
c     but with MHD looking like it kind of, maybe, works.
c
c      version = 'SPRMHD-1.0-MRB-DJP-06-09-05'
c
c   Derived from version = 'SPH3DT-PTMAS-GRAPE-O4-1.0-MRB-17-04-98'
c     but including MHD, gradh terms and radiative transfer.
c
c      version = 'SPH3DTSPGO4P-OPT6-1.0-MRB-29-06-05'
c
c   Derived from version = 'SPH3DT-PTMAS-GRAPE-O4-1.0-MRB-17-04-98'
c     but has been rewritten to have arrays such as x,y,z,pmass,h combined
c     into xyzmh(5,*) etc for better cache useage.  Also, has leapfrog
c     integrator, new tree opening criteria, sink particle boundaries have
c     been removed and the new code has approximately half of the memory foot
c     print of the old code.
c
c   Same as version ='SPH3DT-PTMAS-GRAPE-3.0-MRB-20-02-97' (sph3DTSP_GRAPE3_O4)
c      but includes boundary types 90 and 91 and new grape_main_split.f
c      from sph3DTSP_GRAPE3.  
c      Similar to sph3DTSP_GRAPE3 but includes variable type of 
c      Kernel (3rd and 4th order kernels) from sph3DTSP_GRAPE3_O4.
c      ALSO, can be compiled with either single or double precision.
c
c   Similar to version 'SPH3DT-PTMAS-GRAPE-1.0-MRB-27-04-96'
c
c   Same as version 'SPH3D-TIMESTEPS-PTMAS-5.03-MRB-01-06-95' (infall8)
c      but with TREE/GRAPE switchable computation.  Similar to
c      version 'SPH3D-TIMESTEPS-PTMAS-5.01-MRB-01-02-95'
c
c   Same as version SPH3D-TIMESTEPS-PTMAS-4.01-MRB-30-01-94,
c     but tidied up - all unused stuff thrown away (old boundaries,
c     excessive memory use etc).
c
c   Same as version SPH3D-TIMESTEPS-PTMAS-3.01-MRB-16-05-94,
c     with boundaries as used by Bate,Bonnell,Price 1995
c     but compiled under FORTRAN 90.
c
c
c--Set debug and trace indicator
c
      idebug = 'off'
      itrace = 'off'
c      itrace = 'all'
c
c--Floating point number error checking (DCR C routine)
c
c      CALL floating
c
c--Read main main option for run
c
      IF (imhd.EQ.idim .AND. imhdxsph.NE.imhd) THEN
         WRITE (*,*) 'ERROR - imhd, imhdxsph ',imhd, imhdxsph
         CALL quit
      ELSEIF (XSPH .AND. imhdxsph.NE.idim) THEN
         WRITE (*,*) 'ERROR - XSPH, imhdxsph ',XSPH, imhdxsph
         CALL quit
      ENDIF

#ifdef _OPENMP
      DO i = 1, idim2/10+1
         CALL OMP_INIT_LOCK(forcei_lock(i))
      END DO
      CALL OMP_INIT_LOCK(revtree_lock)
#endif
      nlistinactive = 0

      CALL mainop
c
c--Initialize physical, mathematical and astronomical constants
c
      CALL constan
c
c--Define code's units
c
      CALL unit
c
c--Initialise starting date and time
c  (some versions of the zz file require this)
c
      CALL getdat(iday, imon, iyear)
c
c--Initialise random number generator arrays
c
      iyr = 0
      DO i = 1, NTAB
         iv(i) = 0
      END DO
c
c--Execute main option
c
c  1) define new initial conditions
c
      IF (job(1:7).EQ.'initial') CALL newrun
c
c  2) transfer dumps to another file
c
c      IF (job(1:8).EQ.'transfer') CALL transfd
c
c  3) integrate the system
c
      IF (job(1:9).EQ.'evolution') CALL evol
c
c  4) reduce output files
c
      IF (job(1:6).EQ.'reduce') CALL reduce
c
c  5) extract an output file
c
      IF (job(1:7).EQ.'extract') CALL extract
c
c  6) setup planet calculations
c
      IF (job(1:7).EQ.'planet') CALL newplanet

      CALL quit
      END
