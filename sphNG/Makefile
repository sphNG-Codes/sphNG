##############################################################
###                                                        ###
###           Make File for the sphNG Code                 ###
###      Original version by M. Bate and R. Klessen        ###
###                      March 1996                        ###
###                                                        ###
###           New version by D.Price July 2006,            ###
###        pinched from R.Kurosawa's TORUS Makefile        ###
###                                                        ###
###     This compiles the code using the BINARY TREE       ###
###                                                        ###
### (NOTE: As well as the different Makefile, `igrape'     ###
### must be changed when switching between GRAPE and TREE) ###
###                                                        ###
##############################################################

.KEEP_STATE:

KNOWN_SYSTEM=no
SHELL = /bin/bash

# can comment out the following lines and instead set
# the parameters as environment variables
DEBUG= no
DOUBLEPRECISION= yes
#openmp = yes
#mpi = yes
# endian can be "BIG", "LITTLE" or anything else which has no effect
ENDIAN= LITTLE

##############################################################
# preprocessor options 
# (uncomment as necessary)
#
FPPFLAGS=
FPP_SEP=
#
# This one is for periodic boundary conditions without using
# ghost particles
#
#FPPFLAGS+= -DPERIODIC_NO_GHOSTS 
#
# Add yours here...
#
##############################################################
ifeq ($(mpi), yes)
    mpiall=yes
endif
ifeq ($(mpicopy), yes)
    mpiall=yes
endif

ifeq ($(SYSTEM), ftn)
    FC= $(MPIDIR)ftn
    FFLAGS= -fPIC -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch
    OMPFLAG= -fopenmp
    DBLFLAG= -fdefault-real-8 -fdefault-double-8
    DEBUGFLAG= -g
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fconvert=big-endian
    ENDIANFLAGLITTLE= -fconvert=little-endian
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), ompi)
    FC= $(MPIDIR)mpif90
    FFLAGS= -DMPI -O3 -march=pentium4 -mtune=pentium4 -malign-double -Wno=155
    OMPFLAG= -openmp
    DBLFLAG= -r8
    DEBUGFLAG= -g -ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), intelmac)
    ifeq ($(PROFILE), yes)
      FFLAGS= -pg
    else
      FFLAGS=
    endif
    FC= g95
    DBLFLAG= -r8
    DEBUGFLAG= -g -ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), zen)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
#       FC= ifort -lmpi
    else
       FC= ifort
    endif
    FFLAGS= -O3 -xSSSE3 -mcmodel=medium -shared-intel  \
            -warn uninitialized -warn truncated_source \
            -warn interfaces -nogen-interfaces
# -check
#To use above MPI check need to:
#source /sw/sdev/itac/7.2.0.011/bin/itacvars.csh
#
    OMPFLAG= -openmp #-DZEN
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG= -ipo
#    ZEN_PIN = pin_to_cpu.o
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), skylake)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -xCORE-AVX512 -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check_mpi
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), Cosma8)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -march=core-avx2 -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces -Wl,--no-relax
# -check_mpi
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), DIAL3_Intel)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -mavx2 -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check_mpi
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), DIAL3_AMD)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpif90
    else
       FC= mpifort
    endif
    FFLAGS= -O3 -mcmodel=medium
#            -warn uninitialized -warn truncated_source\
#            -warn interfaces -nogen-interfaces
# -check_mpi
    OMPFLAG= -fopenmp
    DBLFLAG= -fdefault-real-8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), kennedy)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), KNLCam)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -xMIC-AVX512 -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check_mpi
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), isca)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -axAVX -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check_mpi
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), slave)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
    else
       FC= ifort
    endif
    FFLAGS= -O3 -axAVX -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check_mpi
    OMPFLAG= -qopenmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    LDFLAGS= -L/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/

    IPOFLAG=
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), dirac2)
    ifeq ($(mpiall), yes)
       FC= $(MPIDIR)mpiifort
#       FC= ifort -lmpi
    else
       FC= ifort
    endif
    FFLAGS= -O3 -axAVX -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check
#To use above MPI check need to:
#source /sw/sdev/itac/7.2.0.011/bin/itacvars.csh
#
    OMPFLAG= -qopenmp #-DZEN
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG= -ipo
#    ZEN_PIN = pin_to_cpu.o
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), ukaff1a)
    ifeq ($(mpiall), yes)
        FC= $(MPIDIR)mpif77 -qnosave
    else
        FC= xlf_r -qnosave
    endif
#    FFLAGS= -O3 -q64 -qextname -qsmp=noauto 
    FFLAGS= -O3 -q64
    OMPFLAG= -qsmp=noauto
    DBLFLAG= -qrealsize=8
#    DEBUGFLAG= -qflttrap=enable:invalid:zerodivide -g -C -qsigtrap -qfloat=nans 
#    DEBUGFLAG= -C
    KNOWN_SYSTEM=yes
    FPP_SEP=','
    FPP_PREFIX=-d -WF,
#    ZZFILE= zzsun_f77_double1.f zzsun_f77_double2.f
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), apac)
    ifeq ($(mpiall), yes)
       LDFLAGS+= -lmpi
    endif
    FC=ifort
    FFLAGS= -O3 -i-dynamic -WB \
            -w90 -w95 -cm -warn uninitialized -warn unused -warn truncated_source\
            -warn interfaces -gen-interfaces
    OMPFLAG= -openmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fpstkchk
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), msg)
    ifeq ($(mpiall), yes)
       LDFLAGS+= -lmpi
    endif
    FC=ifort
    FFLAGS= -O3 -WB -mcmodel=medium -shared-intel \
            -w90 -w95 -cm -warn uninitialized -warn unused -warn truncated_source\
            -warn interfaces -gen-interfaces
    OMPFLAG= -openmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fpstkchk
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), sunf95)
    FC= sunf95
    FFLAGS= -fast -w4 -errtags -erroff=COMMENT_1582,COMMENT_1744
    DBLFLAG= -xtypemap=real:64,double:64
    OMPFLAG= -xopenmp
    DEBUGFLAG= -g -C -ftrap=%all
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -xfilebyteorder=big16:%all ##:45
    ENDIANFLAGLITTLE= -xfilebyteorder=little16:%all
    ZZFILE= zzsun_f77_double1.f zzsun_f77_double2.f
endif

ifeq ($(SYSTEM), g95)
    FC= g95
    FFLAGS= -O3 -Wall -Wno=112,140 -ffast-math 
    DBLFLAG= -r8
    DEBUGFLAG= -g -ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
# or use setenv G95_ENDIAN='BIG' or 'LITTLE' at runtime
    TORUSLIB= -L/h/mbate/New/sphNG_TORUS/torus -ltorus -L/star/lib `ndf_link` -L/usr/lib/gcc/i586-mandrake-linux-gnu/3.4.1 -lg2c
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), myg95)
    FC= myg95
    FFLAGS= -O3 -Wall -Wno=112,140 -Wextra -ffast-math 
    DBLFLAG= -r8
    DEBUGFLAG= -ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), gfortran)
    FC= gfortran
    ifeq ($(mpi), yes)
       FC= $(MPIDIR)mpif90
    endif
# -I. is needed so files in RT/ can find the include files... (ifort seems to do this implicitly)
    FFLAGS= -m64 -mcmodel=medium -O3 -I. -Wall -Wno-conversion -Wno-unused-dummy-argument -Wno-maybe-uninitialized -Warray-temporaries
    DBLFLAG= -fdefault-real-8 -fdefault-double-8
    OMPFLAG= -fopenmp
    DEBUGFLAG= -g # -fcheck=all
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fconvert=big-endian
    ENDIANFLAGLITTLE= -fconvert=little-endian
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), ifc)
    FC= ifc
    FFLAGS= -O3 -WB -Vaxlib -w90 -w95 -cm
    OMPFLAG= -openmp
    DBLFLAG= -r8
    DEBUGFLAG= -C -g
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), ifort)
    FC= ifort
    FFLAGS= -O3 -xHost -mcmodel=medium -shared-intel  \
            -warn uninitialized -warn truncated_source \
            -warn interfaces -nogen-interfaces 
    OMPFLAG= -openmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    IPOFLAG= -ipo
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

# Set other optional flags depending on settings

ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
endif

ifeq ($(ENDIAN), BIG)
    FFLAGS += ${ENDIANFLAGBIG}
    FPPFLAGS+=${FPP_SEP}-DENDIANBIG
endif

ifeq ($(ENDIAN), LITTLE)
    FFLAGS += ${ENDIANFLAGLITTLE}
    FPPFLAGS+=${FPP_SEP}-DENDIANLITTLE
endif

ifeq ($(openmp), yes)
    FFLAGS += ${OMPFLAG}
endif

ifneq ($(ipo), no)
# No ipo if debugging is on
    ifneq ($(DEBUG), yes)
        FFLAGS += ${IPOFLAG}
    endif
endif

ifeq ($(mpi), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIALL
    FPPFLAGS+=${FPP_SEP}-DMPI
endif

ifeq ($(mpicopy), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIALL
    FPPFLAGS+=${FPP_SEP}-DMPICOPY
endif

ifeq ($(mpimod), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIMOD
endif

ifeq ($(mpidebug), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIDEBUG
endif

ifeq ($(mpidebugt), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIDEBUGT
endif

ifeq ($(mpidebugs), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIDEBUGS
endif

ifeq ($(mpidebugrt), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIDEBUGRT
endif

ifeq ($(dgrid), yes)
    FPPFLAGS+=${FPP_SEP}-DDGRID
endif

ifeq ($(nonideal), yes)
    FPPFLAGS+=${FPP_SEP}-DNONIDEAL
endif

ifeq ($(krome), yes)
   FFLAGS+=-I${HOME}/Codes/krome2/build_complexity
#   KLIB=-L${HOME}/Codes/krome/build
   FPPFLAGS+=${FPP_SEP}-DUSEKROME
   KROMEOBJ=libkrome.so utils/krome_sph_subrtns.F
endif

ifeq ($(kromempi),yes)
   FC=${MPIDIR}mpiifort
   FFLAGS+=-I${HOME}/Codes/krome2/build_complexity
#   KLIB=-L${HOME}/Codes/krome/build
   FPPFLAGS+=${FPP_SEP}-DUSEKROME -DKROMEMPI
   KROMEOBJ=libkrome.so utils/krome_sph.F
endif


# do double precision flag last (append only to FFLAGS)

ZZFFLAGS := ${FFLAGS}
ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAG}
endif

# add preprocessor flags (and prefix required for some compilers)
# the make implicit rule uses CPPFLAGS for both .c and .F files...
CPPFLAGS=${FPP_PREFIX}$(FPPFLAGS)

FFLAGSKEEP := ${FFLAGS}

KERNEL= ktable_QuinticM6.f
OBJKERNEL= $(KERNEL:.f=.o)

# these are the sources common to all compilations

SOURCES= sph.F accrete_P.F addump.f angmom.f \
	bonnorebert.F boundry.f build_kernel_column.f build_kernel_n2.f \
	buildsightlines.f build_stellar_radiation_table.f \
	cartdis.f cartmas.f cartpres.f cartvel.f chanref.f chekopt.f \
	condense.f column_los.F \
	coolcurve.f coolinmo.F cool_func.F \
	coriol.f cyldis.f cylmas.f discboundary.f \
	divBclean_subcycle.F divBsubcycle.F \
	dust_implicit.F dust_drag.F dust_growth.F dust_norm_rhoG.F \
	dragtable.F endrun.f energ.f energ_cooling.f eospg.F error.f \
	evol.f externf.f externBfield.f \
	extract.f file.f gapfactor.f getneigh.f \
	getcv.f get1overmu.f getkappa.f \
	gforsa_P.f gforsn_P.f gforspt_P.f \
	ghostp.f ghostp1.f ghostp2.f ghostp3.F ghostp11.f ghostp100.f \
	ghostp102.f ghostp104.f \
	gptfromgas.f gpti.f header.F hchem.F hillmass.f \
	homexp.f hy09dustevol.f hy09dustfunc.f hy09dustset.f \
	indexx.f indexx2.f indexx3.f indexxi2.f inform.f \
	inopts.f insulate_TREE_P.F integs.F $(KERNEL) \
	kill_particles.F labrun.F \
	makedust.F makeplot.f mesop.f mhdparams_mhd.f \
	modif.f move_particle.F mtree_P.f newrun.f \
	one_fluid_dust_utils.f one_fluid_dustvar_to_dustfrac.f options.f \
	optimise_column_density.f phoenix.f phoenix2.f \
	planetpotential.f planet_setup.f planetesimals.f planck.f \
	potential.f preset.F prout.F polytrope.F Qv.f \
	ran1.f rdump.F rdump_wrapper.F read_dragtable.f reduce.f \
	rho_grain.f \
	maketables/specificheat.f resistivity.f resistivity_implicit.F save.f \
	scaling.f \
	semimajoraxis.f setBfield_mhd.f setpart.F setup_torus.f \
	smoothd_P.f sphericalposition.f \
	sphdis.f sphmas.f spline.F stellar_radiation.f stoppingtime.f \
	thermal_velocity.f thermeq.f toten.f treef_P.F \
	unifdis.F unit.f uset.f \
	wdump.F wdump_wrapper.F wrinsph.F write_dragtable.f zeusread.f

# these are the sources required for the gradh code to compile

SRCGRADH= derivi_P_gradh_MPI.F forcei_P_gradh_MPI.F  \
        densityiterate_gradh_MPI.F divv_gradhsoft.F \
	hcalc_gradh.f dissipative_source.F
	#divBclean_mhd.f divBdiffuse.f divBiterate.f \
	#directsum_poisson_vec_mhd_clean.f 

# these are the extra routines required for the MPI code to compile
ifeq ($(mpiall), yes)
SRCMPI=balance.F classify_MPI.F
endif

# these are the extra routines required for the radiative transfer code to compile

SRCRT= ass.f RT/newgsimpl_MPI.f RT/boundaryparticle.f quartic_gs1t.f \
	RT/getblock.f RT/dust_temperature.f column_density.F
SRCNORT =
RTFPPFLAGS=${FPP_SEP}-DRT

# these are the extra routines required for implicit solve of one-fluid dust

SRCID1= dust1_implicit.F dust1_diffusion_implicit.F quartic_dust1.f \
	quadratic_dust1.f ass_dust1.F

# these are the extra routines required for Monte-Carlo radiative transfer

SRCTORUS= montecarloRT.F

# these are the extra routines required for non-ideal MHD
ifeq ($(nonideal), yes)
SRCNONIDEAL= calc_curlB.F
SRCNI = nicil.F90
SRCNS = nicil_sup.F90
endif

# timestepping routine

SRCRK = step_P.F revtree_P.F
SRCLF = step_P_lf_fast2_gradh.F revtree_P_fast3.f timestep.F

# these are the extra routines required for the stellar evolution code

SRCBAR = RT/Baraffe/1.f RT/Baraffe/gasrad.f RT/Baraffe/sr_all.f \
	RT/Baraffe/2pz.f RT/Baraffe/sr_r.f \
	RT/Baraffe/3libeb.f RT/Baraffe/interNEG.f RT/Baraffe/sr_rx.f \
	RT/Baraffe/blodat0.f RT/Baraffe/kapz2.f RT/Baraffe/sr_x.f \
	RT/Baraffe/new4.f RT/Baraffe/evol.f RT/Baraffe/screenchab99.f
TOBJBAR = $(SRCBAR:.f=.o)
OBJBAR = $(TOBJBAR:.F=.o)
RTBARFPPFLAGS=${RTFPPFLAGS} ${FPP_SEP}-DRTBAR

#FFLAGS := -save ${ZZFFLAGS}

# these routines compiled single precision and with -save

#makebaraffe: $(OBJBAR)
#	$(FC) $(FFLAGS) -o RT/Baraffe/testevol $(OBJBAR)


# construct object files

TOBJECTS = $(SRCMPI:.F=.o)  $(SRCNI:.F90=.o) $(SRCNS:.F90=.o) $(SOURCES:.f=.o) $(ZZFILE:.f=.o) $(SRCNONIDEAL:.F=.o) $(SRCID1:.f=.o)
TOBJGRADHLF = $(TOBJECTS) $(SRCGRADH:.f=.o) $(SRCLF:.f=.o) $(SRCNORT:.f=.o) $(ZEN_PIN)
TOBJGRADHRK = $(TOBJECTS) $(SRCGRADH:.f=.o) $(SRCRK:.f=.o) $(SRCNORT:.f=.o) $(ZEN_PIN)
TOBJGRADHLFRT = $(TOBJECTS) $(SRCGRADH:.f=.o) $(SRCLF:.f=.o) $(SRCRT:.f=.o) $(ZEN_PIN)
TOBJGRADHRKRT = $(TOBJECTS) $(SRCGRADH:.f=.o) $(SRCRK:.f=.o) $(SRCRT:.f=.o) $(ZEN_PIN)
TOBJGRADHRKRTBAR = $(TOBJECTS) $(SRCGRADH:.f=.o) $(SRCRK:.f=.o) $(SRCRT:.f=.o) $(SRCBAR:.f=.o) $(ZEN_PIN)
TOBJTORUS = $(TOBJECTS) $(SRCGRADH:.f=.o) $(SRCRK:.f=.o) $(SRCNORT:.f=.o) $(SRCTORUS:.f=.o) $(ZEN_PIN)

OBJECTS = $(TOBJECTS:.F=.o) 
OBJGRADHLF = $(TOBJGRADHLF:.F=.o) 
OBJGRADHRK = $(TOBJGRADHRK:.F=.o) 
OBJGRADHLFRT = $(TOBJGRADHLFRT:.F=.o)
OBJGRADHRKRT = $(TOBJGRADHRKRT:.F=.o)
OBJGRADHRKRTBAR = $(TOBJGRADHRKRTBAR:.F=.o) 
OBJTORUS = $(TOBJTORUS:.F=.o) 

# to ensure that the .mod file is made
ifeq ($(nonideal), yes)
.SUFFIXES: .F90 
%.o : %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c -o $(SRCNI:.F90=.o) $(SRCNI)
	$(FC) $(FFLAGS) $(CPPFLAGS) -c -o $(SRCNS:.F90=.o) $(SRCNS)
endif

# default action for "make"
default: gradhrk

# specific rules to make MPI object files

#mpif.h :
#	@echo " "
#	@echo "Including mpif.h in the following file............."
#	@echo " "

#%.o: %.f
#	$(FC) $(FFLAGS) -c $<

#sph.o : sph.f.raw
#	ifeq ($(MK_MPI_SRC),yes)
#	   sed 's/C$MPI/     /g' sph.f.raw   > sph.f
#	else
#	   sed '/C$MPI/d' sph.f.raw   > sph.f
#	endif
#	$(FC) $(FFLAGS) -c sph.f

# Recursive make to build Baraffe stellar evo code
baraffe:
	cd RT/Baraffe/ && $(MAKE)

pin_to_cpu.o: pin_to_cpu.c
	icc -c pin_to_cpu.c -openmp

torus: checksystem checkparams sphtorus

gradh: checksystem checkparams sph_tree_lf_gradh 

gradhlf: checksystem checkparams sph_tree_lf_gradh 

gradhrk: checksystem checkparams sph_tree_rk_gradh 

gradhlfrt : FPPFLAGS+=${RTFPPFLAGS}
gradhlfrt: checksystem checkparams sph_tree_lf_gradh_RT 

gradhrkrt : FPPFLAGS+=${RTFPPFLAGS}
gradhrkrt: checksystem checkparams sph_tree_rk_gradh_RT

gradhrkrtbar : FPPFLAGS+=${RTBARFPPFLAGS}
gradhrkrtbar: checksystem checkparams sph_tree_rk_gradh_RT_BAR

all: gradhlf gradhrk

sph_tree_lf_gradh: $(OBJGRADHLF)
	@echo ""
	@echo "This is the gradh, leapfrog version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHLF) $(LDFLAGS)

sph_tree_rk_gradh: $(OBJGRADHRK)
	@echo ""
	@echo "This is the gradh, Runge Kutta version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHRK) $(LDFLAGS)

sph_tree_lf_gradh_RT: $(OBJGRADHLFRT)
	@echo ""
	@echo "This is the gradh, leapfrog version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHLFRT) $(LDFLAGS)

sph_tree_rk_gradh_RT: $(OBJGRADHRKRT)
	@echo ""
	@echo "This is the gradh, Runge Kutta, radiative transfer version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHRKRT) $(LDFLAGS)

sph_tree_rk_gradh_RT_BAR: $(OBJGRADHRKRT) baraffe
	@echo ""
	@echo "This is the gradh, Runge Kutta, radiative transfer version"
	@echo "It includes the Baraffe stellar evolution code"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHRKRT) $(OBJBAR) $(LDFLAGS)

sphtorus: $(OBJTORUS)
	@echo ""
	@echo "This is the TORUS code (Runge Kutta by default)"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJTORUS) -L${TORUS_LIB} $(LDFLAGS)

util: binary convert convertfs convertof ptmass massvsangmom

binary: checksystem checkparams utils/binary.o
	$(FC) $(FFLAGS) -o $@ utils/binary.o

convert: checksystem checkparams utils/convert.o utils/readdump_sphNG.o utils/writedump_ascii.o
	$(FC) $(FFLAGS) -o $@ utils/convert.o utils/writedump_ascii.o utils/readdump_sphNG.o

# convert from scalar dump file to MPI dump file (i.e. multiple blocks)

convertmpi: checksystem checkparams utils/convert_to_mpi.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o 
	$(FC) $(FFLAGS) -o $@ utils/convert_to_mpi.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o 

# convert from MPI dump with type "b" domain decompostion and N domains
#    to N/2 domains (e.g. from 128 domains to 64 domains).

mpi_reduce: checksystem checkparams utils/mpi_reduce.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o balance.o
	$(FC) $(FFLAGS) -o $@ utils/mpi_reduce.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o balance.o

# convert from MPI dump with type "d" domain decomposition to type "b"

mpi_balance: checksystem checkparams utils/mpi_balance.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o balance.o
	$(FC) $(FFLAGS) -o $@ utils/mpi_balance.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o balance.o

# convert from any dump to 'b' with any number of ranks (compile in serial)

mpi_modify: checksystem checkparams utils/mpi_modify.o wdump.o rdump.o indexx3.o zzibm.o error.o file.o
	$(FC) $(FFLAGS) -o $@ utils/mpi_modify.o wdump.o rdump.o indexx3.o zzibm.o error.o file.o

# convert from full to small dump file

convertfs: checksystem checkparams utils/convert_full_small.o wdump.o rdump.o error.o balance.o
	$(FC) $(FFLAGS) -o $@ utils/convert_full_small.o wdump.o rdump.o error.o balance.o

# convert from old dump file to new full dump file

convertof: checksystem checkparams utils/convert_old_full.o wdump.o error.o
	$(FC) $(FFLAGS) -o $@ utils/convert_old_full.o wdump.o error.o

# 3D movie processing

movie3DRG: checksystem checkparams utils/movie3DRG/imcol3DNG_zoom_DT.o rdump.o error.o utils/movie3DRG/grid_P_DT.o utils/movie3DRG/buildcol.o utils/movie3DRG/origin.o
	$(FC) $(FFLAGS) -o $@ utils/movie3DRG/imcol3DNG_zoom_DT.o rdump.o error.o utils/movie3DRG/grid_P_DT.o utils/movie3DRG/buildcol.o utils/movie3DRG/origin.o

movie3DISM: checksystem checkparams utils/movie3DISM/imcol3DNG_zoom_ISM.o rdump.o error.o utils/movie3DISM/grid_P_ISM.o utils/movie3DISM/buildcol.o utils/movie3DISM/origin.o
	$(FC) $(FFLAGS) -o $@ utils/movie3DISM/imcol3DNG_zoom_ISM.o rdump.o error.o utils/movie3DISM/grid_P_ISM.o utils/movie3DISM/buildcol.o utils/movie3DISM/origin.o

# 2D movie processing

imcol3DNG: checksystem checkparams utils/movie_col3DALL/imcol3DNG_follow_iunique.o rdump.o error.o utils/movie_col3DALL/buildcol.o utils/movie_col3DALL/grid_follow.o utils/movie_col3DALL/origin.o utils/movie_col3DALL/value.o
	$(FC) $(FFLAGS) -o $@ utils/movie_col3DALL/imcol3DNG_follow_iunique.o rdump.o error.o utils/movie_col3DALL/buildcol.o utils/movie_col3DALL/grid_follow.o utils/movie_col3DALL/origin.o utils/movie_col3DALL/value.o

# Generate T vs rho chemistry pixel map

pixmap: checksystem checkparams utils/image_pixmap.o rdump.o error.o utils/movie_col3DALL/origin.o unit.o
	$(FC) $(FFLAGS) -o $@ utils/image_pixmap.o rdump.o error.o utils/movie_col3DALL/origin.o unit.o

pixmap_dust: checksystem checkparams utils/imagep3DNG_pixmap/imagep3DNG_pixmap_dust.o rdump.o error.o utils/movie_col3DALL/origin.o unit.o
	$(FC) $(FFLAGS) -o $@ utils/imagep3DNG_pixmap/imagep3DNG_pixmap_dust.o rdump.o error.o utils/movie_col3DALL/origin.o unit.o

pixmap_gas: checksystem checkparams utils/imagep3DNG_pixmap/imagep3DNG_pixmap_gas.o rdump.o error.o utils/movie_col3DALL/origin.o unit.o
	$(FC) $(FFLAGS) -o $@ utils/imagep3DNG_pixmap/imagep3DNG_pixmap_gas.o rdump.o error.o utils/movie_col3DALL/origin.o unit.o

# Disc parameter extraction

disc: checksystem checkparams utils/disc_parameters.o rdump.o wdump.o error.o indexx2.o
	$(FC) $(FFLAGS) -o $@ utils/disc_parameters.o rdump.o wdump.o error.o indexx2.o

# Single disc parameter extraction

disc1: checksystem checkparams utils/discmass.o rdump.o error.o indexx2.o $(ZZFILE:.f=.o)
	$(FC) $(FFLAGS) -o $@ utils/discmass.o rdump.o error.o indexx2.o $(ZZFILE:.f=.o)

# Ptmass evaluate properties of multiple systems

multiples: checksystem checkparams utils/ptmass_multiples.o indexx.o
	$(FC) $(FFLAGS) -o $@ utils/ptmass_multiples.o indexx.o

# Ptmass join multiple P-files together into one P-file

join: checksystem checkparams utils/ptmass_join_Pfiles.o
	$(FC) $(FFLAGS) -o $@ utils/ptmass_join_Pfiles.o

# Take radial cross sections through dump files

profile: checksystem checkparams $(KROMEOBJ) utils/profile.o rdump.o error.o indexx.o utils/movie3DRG/origin.o getcv.o hy09dustset.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $(KROMEOBJ) utils/profile.o rdump.o error.o indexx.o utils/movie3DRG/origin.o getcv.o hy09dustset.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o)

ptmass: utils/ptmassmerge.o
	$(FC) $(FFLAGS) -o $@ utils/ptmassmerge.o

massvsangmom: checksystem checkparams utils/massvsangmom.o utils/readdump_sphNG.o indexx2.o idim
	$(FC) $(FFLAGS) -o $@ utils/massvsangmom.o utils/readdump_sphNG.o indexx2.o

masswithtime: checksystem checkparams utils/masswithtime.o utils/readdump_sphNG.o idim
	$(FC) $(FFLAGS) -o $@ utils/masswithtime.o utils/readdump_sphNG.o

printkernel: checksystem checkparams $(OBJKERNEL) utils/printkernel.o
	$(FC) $(FFLAGS) -o $@ $(OBJKERNEL) utils/printkernel.o

# Additional analysis routines by JHW; see header of each file for description & insrictions
jhw_print_sinkorigin: checksystem checkparams utils/jhw_print_sinkorigin.o error.o rdump.o
	$(FC) $(FFLAGS) -o $@ utils/jhw_print_sinkorigin.o error.o rdump.o

jhw_moddump_CoM: checksystem checkparams utils/jhw_moddump_CoM.o utils/jhw_commons.F indexx.o wdump.o rdump.o error.o file.o options.o chekopt.o
	$(FC) $(FFLAGS) -o $@ utils/jhw_moddump_CoM.o utils/jhw_commons.F indexx.o wdump.o rdump.o error.o file.o options.o chekopt.o

jhw_calc_Mach: checksystem checkparams utils/jhw_calc_Mach.o utils/jhw_commons.o indexx.o error.o rdump.o
	$(FC) $(FFLAGS) -o $@ utils/jhw_calc_Mach.o utils/jhw_commons.o indexx.o error.o rdump.o

jhw_calc_etaNden: checksystem checkparams utils/jhw_calc_etaNden.o utils/jhw_commons.o indexx.o hy09dustset.o rdump.o error.o options.o chekopt.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o)
	$(FC) $(FFLAGS) -o $@ utils/jhw_calc_etaNden.o utils/jhw_commons.o indexx.o hy09dustset.o rdump.o error.o options.o chekopt.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o)

jhw_create_stellarICs: checksystem checkparams utils/jhw_create_stellarICs.o utils/jhw_commons.o getcv.o ran1.o getkappa.o indexx.o wdump.o rdump.o error.o file.o options.o chekopt.o
	$(FC) $(FFLAGS) -o $@ utils/jhw_create_stellarICs.o utils/jhw_commons.o getcv.o ran1.o getkappa.o indexx.o wdump.o rdump.o error.o file.o options.o chekopt.o

jhw_calc_stellarICs: checksystem checkparams utils/jhw_calc_stellarICs.o utils/jhw_commons.o indexx.o error.o rdump.o
	$(FC) $(FFLAGS) -o $@ utils/jhw_calc_stellarICs.o utils/jhw_commons.o indexx.o error.o rdump.o

# set up particle distribution without using entire SPH executable

setup_alone: checksystem checkparams utils/setup_alone.o wdump.o error.o file.o chekopt.o setpart.o $(OBJKERNEL) unifdis.o condense.o cartmas.o cylmas.o sphmas.o cartdis.o cyldis.o sphdis.o cartpres.o cylpres.o sphpres.o eospg.o setBfield_mhd.o inopts.o cartvel.o cylvel.o sphvel.o gapfactor.o uset.o wrinsph.o ran1.o unit.o $(ZZFILE:.f=.o)
	$(FC) $(FFLAGS) -o $@ utils/setup_alone.o wdump.o error.o file.o chekopt.o setpart.o $(OBJKERNEL) unifdis.o condense.o cartmas.o cylmas.o sphmas.o cartdis.o cyldis.o sphdis.o cartpres.o cylpres.o sphpres.o eospg.o setBfield_mhd.o inopts.o cartvel.o cylvel.o sphvel.o gapfactor.o uset.o wrinsph.o ran1.o unit.o $(ZZFILE:.f=.o)

getdeltat: checksystem checkparams utils/getdeltat.o utils/readdump_sphNG.o idim
	$(FC) $(FFLAGS) -o $@ utils/getdeltat.o utils/readdump_sphNG.o

# this routine needs to be compiled with single precision flags only

zzsun_f77_double2.o: zzsun_f77_double2.f
	$(FC) -c $(ZZFFLAGS) zzsun_f77_double2.f -o zzsun_f77_double2.o

krome_sph: checksystem checkparams libkrome.so utils/krome_sph.o rdump.o rdump_wrapper.o wdump.o error.o indexx.o chanref.o scaling.o unit.o $(ZZFILE:.f=.o)
	@echo "cppflags"
	@echo  $(CPPFLAGS)
	@echo "end flags"
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ libkrome.so utils/krome_sph.o rdump.o rdump_wrapper.o wdump.o error.o indexx.o chanref.o scaling.o unit.o $(ZZFILE:.f=.o)

krome_init:
	$(FC) $(FFLAGS) $(CPPFLAGS) libkrome.so -o krome_init.out utils/krome_init.f90

krome_rw: checksystem checkparams libkrome.so utils/krome_rw.o rdump.o wdump.o error.o indexx.o chanref.o scaling.o unit.o $(ZZFILE:.f=.o)
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ libkrome.so utils/krome_rw.o rdump.o wdump.o error.o indexx.o chanref.o  scaling.o unit.o $(ZZFILE:.f=.o)


checksystem:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling sphNG for $(SYSTEM) system..........."
	@echo ""
	@echo "calling creategitversion"
	@./creategitversion
	@echo "...done"
   else
	@echo ""
	@echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=>set the environment variable SYSTEM to one listed "
	@echo "  in the Makefile and try again"
	@echo ""
	quit
   endif

splitdump: checksystem checkparams $(KROMEOBJ) utils/splitdump.F90 rdump.o wdump.o error.o indexx.o options.o chekopt.o labrun.o zzibm.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o) idim igrape
	$(FC) -fixed $(FFLAGS) $(CPPFLAGS) -o $@ $(KROMEOBJ) utils/splitdump.F90 rdump.o wdump.o error.o indexx.o options.o chekopt.o labrun.o zzibm.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o)

mergedump: checksystem checkparams $(KROMEOBJ) utils/mergedump.F90 rdump.o wdump.o error.o indexx.o options.o chekopt.o labrun.o zzibm.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o) idim igrape
	$(FC) -fixed $(FFLAGS) $(CPPFLAGS) -o $@ $(KROMEOBJ) utils/mergedump.F90 rdump.o wdump.o error.o indexx.o options.o chekopt.o labrun.o zzibm.o $(SRCNI:.F90=.o) $(SRCNS:.F90=.o)

checkparams: updatedepends
   ifeq ($(DEBUG), yes)
	@echo "Debugging flags are ON"
   endif
   ifeq ($(DOUBLEPRECISION), yes)
	@echo "Flags set for DOUBLE PRECISION"
   else
	@echo "Flags set for SINGLE PRECISION"
   endif
   ifeq ($(ENDIAN), BIG)
	@echo "Flags set for conversion to BIG endian"
   endif
   ifeq ($(ENDIAN), LITTLE)
	@echo "Flags set for conversion to LITTLE endian"
   endif
   ifeq ($(openmp), yes)
	@echo "openMP flags are ON"
   endif
   ifeq ($(mpi), yes)
	@echo "MPI is ON"
   endif

# this is a way of making all files depend on their INCLUDEd files
# (creates a file called .depends which contain the appropriate dependency lines
#  this file is then included in the Makefile)
#

.depends: updatedepends

updatedepends:
	@if test -e .depends; then mv .depends .dependsold; fi
#       command below greps for INCLUDE in all .f files, uses sed to exclude comment lines, cut to extract text between quotes
#       then writes a dependency line for each .o file in Makefile format to the .depends file
#	@echo updating .depends
	@for x in *.[f,F] RT/*.[f,F]; do echo $${x/.[f,F]/.o}: `grep INCLUDE $$x | sed '/^[c,C,\!]/d' | cut -d\' -f 2` | sed 's/mpif.h//g' | sed 's/..\/COMMONS/COMMONS/g' >> .depends; done
	@diff --brief .dependsold .depends || ${MAKE} ${MAKECMDGOALS}
	@rm -f .dependsold

include .depends

clean:
	\rm -f *.o *.mod RT/Baraffe/*.o RT/Baraffe/*.mod RT/*.o RT/*.mod sph_tree_lf_gradh sph_tree_lf_gradh_RT sph_tree_lf sph_tree_rk_gradh sph_tree_rk_gradh_RT sph_tree_rk_gradh_RT_BAR sph_tree_rk ./utils/*.o ./utils/*/*.o maketables/*.o git_version.h
# DO NOT DELETE
