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
#ENDIAN= default

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

ifeq ($(SYSTEM), ompi)
    FC= mpif90
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
       FC= mpiifort
#       FC= ifort -lmpi
    else
       FC= ifort
    endif
    FFLAGS= -O3 -xSSSE3 -ipo -mcmodel=medium -i-dynamic -WB  \
            -warn nousage -warn uninitialized -warn unused -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check
#To use above MPI check need to:
#source /sw/sdev/itac/7.2.0.011/bin/itacvars.csh
#
    OMPFLAG= -openmp #-DZEN
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
#    ZEN_PIN = pin_to_cpu.o
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), complexity)
    ifeq ($(mpiall), yes)
       FC= mpiifort
#       FC= ifort -lmpi
    else
       FC= ifort
    endif
    FFLAGS= -O3 -axAVX -ipo -mcmodel=medium -WB \
            -warn nousage -warn uninitialized -warn unused -warn truncated_source\
            -warn interfaces -nogen-interfaces
# -check
#To use above MPI check need to:
#source /sw/sdev/itac/7.2.0.011/bin/itacvars.csh
#
    OMPFLAG= -openmp #-DZEN
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
#    ZEN_PIN = pin_to_cpu.o
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), ukaff1a)
    ifeq ($(mpiall), yes)
        FC= mpif77 -qnosave
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
       FC= mpif90
    endif
    FFLAGS= -m64 -mcmodel=medium -O3 -Wall
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
    FFLAGS= -O3 -mcmodel=medium -shared-intel -WB -Vaxlib -w90 -w95 -cm -warn uninitialized -warn unused -warn truncated_source 
    OMPFLAG= -openmp
    DBLFLAG= -r8
    DEBUGFLAG= -check all -traceback -g -fpe0 -fpstkchk
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif

# Set other optional flags depending on settings

ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
endif

ifeq ($(ENDIAN), BIG)
    FFLAGS += ${ENDIANFLAGBIG}
endif

ifeq ($(ENDIAN), LITTLE)
    FFLAGS += ${ENDIANFLAGLITTLE}
endif

ifeq ($(openmp), yes)
    FFLAGS += ${OMPFLAG}
endif

ifeq ($(mpi), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIALL
    FPPFLAGS+=${FPP_SEP}-DMPI
endif
ifeq ($(mpicopy), yes)
    FPPFLAGS+=${FPP_SEP}-DMPIALL
    FPPFLAGS+=${FPP_SEP}-DMPICOPY
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

# do double precision flag last (append only to FFLAGS)

ZZFFLAGS := ${FFLAGS}
ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAG}
endif

# add preprocessor flags (and prefix required for some compilers)
FFLAGS += ${FPP_PREFIX}${FPPFLAGS}

FFLAGSKEEP := ${FFLAGS}

KERNEL= ktable.f
OBJKERNEL= $(KERNEL:.f=.o)

# these are the sources common to all compilations

SOURCES= sph.F accrete_P.F addump.f angmom.f balance.F \
	bonnorebert.f boundry.f buildsightlines.f \
	cartdis.f cartmas.f cartpres.f cartvel.f chanref.f chekopt.f \
	classify_MPI.F condense.f constan.f \
	coolcurve.f coolinmo.F cool_func.F \
	coriol.f cyldis.f cylmas.f cylpres.f cylvel.f discboundary.f \
	divBclean_subcycle.F divBsubcycle.F \
	dust_implicit.F dust_drag.F dust_norm_rhoG.F \
	dragtable.F endrun.f energ.f energ_cooling.f eospg.f error.f \
	evol.f externf.f externBfield.f \
	extract.f file.f gapfactor.f getneigh.f \
	getcd.f getcv.f get1overmu.f getkappa.f \
	gforsa_P.f gforsn_P.f gforspt_P.f \
	ghostp.f ghostp1.f ghostp2.f ghostp3.F ghostp11.f ghostp100.f \
	ghostp102.f ghostp104.f \
	gptall.f gptfromgas.f gpti.f header.F hchem.F hillmass.f \
	homexp.f indexx.f indexx2.f indexx3.f indexxi2.f inform.f \
	inopts.f insulate_TREE_P.F integs.F $(KERNEL) \
	labrun.F lunit_double.f \
	mainop.f makedust.f makeplot.f mesop.f mhdparams_mhd.f \
	modif.f mtree_P.f newrun.f options.f \
	optimise_column_density.f phoenix.f phoenix2.f \
	place.f planetpotential.f planet_setup.f planetesimals.f potential.f \
	preset.F \
	prout.F ran1.f rdump.F rdump_wrapper.F read_dragtable.f reduce.f \
	maketables/specificheat.f resistivity.f resistivity_implicit.F save.f \
	scaling.f secmes.f \
	semimajoraxis.f setBfield_mhd.f setpart.F setup_torus.f \
	smoothd_P.f sphericalposition.f \
	sphdis.f sphmas.f sphpres.f sphvel.f spline.F thermeq.f \
	toten.f treef_P.F unifdis.F unit.f uset.f \
	wdump.F wdump_wrapper.F wrinsph.F write_dragtable.f zeusread.f

# these are the sources required for the gradh code to compile

SRCGRADH= derivi_P_gradh_MPI.F forcei_P_gradh_MPI.F  \
        densityiterate_gradh_MPI.F divv_gradhsoft.F \
	hcalc_gradh.f dissipative_source.F
	#divBclean_mhd.f divBdiffuse.f divBiterate.f \
	#directsum_poisson_vec_mhd_clean.f 

# these are the extra routines required for the radiative transfer code to compile
        
SRCRT= ass.f RT/newgsimpl_MPI.f RT/boundaryparticle.f quartic_gs1t.f \
	RT/getblock.f RT/dust_temperature.f column_density.F 
SRCNORT = ass_dummy.f
RTFPPFLAGS=${FPP_SEP}-DRT

# these are the extra routines required for Monte-Carlo radiative transfer
        
SRCTORUS= montecarloRT.F

# these are the extra routines required for non-ideal MHD
ifeq ($(nonideal), yes)
SRCNONIDEAL= calc_curlB.F
SRCNI = nicil.F90
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

TOBJECTS =  $(SRCNI:.F90=.o) $(SOURCES:.f=.o) $(ZZFILE:.f=.o) $(SRCNONIDEAL:.F=.o)
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
	$(FC) $(FFLAGS) -c -o $(SRCNI:.F90=.o) $<
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

sph_tree_rk_gradh_RT_BAR: $(OBJGRADHRKRT) $(OBJBAR)
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

convertmpi: checksystem checkparams utils/convert_to_mpi.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o lunit_double.o
	$(FC) $(FFLAGS) -o $@ utils/convert_to_mpi.o wdump.o rdump.o error.o classify_MPI.o file.o options.o chekopt.o lunit_double.o

# convert from full to small dump file

convertfs: checksystem checkparams utils/convert_full_small.o wdump.o rdump.o error.o balance.o
	$(FC) $(FFLAGS) -o $@ utils/convert_full_small.o wdump.o rdump.o error.o balance.o

# convert from old dump file to new full dump file

convertof: checksystem checkparams utils/convert_old_full.o wdump.o error.o
	$(FC) $(FFLAGS) -o $@ utils/convert_old_full.o wdump.o error.o

# 3D movie processing

movie3DRG: checksystem checkparams utils/movie3DRG/imcol3DNG_zoom_DT.o rdump.o error.o utils/movie3DRG/grid_P_DT.o utils/movie3DRG/buildcol.o utils/movie3DRG/origin.o
	$(FC) $(FFLAGS) -o $@ utils/movie3DRG/imcol3DNG_zoom_DT.o rdump.o error.o utils/movie3DRG/grid_P_DT.o utils/movie3DRG/buildcol.o utils/movie3DRG/origin.o

# Disc parameter extraction

disc: checksystem checkparams utils/disc_parameters.o rdump.o error.o indexx2.o
	$(FC) $(FFLAGS) -o $@ utils/disc_parameters.o rdump.o error.o indexx2.o

ptmass: utils/ptmassmerge.o
	$(FC) $(FFLAGS) -o $@ utils/ptmassmerge.o

massvsangmom: checksystem checkparams utils/massvsangmom.o utils/readdump_sphNG.o indexx2.o idim
	$(FC) $(FFLAGS) -o $@ utils/massvsangmom.o utils/readdump_sphNG.o indexx2.o

masswithtime: checksystem checkparams utils/masswithtime.o utils/readdump_sphNG.o idim
	$(FC) $(FFLAGS) -o $@ utils/masswithtime.o utils/readdump_sphNG.o

printkernel: checksystem checkparams constan.o $(OBJKERNEL) utils/printkernel.o
	$(FC) $(FFLAGS) -o $@ constan.o $(OBJKERNEL) utils/printkernel.o

# set up particle distribution without using entire SPH executable

setup_alone: checksystem checkparams utils/setup_alone.o wdump.o error.o file.o chekopt.o lunit_double.o setpart.o $(OBJKERNEL) unifdis.o condense.o cartmas.o cylmas.o sphmas.o cartdis.o cyldis.o sphdis.o cartpres.o cylpres.o sphpres.o eospg.o setBfield_mhd.o inopts.o cartvel.o cylvel.o sphvel.o gapfactor.o uset.o wrinsph.o constan.o ran1.o unit.o $(ZZFILE:.f=.o)
	$(FC) $(FFLAGS) -o $@ utils/setup_alone.o wdump.o error.o file.o chekopt.o lunit_double.o setpart.o $(OBJKERNEL) unifdis.o condense.o cartmas.o cylmas.o sphmas.o cartdis.o cyldis.o sphdis.o cartpres.o cylpres.o sphpres.o eospg.o setBfield_mhd.o inopts.o cartvel.o cylvel.o sphvel.o gapfactor.o uset.o wrinsph.o constan.o ran1.o unit.o $(ZZFILE:.f=.o)

getdeltat: checksystem checkparams utils/getdeltat.o utils/readdump_sphNG.o idim
	$(FC) $(FFLAGS) -o $@ utils/getdeltat.o utils/readdump_sphNG.o

# this routine needs to be compiled with single precision flags only

zzsun_f77_double2.o: zzsun_f77_double2.f
	$(FC) -c $(ZZFFLAGS) zzsun_f77_double2.f -o zzsun_f77_double2.o

checksystem:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling sphNG for $(SYSTEM) system..........."
	@echo ""
	@echo "calling createsvnversion"
	@./createsvnversion.sh
	@echo "...done"
   else
	@echo ""
	@echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=>set the environment variable SYSTEM to one listed "
	@echo "  in the Makefile and try again"
	@echo ""
	quit
   endif

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
	@for x in *.[f,F]; do echo $${x/.[f,F]/.o}: `grep INCLUDE $$x | sed '/^[c,C,\!]/d' | cut -d\' -f 2` | sed /mpif.h/d >> .depends; done
	@diff --brief .dependsold .depends || ${MAKE} ${MAKECMDGOALS}
	@rm -f .dependsold

include .depends

clean:
	\rm -f *.o *.mod RT/Baraffe/*.o RT/Baraffe/*.mod RT/*.o RT/*.mod sph_tree_lf_gradh sph_tree_lf_gradh_RT sph_tree_lf sph_tree_rk_gradh sph_tree_rk_gradh_RT sph_tree_rk ./utils/*.o ./utils/*/*.o maketables/*.o svn_version.h
# DO NOT DELETE
