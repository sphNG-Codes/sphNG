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
SHELL = /bin/sh

# can comment out the following lines and instead set
# the parameters as environment variables
debug= yes
doubleprecision= yes
# endian can be "BIG", "LITTLE" or anything else which has no effect
ENDIAN= default

ifeq ($(SYSTEM), ukaff1a)
    FC= xlf_r -qnosave
    FFLAGS= -O3 -q64 -qextname -qsmp=noauto
    DBLFLAG= -qrealsize=8
#    DEBUGFLAG= -qflttrap=enable:invalid:zerodivide -g -C -qsigtrap -qfloat=nans 
    DEBUGFLAG= -C
    KNOWN_SYSTEM=yes
    ZZFILE= zzsun_f77_double1.f zzsun_f77_double2.f
endif

ifeq ($(SYSTEM), sunf95)
    FC= sunf95
    FFLAGS= -fast
    DBLFLAG= -xtypemap=real:64,double:64
    DEBUGFLAG= -g -C -ftrap=common
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
    FFLAGS= -O3 -Wall -mcpu=powerpc -mpowerpc
    DBLFLAG= -fdefault-real-8
    DEBUGFLAG= -g ##-ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
    ZZFILE= zzibm.f
endif

ifeq ($(SYSTEM), ifort)
    FC= ifort
    FFLAGS= -O3 -Vaxlib -w90 -w95 -cm
    DBLFLAG= -r8
    DEBUGFLAG= -C
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
    ZZFILE= zzibm.f
endif


# Set other optional flags depending on settings

ifeq ($(debug), yes)
    FFLAGS += ${DEBUGFLAG}
endif

ifeq ($(ENDIAN), BIG)
    FFLAGS += ${ENDIANFLAGBIG}
endif

ifeq ($(ENDIAN), LITTLE)
    FFLAGS += ${ENDIANFLAGLITTLE}
endif

# do double precision flag last (append only to FFLAGS)

ZZFFLAGS := ${FFLAGS}
ifeq ($(doubleprecision), yes)
    FFLAGS += ${DBLFLAG}
endif

# these are the sources common to all compilations

SOURCES= sph.f accrete_P.f addump.f angmom.f boundry.f \
	cartdis.f cartmas.f cartpres.f cartvel.f chanref.f \
	chekopt.f condense.f constan.f coriol.f cyldis.f cylmas.f cylpres.f \
	cylvel.f endrun.f energ.f \
	eospg.f error.f evol.f externf.f extract.f \
        file.f getneigh.f \
	getcv.f get1overmu.f getkappa.f \
	gforsa_P.f gforsn_P.f gforspt_P.f \
	ghostp1.f ghostp2.f ghostp3.f ghostp11.f ghostp100.f \
        gptall.f hdot_P.f \
	header.f homexp.f indexx.f indexx2.f indexx3.f indexxi2.f inform.f \
	inopts.f insulate_TREE_P.f integs.f ktable.f labrun.f lunit_double.f \
	mainop.f mesop.f mhdparams_mhd.f modif.f mtree_P.f \
	newrun.f options.f phoenix.f phoenix2.f \
	place.f preset.f \
	prout.f ran1.f rdump.f reduce.f save.f \
        scaling.f secmes.f setBfield_mhd.f setpart.f smoothd_P.f \
	sphdis.f sphmas.f sphpres.f sphvel.f \
	timestep.f toten.f treef_P.f unifdis.f unit.f \
	wdump.f wrinsph.f

# these are the sources required for the gradh code to compile

SRCGRADH= densityiterate_gradh.f derivi_P_gradh.f \
	divBclean_mhd.f divBdiffuse.f divBiterate.f \
	directsum_poisson_vec_mhd_clean.f forcei_P_gradh.f \
        hcalc_gradh.f 

# these are the sources required for the non-gradh code to compile
        
SRCNONGRADH= densityi_P.f derivi_P.f forcei_P.f hcalc.f

# these are the extra routines required for the radiative transfer code to compile
        
SRCRT= ass.f newgsimpl_P5_sor.f quartic_gs1t.f derivi_P_RT.f

# these are the extra routines required for Monte-Carlo radiative transfer
        
SRCTORUS= montecarloRT.f

# timestepping routine

SRCRK = step_P.f revtree_P.f
SRCLF = step_P_lf_fast2_gradh.f revtree_P_fast2.f

# construct object files

OBJECTS = $(SOURCES:.f=.o) $(ZZFILE:.f=.o)
OBJLF = $(OBJECTS) $(SRCNONGRADH:.f=.o) $(SRCLF:.f=.o)
OBJRK = $(OBJECTS) $(SRCNONGRADH:.f=.o) $(SRCRK:.f=.o)
OBJRT = $(OBJECTS) $(SRCNONGRADH:.f=.o) $(SRCRK:.f=.o) $(SRCRT:.f=.o)
OBJGRADHLF = $(OBJECTS) $(SRCGRADH:.f=.o) $(SRCLF:.f=.o)
OBJGRADHRK = $(OBJECTS) $(SRCGRADH:.f=.o) $(SRCRK:.f=.o)
OBJTORUS = $(OBJECTS) $(SRCNONGRADH:.f=.o) $(SRCLF:.f=.o) $(SRCTORUS:.f=.o)

lf: checksystem checkparams sph_tree_lf

rk: checksystem checkparams sph_tree_rk

RT: checksystem checkparams sph_tree_RT

torus: checksystem checkparams sphtorus

gradh: checksystem checkparams sph_tree_lf_gradh 

gradhlf: checksystem checkparams sph_tree_lf_gradh 

gradhrk: checksystem checkparams sph_tree_rk_gradh 

all: lf rk RT gradhlf gradhrk

sph_tree_lf_gradh: $(OBJGRADHLF)
	@echo ""
	@echo "This is the gradh, leapfrog version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHLF)

sph_tree_rk_gradh: $(OBJGRADHRK)
	@echo ""
	@echo "This is the gradh, Runge Kutta version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJGRADHRK)

sph_tree_lf: $(OBJLF)
	@echo ""
	@echo "This is the non-gradh, leapfrog version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJLF)

sph_tree_rk: $(OBJRK)
	@echo ""
	@echo "This is the non-gradh, Runge-Kutta version"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJRK)
        
sph_tree_RT: $(OBJRT)
	@echo ""
	@echo "This is the Radiative transfer code (Runge Kutta by default)"
	@echo ""
	$(FC) $(FFLAGS) -o $@ $(OBJRT)

sphtorus: $(OBJTORUS)
	@echo ""
	@echo "This is the TORUS code (Leapfrog by default)"
	@echo ""
	$(FC) $(FFLAGS) $(TORUSLIB) -o $@ $(OBJTORUS)

binary: utils/binary.o
	$(FC) $(FFLAGS) -o $@ utils/binary.o

# this routine needs to be compiled with single precision flags only

zzsun_f77_double2.o: zzsun_f77_double2.f
	$(FC) -c $(ZZFFLAGS) zzsun_f77_double2.f -o zzsun_f77_double2.o
        
checksystem:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling sphNG for $(SYSTEM) system..........."
	@echo ""
   else
	@echo ""
	@echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=>set the environment variable SYSTEM to one listed "
	@echo "  in the Makefile and try again"
	@echo ""
	quit
   endif

checkparams: updatedepends
   ifeq ($(debug), yes)
	@echo "Debugging flags are ON"
   endif
   ifeq ($(doubleprecision), yes)
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
	@for x in *.f; do echo $${x/.f/.o}: `grep INCLUDE $$x | sed '/^[c,C,\!]/d' | cut -d\' -f 2` >> .depends; done
	@diff --brief .dependsold .depends || ${MAKE} ${MAKECMDGOALS}
	@rm -f .dependsold

include .depends
        
clean:
	\rm -f *.o sph_tree_lf_gradh sph_tree_lf sph_tree_rk_gradh sph_tree_rk
# DO NOT DELETE
