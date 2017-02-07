c This header selects whether to load the legacy "mpif.h" header file, or the 
c newer and more comprehensive mpi module. The latter has better interface support, which
c should spot a few bugs at compile time (plus it produces slightly smaller and faster code).
c However, not all combinations of compiler and mpi library support it properly.
c Currently we use the module for ifort, and the header for anything else.
c
c How to use this file:
c  - Put, below any USE statements but above the first INCLUDE or type definition, the following lines
c #ifdef MPIALL
c #include "mpi_sup.h"
c #endif
c
c  - If the subprogramme can use IMPLICIT NONE do the following instead
c #ifdef MPIALL
c #define IMPLICITNONE
c #include "mpi_sup.h"
c #else
c        IMPLICIT NONE
c #endif
c
c Special handling of IMPLICIT NONE is necessary because it must come
c *after* any USE statements but *before* any INCLUDE lines.


#ifndef MODMPI
#ifndef INCMPI
c If not preset above/elsewhere use heuristics
#ifdef __INTEL_COMPILER
c If an intel compiler we assume the module is available and stable
#define MODMPI
#else
c Otherwise we don't.
#define INCMPI
#endif
#endif
#endif

#ifdef MODMPI
      USE mpi
#endif

#ifdef IMPLICITNONE
      IMPLICIT NONE
#endif
#ifdef INCMPI
      INCLUDE "mpif.h"
#endif
