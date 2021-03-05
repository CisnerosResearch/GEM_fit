
#Define PROJECT_ROOT in the main Makefile

SHARED_LIBS=1
DEBUG=1

## Platform options: gfortran ifc
PLATFORM = gfortran
#PLATFORM = ifort

## Path to BLAS and LAPACK libraries
LAPACK=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
BLAS=/usr/lib/x86_64-linux-gnu/blas/libblas.so.3

#############################
## Compiler-specific flags ##
#############################

FPPFLAGS= -DPUBFFT -traditional -I$(PROJECT_ROOT)/include -DDEBUG 
#LDFLAGS= -L$(PROJECT_ROOT)/lib -lGEM -lz
ifeq ($(PLATFORM),ifort)
LDFLAGS= -L$(PROJECT_ROOT)/solib \
    -lxam_io -lxam_math -lxam_strings \
    -lxam_gnonbond \
    -lxam_nonbond -lxam_GEM_hermite\
    -lz -lifcore \
    $(LAPACK) $(BLAS)
endif
## Yes, you need the repeating LDFLAGS
ifeq ($(PLATFORM),gfortran)
LDFLAGS= -L$(PROJECT_ROOT)/solib/ \
    -lxam_gnonbond -lxam_nonbond -lxam_GEM_hermite \
    -lxam_io -lxam_math -lxam_strings \
    -lxam_io -lxam_math -lxam_strings \
    -lxam_gnonbond -lxam_nonbond -lxam_GEM_hermite \
    -lz -lgfortran \
    $(LAPACK) $(BLAS)
endif


FPP=$(PROJECT_ROOT)/tools/xfpp

AR_CREATE=ar rc
AR_EXTRACT=ar x
RANLIB=ranlib

################################################################
# For each platform, set these:
# CC=	        C compiler
# CFLAGS=	Standard C flags
# CPPFLAGS=	C preprocessor flags
# FC=	        FORTAN compiler
# FSTDFLAGS=	Standard F flags
# F90FLAGS=	Flags for .f90 free-format files
# F77FLAGS=	Flags for .f fixed-format files, implicit variables allowed
# FFLAGS =	Unoptimized compile flags
# FOPTFLAGS=	Optimized compile flags
# LD=	        linker FORTRAN based executables
# LDCC=	        linker for C based executables
# RESIDUE=	Additional junk files left by the compiler, etc.
# AR=	        name of the 'ar' static library tool
# RANLIB=	name of the 'ranlib' tool. Use /bin/true" if none.
################################################################


ifeq ($(PLATFORM),gfortran)

CC= gcc -g
CFLAGS=-O2 -fPIC 
CPPFLAGS= -I$(PROJECT_ROOT)/include
FC= gfortran

# FLAGS para gfortran
#FSTDFLAGS= -c -fbounds-check -w -fPIC -fno-range-check
#FSTDFLAGS= -c -w -fPIC -fno-range-check
FSTDFLAGS= -c -w -fPIC -fno-range-check -std=legacy

F90FLAGS=
F77FLAGS=
FFLAGS = $(FSTDFLAGS) # -O0
FOPTFLAGS= $(FSTDFLAGS) -O3 -fPIC

LD= gfortran

LDCC= gcc -g
RESIDUE=

endif

ifeq ($(PLATFORM),ifort)

CC= gcc -g
CFLAGS=-O2 -fPIC
CPPFLAGS= -I$(PROJECT_ROOT)/include
FC= ifort

FSTDFLAGS= -c -w -O3 -tpp6 -align 

F90FLAGS=
F77FLAGS=
FFLAGS = $(FSTDFLAGS) # -O0
FOPTFLAGS= $(FSTDFLAGS) -O3 -fPIC

LD= ifort

LDCC= gcc -g
RESIDUE=

endif

ifeq ($(PLATFORM),ifc)

CC= icc
CFLAGS=-O2
CPPFLAGS=
FC= ifc -u
#FSTDFLAGS= -g -fpp1 -u -w90 -w95 -cm -FR $(FPPFLAGS)
FSTDFLAGS= -fpp1 -u -w90 -w95 -cm -FR $(FPPFLAGS)
FFLAGS = -O0 $(FSTDFLAGS)
FOPTFLAGS= -O3 $(FSTDFLAGS)
LD= ifc 
LDCC= icc 
RESIDUE=ifc?????? work.pc work.pcl

endif


.SUFFIXES:  .f90

.f90.o: $<
	@echo $(FC) -c \$$FPPFLAGS \$$F90FLAGS \$$FOPTFLAGS -o $@ $<
	$(FC) -c $(FPPFLAGS) $(F90FLAGS) $(FOPTFLAGS) -o $@ $<

.f.o:   $<
	@echo $(FC) -c \$$FPPFLAGS \$$F77FLAGS \$$FOPTFLAGS -o $@ $<
	$(FC) -c $(FPPFLAGS) $(F77FLAGS) $(FOPTFLAGS) -o $@ $<

.c.o:
	@echo $(CC) -c \$$CFLAGS \$$CPPFLAGS -o $@ $<
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<
.a.so:
	$(CC) -shared --version-script=exports.scr -o $@ $<


