ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../..
endif
include $(PROJECT_ROOT)/config.mak

NAME = $(notdir $(shell pwd))

# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
F90_SOURCE= ch_coeff_mpoles_new.f90
F_SOURCE=  
C_SOURCE=

include ../libdir.mak
sinclude depend.mak

