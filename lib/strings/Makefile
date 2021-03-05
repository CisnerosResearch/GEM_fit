ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../..
endif
include $(PROJECT_ROOT)/config.mak

NAME = $(notdir $(shell pwd))

# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
#F90_SOURCE = ftype.f90 numfields.f90 strings.f90 tokens.f90
F90_SOURCE = numfields.f90 strings.f90 tokens.f90
F_SOURCE =
C_SOURCE = str_c.c

include ../libdir.mak
sinclude depend.mak

