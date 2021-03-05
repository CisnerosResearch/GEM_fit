ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../..
endif
include $(PROJECT_ROOT)/config.mak

NAME = $(notdir $(shell pwd))

# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
F90_SOURCE = file_io.f90 t_get_next_line.f90
F_SOURCE =
C_SOURCE = direct_io.c

include ../libdir.mak
sinclude depend.mak

