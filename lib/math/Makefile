ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../..
endif
include $(PROJECT_ROOT)/config.mak

NAME = $(notdir $(shell pwd))

#F90_SOURCE = math.f90 spline.f90 func_table.f90 vector_3d.f90 \
#             geometry.f90 utility.f90
# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
F90_SOURCE = math.f90  vector_3d.f90 \
             geometry.f90 utility.f90
F_SOURCE = cubspl.f
C_SOURCE =

include ../libdir.mak
sinclude depend.mak

