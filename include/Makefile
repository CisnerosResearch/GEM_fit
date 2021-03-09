ifndef $(PROJECT_ROOT)
  PROJECT_ROOT=../
endif
include $(PROJECT_ROOT)/config.mak
 
NAME = $(notdir $(shell pwd))
 
# If all three source definitions are blank, libdir.mak
# will automatically include all *.c, *.f, and *.f90 files.
F90_SOURCE= definition.f90 cubes.f90 
F_SOURCE=
C_SOURCE=
 
all:
	$(FC) $(FSTDFLAGS) -o definition.o definition.f90
	$(FC) $(FSTDFLAGS) -o cubes.o cubes.f90
                                                                                
clean:
	/bin/rm -f definition.o definition.mod cubes.o cubes.mod
                                                                                
realclean:
	/bin/rm -f definition.o definition.mod cubes.o cubes.mod


