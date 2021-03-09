# General Makefile directives used by all lib subdirectories.
# Relies on pattern rules for compiling object files from source.

ifdef SHARED_LIBS
all: solib
else
all: lib
endif

LIB=$(PROJECT_ROOT)/lib/libxam_$(NAME).a
lib: $(LIB) depend.mak

# Auto-detect sources only if none are defined:
ifndef F90_SOURCE
  ifndef F_SOURCE
    ifndef C_SOURCE
      F90_SOURCE = $(wildcard *.f90) 
      F_SOURCE = $(wildcard *.f) 
      C_SOURCE = $(wildcard *.c)
      depend: info
    endif
  endif
endif


OBJECTS= $(F90_SOURCE:.f90=.o) $(F_SOURCE:.f=.o) $(C_SOURCE:.c=.o)

MAKEDEPEND=$(PROJECT_ROOT)/tools/GEM_makedep

.PHONY: depend
depend: depend.mak

depend.mak:
	@echo "GENERATING DEPENDENCIES:"
	@echo  C_SOURCE = $(C_SOURCE)
	@echo  F_SOURCE = $(F_SOURCE)
	@echo  F90_SOURCE = $(F90_SOURCE)
	$(MAKEDEPEND) $(FPPFLAGS) $(F90_SOURCE) $(F_SOURCE)
	$(MAKEDEPEND) $(CPPFLAGS) $(C_SOURCE)
	/bin/mv -f rules.mak depend.mak
	@echo "Created depend.mak. You may now re-run make."
	@/bin/false


# The root level Make checks only the directory time,
# so we go ahead and update the lib based on the current dir
# date, even if no sources change.

$(LIB): $(OBJECTS) .
	@echo Creating static library: $(LIB)
	@/bin/rm -f $(LIB)
	@$(AR_CREATE) $(LIB) $(OBJECTS)
	@$(RANLIB) $(LIB)


# An export script in gcc looks like this:
# { global: export1; export2; expat*; local: *; };

SONAME = libxam_$(NAME).so
SOLIB = $(PROJECT_ROOT)/solib/libxam_$(NAME).so
solib: $(LIB) $(SOLIB) depend.mak

$(SOLIB): $(LIB)
	@echo Creating shared library: $(SONAME)
	@if [ -f exports.txt ]; then \
	  awk 'BEGIN{print "{global:"}substr($$1,1,1)!="#"{print $$1"; "$$1"_; "$$1"__;"}END{print "local: *;};"}' \
	      exports.txt > exports.la ; \
	  echo $(FC) -shared -Wl,-soname,$(SONAME),--version-script=exports.la -o $(SOLIB) $(OBJECTS) ; \
	  $(FC) -shared -Wl,-soname,$(SONAME),--version-script=exports.la -o $(SOLIB) $(OBJECTS) ; \
	else \
	  echo $(FC) -shared -Wl,-soname,$(SONAME) -o $(SOLIB) $(OBJECTS) ; \
	  $(FC) -shared -Wl,-soname,$(SONAME) -o $(SOLIB) $(OBJECTS) ; \
	fi

clean:
	/bin/rm -f $(LIB) $(OBJECTS) $(TESTS)
	/bin/rm -f _* *.o *~ $(RESIDUE)

realclean:
	/bin/rm -f depend.mak

