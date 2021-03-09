
include ../config.mak

LIB = libGEM.a
LDFLAGS = -L. -lGEM

DEPDIRS = strings math io nonbond gnonbond GEM_hermite      

DEPLIBS = $(DEPDIRS:%=libxam_%.a)
SOLIBS = $(DEPDIRS:%=../solib/libxam_%.so)

ifdef SHARED_LIBS__
all: $(DEPDIRS) $(SOLIBS)
else
all: $(DEPDIRS) $(LIB)
endif

.PHONY: $(DEPDIRS)

$(DEPDIRS):
	@cd $@; \
	if [ ! -f depend.mak ]; then \
	  echo "Making dependencies in $@:"; \
	  make depend; \
	fi; \
	echo "Making $@ ..."; \
	make

# The only way to merge multiple .a libraries is to
# extract all object files and remake a new lib.
$(LIB): $(DEPLIBS)
	@/bin/rm -rf lib_objs $(LIB)
	@mkdir lib_objs
	@for lib in $(DEPLIBS); do \
          mkdir lib_objs/$$lib ; \
	  cd lib_objs/$$lib ; \
          $(AR_EXTRACT) ../../$$lib ; \
          cd ../.. ; \
        done ;
	cd lib_objs; $(AR_CREATE) ../$(LIB) */*
	@/bin/rm -rf lib_objs
	$(RANLIB) $(LIB)

clean:
	@/bin/rm -rf *.a *.o _*.f lib_objs
	@/bin/rm -rf ../solib/*
	for lib in $(DEPDIRS); do \
          echo "Making clean in $$lib" ; \
	  cd $$lib ; \
          make clean ; \
          cd .. ; \
	done

realclean:
	@/bin/rm -rf *.a *.o _*.f lib_objs
	for lib in $(DEPDIRS); do \
          echo "Making realclean in $$lib" ; \
          cd $$lib ; \
	  make clean ; \
          make realclean ; \
          cd .. ; \
	done

