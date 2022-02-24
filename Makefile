
PROJECT_ROOT = ${PWD}

include config.mak

all: lib include GEM_site_site GEM_site_site2 GEM_fit GEM_calc_coefs

exes: GEM_site_site GEM_site_site2 GEM_numfit GEM_calc_coefs

GEM_site_site: GEM_site_site.o lib/libGEM.a

GEM_site_site2: GEM_site_site2.o lib/libGEM.a

GEM_fit: GEM_fit.o lib/libGEM.a

GEM_calc_coefs: GEM_calc_coefs.o lib/libGEM.a

.o:
	@echo $(FC) -o $@ $< $(LDFLAGS)
	$(FC) -o $@ $< $(LDFLAGS)

#.SILENT:

.PHONY: lib
lib:
	#mkdir solib
	cd include; make
	cd include/source_GRID; cc *.c -lm -o grid_gen.x; mv grid_gen.x ../../
	cd lib; make

depend:
	cd include; make
	cd include/source_GRID; cc *.c -lm -o grid_gen.x; mv grid_gen.x ../../
	cd lib; make depend

clean:
	/bin/rm -f sfpdb
	/bin/rm -f grid_gen.x
	cd include; make clean
	cd lib; make clean
	/bin/rm -f *~ *.o $(RESIDUE)

realclean:
	/bin/rm -f sfpdb
	/bin/rm -f grid_gen.x GEM_site_site GEM_site_site2 GEM_fit GEM_calc_coefs
	cd include; make realclean
	cd lib; make realclean
	/bin/rm -f solib/*
	/bin/rm -f *~ *.o $(RESIDUE)

