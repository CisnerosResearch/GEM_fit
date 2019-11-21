# GEM_fit

Program to calculate GEM Hermites and multipoles with either a numerical or
analytical fit from a Gaussian formatted checkpoint file.

To install:

Uncompress the archive:

```
tar xvfj GEM_fit_ana-num.tar.bz2 
```

After uncompressing go into the main directory (GEM_fit) and:

-  modify "setup" and "setup.profile" so that it points to the correct
directory

-  source the corresponding file (setup for tcsh, setup.profile for bash)

-  modify "config.mak" to choose compiler, either ifort or gfortran. 

-  modify "config.mak" to change location for BLAS and LAPACK, note that
 the program will not generate any executable if it cannot link to both.

-  make

After making you should see three executables: "GEM_fit, GEM_calc_coefs,
GEM_site_site and GEM_site_site2". The first is for the numerical fitting,
the second for analytical fitting and the last two to calculate 
intermolecular interactions with the hermites or multipoles.

______________________________________________________________________________

