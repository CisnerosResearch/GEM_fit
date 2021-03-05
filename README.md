# GEM_fit

Program to calculate GEM Hermites and multipoles with either a numerical or
analytical fit from a Gaussian formatted checkpoint file.

## Installation

Download this repository:

```bash
$ git clone https://github.com/CisnerosResearch/GEM_fit.git
```

Go into the main directory (GEM_fit) and:

-  Modify "setup" and "setup.profile" so that it points to the correct
directory

-  Source the corresponding file based on the shell you use.
(If you don't know what this means do `echo $SHELL`)

   For tcsh:
   ```sh
   $ source setup
   ```
   For bash:
   ```bash
   $ source setup.profile
   ```

-  Modify `config.mak` to choose compiler, either ifort or gfortran.
Note: this selection only needs to be made in the `PLATFORM = ifort`
definition.
Do not do a global string replacement.

-  Modify `config.mak` to change location for BLAS and LAPACK.
Note: the program will not generate any executable if it cannot link to both.
To find these on your machine, you can use `locate`:
   ```bash
   $ locate liblapack.so.3
   /usr/lib/x86_64-linux-gnu/liblapack.so.3
   /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
   /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

   $ locate libblas.so.3
   /usr/lib/x86_64-linux-gnu/libblas.so.3
   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3
   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
   ```

-  You can then use `make` in the top-level directory.
There's no need to do a `./config.mak`.
   ```bash
   $ make
   ```

After making you should see four executables:
- `GEM_fit`: for the numerical fitting
- `GEM_calc_coefs`: for analytical fitting
- `GEM_site_site`: calculates intermolecular interactions with the hermites or
multipoles.
- `GEM_site_site2`: calculates intermolecular interactions with the hermites or
multipoles.

