/*
HPAM: Hirshfeld Partitioned Atomic Multipoles
Copyright (C) 2011  Dennis M. Elking

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Struct.h"

void Get_Grid (char *initial, struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid_hi, int DEBUG, FILE *fp_debug, FILE *fp_grd, FILE *fp_grd2);
/*This function will make the molecular grid (3D).  It will then read in the wavefunction, density matrix, and ab initio (QM)
  basis functions from the formatted checkpoint file for the molecule of interest.  It will then evaluate the ab initio (QM) density
  on each point of the molecular grid.  This function will also make atomic grids (3D) for each element.  It will then read in the wavefunction, 
  density matrix, and ab initio (QM) basis functions from the formatted checkpoint file of each element and each atomic ion defined for that element.
  It will evaluate the ab initio (QM) density on the (3D) atomic grids, and then average over angular shells to get spherically symmetric
  atomic densities which are function of radial distance only.

  The molecular grid (3D), atomic radial grids (1D), molecular charge densities, and atomic charge densities are kept and store for later use.
  The ab initio wave-function, density matrix, and basis sets are no longer needed after this function is called.
*/


