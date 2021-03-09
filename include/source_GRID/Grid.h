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


void Initialize_Grid (char *grid_data_file, FILE *fp_debug, int DEBUG);
/*Reads in grid_data_file containing grid parameters: elements, Slater-Bragg radii for each element,
  Lebedev order, and radial order.  It then creates a single Levedev angular grid and a single Gauss Chebyshev radial grid.
*/

void Debug_Make_Grid (int DEBUG, int PRINT_GRID, FILE *fp_debug);
/*assuming Initialize_Grid is already set*/

void Make_Molecular_Grid (int n_atom, double **mol_r2d, char **mol_element, 
	int *tot_n_pt, double **tot_weight, double ***tot_r2d, FILE *fp_debug, int DEBUG);
/*Once the angular and radial grids have been created in Initialize_Grid, this function will create a molecular grid
  using the Slater-Bragg radii and Becke's smoothing function method.

   INPUT:  n_atom, mol_r2d[n_atom][3], mol_element[n_atom]
   OUTPUT: tot_n_pt, tot_weight[tot_n_pt], tot_r2d[tot_n_pt][3]
*/


void Make_Radial_Atomic_Grid (char *element, int *n_pt_rad, double **radius, double **weight_radius);
/* INPUT:  element
   OUTPUT: n_pt_rad, radius[n_pt_rad], weight_radius[n_pt_rad]

*/


