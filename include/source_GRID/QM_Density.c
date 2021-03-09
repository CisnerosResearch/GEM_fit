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

#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include <math.h> 

#include "Read_Input_Library.h"
#include "Allocate.h"
#include "Math_Functions.h"

#include "Struct.h"
#include "Grid.h"
#include "Gauss_Density_Struct.h"

#define bohr_2_A 0.5291772083    
#define eA_2_bohr (1.0/0.5291772083)       
#define Pi 3.1415926535897932384626433832

/*Grid.c*/
extern int GRID_n_ang_pt;
extern double *GRID_ang_weight;
extern double **GRID_ang_r2d;

void Get_Atom_File_Name (char *initial, char *element, int charge);

struct QM_Data 
/*Stores the ab initio data defined in Gauss_Density_Struct.h*/
{
	int n_AO;
	struct Atomic_Orbital_Data	*AO;
	struct Atom_Data			ATOM;
	struct Wave_Function_Data	WF;
};

void Get_Grid (char *initial, struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid_hi, int DEBUG, FILE *fp_debug, FILE *fp_grd, FILE *fp_grd2 )
/*This function will make the molecular grid (3D).  It will then read in the wavefunction, density matrix, and ab initio (QM)
  basis function from the formatted checkpoint file for the molecule of interest.  It will then evaluate the ab initio (QM) density
  on each point of the molecular grid.  This function will also make atomic grids (3D) for each element.  It will then read in the wavefunction, 
  density matrix, and ab initio (QM) basis function from the formatted checkpoint file of each element and each atomic ion defined for that element.
  It will evaluate the ab initio (QM) density on the (3D) atomic grids, and then average over angular shells to get spherically symmetric
  atomic densities which are function of radial distance only.

  The molecular grid (3D), atomic radial grids (1D), molecular charge densities, and atomic charge densities are kept and store for later use.
  The ab initio wave-function, density matrix, and basis sets are no longer needed after this function is called.
*/
{
	struct QM_Data qm_mol;
	struct QM_Data **qm_ele;		/*[n_element][n_ion_state]*/


	char fchk_file[1000];
	int n, m, i, j, p;
	double sum, norm;
	double *d_temp;
	double **r2d_temp;
	char file_initial[1000];

/*Molecular Grid*/
		if (DEBUG)	printf("\tMake_Molecular_Grid\n");
	Make_Molecular_Grid ( (*mol).n_atom, (*mol).r2d, (*mol).element, 
		&( (*grid_hi).n_grid_pts), &( (*grid_hi).weight), &( (*grid_hi).r2d), fp_debug, 0);

	/*Grid is in terms of angstroms, convert to bohr*/
	for (i = 0; i < (*grid_hi).n_grid_pts; i++)
	{
	/*	for (p = 0; p < 3; p++)
			(*grid_hi).r2d[i][p] *= eA_2_bohr;*/
                        if ((*grid_hi).weight[i] != 0.0)
                        {
                        fprintf(fp_grd, "%15.8f %15.8f %15.8f %15.8f\n", (*grid_hi).r2d[i][0], (*grid_hi).r2d[i][1],
                                        (*grid_hi).r2d[i][2], (*grid_hi).weight[i] );
                        }
                        if ((*grid_hi).weight[i] != 0.0)
                        {
                        fprintf(fp_grd2, "%15.8f %15.8f %15.8f \n", (*grid_hi).r2d[i][0], (*grid_hi).r2d[i][1],
                                        (*grid_hi).r2d[i][2]  );
                        }
	}
	(*grid_hi).p_QM = D_Allocate_1D_Matrix ((*grid_hi).n_grid_pts);

}



