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
#include <time.h>

#include "Read_Input_Library.h"
#include "Allocate.h"
#include "Math_Functions.h"
#include "Struct.h"
#include "HC_Data_Format.h"
#include "QM_Density.h"
#include "Grid.h"
#include "Multipole.h"


#define GRID_PARM      "Grid_Parm.dat"
#define LMAX_MIN           4


int main( int argc, char *argv[])
{
	struct ELE_Data  ELE;
	struct MOL_Data  MOL;
	struct GRID_HI_Data GRID_HI;

	int LOG_FILE_PRESENT;
	char INPUT_FILE[1000];
	char initial[1000];
	char output[1000];
	char DEBUG_FILE[1000];
	char GRD_FILE[1000];
	char GRD_FILE2[1000];
	clock_t time1, time2;
	double cpu_time;
	int LMAX;

/*DEBUG*/
	int DEBUG, PRINT_FIELD;
	FILE *fp_debug;
	FILE *fp_grd;
	FILE *fp_grd2;


	DEBUG = 1;
	PRINT_FIELD = 0;


	if (argc == 2)
	{
		sscanf(argv[1], "%s", INPUT_FILE);
		LOG_FILE_PRESENT = 0;
	}
	else
	{
		printf("%s .input\n", argv[0]);
		exit(0);
	}

	
	Remove_Extension (INPUT_FILE, initial);
	sprintf(DEBUG_FILE, "%s.debug", initial);
	sprintf(GRD_FILE, "grd_pts.xyz", initial);
	sprintf(GRD_FILE2, "grd_pts2.xyz", initial);
	sprintf(MOL.inputfile, "%s", INPUT_FILE);
	time1 = clock();


		fp_grd = Open_File (GRD_FILE, "w");
		fp_grd2 = Open_File (GRD_FILE2, "w");
		/* fprintf(fp_grd, "\n\nTEST");*/
	if (DEBUG)
		fp_debug = Open_File (DEBUG_FILE, "w");

		printf("Reading %s\n", INPUT_FILE);
	Read_Input_File (INPUT_FILE, &ELE, &MOL, DEBUG, fp_debug);
	/*Defined in HC_Data_Format.c.  Reads INPUT_FILE and called from main().  This function reads the geometry of the 
	  molecule of interest (mol data).  It then reads elements, atomic nuclear charges (Z), and ionic states with multiplicities 
	  (ele data).  For each atom (element) in the molecule of interest, it indexes the element number for that atom.
	*/

	LMAX = MOL.lmax_atom;

	if (LMAX < LMAX_MIN)
		LMAX = LMAX_MIN;

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nLMAX = %d\n", LMAX);
	}



		if (DEBUG)	printf("Initialize_Grid %s\n", GRID_PARM);
	Initialize_Grid (GRID_PARM, fp_debug, DEBUG);
	/*Defined in Grid.c.  Reads in grid_data_file containing grid parameters: elements, Slater-Bragg radii for each element,
	  Lebedev order, and radial order.  It then creates a single Levedev angular grid and a single Gauss Chebyshev radial grid.
	*/	


		if (DEBUG)	printf("Get_QM_Density_Grid\n");
	Get_Grid (initial,&ELE,&MOL,&GRID_HI,DEBUG,fp_debug,fp_grd,fp_grd2);
	/*This function will make the molecular grid (3D).  It will then read in the wavefunction, density matrix, and ab initio (QM)
	  basis function from the formatted checkpoint file for the molecule of interest.  It will then evaluate the ab initio (QM) density
	  on each point of the molecular grid.  This function will also make atomic grids (3D) for each element.  It will then read in the wavefunction, 
	  density matrix, and ab initio (QM) basis function from the formatted checkpoint file of each element and each atomic ion defined for that element.
	  It will evaluate the ab initio (QM) density on the (3D) atomic grids, and then average over angular shells to get spherically symmetric
	  atomic densities which are function of radial distance only.

	  The molecular grid (3D), atomic radial grids (1D), molecular charge densities, and atomic charge densities are kept and store for later use.
	  The ab initio wave-function, density matrix, and basis sets are no longer needed after this function is called.
	*/


	if (DEBUG)
		fclose(fp_debug);

}


	
