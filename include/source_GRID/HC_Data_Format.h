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

void Read_Input_File (char *filename, struct ELE_Data *ele, struct MOL_Data *mol, int DEBUG, FILE *fp_debug);
/*Reads INPUT_FILE and called from main().  This function reads the geometry of the molecule of interest (mol data).  
  It then reads elements, atomic nuclear charges (Z), and ionic states with multiplicities (ele data).  For each atom (element) in
  the molecule of interest, it indexes the element number for that atom.
*/

void Read_Gauss_Input (char *initial, struct ELE_Data *ele, struct MOL_Data *mol, int DEBUG, FILE *fp_debug);
/*This function will first check to see if all necessary Gaussian 03 or 09 output files (.log and .fchk) are available in the working directory.
  If the Gaussian output files are available, this function will read in the ChelpG charges (if (*mol).DO_CHELPG == 1) and the molecular multipoles.
  If the Gaussian output files are not available, this function will create the necessary input files (.com) for Gaussian 03 or 09 for
  the molecule of interest and the atomic ions.  This function will then exit and wait for the user to execute Gaussian on all of the created .com files.
  The user will also need to convert the checkpoint (.chk) files created by Gaussian into formatted checkpoint (.fchk) files.
*/


double ***Allocate_Cart_Tensor (int lmax);
/*Allocates Cart[l1][l2][l3], l1+l2+l3 <= lmax*/


void Allocate_Data (struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid_hi,
	int PRINT_FIELD, int DEBUG, FILE *fp_debug);
/*This function will allocate data structures needed in the evaluation of Hirshfeld & Hirshfeld-Iterated charges,
  the evaluation of atomic multipoles, and the evaluation of the molecular multipoles
*/

void Print_Output (double cpu_time, char *filename, struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid_hi,
	int PRINT_FIELD, int DEBUG, FILE *fp_debug);
/*This function prints the output file containing the molecular geometry, Hirshfeld & Hirshfeld-Iterated atomic charges, atomic multipoles (both
  spherical tensor and traceless Cartesian, the QM molecular multipoles, and the Hirshfeld, Hirshfeld-Iterated molecular multipoles for different
  atomic multipole ranks*/



