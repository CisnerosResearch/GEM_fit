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

void Initialize_Atomic_Multipole (int lmax, int DEBUG, FILE *fp_debug);
/*Calculates the constant defined by A[l][m] = sqrt( (l+m)!*(l-m)! )*/

void Calculate_Multipoles (struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid,
	int DEBUG, FILE *fp_debug);
/*This function converts the ab initio (QM) 'pure' Cartesian molecular multipoles obtained from the Gaussian (.log) file
  into traceless Cartesian molecular multipoles.  It then calculates the atomic multipoles from the Hirshfeld and 
  Hirshfeld-Iterated atomic charge densities.  It also calculates the molecular multipoles from the Hirshfeld and 
  Hirshfeld-Iterated atomic multipoles and compares to the QM molecular multipoles
*/


