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

#ifndef GAUSS_DENS_STRUCTURE_H
#define GAUSS_DENS_STRUCTURE_H

struct Atom_Data
/*contains information on the molecule obtained from Gaussian .fchk files.
  atomic units (bohr) are used.
*/
{
	int n_atom;
	int n_electron;
	double **r;		/*x,y,z of nuclie*/
	double *Z;		/*nuclear charge*/
	int charge;
	int multi;
};


struct Atomic_Orbital_Data
/*stores information on the atomic orbitals and QM basis sets.  atomic units (bohr) are used.*/
{
	int atom_num;
	double *r;		/*x,y,z of center*/
	int n_prim;
	double *exp;	/*[n_prim]*/
	double *coef;	/*[n_prim]*/
	int type;		/*0 = s, 1 = p, 2 = d, 3 = f, 4 = g*/

	int n_ang;				/*1 for s, 3 for p, 6 for d, 10 for f*/
	char **ang_descriptor;	/*[n_ang]*/
	int **nxyz;				/*[n_ang][3]*/
	/*nx[i] = nxyz[i][0], ny[i] = nxyz[i][1], nz[i] = nxyz[i][2]*/
	double **normal_const;	/*[n_prim][n_ang]*/

	int ***nxyz_2_ang_mom_num;	/*gives ang_mom_num for [nx][ny][nz]*/
};


struct Wave_Function_Data
/*stores molecular orbitals and one electron density matrix.  the molecular orbtials are tested by 
  calculating the one electron density SCF_P[n_basis][n_basis] and comparing to Gaussian's reference value.  
  In addition, correlated one electron density matrices (e.g. MP2 or CCSD) (if present in the .fchk file)
  are read in and stored in P[n_basis][n_basis]. 
*/
{
	int n_basis;	/*sum of n_ang over atomic orbitals*/
	int n_alpha;	/*# of spin up electrons*/
	int n_beta;		/*# of spin down electrons*/

	double **SCF_P;				/*[n_basis][n_basis] SCF density matrix*/
	double **P;					/*[n_basis][n_basis] total density matrix*/
	int n_orb;

	int IS_CLOSED_SHELL;
	double *alpha_MO_energy;		/*[n_orb]*/
	double **alpha_MO_Coef;			/*[n_basis][n_orb]*/

	double *beta_MO_energy;			/*[n_orb]*/
	double **beta_MO_Coef;			/*[n_basis][n_orb]*/

	int *ao_num;	/*[n_basis]*/
	int *ang_num;	/*[n_basis]*/

	int **ao_ang_num_2_basis_num;	/*[n_AO][n_ang] -> gives basis number*/
};

#endif


