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

#ifndef STRUCTURE_H
#define STRUCTURE_H


struct ELE_Data 
/*contains information on the elements such as names and nuclear charges.  also contains number of ionic states,
  charges of ionic states, spin multiplicities of ionic states, and number of electrons of the ionic state.
  for each element, it contains a radial grid with quadrature weight.  it also contains the ionic (radially 
  symmetric) charge densities.  Angstroms are used.
*/
{ 
/*Number of elements, element names, and nuclear charges*/
	int n_element;
	char **element;			/*[n_element] = name of element*/
	double *Z;				/*[n_element] = (double) nuclear charge of element*/
	int *Z_i;				/*[n_element] = (int)    nuclear charge of element*/

/*Number of ionic states (per element), ionic charges, ionic spin multiplicities, and number of electrons*/
	int *n_ion_state;		/*[n_element]*/
	int **charge;			/*[n_element][n_ion_state]*/
	int **multi;			/*[n_element][n_ion_state]*/
	int **n_electron_i;		/*[n_element][n_ion_state]*/

/*Radial grids for each element with radius and radius weight*/
	int *N_radial_pt;		/*[n_element]*/
	double **radius;		/*[n_element][N_radial_pt]*/
	double **radius_weight;	/*[n_element][N_radial_pt]*/

/*Spherically symmetric atomic ion charge densities for each element and ionic state*/
	double ***density;		/*[n_element][n_ion_state][N_radial_pt]*/



};

struct MOL_Data 
/*contains information on the molecule of interest, such as number of atoms, geometry, elements, and element_num which 
  indexes the element number in the ELE_Data structure.  It contains the molecular spin multiplicity, molecular charge,
  Gaussion 03/09 I/O options such as theory (e.g. MP2/aug-cc-pVTZ), nproc, mem, and DO_CHELPG (0 or 1).  The Hirshfeld,
  atomic charges, Hirshfeld-Iterated atomic charges, number of electrons, and ChelpG charges are stored for each atom.
  The QM 'pure' cartesian molecular multipoles are read in from Gaussian 03/09 output and stored in QM_cart_multi.  The
  'pure' cartesian molecular multipoles are converted to spherical tensor molecular multipoles and then to traceless Cartesian
  molecular multipoles.  The Hirshfeld and Hirshfeld-Iterated atomic multipoles are stored here.  Molecular multipoles
  are calculated from the atomic multipoles and compared with their reference QM values.  Angstroms are used.
*/
{ 
	int n_atom;
	double **r2d;			/*[n_atom][3]*/
	char **element;			/*[n_atom]*/
	int *element_num;		/*[n_atom]*/

	int multi;
	int charge;
	char theory[500];
	int nproc;
	char memory[500];
	char inputfile[1000];
	int DO_CHELPG;


	double *HI_charge;		/*[n_atom]*/
	double *H_charge;		/*[n_atom]*/
	double *n_electron;		/*[n_atom]*/
	double *atom_charge;

	double *chelpg_charge;	/*[n_atom]*/

/*QM molecular multipoles obtained from Gaussian 03/09*/
	double ***QM_cart_multi;	/*[lx][ly][lz], lx+ly+lz <= 4 (for molecular hexadecapole)*/
	double ***QM_cart_multi_TR;	/*[lx][ly][lz], lx+ly+lz <= 4 (for molecular hexadecapole)*/
	double **Q_QM_mol_lm_r;		/*[lmax][2l+1], lmax = 4*/
	double **Q_QM_mol_lm_i;		/*[lmax][2l+1], lmax = 4*/


	double ***cart_TR_atom;		/*temporary array for converting spherical tensor atomic multipoles to traceless cartesian atomic multipoles*/

/*Complex spherical tensor atomic multipoles from Hirshfeld atomic densities*/
	int lmax_atom;
	double ***Qa_H_lm_r;		/*[n_atom][lmax][2l+1]*/
	double ***Qa_H_lm_i;		/*[n_atom][lmax][2l+1]*/

/*Complex spherical tensor atomic multipoles from Hirshfeld-Iterated atomic densities*/
	double ***Qa_HI_lm_r;		/*[n_atom][lmax][2l+1]*/
	double ***Qa_HI_lm_i;		/*[n_atom][lmax][2l+1]*/

/*temporary arrays*/
	double **Q_temp_lm_r;	/*[lmax][2l+1]*/
	double **Q_temp_lm_i;	/*[lmax][2l+1]*/
	double **C_lm_r;
	double **C_lm_i;

/*Molecular multipoles calculated by Hirshfeld atomic multipoles.  
  {Q_H_mol_lm_r[0][][], Q_H_mol_lm_i[0][][] = are the molecular multipoles calculated by Hirshfeld atomic charges.
  {Q_H_mol_lm_r[1][][], Q_H_mol_lm_i[1][][] = are the molecular multipoles calculated by Hirshfeld atomic dipoles.
  {Q_H_mol_lm_r[2][][], Q_H_mol_lm_i[2][][] = are the molecular multipoles calculated by Hirshfeld atomic quadrupoles.
  etc.
*/
	double ***Q_H_mol_lm_r;	/*[lmax] = 0,1,2.. for atomic monopoles, dipoles, etc.   [lmax][2l+1]*/
	double ***Q_H_mol_lm_i;	/*[lmax][lmax][2l+1]*/

/*Molecular multipoles calculated by Hirshfeld-Iterated atomic multipoles.  
  {Q_H_mol_lm_r[0][][], Q_H_mol_lm_i[0][][] = are the molecular multipoles calculated by Hirshfeld atomic charges.
  {Q_H_mol_lm_r[1][][], Q_H_mol_lm_i[1][][] = are the molecular multipoles calculated by Hirshfeld atomic dipoles.
  {Q_H_mol_lm_r[2][][], Q_H_mol_lm_i[2][][] = are the molecular multipoles calculated by Hirshfeld atomic quadrupoles.
  etc.
*/
	double ***Q_HI_mol_lm_r;	/*[lmax] = 0,1,2.. for atomic monopoles, dipoles, etc.   [lmax][2l+1]*/
	double ***Q_HI_mol_lm_i;	/*[lmax][lmax][2l+1]*/

/*Traceless Cartesian molecular multipoles calculated by Hirshfeld atomic multipoles
   MM_H_cart_multi[0][lx][ly][lz] = molecular multipoles calculated by Hirshfeld atomic charges for lx+ly+lz <= 4
   MM_H_cart_multi[1][lx][ly][lz] = molecular multipoles calculated by Hirshfeld atomic dipoles for lx+ly+lz <= 4
   MM_H_cart_multi[2][lx][ly][lz] = molecular multipoles calculated by Hirshfeld atomic quadrupoles for lx+ly+lz <= 4
   etc.
*/
	double ****MM_H_cart_multi;	

/*Traceless Cartesian molecular multipoles calculated by Hirshfeld-Iterated atomic multipoles
   MM_HI_cart_multi[0][lx][ly][lz] = molecular multipoles calculated by Hirshfeld atomic charges for lx+ly+lz <= 4
   MM_HI_cart_multi[1][lx][ly][lz] = molecular multipoles calculated by Hirshfeld atomic dipoles for lx+ly+lz <= 4
   MM_HI_cart_multi[2][lx][ly][lz] = molecular multipoles calculated by Hirshfeld atomic quadrupoles for lx+ly+lz <= 4
   etc.
*/
	double ****MM_HI_cart_multi;	

};

struct GRID_HI_Data 
/*contains molecular grid, the QM molecular charge density on the grid, the proto atomic charge densities on the 
  molecular grid, and the partitioned atomic charge densities on the grid.  Angstroms are used.
*/
{ 
/*molecular grid*/
	int n_grid_pts;	
	double **r2d;	/*[n_grid_pts][3]*/
	double *weight;	/*[n_grid_pts]*/

/*QM molecular charge density*/
	double *p_QM;			/*[n_grid_pts]*/

/*proto atomic charge densities for Hirshfeld and Hirshfeld-Iterated models*/
	double *p0_H_tot;		/*[n_grid_pts]*/
	double *p0_HI_tot;		/*[n_grid_pts]*/

/*total proto atomic charge density from either Hirshfeld or Hirshfeld-Iterated model*/
	double *p0_tot;			/*[n_grid_pts]*/

/*QM atomic ion charge densities on the molecular grid*/
	double ***p_QM_atom;	/*[n_atom][n_ion_state][n_grid_pts]*/

/*present partioned atomic charge density*/
	double **p_MM_atom;		/*[n_atom][n_grid_pts]*/

/*saved Hirshfeld partitioned atomic charge density*/
	double **p_MM_atom_H;	/*[n_atom][n_grid_pts]*/

/*saved Hirshfeld-Iterated partitioned atomic charge density*/
	double **p_MM_atom_HI;	/*[n_atom][n_grid_pts]*/
};

#endif


