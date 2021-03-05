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

/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/*******************PARAMETER READING FUNCTIONS*****************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/

void Read_File_List (char *path, int *n_files, char ***filenames, char *filelist, int DEBUG, FILE *fp_debug);
/*
/home/denny/temp
	num_files 3
		temp1.dat
		temp2.dat
		temp3.dat
*/			


void Read_Atom_Types (int *n_atom_types, char ***atom_types, char *filename, 
	char *flag, FILE *fp_debug, int DEBUG);
/*
	num_atom_types 2
		H  
		HO 
*/

void Read_AT_Double (double **real_value, int n_atom_types, char **atom_types, 
	char *filename, char *flag, FILE *fp_debug, int DEBUG);
/*
	n_polar 2
		H  	0.21
		HO 	0.25
*/

void Read_AT2_Double (int *n_real, double **real_value, int ***ab_num_2_atom_type, int n_atom_types, 
	char **atom_types, char *filename, char *flag, FILE *fp_debug, int DEBUG);
/*
	n_parm 2
		H  	HO	0.21
		HO 	CT	0.25
*/

void Read_AT_Double_2 (double **real_value1, double **real_value2, int n_atom_types, 
	char **atom_types, char *filename, char *flag, FILE *fp_debug, int DEBUG);
/*
	n_vdw 2
		H  	1.21 0.01
		HO 	0.50 0.01
*/

void Read_AT2_Double2 (int *n_real, double **real_value1, double **real_value2, 
	int ***ab_num_2_atom_type, int n_atom_types, char **atom_types, char *filename, 
	char *flag, FILE *fp_debug, int DEBUG);
/*
	n_parm 3
		HC CT	0.00	0.00
		OH HO	2.22	0.0448
		OH OH	3.44 	0.2104
*/

void Read_AT_String (int n_atom_types, char **atom_type_names, 
		char ***string, char *filename, char *flag, int DEBUG, FILE *fp_debug);
/*
@Elements
	num_elements 18
		H		H
		HO		H
		HC		H
*/

void Read_Opt_AT_Parameters (int *n_opt, int **opt_num_2_at_num, int n_atom_types, 
	char **atom_types, char *filename, char *flag, FILE *fp_debug, int DEBUG);
/*
	num_atom_types 2
		H  
		HO 
*/

void Read_Opt_AT_Parameters2 (int *n_opt2, int **opt_num_2_ab_num, int n_ab,
	int **ab_num_2_atom_type, int n_atom_types, char **atom_types, char *filename, 
	char *flag, FILE *fp_debug, int DEBUG);
/*
	num_atom_types 2
		H  	H
		HO 	CT
*/

void Read_Bool_String_Pair ( int **is_value, char ***s_type, int n_atom_types, 
	char **atom_type_names, char *filename, char *flag, FILE *fp_debug, int DEBUG);
/*
@EP_Definition
	num_ep 4
		N2_EP 3_sp3
		OH_EP 2_sp3
		CD_EP 3_sp2
		O_EP  1_sp2
*/

/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/*******************FILE READING FUNCTIONS**********************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/

void Read_ml2_File_W_CE (int *num_atoms, int **atom_type, double **r_1D_atoms, 
		double ***r_2D_atoms, double **charges, int **n_neighbors, int ***neighbors,
        char ***element, int num_atom_types, char **atom_type_names, 
		char *filename, FILE *fp_debug, int DEBUG);
/*Reads in (modified) ml2 file with charges and elements*/

void Read_Conn_Files (int *n_atoms, int **atom_types, double **charges,
		double **polar, int **n_neighbors, int ***neighbors, int n_atom_types, 
		char **atom_type_names, int s_atom_type_size, char *filename, 
		FILE *fp_debug, int DEBUG);
/*Reads Connectivity File in Following Format :
n_atoms = 6
           #  at   chg  n_neigh    neigh
           0   CT    0.24417    0.73630      4             1    2    3    4
           1   H1    0.00000    0.16570      1             0
*/

void Read_EG_File ( int *n_monomer, int **n_atom, double ****R2D_monomer,
	char ***conn_files, int *n_clusters, double *****R2D_cluster, double **QM_Energy,
	double *****QM_Force, char *filename, int DEBUG, FILE *fp_debug);


void Read_XYZ_File (char *filename, char *flag, int *n_atoms, int **mol_atom_types, 
		double ***mol_r2d, int **n_neighbors, int ***neighbors, int n_atom_types, 
		char **atom_type_names, FILE *fp_debug, int DEBUG);
/*
flag
     8
     1 CT  C   -1.594800    0.078100    0.000000     4     2     3     4     5 
     2 CT  C   -0.060500    0.215100   -0.000000     4     1     6     7     8 
*/

void Read_Pdb_File (int *n_atoms, double ***R_2D, char ***atom_type, 
		char *filename, FILE *fp_debug, int DEBUG);
/*
ATOM      1  H1  OME     1       5.322  16.919  10.213
ATOM      2  CH3 OME     1       5.853  15.985  10.410
*/

void Read_Monomer_Cluster_File (char *filename, int *n_monomers, int **n_atoms, 
		char ****element, double ****r2d, char ***connfilename, 
		int **charge, int **multi, int DEBUG, FILE *fp_debug);
/*
@Monomers
	num_molecules 2
		molecule 0	
			num_atoms 5	charge 0	mult. 1		10_form_acid.conn
 				C                 -0.72903000   -0.21501100    0.02393300
 				O                 -0.16249000    0.85642000   -0.00391100
 				O                 -2.06088100   -0.37061600   -0.00185100
 				H                 -0.22075400   -1.18584000    0.07380900
 				H                 -2.46657300    0.51437300   -0.04824200
			
		molecule 1
			num_atoms 3	charge 0	mult. 1		10_water.conn
  				H		     1.76695800    0.36150500	-0.02143800
 				O		     2.47862100   -0.30158700	-0.05407600
  				H		     3.25255000    0.12628900	 0.33098100
*/
void Read_AT2_XYZ_File (char *filename, int *n_atoms, int **atom_types1, int **atom_types2,
		double ***mol_r2d, int **n_neighbors, int ***neighbors, int n_atom_types, 
		char **atom_type_names, FILE *fp_debug, int DEBUG);

/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/*******************PRINTING FILE FUNCTIONS*********************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/

void Print_Conn_File (int n_atoms, int *element_num, char **AT_Elements, 
	double *charge, double *polar, int *n_neighbors, int **neighbors, char *filename);
/*
n_atoms = 6
		    #  at   chg  n_neigh    neigh
           0   CT    0.24417    0.73630      4             1    2    3    4
           1   H1    0.00000    0.16570      1             0
           2   H1    0.00000    0.16570      1             0
           3   HA    0.00000    0.16570      1             0
           4   OH   -0.62388    0.41730      2             0    5
           5   HO    0.37971    0.11040      1             4
*/

void Print_XYZ_Format (int num_atoms, double **r_2D_atoms, 
		int *n_neighbors, int **neighbors, int *atom_type,
		int n_atom_types, char **atom_type_names, char *filename, 
		FILE *fp_debug, int DEBUG);
/*
  1112
     0 CT     -2.852565    7.129814    7.213353     1     	1     2    3    4
     1 HC     -2.010039    6.775897    7.848916     3     	0
     2 HC     -3.744861    7.281564    7.792989     4     	0
*/
void Make_Pdb_File (int num_atoms, char **element, double **r_2d, char *Filename);

void Write_Coordinates (FILE *fp_crd, int num_atoms, double *r_crd);

void Print_Geometry (char *title, FILE *fp, int n_atom, char **atom_names,
	double **r2d, int *n_neigh, int **neigh);
/*
@Geometry
	n_atoms 4
		H1		0.28827   0.0000       -2.0037 		1	1
		O2		0.87274   0.0000       -1.24675		2	0 2
		H3		0.28827   0.0000       -0.4898 		1	1
		X		0.28827   0.0000       -1.24675		0
*/
void Print_Double (char *title, FILE *fp, int n_value, char **atom_names, double *Y);
/*
@Z
	n_atoms 4
		H1	1.0
		O2	6.0
		H3	1.0
		X	0.0
*/

/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/*******************MISC FUNCTIONS******************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
void IO_Error (char *filename);
FILE *Open_File (char *filename, char *action);
int string_compare(char *s1, char *s2, int N);
int Is_New_Line(char *str);
void Get_Next_Line (FILE *fp, char *str, char *filename);
void Go_To_Line (FILE *fp, char *flag, char *filename);
/*skips to line where the flag char string is found*/
void Read_Value (FILE *fp, char *flag, char *filename, char *scan_str, char *return_str);
/*looks for: starting at flag,  if flag == NULL, starts at begining of file.
scan_str return_str
if it finds scan_str, then it returns return_str, if not return_str = NULL
*/
int Find_Atom_Type (char *atom_type_s, int Num_atom_types, char **Atom_Type_Names);
void Get_Current_Directory (char *directory);
void Change_File_Extension (const char *s1, const char *s2, char *s3);
/*attaches s2 to s1 beginning with at the '.' in s1, returns s3*/
void Remove_Extension (char *s1, char *s2);
/*
        s1 = test.dat => s2 = test
*/
int Does_File_Exists (char *filename);
void Create_Sub_Directory (char *directory_name);

void Get_Extension (char *s1, char *s2);
/*
        s1 = test.dat => s2 = dat
*/


