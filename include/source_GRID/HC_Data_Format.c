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


#define RAD_TOL 1.0E-6
#define MAX_ITERATION  500
#define eA_2_DEBYE     (2.541765/0.529177249)



void Read_Geometry (char *filename, int *n_atom, double ***r2d, char ***element,
	int *multi, int *charge, char *theory, int *nproc, char *memory, int *lmax, int *DO_CHELPG,
	FILE *fp_debug, int DEBUG);
void Read_Atom_Charge_Info (char *filename, struct ELE_Data *ele, FILE *fp_debug, int DEBUG);
void Check_Keyword (char *str1, char *str2, char *filename);
void Calc_Molecular_Dipole (double *dip, int n_atom, double **r2d, double *charge);

/*Gaussian IO*/
void Create_Gauss_Input_CHARGE (char *com_file, int n_atoms, char **element, 
	double **r2d, int charge, int multi, char **line, int max_lines, 
	char **keyword, int max_keywords);
int Check_Complete_Log_File (char *log_file);
void Read_Gauss_Multipoles (char *filename, double ***cart, FILE *fp_debug, int DEBUG);
void Read_ChelpG_Charges (char *filename, int n_atom, double **chelpg_charge, int DEBUG, FILE *fp_debug);
double ***Allocate_Cart_Tensor (int lmax);
void Nxyz_Index_2_String (char *str, int n1, int n2, int n3);
void Get_Atom_File_Name (char *initial, char *element, int charge);

double **Allocate_Spherical_Array (int l);
/*Allocates array spherical array of the form:
A[0][0]
A[1][0], A[1][1], A[1][-1]
A[2][0], A[2][1], A[2][2], A[2][-1], A[2][-2]
..

A[l][0], A[l][1],... A[l][l], A[l][-1], .. A[l][-l]
*/
void Convert_Spherical_Multipole_2_TR_Cart_Multipole (int lmax, double ***CART_multi, double **Q_lm_r, 
	double **Q_lm_i, int DEBUG, FILE *fp_debug);
/*Converts Spherical Multipoles into a Traceless Cartesian multipole.  Does not output CART_multi in Buckingham convention, 
  to convert CART_multi to the Buckingham convention, multiply by (2l-1)!!*/


void Read_Input_File (char *filename, struct ELE_Data *ele, struct MOL_Data *mol, int DEBUG, FILE *fp_debug)
/*Reads INPUT_FILE and called from main().  This function reads the geometry of the molecule of interest (mol data).  
  It then reads elements, atomic nuclear charges (Z), and ionic states with multiplicities (ele data).  For each atom (element) in
  the molecule of interest, it indexes the element number for that atom.
*/
{
	int i;
	double x;

/*Reads geometry of the molecule of interest*/
	Read_Geometry (filename, &( (*mol).n_atom), &( (*mol).r2d), &( (*mol).element), 
		&( (*mol).multi), &( (*mol).charge), (*mol).theory, &( (*mol).nproc), (*mol).memory, &( (*mol).lmax_atom), 
		&( (*mol).DO_CHELPG), fp_debug, DEBUG);

/*Reads element definitions*/
	Read_Atom_Types ( &( (*ele).n_element), &( (*ele).element), filename, "@Element", fp_debug, DEBUG);

/*Reads nuclear charge definitions as a real (double) number*/
	Read_AT_Double (&( (*ele).Z), (*ele).n_element, (*ele).element, 
		filename, "@Z", fp_debug, DEBUG);

/*Also stores nuclear charges as integers*/
	(*ele).Z_i = I_Allocate_1D_Matrix ((*ele).n_element);
	for (i = 0; i < (*ele).n_element; i++)
	{
		(*ele).Z_i[i] = (int) ((*ele).Z[i] + 1.0E-6);
		x = (double) (*ele).Z_i[i];
		
		if (fabs(x - (*ele).Z[i]) > 1.0E-8)
		{
			printf("ERROR in (double) Z to (int) Z for %s    Z = %20.15f Z_i = %d\n", (*ele).element[i],
				(*ele).Z[i], (*ele).Z_i[i]);
			exit(0);
		}
	}

/*Reads atomic ion charges with multiplicities*/
/*	Read_Atom_Charge_Info (filename, ele, fp_debug, DEBUG); */

/*indexes the element number each atom for the molecule of interest*/
	(*mol).element_num = I_Allocate_1D_Matrix ((*mol).n_atom);
	for (i = 0; i < (*mol).n_atom; i++)
		(*mol).element_num[i] = Find_Atom_Type ((*mol).element[i], (*ele).n_element, (*ele).element);


}

void Read_Atom_Charge_Info (char *filename, struct ELE_Data *ele, FILE *fp_debug, int DEBUG)
/*Reads atomic ion charges with multiplicities*/
{
	FILE *fp;
	char str[2000];
	int i, j, i_temp, i_temp2, at_num;
	char s_temp[500];
	char s_temp2[500];
	char s_temp3[500];

	fp = Open_File (filename, "r");
	Go_To_Line (fp, "@Ion_States", filename);

	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &i_temp);
	if (i_temp != (*ele).n_element)
	{
		printf("ERROR reading @Ion_States from %s.. n_element do not match\n", filename);
		exit(0);
	}

	(*ele).n_ion_state  = (int *)  malloc ( (*ele).n_element*sizeof(int) );
	(*ele).charge       = (int **) malloc ( (*ele).n_element*sizeof(int *) );
	(*ele).multi        = (int **) malloc ( (*ele).n_element*sizeof(int *) );
	(*ele).n_electron_i = (int **) malloc ( (*ele).n_element*sizeof(int *) );


	for (i = 0; i < (*ele).n_element; i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %s %d %s %d", s_temp, s_temp2, &i_temp, s_temp3, &(i_temp2) );
		Check_Keyword (s_temp2, "multi=",       filename);
		Check_Keyword (s_temp3, "n_ion_state=", filename);

		at_num = Find_Atom_Type (s_temp, (*ele).n_element, (*ele).element);

		/*count neutral atom too*/
		(*ele).n_ion_state[at_num]  = i_temp2+1;
		(*ele).charge[at_num]       = I_Allocate_1D_Matrix ((*ele).n_ion_state[at_num]);
		(*ele).multi[at_num]        = I_Allocate_1D_Matrix ((*ele).n_ion_state[at_num]);
		(*ele).n_electron_i[at_num] = I_Allocate_1D_Matrix ((*ele).n_ion_state[at_num]);

		(*ele).charge[at_num][0]       = 0;
		(*ele).multi[at_num][0]        = i_temp;
		(*ele).n_electron_i[at_num][0] = (*ele).Z_i[at_num] - (*ele).charge[at_num][0];

		for (j = 1; j < (*ele).n_ion_state[at_num]; j++)
		{
			Get_Next_Line (fp, str, filename);
			sscanf(str, "%s %d %s %d", s_temp2, &((*ele).charge[at_num][j]), s_temp3, &((*ele).multi[at_num][j]) );
			(*ele).n_electron_i[at_num][j] = (*ele).Z_i[at_num] - (*ele).charge[at_num][j];
			Check_Keyword (s_temp2, "charge=", filename);
			Check_Keyword (s_temp3, "multi=",  filename);
		}
	}


	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n@Ion_States from %s\n", filename);
		fprintf(fp_debug, "\tn_element %d\n", (*ele).n_element);
		for (i = 0; i < (*ele).n_element; i++)
		{
			fprintf(fp_debug, "\t\t%4s charge= %d n_electron = %d  Z_i = %d multi= %d n_ion_state= %d\n", (*ele).element[i], (*ele).charge[i][0], 
				(*ele).n_electron_i[i][0], (*ele).Z_i[i], (*ele).multi[i][0], (*ele).n_ion_state[i]);
			for (j = 1; j < (*ele).n_ion_state[i]; j++)
				fprintf(fp_debug, "\t\t\tcharge= %d multi= %d n_electron = %d  Z_i = %d\n", (*ele).charge[i][j], (*ele).multi[i][j],
					(*ele).n_electron_i[i][j], (*ele).Z_i[i]);
		}
	}


}

void Check_Keyword (char *str1, char *str2, char *filename)
{
	if (strcmp(str1, str2) != 0)
	{
		printf("ERROR reading %s.. expecting %s found %s\n", filename, str2, str1);
		exit(0);
	}
}



void Read_Geometry (char *filename, int *n_atom, double ***r2d, char ***element,
	int *multi, int *charge, char *theory, int *nproc, char *memory, int *lmax, int *DO_CHELPG,
	FILE *fp_debug, int DEBUG)
/*Reads molecular geometry, molecular spin multiplicity, molecular charge, and also IO options for Gaussian (nproc, memory, DO_CHELPG)*/
{
	int i, j;
	char str[2000];
	FILE *fp;
	int i_temp[5];
	char s_temp[500];
	char s_temp2[500];
	
	fp = Open_File (filename, "r");


	Go_To_Line (fp, "@Geometry", filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &( (*n_atom) ) );
	(*r2d)     = D_Allocate_2D_Matrix ((*n_atom), 3);
	(*element) = C_Allocate_2D_Matrix ((*n_atom), 10);



	for (i = 0; i < (*n_atom); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %lf %lf %lf", (*element)[i], 
			&( (*r2d)[i][0]), &( (*r2d)[i][1]), &( (*r2d)[i][2]) );
	}


	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nReading %s\n", filename);
		fprintf(fp_debug, "\tn_atoms = %d\n", (*n_atom) );
		for (i = 0; i < (*n_atom); i++)
		{
			fprintf(fp_debug, "\t\t%4s %10.6f %10.6f %10.6f\n", (*element)[i],
				(*r2d)[i][0], (*r2d)[i][1], (*r2d)[i][2]);
		}
		fflush(fp_debug);
	}
}










/**************************************************************/
/**************************************************************/
/******************GAUSSIAN IO*********************************/
/**************************************************************/
/**************************************************************/

void Read_Gauss_Input (char *initial, struct ELE_Data *ele, struct MOL_Data *mol, int DEBUG, FILE *fp_debug)
/*This function will first check to see if all necessary Gaussian 03 or 09 output files (.log and .fchk) are available in the working directory.
  If the Gaussian output files are available, this function will read in the ChelpG charges (if (*mol).DO_CHELPG == 1) and the molecular multipoles.
  If the Gaussian output files are not available, this function will create the necessary input files (.com) for Gaussian 03 or 09 for
  the molecule of interest and the atomic ions.  This function will then exit and wait for the user to execute Gaussian on all of the created .com files.
  The user will also need to convert the checkpoint (.chk) files created by Gaussian into formatted checkpoint (.fchk) files.
*/
{
	int i, j;

	char com_file[1000];
	char log_file[1000];
	char chk_file[1000];
	char fchk_file[1000];
	int n_keyword;
	char **keyword;
	int n_line;
	char **line;
	int DONE;
	double **r2d_atom;
	char file_initial[1000];


	DONE = 1;

	sprintf(com_file, "%s.com", initial);
	sprintf(log_file, "%s.log", initial);

	n_keyword = 5;
	keyword   = C_Allocate_2D_Matrix (n_keyword, 200);
	sprintf(keyword[0], "%s", (*mol).theory);
	sprintf(keyword[1], "scf=tight");
	sprintf(keyword[2], "Nosymm");
	sprintf(keyword[3], "density=current");
	sprintf(keyword[4], "pop=chelpg");
	if ( (*mol).DO_CHELPG == 0)
		n_keyword = 4;

	n_line = 3;
	line   = C_Allocate_2D_Matrix (n_line, 200);
	sprintf(line[0], "%%chk=%s.chk", initial);
	if ( (*mol).nproc != -1)
	{
		sprintf(line[1], "%%NPROC=%d", (*mol).nproc);
		sprintf(line[2], "%%MEM=%s", (*mol).memory);
	}
	else
		n_line = 1;


	if (Check_Complete_Log_File (log_file) == 0)
	{
		printf("Creating %s\n", com_file);
		Create_Gauss_Input_CHARGE (com_file, (*mol).n_atom, (*mol).element, 
			(*mol).r2d, (*mol).charge, (*mol).multi, line, n_line, 
			keyword, n_keyword);
		DONE = 0;
	}
	else
	{
		(*mol).QM_cart_multi = Allocate_Cart_Tensor (4);
		Read_Gauss_Multipoles (log_file, (*mol).QM_cart_multi, fp_debug, DEBUG);
		if ( (*mol).DO_CHELPG)
			Read_ChelpG_Charges (log_file, (*mol).n_atom, &( (*mol).chelpg_charge), DEBUG, fp_debug);
		else
			(*mol).chelpg_charge = D_Allocate_1D_Matrix ((*mol).n_atom);

		sprintf(chk_file, "%s.fchk", initial);
		if (Does_File_Exists (chk_file) == 0)
		{
			printf("Waiting for %s.fchk.. try\n", initial);
			printf("formchk %s.chk\n", initial);
			DONE = 0;
		}
	}
	
/*Get Atomic .fchk files*/
	r2d_atom = D_Allocate_2D_Matrix (1, 3);
	r2d_atom[0][0] = 0.0;		r2d_atom[0][1] = 0.0;		r2d_atom[0][2] = 0.0;


	n_keyword = 4;
	keyword   = C_Allocate_2D_Matrix (n_keyword, 200);
	sprintf(keyword[0], "%s", (*mol).theory);
	sprintf(keyword[1], "scf=(qc,tight)");
	sprintf(keyword[2], "Nosymm");
	sprintf(keyword[3], "density=current");


	for (i = 0; i < (*ele).n_element; i++)
		for (j = 0; j < (*ele).n_ion_state[i]; j++)
			if ((*ele).multi[i][j] >= 0)
			{
				Get_Atom_File_Name (file_initial, (*ele).element[i], (*ele).charge[i][j]);
				sprintf(com_file,  "%s.com",  file_initial);
				sprintf(log_file,  "%s.log",  file_initial);
				sprintf(chk_file,  "%s.chk",  file_initial);
				sprintf(fchk_file, "%s.fchk", file_initial);

				sprintf(line[0],  "%%chk=%s", chk_file);
				if (Check_Complete_Log_File (log_file) == 0)
				{
					printf("Creating %s\n", com_file);
					Create_Gauss_Input_CHARGE (com_file, 1, &((*ele).element[i]), 
						r2d_atom, (*ele).charge[i][j], (*ele).multi[i][j], line, n_line, 
						keyword, n_keyword);
					DONE = 0;
				}
				else
				{
					if (Does_File_Exists (fchk_file) == 0)
					{
						printf("Waiting for %s.. try\n", fchk_file);
						printf("formchk %s\n", chk_file);
						DONE = 0;
					}
				}
			}


	if (DONE == 0)
	{
		printf("waiting for Gaussian calculations to finish\n");
		if (DEBUG)	fclose(fp_debug);
		exit(0);
	}

}

void Get_Atom_File_Name (char *initial, char *element, int charge)
{
	if (charge < 0)
		sprintf(initial, "%s_neg_%d", element, -charge);
	else
		sprintf(initial, "%s_pos_%d", element, charge);


}

double ***Allocate_Cart_Tensor (int lmax)
/*Allocates Cart[l1][l2][l3], l1+l2+l3 <= lmax*/
{
	double ***A;
	int l1, l2;

	A = (double ***) malloc ( (lmax+1)*sizeof(double **) );
	for (l1 = 0; l1 <= lmax; l1++)
	{
		A[l1] = (double **) malloc ( (lmax-l1 + 1)*sizeof(double *) );
		for (l2 = 0; l2 <= lmax-l1; l2++)
			A[l1][l2] = (double *) malloc ( (lmax-l1-l2 + 1)*sizeof(double ) );
	}
	return A;
}

void Create_Gauss_Input_CHARGE (char *com_file, int n_atoms, char **element, 
	double **r2d, int charge, int multi, char **line, int max_lines, 
	char **keyword, int max_keywords)
/*Creates .com input file for gaussian*/
{
	int i;
	FILE *fp;
	char input_line[1000];
	char str[2000];
	char s_temp[500];

	fp = Open_File (com_file, "w");

/*Insert command lines.. eg. nproc=2, mem=800MB.. etc*/
	for (i = 0; i < max_lines; i++)
		if (line[i][0] != '\0')
			fprintf(fp, "%s\n", line[i]);


/*keyword line*/
	sprintf(input_line, "# ");
	for (i = 0; i < max_keywords; i++)
		if (keyword[i][0] != '\0')
		{
			strncat(input_line, keyword[i], 200);
			strncat(input_line, " ", 200);
		}
	fprintf(fp, "%s\n", input_line);

/*blank line and then title*/
	fprintf(fp, "\ntitle\n\n");

/*charge and multiplicity*/
	fprintf(fp, "%d %d\n", charge, multi);

/*elements and coordinates*/
	for (i = 0; i < n_atoms; i++)
		fprintf(fp, "%s %12.9f %12.9f %12.9f\n", element[i],
			r2d[i][0], r2d[i][1], r2d[i][2]);
	fprintf(fp, "\n");

	fprintf(fp, "\n\n\n");
	fclose(fp);
}

int Check_Complete_Log_File (char *log_file)
/*looks for :
 Normal termination of Gaussian
*/
{
	char str[2000];
	char str1[200], str2[200], str3[200], str4[200];
	FILE *fp;

	fp = fopen( log_file, "r");
	if (fp == NULL) 
		return 0;

	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s %s %s %s", str1, str2, str3, str4);
			if (strcmp(str1, "Normal") == 0)
				if (strcmp(str2, "termination") == 0)
					if (strcmp(str3, "of") == 0)
						if (strcmp(str4, "Gaussian") == 0)
						{
							fclose(fp);
							return 1;
						}
		}
		else
		{
			fclose(fp);
			return 0;
		}
	}
}

void Read_Gauss_Multipoles (char *filename, double ***cart, FILE *fp_debug, int DEBUG)
/*Reads in charge, dipole, quadrupole, octapole, and hexadecapole in D*A^(n-1)
without converting to traceless forms
 Charge=     0.000000000 electrons
 Dipole moment (field-independent basis, Debye):
    X=     0.000697643    Y=     0.000000000    Z=    -2.198893682  Tot=     2.198893793
 Quadrupole moment (field-independent basis, Debye-Ang):
   XX=    -7.205697372   YY=    -4.105517149   ZZ=    -6.000962372
   XY=     0.000000000   XZ=     0.002362056   YZ=     0.000000000
 Traceless Quadrupole moment (field-independent basis, Debye-Ang):
   XX=    -1.434971741   YY=     1.665208482   ZZ=    -0.230236741
   XY=     0.000000000   XZ=     0.002362056   YZ=     0.000000000
 Octapole moment (field-independent basis, Debye-Ang**2):
  XXX=     0.000356476  YYY=     0.000000000  ZZZ=    -1.430326239  XYY=    -0.000060274
  XXY=     0.000000000  XXZ=    -0.387281484  XZZ=     0.000563697  YZZ=     0.000000000
  YYZ=    -1.354783518  XYZ=     0.000000000
 Hexadecapole moment (field-independent basis, Debye-Ang**3):
 XXXX=    -5.184075654 YYYY=    -5.363435711 ZZZZ=    -5.989339181 XXXY=     0.000000000
 XXXZ=     0.001969478 YYYX=     0.000000000 YYYZ=     0.000000000 ZZZX=     0.002041866
 ZZZY=     0.000000000 XXYY=    -2.016868744 XXZZ=    -1.908606236 YYZZ=    -1.581053691
 XXYZ=     0.000000000 YYXZ=     0.000655311 ZZXY=     0.000000000

*/
{
	int l, l1, l2, l3;
	char str[2000];
	FILE *fp;
	char str1[500];
	char str2[500];

	fp = Open_File (filename, "r");

/*CHARGE*/
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s %*lf %s", str1, str2);
			if (strcmp(str1, "Charge=") == 0)
				if (strcmp(str2, "electrons") == 0)
					break;
		}
		else if (feof(fp))
		{
			printf("Could not find:   Charge=     0.000000000 electrons: in %s\n", filename);
			exit(0);
		}
				
	}
	sscanf(str, "%*s %lf", &(cart[0][0][0]) );

/*DIPOLE*/
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf", &(cart[1][0][0]), &(cart[0][1][0]), &(cart[0][0][1]) );

/*QUADRUPOLE*/
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf", &(cart[2][0][0]), &(cart[0][2][0]), &(cart[0][0][2]));
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf", &(cart[1][1][0]), &(cart[1][0][1]), &(cart[0][1][1]));
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);

/*OCTAPOLE*/
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf %*s %lf", &(cart[3][0][0]), &(cart[0][3][0]), &(cart[0][0][3]), &(cart[1][2][0]));
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf %*s %lf", &(cart[2][1][0]), &(cart[2][0][1]), &(cart[1][0][2]), &(cart[0][1][2]));
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf", &(cart[0][2][1]), &(cart[1][1][1]) );

/*HEXADECAPOLE*/
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf %*s %lf", &(cart[4][0][0]), &(cart[0][4][0]), &(cart[0][0][4]), &(cart[3][1][0]));
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf %*s %lf", &(cart[3][0][1]), &(cart[1][3][0]), &(cart[0][3][1]), &(cart[1][0][3]));
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf %*s %lf", &(cart[0][1][3]), &(cart[2][2][0]), &(cart[2][0][2]), &(cart[0][2][2]));
	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %lf %*s %lf %*s %lf",         &(cart[2][1][1]), &(cart[1][2][1]), &(cart[1][1][2]) );

	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nRead_Gauss_Multipoles for %s\n", filename);
		for (l = 0; l <= 4; l++)
			for (l1 = 0; l1 <= l; l1++)
				for (l2 = 0; l2 <= l - l1; l2++)
				{
					l3 = l - l1 - l2;
					Nxyz_Index_2_String (str, l1, l2, l3);
					fprintf(fp_debug, "%4s %12.8f\n", str, cart[l1][l2][l3]);
				}
	}

/*Applequist divides multipoles by l!
	for (l = 0; l <= 4; l++)
		for (l1 = 0; l1 <= l; l1++)
			for (l2 = 0; l2 <= l - l1; l2++)
			{
				l3 = l - l1 - l2;
				cart[l1][l2][l3] /= factorial (l);
			}

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nAfter dividing by l! %s\n", filename);
		for (l = 0; l <= 4; l++)
			for (l1 = 0; l1 <= l; l1++)
				for (l2 = 0; l2 <= l - l1; l2++)
				{
					l3 = l - l1 - l2;
					Nxyz_Index_2_String (str, l1, l2, l3);
					fprintf(fp_debug, "%4s %12.8f\n", str, cart[l1][l2][l3]);
				}
	}
*/
}

void Read_ChelpG_Charges (char *filename, int n_atom, double **chelpg_charge, int DEBUG, FILE *fp_debug)
/*
 Charges from ESP fit, RMS=   0.00388 RRMS=   0.19719:
 Charge=   0.00000 Dipole=     0.3480    -1.8620     0.0262 Tot=     1.8944
              1
     1  O   -0.674036
     2  H    0.337123
     3  H    0.336913
*/
{
	int i;
	FILE *fp;
	char str[2000];
	char s_temp[5][500];
	int i_temp;
	double sum;

	fp = Open_File (filename, "r");
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s %s %s %s", s_temp[0], s_temp[1], s_temp[2], s_temp[3]);
			if (strcmp(s_temp[0], "Charges") == 0)
				if (strcmp(s_temp[1], "from") == 0)
					if (strcmp(s_temp[2], "ESP") == 0)
						if (strcmp(s_temp[3], "fit,") == 0)
							break;
		}
		else
		{
			printf("ERROR reading %s.. could not find ' Charges from ESP fit,'..\n", filename);
			exit(0);
		}
	}
	fgets(str, sizeof(str), fp);
	fgets(str, sizeof(str), fp);

	(*chelpg_charge) = D_Allocate_1D_Matrix (n_atom);
	for (i = 0; i < n_atom; i++)
	{
		fgets(str, sizeof(str), fp);
		sscanf(str, "%d %*s %lf", &i_temp, &((*chelpg_charge)[i]) );
		if (i_temp != i+1)
		{
			printf("ERROR reading %s while reading in charges.. found '%s'\n", filename, str);
			exit(0);
		}
	}

	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nRead_ChelpG_Charges %s\n", filename);
		sum = 0;
		for (i = 0; i < n_atom; i++)
		{
			fprintf(fp_debug, "\t%10.8f\n", (*chelpg_charge)[i]);
			sum += (*chelpg_charge)[i];
		}
		fprintf(fp_debug, "sum = %10.8f\n", sum);
	}




}


void Allocate_Data (struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid_hi,
	int PRINT_FIELD, int DEBUG, FILE *fp_debug)
/*This function will allocate data structures needed in the evaluation of Hirshfeld & Hirshfeld-Iterated charges,
  the evaluation of atomic multipoles, and the evaluation of the molecular multipoles*/
{
	int i, j, l;
	int ele_num;
	int lmax;


	(*grid_hi).p_MM_atom    = D_Allocate_2D_Matrix ((*mol).n_atom, (*grid_hi).n_grid_pts );
	(*grid_hi).p_MM_atom_H  = D_Allocate_2D_Matrix ((*mol).n_atom, (*grid_hi).n_grid_pts );
	(*grid_hi).p_MM_atom_HI = D_Allocate_2D_Matrix ((*mol).n_atom, (*grid_hi).n_grid_pts );

	(*grid_hi).p_QM_atom = (double ***) malloc ( (*mol).n_atom*sizeof(double **) );
	for (i = 0; i < (*mol).n_atom; i++)
	{
		ele_num = (*mol).element_num[i];
		(*grid_hi).p_QM_atom[i] = (double **) malloc ( (*ele).n_ion_state[ele_num]*sizeof(double *) );
		for (j = 0; j < (*ele).n_ion_state[ele_num]; j++)
			(*grid_hi).p_QM_atom[i][j] = D_Allocate_1D_Matrix ( (*grid_hi).n_grid_pts );
	}

	(*grid_hi).p0_tot    = D_Allocate_1D_Matrix ( (*grid_hi).n_grid_pts );
	(*grid_hi).p0_H_tot  = D_Allocate_1D_Matrix ( (*grid_hi).n_grid_pts );
	(*grid_hi).p0_HI_tot = D_Allocate_1D_Matrix ( (*grid_hi).n_grid_pts );

	(*mol).n_electron  = D_Allocate_1D_Matrix ((*mol).n_atom);
	(*mol).H_charge    = D_Allocate_1D_Matrix ((*mol).n_atom);
	(*mol).HI_charge   = D_Allocate_1D_Matrix ((*mol).n_atom);
	(*mol).atom_charge = D_Allocate_1D_Matrix ((*mol).n_atom);


	for (i = 0; i < (*mol).n_atom; i++)
		(*mol).atom_charge[i] = 0.0;

	lmax = (*mol).lmax_atom;
	if (lmax < 4)
		lmax = 4;


	(*mol).cart_TR_atom = Allocate_Cart_Tensor ((*mol).lmax_atom);
	(*mol).QM_cart_multi_TR = Allocate_Cart_Tensor (4);
	(*mol).C_lm_r        = Allocate_Spherical_Array ( lmax);
	(*mol).C_lm_i        = Allocate_Spherical_Array ( lmax);
	(*mol).Q_temp_lm_r   = Allocate_Spherical_Array ( lmax);
	(*mol).Q_temp_lm_i   = Allocate_Spherical_Array ( lmax);
	(*mol).Q_QM_mol_lm_r = Allocate_Spherical_Array ( 4);
	(*mol).Q_QM_mol_lm_i = Allocate_Spherical_Array ( 4);

	(*mol).Qa_H_lm_r  = (double ***) malloc ( (*mol).n_atom*sizeof(double **) );
	(*mol).Qa_H_lm_i  = (double ***) malloc ( (*mol).n_atom*sizeof(double **) );
	(*mol).Qa_HI_lm_r = (double ***) malloc ( (*mol).n_atom*sizeof(double **) );
	(*mol).Qa_HI_lm_i = (double ***) malloc ( (*mol).n_atom*sizeof(double **) );
	for (i = 0; i < (*mol).n_atom; i++)
	{
		(*mol).Qa_H_lm_r[i]  = Allocate_Spherical_Array ( (*mol).lmax_atom);
		(*mol).Qa_H_lm_i[i]  = Allocate_Spherical_Array ( (*mol).lmax_atom);
		(*mol).Qa_HI_lm_r[i] = Allocate_Spherical_Array ( (*mol).lmax_atom);
		(*mol).Qa_HI_lm_i[i] = Allocate_Spherical_Array ( (*mol).lmax_atom);
	}

	(*mol).Q_H_mol_lm_r     = (double ***)  malloc ( ((*mol).lmax_atom+1)*sizeof(double **) );
	(*mol).Q_H_mol_lm_i     = (double ***)  malloc ( ((*mol).lmax_atom+1)*sizeof(double **) );
	(*mol).Q_HI_mol_lm_r    = (double ***)  malloc ( ((*mol).lmax_atom+1)*sizeof(double **) );
	(*mol).Q_HI_mol_lm_i    = (double ***)  malloc ( ((*mol).lmax_atom+1)*sizeof(double **) );

	(*mol).MM_H_cart_multi  = (double ****) malloc ( ((*mol).lmax_atom+1)*sizeof(double ***) );
	(*mol).MM_HI_cart_multi = (double ****) malloc ( ((*mol).lmax_atom+1)*sizeof(double ***) );

	for (l = 0; l <= (*mol).lmax_atom; l++)
	{
		(*mol).Q_H_mol_lm_r[l]     = Allocate_Spherical_Array ( 4);
		(*mol).Q_H_mol_lm_i[l]     = Allocate_Spherical_Array ( 4);
		(*mol).Q_HI_mol_lm_r[l]    = Allocate_Spherical_Array ( 4);
		(*mol).Q_HI_mol_lm_i[l]    = Allocate_Spherical_Array ( 4);

		(*mol).MM_H_cart_multi[l]  = Allocate_Cart_Tensor (4);
		(*mol).MM_HI_cart_multi[l] = Allocate_Cart_Tensor (4);
	}

}



void Print_Output (double cpu_time, char *filename, struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid_hi,
	int PRINT_FIELD, int DEBUG, FILE *fp_debug)
/*This function prints the output file containing the molecular geometry, Hirshfeld & Hirshfeld-Iterated atomic charges, atomic multipoles (both
  spherical tensor and traceless Cartesian, the QM molecular multipoles, and the Hirshfeld, Hirshfeld-Iterated molecular multipoles for different
  atomic multipole ranks*/
{
	int i, p;
	FILE *fp;
	double sum[3];
	char initial[1000];
	double H_dipole[4];
	double HI_dipole[4];
	double chelpg_dipole[4];
	double QM_mag;
	int l, l1, l2, l3, la;
	char str[2000];
	char flag[5][500];
	double RMSD[5][100];
	double QM, MM;
	int count;
	double H_RMSD, HI_RMSD, CG_RMSD;
	double QM_dip[3];
	char TRUE_FALSE[2][50];

	Remove_Extension (filename, initial);

	fp = Open_File (filename, "w");
	
	fprintf(fp, "@Geometry\n");
	fprintf(fp, "\tn_atoms %d\n", (*mol).n_atom);
	for (i = 0; i < (*mol).n_atom; i++)
		fprintf(fp, "\t\t%4s %15.10f %15.10f %15.10f\n", (*mol).element[i], (*mol).r2d[i][0], (*mol).r2d[i][1], (*mol).r2d[i][2]);
	fprintf(fp, "\n\n");

	fprintf(fp, "@Multi\n");
	fprintf(fp, "\tcharge %d\n", (*mol).charge);
	fprintf(fp, "\tmulti  %d\n", (*mol).multi);
	fprintf(fp, "\n\n");

	sprintf(TRUE_FALSE[0], "FALSE");
	sprintf(TRUE_FALSE[1], "TRUE");

	fprintf(fp, "@IO_options\n");
	fprintf(fp, "\ttheory = %s\n", (*mol).theory);
	fprintf(fp, "\tnproc  = %d\n", (*mol).nproc );
	fprintf(fp, "\tmem    = %s\n", (*mol).memory );
	fprintf(fp, "\tlmax   = %d\n", (*mol).lmax_atom );
	fprintf(fp, "\tChelpG = %s\n", TRUE_FALSE[ (*mol).DO_CHELPG ] );
	fprintf(fp, "\n\n");


	fprintf(fp, "@CPU_TIME\n");
	fprintf(fp, "\tcpu_time = %f s\n", cpu_time);
	fprintf(fp, "\n\n");

	fprintf(fp, "@Charge\n");
	fprintf(fp, "\tn_atom %d\n", (*mol).n_atom);
	fprintf(fp, "     element        Hirshfeld      Hirshfeld-Iterated       ChelpG\n");

	sum[0] = 0.0;	sum[1] = 0.0;	sum[2] = 0.0;

	for (i = 0; i < (*mol).n_atom; i++)
	{
		fprintf(fp, "\t\t%4s      %12.8f       %12.8f       %12.8f\n", (*mol).element[i], (*mol).H_charge[i], (*mol).HI_charge[i], (*mol).chelpg_charge[i]);
		sum[0] += (*mol).H_charge[i];
		sum[1] += (*mol).HI_charge[i];
		sum[2] += (*mol).chelpg_charge[i];
	}
	fprintf(fp, "----------------------------------------------------------------\n");
	fprintf(fp, "sum               %12.8f       %12.8f       %12.8f\n", sum[0], sum[1], sum[2]);
	fprintf(fp, "\n\n");

	fprintf(fp, "\n");
	Calc_Molecular_Dipole (&(H_dipole[0]),      (*mol).n_atom, (*mol).r2d, (*mol).H_charge);
	Calc_Molecular_Dipole (&(HI_dipole[0]),     (*mol).n_atom, (*mol).r2d, (*mol).HI_charge);
	Calc_Molecular_Dipole (&(chelpg_dipole[0]), (*mol).n_atom, (*mol).r2d, (*mol).chelpg_charge);

	for (p = 0; p < 3; p++)
	{
		H_dipole[p]      *= eA_2_DEBYE;
		HI_dipole[p]     *= eA_2_DEBYE;
		chelpg_dipole[p] *= eA_2_DEBYE;
	}
	H_dipole[3]       = sqrt(dot_prod(&(H_dipole[0]),      &(H_dipole[0])) );
	HI_dipole[3]      = sqrt(dot_prod(&(HI_dipole[0]),     &(HI_dipole[0])) );
	chelpg_dipole[3]  = sqrt(dot_prod(&(chelpg_dipole[0]), &(chelpg_dipole[0])) );
	QM_mag            = sqrt((*mol).QM_cart_multi[1][0][0]*(*mol).QM_cart_multi[1][0][0] + 
					(*mol).QM_cart_multi[0][1][0]*(*mol).QM_cart_multi[0][1][0] + (*mol).QM_cart_multi[0][0][1]*(*mol).QM_cart_multi[0][0][1]);

	H_RMSD  = 0;
	HI_RMSD = 0;
	CG_RMSD = 0;
	QM_dip[0] = (*mol).QM_cart_multi[1][0][0];
	QM_dip[1] = (*mol).QM_cart_multi[0][1][0];
	QM_dip[2] = (*mol).QM_cart_multi[0][0][1];
	
	for (p = 0; p < 3; p++)
	{
		H_RMSD   += (H_dipole[p] - QM_dip[p])     *(H_dipole[p] - QM_dip[p]);
		HI_RMSD  += (HI_dipole[p] - QM_dip[p])    *(HI_dipole[p] - QM_dip[p]);
		CG_RMSD  += (chelpg_dipole[p] - QM_dip[p])*(chelpg_dipole[p] - QM_dip[p]);
	}

	H_RMSD  /= 3;
	HI_RMSD /= 3;
	CG_RMSD /= 3;

	H_RMSD  = sqrt(H_RMSD);
	HI_RMSD = sqrt(HI_RMSD);
	CG_RMSD = sqrt(CG_RMSD);


	fprintf(fp, "@Dipole\n");
	fprintf(fp, "\tH_dipole       %8.4f %8.4f %8.4f      RMSD = %8.4f\n", H_dipole[0],      H_dipole[1],      H_dipole[2],      H_RMSD);
	fprintf(fp, "\tHI_dipole      %8.4f %8.4f %8.4f      RMSD = %8.4f\n", HI_dipole[0],     HI_dipole[1],     HI_dipole[2],     HI_RMSD);
	fprintf(fp, "\tChelpg_dipole  %8.4f %8.4f %8.4f      RMSD = %8.4f\n", chelpg_dipole[0], chelpg_dipole[1], chelpg_dipole[2], CG_RMSD);
	fprintf(fp, "\tQM_dipole      %8.4f %8.4f %8.4f\n", (*mol).QM_cart_multi[1][0][0], (*mol).QM_cart_multi[0][1][0], (*mol).QM_cart_multi[0][0][1]);

	fprintf(fp, "\n");
	fprintf(fp, "\tH_dipole       %8.4f\n", H_dipole[3]);
	fprintf(fp, "\tHI_dipole      %8.4f\n", HI_dipole[3]);
	fprintf(fp, "\tChelpg_dipole  %8.4f\n", chelpg_dipole[3]);
	fprintf(fp, "\t-----------------------\n");
	fprintf(fp, "\tQM_dipole      %8.4f\n", QM_mag);
	fprintf(fp, "\n\n\n\n\n");



	fprintf(fp, "\n\n@Hirshfeld_Atomic_Multipoles (e*A^l), Qlm = Re(Qlm) + Im(Qlm)\n");
	fprintf(fp, "\tn_atom %d lmax %d\n", (*mol).n_atom, (*mol).lmax_atom);
	for (i = 0; i < (*mol).n_atom; i++)
	{
		fprintf(fp, "\t\t%4d %4s\n", i, (*mol).element[i]);
		for (l = 0; l <= (*mol).lmax_atom; l++)
			for (p = 0; p <= l; p++)
				fprintf(fp, "\t\t\tQ%d%d = %15.12f + %15.12f\n", l, p, (*mol).Qa_H_lm_r[i][l][p], (*mol).Qa_H_lm_i[i][l][p]);
	}

	fprintf(fp, "\n\n@Hirshfeld_I_Atomic_Multipoles (e*A^l), Qlm = Re(Qlm) + Im(Qlm)\n");
	fprintf(fp, "\tn_atom %d lmax %d\n", (*mol).n_atom, (*mol).lmax_atom);
	for (i = 0; i < (*mol).n_atom; i++)
	{
		fprintf(fp, "\t\t%4d %4s\n", i, (*mol).element[i]);
		for (l = 0; l <= (*mol).lmax_atom; l++)
			for (p = 0; p <= l; p++)
				fprintf(fp, "\t\t\tQ%d%d = %15.12f + %15.12f\n", l, p, (*mol).Qa_HI_lm_r[i][l][p], (*mol).Qa_HI_lm_i[i][l][p]);
	}

	fprintf(fp, "\n\n@Hirshfeld_Atomic_Multipoles_Cartesian (e*A^l)\n");
	fprintf(fp, "\tn_atom %d lmax %d\n", (*mol).n_atom, (*mol).lmax_atom);
	for (i = 0; i < (*mol).n_atom; i++)
	{
		fprintf(fp, "\t\t%4d %4s\n", i, (*mol).element[i]);
		Convert_Spherical_Multipole_2_TR_Cart_Multipole ((*mol).lmax_atom, (*mol).cart_TR_atom, (*mol).Qa_H_lm_r[i], 
			(*mol).Qa_H_lm_i[i], DEBUG, fp_debug);

		for (l = 0; l <= (*mol).lmax_atom; l++)
		{
			for (l1 = l; l1 >= 0; l1--)
				for (l2 = l - l1; l2 >= 0; l2--)
				{
					l3 = l - l1 - l2;
					Nxyz_Index_2_String (str, l1, l2, l3);
					if (l != 0)
						fprintf(fp, "\t%10s = %15.12f\n", str, (*mol).cart_TR_atom[l1][l2][l3]*double_factorial (2*l-1));
					else
						fprintf(fp, "\t%10s = %15.12f\n", "MONOPOLE", (*mol).cart_TR_atom[l1][l2][l3]*double_factorial (2*l-1));
				}
		}
	}

	fprintf(fp, "\n\n@Hirshfeld_I_Atomic_Multipoles_Cartesian (e*A^l)\n");
	fprintf(fp, "\tn_atom %d lmax %d\n", (*mol).n_atom, (*mol).lmax_atom);
	for (i = 0; i < (*mol).n_atom; i++)
	{
		fprintf(fp, "\t\t%4d %4s\n", i, (*mol).element[i]);
		Convert_Spherical_Multipole_2_TR_Cart_Multipole ((*mol).lmax_atom, (*mol).cart_TR_atom, (*mol).Qa_HI_lm_r[i], 
			(*mol).Qa_HI_lm_i[i], DEBUG, fp_debug);

		for (l = 0; l <= (*mol).lmax_atom; l++)
		{
			for (l1 = l; l1 >= 0; l1--)
				for (l2 = l - l1; l2 >= 0; l2--)
				{
					l3 = l - l1 - l2;
					Nxyz_Index_2_String (str, l1, l2, l3);
					if (l != 0)
						fprintf(fp, "\t%10s = %15.12f\n", str, (*mol).cart_TR_atom[l1][l2][l3]*double_factorial (2*l-1));
					else
						fprintf(fp, "\t%10s = %15.12f\n", "MONOPOLE", (*mol).cart_TR_atom[l1][l2][l3]*double_factorial (2*l-1));
				}
		}
	}

	sprintf(flag[1], "@Hirshfeld_Molecular_Dipole       (D)");
	sprintf(flag[2], "@Hirshfeld_Molecular_Quadrupole   (D*A)");
	sprintf(flag[3], "@Hirshfeld_Molecular_Octupole     (D*A^2)");
	sprintf(flag[4], "@Hirshfeld_Molecular_Hexadecapole (D*A^3)");

	for (l = 1; l <= 4; l++)
	{
		fprintf(fp, "\n\n\n%s\n", flag[l]);
		fprintf(fp, "       ");
		for (l1 = l; l1 >= 0; l1--)
			for (l2 = l - l1; l2 >= 0; l2--)
			{
				l3 = l - l1 - l2;
				Nxyz_Index_2_String (str, l1, l2, l3);
				fprintf(fp, "%8s ", str);
			}
		fprintf(fp, "%7s ", "RMSD");
		fprintf(fp, "\n");
		
		for (la = 0; la <= (*mol).lmax_atom; la++)
		{
			RMSD[l][la] = 0;
			count = 0;
			fprintf(fp, "    am_%d ", la);
			for (l1 = l; l1 >= 0; l1--)
				for (l2 = l - l1; l2 >= 0; l2--)
				{
					l3 = l - l1 - l2;
					MM = (*mol).MM_H_cart_multi[la][l1][l2][l3]*eA_2_DEBYE*double_factorial (2*l-1);
					QM = (*mol).QM_cart_multi_TR[l1][l2][l3]*double_factorial (2*l-1);
					fprintf(fp, "%8.4f ",  MM);
					RMSD[l][la] += (MM - QM)*(MM - QM);
					count++;
				}
			RMSD[l][la] /= count;
			RMSD[l][la] = sqrt(RMSD[l][la]);
			fprintf(fp, "%8.4f ",  RMSD[l][la]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "----------------");
		for (l1 = l; l1 >= 0; l1--)
			for (l2 = l - l1; l2 >= 0; l2--)
				fprintf(fp, "--------");
		fprintf(fp, "\n");
		fprintf(fp, "    QM   ");
		for (l1 = l; l1 >= 0; l1--)
			for (l2 = l - l1; l2 >= 0; l2--)
			{
				l3 = l - l1 - l2;
				Nxyz_Index_2_String (str, l1, l2, l3);
				fprintf(fp, "%8.4f ",  (*mol).QM_cart_multi_TR[l1][l2][l3]*double_factorial (2*l-1));
			}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n\n");
	for (l = 1; l <= 4; l++)
	{
		for (la = 0; la <= (*mol).lmax_atom; la++)
			fprintf(fp, "L_%d_RMSD_am_%d = %f\n", l, la, RMSD[l][la]);
		fprintf(fp, "\n");
	}

	sprintf(flag[1], "@Hirshfeld_I_Molecular_Dipole       (D)");
	sprintf(flag[2], "@Hirshfeld_I_Molecular_Quadrupole   (D*A)");
	sprintf(flag[3], "@Hirshfeld_I_Molecular_Octupole     (D*A^2)");
	sprintf(flag[4], "@Hirshfeld_I_Molecular_Hexadecapole (D*A^3)");
	for (l = 1; l <= 4; l++)
	{
		fprintf(fp, "\n\n\n%s\n", flag[l]);
		fprintf(fp, "       ");
		for (l1 = l; l1 >= 0; l1--)
			for (l2 = l - l1; l2 >= 0; l2--)
			{
				l3 = l - l1 - l2;
				Nxyz_Index_2_String (str, l1, l2, l3);
				fprintf(fp, "%8s ", str);
			}
		fprintf(fp, "%7s ", "RMSD");
		fprintf(fp, "\n");
		
		for (la = 0; la <= (*mol).lmax_atom; la++)
		{
			RMSD[l][la] = 0;
			count = 0;
			fprintf(fp, "    am_%d ", la);
			for (l1 = l; l1 >= 0; l1--)
				for (l2 = l - l1; l2 >= 0; l2--)
				{
					l3 = l - l1 - l2;
					MM = (*mol).MM_HI_cart_multi[la][l1][l2][l3]*eA_2_DEBYE*double_factorial (2*l-1);
					QM = (*mol).QM_cart_multi_TR[l1][l2][l3]*double_factorial (2*l-1);
					fprintf(fp, "%8.4f ",  MM);
					RMSD[l][la] += (MM - QM)*(MM - QM);
					count++;
				}
			RMSD[l][la] /= count;
			RMSD[l][la] = sqrt(RMSD[l][la]);
			fprintf(fp, "%8.4f ",  RMSD[l][la]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "----------------");
		for (l1 = l; l1 >= 0; l1--)
			for (l2 = l - l1; l2 >= 0; l2--)
				fprintf(fp, "--------");
		fprintf(fp, "\n");
		fprintf(fp, "    QM   ");
		for (l1 = l; l1 >= 0; l1--)
			for (l2 = l - l1; l2 >= 0; l2--)
			{
				l3 = l - l1 - l2;
				Nxyz_Index_2_String (str, l1, l2, l3);
				fprintf(fp, "%8.4f ",  (*mol).QM_cart_multi_TR[l1][l2][l3]*double_factorial (2*l-1));
			}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n\n");
	for (l = 1; l <= 4; l++)
	{
		for (la = 0; la <= (*mol).lmax_atom; la++)
			fprintf(fp, "L_%d_RMSD_am_%d = %f\n", l, la, RMSD[l][la]);
		fprintf(fp, "\n");
	}
		
	fclose(fp);
	


}

void Calc_Molecular_Dipole (double *dip, int n_atom, double **r2d, double *charge)
/*calculates molecular dipole from the molecular geometry and atmic charges*/
{
	int i, p;

	for (p = 0; p < 3; p++)
		dip[p] = 0;

	for (i = 0; i < n_atom; i++)
		for (p = 0; p < 3; p++)
			dip[p] += r2d[i][p]*charge[i];

}


void Nxyz_Index_2_String (char *str, int n1, int n2, int n3)
{
	int i;
	for (i = 0; i < n1; i++)		str[i]       = 'X';
	for (i = 0; i < n2; i++)		str[i+n1]    = 'Y';
	for (i = 0; i < n3; i++)		str[i+n1+n2] = 'Z';
	str[n1+n2+n3] = '\0';
}

