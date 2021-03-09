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


#include "Allocate.h"

#define MAX_LINE_LENGTH 2000
#define S_ATOM_TYPE_LENGTH 10
#define AT_STRING_LENGTH 10
#define S_FILE_LENGTH 2000
#define MAX_NEIGHBORS 6
#define TRUE 1
#define FALSE 0

/*Used in reading .eg file*/
#define E_MIN -20.0
#define E_MAX 70.0


/*MISC FUNCTIONS*/
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

void Read_File_List (char *path, int *n_files, char ***filenames, char *filelist, int DEBUG, FILE *fp_debug)
/*
/home/denny/temp
	num_files 3
		temp1.dat
		temp2.dat
		temp3.dat
*/			

{
	int i;
	char str[2000];
	FILE *fp;

	fp = fopen( filelist, "r");
	if (fp == NULL) 
	{
		printf("Cant open %s\n", filelist);
		exit(EXIT_FAILURE);
	}

	fgets(str, sizeof(str), fp);
	sscanf(str, "%s", path);

	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %d", &(*n_files) );
	(*filenames) = C_Allocate_2D_Matrix ( *n_files, S_FILE_LENGTH);

	for (i = 0; i < (*n_files); i++)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
			sscanf(str, "%s", (*filenames)[i] );
		else
		{
			printf("EOF while reading in %s\n", filelist);
			exit(0);
		}
	}

	fclose(fp);
/*Print Out to Debug File */
	if (DEBUG)
	{
		fprintf(fp_debug, "Reading %s\n", filelist);
		fprintf(fp_debug, "\tpath = %s\n", path);
		fprintf(fp_debug, "\tn_files = %d\n", *n_files);
		for (i = 0; i < (*n_files); i++)
			fprintf(fp_debug, "\t\t%s\n", (*filenames)[i]);
		fflush(fp_debug);
	}
}


void Read_Atom_Types (int *n_atom_types, char ***atom_types, char *filename, 
	char *flag, FILE *fp_debug, int DEBUG)
/*
	num_atom_types 2
		H  
		HO 
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &(*n_atom_types) );
	(*atom_types) = C_Allocate_2D_Matrix ((*n_atom_types), S_ATOM_TYPE_LENGTH);
	for (i = 0; i < (*n_atom_types); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s", (*atom_types)[i]);
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_atom_types = %d\n", (*n_atom_types) );
		for (i = 0; i < (*n_atom_types); i++)
			fprintf(fp_debug, "\t\t\t%4d %4s\n", i, (*atom_types)[i]);
		fflush(fp_debug);
	}
}
	
void Read_AT_Double (double **real_value, int n_atom_types, char **atom_types, 
	char *filename, char *flag, FILE *fp_debug, int DEBUG)
/*
	n_polar 2
		H  	0.21
		HO 	0.25
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	int n_real;
	double value;
	char s_atom_type[S_ATOM_TYPE_LENGTH];
	int at_num;

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &n_real );

	(*real_value) = D_Allocate_1D_Matrix (n_atom_types);
	for (i = 0; i < n_atom_types; i++)
		(*real_value)[i] = -1.0;

	for (i = 0; i < n_real; i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %lf", s_atom_type, &value);
		at_num = Find_Atom_Type (s_atom_type, n_atom_types, atom_types);
		(*real_value)[at_num] = value;
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_reals = %d\n", n_real );
		for (i = 0; i < n_atom_types; i++)
			fprintf(fp_debug, "\t\t\t%4d %4s %10.5f\n", i, atom_types[i], (*real_value)[i]);
		fflush(fp_debug);
	}
}

void Read_AT2_Double (int *n_real, double **real_value, int ***ab_num_2_atom_type, int n_atom_types, 
	char **atom_types, char *filename, char *flag, FILE *fp_debug, int DEBUG)
/*
	n_parm 2
		H  	HO	0.21
		HO 	CT	0.25
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	char s_atom_type1[S_ATOM_TYPE_LENGTH];
	char s_atom_type2[S_ATOM_TYPE_LENGTH];

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &(*n_real) );

	(*real_value)			= D_Allocate_1D_Matrix (*n_real);
	(*ab_num_2_atom_type)	= I_Allocate_2D_Matrix (*n_real, 2);

	for (i = 0; i < (*n_real); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %s %lf", s_atom_type1, s_atom_type2, &( (*real_value)[i]));
		(*ab_num_2_atom_type)[i][0] = Find_Atom_Type (s_atom_type1, n_atom_types, atom_types);
		(*ab_num_2_atom_type)[i][1] = Find_Atom_Type (s_atom_type2, n_atom_types, atom_types);
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_reals = %d\n", (*n_real) );
		for (i = 0; i < (*n_real); i++)
			fprintf(fp_debug, "\t\t\t%4d %4s %4s %10.5f\n", i, 
				atom_types[ (*ab_num_2_atom_type)[i][0] ], 
				atom_types[ (*ab_num_2_atom_type)[i][1] ], (*real_value)[i]);
		fflush(fp_debug);
	}
}


void Read_AT_Double_2 (double **real_value1, double **real_value2, int n_atom_types, 
	char **atom_types, char *filename, char *flag, FILE *fp_debug, int DEBUG)
/*
	n_vdw 2
		H  	1.21 0.01
		HO 	0.50 0.01
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	int n_real;
	double value1, value2;
	char s_atom_type[S_ATOM_TYPE_LENGTH];
	int at_num;

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &n_real );

	(*real_value1) = D_Allocate_1D_Matrix (n_atom_types);
	(*real_value2) = D_Allocate_1D_Matrix (n_atom_types);

	for (i = 0; i < n_atom_types; i++)
	{
		(*real_value1)[i] = -1.0;
		(*real_value2)[i] = -1.0;
	}
	for (i = 0; i < n_real; i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %lf %lf", s_atom_type, &value1, &value2);
		at_num = Find_Atom_Type (s_atom_type, n_atom_types, atom_types);
		(*real_value1)[at_num] = value1;
		(*real_value2)[at_num] = value2;
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_reals = %d\n", n_real );
		for (i = 0; i < n_atom_types; i++)
			fprintf(fp_debug, "\t\t\t%4d %4s %10.5f %10.5f\n", i, atom_types[i], 
				(*real_value1)[i], (*real_value2)[i]);
		fflush(fp_debug);
	}
}

void Read_AT2_Double2 (int *n_real, double **real_value1, double **real_value2, 
	int ***ab_num_2_atom_type, int n_atom_types, char **atom_types, char *filename, 
	char *flag, FILE *fp_debug, int DEBUG)
/*
	n_parm 3
		HC CT	0.00	0.00
		OH HO	2.22	0.0448
		OH OH	3.44 	0.2104
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	char s_atom_type1[S_ATOM_TYPE_LENGTH];
	char s_atom_type2[S_ATOM_TYPE_LENGTH];

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &(*n_real) );

	(*real_value1)			= D_Allocate_1D_Matrix (*n_real);
	(*real_value2)			= D_Allocate_1D_Matrix (*n_real);
	(*ab_num_2_atom_type)	= I_Allocate_2D_Matrix (*n_real, 2);

	for (i = 0; i < (*n_real); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %s %lf %lf", s_atom_type1, s_atom_type2, &( (*real_value1)[i]), &( (*real_value2)[i]) );
		(*ab_num_2_atom_type)[i][0] = Find_Atom_Type (s_atom_type1, n_atom_types, atom_types);
		(*ab_num_2_atom_type)[i][1] = Find_Atom_Type (s_atom_type2, n_atom_types, atom_types);
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_reals = %d\n", (*n_real) );
		for (i = 0; i < (*n_real); i++)
			fprintf(fp_debug, "\t\t\t%4d %4s %4s %10.5f %10.5f\n", i, 
				atom_types[ (*ab_num_2_atom_type)[i][0] ], 
				atom_types[ (*ab_num_2_atom_type)[i][1] ], (*real_value1)[i], (*real_value2)[i]);
		fflush(fp_debug);
	}
}

void Read_AT_String (int n_atom_types, char **atom_type_names, 
		char ***string, char *filename, char *flag, int DEBUG, FILE *fp_debug)
/*
@Elements
	num_elements 18
		H		H
		HO		H
		HC		H
*/
{
	int i;
	int n_string;
	char s_atom_type[20];
	char s_string[AT_STRING_LENGTH];
	int at_num;
	char str[MAX_LINE_LENGTH];
	FILE *fp;

	fp = fopen( filename, "r");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	(*string) = C_Allocate_2D_Matrix (n_atom_types, AT_STRING_LENGTH);
	for (i = 0; i < n_atom_types; i++)
		sprintf( (*string)[i], "NONE");

    Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &n_string );
	for (i = 0; i < n_string; i++)
	{
		fgets(str, sizeof(str), fp);
		sscanf(str, "%s %s", s_atom_type, s_string);
		at_num = Find_Atom_Type (s_atom_type, n_atom_types, atom_type_names);
		strcpy( (*string)[at_num], s_string);
	}

	fclose(fp);
	
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tReading Element Definitions\n");
		fprintf(fp_debug, "\t\t\tn_string = %d\n", n_string);
		for (i = 0; i < n_atom_types; i++)
			fprintf(fp_debug, "\t%4s %4s\n", atom_type_names[i], (*string)[i]);
		fflush(fp_debug);
	}
}



void Read_Opt_AT_Parameters (int *n_opt, int **opt_num_2_at_num, int n_atom_types, 
	char **atom_types, char *filename, char *flag, FILE *fp_debug, int DEBUG)
/*
	num_atom_types 2
		H  
		HO 
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	char s_atom_type[S_ATOM_TYPE_LENGTH];

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &(*n_opt) );

	(*opt_num_2_at_num) = I_Allocate_1D_Matrix ( (*n_opt) );

	for (i = 0; i < (*n_opt); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s", s_atom_type);
		(*opt_num_2_at_num)[i] = Find_Atom_Type (s_atom_type, n_atom_types, atom_types);
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_opt_atom_types = %d\n", (*n_opt) );
		for (i = 0; i < (*n_opt); i++)
			fprintf(fp_debug, "\t\t\t%4d %4s\n", i, atom_types[(*opt_num_2_at_num)[i]]);
		fflush(fp_debug);
	}
}
	

void Read_Opt_AT_Parameters2 (int *n_opt2, int **opt_num_2_ab_num, int n_ab,
	int **ab_num_2_atom_type, int n_atom_types, char **atom_types, char *filename, 
	char *flag, FILE *fp_debug, int DEBUG)
/*
	num_atom_types 2
		H  	H
		HO 	CT
*/
{
	int i, j;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	char s_atom_type1[S_ATOM_TYPE_LENGTH];
	char s_atom_type2[S_ATOM_TYPE_LENGTH];
	int at1, at2, ab_num;

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &(*n_opt2) );

	(*opt_num_2_ab_num) = I_Allocate_1D_Matrix ( (*n_opt2));

	for (i = 0; i < (*n_opt2); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %s", s_atom_type1, s_atom_type2);
		at1 = Find_Atom_Type (s_atom_type1, n_atom_types, atom_types);
		at2 = Find_Atom_Type (s_atom_type2, n_atom_types, atom_types);
		ab_num = -1;
		for (j = 0; j < n_ab; j++)
		{
			if (at1 == ab_num_2_atom_type[j][0])
				if (at2 == ab_num_2_atom_type[j][1])
				{
					ab_num = j;
					break;
				}
			if (at2 == ab_num_2_atom_type[j][0])
				if (at1 == ab_num_2_atom_type[j][1])
				{
					ab_num = j;
					break;
				}
		}
		if (ab_num == -1)
		{
			printf("Couldn't Find %s %s in AB parameters\n", s_atom_type1, s_atom_type2);
			exit(0);
		}
		(*opt_num_2_ab_num)[i] = ab_num;
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\t\tn_opt_atom_types2 = %d\n", (*n_opt2) );
		for (i = 0; i < (*n_opt2); i++)
		{
			ab_num = (*opt_num_2_ab_num)[i];
			at1 = ab_num_2_atom_type[ab_num][0];
			at2 = ab_num_2_atom_type[ab_num][1];
			fprintf(fp_debug, "\t\t\t%4d %4s %4s\n", i, 
				atom_types[at1], atom_types[at2]);
		}
		fflush(fp_debug);
	}
}

void Read_Bool_String_Pair ( int **is_value, char ***s_type, int n_atom_types, 
	char **atom_type_names, char *filename, char *flag, FILE *fp_debug, int DEBUG)
/*
@EP_Definition
	num_ep 4
		N2_EP 3_sp3
		OH_EP 2_sp3
		CD_EP 3_sp2
		O_EP  1_sp2
*/
{
	int i;
	char str[MAX_LINE_LENGTH];
	FILE *fp;
	char s_atom_type[S_ATOM_TYPE_LENGTH];
	char type[S_ATOM_TYPE_LENGTH];
	int n_types;
	int at_num;

	fp = Open_File (filename, "r");
	Go_To_Line (fp, flag, filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %d", &n_types );

	(*is_value)	= I_Allocate_1D_Matrix (n_atom_types);
	(*s_type)	= C_Allocate_2D_Matrix (n_atom_types, S_ATOM_TYPE_LENGTH);

/*initialize to 0 and NULL*/
	for (i = 0; i < n_atom_types; i++)
	{
		(*is_value)[i]    = 0;
		(*s_type)[i][0]   = '\0';
	}

	for (i = 0; i < n_types; i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%s %s", s_atom_type, type);
		at_num = Find_Atom_Type (s_atom_type, n_atom_types, atom_type_names);
		(*is_value)[at_num] = 1;
		sprintf( ( (*s_type)[at_num]), "%s", type);
	}
	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n%s\n", flag);
		fprintf(fp_debug, "\n\n\nInside Read_Bool_String_Pair\n");
		fprintf(fp_debug, "\t\tn_types = %d\n", n_types );
		for (i = 0; i < n_atom_types; i++)
			if ( (*is_value)[i])
				fprintf(fp_debug, "\t\t\t%5s %s\n", atom_type_names[i],
					(*s_type)[i]);
		fflush(fp_debug);
	}
}




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
		char *filename, FILE *fp_debug, int DEBUG)
/*Reads in (modified) ml2 file with charges and elements*/
{
	int i, j;
	int atom1, atom2, num_neigh1, num_neigh2;
	int num_connects;
	int **connects;
	char str[MAX_LINE_LENGTH];
	char s_temp[MAX_LINE_LENGTH];
	char atom_type_s[S_ATOM_TYPE_LENGTH];
	FILE *fp;


	fp = fopen( filename, "r");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

/*Determine Num_atoms*/
    Go_To_Line (fp, "@<TRIPOS>ATOM", filename);
	*num_atoms = 0;
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s", s_temp);
			if (strcmp(s_temp, "@<TRIPOS>BOND") == 0)
				break;
			sscanf(str, "%*d %*s %*lf %*lf %*lf %s %*d %*s %*s ", s_temp); 
			if ( (s_temp[0] != 'L') && (s_temp[1] != 'p') )
				(*num_atoms)++;
		}
		else
			IO_Error (filename);
	}


/*Determine num_connects*/
	num_connects = 0;
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s", s_temp);
			if (strcmp(s_temp, "@<TRIPOS>SUBSTRUCTURE") == 0)
				break;
			num_connects++;
		}
		else
			IO_Error (filename);
	}
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\n\tReading %s\n", filename);
		fprintf(fp_debug, "\t\tnum_atoms = %d, num_connects = %d\n",
			*num_atoms, num_connects);
		fflush(fp_debug);
	}

/*Read in Atom Types, Coordinates, Elements and Charges*/
	rewind(fp);
	Go_To_Line (fp, "@<TRIPOS>ATOM", filename);

	(*atom_type)	= I_Allocate_1D_Matrix ( *num_atoms );
	(*r_1D_atoms)	= D_Allocate_1D_Matrix ( 3*(*num_atoms) );
    (*r_2D_atoms)	= (double **) malloc((*num_atoms)*sizeof(double*) );
	(*charges)		= D_Allocate_1D_Matrix (*num_atoms);
	(*element)		= C_Allocate_2D_Matrix (*num_atoms, 5);

	for (i = 0; i < *num_atoms; i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%*d %s %lf %lf %lf %s %*d %*s %lf", 
			atom_type_s, &( (*r_1D_atoms)[3*i+0] ), &( (*r_1D_atoms)[3*i+1] ),
			&( (*r_1D_atoms)[3*i+2] ), (*element)[i], &( (*charges)[i]) ) ;
		(*atom_type)[i] = Find_Atom_Type (atom_type_s, num_atom_types, atom_type_names);
	}
	Get_Next_Line (fp, str, filename);


/*Make R_2D_atoms point to R_1D_atoms   */
    for (i = 0; i < *num_atoms; i++)
        (*r_2D_atoms)[i] = &( (*r_1D_atoms)[3*i]);
        
/*Read in Connectivities and turn into Neighbors*/
	connects		= I_Allocate_2D_Matrix ( num_connects, 2 );

	for (i = 0; i < num_connects; i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%*d %d %d", &connects[i][0], &connects[i][1]);
		connects[i][0]--;
		connects[i][1]--;
	}
	
	(*n_neighbors)	= I_Allocate_1D_Matrix ( *num_atoms );
	(*neighbors)	= (int **) malloc ( (*num_atoms)*sizeof(int *) );

	for (i = 0; i < *num_atoms; i++)
		(*n_neighbors)[i] = 0;;

	for (i = 0; i < num_connects; i++)
	{
		atom1 = connects[i][0];
		atom2 = connects[i][1];
		(*n_neighbors)[atom1]++;
		(*n_neighbors)[atom2]++;
	}
	
	for (i = 0; i < *num_atoms; i++)
		(*neighbors)[i] = (int *) malloc ( (*n_neighbors)[i]*sizeof(int) );

	for (i = 0; i < *num_atoms; i++)
		(*n_neighbors)[i] = 0;
	
	for (i = 0; i < num_connects; i++)
	{
		atom1 = connects[i][0];
		atom2 = connects[i][1];
		num_neigh1 = (*n_neighbors)[atom1];
		num_neigh2 = (*n_neighbors)[atom2];
		(*neighbors)[atom1][num_neigh1] = atom2;
		(*neighbors)[atom2][num_neigh2] = atom1;
		(*n_neighbors)[atom1]++;
		(*n_neighbors)[atom2]++;
	}

	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "\t\tnum_atoms = %d\n", *num_atoms);
		for (i = 0; i < *num_atoms; i++)
			fprintf(fp_debug, "\t\t\t%4d %4s %4s	%10.5f %10.5f %10.5f\t\t%f\n", i, 
				(*element)[i], atom_type_names[(*atom_type)[i]], (*r_1D_atoms)[3*i+0], 
                (*r_1D_atoms)[3*i+1], (*r_1D_atoms)[3*i+2], (*charges)[i]);
        
        fprintf(fp_debug, "\n\t\tR_2_D_atoms\n");
        for (i = 0; i < *num_atoms; i++)
            fprintf(fp_debug, "\t\t\t%10.5f %10.5f %10.5f\n", (*r_2D_atoms)[i][0], 
                (*r_2D_atoms)[i][1], (*r_2D_atoms)[i][2]);

		fprintf(fp_debug, "\n\n\t\t num_connects = %d\n", num_connects);
		for (i = 0; i < num_connects; i++)
			fprintf(fp_debug, "\t\t\t%4d	%4d %4d\n", i, connects[i][0], connects[i][1]);

		fprintf(fp_debug, "\n\tNeigbhors\n");
		for (i = 0; i < *num_atoms; i++)
		{
			fprintf(fp_debug, "\t\t\t%4s num_neighors = %d		", 
				atom_type_names[(*atom_type)[i]], (*n_neighbors)[i]);
			for (j = 0; j < (*n_neighbors)[i]; j++)
				fprintf(fp_debug, "\t\t\t%4d	", (*neighbors)[i][j]);
			fprintf(fp_debug, "\n");
		}

		for (i = 0; i < *num_atoms; i++)
		{
			fprintf(fp_debug, "\t\t\t%4d %4s	%10.5f %10.5f %10.5f %4d\t\t", i, 
				atom_type_names[(*atom_type)[i]], (*r_1D_atoms)[3*i+0], 
                (*r_1D_atoms)[3*i+1], (*r_1D_atoms)[3*i+2], (*n_neighbors)[i]);
			for (j = 0; j < (*n_neighbors)[i]; j++)
				fprintf(fp_debug, "%4d	", (*neighbors)[i][j]);
			fprintf(fp_debug, "\n");
		}


	}

	I_free_2D (connects, num_connects);
}


void Read_EG_File ( int *n_monomer, int **n_atom, double ****R2D_monomer,
	char ***conn_files, int *n_clusters, double *****R2D_cluster, double **QM_Energy,
	double *****QM_Force, char *filename, int DEBUG, FILE *fp_debug)
{
	int i, j, n;
	char str[3000];
	FILE *fp;
	int i_temp1, i_temp2;

	fp = Open_File (filename, "r");

/*Monomers*/
	Go_To_Line (fp, "@Molecules", filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %*s %d", &(*n_monomer) );

	(*n_atom)		= I_Allocate_1D_Matrix ( (*n_monomer) );
	(*R2D_monomer)	= (double ***) malloc ( (*n_monomer)*sizeof(double **) );
	(*conn_files)	= C_Allocate_2D_Matrix ( (*n_monomer),  2000);


	for (i = 0; i < (*n_monomer); i++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%*s %*d %*s %*s %d %*s %*s %*d %*s %*s %*d %*s %*s %s",
			&( (*n_atom)[i]), (*conn_files)[i]);
		(*R2D_monomer)[i] = D_Allocate_2D_Matrix ( (*n_atom)[i], 3);
		for (j = 0; j < (*n_atom)[i]; j++)
		{
			Get_Next_Line (fp, str, filename);
			sscanf(str, "%*d %*s %lf %lf %lf", &( (*R2D_monomer)[i][j][0]),
				&( (*R2D_monomer)[i][j][1]), &( (*R2D_monomer)[i][j][2]) );
		}
	}

	if (DEBUG)
	{
		fprintf(fp_debug, "\tReading %s\n", filename);
		fprintf(fp_debug, "\tn_monomers = %d\n", (*n_monomer) );
		for (i = 0; i < (*n_monomer); i++)
		{
			fprintf(fp_debug, "\t\tn_atoms = %d, %s\n", (*n_atom)[i], (*conn_files)[i]);
			for (j = 0; j < (*n_atom)[i]; j++)
				fprintf(fp_debug, "\t\t\t%10.5f %10.5f %10.5f\n", (*R2D_monomer)[i][j][0],
				(*R2D_monomer)[i][j][1], (*R2D_monomer)[i][j][2]);
		}
	}

/*Cluster Geometries, Energies, and Forces*/
	rewind(fp);
	Go_To_Line (fp, "@Cluster_Geometries", filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%*s %*s %d", &(*n_clusters) );
	
	(*QM_Energy)	= D_Allocate_1D_Matrix ( (*n_clusters));
	(*R2D_cluster)	= (double ****) malloc ( (*n_clusters)*sizeof(double ***) );
	(*QM_Force)		= (double ****) malloc ( (*n_clusters)*sizeof(double ***) );

	for (n = 0; n < (*n_clusters); n++)
	{
		(*R2D_cluster)[n]	= (double ***) malloc ( (*n_monomer)*sizeof(double **) );
		(*QM_Force)[n]		= (double ***) malloc ( (*n_monomer)*sizeof(double **) );
		for (i = 0; i < (*n_monomer); i++)
		{
			(*R2D_cluster)[n][i]	= D_Allocate_2D_Matrix ( (*n_atom)[i], 3);
			(*QM_Force)[n][i]		= D_Allocate_2D_Matrix ( (*n_atom)[i], 3);
		}
	}

	for (n = 0; n < (*n_clusters); n++)
	{
		Get_Next_Line (fp, str, filename);
		sscanf(str, "%*s %*d %*s %*s %lf", &( (*QM_Energy)[n]));
		if ( ( (*QM_Energy)[n] > E_MAX) || ( (*QM_Energy)[n] < E_MIN) )
			printf("WARNING for Cluster %d E = %f in %s\n", n, (*QM_Energy)[n], filename);

		for (i = 0; i < (*n_monomer); i++)
			for (j = 0; j < (*n_atom)[i]; j++)
			{
				Get_Next_Line (fp, str, filename);
				sscanf(str, "%d %d %lf %lf %lf %lf %lf %lf", &i_temp1, &i_temp2,
					&( (*R2D_cluster)[n][i][j][0]), &( (*R2D_cluster)[n][i][j][1]),
					&( (*R2D_cluster)[n][i][j][2]), &( (*QM_Force)[n][i][j][0]),
					&( (*QM_Force)[n][i][j][1]), &( (*QM_Force)[n][i][j][2]) );

				if ( (i_temp1 != i) || (i_temp2 != j) )
				{
					printf("ERROR reading in %s cluster geometry %d\n", filename, n);
					printf("i_temp1 = %d != %d, itemp2 = %d != %d\n", i_temp1, i, i_temp2, j);
					exit(0);
	
				}
			}
	}
	fclose(fp);
					
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\tn_clusters = %d\n", (*n_clusters) );
		for (n = 0; n < (*n_clusters); n++)
		{
			fprintf(fp_debug, "\t\tE[%4d] = %f\n", n, (*QM_Energy)[n]);
			for (i = 0; i < (*n_monomer); i++)
				for (j = 0; j < (*n_atom)[i]; j++)
					fprintf(fp_debug, "\t\t\t%4d %4d %10.5f %10.5f %10.5f		%10.5f %10.5f %10.5f\n",
						i, j, (*R2D_cluster)[n][i][j][0], (*R2D_cluster)[n][i][j][1],
						(*R2D_cluster)[n][i][j][2], (*QM_Force)[n][i][j][0],
						(*QM_Force)[n][i][j][1], (*QM_Force)[n][i][j][2]);
			fprintf(fp_debug, "\n");
		}
	}
}





void Read_Conn_Files (int *n_atoms, int **atom_types, double **charges,
		double **polar, int **n_neighbors, int ***neighbors, int n_atom_types, 
		char **atom_type_names, int s_atom_type_size, char *filename, 
		FILE *fp_debug, int DEBUG)
/*Reads Connectivity File in Following Format :
n_atoms = 6
           #  at   chg  n_neigh    neigh
           0   CT    0.24417    0.73630      4             1    2    3    4
           1   H1    0.00000    0.16570      1             0
*/
{
	int i, j;
	char str[500];
	FILE *fp;
	int neighbor_temp[5];
	char s_atom_type[10];


	fp = fopen( filename, "r");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	if (fgets(str, sizeof(str), fp) != NULL)
		sscanf(str, "%*s %*s %d", &(*n_atoms) );
	fgets(str, sizeof(str), fp);

	(*atom_types)	= I_Allocate_1D_Matrix ( (*n_atoms) );
	(*charges)		= D_Allocate_1D_Matrix ( (*n_atoms) );
	(*polar)		= D_Allocate_1D_Matrix ( (*n_atoms) );
	(*n_neighbors)	= I_Allocate_1D_Matrix ( (*n_atoms) );
	(*neighbors)	= (int **) malloc ( (*n_atoms)*sizeof(int *) );


	for (i = 0; i < (*n_atoms); i++)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%*d %s %lf %lf %d %d %d %d %d %d", s_atom_type,
				&( (*charges)[i]), &( (*polar)[i]), &( (*n_neighbors)[i] ), 
				&neighbor_temp[0], &neighbor_temp[1], &neighbor_temp[2], 
				&neighbor_temp[3], &neighbor_temp[4]);

			(*neighbors)[i] = I_Allocate_1D_Matrix ( (*n_neighbors)[i] );
			for (j = 0; j < (*n_neighbors)[i]; j++)
				(*neighbors)[i][j] = neighbor_temp[j];
			s_atom_type[s_atom_type_size] = '\0';
			(*atom_types)[i] = Find_Atom_Type (s_atom_type, n_atom_types, atom_type_names);
		}
		else
		{
			printf("Enexpected EOF while reading %s\n", filename);
			exit(0);
		}
	}

	fclose(fp);

	if (DEBUG)
	{
		fprintf(fp_debug, "Reading %s\n", filename);
		fprintf(fp_debug, "\t(*n_atoms) = %d\n", (*n_atoms));
		for (i = 0; i < (*n_atoms); i++)
		{
			fprintf(fp_debug, "\t\t%4d %4s %10.5f %10.5f %4d	", i, 
				atom_type_names[(*atom_types)[i] ], (*charges)[i],
				(*polar)[i], (*n_neighbors)[i] );
			for (j = 0; j < (*n_neighbors)[i]; j++)
				fprintf(fp_debug, "%d ", (*neighbors)[i][j] );
			fprintf(fp_debug, "\n");
		}
	}
}



void Read_XYZ_File (char *filename, char *flag, int *n_atoms, int **mol_atom_types, 
		double ***mol_r2d, int **n_neighbors, int ***neighbors, int n_atom_types, 
		char **atom_type_names, FILE *fp_debug, int DEBUG)
/*
flag
     8
     1 CT  C   -1.594800    0.078100    0.000000     4     2     3     4     5 
     2 CT  C   -0.060500    0.215100   -0.000000     4     1     6     7     8 
*/
{
	int i, j;
	int neighbors_temp[10];
	char s_atom_type[10];
	char str[500];
	FILE *fp;


	fp = fopen( filename, "r");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	Go_To_Line (fp, flag, filename);

	fgets(str, sizeof(str), fp);
	sscanf(str, "%*s %d", &(*n_atoms) );

	
	(*mol_atom_types)	= I_Allocate_1D_Matrix (*n_atoms);
	(*mol_r2d)			= D_Allocate_2D_Matrix (*n_atoms, 3);
	(*n_neighbors)		= I_Allocate_1D_Matrix (*n_atoms);
	(*neighbors)		= (int **) malloc ( (*n_atoms)*sizeof(int *) );
	
	for (i = 0; i < *n_atoms; i++)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%*d %s %*s %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d", s_atom_type,
				&( (*mol_r2d)[i][0]), &( (*mol_r2d)[i][1]), &( (*mol_r2d)[i][2]),
				&( (*n_neighbors)[i]), &neighbors_temp[0], &neighbors_temp[1],
				&neighbors_temp[2], &neighbors_temp[3], &neighbors_temp[4], &neighbors_temp[5],
				&neighbors_temp[6], &neighbors_temp[7], &neighbors_temp[8], &neighbors_temp[9]);
			(*neighbors)[i] = I_Allocate_1D_Matrix ( (*n_neighbors)[i] );
			for (j = 0; j < (*n_neighbors)[i]; j++)
				(*neighbors)[i][j] = neighbors_temp[j];

			(*mol_atom_types)[i] = Find_Atom_Type (s_atom_type, n_atom_types, atom_type_names);
		}
		else
		{
			printf("Unexpected EOF while reading %s\n", filename);
			exit(0);
		}
	}

	fclose(fp);
	
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nReading %s\n", filename);
		fprintf(fp_debug, "\tn_atoms = %d\n", *n_atoms);
		for (i = 0; i < *n_atoms; i++)
		{
			fprintf(fp_debug, "\t\t%4d %4s %10.5f %10.5f %10.5f\t\t%4d ", i, 
				atom_type_names[ (*mol_atom_types)[i] ], (*mol_r2d)[i][0],
				(*mol_r2d)[i][1], (*mol_r2d)[i][2], (*n_neighbors)[i]);

			for (j = 0; j < (*n_neighbors)[i]; j++)
				fprintf(fp_debug, "%4d ", (*neighbors)[i][j] );
			fprintf(fp_debug, "\n");
		}
	}
}


void Read_AT2_XYZ_File (char *filename, int *n_atoms, int **atom_types1, int **atom_types2,
		double ***mol_r2d, int **n_neighbors, int ***neighbors, int n_atom_types, 
		char **atom_type_names, FILE *fp_debug, int DEBUG)
/*
     8
     1 CT	C     -1.594800    0.078100    0.000000     4     2     3     4     5 
     2 CT	C     -0.060500    0.215100   -0.000000     4     1     6     7     8 
*/
{
	int i, j;
	int neighbors_temp[10];
	char s_atom_type1[10], s_atom_type2[10];
	char str[500];
	FILE *fp;


	fp = fopen( filename, "r");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	fgets(str, sizeof(str), fp);
	sscanf(str, "%d", &(*n_atoms) );

	
	(*atom_types1)		= I_Allocate_1D_Matrix (*n_atoms);
	(*atom_types2)		= I_Allocate_1D_Matrix (*n_atoms);
	(*mol_r2d)			= D_Allocate_2D_Matrix (*n_atoms, 3);
	(*n_neighbors)		= I_Allocate_1D_Matrix (*n_atoms);
	(*neighbors)		= (int **) malloc ( (*n_atoms)*sizeof(int *) );
	
	for (i = 0; i < *n_atoms; i++)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%*d %s %s %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d", 
				s_atom_type1, s_atom_type2, 
				&( (*mol_r2d)[i][0]), &( (*mol_r2d)[i][1]), &( (*mol_r2d)[i][2]),
				&( (*n_neighbors)[i]), &neighbors_temp[0], &neighbors_temp[1],
				&neighbors_temp[2], &neighbors_temp[3], &neighbors_temp[4], &neighbors_temp[5],
				&neighbors_temp[6], &neighbors_temp[7], &neighbors_temp[8], &neighbors_temp[9]);
			(*neighbors)[i] = I_Allocate_1D_Matrix ( (*n_neighbors)[i] );
			for (j = 0; j < (*n_neighbors)[i]; j++)
				(*neighbors)[i][j] = neighbors_temp[j];

			(*atom_types1)[i] = Find_Atom_Type (s_atom_type1, n_atom_types, atom_type_names);
			(*atom_types2)[i] = Find_Atom_Type (s_atom_type2, n_atom_types, atom_type_names);
		}
		else
		{
			printf("Unexpected EOF while reading %s\n", filename);
			exit(0);
		}
	}

	fclose(fp);
	
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nReading %s\n", filename);
		fprintf(fp_debug, "\tn_atoms = %d\n", *n_atoms);
		for (i = 0; i < *n_atoms; i++)
		{
			fprintf(fp_debug, "\t\t%4d %4s %4s %10.5f %10.5f %10.5f\t\t%4d ", i, 
				atom_type_names[ (*atom_types1)[i] ], atom_type_names[ (*atom_types2)[i] ],
				(*mol_r2d)[i][0], (*mol_r2d)[i][1], (*mol_r2d)[i][2], (*n_neighbors)[i]);

			for (j = 0; j < (*n_neighbors)[i]; j++)
				fprintf(fp_debug, "%4d ", (*neighbors)[i][j] );
			fprintf(fp_debug, "\n");
		}
	}
}

void Read_Pdb_File (int *n_atoms, double ***R_2D, char ***atom_type, 
		char *filename, FILE *fp_debug, int DEBUG)
/*
ATOM      1  H1  OME     1       5.322  16.919  10.213
ATOM      2  CH3 OME     1       5.853  15.985  10.410
*/
{
	int i;
	FILE *fp;
	char str[1000];
	char s_temp[50];

	fp = fopen( filename, "r");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

/*Find N_atoms*/
	(*n_atoms) = 0;
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s ", s_temp);
			if (strcmp(s_temp, "ATOM") == 0)
				(*n_atoms)++;
		}
		else
			break;
	}

	if (DEBUG)
		fprintf(fp_debug, "n_atoms = %d\n", *n_atoms);
	rewind(fp);

	
	(*R_2D)			= D_Allocate_2D_Matrix (*n_atoms, 3);
	(*atom_type) 	= C_Allocate_2D_Matrix (*n_atoms, 10);

	i = 0;
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s ", s_temp);
			if (strcmp(s_temp, "ATOM") == 0)
			{
				sscanf(str, "%*s %*d %s %*s %*d %lf %lf %lf",
					(*atom_type)[i], &( (*R_2D)[i][0]), &( (*R_2D)[i][1]), &( (*R_2D)[i][2]) );
				i++;
			}
		}
		else
			break;
	}

	fclose(fp);
	
	if (DEBUG)
	{
		fprintf(fp_debug, "\nReading %s\n", filename);
		fprintf(fp_debug, "N_atoms = %d\n", *n_atoms);
		for (i = 0; i < (*n_atoms); i++)
			fprintf(fp_debug, "\t%4d %4s %10.5f %10.5f %10.5f\n", i, (*atom_type)[i], 
				(*R_2D)[i][0], (*R_2D)[i][1], (*R_2D)[i][2]);
	}
}








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
	double *charge, double *polar, int *n_neighbors, int **neighbors, char *filename)
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
{
	int i, j;
	FILE *fp;

	fp = fopen( filename, "w");
	if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "n_atoms = %d\n", n_atoms);
	fprintf(fp, "			#	at		chg		polar	n_neigh    		neigh\n");
	for (i = 0; i < n_atoms; i++)
	{
		fprintf(fp, "\t\t%4d %4s %10.5f %10.5f\t%4d\t\t", i, AT_Elements[ element_num[i] ],
			charge[i], polar[i], n_neighbors[i]);
		for (j = 0; j < n_neighbors[i]; j++)
			fprintf(fp, "\t%3d", neighbors[i][j]);
		fprintf(fp, "\n");
	}

	fclose(fp);
}


void Print_XYZ_Format (int num_atoms, double **r_2D_atoms, 
		int *n_neighbors, int **neighbors, int *atom_type,
		int n_atom_types, char **atom_type_names, char *filename, 
		FILE *fp_debug, int DEBUG)
/*
  1112
     0 CT     -2.852565    7.129814    7.213353     1     	1     2    3    4
     1 HC     -2.010039    6.775897    7.848916     3     	0
     2 HC     -3.744861    7.281564    7.792989     4     	0
*/
{
	int i, j;
	FILE *fp;

	fp = fopen( filename, "w");
    if (fp == NULL) 
	{
		printf("Cant open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "%6d\n", num_atoms);
	for (i = 0; i < num_atoms; i++)
	{
		fprintf(fp, "%6d %2s %13.6f %11.6f %11.6f %5d ", i, atom_type_names[atom_type[i]],
			r_2D_atoms[i][0], r_2D_atoms[i][1], r_2D_atoms[i][2], 
			n_neighbors[i]);
		for (j = 0; j < n_neighbors[i]; j++)
			fprintf(fp, "%5d ", neighbors[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void Make_Pdb_File (int num_atoms, char **element, double **r_2d, char *Filename) 
{ 
	int i; 
    FILE *fp;

	fp = fopen( Filename, "w"); 
    if (fp == NULL)   
	{ 
		printf("Cant open %s\n", Filename); 
		exit(EXIT_FAILURE); 
	} 
 
	for (i = 0; i < num_atoms; i++) 
		fprintf(fp, "ATOM  %5d %2s   XXX     1    %8.3f %7.3f %7.3f\n",  
			    i+1, element[i], r_2d[i][0], r_2d[i][1], r_2d[i][2]); 
	fprintf(fp, "\n\n");
	fclose(fp); 
}


void Write_Coordinates (FILE *fp_crd, int num_atoms, double *r_crd)
{
	int i;
	int n, r, num_coords;

	num_coords = 3*num_atoms;
	r = num_coords%10;
	n = (num_coords - r)/10;

	for (i = 0; i < n; i++)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", 
			r_crd[10*i  ], r_crd[10*i+1], r_crd[10*i+2], r_crd[10*i+3], 
            r_crd[10*i+4], r_crd[10*i+5], r_crd[10*i+6], r_crd[10*i+7], 
            r_crd[10*i+8], r_crd[10*i+9]);

	if (r == 9)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2], r_crd[10*n+3], 
            r_crd[10*n+4], r_crd[10*n+5], r_crd[10*n+6], r_crd[10*n+7], 
            r_crd[10*n+8]);

	else if (r == 8)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2], r_crd[10*n+3], 
            r_crd[10*n+4], r_crd[10*n+5], r_crd[10*n+6], r_crd[10*n+7]);

	else if (r == 7)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2], r_crd[10*n+3], 
            r_crd[10*n+4], r_crd[10*n+5], r_crd[10*n+6]);

	else if (r == 6)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2], r_crd[10*n+3], 
            r_crd[10*n+4], r_crd[10*n+5]);

	else if (r == 5)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f %7.3f \n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2], r_crd[10*n+3], 
            r_crd[10*n+4]);

	else if (r == 4)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f %7.3f \n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2], r_crd[10*n+3]);

	else if (r == 3)
		fprintf(fp_crd, " %7.3f %7.3f %7.3f\n", 
			r_crd[10*n  ], r_crd[10*n+1], r_crd[10*n+2]);

	else if (r == 2)
		fprintf(fp_crd, " %7.3f %7.3f\n", 
			r_crd[10*n  ], r_crd[10*n+1]);

	else if (r == 1)
		fprintf(fp_crd, " %7.3f\n", 
			r_crd[10*n]);
/*
	fprintf(fp_crd, " %7.3f %7.3f %7.3f\n", dim_x, dim_y, dim_z);
*/

}


void Print_Geometry (char *title, FILE *fp, int n_atom, char **atom_names,
	double **r2d, int *n_neigh, int **neigh)
/*
@Geometry
	n_atoms 4
		H1		0.28827   0.0000       -2.0037 		1	1
		O2		0.87274   0.0000       -1.24675		2	0 2
		H3		0.28827   0.0000       -0.4898 		1	1
		X		0.28827   0.0000       -1.24675		0
*/
{
	int i, j;

	fprintf(fp, "%s", title);
	fprintf(fp, "\tn_atoms %d\n", n_atom);
	for (i = 0; i < n_atom; i++)
	{
		fprintf(fp, "\t\t%4s %12.8f %12.8f %12.8f %4d ", atom_names[i],
			r2d[i][0], r2d[i][1], r2d[i][2], n_neigh[i]);
		for (j = 0; j < n_neigh[i]; j++)
			fprintf(fp, "%4d ", neigh[i][j]);
		fprintf(fp, "\n");
	}
	fflush(fp);
}

void Print_Double (char *title, FILE *fp, int n_value, char **atom_names, double *Y)
/*
@Z
	n_atoms 4
		H1	1.0
		O2	6.0
		H3	1.0
		X	0.0
*/
{
	int i;

	fprintf(fp, "%s", title);
	fprintf(fp, "\tn_value %d\n", n_value);
	for (i = 0; i < n_value; i++)
		fprintf(fp, "\t\t%4s %12.8f\n", atom_names[i], Y[i]);
	fflush(fp);
}






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



void IO_Error (char *filename)
{
	printf("ERROR unexpected EOF while reading in %s\n", filename);
	exit(0);
}

FILE *Open_File (char *filename, char *action)
{
	FILE *fp;

	fp = fopen(filename, action); 
	if (fp == NULL)  
	{ 
		printf("Cant open %s\n", filename); 
		exit(EXIT_FAILURE); 
	}

	return fp;
}




int string_compare(char *s1, char *s2, int N)
/*checks if first N characters of 2 strings are the same*/
{
    int i;

    for (i = 0; i < N; i++)	
	{
		if (s1[i] == '\0')  
		{
			if (s2[i] == '\0') 
				return 1;
			else if (s2[i] == ' ')  
                return 1;
            else
				return 0;
		}
		if (s2[i] == '\0')  
		{
			if (s1[i] == '\0') 
				return 1;
            else if (s1[i] == ' ')  
                return 1;
			else 
				return 0;
		}
        if (s1[i] != s2[i])
			return 0;
	}
    return 1;
}

int Is_New_Line(char *str)
{
	int i;
    for (i = 0; ;i++)
    {
        if (str[i] == '\n')
            return 1;
        if ( (str[i] != ' ') && (str[i] != '\t') && ( (int)str[i] != 13) )
            return 0;
		if (str[i] == '\0')
			return 1;
    }
    return 0;
}

void Get_Next_Line (FILE *fp, char *str_return, char *filename)
{
	char str[3000];
    for (;;)
    {
        if (fgets(str, sizeof(str), fp) != NULL)
		{
			if ( (str[0] != '/') &&(Is_New_Line(str) == 0) )
			{
				sprintf(str_return, "%s", str);
				return;
			}
		}
		else
			IO_Error(filename);
    }
}


void Go_To_Line (FILE *fp, char *flag, char *filename)
/*skips to line where the flag char string is found*/
{
    char str[2000], f_str[200];

    for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
	        sscanf(str, "%s", f_str);
    	    if ( strcmp(f_str, flag) == 0 )
    	        return;
		}
		else
		{
			printf("Could not find %s in file %s\n", flag, filename);
			exit(0);
		}
	}
}

void Read_Value (FILE *fp, char *flag, char *filename, char *scan_str, char *return_str)
/*looks for:
scan_str return_str
if it finds scan_str, then it returns return_str, if not return_str = NULL
*/
{
	char str[4000];
	char str_temp[4000];

	return_str[0] = '\0';
	rewind(fp);
	Go_To_Line (fp, flag, filename);
	for (;;)
	{
		if (fgets(str, sizeof(str), fp) != NULL)
		{
			sscanf(str, "%s ", str_temp);
			if (strcmp(str_temp, scan_str) == 0)
			{
				sscanf(str, "%*s %s", return_str);
				return;
			}
		}
		else
			return;
	}
}

int Find_Atom_Type (char *atom_type_s, int Num_atom_types, char **Atom_Type_Names)
{
	int i;
	for (i = 0; i < Num_atom_types; i++)
		if ( strcmp(atom_type_s, Atom_Type_Names[i]) == 0 )
			return i;
	printf("Could not find atom type number for %s\n", atom_type_s);
	exit(0);
}

void Get_Current_Directory (char *directory)
{
	FILE *fp;
	char str[2000];

	system("pwd > directory.txt");
	fp = fopen( "directory.txt", "r");
	if (fp == NULL) 
	{
		printf("Cant open %s\n", "directory.txt");
		exit(EXIT_FAILURE);
	}
	fgets(str, sizeof(str), fp);
	sscanf(str, "%s", directory);
	fclose(fp);	
	system("rm directory.txt");
}

void Change_File_Extension (const char *s1, const char *s2, char *s3)
/*attaches s2 to s1 beginning with at the '.' in s1, returns s3*/
{
	int i, j;

	for (i = 0; ; i++)
	{
		if ( (s1[i] == '.') || (s1[i] == '\0')	)
			break;
		s3[i] = s1[i];
	}
	
	
	for (j = 0; ; j++)
	{	
		s3[j+i] = s2[j];
		if (s2[j] == '\0')
			break;
	}
}		

void Remove_Extension (char *s1, char *s2)
/*
        s1 = test.dat => s2 = test
*/
{
        int i;
        i = 0;
        for (;;)
        {
                if ( (s1[i] == '\0') || (s1[i] == '.') )
                        break;
                s2[i] = s1[i];
                i++;
        }
        s2[i] = '\0';
}


void Get_Extension (char *s1, char *s2)
/*
        s1 = test.dat => s2 = dat
*/
{
	int i, start;

	for (i = 0; ;i++)
	{
		if (s1[i] == '.')
			break;
		if (s1[i] == '\0')
		{
			s2[0] = '\0';
			return;
		}
	}

	start = i+1;

	for (i = 0; ;i++)
	{
		s2[i] = s1[i+start];
		if (s2[i] == '\0')
			break;
	}


}

int Does_File_Exists (char *filename)
{
	FILE *fp;

	fp = fopen(filename, "r"); 
	if (fp == NULL)  
	{ 
		return 0;
	}
	else
	{
		fclose(fp);
		return 1;
	}
}




void Create_Sub_Directory (char *directory_name)
{
	char command[2000];

	if (Does_File_Exists(directory_name) )
	{
		sprintf(command, "rm -r %s", directory_name);
		system(command);
	}

	sprintf(command, "mkdir %s", directory_name);
	system(command);
}







