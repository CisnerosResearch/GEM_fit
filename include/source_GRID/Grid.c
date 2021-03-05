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
#include "Lebedev.h"

#define Pi 3.1415926535897932384626433832
#define N_P_ITERATION 3
#define MAX_RAD 20.0
#define GRID_WEIGHT_TOL 1.0E-15

int GRID_n_ang_pt;
double *GRID_ang_weight;
double **GRID_ang_r2d;

/*Gauss Chebyshev points*/
int GRID_n_GC2_pt;
double *GRID_GC2_xi;
double *GRID_GC2_wi;

struct GRID_Data 
{
	int n_element;
	char **element;
	double *bragg_radii;

	int ang_index1;
	int ang_index2;
	int n_radial_pt;
};
struct GRID_Data GRID;




void Read_Grid_Instruct (char *filename, struct GRID_Data *grid, FILE *fp_debug, int DEBUG);
double Iterate_p_u (double mu, int N);
double Calc_s_u (double mu, int N);
void Set_Gauss_Chebyshev_Grid (int n, double **xi, double **wi, int DEBUG, FILE *fp_debug);
void Make_Atomic_Grid (double rm, double *r_atom, int *n_pt, double ***r2d, double **weight);
/*rm is half of the Bragg radius except for hydrogen*/
double Calc_Atom_Rel_Weight (int atom_i, int n_atom, double **r2d_atom, double *r_grid, double *bragg_radii);
void Sort_Array (int n, double *r, double *w);
/*sorts in order of increasing r*/




void Initialize_Grid (char *grid_data_file, FILE *fp_debug, int DEBUG)
/*Reads in grid_data_file containing grid parameters: elements, Slater-Bragg radii for each element,
  Lebedev order, and radial order.  It then creates a single Levedev angular grid and a single Gauss Chebyshev radial grid.
*/
{
	struct LEB_Data leb;

	Read_Grid_Instruct (grid_data_file, &GRID, fp_debug, DEBUG);
	/*Reads in grid_data_file containing grid parameters: elements, Slater-Bragg radii for each element,
	  Lebedev order, and radial order*/

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nInitialize_Grid\n");
		fprintf(fp_debug, "Setting Angular Grid\n");
	}


	Set_Lebedev_Grid_Parm (GRID.ang_index1, GRID.ang_index2, &leb, DEBUG, fp_debug);
	/*Defines Lebedev Angular Grid Parameters in terms of constants defined in Lebedev's papers*/

	Calc_Lebedev_Grid (&GRID_n_ang_pt, &GRID_ang_r2d, &GRID_ang_weight, &leb);
	/*Calculates Lebedev Angular Grid*/

	GRID_n_GC2_pt = GRID.n_radial_pt;
	/*Initialize Gauss Chebyshev Grid*/

	Set_Gauss_Chebyshev_Grid (GRID_n_GC2_pt, &GRID_GC2_xi, &GRID_GC2_wi, DEBUG, fp_debug);
	/*Calculates Gauss Chebyshev Radial Grid*/
}

void Read_Grid_Instruct (char *filename, struct GRID_Data *grid, FILE *fp_debug, int DEBUG)
/*Reads in grid_data_file containing grid parameters: elements, Slater-Bragg radii for each element,
  Lebedev order, and radial order*/
{
	FILE *fp;
	char str[2000];
	char s_temp[500];

	Read_Atom_Types ( &( (*grid).n_element), &( (*grid).element), filename, "@Element", fp_debug, DEBUG);
	Read_AT_Double  ( &( (*grid).bragg_radii), (*grid).n_element, (*grid).element, filename, "@Bragg_Radii", fp_debug, DEBUG);

	fp = Open_File (filename, "r");
	Go_To_Line (fp, "@Instruct", filename);
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%s %d %d", s_temp, &( (*grid).ang_index1), &( (*grid).ang_index2) );
	if (strcmp(s_temp, "Lebedev_Ang_Type") != 0)
	{
		printf("ERROR reading %s.. expecting Lebedev_Ang_Type\n", filename);
		exit(0);
	}
	Get_Next_Line (fp, str, filename);
	sscanf(str, "%s %d", s_temp, &( (*grid).n_radial_pt) );
	if (strcmp(s_temp, "N_Radial_Pt") != 0)
	{
		printf("ERROR reading %s.. expecting N_Radial_Pt\n", filename);
		exit(0);
	}
	fclose(fp);


	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nReading @Instruct\n");
		fprintf(fp_debug, "\tLebedev_Ang_Type %d %d\n", (*grid).ang_index1, (*grid).ang_index2);
		fprintf(fp_debug, "\tN_Radial_Pt      %d\n",    (*grid).n_radial_pt);
		fflush(fp_debug);
	}

}

void Make_Molecular_Grid (int n_atom, double **mol_r2d, char **mol_element, 
	int *tot_n_pt, double **tot_weight, double ***tot_r2d, FILE *fp_debug, int DEBUG)
/*Once the angular and radial grids have been created in Initialize_Grid, this function will create a molecular grid
  using the Slater-Bragg radii and Becke's smoothing function method.
*/
{
	int i, j, m, p, count;
	double rm;
	double sum;
	double delta_r[3];
	
/*Atomic Grids*/
	int *ag_n_pt;			/*[n_atom]*/
	double **ag_weight;		/*[n_atom][ag_n_pt]*/
	double ***ag_r2d;		/*[n_atom][ag_n_pt][3]*/

/*smooth function weight*/
	double **sm_weight;		/*[n_atom][ag_n_pt]*/

/*local bragg radii*/
	int ele_num;
	double *bragg_radii;


/*Atomic Grids*/
	ag_n_pt   = (int *)      malloc ( n_atom*sizeof(int) );
	ag_weight = (double **)  malloc ( n_atom*sizeof(double *) );
	ag_r2d    = (double ***) malloc ( n_atom*sizeof(double **) );

	bragg_radii = (double *)  malloc ( n_atom*sizeof(double ) );
	for (i = 0; i < n_atom; i++)
	{
		ele_num = Find_Atom_Type (mol_element[i], GRID.n_element, GRID.element);
		bragg_radii[i] = GRID.bragg_radii[ele_num];
	}

/*Make Atomic Grids*/
	for (i = 0; i < n_atom; i++)
	{
		if (strcmp(mol_element[i], "H") == 0)
			rm = bragg_radii[i];
		else
			rm = bragg_radii[i]/2.0;
		Make_Atomic_Grid (rm, mol_r2d[i], &(ag_n_pt[i]), &(ag_r2d[i]), &(ag_weight[i]));
	}

	sm_weight = (double **) malloc (n_atom*sizeof(double *) );
	for (i = 0; i < n_atom; i++)
	{
		sm_weight[i] = D_Allocate_1D_Matrix (ag_n_pt[i]);
		for (m = 0; m < ag_n_pt[i]; m++)
		{
			sm_weight[i][m] = Calc_Atom_Rel_Weight (i, n_atom, mol_r2d, ag_r2d[i][m], bragg_radii);
			/*Calculate relative weight from other atoms*/
			sum = sm_weight[i][m];
			for (j = 0; j < n_atom; j++)
				if (i != j)
					sum += Calc_Atom_Rel_Weight (j, n_atom, mol_r2d, ag_r2d[i][m], bragg_radii);
			sm_weight[i][m] /= sum;
		}
	}


	if (DEBUG)
	{
		fprintf(fp_debug, "sm_weight\n");
		for (i = 0; i < n_atom; i++)
			for (m = 0; m < ag_n_pt[i]; m++)
			{
				for (p = 0; p < 3; p++)
					delta_r[p] = mol_r2d[i][p] - ag_r2d[i][m][p];
				rm = sqrt(dot_prod (&(delta_r[0]), &(delta_r[0]) ) );

				fprintf(fp_debug, "atom %5d grid_pt %5d r = %10.6f  weight = %10.6f\n", i, m, rm, sm_weight[i][m]);
				for (j = 0; j < n_atom; j++)
					if (i != j)
					{
						for (p = 0; p < 3; p++)
							delta_r[p] = mol_r2d[j][p] - ag_r2d[i][m][p];
						rm = sqrt(dot_prod (&(delta_r[0]), &(delta_r[0]) ) );
						fprintf(fp_debug, "\tr_%d = %8.4f\n", j, rm);
					}
				fprintf(fp_debug, "\n");
			}
	}


	(*tot_n_pt) = 0;
	for (i = 0; i < n_atom; i++)
		for (m = 0; m < ag_n_pt[i]; m++)
			if (fabs(sm_weight[i][m]*ag_weight[i][m]) > GRID_WEIGHT_TOL)
				(*tot_n_pt)++;

	(*tot_weight) = D_Allocate_1D_Matrix ((*tot_n_pt));
	(*tot_r2d)    = D_Allocate_2D_Matrix ((*tot_n_pt), 3);

	count = 0;
	for (i = 0; i < n_atom; i++)
		for (m = 0; m < ag_n_pt[i]; m++)
			if (fabs(sm_weight[i][m]*ag_weight[i][m]) > GRID_WEIGHT_TOL)
			{
				for (p = 0; p < 3; p++)
					(*tot_r2d)[count][p] = ag_r2d[i][m][p];
				(*tot_weight)[count] = sm_weight[i][m]*ag_weight[i][m];
				count++;
			}
		
	


/*free atomic grids*/
	for (i = 0; i < n_atom; i++)
	{
		for (j = 0; j < ag_n_pt[i]; j++)
			free(ag_r2d[i][j]);
		free(ag_r2d[i]);
		free(ag_weight[i]);
		free(sm_weight[i]);
	}
	free(ag_n_pt);
	free(ag_weight);
	free(ag_r2d);
	free(sm_weight);
	free(bragg_radii);


}

double Calc_Atom_Rel_Weight (int atom_i, int n_atom, double **r2d_atom, double *r_grid, double *bragg_radii)
{
	double rel_weight;
	int j, p;
	double ri, rj, Rij;
	double mu_ij, chi, u_ij, a_ij, v_ij;
	double delta_r[3];

	rel_weight = 1.0;
	for (j = 0; j < n_atom; j++)
		if (atom_i != j)
		{
			for (p = 0; p < 3; p++)
				delta_r[p] = r2d_atom[atom_i][p] - r_grid[p];
			ri = sqrt(dot_prod (&(delta_r[0]), &(delta_r[0]) ) );

			for (p = 0; p < 3; p++)
				delta_r[p] = r2d_atom[j][p] - r_grid[p];
			rj = sqrt(dot_prod (&(delta_r[0]), &(delta_r[0]) ) );

			for (p = 0; p < 3; p++)
				delta_r[p] = r2d_atom[atom_i][p] - r2d_atom[j][p];
			Rij = sqrt(dot_prod (&(delta_r[0]), &(delta_r[0]) ) );

			mu_ij = (ri-rj)/Rij;
			chi   = bragg_radii[atom_i]/bragg_radii[j];
			u_ij  = (chi-1.0)/(chi+1.0);
			a_ij  = u_ij/(u_ij*u_ij-1.0);
			if (a_ij >  0.5)	a_ij =  0.5;
			if (a_ij < -0.5)	a_ij = -0.5;
			v_ij  = mu_ij + a_ij*(1.0-mu_ij*mu_ij);


			rel_weight *= Calc_s_u (v_ij, N_P_ITERATION);
		}
	return rel_weight;
}


void Sort_Array (int n, double *r, double *w)
/*sorts in order of increasing r*/
{
	int Is_Done, i;
	double temp1, temp2;

	for (;;)
	{
		Is_Done = 1;
		for (i = 0; i < n-1; i++)
		{
			if (r[i+1] < r[i])
			{
				temp1  = r[i];
				temp2  = r[i+1];
				r[i]   = temp2;
				r[i+1] = temp1;

				temp1  = w[i];
				temp2  = w[i+1];
				w[i]   = temp2;
				w[i+1] = temp1;
				Is_Done = 0;
			}
		}
		if (Is_Done)
			return;
	}
}


double Calc_s_u (double mu, int N)
{
	return 0.5*(1.0-Iterate_p_u (mu, N) );
}


double Iterate_p_u (double mu, int N)
{
	double v;

	if (N == 1)
		return 1.5*mu - 0.5*mu*mu*mu;
	else
	{
		v = Iterate_p_u (mu, N-1);
		return 1.5*v - 0.5*v*v*v;
	}
}




void Debug_Make_Grid (int DEBUG, int PRINT_GRID, FILE *fp_debug)
/*assuming Initialize_Grid is already set*/
{
	int i;
	double sum, f;

/*atom grid*/
	double rm;
	double r_atom[3];
	int n_pt;
	double **r2d;
	double *weight;
	double a;
	double x, y, z, r2;


	if (DEBUG == 0)
		return;

/*Gauss Chebyshev Grid of the Second Kind*/
/*check integral of (1-x^2)*3x^2 from -1 to 1 (should be 0.8)*/
	sum = 0;
	for (i = 0; i < GRID_n_GC2_pt; i++)
	{
		x = GRID_GC2_xi[i];
		f = sqrt(1.0-x*x)*3*x*x;
		sum += GRID_GC2_wi[i]*f;
	}
	fprintf(fp_debug, "integral of (1-x^2)*3x^2 from -1 to 1 = %15.12f (NUME)\n", sum);
	fprintf(fp_debug, "integral of (1-x^2)*3x^2 from -1 to 1 = %15.12f (ANAL)\n", 0.8);


	fprintf(fp_debug, "\n\nMake_Atomic_Grid\n");
	rm = 0.5;
	r_atom[0] = 0.0;	r_atom[1] = 0.0;	r_atom[2] = 0.0;
	Make_Atomic_Grid (rm, &(r_atom[0]), &n_pt, &r2d, &weight);

	fprintf(fp_debug, "n_pt = %d\n", n_pt);
	if (PRINT_GRID)
		for (i = 0; i < n_pt; i++)
			fprintf(fp_debug, "\tr2d = %10.6f %10.6f %10.6f    weight = %15.12f\n", r2d[i][0], r2d[i][1], r2d[i][2], weight[i]);

/*test on integral of exp(-ar*r)*(x^2+y^2) over all space*/
	a = 2.5;
	sum = 0;
	for (i = 0; i < n_pt; i++)
	{
		x = r2d[i][0];
		y = r2d[i][1];
		z = r2d[i][2];
		
		r2 = x*x+y*y+z*z;
		sum += weight[i]*exp(-a*r2)*(x*x+y*y);
	}

	fprintf(fp_debug, "\n\ntest on integral of exp(-ar*r)*(x^2+y^2) over all space\n");
	fprintf(fp_debug, "NUME = %15.12f\n", sum);
	fprintf(fp_debug, "ANAL = %15.12f\n", pow(Pi, 1.5)/pow(a, 2.5));


	fprintf(fp_debug, "\n\nIterate_p_u\n");
	fprintf(fp_debug, "               %10d %10d %10d %10d %10d %10d\n", 1, 2, 3, 4, 5, 6);
	for (i = 0; i < 20; i++)
	{
		x = -1 + i*0.1;
		fprintf(fp_debug, "x = %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", x, Iterate_p_u (x, 1), Iterate_p_u (x, 2),
			Iterate_p_u (x, 3), Iterate_p_u (x, 4), Iterate_p_u (x, 5), Iterate_p_u (x, 6) );
	}

	fprintf(fp_debug, "\n\nCalc_s_u\n");
	fprintf(fp_debug, "               %10d %10d %10d %10d %10d %10d %10d %10d %10d\n", 1, 2, 3, 4, 5, 6, 7, 8, 9);
	for (i = 0; i < 20; i++)
	{
		x = -1 + i*0.1;
		fprintf(fp_debug, "x = %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", x, Calc_s_u (x, 1), Calc_s_u (x, 2),
			Calc_s_u (x, 3), Calc_s_u (x, 4), Calc_s_u (x, 5), Calc_s_u (x, 6), Calc_s_u (x, 7), Calc_s_u (x, 8), Calc_s_u (x, 9) );
	}


}



void Set_Gauss_Chebyshev_Grid (int n, double **xi, double **wi, int DEBUG, FILE *fp_debug)
{
	int i;

	(*xi) = D_Allocate_1D_Matrix (n);
	(*wi) = D_Allocate_1D_Matrix (n);

	for (i = 0; i < n; i++)
	{
		(*xi)[n-1-i] = cos( (i+1.0)/(n+1.0)*Pi);
		(*wi)[n-1-i] = sin( (i+1.0)/(n+1.0)*Pi)*sin( (i+1.0)/(n+1.0)*Pi)*Pi/(n+1);
	}

	if (DEBUG)
	{
		fprintf(fp_debug, "Setting Gauss Chebyshev Grid of the Second Kind\n");
		fprintf(fp_debug, "n = %d\n", n);
		for (i = 0; i < n; i++)
			fprintf(fp_debug, "xi = %8.5f wi = %8.5f\n", (*xi)[i], (*wi)[i]);
	}
}

void Make_Atomic_Grid (double rm, double *r_atom, int *n_pt, double ***r2d, double **weight)
/*rm is half of the Bragg radius except for hydrogen*/
{
	int n, m, p, count;
	double x, w_rad, rad;
	int n_rad;

	n_rad = 0;
	for (n = 0; n < GRID_n_GC2_pt; n++)
	{
		x     = GRID_GC2_xi[n];
		rad   = rm*(1.0+x)/(1.0-x);
		if (rad < MAX_RAD)
			n_rad++;
	}

	(*n_pt)   = GRID_n_ang_pt*n_rad;
	(*r2d)    = D_Allocate_2D_Matrix ((*n_pt), 3);
	(*weight) = D_Allocate_1D_Matrix ((*n_pt) );

	count = 0;
	for (n = 0; n < GRID_n_GC2_pt; n++)
	{
		x     = GRID_GC2_xi[n];
		rad   = rm*(1.0+x)/(1.0-x);
		w_rad = GRID_GC2_wi[n]*(rad+rm)*(rad+rm)*(rad+rm)/4.0*pow(rad/rm, 1.5);
		if (rad < MAX_RAD)
			for (m = 0; m < GRID_n_ang_pt; m++)
			{
				for (p = 0; p < 3; p++)
					(*r2d)[count][p] = rad*GRID_ang_r2d[m][p] + r_atom[p];
				(*weight)[count] = w_rad*GRID_ang_weight[m];
				count++;
			}
	}
}

void Make_Radial_Atomic_Grid (char *element, int *n_pt_rad, double **radius, double **weight_radius)
/* INPUT:  element
   OUTPUT: n_pt_rad, radius[n_pt_rad], weight_radius[n_pt_rad]

*/
{
	int n, m, p, count;
	double x, w_rad, rad;
	int n_rad;

	int ele_num;
	double rm;


/*rm is half of the Bragg radius except for hydrogen*/
	ele_num = Find_Atom_Type (element, GRID.n_element, GRID.element);
	rm = GRID.bragg_radii[ele_num];
	if (strcmp(element, "H") != 0)
		rm /= 2;

	(*n_pt_rad) = 0;
	for (n = 0; n < GRID_n_GC2_pt; n++)
	{
		x     = GRID_GC2_xi[n];
		rad   = rm*(1.0+x)/(1.0-x);
		if (rad < MAX_RAD)
			(*n_pt_rad)++;
	}

	(*radius)        = D_Allocate_1D_Matrix ((*n_pt_rad) );
	(*weight_radius) = D_Allocate_1D_Matrix ((*n_pt_rad) );

	count = 0;
	for (n = 0; n < GRID_n_GC2_pt; n++)
	{
		x     = GRID_GC2_xi[n];
		rad   = rm*(1.0+x)/(1.0-x);
		w_rad = GRID_GC2_wi[n]*(rad+rm)*(rad+rm)*(rad+rm)/4.0*pow(rad/rm, 1.5);

		if (rad < MAX_RAD)
		{
			(*radius)[count]        = rad;
			(*weight_radius)[count] = w_rad;
			count++;
		}
	}
}

