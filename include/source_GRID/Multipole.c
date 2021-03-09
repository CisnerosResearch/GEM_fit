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


double **A_LM;

void Translate_Multipole2_To_Multipole1 (double **C_lm_R_r, double **C_lm_R_i,
	int lmax1, double **Q1_lm_r, double **Q1_lm_i, double *R1,
	int lmax2, double **Q2_lm_r, double **Q2_lm_i, double *R2);
/*C_lm_R is a temporary array with the same size lmax2.
  Translates multipole set 2 at R2 to multipole set 1 at R1.
*/
double **Allocate_Spherical_Array (int l);
/*Allocates array spherical array of the form:
A[0][0]
A[1][0], A[1][1], A[1][-1]
A[2][0], A[2][1], A[2][2], A[2][-1], A[2][-2]
..

A[l][0], A[l][1],... A[l][l], A[l][-1], .. A[l][-l]
*/
int Array_Index_2_Sphere_Index (int p, int l);
int Sphere_Index_2_Array_Index (int m, int l);


void Solid_Harm_Recur (double **C_lm_r, double **C_lm_i, int lmax, double x, double y, double z);
/*Given a point x,y,z, this function calculates solid harmonic function with both real C_lm_r[l][m] and 
  imaginary parts C_lm_i[l][m] up to l <= lmax*/
void Convert_Spherical_Multipole_2_TR_Cart_Multipole (int lmax, double ***CART_multi, double **Q_lm_r, 
	double **Q_lm_i, int DEBUG, FILE *fp_debug);
/*Converts Spherical Multipoles into a Traceless Cartesian multipole.  Does not output CART_multi in Buckingham convention, 
  to convert CART_multi to the Buckingham convention, multiply by (2l-1)!!*/
void Convert_Pure_Cart_Multipole_2_Spherical_Multipole (int LMAX, double ***cart_moment, double **Q_lm_r, 
	double **Q_lm_i, int DEBUG, FILE *fp_debug);
/*Converts a pure cartesian multipole moment (cart_moment[lx][ly][lz] lx+ly+lz <= LMAX) into complex spherical
  multipole moments with real Q_lm_r[l][m] and imaginary Q_lm_i[l][m] parts for l <= LMAX.
*/


double Calc_I_n_lx_ly (int lx, int ly, int n);
/*a constant defined by Sum over q ( lx!/(n-q)!(lx-n+q)!*ly!/q!/(ly-q)!*(-1)^(n-q) from q = max(0,n-lx) to q = min(n,ly)*/

void Nxyz_Index_2_String (char *str, int n1, int n2, int n3);

void Initialize_Atomic_Multipole (int lmax, int DEBUG, FILE *fp_debug)
/*Calculates the constant defined by A[l][m] = sqrt( (l+m)!*(l-m)! )*/
{
	int l, m;

	A_LM = Allocate_Spherical_Array (lmax);
	for (l = 0; l <= lmax; l++)
		for (m = 0; m <= l; m++)
		{
			A_LM[l][m]   = sqrt_factorial (l+m)*sqrt_factorial (l-m);
			A_LM[l][m+l] = A_LM[l][m];
		}
}

void Calculate_Multipoles (struct ELE_Data *ele, struct MOL_Data *mol, struct GRID_HI_Data *grid,
	int DEBUG, FILE *fp_debug)
/*This function converts the ab initio (QM) 'pure' Cartesian molecular multipoles obtained from the Gaussian (.log) file
  into traceless Cartesian molecular multipoles.  It then calculates the atomic multipoles from the Hirshfeld and 
  Hirshfeld-Iterated atomic charge densities.  It also calculates the molecular multipoles from the Hirshfeld and 
  Hirshfeld-Iterated atomic multipoles and compares to the QM molecular multipoles
*/
{
	int i, l, p, n, m, l1, l2;
	double r_origin[3];


/*First calculate ab initio (QM) molecular multipoles.  'Convert_Pure_Cart_Multipole_2_Spherical_Multipole'
  converts the 'pure' Cartesian molecular multipole obtained from the Gaussian (.log) output file into
  a complex spherical tensor molecular multipole.  'Convert_Spherical_Multipole_2_TR_Cart_Multipole' converts
  the complex spherical tensor molecular multipole into a traceless Cartesian molecular multipole
*/
	Convert_Pure_Cart_Multipole_2_Spherical_Multipole (4, (*mol).QM_cart_multi, (*mol).Q_QM_mol_lm_r, 
		(*mol).Q_QM_mol_lm_i, DEBUG, fp_debug);
	Convert_Spherical_Multipole_2_TR_Cart_Multipole (4, (*mol).QM_cart_multi_TR, (*mol).Q_QM_mol_lm_r, 
		(*mol).Q_QM_mol_lm_i, DEBUG, fp_debug);


	for (i = 0; i < (*mol).n_atom; i++)
	{
		for (l = 0; l <= (*mol).lmax_atom; l++)
			for (p = 0; p < 2*l+1; p++)
			{
				(*mol).Qa_H_lm_r[i][l][p] = 0;
				(*mol).Qa_H_lm_i[i][l][p] = 0;
				(*mol).Qa_HI_lm_r[i][l][p] = 0;
				(*mol).Qa_HI_lm_i[i][l][p] = 0;
			}
		(*mol).Qa_H_lm_r[i][0][0]  = (*ele).Z[ (*mol).element_num[i] ];
		(*mol).Qa_HI_lm_r[i][0][0] = (*ele).Z[ (*mol).element_num[i] ];

		for (n = 0; n < (*grid).n_grid_pts; n++)
		{
			Solid_Harm_Recur ((*mol).C_lm_r, (*mol).C_lm_i, (*mol).lmax_atom, 
				(*grid).r2d[n][0] - (*mol).r2d[i][0], 
				(*grid).r2d[n][1] - (*mol).r2d[i][1], 
				(*grid).r2d[n][2] - (*mol).r2d[i][2]);
			for (l = 0; l <= (*mol).lmax_atom; l++)
				for (p = 0; p < 2*l+1; p++)
				{
					if (fabs((*grid).p0_H_tot[n]) > 1.0E-12)
					{
						(*mol).Qa_H_lm_r[i][l][p]  -= (*mol).C_lm_r[l][p]*(*grid).p_MM_atom_H[i][n]/(*grid).p0_H_tot[n]*(*grid).p_QM[n]*(*grid).weight[n];
						(*mol).Qa_H_lm_i[i][l][p]  -= (*mol).C_lm_i[l][p]*(*grid).p_MM_atom_H[i][n]/(*grid).p0_H_tot[n]*(*grid).p_QM[n]*(*grid).weight[n];
					}
					if (fabs((*grid).p0_HI_tot[n]) > 1.0E-12)
					{
						(*mol).Qa_HI_lm_r[i][l][p] -= (*mol).C_lm_r[l][p]*(*grid).p_MM_atom_HI[i][n]/(*grid).p0_HI_tot[n]*(*grid).p_QM[n]*(*grid).weight[n];
						(*mol).Qa_HI_lm_i[i][l][p] -= (*mol).C_lm_i[l][p]*(*grid).p_MM_atom_HI[i][n]/(*grid).p0_HI_tot[n]*(*grid).p_QM[n]*(*grid).weight[n];
					}
				}
		}
	}

	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nCalculate_Atomic_Multipoles\n");
		for (i = 0; i < (*mol).n_atom; i++)
		{
			fprintf(fp_debug, "Atom %d\n", i);
			for (l = 0; l <= 2; l++)
				for (p = 0; p < 2*l+1; p++)
				{
					m = Array_Index_2_Sphere_Index (p, l);
					fprintf(fp_debug, "\tQ_H[%d][%d]  = %10.6f + %10.6fi\n", l, m, (*mol).Qa_H_lm_r[i][l][p], (*mol).Qa_H_lm_i[i][l][p]);
					fprintf(fp_debug, "\tQ_HI[%d][%d] = %10.6f + %10.6fi\n", l, m, (*mol).Qa_HI_lm_r[i][l][p], (*mol).Qa_HI_lm_i[i][l][p]);
				}
		}
	}

/*Calculate Molecular Multipole*/
	r_origin[0] = 0;  r_origin[1] = 0;  r_origin[2] = 0;

	for (l1 = 0; l1 <= (*mol).lmax_atom; l1++)
	{
		for (l2 = 0; l2 <= 4; l2++)
			for (p = 0; p < 2*l2+1; p++)
			{
				(*mol).Q_H_mol_lm_r[l1][l2][p] = 0;
				(*mol).Q_H_mol_lm_i[l1][l2][p] = 0;
			}
		for (i = 0; i < (*mol).n_atom; i++)
		{
			Translate_Multipole2_To_Multipole1 ((*mol).C_lm_r, (*mol).C_lm_i,
				4,           (*mol).Q_temp_lm_r,   (*mol).Q_temp_lm_i, &(r_origin[0]),
				l1,          (*mol).Qa_H_lm_r[i],  (*mol).Qa_H_lm_i[i],   (*mol).r2d[i]);

			for (l2 = 0; l2 <= 4; l2++)
				for (p = 0; p < 2*l2+1; p++)
				{
					(*mol).Q_H_mol_lm_r[l1][l2][p] += (*mol).Q_temp_lm_r[l2][p];
					(*mol).Q_H_mol_lm_i[l1][l2][p] += (*mol).Q_temp_lm_i[l2][p];
				}
		}
	}

	for (l1 = 0; l1 <= (*mol).lmax_atom; l1++)
	{
		for (l2 = 0; l2 <= 4; l2++)
			for (p = 0; p < 2*l2+1; p++)
			{
				(*mol).Q_HI_mol_lm_r[l1][l2][p] = 0;
				(*mol).Q_HI_mol_lm_i[l1][l2][p] = 0;
			}
		for (i = 0; i < (*mol).n_atom; i++)
		{
			Translate_Multipole2_To_Multipole1 ((*mol).C_lm_r, (*mol).C_lm_i,
				4,           (*mol).Q_temp_lm_r, (*mol).Q_temp_lm_i, &(r_origin[0]),
				l1,          (*mol).Qa_HI_lm_r[i],  (*mol).Qa_HI_lm_i[i],   (*mol).r2d[i]);

			for (l2 = 0; l2 <= 4; l2++)
				for (p = 0; p < 2*l2+1; p++)
				{
					(*mol).Q_HI_mol_lm_r[l1][l2][p] += (*mol).Q_temp_lm_r[l2][p];
					(*mol).Q_HI_mol_lm_i[l1][l2][p] += (*mol).Q_temp_lm_i[l2][p];
				}
		}
	}
/*
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nMolecular Multipole\n");
		for (l = 0; l <= 2; l++)
			for (p = 0; p < 2*l+1; p++)
			{
				m = Array_Index_2_Sphere_Index (p, l);
				fprintf(fp_debug, "\tQ[%d][%d] = %10.6f + %10.6fi\n", l, m, (*mol).Q_mol_lm_r[l][p], (*mol).Q_mol_lm_i[l][p]);
			}
	}
*/
	for (l1 = 0; l1 <= (*mol).lmax_atom; l1++)
	{
		Convert_Spherical_Multipole_2_TR_Cart_Multipole (4, (*mol).MM_H_cart_multi[l1], 
			(*mol).Q_H_mol_lm_r[l1], (*mol).Q_H_mol_lm_i[l1], DEBUG, fp_debug);
		Convert_Spherical_Multipole_2_TR_Cart_Multipole (4, (*mol).MM_HI_cart_multi[l1], 
			(*mol).Q_HI_mol_lm_r[l1], (*mol).Q_HI_mol_lm_i[l1], DEBUG, fp_debug);
	}



}


void Translate_Multipole2_To_Multipole1 (double **C_lm_R_r, double **C_lm_R_i,
	int lmax1, double **Q1_lm_r, double **Q1_lm_i, double *R1,
	int lmax2, double **Q2_lm_r, double **Q2_lm_i, double *R2)
/*C_lm_R is a temporary array with the same size lmax2.
  Translates multipole set 2 at R2 to multipole set 1 at R1.
*/
{
	int l,  m,  p;
	int l1, m1, p1;
	int l2, m2, p2;
	double C;

	Solid_Harm_Recur (C_lm_R_r, C_lm_R_i, lmax1, R2[0] - R1[0], R2[1] - R1[1], R2[2] - R1[2]);

	for (l = 0; l <= lmax1; l++)
		for (p = 0; p <= 2*l; p++)
		{
			m = Array_Index_2_Sphere_Index (p, l);

			Q1_lm_r[l][p] = 0;
			Q1_lm_i[l][p] = 0;
			for (l1 = 0; l1 <= lmax2; l1++)
				for (p1 = 0; p1 <= 2*l1; p1++)
				{
					m1 = Array_Index_2_Sphere_Index (p1, l1);
					l2 = l - l1;
					if (l2 >= 0)
					{
						m2 = m - m1;
						p2 = Sphere_Index_2_Array_Index (m2, l2);
						if (p2 != -1)
						{
							C = (A_LM[l][p]/A_LM[l1][p1])/A_LM[l2][p2];
							Q1_lm_r[l][p] += C*(Q2_lm_r[l1][p1]*C_lm_R_r[l2][p2] - Q2_lm_i[l1][p1]*C_lm_R_i[l2][p2]);
							Q1_lm_i[l][p] += C*(Q2_lm_r[l1][p1]*C_lm_R_i[l2][p2] + Q2_lm_i[l1][p1]*C_lm_R_r[l2][p2]);
						}
					}		
				}
		}
}

void Solid_Harm_Recur (double **C_lm_r, double **C_lm_i, int lmax, double x, double y, double z)
{
	int l, m;
	double fact, fact1, fact2;
	double r2;

	r2 = x*x + y*y + z*z;


/*Diagonal Recursion*/
	C_lm_r[0][0] = 1;
	C_lm_i[0][0] = 0;

	for (l = 0; l < lmax; l++)
	{
		C_lm_r[l+1][l+1] = -sqrt( (2.0*l+1)/(2.0*l+2) );
		C_lm_i[l+1][l+1] = C_lm_r[l+1][l+1];
		C_lm_r[l+1][l+1] *= (x*C_lm_r[l][l] - y*C_lm_i[l][l]);
		C_lm_i[l+1][l+1] *= (x*C_lm_i[l][l] + y*C_lm_r[l][l]);
	}

/*Downward Recursion*/
	for (l = 0; l < lmax; l++)
	{
		fact = sqrt(2.0*l+1);
		C_lm_r[l+1][l] = fact*z*C_lm_r[l][l];
		C_lm_i[l+1][l] = fact*z*C_lm_i[l][l];
	}

	for (l = 1; l < lmax; l++)
	{
		for (m = 0; m <= l-1; m++)
		{
			fact1 = sqrt( 1.0/(l+m+1.0)/(l-m+1.0) );
			fact2 = sqrt( 1.0*(l+m)*(l-m) );

			C_lm_r[l+1][m] = fact1*( (2.0*l+1)*z*C_lm_r[l][m] - fact2*r2*C_lm_r[l-1][m]);
			C_lm_i[l+1][m] = fact1*( (2.0*l+1)*z*C_lm_i[l][m] - fact2*r2*C_lm_i[l-1][m]);
		}
	}

/*Negative m*/
	for (l = 0; l <= lmax; l++)
		for (m = 1; m <= l; m++)
		{
			C_lm_r[l][l+m] =  C_lm_r[l][m];
			C_lm_i[l][l+m] = -C_lm_i[l][m];
			if (m%2!=0)
			{
				C_lm_r[l][l+m] *= -1;
				C_lm_i[l][l+m] *= -1;
			}
		}

}


double **Allocate_Spherical_Array (int l)
/*Allocates array spherical array of the form:
A[0][0]
A[1][0], A[1][1], A[1][-1]
A[2][0], A[2][1], A[2][2], A[2][-1], A[2][-2]
..

A[l][0], A[l][1],... A[l][l], A[l][-1], .. A[l][-l]
*/
{
	int i;
	double **A;
	A = (double **) malloc( (l+1)*sizeof(double*) );
	for (i = 0; i <= l; i++)
		A[i] = (double *) malloc((2*i+1)*sizeof(double) );
	return A;
}

int Array_Index_2_Sphere_Index (int p, int l)
{
	if ( (p > 2*l) || (p < -2*l) )
	{
		printf("ERROR: Array_Index_2_Sphere_Index p = %d, l = %d\n", p, l);
		exit(0);
	}
	if (p > l)
		return -(p - l);
	else
		return p;
}
int Sphere_Index_2_Array_Index (int m, int l)
{
	if ( (m < -l) || (m > l) )
		return -1;
	if (m < 0)
		return l-m;
	else
		return m;
}

void Convert_Pure_Cart_Multipole_2_Spherical_Multipole (int LMAX, double ***cart_moment, double **Q_lm_r, 
	double **Q_lm_i, int DEBUG, FILE *fp_debug)
/*Converts a pure cartesian multipole moment (cart_moment[lx][ly][lz] lx+ly+lz <= LMAX) into complex spherical
  multipole moments with real Q_lm_r[l][m] and imaginary Q_lm_i[l][m] parts for l <= LMAX.
*/
{

	int l, m, p, q, k;
	double term;

	for (l = 0; l <= LMAX; l++)
		for (m = 0; m <= l; m++)
		{
			Q_lm_r[l][m] = 0;
			Q_lm_i[l][m] = 0;

			term = 0;
			for (k = 0; 2*k <= (l - m); k++)
				for (p = 0; p <= m+k; p++)
					for (q = 0; q <= k; q++)
					{
						term = A_LM[l][m]/int_pow (2.0, m+2*k);
						term /= (factorial (p)*factorial (m+k-p)*factorial (q)*factorial (k-q)*factorial (l-m-2*k) );
						term *= cart_moment[p+q][m+2*k-p-q][l-m-2*k];

						if ( (m+q)%2 != 0)
							term = -term;

						if      ( (m-p-q+2*k)%4 == 0)
							Q_lm_r[l][m] += term;
						else if ( (m-p-q+2*k)%4 == 1)
							Q_lm_i[l][m] += term;
						else if ( (m-p-q+2*k)%4 == 2)
							Q_lm_r[l][m] += -term;
						else if ( (m-p-q+2*k)%4 == 3)
							Q_lm_i[l][m] += -term;
					}

			p = Sphere_Index_2_Array_Index (-m, l);
			if (m%2 == 0)
			{
				Q_lm_r[l][p] =  Q_lm_r[l][m];
				Q_lm_i[l][p] = -Q_lm_i[l][m];
			}
			else
			{
				Q_lm_r[l][p] = -Q_lm_r[l][m];
				Q_lm_i[l][p] =  Q_lm_i[l][m];
			}
	
		}
}

void Convert_Spherical_Multipole_2_TR_Cart_Multipole (int lmax, double ***CART_multi, double **Q_lm_r, 
	double **Q_lm_i, int DEBUG, FILE *fp_debug)
/*Converts Spherical Multipoles into a Traceless Cartesian multipole.  Does not output CART_multi in Buckingham convention, 
  to convert CART_multi to the Buckingham convention, multiply by (2l-1)!!*/
{
	int l, l1, l2, l3, m, m_abs, I_n_lx_ly, n;
	double real, imag, fact;
	char str[100];

	for (l = 0; l <= lmax; l++)
		for (l1 = 0; l1 <= l; l1++)
			for (l2 = 0; l2 <= l - l1; l2++)
			{
				l3 = l - l1 - l2;

				real = 0;
				imag = 0;
				fact = 1.0/factorial (l)/double_factorial (2*l-1)/int_pow (2.0, l1 + l2);
				for (n = 0; n <= l1+l2; n++)
				{
					m     = 2*n - l1 - l2;
					m_abs = abs(m);
					I_n_lx_ly = Calc_I_n_lx_ly (l1, l2, n);
					if (m < 0)
					{
						if (m%2 == 0)
						{
							real += A_LM[l][m_abs]*Q_lm_r[l][m_abs]*I_n_lx_ly;
							imag -= A_LM[l][m_abs]*Q_lm_i[l][m_abs]*I_n_lx_ly;
						}
						else
						{
							real -= A_LM[l][m_abs]*Q_lm_r[l][m_abs]*I_n_lx_ly;
							imag += A_LM[l][m_abs]*Q_lm_i[l][m_abs]*I_n_lx_ly;
						}
					}
					else
					{
						real += A_LM[l][m_abs]*Q_lm_r[l][m_abs]*I_n_lx_ly;
						imag += A_LM[l][m_abs]*Q_lm_i[l][m_abs]*I_n_lx_ly;
					}
				}
				if (l2%4 == 0)
					CART_multi[l1][l2][l3] =  fact*real;
				else if (l2%4 == 1)
					CART_multi[l1][l2][l3] = -fact*imag;
				else if (l2%4 == 2)
					CART_multi[l1][l2][l3] = -fact*real;
				else if (l2%4 == 3)
					CART_multi[l1][l2][l3] =  fact*imag;

				if (l2%2 == 0)
				{
					if (fabs(imag) > 1.0E-8)
						printf("WARNING in Convert_Spherical_Multipole_2_TR_Cart_Multipole.. l2 = %d, imag = %12.10f != 0.0\n", l2, imag);
				}
				else
				{
					if (fabs(real) > 1.0E-8)
						printf("WARNING in Convert_Spherical_Multipole_2_TR_Cart_Multipole.. l2 = %d, real = %12.10f != 0.0\n", l2, real);
				}			
			}


	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nConvert_Spherical_Multipole_2_TR_Cart_Multipole\n");
		fprintf(fp_debug, "\tlmax = %d\n", lmax);
		for (l = 0; l <= lmax; l++)
			for (l1 = 0; l1 <= l; l1++)
				for (l2 = 0; l2 <= l - l1; l2++)
				{
					l3 = l - l1 - l2;
					Nxyz_Index_2_String (str, l1, l2, l3);
					fprintf(fp_debug, "\t%10s = %12.8f\n", str, CART_multi[l1][l2][l3]);
				}
	}


}

double Calc_I_n_lx_ly (int lx, int ly, int n)
/*a constant defined by Sum over q ( lx!/(n-q)!(lx-n+q)!*ly!/q!/(ly-q)!*(-1)^(n-q) from q = max(0,n-lx) to q = min(n,ly)*/
{
	int start, end, q;
	double sum, term;
	
	start = 0;
	if (start < n - lx)
		start = n - lx;

	end = n;
	if (end > ly)
		end = ly;

	sum = 0;
	for (q = start; q <= end; q++)
	{
		term = Factorial_2 (lx, n - q)/factorial (lx - n + q)*Factorial_2 (ly, q)/factorial (ly - q);
		if ( (n-q)%2 != 0)
			term = -term;
		sum += term;
	}
	return sum;
}


