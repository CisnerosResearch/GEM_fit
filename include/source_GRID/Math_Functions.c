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

#define MAX_INVERSE_ITERATION 1000000
#define h 1.0E-6
#define eps_inverse 1.0E-3 

void Print_Matrix (FILE *fp, double **A, int N, int M);

/*External Functions Used in Derivative Debugger
    double f_x (double *x);
    void Calc_df_dx (double *df_dx, double *x);
*/

void Switch_Row (double **A, int row1, int row2, int N);

void fsolve(double **H, double *b, double *x, int N)
{
	int i,j,m;
	double sum;
	double r_ii;
	
	for (i = 0; i < N; i++)
		x[i] = 0.0;

	for (m = 0; m < MAX_INVERSE_ITERATION; m++)
		for (i = 0; i < N; i++)
		{
			sum = 0.0;
			
			for (j = 0; j < N; j++)
				if (i != j)
					sum += H[i][j]*x[j];
			
			r_ii = b[i] - sum;
	
			x[i] = r_ii/H[i][i];
		}
}

void Matrix_Inverse (double **A_original, double **A, double **B, int N, int DEBUG, FILE *fp_debug)
{
    int i, j, k;
    double Aii, Aji;
	double sum;

/*Initialize B to Identity*/
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            if (i == j)
                B[i][j] = 1.0;
            else
                B[i][j] = 0.0;


/*Make A upper triangular with diagonal 1's*/
    for (i = 0; i < N; i++)
    {
/*If A[i][i] is small, switch ith row with a row below it*/
        if (fabs(A[i][i]) < eps_inverse)
        {
            for (j = i+1; j < N; j++)
                if (A[j][j] > eps_inverse)
                {
                    Switch_Row (A, i, j, N);
                    Switch_Row (B, i, j, N);
                    break;
                }
        }

/*Divide row i by A[i][i]*/
        Aii      = A[i][i];
        for (j = 0; j < N; j++)
        {
            A[i][j] /= Aii;
            B[i][j] /= Aii;
        }

/*For each row below i, make A[j][i] = 0, (j > i)*/
        for (j = i+1; j < N; j++)
        {
            Aji = A[j][i]; 
            for (k = 0; k < N; k++)
            {
                A[j][k] += -Aji*A[i][k];
                B[j][k] += -Aji*B[i][k];
            }
        }

    }

/*Now make A's upper triangular go to 0*/
    for (i = N-1; i > -1; i--)
        for (j = i-1; j > -1; j--)
        {
            Aji = A[j][i];
            for (k = 0; k < N; k++)
            {
                A[j][k] -= Aji*A[i][k];
                B[j][k] -= Aji*B[i][k];
            }
        }
	if (DEBUG)
	{
		fprintf(fp_debug, "\n\nA_inv = \n");
		Print_Matrix (fp_debug, B, N, N);
		fprintf(fp_debug, "\n\nA*A_inv = \n");
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				sum = 0;
				for (k = 0; k < N; k++)
					sum += A_original[i][k]*B[k][j];
				fprintf(fp_debug, "%8.4f ", sum);
			}
			fprintf(fp_debug, "\n");
		}
	}


}

void Switch_Row (double **A, int row1, int row2, int N)
{
    int i;

    double temp;

    for (i = 0; i < N; i++)
    {
        temp        = A[row1][i];
        A[row1][i]  = A[row2][i];
        A[row2][i]  = temp;
    }
}









/*************MATRIX OPERATIONS******************/
void Vector_Add (double *a, double *b, double *c, int N)
{
	int i;

	for (i = 0; i < N; i++)
		a[i] = b[i] + c[i];
}

void Matrix_Prod (double **A, double **B, double **C, int N, int M, int P)
{
    int i, j, k;
    
    for (i = 0; i < N; i++)
        for (j = 0; j < P; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < M; k++)
                C[i][j] += A[i][k]*B[k][j];
        }
}

void Print_Matrix (FILE *fp, double **A, int N, int M)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < M; j++)
			fprintf(fp, "%12.8f ", A[i][j]);
		fprintf(fp, "\n");
	}
}

void Print_Vector (FILE *fp, double *A, int N)
{
	int i;

	for (i = 0; i < N; i++)
		fprintf(fp, "%15.12f ", A[i]);
	fprintf(fp, "\n");
}


void Copy_Matrix (double **A, double **B, int N, int M)
/*Copies A into B, Both of size NxM*/
{
	int i, j;

	for (i = 0; i < N; i++)
		for (j = 0; j < M; j++)
			B[i][j] = A[i][j];
}

double N_dot_prod (double *a, double *b, int N)
{
	int i;
	double sum;

	sum = 0;
	for (i = 0; i < N; i++)
		sum += a[i]*b[i];
	return sum;
}

void Mat_Vec_Product (double *y, double **A, double *x, double Constant, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        y[i] = 0;
        for (j = 0; j < N; j++)
            y[i] += A[i][j]*x[j];
        y[i] *= Constant;
    }
}


double int_pow (double x, int n)
{
	int i;
	double prod;
	prod = 1;
	if (n > 0)
		for (i = 0; i < n; i++)
			prod *= x;	
	else
	{
		n *= -1;
		for (i = 0; i < n; i++)
			prod /= x;	
	}
	return prod;
}

double dot_prod (double *a, double *b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void normalize (double *a)
{
	int p;
	double d;

	d = sqrt( dot_prod (a,a) );
	for (p = 0; p < 3; p++)
		a[p] /= d;
}	

double double_factorial (int n)
{
	int i;
	double fact;

	if (n < 2)
		return 1;

	fact = 1;
	for (;;)
	{
		fact *= n;
		n -= 2;
		if (n < 2)
			return fact;
	}
}

double sqrt_Factorial_2 (int n, int m)
/*= sqrt(n!/(m)! )*/
{
	int i;
	double fact;

	fact = 1.0;

	if (n > m)
	{
		for (i = n; i > m; i--)
			fact *= sqrt(i);
	}
	else if (n < m)
	{
		for (i = m; i > n; i--)
			fact /= sqrt(i);
	}

	return fact;
}


double Factorial_2 (int n, int m)
/*= n!/(m)!*/
{
	int i;
	double fact;

	fact = 1.0;

	if (n > m)
	{
		for (i = n; i > m; i--)
			fact *= i;
	}
	else if (n < m)
	{
		for (i = m; i > n; i--)
			fact /= i;
	}

	return fact;
}

double factorial (int n)
{
	int i;
	double fact;

	fact = 1;
	for (i = 1; i < n+1; i++)
		fact *= i;
	return fact;
}

double sqrt_factorial (int n)
{
	int i;
	double fact;

	fact = 1.0;
	for (i = 1; i < n+1; i++)
		fact *= sqrt(i);
	return fact;
}

int Is_N_Perm (int *I, int *J, int N)
/*each element of I must have same number of occurences in J*/
{
	int i, j;
	int n_occur_I, n_occur_J;

	for (i = 0; i < N; i++)
	{
		n_occur_I = 0;
		for (j = 0; j < N; j++)
			if (I[i] == I[j])
				n_occur_I++;
		
		n_occur_J = 0;
		for (j = 0; j < N; j++)
			if (I[i] == J[j])
				n_occur_J++;

		if (n_occur_I != n_occur_J)
			return 0;
	}
	return 1;
}


int Integer_Min (int a, int b)
{
	if (a > b)
		return b;
	else
		return a;
}

int Integer_Max (int a, int b)
{
	if (a > b)
		return a;
	else
		return b;
}
double Get_Double_Rand (double a, double b)
/*Gets a random number between a and b*/
{
    return a + (b - a)* ( (double) rand()/ (double)RAND_MAX );
}


