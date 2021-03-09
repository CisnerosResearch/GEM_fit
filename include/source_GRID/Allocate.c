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

#include "Allocate.h"



/********************DOUBLE *****************************/
double *D_Allocate_1D_Matrix (int N)
{
	double *A;
	A = (double *) malloc(N*sizeof(double) );
	return A;
}

double **D_Allocate_2D_Matrix (int N, int M)
{
	int i;
	double **A;

	A = (double **) malloc(N*sizeof(double*) );
	for (i = 0; i < N; i++)
		A[i] = (double *) malloc(M*sizeof(double) );
	return A;
}

double ***D_Allocate_3D_Matrix (int N, int M, int L)
{
	int i, j;
	double ***A;

	A = (double ***) malloc(N*sizeof(double**) );
	for (i = 0; i < N; i++)
	{
		A[i] = (double **) malloc(M*sizeof(double*) );
		for (j = 0; j < M; j++)
			A[i][j] = (double *) malloc(L*sizeof(double) );
	}
	return A;
}

double ****D_Allocate_4D_Matrix (int N, int M, int L, int K)
{
	int i;
	double ****A;

	A = (double ****) malloc(N*sizeof(double***) );
	for (i = 0; i < N; i++)
		A[i] = D_Allocate_3D_Matrix (M, L, K);
	return A;
}

double *****D_Allocate_5D_Matrix (int N, int M, int L, int K, int P)
{
	int i;
	double *****A;

	A = (double *****) malloc(N*sizeof(double****) );
	for (i = 0; i < N; i++)
		A[i] = D_Allocate_4D_Matrix (M, L, K, P);

	return A;
}

void D_free_2D (double **A, int N)
{
	int i;
	for (i = 0; i < N; i++)
		free(A[i]);
	free(A);
}

void D_free_3D (double ***A, int N, int M)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < M; j++)
			free(A[i][j]);
		free(A[i]);
	}
	free(A);
}

void D_free_4D (double ****A, int N, int M, int L)
{
	int i;

	for (i = 0; i < N; i++)
		D_free_3D (A[i], M, L);

	free(A);
}

void D_free_5D (double *****A, int N, int M, int L, int K)
{
	int i;

	for (i = 0; i < N; i++)
		D_free_4D (A[i], M, L, K);

	free(A);
}


/********************INT *****************************/
int *I_Allocate_1D_Matrix (int N)
{
	int *A;
	A = (int *) malloc(N*sizeof(int) );
	return A;
}

int **I_Allocate_2D_Matrix (int N, int M)
{
	int i;
	int **A;

	A = (int **) malloc(N*sizeof(int*) );
	for (i = 0; i < N; i++)
		A[i] = (int *) malloc(M*sizeof(int) );
	return A;
}

int ***I_Allocate_3D_Matrix (int N, int M, int L)
{
	int i, j;
	int ***A;

	A = (int ***) malloc(N*sizeof(int**) );
	for (i = 0; i < N; i++)
	{
		A[i] = (int **) malloc(M*sizeof(int*) );
		for (j = 0; j < M; j++)
			A[i][j] = (int *) malloc(L*sizeof(int) );
	}
	return A;
}

void I_free_2D (int **A, int N)
{
	int i;
	for (i = 0; i < N; i++)
		free(A[i]);
	free(A);
}

void I_free_3D (int ***A, int N, int M)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < M; j++)
			free(A[i][j]);
		free(A[i]);
	}
	free(A);
}





/********************CHAR *****************************/
char *C_Allocate_1D_Matrix (int N)
{
	char *A;
	A = (char *) malloc(N*sizeof(char) );
	return A;
}

char **C_Allocate_2D_Matrix (int N, int M)
{
	int i;
	char **A;

	A = (char **) malloc(N*sizeof(char*) );
	for (i = 0; i < N; i++)
		A[i] = (char *) malloc(M*sizeof(char) );
	return A;
}

char ***C_Allocate_3D_Matrix (int N, int M, int L)
{
	int i, j;
	char ***A;

	A = (char ***) malloc(N*sizeof(char**) );
	for (i = 0; i < N; i++)
	{
		A[i] = (char **) malloc(M*sizeof(char*) );
		for (j = 0; j < M; j++)
			A[i][j] = (char *) malloc(L*sizeof(char) );
	}
	return A;
}


void C_free_2D (char **A, int N)
{
	int i;
	for (i = 0; i < N; i++)
		free(A[i]);
	free(A);
}

void C_free_3D (char ***A, int N, int M)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < M; j++)
			free(A[i][j]);
		free(A[i]);
	}
	free(A);
}

