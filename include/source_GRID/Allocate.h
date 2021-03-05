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

double *D_Allocate_1D_Matrix (int N);
double **D_Allocate_2D_Matrix (int N, int M);
double ***D_Allocate_3D_Matrix (int N, int M, int L);
double ****D_Allocate_4D_Matrix (int N, int M, int L, int K);
double *****D_Allocate_5D_Matrix (int N, int M, int L, int K, int P);

void D_free_2D (double **A, int N);
void D_free_3D (double ***A, int N, int M);
void D_free_4D (double ****A, int N, int M, int L);
void D_free_5D (double *****A, int N, int M, int L, int K);

int *I_Allocate_1D_Matrix (int N);
int **I_Allocate_2D_Matrix (int N, int M);
int ***I_Allocate_3D_Matrix (int N, int M, int L);
void I_free_2D (int **A, int N);
void I_free_3D (int ***A, int N, int M);


char *C_Allocate_1D_Matrix (int N);
char **C_Allocate_2D_Matrix (int N, int M);
char ***C_Allocate_3D_Matrix (int N, int M, int L);
void C_free_2D (char **A, int N);
void C_free_3D (char ***A, int N, int M);
