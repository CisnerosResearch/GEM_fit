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

/************LINEAR EQUATION SOLVER*************/
void fsolve(double **H, double *b, double *x, int N);

/******ELEMENTARY MATRIX OPERATIONS************/
void Matrix_Inverse (double **A_original, double **A, double **B, int N, int DEBUG, FILE *fp_debug); /*A = matrix to be inversed*/
void Switch_Row (double **A, int row1, int row2, int N);
void Vector_Add (double *a, double *b, double *c, int N);
void Print_Matrix (FILE *fp, double **A, int N, int M);
void Print_Vector (FILE *fp, double *A, int N);
void Copy_Matrix (double **A, double **B, int N, int M);
double N_dot_prod (double *a, double *b, int N);
void Mat_Vec_Product (double *y, double **A, double *x, double Constant, int N);
void Matrix_Prod (double **A, double **B, double **C, int N, int M, int P);

/******************MISC FUNCTIONS***************/
double int_pow (double x, int n);
double dot_prod (double *a, double *b);
void normalize (double *a);
double double_factorial (int n);

double sqrt_factorial (int n);
double Factorial_2 (int n, int m);
/*= n!/(m)!*/
double sqrt_Factorial_2 (int n, int m);
/*= sqrt(n!/(m)! )*/
double factorial (int n);
int Integer_Min (int a, int b);
int Integer_Max (int a, int b);
int Is_N_Perm (int *I, int *J, int N);

double Get_Double_Rand (double a, double b);
/*Gets a random number between a and b*/
