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

struct LEB_Data 
{
	double A1;
	double A2;
	double A3;

	int N1;
	double *Bk;	/*[N1]*/
	double *mk;	/*[N1]*/

	int N2;
	double *C1;
	double *p;

	int N3;
	double *D1;
	double *r;
	double *u;

	int N_check;
};

void Debug_Lebedev_Grid (int DEBUG, FILE *fp_debug);

void Set_Lebedev_Grid_Parm (int index1, int index2, struct LEB_Data *leb, int DEBUG, FILE *fp_debug);
/*Defines Lebedev Angular Grid Parameters in terms of constants defined in Lebedev's papers*/

void Calc_Lebedev_Grid (int *n_pt, double ***r2d, double **weight, struct LEB_Data *leb);
/*Calculates Lebedev Angular Grid*/


