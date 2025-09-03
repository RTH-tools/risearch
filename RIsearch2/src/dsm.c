/***********************************************************
  Copyright 2012 Anne Wenzel <wenzel@rth.dk>

  This file is part of RIsearch.

  RIsearch is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.

***********************************************************/

#include <stdio.h>		/* printf, NULL */
#include <stdlib.h>		/* strtof */
#include <string.h>		/* strlen */
#include <math.h>		/* round */
#include <assert.h>
#include <unistd.h>

#include "main.h"
#include "dsm.h"

/*equals dsm for extpen=0 */
/* based on Turner 1999 parameters. */
/* Compared to RIsearch1.1, here ending with a bulge has the same cost as continuing a bulge.*/
/* ../tables/risearch2/dsm_t04.1.1.tsv */
/* ../tables/risearch2/dsm_t04.1.2.tsv */


/* based on Turner 2004 parameters Compared to RIsearch1.1, here ending with a
 * bulge has the same cost as continuing a bulge 
 * ../tables/risearch2/dsm_t99.1.1.tsv 
 * ../tables/risearch2/dsm_t99.1.2.tsv */

/* DNA_DNA matrix based on Santa Lucia and Hicks 2004 parameters. G-T wobble is
 * not allowed Matrix converted to RIsearch1 format from RIsearch2.1. Compared
 * to RIsearch2.1, here ending with a bulge has the same cost as continuing a
 * bulge.*/

/* ../tables/risearch2/dsm_slh04_noGU.1.1.tsv */
/* ../tables/risearch2/dsm_slh04_noGU.1.2.tsv */

/* RNA_DNA matrix based on Sugimoto 95 stacking parameter. Bulges and internal
 * loops are estimated as the mean of t04 and slh04 initiation or extension
 * parameters. G-U wobble is allowed. Matrix converted to RIsearch1 format from
 * RIsearch2.1. Compared to RIsearch2.1, here ending with a bulge has the same
 * cost as continuing a bulge. */

/* ../tables/risearch2/dsm_su95_rev_wGU.1.1.tsv 
 * ../tables/risearch2/dsm_su95_rev_wGU.1.2.tsv */

/* RNA_DNA matrix based on Sugimoto 95 stacking parameter. Bulges and internal
 * loops are estimated as the mean of t04 and slh04 initiation or extension
 * parameters. G-U wobble is not allowed in the seed Matrix converted to
 * RIsearch1 format from RIsearch2.1. Compared to RIsearch2.1, here ending with
 * a bulge has the same cost as continuing a bulge.*/
/* ../tables/risearch2/dsm_su95_rev_noGU.1.2.tsv */
/* ../tables/risearch2/dsm_su95_rev_noGU.1.1.tsv */

int dsm_offset;
dsm_t dsm_T;
/*#nt added, to be multiplied by extpen */
int dsm_extend[DSM_N][DSM_N][DSM_N][DSM_N] = {\
{ \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} }, \
} , {
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} }, \
} , {
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} }, \
} , {
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} }, \
} , {
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} }, \
{ {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} }, \
} , {
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {0,0,0,2,0,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {0,0,2,0,0,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {0,2,0,2,0,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,0,2,0,0,1} }, \
{ {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {2,2,2,2,2,1} , {0,0,0,0,0,1} }, \
{ {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} , {1,1,1,1,1,0} } \
} };
/*0s in last block to have init score clean!*/


/* transform matrix into Risearch1 layout*/
static void to_version(int version, int *i0, int *i1, int *i2, int *i3)
{
	int j0;
	int j1;
	int j2;
	int j3;
	if (version == 2) {
		int to2[DSM_N] = { 1, 3, 2, 4, 5, 0 };
		*i0 = to2[*i0];
		*i1 = to2[*i1];
		*i2 = to2[*i2];
		*i3 = to2[*i3];
		j0 = *i0;
		j1 = *i1;
		j2 = *i2;
		j3 = *i3;
		*i0 = j3;
		*i1 = j2;
		*i2 = j1;
		*i3 = j0;
		return;
	} else if (version == 1) {
		j0 = *i0;
		j1 = *i1;
		j2 = *i2;
		j3 = *i3;
		*i0 = j3;
		*i1 = j2;
		*i2 = j1;
		*i3 = j0;
		//int to1[DSM_N] = {5,0,2,1,3,4};
		return;
	}
	assert(version > 0 && version < 3);
}

/* load mat from tsv file, first line is off-set*/
static int load_mat_tsv(dsmf_t newmat, float *offset, const char *fn, int version)
{
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	int i = 0;
	fp = fopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "failed to open '%s'\n", fn);
		return (1);
	}
	if (verbose > 0)
		fprintf(stderr, "mapping to version %d\n", version);
	while ((read = getline(&line, &len, fp)) != -1) {
		char *pEnd = 0;
		char *pSt = line;
		float f1;
		while (i < DSM_DIM + DSM_N_F) {
			f1 = strtof(pSt, &pEnd);
			if (pEnd == pSt) {
				break;
			}
			if (i == 0) {
				*offset = f1;
			} else {
				int i0, i1, i2, i3;
				i0 = (i - 1) % DSM_N;
				i1 = (i - 1 - i0) / DSM_N % DSM_N;
				i2 = (i - 1 - i0 - DSM_N * i1) / DSM_N / DSM_N % DSM_N;
				i3 = (i - 1 - i0 - DSM_N * i1 - DSM_N * DSM_N * i2) / DSM_N / DSM_N / DSM_N % DSM_N;
				to_version(version, &i0, &i1, &i2, &i3);
				newmat[i0][i1][i2][i3] = -f1;
			}
			pSt = pEnd;
			i++;
		}
	}
	if (verbose > 0)
		fprintf(stderr, "Read %d floats, expected %d floats from %s\n", i, DSM_DIM + DSM_N_F, fn);
	if (line)
		free(line);
	fclose(fp);
	return (0);
}

#define TO_T(E1,E2,T) round(SCALEFACTOR*((T[0]-T[2])/(T[1]-T[2])*(E1 - E2)+E2))

static void scale_to_T(dsmf_t dsmf_1, dsmf_t dsmf_2, float *T)
{
	int i0, i1, i2, i3;
	if (abs(T[1] - T[2]) < 0.001) {
#if VERBOSE > 1
		fprintf(stderr, "T1 == T2, so no scaling\n");
#endif
		return;
	}
	for (i0 = 0; i0 < DSM_N; i0++) {
		for (i1 = 0; i1 < DSM_N; i1++) {
			for (i2 = 0; i2 < DSM_N; i2++) {
				for (i3 = 0; i3 < DSM_N; i3++) {
					dsm_T[i0][i1][i2][i3] = TO_T(dsmf_1[i0][i1][i2][i3], dsmf_2[i0][i1][i2][i3], T);
#if VERBOSE > 1
					fprintf(stderr, "%d %d %d %d %d\n", dsm_T[i0][i1][i2][i3], i0, i1, i2, i3);
#endif
				}
			}
		}
	}
}

void load_mats(const char *matpath, const char *matname, const char *matname2, float *T, int matversion)
{
	float offset_1;
	float offset_2;
	dsmf_t dsmf_1;		/*  -deltaG table at 37 deg celcius, that is 301.15K */
	dsmf_t dsmf_2;		/*  -deltaG table at  0 deg kelvin, that is deltaH */
	char *mat1name = NULL, *mat2name = NULL;

	int N1 = strlen(matpath) + strlen(matname) + 10;
	mat1name = (char *)malloc(N1 * sizeof(char));
	sprintf(mat1name, "%s/%s.tsv", matpath, matname);

	if (matname2 == NULL) {
		mat2name = (char *)malloc(N1 * sizeof(char));
		sprintf(mat2name, "%s/%s.2.tsv", matpath, matname);
	} else {
		int N2 = strlen(matpath) + strlen(matname2) + 10;
		mat2name = (char *)malloc(N2 * sizeof(char));
		sprintf(mat2name, "%s/%s.tsv", matpath, matname2);
	}
	if (access(mat1name, R_OK) != 0) {
		fprintf(stderr, "%s not found or not readable\n", mat1name);
		exit(1);
	}
	if (access(mat2name, R_OK) != 0) {
		fprintf(stderr, "%s not found or not readable\n", mat2name);
		mat2name = mat1name;
	}

	if (strcmp(mat2name, mat1name) == 0) {
		fprintf(stderr, "Temperature ignored, using %s directly.\n", mat1name);
	} else {
		fprintf(stderr, "Temperature Scaling of T=%.3f, using T1=%.3f with %s and T2=%.3f with %s\n", T[0], T[1], mat1name, T[2], mat2name);
	}

	if (load_mat_tsv(dsmf_1, &offset_1, mat1name, matversion) != 0) {
		if (mat2name != mat1name)
			free(mat2name);
		free(mat1name);
		exit(1);
	}

	if (load_mat_tsv(dsmf_2, &offset_2, mat2name, matversion) != 0) {
		if (mat2name != mat1name)
			free(mat2name);
		free(mat1name);
		exit(1);
	}

	if (mat2name != mat1name)
		free(mat2name);
	free(mat1name);
	scale_to_T(dsmf_1, dsmf_2, T);
	dsm_offset = TO_T(offset_1, offset_2, T);
}

void check_mat_settings(const char *matname, const char *matname2, float *T)
{
	int i;
	for (i = 0; i < 3; i++) {
		if (T[i] < 0 || T[i] > 410.15) {
			fprintf(stderr, "Temperatures must be in the range 0 - 410.15K\n");
			exit(1);
		}
	}
	if (matname2 == NULL) {
		if (abs(T[1] - TEMPERATURE1) > 0.001 || abs(T[2] - TEMPERATURE2) > 0.001) {
			fprintf(stderr, "Non standard tempatures given, needs --matrix2 parameter\n");
			exit(1);
		}
	}
}

/** @brief Initializes RIsearch energy scoring matrix incl extension
 *  @param matname The name of the dsm matrix to use
 *  @param bA_nu The address to which the combined matrix will be written
 */
void getMat(const char *matname, const char *matname2, const char *matpath, float *T, int *bA_nu)
{
	const int *bA_bas, *bA_ext;
	int i;

	check_mat_settings(matname, matname2, T);
	load_mats(matpath, matname, matname2, T, MATVERSION);

	bA_bas = &dsm_T[0][0][0][0];
	bA_ext = &dsm_extend[0][0][0][0];

	/* create dsm from   dsm_base - d * dsm_extend   */
	for (i = 0; i < 1296; i++) {
		*(bA_nu + i) = *(bA_bas + i) - extPen * *(bA_ext + i);	/* bA_nu[i] =  ... also works */
	}
}
