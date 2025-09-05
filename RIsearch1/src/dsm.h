/*********************************************************************

  Copyright (c) 2016 by the contributors (see AUTHORS file)

  This file is part of RIsearch2.

  RIsearch2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
#ifndef __DSM_H__
#define __DSM_H__

#if RISVERSION == 2
 /* position of '-' in alphabet, not as define if read from matrix... */
#define GAP 0
#define MATVERSION 2
#elif RISVERSION == 1
/* position of '-' in alphabet, not as define if read from matrix... */
#define GAP 5
#define MATVERSION 1
#endif

/* dsm size */
#define DSM_N 6
#define DSM_DIM 6*6*6*6
#define DSM_N_F 1
#define SCALEFACTOR 10000.0

#define TEMPERATURE1 310.15
#define TEMPERATURE2 273.15

typedef int dsm_t[DSM_N][DSM_N][DSM_N][DSM_N];
typedef float dsmf_t[DSM_N][DSM_N][DSM_N][DSM_N];

typedef enum qt_type_enum {
	RRNA,
	DDNA,
	RDNA,
	DRNA,
	NNDNA,
}Qt_type;

extern char Qt_path[5][10]; 

extern int dsm_offset;
extern dsm_t dsm_T;
extern dsm_t dsm_extend;
extern void getMat(Qt_type qt_type, const char *matname, const char *matname2, const char *matpath, float *T, int *bA_nu);
extern int extPen;
#endif
