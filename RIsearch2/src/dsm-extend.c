/*********************************************************************
  RIsearch2   --   RNA-RNA interaction search

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
#include "dsm.h"

/*#nt added, to be multiplied by extpen */
const dsm_t dsm_extend = {\
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

