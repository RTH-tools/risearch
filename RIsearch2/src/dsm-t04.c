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

/* based on Turner 2004 parameters */
const dsm_t dsm_t04 = {\
{ \
{ {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,30,-22,-71} , {-230,-230,-150,90,-230,-285} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-45,-22,-40} },\
{ {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,100,-22,-22,-71} , {-230,-230,220,-230,-230,-285} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-40} },\
{ {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,100,-22,30,-22,-71} , {-130,210,-110,60,-230,-285} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-45,-22,-40} },\
{ {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {30,-22,30,-22,-22,-71} , {110,-230,140,-160,-230,-285} , {-70,-22,-70,-22,-22,-71} , {-45,-22,-45,-22,-22,-40} },\
{ {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-230,-230,-230,-230,-230,-285} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-40} },\
{ {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-285,-285,-285,-285,-285,-45} , {-71,-71,-71,-71,-71,-2000} , {-123,-123,-123,-123,-123,-123} },\
} , {
{ {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-160,-160,-80,210,-160,-240} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-45,-22,-40} },\
{ {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-160,-160,330,-160,-160,-240} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-40} },\
{ {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-60,240,-40,140,-160,-240} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-45,-22,-40} },\
{ {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {210,-160,210,-90,-160,-240} , {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {-45,-22,-45,-22,-22,-40} },\
{ {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-160,-160,-160,-160,-160,-240} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-40} },\
{ {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-240,-240,-240,-240,-240,0} , {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-123,-123,-123,-123,-123,-123} },\
} , {
{ {-22,-22,-22,10,-22,-71} , {-160,-160,-80,240,-160,-240} , {-22,-22,-22,50,-22,-71} , {-230,-230,-150,130,-230,-285} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-45,-22,-40} },\
{ {-22,-22,80,-22,-22,-71} , {-160,-160,340,-160,-160,-240} , {-22,-22,120,-22,-22,-71} , {-230,-230,250,-230,-230,-285} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-40} },\
{ {-22,80,-22,10,-22,-71} , {-60,330,-40,150,-160,-240} , {-22,120,-22,50,-22,-71} , {-130,210,-110,50,-230,-285} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-45,-22,-40} },\
{ {10,-22,10,-22,-22,-71} , {220,-160,250,-90,-160,-240} , {50,-22,50,-22,-22,-71} , {140,-230,-130,-160,-230,-285} , {-70,-22,-70,-22,-22,-71} , {-45,-22,-45,-22,-22,-40} },\
{ {-22,-22,-22,-22,-22,-71} , {-160,-160,-160,-160,-160,-240} , {-22,-22,-22,-22,-22,-71} , {-230,-230,-230,-230,-230,-285} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-40} },\
{ {-71,-71,-71,-71,-71,-2000} , {-240,-240,-240,-240,-240,0} , {-71,-71,-71,-71,-71,-2000} , {-285,-285,-285,-285,-285,-45} , {-71,-71,-71,-71,-71,-2000} , {-123,-123,-123,-123,-123,-123} },\
} , {
{ {-230,-230,-150,130,-230,-285} , {-22,-22,-22,-70,-22,-71} , {-230,-230,-150,100,-230,-285} , {-22,-22,-22,0,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-45,-22,-40} },\
{ {-230,-230,240,-230,-230,-285} , {-22,-22,0,-22,-22,-71} , {-230,-230,150,-230,-230,-285} , {-22,-22,70,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-40} },\
{ {-130,210,-110,100,-230,-285} , {-22,0,-22,-70,-22,-71} , {-130,140,-110,-30,-230,-285} , {-22,70,-22,0,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-45,-22,-40} },\
{ {90,-230,130,-160,-230,-285} , {-70,-22,-70,-22,-22,-71} , {60,-230,50,-160,-230,-285} , {0,-22,0,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {-45,-22,-45,-22,-22,-40} },\
{ {-230,-230,-230,-230,-230,-285} , {-22,-22,-22,-22,-22,-71} , {-230,-230,-230,-230,-230,-285} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-40} },\
{ {-285,-285,-285,-285,-285,-45} , {-71,-71,-71,-71,-71,-2000} , {-285,-285,-285,-285,-285,-45} , {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-123,-123,-123,-123,-123,-123} },\
} , {
{ {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-70,-22,-71} , {-22,-22,-22,-45,-22,-40} },\
{ {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-71} , {-22,-22,0,-22,-22,-40} },\
{ {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-70,-22,-71} , {-22,0,-22,-45,-22,-40} },\
{ {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {-70,-22,-70,-22,-22,-71} , {-45,-22,-45,-22,-22,-40} },\
{ {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-71} , {-22,-22,-22,-22,-22,-40} },\
{ {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-71,-71,-71,-71,-71,-2000} , {-123,-123,-123,-123,-123,-123} },\
} , {
{ {-22,-22,-22,-45,-22,-123} , {-22,-22,-22,-45,-22,-123} , {-22,-22,-22,-45,-22,-123} , {-22,-22,-22,-45,-22,-123} , {-22,-22,-22,-45,-22,-123} , {0,0,0,105,0,-123} },\
{ {-22,-22,0,-22,-22,-123} , {-22,-22,0,-22,-22,-123} , {-22,-22,0,-22,-22,-123} , {-22,-22,0,-22,-22,-123} , {-22,-22,0,-22,-22,-123} , {0,0,150,0,0,-123} },\
{ {-22,0,-22,-45,-22,-123} , {-22,0,-22,-45,-22,-123} , {-22,0,-22,-45,-22,-123} , {-22,0,-22,-45,-22,-123} , {-22,0,-22,-45,-22,-123} , {0,150,0,105,0,-123} },\
{ {-45,-22,-45,-22,-22,-123} , {-45,-22,-45,-22,-22,-123} , {-45,-22,-45,-22,-22,-123} , {-45,-22,-45,-22,-22,-123} , {-45,-22,-45,-22,-22,-123} , {105,0,105,0,0,-123} },\
{ {-22,-22,-22,-22,-22,-123} , {-22,-22,-22,-22,-22,-123} , {-22,-22,-22,-22,-22,-123} , {-22,-22,-22,-22,-22,-123} , {-22,-22,-22,-22,-22,-123} , {0,0,0,0,0,-123} },\
{ {-40,-40,-40,-40,-40,-123} , {-40,-40,-40,-40,-40,-123} , {-40,-40,-40,-40,-40,-123} , {-40,-40,-40,-40,-40,-123} , {-40,-40,-40,-40,-40,-123} , {-123,-123,-123,-123,-123,-2000} },\
} };

