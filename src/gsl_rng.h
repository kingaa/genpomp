/* rng/gsl_rng.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_RNG_H__
#define __GSL_RNG_H__

#include <stdlib.h>
#include <stdio.h>
#include "RngStream.h"
#include "constants.h"
//#include <gsl/gsl_types.h>
//#include <gsl/gsl_errno.h>

//typedef struct RngStream_InfoState gsl_rng;
typedef RngStream gsl_rng;

unsigned long int gsl_rng_get (gsl_rng * r);
double gsl_rng_uniform (gsl_rng * r);
double gsl_rng_uniform_pos (gsl_rng * r);
unsigned long int gsl_rng_uniform_int (gsl_rng * r, unsigned long int n);

#endif /* __GSL_RNG_H__ */



