
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of WIGXJPF.
 *
 *  WIGXJPF is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  WIGXJPF is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with WIGXJPF.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#ifndef __FPSIMPLEXJ_H__
#define __FPSIMPLEXJ_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

void fpsimplexj_setup_factorial_precalc();

void fpsimple_3j(double *result,
		 int two_j1, int two_j2, int two_j3,
		 int two_m1, int two_m2, int two_m3);
void fpsimple_6j(double *result,
		 int two_j1, int two_j2, int two_j3,
		 int two_j4, int two_j5, int two_j6);
void fpsimple_9j(double *result,
		 int two_j1, int two_j2, int two_j3,
		 int two_j4, int two_j5, int two_j6,
		 int two_j7, int two_j8, int two_j9);

void fpsimple_log_3j(double *result,
		     int two_j1, int two_j2, int two_j3,
		     int two_m1, int two_m2, int two_m3);
void fpsimple_log_6j(double *result,
		     int two_j1, int two_j2, int two_j3,
		     int two_j4, int two_j5, int two_j6);
void fpsimple_log_9j(double *result,
		     int two_j1, int two_j2, int two_j3,
		     int two_j4, int two_j5, int two_j6,
		     int two_j7, int two_j8, int two_j9);

#ifdef __cplusplus
}
#endif

#endif/*__FPSIMPLEXJ_H__*/
