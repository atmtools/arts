
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

#ifndef __WIGXJPF_H__
#define __WIGXJPF_H__

#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************************
 *
 * First a simplified interface (user need not care about the
 * temporary array pointer).  Uses functions from the normal interface
 * below.
 */

void wig_table_init(int max_two_j, int wigner_type);
void wig_table_free();
void wig_temp_init(int max_two_j);
void wig_temp_free();

/* Use the instead of wig_temp_init to ensure that __thread support
 * is compiled into the library. */
void wig_thread_temp_init(int max_two_j);

double wig3jj(int two_j1, int two_j2, int two_j3,
	      int two_m1, int two_m2, int two_m3);
double wig6jj(int two_j1, int two_j2, int two_j3,
	      int two_j4, int two_j5, int two_j6);
double wig9jj(int two_j1, int two_j2, int two_j3,
	      int two_j4, int two_j5, int two_j6,
	      int two_j7, int two_j8, int two_j9);
  
/***********************************************************************
 *
 * And then the normal interface.
 */

#include <stdlib.h>

/* When calculating symbols up to a maximum j, the largest
 * factorial and iteration needed are:
 *
 * 3j:  max_fact: 3*j + 1  max_iter:  j + 1
 * 6j:  max_fact: 4*j + 1  max_iter:  j + 1
 * 9j:  max_fact: 5*j + 1  max_iter:  j + 1
 *
 * For odd j, the values can be rounded down.  Note, the above table
 * is in j, not two_j.
 */

/* Only call this once, used by all threads.  Returns allocated size,
 * useful for statistics purposes. */
size_t wigxjpf_fill_factors(int max_factorial); /* 0 to free memory */

struct wigxjpf_temp;

/* One temporary array is needed per thread (when used simultaneously). */
struct wigxjpf_temp *wigxjpf_temp_alloc(int max_iter);

void wigxjpf_temp_free(struct wigxjpf_temp *temp);

size_t wigxjpf_temp_size(struct wigxjpf_temp *temp);

int trivial_zero_3j(int two_j1, int two_j2, int two_j3,
		    int two_m1, int two_m2, int two_m3);
int trivial_zero_6j(int two_j1, int two_j2, int two_j3,
		    int two_j4, int two_j5, int two_j6);
int trivial_zero_9j(int two_j1, int two_j2, int two_j3,
		    int two_j4, int two_j5, int two_j6,
		    int two_j7, int two_j8, int two_j9);

void calc_3j_double(double *result,
		    int two_j1, int two_j2, int two_j3,
		    int two_m1, int two_m2, int two_m3,
		    struct wigxjpf_temp *temp);
void calc_6j_double(double *result,
		    int two_j1, int two_j2, int two_j3,
		    int two_j4, int two_j5, int two_j6,
		    struct wigxjpf_temp *temp);
void calc_9j_double(double *result,
		    int two_j1, int two_j2, int two_j3,
		    int two_j4, int two_j5, int two_j6,
		    int two_j7, int two_j8, int two_j9,
		    struct wigxjpf_temp *temp);

void calc_3j_long_double(long double *result,
			 int two_j1, int two_j2, int two_j3,
			 int two_m1, int two_m2, int two_m3,
			 struct wigxjpf_temp *temp);
void calc_6j_long_double(long double *result,
			 int two_j1, int two_j2, int two_j3,
			 int two_j4, int two_j5, int two_j6,
			 struct wigxjpf_temp *temp);
void calc_9j_long_double(long double *result,
			 int two_j1, int two_j2, int two_j3,
			 int two_j4, int two_j5, int two_j6,
			 int two_j7, int two_j8, int two_j9,
			 struct wigxjpf_temp *temp);

/**********************************************************************/

#ifdef __cplusplus
}
#endif

#endif/*__WIGXJPF_H__*/
