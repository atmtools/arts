
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of FASTWIGXJ.
 *
 *  FASTWIGXJ is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  FASTWIGXJ is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with FASTWIGXJ.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#ifndef __FASTWIGXJ_H__
#define __FASTWIGXJ_H__

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct fastwigxj_header;

size_t fastwigxj_load(const char *filename, int type,
		      struct fastwigxj_header *header);
void fastwigxj_unload(int type);

size_t fastwigxj_dyn_init(int type, size_t entries);

void fastwigxj_dyn_free(int type);

void fastwigxj_print_stats();

/* 3j */

static void fw3jj_canon(const int *two_jv, uint64_t *rx);

static void fw3jj_prefetch(uint64_t x);

static double fw3jj_get(const int *two_jv, uint64_t x);

static double fw3jjl(const int *two_jv);

static double fw3jja(int two_j1, int two_j2, int two_j3,
		     int two_m1, int two_m2);

/* TODO: remove me? */
static double fw3jja6(int two_j1, int two_j2, int two_j3,
		      int two_m1, int two_m2, int two_m3);

/* 6j */

static void fw6jj_canon(const int *two_jv, uint64_t *rx);

static void fw6jj_prefetch(uint64_t x);

static double fw6jj_get(const int *two_jv, uint64_t x);

static double fw6jjl(const int *two_jv);

static double fw6jja(int two_j1, int two_j2, int two_j3,
		     int two_j4, int two_j5, int two_j6);

/* 9j */

static void fw9jj_canon(const int *two_jv, uint64_t *rkey, uint64_t *rx);

static void fw9jj_prefetch(uint64_t x);

static double fw9jj_get(const int *two_jv, uint64_t key, uint64_t x);

static double fw9jjl(const int *two_jv);

static double fw9jja(int two_j1, int two_j2, int two_j3,
		     int two_j4, int two_j5, int two_j6,
		     int two_j7, int two_j8, int two_j9);

#ifdef __cplusplus
}
#endif

/* Include inline function declarations. */

#include "fastwigxj_inc.h"

#endif/*__FASTWIGXJ_H__*/
