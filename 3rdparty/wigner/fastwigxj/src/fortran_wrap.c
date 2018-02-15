
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

#include "fastwigxj.h"
#include "fastwigxj_config.h"

#include <string.h>

void ffastwigxj_load_(const char *filename, int32_t *type,
		      int32_t filename_len)
{
  char *filename_arg;

  filename_arg = malloc((size_t) filename_len + 1);

  if (filename_arg == NULL)
    {
      fprintf (stderr,
	       "Failure allocating memory "
	       "for temporary filename (%zd bytes).\n",
	       (size_t) filename_len);
      exit(1);
    }
  memcpy(filename_arg, filename, (size_t) filename_len);
  filename_arg[filename_len] = 0;
  
  fastwigxj_load(filename_arg, *type, NULL);

  free(filename_arg);
}

void ffastwigxj_unload_(int32_t *type)
{
  fastwigxj_unload(*type);
}

void ffastwigxj_print_stats_()
{
  fastwigxj_print_stats();
}

/* 3j */

void ffw3jj_canon_(const int32_t *two_jv, int64_t *rx)
{
  fw3jj_canon(two_jv, (uint64_t *) rx);
}

void ffw3jj_prefetch_(uint64_t x)
{
  fw3jj_prefetch((uint64_t) x);
}

double ffw3jj_get_(const int32_t *two_jv, uint64_t x)
{
  return fw3jj_get((const int *) two_jv, (uint64_t) x);
}

double ffw3jjl_(const int32_t *two_jv)
{
  return fw3jjl((const int *) two_jv);
}

double ffw3jja_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		int32_t *two_m1, int32_t *two_m2)
{
  return fw3jja(*two_j1, *two_j2, *two_j3,
		*two_m1, *two_m2);
}

double ffw3jja6_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		 int32_t *two_m1, int32_t *two_m2, int32_t *two_m3)
{
  return fw3jja6(*two_j1, *two_j2, *two_j3,
		 *two_m1, *two_m2, *two_m3);
}

/* 6j */

void ffw6jj_canon_(const int32_t *two_jv, int64_t *rx)
{
  fw6jj_canon(two_jv, (uint64_t *) rx);
}

void ffw6jj_prefetch_(uint64_t x)
{
  fw6jj_prefetch((uint64_t) x);
}

double ffw6jj_get_(const int32_t *two_jv, uint64_t x)
{
  return fw6jj_get((const int *) two_jv, (uint64_t) x);
}

double ffw6jjl_(const int32_t *two_jv)
{
  return fw6jjl((const int *) two_jv);
}

double ffw6jja_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		int32_t *two_j4, int32_t *two_j5, int32_t *two_j6)
{
  return fw6jja(*two_j1, *two_j2, *two_j3,
		*two_j4, *two_j5, *two_j6);
}

/* 9j */

void ffw9jj_canon_(const int32_t *two_jv, int64_t *rkey, int64_t *rx)
{
  fw9jj_canon(two_jv, (uint64_t *) rkey, (uint64_t *) rx);
}

void ffw9jj_prefetch_(int64_t *x)
{
  fw9jj_prefetch((uint64_t) x);
}

double ffw9jj_get_(const int32_t *two_jv, uint64_t key, uint64_t x)
{
  return fw9jj_get((const int *) two_jv, (uint64_t) key, (uint64_t) x);
}

double ffw9jjl_(const int32_t *two_jv)
{
  return fw9jjl((const int *) two_jv);
}

double ffw9jja_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		int32_t *two_j4, int32_t *two_j5, int32_t *two_j6,
		int32_t *two_j7, int32_t *two_j8, int32_t *two_j9)
{
  return fw9jja(*two_j1, *two_j2, *two_j3,
		*two_j4, *two_j5, *two_j6,
		*two_j7, *two_j8, *two_j9);
}
