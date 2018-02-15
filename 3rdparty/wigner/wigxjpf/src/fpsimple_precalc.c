
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

#include "multi_word_int.h"

double factorial_precalc[FPSIMPLE_MAX_FACTORIAL+1];
int max_factorial_precalc = -1;

double factorial_log_precalc[FPSIMPLE_MAX_FACTORIAL+1];
int max_factorial_log_precalc = -1;

void fpsimplexj_setup_factorial_precalc()
{
  int i;

  factorial_precalc[0] = 1;
  max_factorial_precalc = 0;

  struct multi_word_int big_fact;

  mwi_alloc(&big_fact);
  mwi_set_one_word(&big_fact, 1);

  for (i = 1; i <= FPSIMPLE_MAX_FACTORIAL; i++)
    {
      mwi_mul_plain(&big_fact, (uint32_t) i);

      long double value;
      int expo;

      mwi_to_long_double(&value, &expo, &big_fact);

      factorial_precalc[i] = ldexp((double) value, expo);
      factorial_log_precalc[i] = (double) logl(ldexpl(value, expo));

      if (isfinite(factorial_precalc[i]))
	max_factorial_precalc = i;

      max_factorial_log_precalc = i;

      /* printf ("%d  %.1f\n", i, factorial_precalc[i]); */
    }
  mwi_free(&big_fact);
}
