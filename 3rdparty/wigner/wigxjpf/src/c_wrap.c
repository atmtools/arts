
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

#include <stdio.h>

#include "c_wrap.h"
#include "wigxjpf_error.h"
#include "calc.h"

#if WIGXJPF_HAVE_THREAD
__thread struct wigxjpf_temp *wigxjpf_global_temp = NULL;
#else
struct wigxjpf_temp *wigxjpf_global_temp = NULL;
#endif

void wig_table_init(int max_two_j, int wigner_type)
{
  int max_factorial;

  PYWIGXJPF_ERROR_SETUP_void;

  if (max_two_j < 0) {
    fprintf (stderr,
	     "wigxjpf: Negative max_two_j in wig_table_init.\n");
    wigxjpf_error();
  }

  switch (wigner_type)
    {
    case 0: max_factorial = max_two_j; break;
    case 3: max_factorial = (3 * max_two_j) / 2 + 1; break;
    case 6: max_factorial = (4 * max_two_j) / 2 + 1; break;
    case 9: max_factorial = (5 * max_two_j) / 2 + 1; break;
    default:
      fprintf (stderr,
	       "wigxjpf: Bad wigner_type (%d) in wig_table_init.\n",
	       wigner_type);
      wigxjpf_error();
      /* A cleaner solution to the compiler warning about later
       * uninitialised use of max_factorial would be to have a
       * noreturn attribute on wigxjpf_error().  But then we may end
       * up with compiler issues if they do not understand attributes?
       */
      return; /* In fact, we never return... */
    }

  wigxjpf_fill_factors(max_factorial);
}

void wig_table_free(void)
{
  PYWIGXJPF_ERROR_SETUP_void;

  wigxjpf_fill_factors(0);
}

void wig_temp_init(int max_two_j)
{
  int max_iter = (max_two_j / 2) + 1;

  PYWIGXJPF_ERROR_SETUP_void;

  if (max_two_j < 0) {
    fprintf (stderr, "wigxjpf: Negative max_two_j in wig_temp_init.\n");
    wigxjpf_error();
  }

  wigxjpf_global_temp = wigxjpf_temp_alloc(max_iter);
}

#if WIGXJPF_HAVE_THREAD
void wig_thread_temp_init(int max_two_j)
{
  PYWIGXJPF_ERROR_SETUP_void;

  wig_temp_init(max_two_j);
}
#endif

void wig_temp_free(void)
{
  PYWIGXJPF_ERROR_SETUP_void;

  wigxjpf_temp_free(wigxjpf_global_temp);
  wigxjpf_global_temp = NULL;
}

#if PYWIGXJPF_ERROR_HANDLING
/* This function is only called from the error handler, to allow reuse
 * of the temp array.
 */
void wigxjpf_drop_temp(void)
{
  if (wigxjpf_global_temp)
    wigxjpf_global_temp->inuse = 0;
}
#endif

double wig3jj(int two_j1, int two_j2, int two_j3,
	      int two_m1, int two_m2, int two_m3)
{
  double result;

  PYWIGXJPF_ERROR_SETUP_NaN;

  calc_3j_double(&result,
		 two_j1, two_j2, two_j3,
		 two_m1, two_m2, two_m3,
		 wigxjpf_global_temp);

  return result;
}

double wig6jj(int two_j1, int two_j2, int two_j3,
	      int two_j4, int two_j5, int two_j6)
{
  double result;

  PYWIGXJPF_ERROR_SETUP_NaN;

  calc_6j_double(&result,
		 two_j1, two_j2, two_j3,
		 two_j4, two_j5, two_j6,
		 wigxjpf_global_temp);

  return result;
}

double wig9jj(int two_j1, int two_j2, int two_j3,
	      int two_j4, int two_j5, int two_j6,
	      int two_j7, int two_j8, int two_j9)
{
  double result;

  PYWIGXJPF_ERROR_SETUP_NaN;

  calc_9j_double(&result,
		 two_j1, two_j2, two_j3,
		 two_j4, two_j5, two_j6,
		 two_j7, two_j8, two_j9,
		 wigxjpf_global_temp);

  return result;
}
