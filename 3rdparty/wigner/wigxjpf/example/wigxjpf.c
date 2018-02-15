
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

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <inttypes.h>
#include <sys/time.h>

#include "wigxjpf.h"

#if 1
# define DOUBLE_TYPE         double
# define DBL_FCN_POSTFIX(x)  x##_double
# define FMT_DBL_g           "g"
#else
# define DOUBLE_TYPE         long double
# define DBL_FCN_POSTFIX(x)  x##_long_double
# define FMT_DBL_g           "Lg"
#endif

struct wigxjpf_temp *temp = NULL;
size_t factors_size;
int showstat = 0;

void alloc_temps(int maxfact, int maxiter)
{
  factors_size = wigxjpf_fill_factors(maxfact);
  temp = wigxjpf_temp_alloc(maxiter);
}

void free_temps()
{
  if (!temp)
    return;

  if (showstat)
    {
      size_t temp_size = wigxjpf_temp_size(temp);

      printf ("Memory used (bytes): %zd + %zd = %zd\n",
	      factors_size, temp_size, factors_size + temp_size);
    }

  wigxjpf_fill_factors(0);
  wigxjpf_temp_free(temp);  temp = NULL;
}

#define DBL_TO_TWO(j) do {				\
    two_##j = (int) ((j) * 2.);				\
    if (two_##j != (j) * 2.) {				\
      fprintf (stderr, "Given j (%f) in option '%s' "	\
	       "is neither integer nor "		\
	       "half-integer.\n", j, argv[i]);		\
      exit(1);						\
    }							\
    if (two_##j > max_two_j) {				\
      max_two_j = two_##j;				\
    }							\
  }							\
  while (0)

void usage(const char *argv0)
{
  printf ("\n");
  printf ("WIGXJPF - Wigner 3j, 6j, 9j evaluator. \n");
  printf ("\n");
  printf ("Usage: %s [OPTIONS]\n", argv0);
  printf ("\n");
  printf ("  --3j=j1,j2,j3,m1,m2,m3           Evaluate Wigner 3j symbol.\n");
  printf ("  --6j=j1,j2,j3,j4,j5,j6           Evaluate Wigner 6j symbol.\n");
  printf ("  --9j=j1,j2,j3,j4,j5,j6,j7,j8,j9  Evaluate Wigner 9j symbol.\n");
  printf ("  --stats    Print memory and time usage after each evaluation.\n");
  printf ("  --fact=N   Set up dummy prime factorised factorial tables.\n");
  printf ("  --help     Print this usage information and quit.\n");
  printf ("\n");
}

int main(int argc, char *argv[])
{
  int i;
  int badopt;
  struct timeval t1;
  struct timeval t2;
  double elapsed;

  gettimeofday(&t1, NULL);

  for (i = 1; i < argc; i++)
    {
      badopt = 0;

      gettimeofday(&t1, NULL);

      if (strcmp(argv[i],"--help") == 0)
	{
	  usage(argv[0]);
	  exit(0);
	}
      else if (strcmp(argv[i],"--stats") == 0)
	{
	  showstat = 2;

	  for ( ; ; )
	    {
	      gettimeofday(&t2, NULL);
	      elapsed =
		(double) (t2.tv_sec - t1.tv_sec) +
		1.e-6 * (double) (t2.tv_usec - t1.tv_usec);
	      if (elapsed < 0 || elapsed > 0.1)
		break;
	    }
	}
      else if (strncmp(argv[i],"--fact=",7) == 0)
	{
	  int max_fact = atoi(argv[i]+7);

	  alloc_temps(max_fact,0);

	  free_temps();
	}
      else if (strncmp(argv[i],"--3j=",5) == 0)
	{
	  double j1, j2, j3, m1, m2, m3;
	  int two_j1, two_j2, two_j3,
	    two_m1, two_m2, two_m3;
	  int max_two_j = 0;
	  int n;

	  n = sscanf(argv[i]+5,"%lf,%lf,%lf,%lf,%lf,%lf",
		     &j1, &j2, &j3, &m1, &m2, &m3);

	  if (n == 6)
	    {
	      DBL_TO_TWO(j1);
	      DBL_TO_TWO(j2);
	      DBL_TO_TWO(j3);
	      DBL_TO_TWO(m1);
	      DBL_TO_TWO(m2);
	      DBL_TO_TWO(m3);

	      printf ("3j(%.1f %.1f %.1f   %.1f %.1f %.1f) = ",
		      0.5 * two_j1, 0.5 * two_j2, 0.5 * two_j3,
		      0.5 * two_m1, 0.5 * two_m2, 0.5 * two_m3);

	      if (trivial_zero_3j(two_j1, two_j2, two_j3,
				  two_m1, two_m2, two_m3))
		printf ("trivially 0\n");
	      else
		{
		  DOUBLE_TYPE result;

		  alloc_temps((3 * max_two_j) / 2 + 1,
			      max_two_j / 2 + 1);

		  DBL_FCN_POSTFIX(calc_3j)(&result,
					   two_j1, two_j2, two_j3,
					   two_m1, two_m2, two_m3,
					   temp);

		  if (result == 0.0)
		    printf ("0\n");
		  else
		    printf ("%#.15" FMT_DBL_g "\n", result);

		  free_temps();
		}
	    }
	  else
	    badopt = 1;
	}
      else if (strncmp(argv[i],"--6j=",5) == 0)
	{
	  double j1, j2, j3, j4, j5, j6;
	  int two_j1, two_j2, two_j3,
	    two_j4, two_j5, two_j6;
	  int max_two_j = 0;
	  int n;

	  n = sscanf(argv[i]+5,"%lf,%lf,%lf,%lf,%lf,%lf",
		     &j1, &j2, &j3, &j4, &j5, &j6);

	  if (n == 6)
	    {
	      DBL_TO_TWO(j1);
	      DBL_TO_TWO(j2);
	      DBL_TO_TWO(j3);
	      DBL_TO_TWO(j4);
	      DBL_TO_TWO(j5);
	      DBL_TO_TWO(j6);

	      printf ("6j(%.1f %.1f %.1f   %.1f %.1f %.1f) = ",
		      0.5 * two_j1, 0.5 * two_j2, 0.5 * two_j3,
		      0.5 * two_j4, 0.5 * two_j5, 0.5 * two_j6);

	      if (trivial_zero_6j(two_j1, two_j2, two_j3,
				  two_j4, two_j5, two_j6))
		printf ("trivially 0\n");
	      else
		{
		  DOUBLE_TYPE result;

		  alloc_temps((4 * max_two_j) / 2 + 1,
			      max_two_j / 2 + 1);

		  DBL_FCN_POSTFIX(calc_6j)(&result,
					   two_j1, two_j2, two_j3,
					   two_j4, two_j5, two_j6,
					   temp);

		  if (result == 0.0)
		    printf ("0\n");
		  else
		    printf ("%#.15" FMT_DBL_g "\n", result);

		  free_temps();
		}
	    }
	  else
	    badopt = 1;
	}
      else if (strncmp(argv[i],"--9j=",5) == 0)
	{
	  double j1, j2, j3, j4, j5, j6, j7, j8, j9;
	  int two_j1, two_j2, two_j3,
	    two_j4, two_j5, two_j6,
	    two_j7, two_j8, two_j9;
	  int max_two_j = 0;
	  int n;

	  n = sscanf(argv[i]+5,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		     &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9);

	  if (n == 9)
	    {
	      DBL_TO_TWO(j1);
	      DBL_TO_TWO(j2);
	      DBL_TO_TWO(j3);
	      DBL_TO_TWO(j4);
	      DBL_TO_TWO(j5);
	      DBL_TO_TWO(j6);
	      DBL_TO_TWO(j7);
	      DBL_TO_TWO(j8);
	      DBL_TO_TWO(j9);

	      printf ("9j(%.1f %.1f %.1f   %.1f %.1f %.1f   "
		      "%.1f %.1f %.1f) = ",
		      0.5 * two_j1, 0.5 * two_j2, 0.5 * two_j3,
		      0.5 * two_j4, 0.5 * two_j5, 0.5 * two_j6,
		      0.5 * two_j7, 0.5 * two_j8, 0.5 * two_j9);

	      if (trivial_zero_9j(two_j1, two_j2, two_j3,
				  two_j4, two_j5, two_j6,
				  two_j7, two_j8, two_j9))
		printf ("trivially 0\n");
	      else
		{
		  DOUBLE_TYPE result;

		  alloc_temps((5 * max_two_j) / 2 + 1,
			      max_two_j / 2 + 1);

		  DBL_FCN_POSTFIX(calc_9j)(&result,
					   two_j1, two_j2, two_j3,
					   two_j4, two_j5, two_j6,
					   two_j7, two_j8, two_j9,
					   temp);

		  if (result == 0.0)
		    printf ("0\n");
		  else
		    printf ("%#.15" FMT_DBL_g "\n", result);

		  free_temps();
		}
	    }
	  else
	    badopt = 1;
	}

      gettimeofday(&t2, NULL);

      if (showstat == 1)
	{
	  elapsed =
	    (double) (t2.tv_sec - t1.tv_sec) +
	    1.e-6 * (double) (t2.tv_usec - t1.tv_usec);

	  printf ("Elapsed time: %.6f s.\n", elapsed);
	}
      showstat = !!showstat;

      if (badopt)
	{
	  fprintf (stderr, "Unknown option '%s'.\n", argv[i]);
	  exit (1);
	}
    }

  if (argc <= 1)
    usage(argv[0]);
  exit(0);
}
