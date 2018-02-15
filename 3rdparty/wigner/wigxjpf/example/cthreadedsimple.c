#include <stdio.h>
#include <time.h>
#include <math.h>

#include "wigxjpf.h"

/* This program tests to use the wigner calculation functions from
 * multiple threads.  The threads are forked using openmp.
 *
 * Note that wig_thread_temp_init() and wig_temp_free() is called from
 * each thread.  (wig_thread_temp_init() is the same as wig_temp_init(),
 * but ensures that __thread was applied during library compilation.)
 *
 * wig_table_init() and wig_table_free() must be called globally,
 * i.e. only for one thread.
 *
 * For each iteration, compare the result of two calculations which
 * shall yield the same result.
 */

int main()
{
  wig_table_init(2*100,9);

#pragma omp parallel
  {
    time_t start = time(NULL);

    int two_j1, two_j2, two_j3;
    int two_m1, two_m2, two_m3;

    double v1, v2, d;

    wig_thread_temp_init(2*100);

    for (two_j1 = 0; two_j1 < 25; two_j1++)
      for (two_j2 = 0; two_j2 < 25; two_j2++)
	for (two_j3 = abs(two_j1 - two_j2);
	     two_j3 <= two_j1 + two_j2; two_j3 += 2)
	  {
	    for (two_m1 = -two_j1; two_m1 <= two_j1; two_m1 += 2)
	      for (two_m2 = -two_j2; two_m2 <= two_j2; two_m2 += 2)
		{
		  two_m3 = -(two_m1 + two_m2);

		  v1 = wig3jj(two_j1, two_j2, two_j3,
			      two_m1, two_m2, two_m3);

		  v2 = wig3jj(two_j2, two_j3, two_j1,
			      two_m2, two_m3, two_m1);

		  d = fabs(v1 - v2);

		  if (d > fabs(v1 + v2) * 1e-13)
		    {
		      fprintf (stderr, "DIFF: %.20g %.20g\n", v1, v2);
		      exit(1);
		    }
		}

	    time_t end = time(NULL);

	    if (end < start || end > start + 1)
	      break;
	  }

    wig_temp_free();
  }
  wig_table_free();

  printf ("WIGXJPF openmp threaded C test program completed.\n");

  return 0;
}
